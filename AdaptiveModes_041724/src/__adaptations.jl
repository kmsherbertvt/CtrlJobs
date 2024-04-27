module Adaptations
    import Random
    import LinearAlgebra

    import CtrlVQE
    include("ModalSignals.jl")
    Harmonic = HarmonicSignals.Harmonic

    import ..StateVars
    import ..Float
    import ..require_work

    #= Candidate generation.

    Each of these functions takes the active vars,
        and returns a list of candidate split times,
        represented with the tuple `(i,s)` where i indexes the drive.

    Note these functions are responsible for filtering out candidates
        which would violate fMAX constraint.
    =#

    function makepool_harmonics(vars)
        T = vars.setup.T
        fMAX = vars.setup.fMAX
        nMAX = floor(Int, fMAX * 2T)
        return [CtrlVQE.ConstrainedSignal(Harmonic(0.0, 0.0, n, T), :T) for n in 1:nMAX]
    end






    #= Candidate selection.

    Each of these functions takes a matrix of scores,
        and returns one candidate or a list of candidates.

    =#

    function select_one(scores)
        i, n = Tuple(argmax(scores))
        return i => n
    end

    function select_oneperpulse(scores)
        return [(
            (i, n) = Tuple(modes);
            i => n
        ) for modes in vec(argmax(scores; dims=2))]
    end




    #= Hessian update.

    Each of these functions takes the current Hk and an index map
        identifying which new parameters are taken from which old parameters,
        and they return a new Hk on the new parameter space.
    =#

    """ Hessian adaptation: Use identity on whole space. """
    function upHk_naive(Hk, imap)
        L = length(imap)
        return Matrix{eltype(Hk)}(I, L, L)
    end

    """ Hessian adaptation: Recylce but identity out the subspace of split parameters. """
    function upHk_slate(Hk, imap)
        # FLAG ANY INDICES IN `imap` WHICH ARE DUPLICATED
        imap = deepcopy(imap)
        for ix in imap
            ix == 0 && continue                         # ALREADY FLAGGED
            dependents = findall(i -> i == ix, imap)
            length(dependents) == 1 && continue         # PARAMETER NOT SPLIT
            imap[dependents] .= 0
        end

        # RECYCLE THE HESSIAN, BUT USE IDENTITY ON THE SUBSPACE OF SPLIT PARAMETERS
        new_Hk = LinearAlgebra.diagm(ones(eltype(Hk), length(imap)))
        for ix in CartesianIndices(Hk); i,j = Tuple(ix)
            i0 = imap[i]; i0 == 0 && continue
            j0 = imap[j]; j0 == 0 && continue
            new_Hk[i,j] = Hk[i0,j0]
        end

        return new_Hk
    end

    function upHk_recycle(Hk, imap)
        #= TODO: I don't really know.
            It should GUARANTEE positive definiteness,
                and it should USE information obtained on all parameters.
            Somehow.
            Mafalda may have insights.
            Not a high priority.
        =#
        error("Dunno how to do this yet.")
    end






    #= Stock behavior. =#

    """ Prepare a candidate device which includes an extra mode.

    Also return an index map relating new parameters to old.
    """
    function prepare_candidate(vars, mode::Pair{Int,Int})
        require_work(vars)
        pool = vars.work.pool
        i, n = mode

        # ADD IN THE NEW MODE
        candidate = deepcopy(vars.work.device)
        push!(
            CtrlVQE.drivesignal(candidate, i).components,
            deepcopy(pool[n]),
        )

        # CONSTRUCT AN INDEX MAP WITH ZEROS FOR THE NEW PARAMETERS
        Ω = deepcopy(vars.state.Ω)
        append!(Ω[i], zeros(CtrlVQE.Parameters.count(pool[n])))
        imap = vcat(Ω..., vars.state.ν)

        return candidate, imap
    end

    function prepare_candidate(vars, pool, modes::Vector{Tuple{Int,Int}})
        #= TODO: Add multiple modes simultaneously.
                The tricky part with square windows was getting the correct index map,
                    but it should be easier with the modal version.
        =#
    end

    """ Calculate the total gradient norm for each candidate in a pool. """
    function score_candidates(vars, ϕ)
        require_work(vars)
        pool = vars.work.pool
        lossfn = vars.work.lossfn

        nD = CtrlVQE.ndrives(vars.work.device)
        scores = Matrix{Float}(undef, nD, length(pool))
        for i in 1:nD; for n in eachindex(pool)
            # PREPARE CANDIDATE DEVICE AND EXTENDED PARAMETER VECTOR
            candidate, _ = prepare_candidate(vars, i => n)
            x = CtrlVQE.Parameters.values(candidate)

            # CALCULATE GRADIENT
            g = CtrlVQE.Devices.gradient(candidate, vars.work.grid, ϕ)
            for l in 1:length(lossfn.penalties)
                # We need to extend the penalty function to work on more parameters.
                # HACK: Construct a new fn using each field except one called "device".
                protofn = lossfn.penaltyfns[l]
                fntype = Base.typename(typeof(protofn)).wrapper
                    # This "typename" "wrapper" strips away parametric args in the type,
                    #   so that Julia can find the generic constructor.
                penaltyfn = fntype((
                    field == :device ? candidate : getfield(protofn, field)
                        for field in fieldnames(fntype)
                )...)

                gd = CtrlVQE.grad_function(penaltyfn)

                # MODIFY THE GRADIENT WITH PENALTY CONTRIBUTION
                g .+= lossfn.weights[l] * gd(x)
            end

            # REGISTER THE GRADIENT NORM
            scores[i,n] = LinearAlgebra.norm(g)
        end; end
        return scores
    end

    #= Update variables with given segment(s). =#

    """ Update trace object with adaptive information. """
    function adapt_trace!(vars, scores)
        G_max = isempty(scores) ? 0.0 : maximum(scores)
        push!(vars.trace.adaptations, last(vars.trace.iterations))
        push!(vars.trace.poolsize, length(scores))
        push!(vars.trace.G_max, G_max)
        push!(vars.trace.parameters, [length(Ω) for Ω in vars.state.Ω])
    end

    """ Update state objects with one or more new segments.

    The `modes` argument may be a tuple or a vector of tuples `(i,n)`,
        where `i` is the drive index and `n` is the mode index of a new signal.
    """
    function adapt_state!(vars, modes, upHk)
        require_work(vars)

        # PREPARE A CANDIDATE DEVICE (which knows all about its own parameters)
        candidate, imap = prepare_candidate(vars, modes)

        # FETCH THE NEW PARAMETER VECTOR
        x = CtrlVQE.Parameters.values(candidate)
        L = CtrlVQE.Parameters.count(candidate)

        # UPDATE THE INVERSE HESSIAN
        Hk = upHk(vars.state.Hk, imap)

        # FETCH WINDOW START TIMES
        n = deepcopy(vars.state.n)
        if modes isa Pair{Int,Int}
            push!(n[modes.first], modes.second)
        else # It's a Vector of such items.
            for (i, n_) in modes
                push!(n[i], n_)
            end
        end

        # PREPARE SIGNAL INDEXING VECTORS
        nD = CtrlVQE.ndrives(candidate)
        signals = [CtrlVQE.drivesignal(candidate, i) for i in 1:nD]
        LΩ = [CtrlVQE.Parameters.count(signal) for signal in signals]
        Ω = [collect(1+sum(LΩ[1:i]):sum(LΩ[1:i+1])) for i in 0:nD-1]

        # PREPARE FREQUENCY INDEXING VECTORS -- Implementation may vary with type someday?
        ν = vars.state.ν

        vars.state = StateVars(x, Hk, n, Ω, ν)
    end

    """ Update work objects to match the current state. """
    function adapt_work!(vars)
        require_work(vars)
        pool = vars.work.pool

        # MUTATE THE DEVICE TO ACCOUNT FOR NEW PULSE SHAPES
        for i in 1:CtrlVQE.ndrives(vars.work.device)
            components = [deepcopy(pool[n]) for n in vars.state.n[i]]
            signal = CtrlVQE.CompositeSignal(deepcopy(vars.work.protopulse), components...)
            CtrlVQE.set_drivesignal(vars.work.device, i, signal)
        end
        CtrlVQE.Parameters.bind(vars.work.device, vars.state.x)
    end

    #= Check trace against metavariables for convergence. =#

    """ Inspect last record of the trace to see if it meets the convergence conditions. """
    function adapt_is_converged(vars)
        isnothing(vars.trace) && return false

        isnothing(vars.trace) && return false
        isempty(vars.trace.adaptations) && return false
        last(vars.trace.G_max) ≤ vars.meta.G_tol && return true

        return false
    end

    """ Inspect last record of the trace to see if it meets the termination conditions. """
    function adapt_is_terminated(vars)
        adapt_is_converged(vars) && return true

        isnothing(vars.trace) && return false
        isempty(vars.trace.adaptations) && return false
        length(vars.trace.adaptations) >= vars.meta.maxadapt && return true

        return false
    end

end