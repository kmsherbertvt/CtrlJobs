module Adaptations
    import Random
    import LinearAlgebra

    import CtrlVQE

    import ..StateVars
    import ..Float
    import ..require_work

    #= Candidate generation.

    Each of these functions takes the active vars and the gradient signals,
        and returns a list of candidate split times,
        represented with the tuple `(i,s)` where i indexes the drive.

    Note these functions are responsible for filtering out candidates
        which would violate fMAX constraint.
    =#

    """ Split existing windows in half. """
    function makepool_bisection(vars, ϕ)
        T = vars.setup.T
        sMIN = 1 / vars.setup.fMAX
        s = vars.state.s

        pool = Tuple{Int,Float}[]
        for i in eachindex(s)
            Δs = diff([s[i]; T]) ./ 2
            candidates = s[i] .+ Δs

            # FILTER OUT BISECTIONS RESULTING IN TOO SMALL A WINDOW
            toosmall = findall(Δs .< sMIN)
            isnothing(toosmall) || deleteat!(candidates, toosmall)
            append!(pool, (i, candidate) for candidate in candidates)
        end

        return frequencyfilter!(vars, pool)
    end

    # """ Pick a random time for each pulse, but ensure it is valid under fMAX. """
    # function makepool_random(vars, ϕ)
    #     s̄ = vars.state.s
    #     ΔsMIN = 1 / vars.setup.fMAX
    #     T = vars.setup.T

    #     #= This is a little weird.
    #         I want to pick a single random time for each drive,
    #         sampled uniformly from the parts of the number line within ±ΔsMIN from each s.
    #     =#
    #     pool = Tuple{Int,Float}[]
    #     for i in 1:CtrlVQE.Devices.ndrives(vars.work.device)
    #         # IDENTIFY VALID TIME RANGES
    #         intervals = Pair{Float,Float}[]
    #         for k in eachindex(s̄[i])
    #             first = s̄[i][k] + ΔsMIN
    #             second = (k == length(s̄[i]) ? T : s̄[i][k+1]) - ΔsMIN
    #             first >= second && continue
    #             push!(intervals, first => second)
    #         end
    #         isempty(intervals) && continue

    #         # MAP VALID TIME RANGES ONTO A STANDARD INTERVAL [0,t)
    #         Δs_avail = [interval.second - interval.first for interval in intervals]
    #         Δs_cum = cumsum(Δs_avail)

    #         # GENERATE A RANDOM NUMBER - note we rely on script to control the seed
    #         r = rand(Float) * last(Δs_cum)

    #         # FIND WHICH SUB-INTERVAL THIS RANDOM NUMBER FALLS
    #         k = findfirst(Δs -> Δs > r, Δs_cum)

    #         # SLIDE RANDOM NUMBER OVER SO IT IS THE Δ WITH RESPECT TO INTERVAL START TIME
    #         k > 1 && (r -= Δs_cum[k-1])

    #         # ADD THE SELECTION TO THE POOL
    #         push!(pool, (i, intervals[k].first + r))
    #     end
    #     return pool
    # end

    """ Pick a random time for each window long enough to be split. """
    function makepool_random(vars, ϕ)
        s̄ = vars.state.s
        ΔsMIN = 1 / vars.setup.fMAX
        T = vars.setup.T

        pool = Tuple{Int,Float}[]
        for i in 1:CtrlVQE.Devices.ndrives(vars.work.device)
            # IDENTIFY VALID TIME RANGES
            intervals = Pair{Float,Float}[]
            for k in eachindex(s̄[i])
                first = s̄[i][k] + ΔsMIN
                second = (k == length(s̄[i]) ? T : s̄[i][k+1]) - ΔsMIN
                first >= second && continue
                push!(intervals, first => second)
            end
            isempty(intervals) && continue

            # PICK A RANDOM TIME IN THE MARGINED INTERVAL
            for interval in intervals
                r = rand(Float)
                range = interval.second - interval.first
                s = (r * range) + interval.first
                push!(pool, (i, s))
            end
        end
        return pool
    end

    """ Identify nodes in each gradient signal. """
    function makepool_nodes(vars, ϕ)
        t = CtrlVQE.lattice(vars.work.grid)

        pool = Tuple{Int,Float}[]
        for i in 1:CtrlVQE.Devices.ndrives(vars.work.device)
            j = 2i - 1
            ϕα = @view(ϕ[:,j])
            ϕβ = @view(ϕ[:,j+1])

            append!(pool, (i,s) for s in identifynodes(t, ϕα))
            append!(pool, (i,s) for s in identifynodes(t, ϕβ))
        end

        return frequencyfilter!(vars, pool)
    end

    """ Choose an optimal split-point for each eligible window,
            based on combination of gradients.

    Implementation assumes each window has two adjacent parameters α and β!

    A more sophisticated version might do some sort of uni-variate optimization
        over a loss function defined in terms of continuous-valued ϕα and ϕβ
        inferred from an interpolation over the relevant columns in ϕ.
    That sounds really fun.
    But for now, let's stick to the time grid used in simulation,
        and just exhaustively search for which time gives the best integrals. :)

    TODO:
        Brute force consistently picks one time-step later.
        So maybe the more sophisticated version is worthwhile.

    """
    function makepool_optimal(vars, ϕ)
        s̄ = vars.state.s
        ΔsMIN = 1 / vars.setup.fMAX
        T = vars.setup.T

        t = CtrlVQE.lattice(vars.work.grid)
        τ = CtrlVQE.stepsize(vars.work.grid)

        pool = Tuple{Int,Float}[]
        for i in 1:CtrlVQE.Devices.ndrives(vars.work.device)
            j = 2i - 1
            ϕα = @view(ϕ[:,j])
            ϕβ = @view(ϕ[:,j+1])

            # INITIALIZE INTEGRALS FOR EACH TIME CANDIDATE
            L = zero(t)

            # IDENTIFY VALID TIME RANGES
            intervals = Pair{Float,Float}[]
            for k in eachindex(s̄[i])
                first = s̄[i][k]
                second = (k == length(s̄[i]) ? T : s̄[i][k+1])
                push!(intervals, first => second)
            end

            # IDENTIFY THE OPTIMAL TIME WITHIN EACH TIME RANGE
            for interval in intervals

                # START AND END INDICES, FOR INTEGRATION
                i1 = findfirst(t_ -> t_ >= interval.first, t)
                i2 = findlast( t_ -> t_ < interval.second, t)

                # COMPUTE INTEGRALS FOR ALL t IN INTERVAL
                I_α = cumsum(ϕα[i1:i2]) .* τ
                I_β = cumsum(ϕβ[i1:i2]) .* τ
                L[i1:i2] .= (I_α .^ 2) .+ (I_β .^ 2)

                # ERASE L FOR ANY t VIOLATING BANDWIDTH CONSTRAINTS
                iL = findfirst(t_ -> t_ >= interval.first + ΔsMIN, t)
                iR = findlast( t_ -> t_ < interval.second - ΔsMIN, t)
                L[i1:iL] .= 0
                L[iR:i2] .= 0

                # IDENTIFY THE TIME YIELDING THE MAXIMUM LOSS FUNCTION
                LMAX, iO = findmax(L[i1:i2])
                √(2LMAX) ≤ vars.setup.g_tol && continue         # TOO SMALL TO BOTHER
                    # Factor of 2 is the contribution from the other half of the window.
                    # √(2LMAX) is the 2-norm of the gradient, basically.
                    #   g_gol is an ∞-norm. Not *exactly* comparable but meh.
                tO = t[i1:i2][iO]

                # REGISTER THIS CANDIDATE
                push!(pool, (i, tO))
            end
        end

        return pool
    end

    """ Literally every time in the time grid. >_> """
    function makepool_brute(vars, ϕ)
        t = CtrlVQE.lattice(vars.work.grid)
        pool = Tuple{Int,Float}[]
        for i in 1:CtrlVQE.Devices.ndrives(vars.work.device)
            append!(pool, (i, t_) for t_ in t)
        end
        return frequencyfilter!(vars, pool)
    end



    """ Helper method to remove segments violating bandwidth constraints. """
    function frequencyfilter!(vars, pool)
        s̄ = vars.state.s
        ΔsMIN = 1 / vars.setup.fMAX
        T = vars.setup.T

        toosmall = Int[]
        for ix in eachindex(pool)
            i, sC = pool[ix]

            # FIND CURRENT WINDOW SEGMENTS BEFORE AND AFTER THE CANDIDATE
            iL = findlast(s -> s <= sC, s̄[i])
            isnothing(iL) && (push!(toosmall, ix); continue)
            sL = s̄[i][iL]
            sR = iL == length(s̄[i]) ? T : s̄[i][iL+1]

            # FLAG CANDIDATE AS TO-BE-DISCARDED IF IT ISN'T FAR ENOUGH AWAY
            ((sC - sL) < ΔsMIN || (sR - sC) < ΔsMIN) && (push!(toosmall, ix))
        end
        return deleteat!(pool, toosmall)
    end

    """ Helper method to linearly extrapolate nodes in any time series. """
    function identifynodes(t̄, ϕ)
        s̄ = Float[]
        isempty(ϕ) && return s̄

        phase = sign(first(ϕ))          # TRACK IF WE ARE ABOVE OR BELOW ϕ=0 LINE
        for (i, ϕi) in enumerate(ϕ)
            sign(ϕi) == phase && continue   # SKIPS FIRST ITERATION, SO ϕ̄[i-1] IS SAFE

            # LINEAR INTERPOLATION TO IDENTIFY TIME OF PHASE FLIP
            t = t̄[i]; t0 = t̄[i-1]; ϕ0 = ϕ[i-1]
            m = abs(ϕi-ϕ0)/(t-t0)           #     SLOPE OF LINE
            s = t0 + m*ϕ0                   # INTERCEPT OF LINE
            push!(s̄, s)

            phase = sign(ϕi)                # UPDATE WHICH SIDE OF ϕ=0 WE ARE ON
        end

        return s̄
    end



    #= Candidate selection.

    Each of these functions takes a list of candidates and a list of scores,
        and returns one candidate or a list of candidates.

    =#

    function select_one(vars, pool, scores)
        return pool[argmax(scores)]
    end

    function select_oneperpulse(vars, pool, scores)
        segments = Dict{Int,Int}()  # MAPS DRIVE INDEX TO POOL INDEX
        for ix in eachindex(pool)
            # FETCH THE DRIVE INDEX FOR THIS SEGMENT
            i, _ = pool[ix]
            # REGISTER INDEX IF IT HAS THE BEST SCORE SO FAR FOR THIS DRIVE
            if !(i in keys(segments)) || scores[segments[i]] < scores[ix]
                segments[i] = ix
            end
        end
        return [pool[ix] for ix in values(segments)]
    end

    function select_oneperwindow(vars, pool, scores)
        s̄ = vars.state.s

        segments = Dict{Tuple{Int,Int},Int}()  # MAPS (DRIVE,WINDOW) TO POOL INDEX
        for ix in eachindex(pool)
            # FETCH THE DRIVE INDEX FOR THIS SEGMENT
            i, s = pool[ix]
            # FETCH THE WINDOW INDEX FOR THIS SEGMENT
            k = findlast(s̄[i] .≤ s)
            isnothing(k) && continue    # IGNORE IF CANDIDATE IS NEGATIVE
                # NOTE: Only relevant in an annoying edge case of linear extrapolation.

            if !((i,k) in keys(segments)) || scores[segments[(i,k)]] < scores[ix]
                segments[i,k] = ix
            end
        end

        return [pool[ix] for ix in values(segments)]
    end




    #= Hessian update.

    Each of these functions takes the current Hk and an index map
        identifying which new parameters are taken from which old parameters,
        and they return a new Hk on the new parameter space.
    =#

    """ Hessian adaptation: Use identity on whole space. """
    function upHk_naive(Hk, imap)
        L = length(imap)
        return Matrix{eltype(Hk)}(LinearAlgebra.I, L, L)
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
        # TODO: I don't really know.
        #   It should GUARANTEE positive definiteness,
        #       and it should USE information obtained on all parameters.
        #   Somehow.
        error("Dunno how to do this yet.")
    end






    #= Stock behavior. =#

    """ Insert an extra window into a windowed signal, duplicating the sub-signal.

    Also return an index map relating new parameters to old.
    """
    function insert_window(signal, s::Float)
        # IDENTIFY THE WINDOW WHERE s CURRENTLY FALLS
        k = findlast(s0 -> s0 < s, signal.starttimes)
        isnothing(k) && error("`s` does not fall in any window.")
        window = deepcopy(signal.windows[k]) # DUPLICATE PARAMETERS SO Ω(t) IS UNCHANGED

        # CREATE THE NEW PULSE
        windows = deepcopy(signal.windows);         insert!(windows, k+1, window)
        starttimes = deepcopy(signal.starttimes);   insert!(starttimes, k+1, s)

        # MAP NEW PARAMETERS TO OLD PARAMETERS
        Ll = sum(
            CtrlVQE.Parameters.count(window) for window in signal.windows[1:k-1];
            init=0,
        )
        Lk = CtrlVQE.Parameters.count(signal.windows[k])
        Lr = sum(
            CtrlVQE.Parameters.count(window) for window in signal.windows[1+k:end];
            init=0,
        )
        imap = vcat(1:Ll..., 1+Ll:Lk+Ll..., 1+Ll:Lk+Ll..., 1+Ll+Lk:Lr+Ll+Lk)

        return CtrlVQE.WindowedSignal(windows, starttimes), imap
    end

    function insert_window(signal, s̄::Vector{Float})
        latestsignal = Ref(signal)
        latestimap = Ref(collect(1:CtrlVQE.Parameters.count(signal)))
        for s in s̄
            latestsignal[], imap = insert_window(latestsignal[], s)
            latestimap[] = latestimap[][imap]
        end
        return latestsignal[], latestimap[]
    end

    """ Prepare a candidate device which includes an extra window.

    Also return an index map relating new parameters to old.
    """
    function prepare_candidate(vars, segment::Tuple{Int,Float})
        i, s = segment
        candidate = deepcopy(vars.work.device)
        signal = CtrlVQE.drivesignal(candidate, i)
        newsignal, pulse_imap = insert_window(signal, s)
        CtrlVQE.set_drivesignal(candidate, i, newsignal)

        nD = CtrlVQE.ndrives(candidate)
        Ω = [i == i_ ? vars.state.Ω[i][pulse_imap] : vars.state.Ω[i_] for i_ in 1:nD]
        imap = vcat(Ω..., vars.state.ν)

        return candidate, imap
    end

    function prepare_candidate(vars, segments::Vector{Tuple{Int,Float}})
        candidate = deepcopy(vars.work.device)
        Ω = Vector{Int}[]
        for i in 1:CtrlVQE.ndrives(candidate)
            s̄ = [segment[2] for segment in segments if segment[1] == i]
            signal = CtrlVQE.drivesignal(candidate, i)
            newsignal, pulse_imap = insert_window(signal, s̄)
            CtrlVQE.set_drivesignal(candidate, i, newsignal)
        push!(Ω, vars.state.Ω[i][pulse_imap])
        end
        imap = vcat(Ω..., vars.state.ν)
        return candidate, imap
    end

    """ Calculate the total gradient norm for each candidate in a pool. """
    function score_candidates(vars, pool, ϕ)
        lossfn = vars.work.lossfn

        scores = Vector{Float}(undef, length(pool))
        for (ix, segment) in enumerate(pool)
            # PREPARE CANDIDATE DEVICE AND EXTENDED PARAMETER VECTOR
            candidate, _ = prepare_candidate(vars, segment)
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
            scores[ix] = LinearAlgebra.norm(g)
                #= NOTE: We get different behavior if we use ∞-norm.
                    The 2-norm makes more sense, but it's not necessarily better. =#
        end
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
        push!(vars.trace.Δs_min, [minimum(diff([s; vars.setup.T])) for s in vars.state.s])
    end

    """ Update state objects with one or more new segments.

    The `segments` argument may be a tuple or a vector of tuples `(i,s)`,
        where `i` is the drive index and `s` is the starttime of a new window.
    """
    function adapt_state!(vars, segments, upHk)
        # PREPARE A CANDIDATE DEVICE (which knows all about its own parameters)
        candidate, imap = prepare_candidate(vars, segments)

        # FETCH THE NEW PARAMETER VECTOR
        x = CtrlVQE.Parameters.values(candidate)
        L = CtrlVQE.Parameters.count(candidate)

        # UPDATE THE INVERSE HESSIAN
        Hk = upHk(vars.state.Hk, imap)

        # FETCH WINDOW START TIMES
        n = CtrlVQE.ndrives(candidate)
        signals = [CtrlVQE.drivesignal(candidate, i) for i in 1:n]
        s = [signal.starttimes for signal in signals]

        # PREPARE SIGNAL INDEXING VECTORS
        LΩ = [CtrlVQE.Parameters.count(signal) for signal in signals]
        Ω = [collect(1+sum(LΩ[1:i]):sum(LΩ[1:i+1])) for i in 0:n-1]

        # PREPARE FREQUENCY INDEXING VECTORS -- Implementation may vary with type someday?
        ν = vars.state.ν

        vars.state = StateVars(x, Hk, s, Ω, ν)
    end

    """ Update work objects to match the current state. """
    function adapt_work!(vars)
        require_work(vars)

        # MUTATE THE DEVICE TO ACCOUNT FOR NEW PULSE SHAPES
        for i in 1:CtrlVQE.ndrives(vars.work.device)
            windows = [deepcopy(vars.work.protopulse) for _ in vars.state.s[i]]
            signal = CtrlVQE.WindowedSignal(windows, vars.state.s[i])
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