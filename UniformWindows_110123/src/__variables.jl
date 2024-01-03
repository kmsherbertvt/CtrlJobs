module Variables
    import Random
    import LinearAlgebra

    import Parameters: @with_kw

    import CtrlVQE
    import CtrlVQE: SignalType, IntegrationType, DeviceType
    import CtrlVQE: Normalization, ConstrainedEnergyFunction

    import ..System
    import ..Float

    """ Uniquely define a trajectory. """
    @with_kw mutable struct SetupVars
        # SYSTEM PARAMETERS
        code::String = "H215"
        T::Float = 10.0                     # DURATION OF PULSE
        W::Int = 3                          # NUMBER OF WINDOWS PER PULSE

        # SIMULATION PARAMETERS
        r::Int = round(Int,20T)             # NUMBER OF TROTTER STEPS
        m::Int = 2                          # NUMBER OF LEVELS PER TRANSMON

        # HARDWARE BOUNDS
        ΩMAX::Float = 2π * 0.02 # GHz       # LARGEST PULSE AMPLITUDE
        ΔMAX::Float = 2π * 1.00 # GHz       # LARGEST PULSE DETUNING

        # PARAMETER INITIALIZATION
        seed_Ω::Int = 0000                  # DEFINE A PSEUDO-RANDOM NUMBER SEQUENCE
        seed_φ::Int = 1111
        seed_Δ::Int = 2222
        kick_Ω::Float = 0.0                 # STRENGTH OF INITIAL PERTURBATION (±)
        kick_φ::Float = 0.0
        kick_Δ::Float = 0.0

        # PENALTY VARIABLES
        λΩ::Float = 1.0 # Ha                # STRENGTH OF PENALTY TERMS
        λΔ::Float = 1.0 # Ha
        σΩ::Float = ΩMAX                    # SMOOTHNESS OF PENALTY TERMS
        σΔ::Float = ΔMAX
    end

    """ How scripts interact with a (fixed) trajectory. Mostly, convergence criteria. """
    @with_kw mutable struct MetaVars
        f_tol::Float = eps(Float)
        g_tol::Float = 1e-6
        maxiter::Int = 10000
        fnRATIO::Float = 10.0               # NEVER EXCEED RATIO FN CALLS / ITERATIONS
        update::Int = 500                   # CALL `report` EVERY SO MANY ITERATIONS
    end

    """ Store plottable data over the course of a trajectory. """
    @with_kw mutable struct TraceVars
        iterations::Vector{Int}
        f_calls::Vector{Int}
        g_calls::Vector{Int}

        fn::Vector{Float}
        gd::Vector{Float}
        energy::Vector{Float}
        penalties::Vector{Vector{Float}}

        Ωmax::Vector{Vector{Float}}
        ν::Vector{Vector{Float}}
    end

    """ All the dynamic data needed to resume a trajectory. """
    @with_kw mutable struct StateVars
        x::Vector{Float}
        Ω::Vector{Int}              # INDICES OF x REPRESENTING AMPLITUDES
        φ::Vector{Int}              # INDICES OF x REPRESENTING PHASES
        ν::Vector{Int}              # INDICES OF x REPRESENTING FREQUENCIES
        Hk::Matrix{Float}           # INVERSE HESSIAN APPROXIMATION
    end

    """ All the non-scalar objects needed to calculate a trajectory. """
    @with_kw mutable struct WorkVars
        system::System
        grid::IntegrationType
        protopulse::SignalType
        device::DeviceType
        normfn::Normalization
        lossfn::ConstrainedEnergyFunction
    end

    """ Bundle all of the above. """
    @with_kw mutable struct Vars
        outdir::String = "jobs/job"

        setup::SetupVars = SetupVars()
        meta::MetaVars = MetaVars()
        trace::Union{Nothing,TraceVars} = nothing
        state::Union{Nothing,StateVars} = nothing
        work::Union{Nothing,WorkVars} = nothing

        run::Bool = true    # FLAG THAT SCRIPT SHOULD RUN LENGTHY CALCULATIONS
    end


    const require_work(vars) = isnothing(vars.work) && error("Set work variables first!")


    """ Prepare a standard state suitable for starting a trajectory defined by vars. """
    function initial_state(vars)
        require_work(vars)
        # ALIASES FOR CONVENIENCE
        setup = vars.setup
        work = vars.work

        # PREPARE INDEXING VECTORS
        n = CtrlVQE.ndrives(work.device)
        Ω, φ, ν = indices(vars, work.protopulse, work.device)
        L = sum(length, (Ω, φ, ν))

        # PREPARE INITIAL PARAMETER VECTOR
        x = zeros(Float, L)
        for i in eachindex(ν)
            q = CtrlVQE.drivequbit(work.device, i)
            x[ν[i]] = CtrlVQE.resonancefrequency(work.device, q)
        end
        Random.seed!(setup.seed_Ω); x[Ω] .+= setup.kick_Ω .* (2 .* rand(length(Ω)) .- 1)
        Random.seed!(setup.seed_φ); x[φ] .+= setup.kick_φ .* (2 .* rand(length(φ)) .- 1)
        Random.seed!(setup.seed_Δ); x[ν] .+= setup.kick_Δ .* (2 .* rand(length(ν)) .- 1)

        # PREPARE INITIAL INVERSE HESSIAN APPROXIMATION
        Hk = LinearAlgebra.diagm(ones(Float, L))

        return StateVars(x, Ω, φ, ν, Hk)
    end

    """ Prepare a standard trace suitable for describing a trajectory defined by vars. """
    function initial_trace(vars)
        return TraceVars(
            iterations = Int[],
            f_calls = Int[],
            g_calls = Int[],

            fn = Float[],
            gd = Float[],
            energy = Float[],
            penalties = Vector{Float}[],

            Ωmax = Vector{Float}[],
            ν = Vector{Float}[],
        )
    end

    """ Update trace variables and check whether trajectory should terminate.

    For this job,
        we need the iteration count, current function value, and current gradient norm.
    The rest can be found by inspecting work variables.

    """
    function update_trace!(vars, iteration, fn, gd)
        require_work(vars)
        trace = vars.trace
        work  = vars.work

        push!(trace.iterations, iteration)
        push!(trace.f_calls, work.lossfn.f_counter[])
        push!(trace.g_calls, work.lossfn.g_counter[])

        push!(trace.fn, fn)
        push!(trace.gd, gd)

        push!(trace.energy, work.lossfn.energy[])
        push!(trace.penalties, copy(work.lossfn.penalties))

        nD = CtrlVQE.ndrives(work.device)
        Ωmax = Vector{Float}(undef, nD)
        ν    = Vector{Float}(undef, nD)
        for i in 1:nD
            # CALCULATE LARGEST Ω OVER ALL WINDOWS IN PULSE
            pulse = CtrlVQE.drivesignal(work.device, i)
            # TODO: Who signed off on this?! Needs to be type-agnostic.
            # Ωmax[i] = maximum(s -> abs(pulse(s)), pulse.starttimes)
            # TODO: This is a hack solution but we should use work variables more wisely.
            Ωmax[i] = maximum(abs.(CtrlVQE.valueat(pulse, CtrlVQE.lattice(work.grid))))

            # FETCH DRIVE FREQUENCY OF PULSE
            ν[i] = CtrlVQE.drivefrequency(vars.work.device, i)
        end
        push!(trace.Ωmax, Ωmax)
        push!(trace.ν, ν)
    end

    """ Inspect last record of the trace to see if it meets the convergence conditions. """
    function is_converged(vars)
        isnothing(vars.trace) && return false

        # SINGLE-ITERATION CHECKS
        isempty(vars.trace.iterations) && return false
        last(vars.trace.gd) ≤ vars.meta.g_tol && return true

        # TWO-ITERATION CHECKS
        length(vars.trace.iterations) < 2 && return false
        abs(vars.trace.fn[end-1] - vars.trace.fn[end]) ≤ vars.meta.f_tol && return true

        return false
    end

    """ Inspect last record of the trace to see if it meets the termination conditions. """
    function is_terminated(vars)
        isnothing(vars.trace) && return false
        isempty(vars.trace.iterations) && return false

        iterations = last(vars.trace.iterations)
        f_calls = last(vars.trace.f_calls)

        return any((
            is_converged(vars),
            # TOO MANY ITERATIONS - You can just run again from the final state.
            iterations >= vars.meta.maxiter,
            # TOO MANY FUNCTION EVALUATIONS - The linesearch seems to be unstable.
            iterations >= 10 && (f_calls > vars.meta.fnRATIO * iterations),
        ))
    end

    """ Update work objects to match the current state.

    TODO Two things:
    1. This method is script-dependent and therefore does not belong here.
    2. This method is entirely unnecessary for this job which has a fixed ansatz.
    I am leaving it here temporarily, so that it's easy to grab when implementing ADAPT job.

    """
    function update_work!(vars)
        require_work(vars)
        pulsetype = typeof(vars.work.protopulse)

        cursor = 1
        for i in 1:CtrlVQE.ndrives(vars.work.device)
            windows = pulsetype[]
            for (k, s) in enumerate(vars.state.s[i])
                push!(windows, pulsetype(
                    _!.state.x[cursor],
                    _!.state.x[cursor+1],
                ))
                cursor += 2
            end
            signal = CtrlVQE.WindowedSignal(windows, _!.state.s[i])
            CtrlVQE.set_drivesignal(device, i, signal)
        end
        CtrlVQE.Parameters.bind(device, _!.state.x)
    end



    ######################################################################################
    #= PARAMETER INDEXING - dispatched on pulse type and device type =#

    function indices(vars,
        ::CtrlVQE.SignalType,
        ::CtrlVQE.DeviceType,
    )
        require_work(vars)
        n = CtrlVQE.ndrives(vars.work.device)
        p = CtrlVQE.Parameters.count(vars.work.protopulse)
        L = n * (p * vars.setup.W + 1)
        return (Ω=collect(1:L-n), φ=Int[], ν=collect(1+L-n:L))
    end

    function indices(vars,
        ::CtrlVQE.SignalType,
        ::CtrlVQE.FixedFrequencyTransmonDevice,
    )
        require_work(vars)
        n = CtrlVQE.ndrives(vars.work.device)
        p = CtrlVQE.Parameters.count(vars.work.protopulse)
        L = n * p * vars.setup.W
        return (Ω=collect(1:L), φ=Int[], ν=Int[])
    end

    function indices(vars,
        ::CtrlVQE.PolarComplexConstant,
        ::CtrlVQE.DeviceType,
    )
        require_work(vars)
        n = CtrlVQE.ndrives(vars.work.device)
        L = n * (2 * vars.setup.W + 1)
        return (Ω=collect(1:2:L-n), φ=collect(2:2:L-n), ν=collect(1+L-n:L))
    end

    function indices(vars,
        ::CtrlVQE.PolarComplexConstant,
        ::CtrlVQE.FixedFrequencyTransmonDevice,
    )
        require_work(vars)
        n = CtrlVQE.ndrives(vars.work.device)
        L = n * 2 * vars.setup.W
        return (Ω=collect(1:2:L), φ=collect(2:2:L), ν=Int[])
    end
end