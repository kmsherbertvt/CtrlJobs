module Variables
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

        # SIMULATION PARAMETERS
        r::Int = round(Int,20T)             # NUMBER OF TROTTER STEPS
        m::Int = 2                          # NUMBER OF LEVELS PER TRANSMON
        seed::Int = 0                       # RANDOM SEED TO BE USED BY SCRIPTS

        # HARDWARE BOUNDS
        ΩMAX::Float = 2π * 0.02 # GHz       # LARGEST PULSE AMPLITUDE
        ΔMAX::Float = 2π * 1.00 # GHz       # LARGEST PULSE DETUNING
        fMAX::Float =      0.17 # GHz       # LARGEST PULSE FREQUENCY (ish)

        # PENALTY VARIABLES
        λΩ::Float = 1.0 # Ha                # STRENGTH OF PENALTY TERMS
        λΔ::Float = 1.0 # Ha
        σΩ::Float = ΩMAX                    # SMOOTHNESS OF PENALTY TERMS
        σΔ::Float = ΔMAX

        # OPTIMIZATION CONVERGENCE
        f_tol::Float = eps(Float)
        g_tol::Float = 1e-6
        maxiter::Int = 10000
        fnRATIO::Float = 10.0               # NEVER EXCEED RATIO FN CALLS / ITERATIONS
    end

    """ How scripts interact with a (fixed) trajectory. Mostly, convergence criteria. """
    @with_kw mutable struct MetaVars
        G_tol::Float = 1e-3                 # ADAPT CONVERGENCE CRITERION
        maxadapt::Int = 200                 # NUMBER OF ADAPT ITERATIONS
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

        adaptations::Vector{Int}            # INDICES WHERE VQE IS OPTIMIZED
        poolsize::Vector{Int}               # NUMBER OF OPERATORS IN THE CURRENT POOL
        G_max::Vector{Float}                # LARGEST GRADIENT IN THE POOL
        parameters::Vector{Vector{Int}}     # COUNTS PARAMETERS FOR EACH DRIVE
    end

    """ All the dynamic data needed to resume a trajectory. """
    @with_kw mutable struct StateVars
        x::Vector{Float}
        Hk::Matrix{Float}           # INVERSE HESSIAN APPROXIMATION
        n::Vector{Vector{Int}}      # LIST OF MODES FOR EACH DRIVE
        Ω::Vector{Vector{Int}}      # INDICES OF x FOR EACH DRIVE
        ν::Vector{Int}              # INDICES OF x REPRESENTING FREQUENCIES
    end

    """ All the non-scalar objects needed to calculate a trajectory. """
    @with_kw mutable struct WorkVars
        system::System
        grid::IntegrationType
        protopulse::SignalType
        device::DeviceType
        normfn::Normalization
        lossfn::ConstrainedEnergyFunction
        pool::Vector{<:SignalType}
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
        # ALIAS FOR CONVENIENCE
        work = vars.work

        # PREPARE SIGNAL INDEXING VECTORS
        nD = CtrlVQE.ndrives(work.device)
        Ω = [Int[] for i in 1:nD]
        L = sum(length, Ω)

        # PREPARE FREQUENCY INDEXING VECTORS -- Implementation may vary with type someday?
        ν = collect(L+1:L+nD)
        work.device isa CtrlVQE.FixedFrequencyTransmonDevice && empty!(ν)
        L += length(ν)

        # PREPARE INITIAL PARAMETER VECTOR
        x = zeros(Float, L)
        for i in eachindex(ν)
            q = CtrlVQE.drivequbit(work.device, i)
            x[ν[i]] = CtrlVQE.resonancefrequency(work.device, q)
        end

        # PREPARE INITIAL INVERSE HESSIAN APPROXIMATION
        Hk = LinearAlgebra.diagm(ones(Float, L))

        # PREPARE INITIAL STARTTIMES
        n = [Int[] for _ in 1:nD]

        return StateVars(x, Hk, n, Ω, ν)
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

            adaptations = Int[],
            poolsize = Int[],
            G_max = Float[],
            parameters = Vector{Int}[],
        )
    end

end