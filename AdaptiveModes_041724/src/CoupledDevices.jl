#= The ultimate strategy for superposition is
    a linesearch along the gradient signal at each adaptation.

This means associating one parameter with pulses over many drives.

So, we need a DeviceType that can handle that.

=#

module CoupledDevices
    import CtrlVQE: Parameters, Devices
    export CoupledDevice

    import CtrlVQE

    """

    A generic, extensible version of this idea is under way in CtrlVQE.jl,
        but it is not likely to be finished anytime soon.

    So we are making one here, with simplifying assumptions.
    - Each drivesignal is a WeightedCompositeSignal of `LΩ` parameter-less signals.
    - The first LΩ * nD parameters in base are just the weights for each drivesignal.

    """
    struct CoupledDevice{F,FΩ} <: Devices.LocallyDrivenDevice{F,FΩ}
        base::Devices.DeviceType{F,FΩ}
        x::Vector{F}
    end

    #= Some convenient functions to manage the particulars of this use case.
        They don't necessarily belong here. =#

    """ Convert a vector with a time grid into a WindowedSignal. """
    function as_trotterized_signal(α, β, grid)
        t̄ = CtrlVQE.lattice(grid)
        starttimes = t̄ .- [1; diff(t̄) ./ 2] # Slight offset so Ω(t) is never vague.

        component(i) = CtrlVQE.ConstrainedSignal(
            CtrlVQE.ComplexConstant(α[i], β[i]),
            :A, :B
        )

        return CtrlVQE.WindowedSignal(
            [component(i) for i in eachindex(starttimes)],
            starttimes,
        )
    end

    """ Completely constrain a pulse to the given values. """
    function constrained(pulse::CtrlVQE.ParametricSignal)
        return CtrlVQE.ConstrainedSignal(pulse, CtrlVQE.parameters(pulse)...)
    end

    constrained(pulse::CtrlVQE.ConstrainedSignal) = constrained(pulse.constrained)

    #= TODO: Generalize for the various composite pulses?
        I mean, there's no reason to...but it would be "robust". =#

    """ Make a superposition of constrained single-parameter pulses. """
    function as_composite_signal(basis, weights)
        pulses = CtrlVQE.CompositeSignal((deepcopy(pulse) for pulse in basis)...)
        CtrlVQE.Parameters.bind(pulses, weights)

        for i in eachindex(pulses.components)
            pulses.components[i] = constrained(pulses.components[i])
        end
        return pulses
    end

    #= TODO: Some convenience function for taking a basis set and a vector,
            binding each element of the vector to the element of the basis set,
            and then "freezing" each basis set,
            and combining the frozen ones into a CompositeSignal.
        Clearly depends on the type of the basis set.
        But if its type is a ParametricSignal or a ConstrainedSignal,
            it can be done automatically.
        Still, such a method clearly doesn't belong here. :)
    =#

    #= Special functions belonging to the `CoupledDevice` type. =#

    """ Constructor initializing a CoupledDevice with no parameters.
    Note that the drive signals in `base` will be overwritten,
        but they must still be chosen with care to fix the type parameters.
    """
    function CoupledDevice(base::Devices.DeviceType{F,FΩ}) where {F,FΩ}
        base = deepcopy(base)

        # Initialize drivesignals of base to empty WeightedCompositeSignals.
        initpulse = CtrlVQE.WeightedCompositeSignal(CtrlVQE.SignalType{F,FΩ}[], F[])
        for i in 1:Devices.ndrives(base)
            Devices.set_drivesignal(base, i, deepcopy(initpulse))
        end

        return CoupledDevice(base, F[])
    end

    """ Fetch the number of parameters belonging to the drive parameters
        (which is presumably the total number of parameters,
            less the number controlling drive frequencies). """
    function ndriveparams(device::CoupledDevice)
        signals = Devices.__get__drivesignals(device.base)
        return isempty(signals) ? 0 : length(first(signals).weights)
        #= TODO: I'd prefer `Parameter.count(first(signals))`
            but `WeightedCompositeSignal` is not very robust evidently. =#
    end

    """ Remove all drive parameters. """
    function empty_signals!(device::CoupledDevice)
        drivesignals = Devices.__get__drivesignals(device)
        for i in eachindex(drivesignals)
            empty!(drivesignals[i].components)
            empty!(drivesignals[i].weights)
        end

        empty!(device.x)
        return device.x
    end

    """ Add a new drive parameter. """
    function add_signals!(device::CoupledDevice, newsignals)
        drivesignals = Devices.__get__drivesignals(device)
        @assert length(newsignals) == length(drivesignals)

        for i in eachindex(drivesignals)
            push!(drivesignals[i].components, newsignals[i])
            push!(drivesignals[i].weights, zero(eltype(device.x)))
        end

        push!(device.x, zero(eltype(device.x)))

        return device
    end

    #= Methods to juggle parameters. =#

    function Parameters.count(device::CoupledDevice)
        return length(device.x)
    end

    function Parameters.names(device::CoupledDevice)
        basenames = Parameters.names(device.base)

        LΩ = ndriveparams(device)
        return [("c$l" for l in 1:LΩ)...; basenames[1+LΩ:end]]
    end

    function Parameters.values(device::CoupledDevice)
        return device.x
    end

    function Parameters.bind(device::CoupledDevice, x::AbstractVector)
        device.x .= x

        isempty(device.x) && return
        #= TODO: Shouldn't be necessary,
            and also I think even this breaks if you have frequency parameters.
            `WeightedCompositeSignal` is not very robust, evidently.

            (It needs to handle cases where `components` is empty.)
            =#

        LΩ = ndriveparams(device)
        nD = Devices.ndrives(device)
        expanded = [repeat(x[1:LΩ], nD); x[1+LΩ:end]]
        Parameters.bind(device.base, expanded)
    end

    function Devices.gradient(device::CoupledDevice,
        grid::CtrlVQE.Integrations.IntegrationType,
        ϕ̄::AbstractMatrix;
        result=nothing,
    )
        isempty(device.x) && return similar(device.x)
        #= TODO: Shouldn't be necessary,
            and also I think even this breaks if you have frequency parameters.
            `WeightedCompositeSignal` is not very robust, evidently.

            (It needs to handle cases where `components` is empty.)
            =#

        basegradient = Devices.gradient(device.base, grid, ϕ̄; result=result)

        LΩ = ndriveparams(device)
        nD = Devices.ndrives(device)
        g = similar(device.x)
        for l in 1:LΩ
            g[l] = sum(basegradient[l:LΩ:LΩ*nD])
        end
        g[1+LΩ:end] .= basegradient[1+LΩ*nD:end]

        return g
    end


    #= Delegate all remaining functions to `base`. =#

    Devices.nqubits(device::CoupledDevice) = Devices.nqubits(device.base)
    Devices.nlevels(device::CoupledDevice) = Devices.nlevels(device.base)
    Devices.ndrives(device::CoupledDevice) = Devices.ndrives(device.base)
    Devices.ngrades(device::CoupledDevice) = Devices.ngrades(device.base)

    Devices.localloweringoperator(device::CoupledDevice; result=nothing) =
        Devices.localloweringoperator(device.base; result=result)

    MatrixList = CtrlVQE.LinearAlgebraTools.MatrixList
    Devices.qubithamiltonian(device::CoupledDevice,
            ā::MatrixList, q::Int; result=nothing) =
        Devices.qubithamiltonian(device.base, ā, q; result=result)
    Devices.staticcoupling(device::CoupledDevice,
            ā::MatrixList; result=nothing) =
        Devices.staticcoupling(device.base, ā; result=result)
    Devices.driveoperator(device::CoupledDevice,
            ā::MatrixList, i::Int, t::Real; result=nothing) =
        Devices.driveoperator(device.base, ā, i, t; result=result)
    Devices.gradeoperator(device::CoupledDevice,
            ā::MatrixList, j::Int, t::Real; result=nothing) =
        Devices.gradeoperator(device.base, ā, j, t; result=result)

    Devices.eltype_localloweringoperator(device::CoupledDevice) =
        Devices.eltype_localloweringoperator(device.base)
    Devices.eltype_qubithamiltonian(device::CoupledDevice) =
        Devices.eltype_qubithamiltonian(device.base)
    Devices.eltype_staticcoupling(device::CoupledDevice) =
        Devices.eltype_staticcoupling(device.base)
    Devices.eltype_driveoperator(device::CoupledDevice) =
        Devices.eltype_driveoperator(device.base)
    Devices.eltype_gradeoperator(device::CoupledDevice) =
        Devices.eltype_gradeoperator(device.base)

    Devices.resonancefrequency(device::CoupledDevice, q::Int) =
        Devices.resonancefrequency(device.base, q)
    Devices.drivefrequency(device::CoupledDevice, i::Int) =
        Devices.drivefrequency(device.base, i)

    Devices.__get__drivesignals(device::CoupledDevice) =
        Devices.__get__drivesignals(device.base)

    Devices.drivequbit(device::CoupledDevice, i::Int) =
        Devices.drivequbit(device.base, i)
    Devices.gradequbit(device::CoupledDevice, j::Int) =
        Devices.gradequbit(device.base, j)









    import CtrlVQE: CostFunctions
    export CoupledGlobalAmplitudeBound, CoupledGlobalFrequencyBound

    import CtrlVQE: Parameters, Integrations, Devices, Signals

    # NOTE: Implicitly use smooth bounding function.
    wall(u) = exp(u - 1/u)
    grad(u) = exp(u - 1/u) * (1 + 1/u^2)

    """ Identical to GlobalAmplitudeBound except grad_function assumes CoupledDevice. """
    struct CoupledGlobalAmplitudeBound{F,FΩ} <: CostFunctions.CostFunctionType{F}
        device::CoupledDevice{F,FΩ}
        grid::Integrations.IntegrationType{F}
        ΩMAX::F             # MAXIMUM PERMISSIBLE AMPLITUDE
        λ::F                # STRENGTH OF BOUND
        σ::F                # STEEPNESS OF BOUND

        function CoupledGlobalAmplitudeBound(
            device::CoupledDevice{DF,FΩ},
            grid::Integrations.IntegrationType{IF},
            ΩMAX::Real,
            λ::Real,
            σ::Real,
        ) where {DF,FΩ,IF}
            F = promote_type(Float16, DF, real(FΩ), IF, eltype(ΩMAX), eltype(λ), eltype(σ))
            return new{F,FΩ}(device, grid, ΩMAX, λ, σ)
        end
    end

    Base.length(fn::CoupledGlobalAmplitudeBound) = Parameters.count(fn.device)

    function CostFunctions.cost_function(fn::CoupledGlobalAmplitudeBound{F,FΩ}) where {F,FΩ}
        t̄ = Integrations.lattice(fn.grid)                       # CACHED, THEREFORE FREE
        Ω̄ = Vector{FΩ}(undef, length(t̄))                        # TO FILL, FOR EACH DRIVE

        Φ(t, Ω) = (
            u = (abs(Ω) - fn.ΩMAX) / fn.σ;
            u ≤ 0 ? zero(u) : fn.λ * wall(u)
        )

        return (x̄) -> (
            Parameters.bind(fn.device, x̄);
            total = zero(F);
            for i in 1:Devices.ndrives(fn.device);
                signal = Devices.drivesignal(fn.device, i);
                Ω̄ = Signals.valueat(signal, t̄; result=Ω̄);
                total += Integrations.integrate(fn.grid, Φ, Ω̄)
            end;
            total
        )
    end

    function CostFunctions.grad_function_inplace(
        fn::CoupledGlobalAmplitudeBound{F,FΩ},
    ) where {F,FΩ}
        t̄ = Integrations.lattice(fn.grid)                       # CACHED, THEREFORE FREE
        Ω̄ = Vector{FΩ}(undef, length(t̄))                        # TO FILL, FOR EACH DRIVE
        ∂̄ = Vector{FΩ}(undef, length(t̄))                        # TO FILL, FOR EACH PARAMETER

        Φ(t, Ω, ∂) = (
            u = (abs(Ω) - fn.ΩMAX) / fn.σ;
            u ≤ 0 ? zero(u) : fn.λ * grad(u) * real(conj(Ω)*∂) / (abs(Ω)*fn.σ)
        )

        return (∇f̄, x̄) -> (
            Parameters.bind(fn.device, x̄);
            ∇f̄ .= 0;
            for i in 1:Devices.ndrives(fn.device);
                signal = Devices.drivesignal(fn.device, i);
                Ω̄ = Signals.valueat(signal, t̄; result=Ω̄);
                for l in 1:ndriveparams(fn.device);
                    ∂̄ = Signals.partial(l, signal, t̄; result=∂̄);
                    ∇f̄[l] += Integrations.integrate(fn.grid, Φ, Ω̄, ∂̄);
                end;
            end;
            ∇f̄
        )
    end



    """ Identical to GlobalFrequencyBound except grad_function assumes CoupledDevice. """
    struct CoupledGlobalFrequencyBound{F,FΩ} <: CostFunctions.CostFunctionType{F}
        device::CoupledDevice{F,FΩ}
        grid::Integrations.IntegrationType{F}
        ΔMAX::F             # MAXIMUM PERMISSIBLE DETUNING
        λ::F                # STRENGTH OF BOUND
        σ::F                # STEEPNESS OF BOUND

        function CoupledGlobalFrequencyBound(
            device::CoupledDevice{DF,FΩ},
            grid::Integrations.IntegrationType{IF},
            ΔMAX::Real,
            λ::Real,
            σ::Real,
        ) where {DF,FΩ,IF}
            F = promote_type(Float16, DF, real(FΩ), IF, eltype(ΔMAX), eltype(λ), eltype(σ))
            return new{F,FΩ}(device, grid, ΔMAX, λ, σ)
        end
    end

    Base.length(fn::CoupledGlobalFrequencyBound) = Parameters.count(fn.device)

    function CostFunctions.cost_function(fn::CoupledGlobalFrequencyBound{F,FΩ}) where {F,FΩ}
        return (x̄) -> (
            Parameters.bind(fn.device, x̄);
            total = zero(F);
            for i in 1:Devices.ndrives(fn.device);
                q = Devices.drivequbit(fn.device, i);
                Δ = Devices.detuningfrequency(fn.device, i, q);
                u = (abs(Δ) - fn.ΔMAX) / fn.σ;
                total += u ≤ 0 ? zero(u) : fn.λ * wall(u);
            end;
            total
        )
    end

    function CostFunctions.grad_function_inplace(
        fn::CoupledGlobalFrequencyBound{F,FΩ},
    ) where {F,FΩ}
        # INFER WHICH PARAMETERS REFER TO DRIVE FREQUENCIES
        nD = Devices.ndrives(fn.device)
        LΩ = ndriveparams(fn.device)

        LΩ == length(fn) && return (∇f̄, x̄) -> (∇f̄ .= 0; ∇f̄)     # NO FREQUENCIES
        LΩ + nD == length(fn) || error("Ill-defined number of frequency parameters.")

        # AT THIS POINT (for now) ASSUME x[offset+i] == ith frequency
        return (∇f̄, x̄) -> (
            Parameters.bind(fn.device, x̄);
            ∇f̄ .= 0;
            for i in 1:nD;
                q = Devices.drivequbit(fn.device, i);
                Δ = Devices.detuningfrequency(fn.device, i, q);
                u = (abs(Δ) - fn.ΔMAX) / fn.σ;
                ∇f̄[LΩ+i] += u ≤ 0 ? zero(u) : fn.λ * grad(u) * sign(Δ) / fn.σ;
            end;
            ∇f̄
        )
    end

end