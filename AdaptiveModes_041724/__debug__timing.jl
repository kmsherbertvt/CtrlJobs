#= `gradients.harmonics` is taking WAY too long.

Acknowledged, `adapt_work!` is sillily inefficient;
    we should just be saving a matrix of coefficients in the trace,
    rather than re-calculating the gradient signals at every adaptation,
    at every adaptation.
I can fake this by manually adding in a `coefficients.npy` file via the script.

But even individual function evaluations are taking longer than the really should,
    for a 4-qubit system.
I think the composite signals are just too much for Julia.
Too many nested function calls, too much stack trace.
So let us instead customize a new SignalType that has everything all in one.

I'm afraid the real killer is the penalty function, which this alone won't address.
But let's see how much this helps.

=#

module ModalHarmonics
    import CtrlVQE: Parameters, Signals, SignalType
    export ModalHarmonic

    import LinearAlgebra: dot

    struct ModalHarmonic{F} <: SignalType{F,Complex{F}}
        A::Vector{Vector{F}}    # AMPLITUDE - REAL PART     (A[n][a])
        B::Vector{Vector{F}}    # AMPLITUDE - IMAG PART     (B[n][a])
        C::Vector{F}            # VARIATIONAL PARAMETERS    (C[a])
        T::F                    # FUNDAMENTAL PERIOD (aka T such that ν=n/2T)
    end

    #= Utilities =#

    nmodes(signal::ModalHarmonic) = length(signal.A)
    eachmode(signal::ModalHarmonic) = 1:nmodes(signal)

    get_Ωn(signal::ModalHarmonic, n::Int) = Complex(
        dot(signal.A[n], signal.C),
        dot(signal.B[n], signal.C),
    )

    get_νn(signal::ModalHarmonic, n::Int) = n * (π / signal.T)


    #= Constructors =#

    """ Initialize an empty ModalHarmonic signal with nMAX modes. """
    function ModalHarmonic(nMAX::Int, T::F) where F
        return ModalHarmonic(
            [F[] for n in 1:nMAX],
            [F[] for n in 1:nMAX],
            F[], T,
        )
    end

    """ Initialize a ModalHarmonic signal from a coefficient matrix. """
    function ModalHarmonic(AB::Matrix{F}, T::F) where F
        return ModalHarmonic(
            Vector{F}[eachcol(real.(AB))...],
            Vector{F}[eachcol(imag.(AB))...],
            zeros(F, size(AB, 1)),
            T,
        )
    end

    #= Mutators =#

    """ Add on a new row, with coefficient initialized to zero. """
    function Base.append!(signal::ModalHarmonic{F}, A::Vector{F}, B::Vector{F}) where F
        @assert length(A) == length(B)
        for n in eachmode(signal)
            push!(signal.A[n], A[n])
            push!(signal.B[n], B[n])
        end
        push!(signal.C, 0)
        return signal
    end

    """ Add on many rows, with coefficients initialized to zero. """
    function Base.append!(signal::ModalHarmonic{F}, AB::Matrix{Complex{F}}) where F
        for n in eachmode(signal)
            append!(signal.A[n], real.(AB[:,n]))
            append!(signal.B[n], imag.(AB[:,n]))
        end
        append!(signal.C, zeros(F, size(AB,1)))
        return signal
    end

    """ Empty out all rows. """
    function empty!(signal::ModalHarmonic)
        for n in eachmode(signal)
            empty!(signal.A[n])
            empty!(signal.B[n])
        end
        empty!(signal.C)
        return signal
    end

    #= SignalType interface =#

    function Signals.valueat(
        signal::ModalHarmonic, t::Real;
        Ω=nothing, ν=nothing,
    )
        isnothing(Ω) && (Ω = get_Ωn.(Ref(signal), eachmode(signal)))
        isnothing(ν) && (ν = get_νn.(Ref(signal), eachmode(signal)))
        return sum(Ω .* sin.(ν.*t))
    end

    """ Same vectorization as default but only calculate Ω, ν once. """
    function Signals.valueat(
        signal::ModalHarmonic{F}, t̄::AbstractVector{<:Real};
        Ω=nothing, ν=nothing, result=nothing,
    ) where {F}
        isnothing(Ω) && (Ω = get_Ωn.(Ref(signal), eachmode(signal)))
        isnothing(ν) && (ν = get_νn.(Ref(signal), eachmode(signal)))
        isnothing(result) && return Signals.valueat(
            signal, t̄;
            Ω=Ω, ν=ν, result=Vector{Complex{F}}(undef, size(t̄)),
        )
        for (i, t) in enumerate(t̄)
            result[i] = Signals.valueat(signal, t; Ω=Ω, ν=ν)
        end
        return result
    end

    function Signals.partial(
        k::Int, signal::ModalHarmonic{F}, t::Real;
        Ω=nothing, ν=nothing,
    ) where F
        isnothing(Ω) && (Ω = Complex.(
            [signal.A[n][k] for n in eachmode(signal)],
            [signal.B[n][k] for n in eachmode(signal)],
        ))
        isnothing(ν) && (ν = get_νn.(Ref(signal), eachmode(signal)))

        part = zero(F)
        for n in eachmode(signal)
            part += Ω[n] * sin(ν[n]*t)
        end
        return part
    end

    """ Same vectorization as default but only calculate Ω, ν once. """
    function Signals.partial(
        k::Int, signal::ModalHarmonic{F}, t̄::AbstractVector{<:Real};
        Ω=nothing, ν=nothing, result=nothing,
    ) where {F}
        isnothing(Ω) && (Ω = Complex.(
            [signal.A[n][k] for n in eachmode(signal)],
            [signal.B[n][k] for n in eachmode(signal)],
        ))
        isnothing(ν) && (ν = get_νn.(Ref(signal), eachmode(signal)))
        isnothing(result) && return Signals.partial(
            n, signal, t̄;
            Ω=Ω, ν=ν, result=Vector{Complex{F}}(undef, size(t̄)),
        )
        for (i, t) in enumerate(t̄)
            result[i] = Signals.partial(k, signal, t; Ω=Ω, ν=ν)
        end
        return result
    end

    function Base.string(signal::ModalHarmonic, names::AbstractVector{String})
        A = signal.A
        B = signal.B
        C = names
        T = signal.T

        return join(("$(C[n]) ($A+i$B) sin($(n)π⋅t/$T)" for n in eachindex(C)), " + ")
    end

    #= Parameters interface. =#

    function Parameters.count(signal::ModalHarmonic)
        return length(signal.C)
    end

    Parameters.count(signal::ModalHarmonic) = length(signal.C)
    Parameters.names(signal::ModalHarmonic) = ["c$n" for n in eachmode(signal)]
    Parameters.values(signal::ModalHarmonic) = signal.C
    function Parameters.bind(signal::ModalHarmonic{F}, x̄::Vector{F}) where {F}
        signal.C .= x̄
        return signal
    end









    import CtrlVQE: Parameters, Devices, TransmonDevices
    export ModalHarmonicTransmonDevice

    import CtrlVQE
    import CtrlVQE: Integrations
    import CtrlVQE: FixedFrequencyTransmonDevice
    import CtrlVQE.TransmonDevices: AbstractTransmonDevice
    import CtrlVQE.LinearAlgebraTools: MatrixList

    import CtrlVQE.TempArrays: array
    const LABEL = Symbol(@__MODULE__)

    """ Highly customized version of a FixedFrequencyTransmonDevice,
            where every drive signal is a ModalHarmonic.

    Not meant to be extensible in the slightest.

    Uses a "template" device, which is NOT copied.
    Its variational parameters are never used,
        though technically the frequencies and such could still be mutated.

    Throughout, we interact with vectors represented in a multi-index basis [ n ← αβ, i].
        That is, each column is associated with drive i.
        Within that, the first half are for real modes and the latter for imaginary.
        Within that, each index corresponds to the harmonic mode n.
    Collections of such vectors indexed by a are represented as [ a, n ← αβ, i].

    """
    struct ModalHarmonicTransmonDevice{F,FΩ} <: AbstractTransmonDevice{F,FΩ}
        template::FixedFrequencyTransmonDevice{F,FΩ}
        signals::Vector{ModalHarmonic{F}}
    end

    """ Constructor initializing a ModalHarmonicTransmonDevice with no parameters. """
    function ModalHarmonicTransmonDevice(
        template::FixedFrequencyTransmonDevice{F,FΩ},
        nMAX::Int,
        T::F,
    ) where {F,FΩ}
        signals = [ModalHarmonic(nMAX,T) for _ in 1:CtrlVQE.ndrives(template)]
        return ModalHarmonicTransmonDevice(template, signals)
    end

    """ Constructor initializing a ModalHarmonicTransmonDevice from a vector space. """
    function ModalHarmonicTransmonDevice(
        template::FixedFrequencyTransmonDevice{F,FΩ},
        space::Array{F,3},  # [ a, n ← σ, i ]
        T::F,
    ) where {F,FΩ}
        nMAX = size(space,2) >> 1   # Divide by 2 since each mode has α and β component.
        device = ModalHarmonicTransmonDevice(template, nMAX, T)
        append!(device, space)
        return device
    end

    """ Designate `device.signals` as the array of drive signals. """
    Devices.__get__drivesignals(device::ModalHarmonicTransmonDevice) = device.signals

    #= Special functions belonging to the `ModalHarmonicTransmonDevice` type. =#

    nmodes(device::ModalHarmonicTransmonDevice) = (
        isempty(device.signals) ? 0 : nmodes(first(device.signals))
    )
    eachmode(device::ModalHarmonicTransmonDevice) = 1:nmodes(device)
    eachdrive(device::ModalHarmonicTransmonDevice) = 1:CtrlVQE.ndrives(device)
        # TODO: This last, along with `eachqubit`, `eachstate`, etc. should be in Devices.

    """ Add on a new parameter along the given direction. """
    function Base.append!(
        device::ModalHarmonicTransmonDevice{F},
        vector::Array{F,2},  # [ n ← σ, i]
    ) where F
        @assert size(space,2) == 2 * nmodes(device)
        @assert size(space,3) == CtrlVQE.ndrives(device)
        V = reshape(vector, nmodes(device), 2, CtrlVQE.ndrives(device))     # V[n,σ,i]
        for i in eachdrive(device)
            append!(device.signals[i], @view(V[:,1,i]), @view(V[:,2,i]))
        end
        return device
    end

    """ Add on a whole vector space all at once, with coefficients initialized to zero. """
    function Base.append!(
        device::ModalHarmonicTransmonDevice{F},
        space::Array{F,3},  # [ a, n ← σ, i ]
    ) where F
        # for a in axes(space, 2)
        #     append!(device, @view(space[a,:,:]))
        # end
        # return device

        @assert size(space,2) == 2 * nmodes(device)
        @assert size(space,3) == CtrlVQE.ndrives(device)
        B = reshape(space, :, nmodes(device), 2, CtrlVQE.ndrives(device)) # B[a,n,σ,i]
        for i in eachdrive(device)
            AB = Complex.(B[:,:,1,i], B[:,:,2,i])
            append!(device.signals[i], AB)
        end
        return device
    end

    """ Empty out all variational parameters. """
    function empty!(device::ModalHarmonicTransmonDevice)
        for i in eachdrive(device)
            empty!(device.signals[i])
        end
        return device
    end


    #= Methods to juggle parameters. =#

    function Parameters.count(device::ModalHarmonicTransmonDevice)
        isempty(device.signals) && return 0
        return Parameters.count(first(device.signals))
    end

    function Parameters.names(device::ModalHarmonicTransmonDevice)
        isempty(device.signals) && return String[]
        return Parameters.names(first(device.signals))
    end

    function Parameters.values(device::ModalHarmonicTransmonDevice{F}) where {F}
        isempty(device.signals) && return F[]
        return Parameters.values(first(device.signals))
    end

    function Parameters.bind(device::ModalHarmonicTransmonDevice, x::AbstractVector)
        for i in eachdrive(device)
            Parameters.bind(device.signals[i], x)
        end
        return device
    end

    function Devices.gradient(device::ModalHarmonicTransmonDevice{F,FΩ},
        grid::Integrations.IntegrationType,
        ϕ̄::AbstractMatrix;
        result=nothing,
    ) where {F,FΩ}
        isnothing(result) && return Devices.gradient(
            device, grid, ϕ̄;
            result=Vector{F}(undef, Parameters.count(device)),
        )

        # CALCULATE GRADIENT FOR SIGNAL PARAMETERS
        t̄ = Integrations.lattice(grid)
        ∂̄ = array(FΩ, size(t̄), LABEL)
        Φ = (t, ∂, ϕα, ϕβ) -> (real(∂)*ϕα + imag(∂)*ϕβ)

        result .= 0
        for i in eachdrive(device)
            j = 2i - 1
            ϕ̄α = @view(ϕ̄[:,j])
            ϕ̄β = @view(ϕ̄[:,j+1])

            signal = Devices.drivesignal(device, i)
            for k in 1:CtrlVQE.Parameters.count(signal)
                ∂̄ = Signals.partial(k, signal, t̄; result=∂̄)
                result[k] += Integrations.integrate(grid, Φ, ∂̄, ϕ̄α, ϕ̄β)
            end
        end

        return result
    end

    #= Delegate all remaining transmon behavior to `template`. =#

    Devices.nlevels(device::ModalHarmonicTransmonDevice) =
        Devices.nlevels(device.template)
    Devices.nqubits(device::ModalHarmonicTransmonDevice) =
        Devices.nqubits(device.template)
    Devices.resonancefrequency(device::ModalHarmonicTransmonDevice, q::Int) =
        Devices.resonancefrequency(device.template, q)
    TransmonDevices.anharmonicity(device::ModalHarmonicTransmonDevice, q::Int) =
        TransmonDevices.anharmonicity(device.template, q)

    TransmonDevices.ncouplings(device::ModalHarmonicTransmonDevice) =
        TransmonDevices.ncouplings(device.template)
    TransmonDevices.couplingpair(device::ModalHarmonicTransmonDevice, k::Int) =
        TransmonDevices.couplingpair(device.template, k)
    TransmonDevices.couplingstrength(device::ModalHarmonicTransmonDevice, k::Int) =
        TransmonDevices.couplingstrength(device.template, k)

    Devices.ndrives(device::ModalHarmonicTransmonDevice) =
        Devices.ndrives(device.template)
    Devices.drivequbit(device::ModalHarmonicTransmonDevice, i::Int) =
        Devices.drivequbit(device.template, i)
    Devices.drivefrequency(device::ModalHarmonicTransmonDevice, i::Int) =
        Devices.drivefrequency(device.template, i)

    # NOTE: `__get__drivesignals` implemented above
    # NOTE: `bindfrequencies` not implemented






    import CtrlVQE: CostFunctions
    export ModalHarmonicAmplitudeBound

    import CtrlVQE: Parameters, Integrations, Devices, Signals

    # NOTE: Implicitly use smooth bounding function.
    wall(u) = exp(u - 1/u)
    grad(u) = exp(u - 1/u) * (1 + 1/u^2)

    """ Simplified GlobalAmplitudeBound restricted to `ModalHarmonicTransmonDevice`. """
    struct ModalHarmonicAmplitudeBound{F,FΩ} <: CostFunctions.CostFunctionType{F}
        device::ModalHarmonicTransmonDevice{F,FΩ}
        grid::Integrations.IntegrationType{F}
        ΩMAX::F             # MAXIMUM PERMISSIBLE AMPLITUDE
        λ::F                # STRENGTH OF BOUND
        σ::F                # STEEPNESS OF BOUND

        function ModalHarmonicAmplitudeBound(
            device::ModalHarmonicTransmonDevice{DF,FΩ},
            grid::Integrations.IntegrationType{IF},
            ΩMAX::Real,
            λ::Real,
            σ::Real,
        ) where {DF,FΩ,IF}
            F = promote_type(Float16, DF, real(FΩ), IF, eltype(ΩMAX), eltype(λ), eltype(σ))
            return new{F,FΩ}(device, grid, ΩMAX, λ, σ)
        end
    end

    Base.length(fn::ModalHarmonicAmplitudeBound) = Parameters.count(fn.device)

    function CostFunctions.cost_function(fn::ModalHarmonicAmplitudeBound{F,FΩ}) where {F,FΩ}
        t̄ = Integrations.lattice(fn.grid)                       # CACHED, THEREFORE FREE
        Ω̄ = Vector{FΩ}(undef, length(t̄))                        # TO FILL, FOR EACH DRIVE

        Φ(t, Ω) = (
            u = (abs(Ω) - fn.ΩMAX) / fn.σ;
            u ≤ 0 ? zero(u) : fn.λ * wall(u)
        )

        return (x̄) -> (
            Parameters.bind(fn.device, x̄);
            total = zero(F);
            for i in eachdrive(fn.device);
                signal = Devices.drivesignal(fn.device, i);
                Ω̄ = Signals.valueat(signal, t̄; result=Ω̄);
                total += Integrations.integrate(fn.grid, Φ, Ω̄)
            end;
            total
        )
    end

    function CostFunctions.grad_function_inplace(
        fn::ModalHarmonicAmplitudeBound{F,FΩ},
    ) where {F,FΩ}
        t̄ = Integrations.lattice(fn.grid)               # CACHED, THEREFORE FREE
        Ω̄ = Vector{FΩ}(undef, length(t̄))                # TO FILL, FOR EACH DRIVE
        ∂̄ = Vector{FΩ}(undef, length(t̄))                # TO FILL, FOR EACH PARAMETER

        Φ(t, Ω, ∂) = (
            u = (abs(Ω) - fn.ΩMAX) / fn.σ;
            u ≤ 0 ? zero(u) : fn.λ * grad(u) * real(conj(Ω)*∂) / (abs(Ω)*fn.σ)
        )

        return (∇f̄, x̄) -> (
            Parameters.bind(fn.device, x̄);
            ∇f̄ .= 0;
            for i in eachdrive(fn.device);
                signal = Devices.drivesignal(fn.device, i);
                Ω̄ = Signals.valueat(signal, t̄; result=Ω̄);
                for k in 1:CtrlVQE.Parameters.count(signal);
                    ∂̄ = Signals.partial(k, signal, t̄; result=∂̄);
                    ∇f̄[k] += Integrations.integrate(fn.grid, Φ, Ω̄, ∂̄);
                end;
            end;
            ∇f̄
        )
    end

end

##########################################################################################
#= SET UP WORK VARIABLES =#

import CtrlVQE

code = "lih30"
T = 20.0
nadapt = 20

setup=JOB.SetupVars(code=code, T=T)
vars = JOB.Vars(setup=setup)
system = JOB.System(code)

# PREP THE PULSES

r = 400
grid = CtrlVQE.TemporalLattice(T, r)

pool = JOB.makepool_harmonics(vars)

import Random; Random.seed!(0);
scores = rand(nadapt, length(pool), system.n)
    # scores[a,n←σ,i]

# INITIALIZE THE DEVICES


device_template = CtrlVQE.Systematic(
    CtrlVQE.FixedFrequencyTransmonDevice,
    system.n,
    CtrlVQE.ComplexConstant(zero(JOB.Float), zero(JOB.Float)),
)

device_modal = ModalHarmonics.ModalHarmonicTransmonDevice(device_template, scores, T)

device_coupled = JOB.CoupledDevice(device_template)
for a in 1:nadapt
    newsignals = [
        JOB.CoupledDevices.as_composite_signal(pool, scores[a,:,i])
            for i in 1:CtrlVQE.ndrives(device_coupled)
    ]
    JOB.CoupledDevices.add_signals!(device_coupled, newsignals)
end

# REMAINING WORK VARIABLES

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC

ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, device_template)
O0 = CtrlVQE.QubitOperators.project(system.model.H, device_template)

lossfn_modal = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device_modal, basis, frame, grid, ψ0, O0),
    ModalHarmonics.ModalHarmonicAmplitudeBound(device_modal,
        grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
)

lossfn_coupled = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device_coupled, basis, frame, grid, ψ0, O0),
    JOB.CoupledGlobalAmplitudeBound(device_coupled, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
)


##########################################################################################
#= TIME SOME THINGS =#

# A single evolution

ψ_modal = similar(ψ0)
CtrlVQE.evolve(evolution, device_modal, basis, grid, ψ0; result=ψ_modal)
println("Modal Evolution:")
@time CtrlVQE.evolve(evolution, device_modal, basis, grid, ψ0; result=ψ_modal)


ψ_coupled = similar(ψ0)
CtrlVQE.evolve(evolution, device_coupled, basis, grid, ψ0; result=ψ_coupled)
println("Composite Evolution:")
@time CtrlVQE.evolve(evolution, device_coupled, basis, grid, ψ0; result=ψ_coupled)

println("Fidelity: $(abs2(ψ_modal'*ψ_coupled))")
println()

# A lossfn evaluation

x0 = zeros(nadapt)

fn_modal = CtrlVQE.cost_function(lossfn_modal)
gd!_modal = CtrlVQE.grad_function_inplace(lossfn_modal)

g_modal = similar(x0)
println("Modal Cost/Gradient:")
fn_modal(x0); @time f_modal = fn_modal(x0)
gd!_modal(g_modal, x0); @time gd!_modal(g_modal, x0);


fn_coupled = CtrlVQE.cost_function(lossfn_coupled)
gd!_coupled = CtrlVQE.grad_function_inplace(lossfn_coupled)

g_coupled = similar(x0)
println("Composite Cost/Gradient:")
fn_coupled(x0); @time f_coupled = fn_coupled(x0)
gd!_coupled(g_coupled, x0); @time gd!_coupled(g_coupled, x0);

println("Energy error: $(f_modal - f_coupled)")
println("Gradient RMS Error: $(sqrt(sum((g_modal .- g_coupled).^2))/nadapt)")
println()


# Focus on energy evaluation

x0 = zeros(nadapt)

fn_modal = CtrlVQE.cost_function(lossfn_modal.energyfn)
gd!_modal = CtrlVQE.grad_function_inplace(lossfn_modal.energyfn)

g_modal = similar(x0)
println("Modal Energy Cost/Gradient:")
fn_modal(x0); @time f_modal = fn_modal(x0)
gd!_modal(g_modal, x0); @time gd!_modal(g_modal, x0);


fn_coupled = CtrlVQE.cost_function(lossfn_coupled.energyfn)
gd!_coupled = CtrlVQE.grad_function_inplace(lossfn_coupled.energyfn)

g_coupled = similar(x0)
println("Composite Energy Cost/Gradient:")
fn_coupled(x0); @time f_coupled = fn_coupled(x0)
gd!_coupled(g_coupled, x0); @time gd!_coupled(g_coupled, x0);

println("Energy error: $(f_modal - f_coupled)")
println("Gradient RMS Error: $(sqrt(sum((g_modal .- g_coupled).^2))/nadapt)")
println()

# Focus on grad evaluation

x0 = ones(nadapt) * vars.setup.ΩMAX / 2system.n / √2

fn_modal = CtrlVQE.cost_function(lossfn_modal.penaltyfns[1])
gd!_modal = CtrlVQE.grad_function_inplace(lossfn_modal.penaltyfns[1])

g_modal = similar(x0)
println("Modal Penalty Cost/Gradient:")
fn_modal(x0); @time f_modal = fn_modal(x0)
gd!_modal(g_modal, x0); @time gd!_modal(g_modal, x0);


fn_coupled = CtrlVQE.cost_function(lossfn_coupled.penaltyfns[1])
gd!_coupled = CtrlVQE.grad_function_inplace(lossfn_coupled.penaltyfns[1])

g_coupled = similar(x0)
println("Composite Penalty Cost/Gradient:")
fn_coupled(x0); @time f_coupled = fn_coupled(x0)
gd!_coupled(g_coupled, x0); @time gd!_coupled(g_coupled, x0);

println("Energy error: $(f_modal - f_coupled)")
println("Gradient RMS Error: $(sqrt(sum((g_modal .- g_coupled).^2))/nadapt)")
println()
