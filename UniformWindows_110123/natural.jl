#= This is a prototype for a device type which has its pulses "baked in".
    I suppose they can be edited post-hoc, eg. in an adaptive protocol,
    but all the parameters are just weights, and maybe (but not for now) detunings. =#

import CtrlVQE.TransmonDevices: AbstractTransmonDevice
import CtrlVQE: Parameters, Devices, TransmonDevices
import CtrlVQE: Quple, Integrations, Signals
import CtrlVQE: SignalType, CompositeSignal, ModulatedSignal, Constant, ComplexConstant

import CtrlVQE.TempArrays: array
const LABEL = Symbol(@__MODULE__)



struct NormalTransmonDevice{F,FΩ} <: AbstractTransmonDevice{F,FΩ}
    # QUBIT LISTS
    ω̄::Vector{F}
    δ̄::Vector{F}
    # COUPLING LISTS
    ḡ::Vector{F}
    quples::Vector{Quple}
    # DRIVE LISTS
    q̄::Vector{Int}
    ν̄::Vector{F}
    Ω̄::Vector{SignalType{F,FΩ}}
    eigenpulses::Vector{Vector{SignalType{F,FΩ}}}
    weights::Vector{F}
    # OTHER PARAMETERS
    m::Int

    function NormalTransmonDevice(
        ω̄::AbstractVector{<:Real},
        δ̄::AbstractVector{<:Real},
        ḡ::AbstractVector{<:Real},
        quples::AbstractVector{Quple},
        q̄::AbstractVector{Int},
        ν̄::AbstractVector{<:AbstractFloat},
        pulsesets::AbstractVector{<:AbstractVector{SignalType{F,FΩ}}},
        weights::AbstractVector{<:Real},
        m::Int,
    ) where {F,FΩ}
        # VALIDATE PARALLEL LISTS ARE CONSISTENT SIZE
        @assert length(ω̄) == length(δ̄) ≥ 1              # NUMBER OF QUBITS
        @assert length(ḡ) == length(quples)             # NUMBER OF COUPLINGS
        @assert length(q̄) == length(ν̄)                  # NUMBER OF DRIVES

        # VALIDATE QUBIT INDICES
        for (p,q) in quples
            @assert 1 <= p <= length(ω̄)
            @assert 1 <= q <= length(ω̄)
        end
        for q in q̄
            @assert 1 <= q <= length(ω̄)
        end

        # VALIDATE THAT THE HILBERT SPACE HAS SOME VOLUME...
        @assert m ≥ 2

        # CONSTRUCT THE ACTUAL SIGNALS
        @assert all(length(q̄) == length(pulseset) for pulseset in pulsesets)
        @assert length(weights) == length(pulsesets)
        Ω̄ = SignalType{F,FΩ}[]
        for i in 1:length(q̄)            # Initialize each drive as a CompositeSignal.
            Ω = CompositeSignal(SignalType{F,FΩ}[])
            for l in eachindex(pulsesets)   # Fill each drive from each pulseset.
                push!(Ω.components, ModulatedSignal(
                    ComplexConstant(weights[l], 0.0),   # TODO: Generalize
                    # Constant(weights[l]),   # TODO: Generalize
                    pulsesets[l][i],
                ))
            end
            push!(Ω̄, Ω)
        end

        # STANDARDIZE TYPING
        return new{F,FΩ}(
            convert(Vector{F}, ω̄),
            convert(Vector{F}, δ̄),
            convert(Vector{F}, ḡ),
            quples,
            q̄,
            convert(Vector{F}, ν̄),
            Ω̄,
            pulsesets, # TODO: pry need to convert this somehow?
            convert(Vector{F}, weights),
            m,
        )
    end
end

Devices.nlevels(device::NormalTransmonDevice) = device.m

Devices.nqubits(device::NormalTransmonDevice) = length(device.ω̄)
Devices.resonancefrequency(device::NormalTransmonDevice, q::Int) = device.ω̄[q]
TransmonDevices.anharmonicity(device::NormalTransmonDevice, q::Int) = device.δ̄[q]

TransmonDevices.ncouplings(device::NormalTransmonDevice) = length(device.quples)
TransmonDevices.couplingpair(device::NormalTransmonDevice, k::Int) = device.quples[k]
TransmonDevices.couplingstrength(device::NormalTransmonDevice, k::Int) = device.ḡ[k]

Devices.ndrives(device::NormalTransmonDevice) = length(device.q̄)
Devices.drivequbit(device::NormalTransmonDevice, i::Int)=device.q̄[i]
Devices.drivefrequency(device::NormalTransmonDevice, i::Int) = device.ν̄[i]
Devices.__get__drivesignals(device::NormalTransmonDevice) = device.Ω̄

TransmonDevices.bindfrequencies(device::NormalTransmonDevice, ν̄::AbstractVector) = nothing


function Parameters.count(device::NormalTransmonDevice)
    return length(device.weights)
end

function Parameters.names(device::NormalTransmonDevice)
    return ["c$l" for l in eachindex(device.weights)]
end

function Parameters.values(device::NormalTransmonDevice{F,FΩ}) where {F,FΩ}
    return copy(device.weights)
end

function Parameters.bind(
    device::NormalTransmonDevice,
    x̄::AbstractVector{F},
) where {F}
    device.weights .= x̄
    for i in 1:Devices.ndrives(device)
        for l in 1:length(device.weights)
            device.Ω̄[i].components[l].components[1].A = device.weights[l]
        end
    end
end

function Devices.gradient(
    device::NormalTransmonDevice{F,FΩ},
    grid::Integrations.IntegrationType,
    ϕ̄::AbstractMatrix;
    result=nothing,
) where {F,FΩ}
    L = Parameters.count(device)::Int
    isnothing(result) && return Devices.gradient(
        device, grid, ϕ̄;
        result=Vector{F}(undef, L),
    )

    # CALCULATE GRADIENT FOR SIGNAL PARAMETERS
    t̄ = Integrations.lattice(grid)
    ∂̄ = array(FΩ, size(t̄), LABEL)
    Φ = (t, ∂, ϕα, ϕβ) -> (real(∂)*ϕα + imag(∂)*ϕβ)

    # The loop will add rather than set, so tare the result vector.
    result .= 0

    for i in 1:Devices.ndrives(device)
        j = 2i - 1
        ϕ̄α = @view(ϕ̄[:,j])
        ϕ̄β = @view(ϕ̄[:,j+1])

        signal = Devices.drivesignal(device, i)

        for k in 1:length(device.weights)
            # Normally fill ∂̄ with its partial derivative; here that is the normal pulse.
            normalpulse = signal.components[k].components[2]
            ∂̄ = Signals.valueat(normalpulse, t̄; result=∂̄)

            result[k] += Integrations.integrate(grid, Φ, ∂̄, ϕ̄α, ϕ̄β)
        end
    end

    return result
end




import CtrlVQE
function CtrlVQE.Systematic(
    TransmonDeviceType::Type{<:NormalTransmonDevice},
    n::Int,
    pulsesets;
    m=2,
    F=Float64,
)
    # DEFINE STANDARDIZED PARAMETERS
    ω0 = F(2π * 4.80)
    Δω = F(2π * 0.02)
    δ0 = F(2π * 0.30)
    g0 = F(2π * 0.02)

    # ASSEMBLE THE DEVICE
    ω̄ = collect(ω0 .+ (Δω * (1:n)))
    δ̄ = fill(δ0, n)
    ḡ = fill(g0, n-1)
    quples = [Quple(q,q+1) for q in 1:n-1]
    q̄ = 1:n
    ν̄ = copy(ω̄)
    weights = zeros(F, length(pulsesets))
    return TransmonDeviceType(ω̄, δ̄, ḡ, quples, q̄, ν̄, pulsesets, weights, m)
end


################################################################################
#= Tools to construct "natural" (ie. from Jacobian) pulse sets. =#

import CtrlVQE
import LinearAlgebra
import FiniteDifferences

function construct__staticevolver(evolution, device, basis, grid, ψ0)
    #= We need a function evolving on a device with a certain parameter set,
        but we can't mutate the device.
        So this function just un-mutates the parameter bind. =#
    return x ->(
        x0 = CtrlVQE.Parameters.values(device);
        CtrlVQE.Parameters.bind(device, x);
        ψ = CtrlVQE.evolve(evolution, device, basis, grid, ψ0);
        CtrlVQE.Parameters.bind(device, x0);
        ψ
    )
end

function construct__jacobian(evolver, device)
    cfd = FiniteDifferences.central_fdm(5, 1)
    x0 = CtrlVQE.Parameters.values(device)
    jac = FiniteDifferences.jacobian(cfd, evolver, x0)[1]
end

function construct__naturalpulsesets(device, jac)
    USV = LinearAlgebra.svd(jac)
    nS = count(x -> x > 1e-8, USV.S)
    nD = CtrlVQE.ndrives(device)

    pulsesets = Vector{SignalType{Float64,ComplexF64}}[] # TODO: generalize
    for l in 1:nS
        # DEEP COPY THE DEVICE, USING THE NATURAL PARAMETERS
        candidate = deepcopy(device)
        CtrlVQE.Parameters.bind(candidate, USV.V[:,l])

        # EXTRACT THE SIGNALS
        Ω̄ = [CtrlVQE.drivesignal(candidate, i) for i in 1:nD]
        push!(pulsesets, Ω̄)
    end
    return pulsesets
end



################################################################################
#= Proof-of-concept optimization. =#

# 1. Construct a regular device with reasonably-windowed pulse.
# 2. Construct the natural pulsesets.
# 3. Construct a normal device.    ...words are weird. This one is normal but not regular. ^_^
# 4. Construct all the other nonsense and run the optimization.




#= This is the standard workhorse script: complex pulses, fixed frequency. =#
import CtrlVQE

import UniformWindows_110123 as JOB
import UniformWindows_110123: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
protopulse = CtrlVQE.ComplexConstant(zero(Float),zero(Float))
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

##########################################################################################
#= Construct the normal device. =#

grid = CtrlVQE.TemporalLattice(setup.T, setup.r)
system = JOB.System(setup.code)

protodevice = CtrlVQE.Systematic(
    devicetype,
    system.n,
    CtrlVQE.UniformWindowed(protopulse, setup.T, setup.W),
    m = setup.m,
)

ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, protodevice)

fun = construct__staticevolver(evolution, protodevice, basis, grid, ψ0)
@time jac = construct__jacobian(fun, protodevice)
USV = LinearAlgebra.svd(jac)
pulsesets = construct__naturalpulsesets(protodevice, jac)

device = CtrlVQE.Systematic(NormalTransmonDevice, system.n, pulsesets; m=setup.m)

##########################################################################################
#= Construct the rest of the work variables. =#

O0 = CtrlVQE.QubitOperators.project(system.model.H, device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, grid, ψ0, O0),

    CtrlVQE.CostFunctionType{Float}[], Float[],
    # CtrlVQE.GlobalAmplitudeBound(device, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
    # TODO: This type assumes all signal parameters belong to signals.
    # CtrlVQE.GlobalFrequencyBound(device, grid, setup.ΔMAX, setup.λΔ, setup.σΔ),
    # TODO: This type assumes all device parameters belong to signals or to frequency.
)








##########################################################################################
#= Finalize the variable objects. =#

work = _!.work = JOB.WorkVars(
    system,
    grid,
    protopulse,
    device,
    normfn,
    lossfn,
)

# state = _!.state = isnothing(_!.state) ? JOB.initial_state() : _!.state
# TODO: hack
initstate = JOB.StateVars(
    zeros(Float, length(pulsesets)),
    Int[], Int[], Int[],
    LinearAlgebra.diagm(ones(Float, length(pulsesets))),
)
state = _!.state = isnothing(_!.state) ? initstate : _!.state
trace = _!.trace = isnothing(_!.trace) ? JOB.initial_trace() : _!.trace

CtrlVQE.Parameters.bind(device, _!.state.x)
lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)


!_!.run && error("Setup Finished")

##########################################################################################
#= Additional setup. =#

import Optim, LineSearches, Printf

# SELECT THE OPTIMIZATION METHOD
linesearch = LineSearches.MoreThuente()
bfgs_optimizer = Optim.BFGS(linesearch=linesearch)

# INITIALIZE THE OBJECTIVE FUNCTION
f = CtrlVQE.cost_function(lossfn)
g! = CtrlVQE.grad_function_inplace(lossfn)
bfgs_objective = Optim.OnceDifferentiable(f, g!, _!.state.x)

# INITIALIZE THE OPTIMIZATION STATE
bfgs_state = Optim.initial_state(
    bfgs_optimizer,
    Optim.Options(),
    bfgs_objective,
    _!.state.x,
)
    #= NOTE: Optim.Options() is a dummy argument, unused for BFGS state initialization. =#
bfgs_state.invH .= _!.state.Hk

iteration = isempty(_!.trace.iterations) ? Ref(0) : Ref(last(_!.trace.iterations))

# FORMALIZE THE CALLBACK
function bfgs_callback(callback_state)
    iteration[] += 1

    # PRINT ITERATION LINE
    println(Printf.format(
        Printf.Format("%8s    " * "%13.8g    "^3),
        callback_state.iteration,
        callback_state.value,
        callback_state.g_norm,
        callback_state.metadata["time"],
    ))

    # UPDATE TRACES
    JOB.update_trace!(iteration[], callback_state.value, callback_state.g_norm)

    # UPDATE DATA
    _!.state.x  .= bfgs_state.x
    _!.state.Hk .= bfgs_state.invH
    #= TODO: BFGS state also has x_previous, g_previous, etc.
        How much is needed to have identical restart?
        I think nothing, right? They are just used to calculate Hk for this iteration?
        Better check sometime...
    =#
    JOB.save()

    # MORE EXTENSIVE UPDATES EVERY SO OFTEN
    if iteration[] % _!.meta.update == 0
        JOB.report()
        JOB.plot_trace("")
        JOB.plot_pulse("", :amplitudes)
        JOB.archive(JOB.iterid(iteration[]))
    end

    return JOB.is_terminated()
end

# FORMALIZE THE OPTIONS OBJECT
bfgs_options = Optim.Options(
    f_tol = _!.meta.f_tol,
    g_tol = _!.meta.g_tol,
    iterations = _!.meta.maxiter,
    callback = bfgs_callback,
)

function do_optimization()
    # RUN OPTIMIZATION
    global bfgs_result = Optim.optimize(
        bfgs_objective,
        _!.state.x,
        bfgs_optimizer,
        bfgs_options,
        bfgs_state,
    )

    println(bfgs_result)

    return JOB.is_converged()
end

##########################################################################################
#= RUN OPTIMIZATION =#

!do_optimization() && println("Optimization did not converge. Run again to resume.")

JOB.save()
JOB.report()

JOB.plot_trace("")
JOB.plot_pulse("", :amplitudes)

JOB.archive("final")
JOB.plot_pulse("final")


