#=

There are a LOT of things you can do in GMQC by attaching a "second line" of qubits.

    o-o-o-o     ->      o-o-o-o
                        | | | |
                        o-o-o-o

Error correction, virtual distillation, controlled multi-qubit operations, etc.

Can it also serve a purpose in ctrl-VQE?

Two ideas of how it might:
1. Larger space allowing "hyperspace" travel to the ground state, just like leakage.
2. Normally entanglement with ancilla will be implicitly penalized by an unnormalized energy. But this is not so in gapless systems where the ground state is degenerate. So...perhaps purity of an optimized solution can witness tightly gapped systems! Optimization pry still runs into serious issues here, though. Need to think about this more.

Anyway, we'll try the 2xN architecture in the `ancilla` script,
    and in `ancilla.linear` we'll use a "control" where the ancilla are just linearly appended to the system register.

In THIS script, we'll try to get some hyperspace action by letting coupled qubits be a potentially more direct path. So for now couple one ancilla with all data qubits. But project out the |1⟩ states on the ancilla. Actually that's what we should be doing in all of these scripts probably...

And it turns out I've been misunderstanding something major about how the speed-up works lol; we need to re-normalize. Which means we'll need to manually construct a new CostFunction which calculates the "projected energy" and the "normalization" just like "NormalizedEnergy" but with a subset of qubits as the computational space. Whew! TODO

=#

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
protopulse = CtrlVQE.ComplexConstant(zero(Float), zero(Float))
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

##########################################################################################
#= Define a specialized cost-function projecting and normalizing out ancillae. =#

import CtrlVQE: CostFunctions

import CtrlVQE: LinearAlgebraTools, QubitOperators
import CtrlVQE: Parameters, Integrations, Devices, Evolutions
import CtrlVQE: Bases, Operators

import CtrlVQE.TrapezoidalIntegrations: TrapezoidalIntegration

struct AncillaNormalizedEnergy{F} <: CostFunctions.EnergyFunction{F}
    evolution::Evolutions.EvolutionType
    device::Devices.DeviceType
    basis::Bases.BasisType
    frame::Operators.StaticOperator
    grid::TrapezoidalIntegration
    ψ0::Vector{Complex{F}}
    O0::Matrix{Complex{F}}
    nH::Int

    function AncillaNormalizedEnergy(
        evolution::Evolutions.EvolutionType,
        device::Devices.DeviceType,
        basis::Bases.BasisType,
        frame::Operators.StaticOperator,
        grid::TrapezoidalIntegration,
        ψ0::AbstractVector,
        O0::AbstractMatrix,
        nH::Integer,
    )
        # INFER FLOAT TYPE AND CONVERT ARGUMENTS
        F = real(promote_type(Float16, eltype(O0), eltype(ψ0), eltype(grid)))

        # CREATE OBJECT
        return new{F}(
            evolution, device, basis, frame, grid,
            convert(Array{Complex{F}}, ψ0),
            convert(Array{Complex{F}}, O0),
            nH
        )
    end
end

Base.length(fn::AncillaNormalizedEnergy) = Parameters.count(fn.device)

function CostFunctions.trajectory_callback(
    fn::AncillaNormalizedEnergy,
    En::AbstractVector;
    callback=nothing
)
    workbasis = Evolutions.workbasis(fn.evolution)      # BASIS OF CALLBACK ψ
    U = Devices.basisrotation(fn.basis, workbasis, fn.device)
    π̄ = QubitOperators.localqubitprojectors(fn.device)
    #= --- SPECIAL TO ANCILLAE --- =#
    n = size(π̄,3); for q in 1+n-nH:n; π̄[2,2,q] = 0; end
        # Qubit projectors onto each ancilla will use |0⟩⟨0| instead of |0⟩⟨0|+|1⟩⟨1|
    ψ_ = similar(fn.ψ0)

    return (i, t, ψ) -> (
        ψ_ .= ψ;
        LinearAlgebraTools.rotate!(U, ψ_);  # ψ_ IS NOW IN MEASUREMENT BASIS
        LinearAlgebraTools.rotate!(π̄, ψ_);  # ψ_ IS NOW "MEASURED"
        # APPLY FRAME ROTATION TO STATE RATHER THAN OBSERVABLE
        Devices.evolve!(fn.frame, fn.device, fn.basis, -t, ψ_);
            # NOTE: Rotating observable only makes sense when t is always the same.
        E = real(LinearAlgebraTools.expectation(fn.O0, ψ_));
        F = real(LinearAlgebraTools.expectation(π̄, ψ_));
        En[i] = E / F;
        !isnothing(callback) && callback(i, t, ψ)
    )
end

function CostFunctions.cost_function(fn::AncillaNormalizedEnergy; callback=nothing)
    # THE PROJECTION OPERATOR
    π̄ = QubitOperators.localqubitprojectors(fn.device)
    #= --- SPECIAL TO ANCILLAE --- =#
    n = size(π̄,3); for q in 1+n-nH:n; π̄[2,2,q] = 0; end
        # Qubit projectors onto each ancilla will use |0⟩⟨0| instead of |0⟩⟨0|+|1⟩⟨1|
    # DYNAMICALLY UPDATED STATEVECTOR
    ψ = copy(fn.ψ0)
    # OBSERVABLE, IN MEASUREMENT FRAME
    T = Integrations.endtime(fn.grid)
    OT = copy(fn.O0); Devices.evolve!(fn.frame, fn.device, fn.basis, T, OT)
    # INCLUDE PROJECTION ONTO COMPUTATIONAL SUBSPACE IN THE MEASUREMENT
    LinearAlgebraTools.rotate!(π̄, OT)

    return (x̄) -> (
        Parameters.bind(fn.device, x̄);
        Evolutions.evolve(
            fn.evolution,
            fn.device,
            fn.basis,
            fn.grid,
            fn.ψ0;
            result=ψ,
            callback=callback,
        );
        E = real(LinearAlgebraTools.expectation(OT, ψ));
        F = real(LinearAlgebraTools.expectation(π̄, ψ));
        E / F
    )
end

function CostFunctions.grad_function_inplace(fn::AncillaNormalizedEnergy{F}; ϕ=nothing) where {F}
    r = Integrations.nsteps(fn.grid)

    if isnothing(ϕ)
        return CostFunctions.grad_function_inplace(
            fn;
            ϕ=Array{F}(undef, r+1, Devices.ngrades(fn.device), 2)
        )
    end

    # THE PROJECTION OPERATOR, FOR COMPONENT COST FUNCTION EVALUATIONS
    π̄ = QubitOperators.localqubitprojectors(fn.device)
    #= --- SPECIAL TO ANCILLAE --- =#
    n = size(π̄,3); for q in 1+n-nH:n; π̄[2,2,q] = 0; end
        # Qubit projectors onto each ancilla will use |0⟩⟨0| instead of |0⟩⟨0|+|1⟩⟨1|

    # DYNAMICALLY UPDATED STATEVECTOR
    ψ = copy(fn.ψ0)

    # THE "MATRIX LIST" (A 3D ARRAY), FOR EACH GRADIENT SIGNAL
    Ō = Array{eltype(ψ)}(undef, (size(fn.O0)..., 2))
    # FIRST MATRIX: THE OBSERVABLE, IN MEASUREMENT FRAME
    OT = @view(Ō[:,:,1])
    T = Integrations.endtime(fn.grid)
    OT .= fn.O0; Devices.evolve!(fn.frame, fn.device, fn.basis, T, OT)
    # INCLUDE PROJECTION ONTO COMPUTATIONAL SUBSPACE IN THE MEASUREMENT
    LinearAlgebraTools.rotate!(π̄, OT)
    # SECOND MATRIX: PROJECTION OPERATOR, AS A GLOBAL OPERATOR
    LinearAlgebraTools.kron(π̄; result=@view(Ō[:,:,2]))

    # GRADIENT VECTORS
    ∂E = Array{F}(undef, length(fn))
    ∂N = Array{F}(undef, length(fn))

    return (∇f̄, x̄) -> (
        Parameters.bind(fn.device, x̄);
        Evolutions.evolve(
            fn.evolution,
            fn.device,
            fn.basis,
            fn.grid,
            fn.ψ0;
            result=ψ,
        );
        E = real(LinearAlgebraTools.expectation(OT, ψ));
        N = real(LinearAlgebraTools.expectation(π̄, ψ));

        Parameters.bind(fn.device, x̄);
        Evolutions.gradientsignals(
            fn.evolution,
            fn.device,
            fn.basis,
            fn.grid,
            fn.ψ0,
            Ō;
            result=ϕ,   # NOTE: This writes the gradient signal as needed.
        );
        ∂E .= Devices.gradient(fn.device, fn.grid, @view(ϕ[:,:,1]));
        ∂N .= Devices.gradient(fn.device, fn.grid, @view(ϕ[:,:,2]));

        ∇f̄ .= (∂E./N) .- (E/N) .* (∂N./N)
    )
end

##########################################################################################
#= Prepare all the work variables. =#

system = JOB.System(setup.code)
grid = CtrlVQE.TemporalLattice(setup.T, setup.r)

device = CtrlVQE.Systematic(
    devicetype,
    1+system.n,                         # NOTE: Use one extra qubit.
    CtrlVQE.UniformWindowed(protopulse, setup.T, setup.W),
    m = setup.m,
)
append!(device.quples, CtrlVQE.Quple(q,1+system.n) for q in 1:system.n)
append!(device.ḡ, device.ḡ[1:system.n])
    # Attaches all data qubits to an ancilla.

nH = CtrlVQE.nqubits(device) - system.n

import LinearAlgebra
ψ0 = CtrlVQE.QubitOperators.project(CtrlVQE.kron(       # Extend statevector on ancillae.
    system.ψ0,
    CtrlVQE.LinearAlgebraTools.basisvector(setup.m, 1),
), device)
O0 = CtrlVQE.QubitOperators.project(CtrlVQE.kron(       # Extend Hamiltonian on ancillae.
    system.model.H,
    Matrix{Bool}(LinearAlgebra.I, setup.m, setup.m),
), device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    AncillaNormalizedEnergy(evolution, device, basis, frame, grid, ψ0, O0, nH),
    CtrlVQE.GlobalAmplitudeBound(device, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
    CtrlVQE.GlobalFrequencyBound(device, grid, setup.ΔMAX, setup.λΔ, setup.σΔ),
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

state = _!.state = isnothing(_!.state) ? JOB.initial_state() : _!.state
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
JOB.plot_pulse("final", :moduli, :phases, :trajectory)