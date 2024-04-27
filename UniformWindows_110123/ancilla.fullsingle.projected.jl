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

ψ0 = CtrlVQE.QubitOperators.project(CtrlVQE.kron(       # Extend statevector on ancillae.
    system.ψ0,
    CtrlVQE.LinearAlgebraTools.basisvector(2, 1),
), device)
O0 = CtrlVQE.QubitOperators.project(CtrlVQE.kron(       # Extend Hamiltonian on ancillae.
    system.model.H,
    [1 0; 0 0],             # NOTE: Project out the ancilla.
), device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, grid, ψ0, O0),
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
JOB.plot_pulse("final")