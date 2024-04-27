#= Match the device parameters provided by AQT folks in mid-March 2024. =#

import CtrlVQE

import UniformWindows_110123 as JOB
import UniformWindows_110123: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Define a "master" device, from which will be constructed the n-sized device. =#
doi_10_1103 = CtrlVQE.TransmonDevice(
    2π .* [5.230708, 5.297662, 5.459108, 5.633493],
    2π .* [0.27366, 0.27332, 0.27070, 0.26745],
        # TODO: Verify AQT's definition of anharmonicity. There's a question of -1 and 2.
    2π .* [0.0025, 0.00273, 0.00311],
    [CtrlVQE.Quple(1,2), CtrlVQE.Quple(2,3), CtrlVQE.Quple(3,4)],
    [1, 2, 3, 4],
    2π .* [5.230708, 5.297662, 5.459108, 5.633493], # Initial drive matches resonance
    [CtrlVQE.Constant(0.0) for _ in 1:4],    # Just a placeholder.
    2,                                              # Just a placeholder.
)

"""
    subdevice(devicetype, master, pulses, qubits)

Construct one transmon-like device from another, selecting a subset of qubits.

# Arguments
- `devicetype`: the type of transmon device
- `template`: the original device from which to take qubits
- `qubits`: the list of qubits to keep from template

# Keyword Arguments
- `pulses`: pulses to use for the new device
- `m`: truncation level

"""
function subdevice(devicetype, template, qubits; pulses, m)
    # IDENTIFY WHICH COUPLINGS ARE RETAINED
    qupleindices = [
        i for i in eachindex(template.quples)
            if template.quples[i].q1 in qubits && template.quples[i].q2 in qubits
    ]
    # RELABEL COUPLINGS WITH THE NEW QUBIT ORDERING
    quples = [
        CtrlVQE.Quple(findfirst(==(quple.q1), qubits), findfirst(==(quple.q2), qubits))
            for quple in template.quples[qupleindices]
    ]
    # CONSTRUCT THE NEW DEVICE
    return devicetype(
        template.ω̄[qubits],
        template.δ̄[qubits],
        template.ḡ[qupleindices],
        quples,
        template.q̄[qubits],
        template.ν̄[qubits],
        pulses,
        m,
    )
end

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
protopulse = CtrlVQE.ComplexConstant(zero(Float), zero(Float))
devicetype = CtrlVQE.TransmonDevice

##########################################################################################
#= Prepare all the work variables. =#

system = JOB.System(setup.code)
grid = CtrlVQE.TemporalLattice(setup.T, setup.r)

pulses =
device = subdevice(
    devicetype, doi_10_1103, 1:system.n;
    pulses=[
        CtrlVQE.CtrlVQE.UniformWindowed(deepcopy(protopulse), setup.T, setup.W)
            for _ in 1:system.n
    ],
    m=setup.m,
)

ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, device)
O0 = CtrlVQE.QubitOperators.project(system.model.H, device)
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