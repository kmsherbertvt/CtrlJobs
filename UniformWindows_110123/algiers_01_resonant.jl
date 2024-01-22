#= Match the device parameters Hisham pulled for Algiers qubits 0 and 1, on 1/11/24. =#

import CtrlVQE

import UniformWindows_110123 as JOB
import UniformWindows_110123: _!, Float

#= Generally setup is defined wholly separately,
    but for convenience and posterity we'll write this in.
    Usually it'll be commented out though.
=#
_!.setup.W = 3
_!.setup.ΩMAX = 2π*0.1275

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

device = CtrlVQE.FixedFrequencyTransmonDevice(
    2π .* [4.9462424483428675, 4.8360305756841],
    2π .* [0.3441379643980549, 0.34870027897863465],
    2π .* [0.0018554933031236817],
    [CtrlVQE.Quple(1,2)],
    [1, 2],
    2π .* [4.9462424483428675, 4.8360305756841], # drive matches resonance
    [CtrlVQE.CtrlVQE.UniformWindowed(deepcopy(protopulse), setup.T, setup.W) for _ in 1:2],
    setup.m,
)
#= The actual amplitude caps: 2π .* [0.1276456365439478, 0.13370758447967016]
    So, we'll be sure to set setup.Ω0 to, say, 2π*0.1275.
=#


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