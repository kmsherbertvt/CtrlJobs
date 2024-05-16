#=

- Candidates are harmonic modes with nodes at the start and end (aka sine waves).
- Add one mode at a time, including real and imaginary parameters.
- Hessian is identity on the updated parameters only.

=#

import CtrlVQE

import AdaptiveModes_041724 as JOB
import AdaptiveModes_041724: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_complexharmonics
select = JOB.select_one
upHk = JOB.upHk_slate

##########################################################################################
#= Prepare all the work variables. =#

system = JOB.System(setup.code)
grid = CtrlVQE.TemporalLattice(setup.T, setup.r)

protopulse = CtrlVQE.ConstrainedSignal(
    CtrlVQE.ComplexConstant(zero(Float), zero(Float)),
    :A, :B,
)

device = CtrlVQE.Systematic(
    devicetype,
    system.n,
    CtrlVQE.CompositeSignal(protopulse),
    m = setup.m,
)

ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, device)
O0 = CtrlVQE.QubitOperators.project(system.model.H, device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, grid, ψ0, O0),
    CtrlVQE.GlobalAmplitudeBound(device, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
    CtrlVQE.GlobalFrequencyBound(device, grid, setup.ΔMAX, setup.λΔ, setup.σΔ),
)

pool = makepool(_!)

##########################################################################################
#= Finalize the variable objects. =#

work = _!.work = JOB.WorkVars(
    system,
    grid,
    protopulse,
    device,
    normfn,
    lossfn,
    pool,
)

state = _!.state = isnothing(_!.state) ? JOB.initial_state() : _!.state
trace = _!.trace = isnothing(_!.trace) ? JOB.initial_trace() : _!.trace

CtrlVQE.Parameters.bind(device, _!.state.x)
lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)

!_!.run && error("Setup Finished")

##########################################################################################
#= Additional setup. =#

import Optim, LineSearches

function do_optimization()
    # ENSURE WORK OBJECTS ARE CONSISTENT WITH PARAMETERS
    JOB.adapt_work!()

    # SEED THE VERY FIRST ITERATION WITH INERT NUMBERS
    if isempty(_!.trace.iterations)
        fn = _!.work.lossfn(Float[])
        gd = zero(Float)
        JOB.update_trace!(1, fn, gd)
    end

    # DON'T OPTIMIZE AN EMPTY PARAMETER VECTOR
    #= Note this ALSO happens for the very first iteration,
        but it needs to be a separate check to be robust against restarts. =#
    if isempty(_!.state.x)
        return true
    end

    # SPECIFY OPTIMIZATION ROUTINE
    global linesearch = LineSearches.MoreThuente()
    global bfgs_optimizer = Optim.BFGS(linesearch=linesearch)

    # PREPARE OPTIMIZATION WORK VARIABLES
    global bfgs_objective = JOB.make_objective()
    global bfgs_state = JOB.make_state(bfgs_optimizer, bfgs_objective)
    global bfgs_callback = JOB.make_callback(bfgs_state)
    global bfgs_options = JOB.make_options(bfgs_callback)

    # RUN OPTIMIZATION
    global bfgs_result = Optim.optimize(
        bfgs_objective,
        _!.state.x,
        bfgs_optimizer,
        bfgs_options,
        bfgs_state,
    )

    # REPORT RESULTS
    println(bfgs_result)

    return JOB.optimization_is_converged()
end

function trace_is_optimized()
    isnothing(_!.trace) && return false
    isempty(_!.trace.iterations) && return false
    isempty(_!.trace.adaptations) && return false
    return last(_!.trace.iterations) == last(_!.trace.adaptations)
end

##########################################################################################
#= RUN ADAPT =#

loaded_converged = trace_is_optimized()

while loaded_converged || do_optimization()
    JOB.save()
    JOB.report()

    # ARCHIVE THE OPTIMIZED STATE
    adaptiter = 1 + length(_!.trace.adaptations)
    name = JOB.adaptid(adaptiter)
    JOB.archive(name)
    JOB.plot_pulse(name, :amplitudes)

    # CALCULATE GRADIENT SIGNALS
    ϕ, ϕα, ϕβ = JOB.make_gradientsignals()

    # CONSTRUCT POOL AND CALCULATE SCORES
    scores = JOB.score_candidates(ϕ)

    # UPDATE TRACE
    if loaded_converged
        global loaded_converged = false
    else
        JOB.adapt_trace!(_!, scores)
    end
    JOB.adapt_is_terminated() && break
        # NOTE: If pool is empty, update sets G_max to 0.0, which always terminates.

    # PERFORM ADAPTATION
    modes = select(scores)
    JOB.adapt_state!(modes, upHk)
    JOB.adapt_work!()

end
JOB.save()
JOB.report()

JOB.plot_trace("")
JOB.plot_pulse("", :amplitudes)

if trace_is_optimized()
    JOB.archive("final")
    JOB.plot_pulse("final")
    println("ctrl-ADAPT-VQE has terminated.")
else
    println("Optimization did not converge. Run again to resume.")
end
