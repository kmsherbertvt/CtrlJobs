#=

For the sake of comparison, define a not-exactly-adaptive trajectory
    which adds one window to each pulse at each "adaptation"
    (so it is directly comparable to "each" runs),
    by simply re-dividing each pulse into so-many chunks.
Keep going until so-many chunks require windows less than the usual ΔsMIN.

=#

import CtrlVQE

import AdaptiveWindows_040624 as JOB
import AdaptiveWindows_040624: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
protopulse = CtrlVQE.ComplexConstant(zero(Float), zero(Float))
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_bisection
select = JOB.select_one
upHk = JOB.upHk_slate

##########################################################################################
#= Prepare all the work variables. =#

system = JOB.System(setup.code)
grid = CtrlVQE.TemporalLattice(setup.T, setup.r)

device = CtrlVQE.Systematic(
    devicetype,
    system.n,
    CtrlVQE.UniformWindowed(protopulse, setup.T, 1),
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

JOB.adapt_work!()       # SYNC WORK OBJECTS WITH THE CURRENT STATE

!_!.run && error("Setup Finished")

##########################################################################################
#= Additional setup. =#

import Optim, LineSearches

function do_optimization()
    # ENSURE WORK OBJECTS ARE CONSISTENT WITH PARAMETERS
    JOB.adapt_work!()

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

""" Update trace object with adaptive information.

Poolsize and G_max do not have their usual meaning in this trajectory,
    so we'll use poolsize for the total number of adaptations (constant throughout trajectory)
    and G_max for the number of adaptations remaining,
        so that G_max == 0.0 at convergence.

"""
function my_adapt_trace!(vars)
    # COUNT HOW MANY ADAPTATIONS ARE LEFT TO GO
    ΔsMIN = 1 / vars.setup.fMAX
    WMAX = floor(Int, vars.setup.T / ΔsMIN)
    WNOW = length(first(vars.state.s))
    cnt = WMAX - WNOW

    # UPDATE ADAPTIVE TRACE VARS
    push!(vars.trace.adaptations, last(vars.trace.iterations))
    push!(vars.trace.poolsize, WMAX)
    push!(vars.trace.G_max, Float(cnt))
    push!(vars.trace.parameters, [length(Ω) for Ω in vars.state.Ω])
    push!(vars.trace.Δs_min, [minimum(diff([s; vars.setup.T])) for s in vars.state.s])
end

""" Update state objects with one or more new segments.

The `segments` argument may be a tuple or a vector of tuples `(i,s)`,
    where `i` is the drive index and `s` is the starttime of a new window.
"""
function my_adapt_state!(vars)
    # PREPARE A CANDIDATE DEVICE (which knows all about its own parameters)
    WNOW = length(first(vars.state.s))
    candidate = CtrlVQE.Systematic(
        devicetype,
        vars.work.system.n,
        CtrlVQE.UniformWindowed(vars.work.protopulse, setup.T, 1+WNOW),
        m = setup.m,
    )

    # FETCH THE NEW PARAMETER VECTOR
    x = CtrlVQE.Parameters.values(candidate)    # All zeros.
    L = CtrlVQE.Parameters.count(candidate)

    # UPDATE THE INVERSE HESSIAN
    imap = 1:L
    Hk = JOB.upHk_naive(vars.state.Hk, imap)

    # FETCH WINDOW START TIMES
    n = CtrlVQE.ndrives(candidate)
    signals = [CtrlVQE.drivesignal(candidate, i) for i in 1:n]
    s = [signal.starttimes for signal in signals]

    # PREPARE SIGNAL INDEXING VECTORS
    LΩ = [CtrlVQE.Parameters.count(signal) for signal in signals]
    Ω = [collect(1+sum(LΩ[1:i]):sum(LΩ[1:i+1])) for i in 0:n-1]

    # PREPARE FREQUENCY INDEXING VECTORS -- Implementation may vary with type someday?
    ν = vars.state.ν

    vars.state = JOB.StateVars(x, Hk, s, Ω, ν)
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

    # UPDATE TRACE
    if loaded_converged
        global loaded_converged = false
    else
        my_adapt_trace!(_!)
    end
    JOB.adapt_is_terminated() && break
        # NOTE: If pool is empty, update sets G_max to 0.0, which always terminates.

    # PERFORM ADAPTATION
    my_adapt_state!(_!)
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
