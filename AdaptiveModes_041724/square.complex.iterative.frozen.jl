#=

- Candidates are square pulses with n windows.
- Add modes iteratively, so that in the end, the ansatz consists of all modes.
- Freeze parameters at each optimization.

=#

import CtrlVQE

import AdaptiveModes_041724 as JOB
import AdaptiveModes_041724: _!, Float

setup = _!.setup
meta  = _!.meta



##########################################################################################
#= Override adaptive functionality for frozen parameters. =#

""" Totally different and a bit simpler, since we mostly discard previous parameters. """
function frozen_adapt_state!(vars, modes, upHk)
    JOB.require_work(vars)

    # PARSE MODES ARGUMENT
    nD = CtrlVQE.ndrives(vars.work.device)
    n = [Int[] for i in 1:CtrlVQE.ndrives(vars.work.device)]
    LΩ = zeros(Int, CtrlVQE.ndrives(vars.work.device))
    for (i, n_) in (modes isa Pair{Int,Int} ? [modes] : modes)
        # Register the mode.
        push!(n[i], n_)
        # Tabulate the number of parameters for each drive.
        LΩ[i] += CtrlVQE.Parameters.count(vars.work.pool[n_])
    end
    L = sum(LΩ)     # Anachronistically, does NOT include frequency counts.

    # CONSTRUCT THE VECTORS OF ZEROS (plus carrying over frequency)
    x = [zeros(Float, L); vars.state.x[vars.state.ν]]
    imap = [zeros(Int, L); vars.state.ν]

    # "UPDATE" THE INVERSE HESSIAN
    Hk = upHk(vars.state.Hk, imap)

    # UPDATE THE INDEXING VECTORS
    Ω = [collect(1+sum(LΩ[1:i]):sum(LΩ[1:i+1])) for i in 0:nD-1]
    ν = isempty(vars.state.ν) ? Int[] : collect(1+L:1+length(vars.state.ν))

    vars.state = JOB.StateVars(x, Hk, n, Ω, ν)
end
frozen_adapt_state!(modes, upHk) = frozen_adapt_state!(_!, modes, upHk)

function frozen_adapt_work!(vars)
    JOB.require_work(vars)
    pool = vars.work.pool

    # MUTATE THE DEVICE TO ACCOUNT FOR NEW PULSE SHAPES
    for i in 1:CtrlVQE.ndrives(vars.work.device)
        components = [deepcopy(pool[n]) for n in vars.state.n[i]]
        signal = CtrlVQE.CompositeSignal(deepcopy(vars.work.protopulse), components...)
        CtrlVQE.set_drivesignal(vars.work.device, i, signal)
    end
    CtrlVQE.Parameters.bind(vars.work.device, vars.state.x)
    #= NOTE: This line handles frequencies. The next section SHOULD never touch them. =#

    # NOW ADD IN ALL THE PULSES FROM PREVIOUS RUNS, FROZEN
    currentstate = vars.state
    for a in 1:length(vars.trace.adaptations)
        # LOAD THE PREVIOUS STATE FROM A FILE
        unarchive!(vars, JOB.adaptid(a))

        # ASSUME PARAMETERS ARE ORDERED drive -> mode -> signal parameters
        offset = 0
        for i in 1:CtrlVQE.ndrives(vars.work.device)
            for n in vars.state.n[i]
                # Create the mode with the optimized parameter(s).
                signal = deepcopy(pool[n])
                L = CtrlVQE.Parameters.count(signal)
                CtrlVQE.Parameters.bind(signal, vars.state.x[1+offset:L+offset])
                offset += L

                # Re-create this signal as a totally constrained harmonic.
                #= NOTE: This assumes each pool element is a WindowedSignal of constants.
                        If we adapt this script to other pools,
                            I think we'll have no recourse but to edit this line. :/ =#
                frozensignal = CtrlVQE.WindowedSignal([
                    CtrlVQE.ConstrainedSignal(window, :A, :B)
                        for window in signal.windows
                ], signal.starttimes)

                # Add this frozen signal to the active device.
                push!(CtrlVQE.drivesignal(device, i).components, frozensignal)
            end
        end
    end
    vars.state = currentstate
end
frozen_adapt_work!() = frozen_adapt_work!(_!)

""" Count parameters cumulatively, since the instantaneus count is uninformative. """
function frozen_adapt_trace!(vars, scores)
    G_max = isempty(scores) ? 0.0 : maximum(scores)
    push!(vars.trace.adaptations, last(vars.trace.iterations))
    push!(vars.trace.poolsize, length(scores))
    push!(vars.trace.G_max, G_max)

    #= Let's define, for frozen runs, parameters should count CUMULATIVE parameters,
        rather than the current number (which would just always be the same). =#
    LΩ = [length(Ω) for Ω in vars.state.Ω]
    currentstate = vars.state
    for a in 1:length(vars.trace.adaptations)-1     # NOTHING HAPPENS IN 1st/2nd ADAPTs
        unarchive!(vars, JOB.adaptid(a))
        LΩ .+= [length(Ω) for Ω in vars.state.Ω]
    end
    vars.state = currentstate
    push!(vars.trace.parameters, LΩ)
end
frozen_adapt_trace!(scores) = frozen_adapt_trace!(_!, scores)

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_complexsquare
select = JOB.select_iterative
    # HACK: Alas, `select_iterative` needs a hack at the end of the adapt loop also. =#
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
    frozen_adapt_work!()

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
        frozen_adapt_trace!(_!, scores)
    end
    JOB.adapt_is_terminated() && break
        # NOTE: If pool is empty, update sets G_max to 0.0, which always terminates.

    # PERFORM ADAPTATION
    modes = select(scores)
    JOB.adapt_is_terminated() && break
        # HACK: select_iterative hacks in the last G_max to 0.0, so *now* it terminates.
    frozen_adapt_state!(modes, upHk)
    frozen_adapt_work!()

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
