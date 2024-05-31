#=

In each adaptation, add the current gradient signal to the basis.

This version does not require "scoring" anything.
You just make the gradient signals,
    then convert them into (constrained) windowed pulses,
    then add it into the appropriate drive signal.

Robustly "adapting" the device requires loading each adaptation,
    measuring the gradient signals, and adding them to the appropriate drive signal.

=#

import CtrlVQE

import AdaptiveModes_041724 as JOB
import AdaptiveModes_041724: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Some custom adapt methods. =#

import LinearAlgebra
""" Exact gradients don't actually need scores at all.
    The "pool" has exactly one element in it, so we don't need to "select" anything.

But, we _do_ need to know when to _stop_.

So, we should compute a single-element score vector
    which is the norm of the initial gradient when adding in the new direction.

For the energy alone, one can simply compute ∑⟨ϕ,ϕ⟩, I think.

But we have to account for the bloody penalties.

So actually we're going to have to construct a candidate device...

"""
function gradients_score_candidates(vars, ϕ)
    JOB.require_work(vars)
    lossfn = vars.work.lossfn

    # SPLIT ϕ UP INTO ϕα, ϕβ
    nD = CtrlVQE.ndrives(vars.work.device)
    ϕα = Array{Float}(undef, vars.setup.r+1, nD)
    ϕβ = Array{Float}(undef, vars.setup.r+1, nD)
    for i in 1:nD
        j = 2(i-1) + 1
        ϕα[:,i] .= ϕ[:,j  ]
        ϕβ[:,i] .= ϕ[:,j+1]
    end

    # PREPARE THE CANDIDATE DEVICE
    candidate = deepcopy(vars.work.device)
    JOB.CoupledDevices.add_signals!(candidate, [
        JOB.CoupledDevices.as_trotterized_signal(ϕα[:,i], ϕβ[:,i], vars.work.grid)
            for i in 1:CtrlVQE.ndrives(candidate)
    ])
    x = CtrlVQE.Parameters.values(candidate)



    scores = Matrix{Float}(undef, nD, 1)

    # CALCULATE GRADIENT
    g = CtrlVQE.Devices.gradient(candidate, vars.work.grid, ϕ)
    for l in 1:length(lossfn.penalties)
        # We need to extend the penalty function to work on more parameters.
        # HACK: Construct a new fn using each field except one called "device".
        protofn = lossfn.penaltyfns[l]
        fntype = Base.typename(typeof(protofn)).wrapper
            # This "typename" "wrapper" strips away parametric args in the type,
            #   so that Julia can find the generic constructor.
        penaltyfn = fntype((
            field == :device ? candidate : getfield(protofn, field)
                for field in fieldnames(fntype)
        )...)

        gd = CtrlVQE.grad_function(penaltyfn)

        # MODIFY THE GRADIENT WITH PENALTY CONTRIBUTION
        g .+= lossfn.weights[l] * gd(x)
    end

    scores .= LinearAlgebra.norm(g)

    return scores
end
gradients_score_candidates(args...) = gradients_score_candidates(_!, args...)

""" Just inserts one more Ω parameter, starting from zero. """
function gradients_adapt_state!(vars)
    JOB.require_work(vars)

    imap = collect(1:length(_!.state.x))

    LΩ = JOB.CoupledDevices.ndriveparams(vars.work.device)

    x = deepcopy(_!.state.x)
    insert!(x, 1+LΩ, zero(eltype(_!.state.x)))
    insert!(imap, 1+LΩ, 0)

    Hk = upHk(vars.state.Hk, imap)

    n = _!.state.n  # NOT USED

    Ω = deepcopy(_!.state.Ω); foreach(Ωi -> push!(Ωi, LΩ), Ω)
    ν = 1 .+ _!.state.ν

    vars.state = JOB.StateVars(x, Hk, n, Ω, ν)
end
gradients_adapt_state!(args...) = gradients_adapt_state!(_!, args...)

function gradients_adapt_work!(vars)
    JOB.require_work(vars)
    device = vars.work.device

    # ITERATIVELY ADD IN EACH GRADIENT SIGNAL FROM EACH ADAPTATION
    #= NOTE: This is more than is necessary when running a script;
        one SHOULD do this once at the beginning,
        and from then on just add the one ϕ as you finish each optimization.
        But, it's easier to just write this one function that does it all each time. ^_^
    =#
    currentstate = vars.state
    JOB.CoupledDevices.empty_signals!(device)
    # for a in 1:length(vars.trace.adaptations)   # TODO: actually just up to length of x?
    for a in 1:length(vars.state.x)
        # LOAD THE PREVIOUS STATE FROM A FILE
        unarchive!(vars, JOB.adaptid(a))
        CtrlVQE.Parameters.bind(device, vars.state.x)

        # CALCULATE THE GRADIENT SIGNAL
        ϕ, ϕα, ϕβ = JOB.make_gradientsignals()

        # NORMALIZE SO THAT MAX IS 1
        A = maximum(abs.(ϕ))
        ϕα ./= A
        ϕβ ./= A

        # ADD IN THE GRADIENT SIGNALS AS BASIS PULSES TO THE DEVICE
        JOB.CoupledDevices.add_signals!(device, [
            JOB.CoupledDevices.as_trotterized_signal(ϕα[:,i], ϕβ[:,i], vars.work.grid)
                for i in 1:CtrlVQE.ndrives(device)
        ])
    end
    vars.state = currentstate
    CtrlVQE.Parameters.bind(device, vars.state.x)
end
gradients_adapt_work!(args...) = gradients_adapt_work!(_!, args...)

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_harmonics
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

protodevice = CtrlVQE.Systematic(
    devicetype,
    system.n,
    protopulse,
    m = setup.m,
)

device = JOB.CoupledDevice(protodevice)

ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, device)
O0 = CtrlVQE.QubitOperators.project(system.model.H, device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, grid, ψ0, O0),
    JOB.CoupledGlobalAmplitudeBound(device, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
    JOB.CoupledGlobalFrequencyBound(device, grid, setup.ΔMAX, setup.λΔ, setup.σΔ),
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

gradients_adapt_work!()
lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)

!_!.run && error("Setup Finished")

##########################################################################################
#= Additional setup. =#

import Optim, LineSearches

function do_optimization()
    # ENSURE WORK OBJECTS ARE CONSISTENT WITH PARAMETERS
    gradients_adapt_work!()

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
    # global linesearch = LineSearches.MoreThuente()
    global linesearch = LineSearches.BackTracking()
        #= TODO: MoreThuente's H215 has ..* infinities after normalization.
        Using BackTracking, we don't run into infinities.
        But the linesearch does blow up at some point.
        We need to load those points and double-check the gradient.
        There's pry something wrong with the CoupledGlobalAmplitudeBound.

        Otherwise, maybe go back to MoreThuente and try with some other InitialGuess..?
        =#
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
    JOB.plot_trace("");
    flush(stdout)

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
loaded_converged && gradients_adapt_work!()

while loaded_converged || do_optimization()
    JOB.save()
    JOB.report()

    # ARCHIVE THE OPTIMIZED STATE
    adaptiter = 1 + length(_!.trace.adaptations)
    name = JOB.adaptid(adaptiter)
    JOB.archive(name)
    JOB.plot_pulse(name, :amplitudes, :gradients)

    # CALCULATE GRADIENT SIGNALS
    ϕ, ϕα, ϕβ = JOB.make_gradientsignals()

    # CONSTRUCT POOL AND CALCULATE SCORES
    scores = gradients_score_candidates(ϕ)

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
    gradients_adapt_state!()
    gradients_adapt_work!()

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
