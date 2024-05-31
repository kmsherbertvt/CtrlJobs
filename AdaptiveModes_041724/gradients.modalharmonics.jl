#=

In each adaptation, add a linear combination of all basis vectors, weighted by score.
NOTE: score here should NOT be |⋅|'d.

This version does not require "scoring" anything.
You just make the gradient signals,
    then convert them into (constrained) windowed pulses,
    then add it into the appropriate drive signal.

This file is distinct from `gradients.harmonics`
    because we are using the super-customized ModalHarmonics
    rather than the hypothetically-extensible CoupledDevices.
They are explicitly fixed frequency, for now.

Also, at each adaptation, we'll manually store the array of coefficients with NPZ.
That way, syncing to a particular adaptation
    just requires taking the first so many rows of the saved matrix,
    rather than loading the state of each adaptation in turn
    and recomputing gradient signals.

=#

import NPZ

import CtrlVQE

import AdaptiveModes_041724 as JOB
import AdaptiveModes_041724: _!, Float

setup = _!.setup
meta  = _!.meta

##########################################################################################
#= Some custom adapt methods. =#

function prepare_candidate(vars, mode::Pair{Int,Int})
    JOB.require_work(vars)
    pool = vars.work.pool
    i, n = mode

    nD = CtrlVQE.ndrives(vars.work.device)

    scores = zeros(JOB.Float, length(pool), nD)
    scores[n,i] = 1

    candidate = deepcopy(vars.work.device)
    append!(candidate, scores)

    return candidate, nothing
end

""" Calculate the total gradient norm for each candidate in a pool.

The major differences here are:
1. This script only makes sense when each basis pulse has a single parameter.
2. We actually need the scores to be unnormalized, ie. retain sign.
3. Actually turns out better to have index [n,i] not [i,n]. Obvious in retrospect.

Thus, everything proceeds normally to compute `g`,
    the gradient vector of the new device with respect to a given pulse,
    except that at the end, we take `g[end]` rather than `LinearAlgebra.norm(g)`.

"""
function score_candidates(vars, ϕ)
    JOB.require_work(vars)
    pool = vars.work.pool
    lossfn = vars.work.lossfn

    nD = CtrlVQE.ndrives(vars.work.device)
    scores = Matrix{Float}(undef, length(pool), nD)
    for i in 1:nD; for n in eachindex(pool)
        # PREPARE CANDIDATE DEVICE AND EXTENDED PARAMETER VECTOR
        candidate, _ = prepare_candidate(vars, i => n)
        x = CtrlVQE.Parameters.values(candidate)

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

        # REGISTER THE GRADIENT NORM
        scores[n,i] = g[end]
    end; end
    return scores
end
score_candidates(args...) = score_candidates(_!, args...)

""" Just inserts one more Ω parameter, starting from zero. """
function adapt_state!(vars, upHk)
    JOB.require_work(vars)

    imap = collect(1:length(_!.state.x))

    LΩ = CtrlVQE.Parameters.count(vars.work.device)
    #= NOTE: Semantically this should explicitly omit frequency parameters,
        but this script is only currently supported for fixed frequency devices. =#

    x = deepcopy(_!.state.x)
    insert!(x, 1+LΩ, zero(eltype(_!.state.x)))
    insert!(imap, 1+LΩ, 0)

    Hk = upHk(vars.state.Hk, imap)

    n = _!.state.n  # NOT USED

    Ω = deepcopy(_!.state.Ω); foreach(Ωi -> push!(Ωi, LΩ), Ω)
    ν = 1 .+ _!.state.ν

    vars.state = JOB.StateVars(x, Hk, n, Ω, ν)
end
adapt_state!(args...) = adapt_state!(_!, args...)

function adapt_work!(vars)
    JOB.require_work(vars)
    device = vars.work.device
    x = vars.state.x

    # LOAD THE CURRENT VECTOR SPACE
    space = NPZ.npzread("$(vars.outdir)/vectorspace.npy")    # [ a, n ← σ, i ]

    # TRUNCATE THE SPACE TO MATCH THE SIZE OF THE CURRENT PARAMETER VECTOR
    space = space[1:length(x),:,:]

    # SYNC UP ALL THE MODES
    empty!(device)
    append!(device, space)

    # BIND THE CURRENT PARAMETERS
    CtrlVQE.Parameters.bind(device, x)
end
adapt_work!(args...) = adapt_work!(_!, args...)

##########################################################################################
#= Some script-specific settings. =#

evolution = CtrlVQE.TOGGLE
basis = CtrlVQE.DRESSED
frame = CtrlVQE.STATIC
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_harmonics
upHk = JOB.upHk_slate

##########################################################################################
#= Prepare all the work variables. =#


pool = makepool(_!)

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

device = JOB.ModalHarmonicTransmonDevice(protodevice, length(pool) >> 1, setup.T)


ψ0 = CtrlVQE.QubitOperators.project(system.ψ0, device)
O0 = CtrlVQE.QubitOperators.project(system.model.H, device)
normfn = CtrlVQE.Normalization(evolution, device, basis, grid, ψ0)

lossfn = CtrlVQE.ConstrainedEnergyFunction(
    CtrlVQE.ProjectedEnergy(evolution, device, basis, frame, grid, ψ0, O0),
    JOB.ModalHarmonicAmplitudeBound(device, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
    # NOTE: This script is only compatible with fixed frequency devices.
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
    pool,
)

state = _!.state = isnothing(_!.state) ? (
    # Initialize `vectorspace` file if we are starting from scratch.
    NPZ.npzwrite("$(_!.outdir)/vectorspace.npy", Array{JOB.Float}(
        undef, 0, length(pool), CtrlVQE.ndrives(device))
    );
     JOB.initial_state()
) : _!.state
trace = _!.trace = isnothing(_!.trace) ? JOB.initial_trace() : _!.trace

adapt_work!()
lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)

!_!.run && error("Setup Finished")

##########################################################################################
#= Additional setup. =#

import Optim, LineSearches

function do_optimization()
    # ENSURE WORK OBJECTS ARE CONSISTENT WITH PARAMETERS
    adapt_work!()

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
loaded_converged && adapt_work!()

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
    scores = score_candidates(ϕ)

    # UPDATE THE `vectorspace` ARRAY
    space = NPZ.npzread("$(_!.outdir)/vectorspace.npy")
    newscores = reshape(
        scores ./ maximum(abs.(scores)),    # NORMALIZE SO PARAMETER UNITS ARE CONSISTENT
        1, size(scores)...,                 # RESIZE FOR CONCATENATION INTO VECTOR SPACE
    )
    NPZ.npzwrite("$(_!.outdir)/vectorspace.npy", vcat(space, newscores))

    # UPDATE TRACE
    if loaded_converged
        global loaded_converged = false
    else
        JOB.adapt_trace!(_!, abs.(scores))  # NOTE: adapt_trace! expects normalized scores
    end
    JOB.adapt_is_terminated() && break
        # NOTE: If pool is empty, update sets G_max to 0.0, which always terminates.

    # PERFORM ADAPTATION
    adapt_state!(upHk)
    adapt_work!()

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
