#= Sandbox a version of adaptive windows where real and imaginary components
    are added separately.

Particularly suited to "nodes" pool,
    though the pool will need to be organized a bit differently.

I guess the best way to do this is to add more _drives_, right?
    An X drive and a Y drive?
The drives will be WindowedSignals of Constrained ComplexConstants,
    half fixing B and half fixing A.
We'll need to make a protodevice from systematic as usual,
    but then copy it ala how we did ancillae?

Finally, we'll need to manually drop the α nodes from Y drives,
    and the β nodes from X drives.

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
protopulse = CtrlVQE.ComplexConstant(zero(Float), zero(Float)) #= CtrlVQE.Constant(zero(Float)) =#
devicetype = CtrlVQE.FixedFrequencyTransmonDevice

makepool = JOB.makepool_nodes
select = JOB.select_oneperwindow
upHk = JOB.upHk_slate

##########################################################################################
#= Prepare all the work variables. =#

system = JOB.System(setup.code)
grid = CtrlVQE.TemporalLattice(setup.T, setup.r)

protodevice = CtrlVQE.Systematic(
    devicetype,
    system.n,
    CtrlVQE.UniformWindowed(protopulse, setup.T, 1),
    m = setup.m,
)

pulses = [
    CtrlVQE.UniformWindowed(
        CtrlVQE.ConstrainedSignal(
            protopulse,
            i <= system.n ? :B : :A,
        ),
        setup.T, 1,
    ) for i in 1:2system.n
]

""" Specialized adapt_work! function taking the ConstrainedSignal layer into account. """
function adapt_work!(vars)
    JOB.require_work(vars)

    # MUTATE THE DEVICE TO ACCOUNT FOR NEW PULSE SHAPES
    for i in 1:CtrlVQE.ndrives(vars.work.device)
        signal = CtrlVQE.WindowedSignal([
            CtrlVQE.ConstrainedSignal(
                deepcopy(protopulse),
                i <= vars.work.system.n ? :B : :A,
            ) for _ in vars.state.s[i]],
            vars.state.s[i],
        )
        CtrlVQE.set_drivesignal(vars.work.device, i, signal)
    end
    CtrlVQE.Parameters.bind(vars.work.device, vars.state.x)
end

import LinearAlgebra
""" Specialized initial_state function taking the ConstrainedSignal layer into account. """
function initial_state(vars)
    JOB.require_work(vars)
    # ALIAS FOR CONVENIENCE
    work = vars.work

    # PREPARE SIGNAL INDEXING VECTORS
    p = CtrlVQE.Parameters.count(work.protopulse) ÷ 2
    n = CtrlVQE.ndrives(work.device)
    Ω = [collect(1+(i-1)*p:i*p) for i in 1:n]
    L = sum(length, Ω)

    # PREPARE FREQUENCY INDEXING VECTORS -- Implementation may vary with type someday?
    ν = collect(L+1:L+n)
    work.device isa CtrlVQE.FixedFrequencyTransmonDevice && empty!(ν)
    L += length(ν)

    # PREPARE INITIAL PARAMETER VECTOR
    x = zeros(Float, L)
    for i in eachindex(ν)
        q = CtrlVQE.drivequbit(work.device, i)
        x[ν[i]] = CtrlVQE.resonancefrequency(work.device, q)
    end

    # PREPARE INITIAL INVERSE HESSIAN APPROXIMATION
    Hk = LinearAlgebra.diagm(ones(Float, L))

    # PREPARE INITIAL STARTTIMES
    s = fill([0.0], n)

    return JOB.StateVars(x, Hk, s, Ω, ν)
end



device = devicetype(
    protodevice.ω̄,
    protodevice.δ̄,
    protodevice.ḡ,
    protodevice.quples,
    vcat(protodevice.q̄, protodevice.q̄),
    vcat(protodevice.ν̄, protodevice.ν̄),
    pulses,
    protodevice.m,
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

state = _!.state = isnothing(_!.state) ? initial_state(_!) : _!.state # NOTE: use specialized initial_state
trace = _!.trace = isnothing(_!.trace) ? JOB.initial_trace() : _!.trace

CtrlVQE.Parameters.bind(device, _!.state.x)
lossfn.f_counter[] = isempty(_!.trace.f_calls) ? 0 : last(_!.trace.f_calls)
lossfn.g_counter[] = isempty(_!.trace.g_calls) ? 0 : last(_!.trace.g_calls)

adapt_work!(_!)         # NOTE: Calling specialized adapt_work!


##########################################################################################
#= Additional setup. =#

import Optim, LineSearches

function do_optimization()
    # ENSURE WORK OBJECTS ARE CONSISTENT WITH PARAMETERS
    adapt_work!(_!)         # NOTE: Calling specialized adapt_work!

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


# ##########################################################################################
# #= TEMP: Set up debugging. =#

# protolossfn = CtrlVQE.ConstrainedEnergyFunction(
#     CtrlVQE.ProjectedEnergy(evolution, protodevice, basis, frame, grid, ψ0, O0),
#     CtrlVQE.GlobalAmplitudeBound(protodevice, grid, setup.ΩMAX, setup.λΩ, setup.σΩ),
#     CtrlVQE.GlobalFrequencyBound(protodevice, grid, setup.ΔMAX, setup.λΔ, setup.σΔ),
# )
# protovars = deepcopy(_!)
# protovars.work.device = protodevice
# protovars.work.lossfn = protolossfn
# protoobj = JOB.make_objective(protovars)

# # Check evolution matches for initial parameters.
# ψ  = CtrlVQE.evolve(evolution, protodevice, basis, grid, ψ0)
# ψ_ = CtrlVQE.evolve(evolution, device, basis, grid, ψ0)
# display([ψ;; ψ_])
# display([abs2.(ψ);; abs2.(ψ_)])

# # Check evolution matches for other parameters.
# x0 = CtrlVQE.Parameters.values(device)
# protox = [0.01, -0.02, -0.01, 0.02]
# x = protox[[1,3,2,4]]
# CtrlVQE.Parameters.bind(protodevice, protox)
# CtrlVQE.Parameters.bind(device, x)
# ψ  = CtrlVQE.evolve(evolution, protodevice, basis, grid, ψ0)
# ψ_ = CtrlVQE.evolve(evolution, device, basis, grid, ψ0)
# display([ψ;; ψ_])
# display([abs2.(ψ);; abs2.(ψ_)])

##########################################################################################
#= RUN ADAPT =#

!_!.run && error("Setup Finished")

loaded_converged = trace_is_optimized()

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

    # NOTE: Special for this sandbox: Zero out the even indices in the first half,
    #       and the odd indices in the second half, for node identification.
    ϕ_ = deepcopy(ϕ)
    ϕ_[:,2:2:2system.n] .= 0        # ZERO OUT ϕβ ON α DRIVES
    ϕ_[:,2system.n+1:2:end] .= 0    # ZERO OUT ϕα ON β DRIVES

    # CONSTRUCT POOL AND CALCULATE SCORES
    pool = makepool(_!, ϕ_)
    scores = JOB.score_candidates(pool, ϕ)

    # UPDATE TRACE
    if loaded_converged
        global loaded_converged = false
    else
        JOB.adapt_trace!(_!, scores)
    end
    JOB.adapt_is_terminated() && break
        # NOTE: If pool is empty, update sets G_max to 0.0, which always terminates.

    # PERFORM ADAPTATION
    segments = select(pool, scores)
    JOB.adapt_state!(segments, upHk)
    adapt_work!(_!)         # NOTE: Calling specialized adapt_work!

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
