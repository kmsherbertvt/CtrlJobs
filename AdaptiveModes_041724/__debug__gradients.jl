#= Quick dirty script to compare finite difference and analytical gradient
    at some point in some trajectory.

Specifically checking the first point in `gradients.exact`
    which "converges" by Δf rather than ∇f,
    since it feels to me like that's the fault of the linesearch,
    which is often catastrophic when the analytical gradient is wrong.

    =#

import CtrlVQE
import FiniteDifferences; cfd = FiniteDifferences.central_fdm(5,1)

##########################################################################################
#= Scoping issue with `try` means we have to define custom auxiliary functions here. =#

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
#= Load job and the particular trial point. =#

jobdir = "jobs/lih30_gradients.exact/T20.0"
load!(jobdir; run=false)            # Load data into `_!`
try JOB.runlocal() catch end        # Initialize work variables and auxiliary functions.

adaptiter = 12
unarchive!(JOB.adaptid(adaptiter))  # Load a particular state into `_!`.
gradients_adapt_work!()             # Update work variables to the particular state.

##########################################################################################
#= Construct cost functions and grad functions for all the pieces. =#

fn = CtrlVQE.cost_function(_!.work.lossfn)
fn_E = CtrlVQE.cost_function(_!.work.lossfn.energyfn)
fn_Λ1 = CtrlVQE.cost_function(_!.work.lossfn.penaltyfns[1])
fn_Λ2 = CtrlVQE.cost_function(_!.work.lossfn.penaltyfns[2])

gd = CtrlVQE.grad_function(_!.work.lossfn)
gd_E = CtrlVQE.grad_function(_!.work.lossfn.energyfn)
gd_Λ1 = CtrlVQE.grad_function(_!.work.lossfn.penaltyfns[1])
gd_Λ2 = CtrlVQE.grad_function(_!.work.lossfn.penaltyfns[2])

##########################################################################################
#= Calculate gradients. =#

x0 = _!.state.x     # Alias the current state

# Analytical gradients:
g0 = gd(x0)
g0_E = gd_E(x0)
g0_Λ1 = gd_Λ1(x0)
g0_Λ2 = gd_Λ2(x0)

# Finite differences:
gΔ = FiniteDifferences.grad(cfd, fn, x0)[1]
gΔ_E = FiniteDifferences.grad(cfd, fn_E, x0)[1]
gΔ_Λ1 = FiniteDifferences.grad(cfd, fn_Λ1, x0)[1]
gΔ_Λ2 = FiniteDifferences.grad(cfd, fn_Λ2, x0)[1]

##########################################################################################
#= Display results. =#

function contrast(g0, gΔ)
    ε = gΔ .- g0

    parse(f) = rpad("$(round(f; digits=8))", 12)    # NOTE: Not robust!

    println("gΔ            g0               gΔ - g0     ")
    println("-------------------------------------------")
    for i in eachindex(ε)
        println("$(parse(gΔ[i]))  $(parse(g0[i]))  |  $(parse(ε[i]))")
    end
    println("-------------------------------------------")
    rms = sqrt(sum(ε.^2)/length(ε))
    println("RMS Error: $rms")
end

println("Total gradient:")
contrast(g0, gΔ)
println()

println("Energy:")
contrast(g0_E, gΔ_E)
println()

println("Amplitude penalty:")
contrast(g0_Λ1, gΔ_Λ1)
println()

println("Frequency penalty:")
contrast(g0_Λ2, gΔ_Λ2)
println()

##########################################################################################
#= Manually compute the amplitude penalty contribution as a function of time. =#

ΩMAX = _!.setup.ΩMAX
t = CtrlVQE.lattice(_!.work.grid)
Ω, α, β = JOB.make_pulses()
ΩL2 = abs.(Ω)

wall(u) = exp(u - 1/u)
grad(u) = exp(u - 1/u) * (1 + 1/u^2)
fΦ(Ω) = (
    u = (abs(Ω) - ΩMAX) / ΩMAX;
    u ≤ 0 ? zero(u) : wall(u)
)
gΦ(Ω, ∂) = (
    u = (abs(Ω) - ΩMAX) / ΩMAX;
    u ≤ 0 ? zero(u) : grad(u) * real(conj(Ω)*∂) / (abs(Ω)*ΩMAX)
)

display(fΦ.(Ω))

# Evidently, only the third pulse has quite reached over the boundary, so focus on it.
L = JOB.CoupledDevices.ndriveparams(_!.work.device)
Ω3 = Ω[:,3]
signal_3 = CtrlVQE.Devices.drivesignal(_!.work.device, 3)
∂3 = Matrix{Complex{JOB.Float}}(undef, length(t), L)
gΦ3 = Matrix{Complex{JOB.Float}}(undef, length(t), L)
for l in 1:L
    CtrlVQE.Signals.partial(l, signal_3, t; result=@view(∂3[:,l]))
    gΦ3[:,l] .= gΦ.(Ω3, @view(∂3[:,l]))
end

display(gΦ3)