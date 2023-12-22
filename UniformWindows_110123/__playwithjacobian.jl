#= Consider an optimized state ψ(θ).

Construct the Jacobian (each column lives in Hilbert space, ∂|ψ⟩ wrt one parameter.
What can we learn from this?
- Immediately, rank identifies the dimension actually explored by our optimization.
    (*Technically it only lower-bounds it, but I think the optimized point
        is almost certainly non-singular? Need to think that through...)
- Identifying which columns are not linearly independent with the preceding set,
    one can immediately identify redundant parameters.

My sense is that doing an SVD should go even further:
- The parameter-wise eigenvectors of the largest singular values
    should identify exactly which parameters are the most "useful".
- The condition (inverse of smallest non-zero singular value)
    should indicate that the parameterization,
    while technically full rank, mayn't be as effective as another.
- The parameter-wise eigenvectors might just identify the constraints (to first-order)
    needed for an "optimal" parameterization.

So, let's load up a two-qubit system, calculate the Jacobian,
    and see if we can't work out an ideal pulse parameterization.


=#

import CtrlVQE

using Random
using FiniteDifferences
using LinearAlgebra

##########################################################################################
#= PARAMETERS =#

jobdirs = (
    #############
    # H215
    #############
    # detuned_success="jobs/H215_detuned_ΔsMAX3.0/T6.0_W2",
    # detuned_failed="jobs/H215_detuned_ΔsMAX3.0/T5.0_W1",
    # detuned_overdone="jobs/H215_detuned_ΔsMAX3.0/T12.0_W4",

    # resonant_success="jobs/H215_resonant_ΔsMAX3.0/T7.0_W2",
    # resonant_transition="jobs/H215_resonant_ΔsMAX3.0/T6.0_W2",
    # resonant_failed="jobs/H215_resonant_ΔsMAX3.0/T5.0_W1",
    # resonant_overdone="jobs/H215_resonant_ΔsMAX3.0/T12.0_W4",

    # smooth_success="jobs/H215_resonant_perns20/T6.0_W120",
    # smooth_failed="jobs/H215_resonant_perns20/T5.0_W100",
    # smooth_overdone="jobs/H215_resonant_perns20/T12.0_W240",

    #############
    # lih30
    #############
    R   = "jobs/lih30_resonant.real_ΔsMAX1.5/T30.0_W20",
    RF  = "jobs/lih30_detuned.real_ΔsMAX1.5/T30.0_W20",
    RI  = "jobs/lih30_resonant_ΔsMAX3.0/T30.0_W10",
    RIF = "jobs/lih30_detuned_ΔsMAX3.0/T30.0_W10",
    MP  = "jobs/lih30_resonant.polar_ΔsMAX3.0/T30.0_W10",
    MPF = "jobs/lih30_detuned.polar_ΔsMAX3.0/T30.0_W10",
)

##########################################################################################
#= SETUP =#

cfd = central_fdm(5, 1)

function evolve_function(x)
    device = _!.work.device
    energyfn = _!.work.lossfn.energyfn

    x0 = CtrlVQE.Parameters.values(device)
    CtrlVQE.Parameters.bind(device, x)
    ψ = CtrlVQE.evolve(
        energyfn.evolution,
        energyfn.device,
        energyfn.basis,
        energyfn.grid,
        energyfn.ψ0,
    )
    CtrlVQE.Parameters.bind(device, x0)

    return ψ
end

function calculate_initjacobian!(jobdir)
    println("REF $jobdir...")
    JOB.load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl"); catch end
    return jacobian(cfd, evolve_function, JOB.initial_state().x)[1]
end

function perturb_initjacobian!(jobdir, seed, ΩMAXatNσ)
    println("PTB $seed $ΩMAXatNσ $jobdir...")
    JOB.load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl"); catch end
    x0 = JOB.initial_state().x
    Random.seed!(seed)
    x = x0 .+ randn(length(x0)) .* (_!.setup.ΩMAX / ΩMAXatNσ)
    return jacobian(cfd, evolve_function, x)[1]
end

function calculate_jacobian!(jobdir)
    println("OPT $jobdir...")
    JOB.load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl"); catch end
    return jacobian(cfd, evolve_function, _!.state.x)[1]
end

function calculate_perturbation!(jobdir)
    JOB.load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl"); catch end
    return _!.state.x .- JOB.initial_state().x
end

##########################################################################################
#= CALCULATIONS =#

jacobians = NamedTuple(job => calculate_jacobian!(dir) for (job, dir) in pairs(jobdirs))
numparams = NamedTuple(job => size(J,2) for (job, J) in pairs(jacobians))
perturbs = NamedTuple(job => calculate_perturbation!(dir) for (job, dir) in pairs(jobdirs))
numhilbert = NamedTuple(job => size(J,1) for (job, J) in pairs(jacobians))
ranks = NamedTuple(job => rank(J) for (job, J) in pairs(jacobians))
svds = NamedTuple(job => svd(J) for (job, J) in pairs(jacobians))
conditions = NamedTuple(job => cond(J) for (job, J) in pairs(jacobians))
weights = NamedTuple(job => diag(J'*J) for (job, J) in pairs(jacobians))
∂ψs = NamedTuple(job => J[1:2:end,:] .+ im.*J[2:2:end,:] for (job, J) in pairs(jacobians))
n∂ψs = NamedTuple(job => hcat([∂ψ[:,i]./norm(∂ψ[:,i]) for i in axes(∂ψ,2)]...) for (job, ∂ψ) in pairs(∂ψs))
volumes = NamedTuple(job => sum(USV.S) for (job, USV) in pairs(svds))

initjacobians = NamedTuple(job => calculate_initjacobian!(dir) for (job, dir) in pairs(jobdirs))
initranks = NamedTuple(job => rank(J) for (job, J) in pairs(initjacobians))
initsvds = NamedTuple(job => svd(J) for (job, J) in pairs(initjacobians))
initconditions = NamedTuple(job => cond(J) for (job, J) in pairs(initjacobians))
initweights = NamedTuple(job => diag(J'*J) for (job, J) in pairs(initjacobians))
init∂ψs = NamedTuple(job => J[1:2:end,:] .+ im.*J[2:2:end,:] for (job, J) in pairs(initjacobians))
initvolumes = NamedTuple(job => sum(USV.S) for (job, USV) in pairs(initsvds))
initn∂ψs = NamedTuple(job => hcat([∂ψ[:,i]./norm(∂ψ[:,i]) for i in axes(∂ψ,2)]...) for (job, ∂ψ) in pairs(init∂ψs))

##########################################################################################
#= ANALYSIS =#

# Plug columns of V into a pulse, and then plot Ω(t) in both cartesian and polar coordinates.

""" Plot ith eigenpulse of the given svd. Assumes the job is already loaded. """
function plot_eigenpulses!(label, svd, i)
    # Load the eigenparameters.
    _!.state.x = svd.V[:,i]

    # We only care about the shape, so normalize the amplitudes to just saturate bounds.
    Ω, α, β = JOB.Calculations.make_pulses()
    Ωmax = maximum(abs.(Ω))
    _!.state.x[_!.state.Ω] .*= _!.setup.ΩMAX / Ωmax

    # Make the plots.
    prefix = "$(label)/$(i).δ$(round(svd.S[i], digits=2))"
    JOB.plot_pulse(prefix, :moduli, :phases, :amplitudes)
end

""" Plot the eigenpulses for each singular value, at initial and optimized points. """
function plot_eigenpulses!(job)
    JOB.load!(jobdirs[job]; run=false)
    try include("$(_!.outdir)/script.jl"); catch end

    isdir("$(_!.outdir)/jac_opt") || mkdir("$(_!.outdir)/jac_opt")
    isdir("$(_!.outdir)/jac_ref") || mkdir("$(_!.outdir)/jac_ref")

    for i in eachindex(svds[job].S)
        plot_eigenpulses!("jac_opt", svds[job], i)
        plot_eigenpulses!("jac_ref", initsvds[job], i)
    end
end

#= TODO:
What can we argue from the Jacobian?

This is an interesting analysis:
1. Collapse J back into Hilbert space, each column ∂ψ is J[1:2:end] + im J[2:2:end]
2. Normalize each column ∂ψ
3. Calculate revised J' * J, to give ⟨∂iψ|∂jψ⟩
4. Calculate |⋅|² to give overlap |⟨∂iψ|∂jψ⟩|²
5. Subtract off the diagonal, which has no information.
6. Norm of resulting matrix is a measure of parameter redundancy:
    perturbations of different parameters push in the same direction

Meanwhile, omitting the normalization perhaps gives a sense of "volume",
    at least the volume accessible locally from the point where the Jacobian is taken.
Naw this doesn't really work. ∂[ν]ψREF = 0 always, so it doesn't tell anything.
And generally ∂[ν]ψ is actually quite large (which makes sense; it affects full duration)
    so this isn't what you are after.
Naw the most straightforward explanation is that ∂ν lets you push in a direction
    which you can't get from ∂R, but which you can from ∂I.
So:
1. Consider each normalized ∂ν.
2. Subtract projections from non-∂ν.
3. What is norm of the remainder?
    I expect to find norm much larger for RF than for RIF or MPF.
Ngh, have to renormalize at each step. So I guess, keep track of the norm at each step.
Noooo, this still isn't going to work. We KNOW there won't be anything left, becasue (eventually) the resonant real pulse DOES optimize. The DIMENSION isn't the issue. It's the VOLUME. How do we measure the VOLUME? We really need to INTEGRATE over something, like the path trajectory. That is too much!
    And yet, surely, a redundant parameter which is more orthogonal than the others would surely increase the volume? Cross products. Cross products. Yes, that's it. Cross products. A⨯B is small if A and B point in similar directions.
Gar can we just ask what is the total projection into each basis from all real parameters, all complex parameters, and all frequency parameters?



=#