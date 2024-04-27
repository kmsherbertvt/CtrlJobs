#= Run the current pulse (in `_!.state`) with a callback that records ψ(t) at each timestep,
    and then save it with a meaningful name. =#

# This script expects to be in a ./juliarepl, using the project module aliased as `JOB`.

import CtrlVQE
import NPZ
import Plots
import ColorSchemes

"""
Initialize an array to store statevectors at every iteration of a time-evolution,
    and the callback to pass to the evolution method which fills it.

Note that Ψ will be filled in the evolution algorithm's work basis, and in the lab frame.
"""
function statevector_callback(vars)
    r = CtrlVQE.nsteps(vars.work.grid)
    N = CtrlVQE.nstates(vars.work.device)

    Ψ = Matrix{ComplexF64}(undef, N, r+1)
    return Ψ, (i, t, ψ) -> (Ψ[:,i] .= ψ)
end

function rotated_statevector_callback(vars)
    r = CtrlVQE.nsteps(vars.work.grid)
    N = CtrlVQE.nstates(vars.work.device)

    basis = vars.work.lossfn.energyfn.basis
    frame = vars.work.lossfn.energyfn.frame
    workbasis = CtrlVQE.Evolutions.workbasis(vars.work.lossfn.energyfn.evolution)
    U = CtrlVQE.Devices.basisrotation(basis, workbasis, vars.work.device)
    ψ_ = Vector{ComplexF64}(undef, N)

    Ψ = Matrix{ComplexF64}(undef, N, r+1)
    return Ψ, (i, t, ψ) -> (
        ψ_ .= ψ;
        CtrlVQE.LinearAlgebraTools.rotate!(U, ψ_);
        # TODO: If we decide to generalize this, account for projection.
        CtrlVQE.Devices.evolve!(frame, vars.work.device, basis, -t, ψ_);
        Ψ[:,i] .= ψ_
    )
end

"""
Load a job, run its initialization script, and run the pulse
    to generate a matrix of the statevectors for each timestep.

# Keyword Arguments
- archive: pass in a string to unarchive a particular state rather than the one in `state`
- save: If false, doesn't write statevector matrix to a file.
        If true, serializes statevector matrix to a file called "statevectors".
        If a string, serializes statevector matrix to a file with that name.
"""
function collect_statevector!(jobdir; archive=nothing, save=true)
    save === true && (save = "statevectors")
    filename = "$jobdir/$save.npy"
    save isa String && isfile(filename) && return NPZ.npzread(filename)

    # SETUP
    load!(jobdir; run=false)                            # LOAD THE FILE
    try include("$jobdir/script.jl") catch end          # INITIALIZE WORK VARIABLES
    isnothing(archive) || unarchive!(archive)           # LOAD A PARTICULAR STATE
    Ψ, callback = rotated_statevector_callback(_!)              # PREP THE STATEVECTOR MATRIX

    # RUN THE TARGET PULSE
    CtrlVQE.cost_function(_!.work.lossfn.energyfn; callback=callback)(_!.state.x)

    # SAVE IF DESIRED
    save isa String && NPZ.npzwrite(filename, Ψ)
    return Ψ
end


function fidelity_plot!(job1, job2)
    Ψ1 = collect_statevector!(job1)
    Ψ2 = collect_statevector!(job2)

    matrix = abs2.(Ψ2' * Ψ1)
    return Plots.heatmap(matrix;
        ylabel="Ψ2", xlabel="Ψ1", color=:hot,
        aspect_ratio=:equal, #size=(500,500),
        xlims=[0.5, 0.5+size(Ψ1,2)],
        ylims=[0.5, 0.5+size(Ψ2,2)],
    )
end

function eigenbasis_plot!(job)
    Ψ = collect_statevector!(job)

    vars = load(job)
    system = JOB.System(vars.setup.code)
    matrix = abs2.(system.model.U' * Ψ)
    return Plots.heatmap(matrix; ylabel="Ψ", xlabel="t", color=:hot)
end

function weight(z)
    return z == 0 ? 0 : weight(z>>1) + z & 1
end

"""
This is a bit confusing...

We like the eigenbasis plot.
The vanilla ordering puts higher-energy states at the top,
    so bright bands there means we are exploring the high-energy space.
    Very interesting.

But, we can consider a different ordering which is "accessibility via pulse".
So, bright bands at the top would mean the pulse is doing "more work" to make it happen.

Personally I think basis state occupation is more natural in this ranking,
    rather than the eigenbasis occupation.
We'll have to try it both ways...


"""
function perturbasis_plot!(job)
    Ψ = collect_statevector!(job)

    # FETCH THE BITSTRING CORRESPONDING TO OUR REFERENCE STATE
    zREF = argmax(abs2.(Ψ[:,1])) - 1

    # RANK ALL THE BASIS VECTORS ACCORDING TO HAMMING DISTANCE FROM REFERENCE
    r = [weight(z ⊻ zREF) for z in 0:(size(Ψ,1)-1)]
    σ = sortperm(r)
    rσ = r[σ]                                   # SORTED WEIGHTS
    increments = findall(Bool.(diff(rσ)))       # INDICES WHERE WEIGHT CHANGES

    matrix = abs2.(Ψ[σ,:])
    return Plots.heatmap(
        matrix;
        yticks=(increments, rσ[increments]),
        xlabel="Time Index",
    )
end


# GOOD RESULTS FROM SHOTGUN SURVEY
#= Best run for Δs ≥ 3.0ns was 44. Roughly similar pulses: 66, 73, 57, 62...
    91 and 49 are also pretty similar, but distinct enough to check.
        (but it turns out, yeah, they're pretty darned similar...).

    So it looks like 44 and 90 are unique paths. =#
T20_Δs3_44 = "jobs/lih30_resonant_ΔsMAX3.0/T20.0_W6/shotgun_44"
T20_Δs3_91 = "jobs/lih30_resonant_ΔsMAX3.0/T20.0_W6/shotgun_91"     # Sim. path to 44.
T20_Δs3_90 = "jobs/lih30_resonant_ΔsMAX3.0/T20.0_W6/shotgun_90"
T20_Δs3_49 = "jobs/lih30_resonant_ΔsMAX3.0/T20.0_W6/shotgun_49"     # Sim. path to 44.

#= All the top ten runs for Δs ≥ 1.5ns look different.
    Correlations do reveal similar paths, thank goodness.

    So it looks like 14, 31, 49, 39, and 19 are unique paths. =#
T18_Δs15_14 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_14"
T18_Δs15_31 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_31"
T18_Δs15_2 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_2"     # Sim. path to 14.
T18_Δs15_36 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_36"   # Sim. path to 14.
T18_Δs15_49 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_49"
T18_Δs15_39 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_39"
T18_Δs15_6 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_6"     # Sim. path to 49.
T18_Δs15_19 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_19"
T18_Δs15_40 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_40"   # Sim. path to 14.
T18_Δs15_4 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_4"     # Sim. path to 39.
T18_Δs15_15 = "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12/shotgun_15"   # Sim. path to 39.

#= Successful runs for continuous trials. [0], 3, 5, 17
    The moduli are all close to saturated, but phases are a little different.

    Correlations reveal from-zero and seed 3 are essentially identical,
        and 5 is very similar but not totally so.
    17 pulse looks identical to 5, so no need to check it.
    =#
T17_Δs0 = "jobs/lih30_resonant.limited_perns20/T17.0_W340"
T17_Δs0_3 = "jobs/lih30_resonant.limited_perns20/T17.0_W340/shotgun_3"
T17_Δs0_5 = "jobs/lih30_resonant.limited_perns20/T17.0_W340/shotgun_5"

#= Just for fun, take the best run from Δs≥3.0ns, at T=19.0ns.
    This run is one I consider a *failure*,
        but let's check if the path resembles successes at shorter times. =#
T19_Δs3_26 = "jobs/lih30_resonant_ΔsMAX3.0/T19.0_W6/shotgun_26"



pert = "figs/pathperturbasis"
eign = "figs/patheigenbasis"
corr = "figs/pathcorrelation"


# ##########################################################################################
# #= PLOTTING TRAJECTORIES, FOR EACH PARTICULARLY UNIQUE PLOT =#


# Plots.savefig(fidelity_plot!(T20_Δs3_44, T20_Δs3_44), "$corr/auto.T20_Δs3_44.pdf")
# Plots.savefig(eigenbasis_plot!(T20_Δs3_44), "$eign/T20_Δs3_44.pdf")
# Plots.savefig(perturbasis_plot!(T20_Δs3_44), "$pert/T20_Δs3_44.pdf")

# Plots.savefig(fidelity_plot!(T20_Δs3_90, T20_Δs3_90), "$corr/auto.T20_Δs3_90.pdf")
# Plots.savefig(eigenbasis_plot!(T20_Δs3_90), "$eign/T20_Δs3_90.pdf")
# Plots.savefig(perturbasis_plot!(T20_Δs3_90), "$pert/T20_Δs3_90.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_14), "$corr/auto.T18_Δs15_14.pdf")
# Plots.savefig(eigenbasis_plot!(T18_Δs15_14), "$eign/T18_Δs15_14.pdf")
# Plots.savefig(perturbasis_plot!(T18_Δs15_14), "$pert/T18_Δs15_14.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_31), "$corr/auto.T18_Δs15_31.pdf")
# Plots.savefig(eigenbasis_plot!(T18_Δs15_31), "$eign/T18_Δs15_31.pdf")
# Plots.savefig(perturbasis_plot!(T18_Δs15_31), "$pert/T18_Δs15_31.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_49), "$corr/auto.T18_Δs15_49.pdf")
# Plots.savefig(eigenbasis_plot!(T18_Δs15_49), "$eign/T18_Δs15_49.pdf")
# Plots.savefig(perturbasis_plot!(T18_Δs15_49), "$pert/T18_Δs15_49.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_39, T18_Δs15_39), "$corr/auto.T18_Δs15_39.pdf")
# Plots.savefig(eigenbasis_plot!(T18_Δs15_39), "$eign/T18_Δs15_39.pdf")
# Plots.savefig(perturbasis_plot!(T18_Δs15_39), "$pert/T18_Δs15_39.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_19, T18_Δs15_19), "$corr/auto.T18_Δs15_19.pdf")
# Plots.savefig(eigenbasis_plot!(T18_Δs15_19), "$eign/T18_Δs15_19.pdf")
# Plots.savefig(perturbasis_plot!(T18_Δs15_19), "$pert/T18_Δs15_19.pdf")

# Plots.savefig(fidelity_plot!(T17_Δs0, T17_Δs0), "$corr/auto.T17_Δs0.pdf")
# Plots.savefig(eigenbasis_plot!(T17_Δs0), "$eign/T17_Δs0.pdf")
# Plots.savefig(perturbasis_plot!(T17_Δs0), "$pert/T17_Δs0.pdf")

# Plots.savefig(fidelity_plot!(T17_Δs0_5, T17_Δs0_5), "$corr/auto.T17_Δs0_5.pdf")
# Plots.savefig(eigenbasis_plot!(T17_Δs0_5), "$eign/T17_Δs0_5.pdf")
# Plots.savefig(perturbasis_plot!(T17_Δs0_5), "$pert/T17_Δs0_5.pdf")

# Plots.savefig(fidelity_plot!(T19_Δs3_26, T19_Δs3_26), "$corr/auto.T19_Δs3_26.pdf")
# Plots.savefig(eigenbasis_plot!(T19_Δs3_26), "$eign/T19_Δs3_26.pdf")
# Plots.savefig(perturbasis_plot!(T19_Δs3_26), "$pert/T19_Δs3_26.pdf")

# ##########################################################################################
# #= CORRELATING TRAJECTORIES =#

# Plots.savefig(fidelity_plot!(T20_Δs3_44, T20_Δs3_91), "$corr/T20_Δs3_44.T20_Δs3_91.pdf")
# Plots.savefig(fidelity_plot!(T20_Δs3_44, T20_Δs3_90), "$corr/T20_Δs3_44.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T20_Δs3_44, T20_Δs3_49), "$corr/T20_Δs3_44.T20_Δs3_49.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_31), "$corr/T18_Δs15_14.T18_Δs15_31.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_2), "$corr/T18_Δs15_14.T18_Δs15_2.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_36), "$corr/T18_Δs15_14.T18_Δs15_36.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_49), "$corr/T18_Δs15_14.T18_Δs15_49.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_39), "$corr/T18_Δs15_14.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_6), "$corr/T18_Δs15_14.T18_Δs15_6.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_19), "$corr/T18_Δs15_14.T18_Δs15_19.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_40), "$corr/T18_Δs15_14.T18_Δs15_40.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_4), "$corr/T18_Δs15_14.T18_Δs15_4.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_14, T18_Δs15_15), "$corr/T18_Δs15_14.T18_Δs15_15.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_49), "$corr/T18_Δs15_31.T18_Δs15_49.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_39), "$corr/T18_Δs15_31.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_6), "$corr/T18_Δs15_31.T18_Δs15_6.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_19), "$corr/T18_Δs15_31.T18_Δs15_19.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_4), "$corr/T18_Δs15_31.T18_Δs15_4.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T18_Δs15_15), "$corr/T18_Δs15_31.T18_Δs15_15.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_39), "$corr/T18_Δs15_49.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_6), "$corr/T18_Δs15_49.T18_Δs15_6.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_19), "$corr/T18_Δs15_49.T18_Δs15_19.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_4), "$corr/T18_Δs15_49.T18_Δs15_4.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T18_Δs15_15), "$corr/T18_Δs15_49.T18_Δs15_15.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_39, T18_Δs15_19), "$corr/T18_Δs15_39.T18_Δs15_19.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_39, T18_Δs15_4), "$corr/T18_Δs15_39.T18_Δs15_4.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_39, T18_Δs15_15), "$corr/T18_Δs15_39.T18_Δs15_15.pdf")

# Plots.savefig(fidelity_plot!(T17_Δs0, T17_Δs0_3), "$corr/T17_Δs0.T17_Δs0_3.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T17_Δs0_5), "$corr/T17_Δs0.T17_Δs0_5.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_3, T17_Δs0_5), "$corr/T17_Δs0_3.T17_Δs0_5.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_14, T20_Δs3_44), "$corr/T18_Δs15_14.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T20_Δs3_44), "$corr/T18_Δs15_31.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T20_Δs3_44), "$corr/T18_Δs15_49.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_39, T20_Δs3_44), "$corr/T18_Δs15_39.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_19, T20_Δs3_44), "$corr/T18_Δs15_19.T20_Δs3_44.pdf")

# Plots.savefig(fidelity_plot!(T18_Δs15_14, T20_Δs3_90), "$corr/T18_Δs15_14.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_31, T20_Δs3_90), "$corr/T18_Δs15_31.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_49, T20_Δs3_90), "$corr/T18_Δs15_49.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_39, T20_Δs3_90), "$corr/T18_Δs15_39.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T18_Δs15_19, T20_Δs3_90), "$corr/T18_Δs15_19.T20_Δs3_90.pdf")

# Plots.savefig(fidelity_plot!(T17_Δs0, T20_Δs3_44), "$corr/T17_Δs0.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T20_Δs3_90), "$corr/T17_Δs0.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T18_Δs15_14), "$corr/T17_Δs0.T18_Δs15_14.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T18_Δs15_31), "$corr/T17_Δs0.T18_Δs15_31.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T18_Δs15_49), "$corr/T17_Δs0.T18_Δs15_49.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T18_Δs15_39), "$corr/T17_Δs0.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0, T18_Δs15_19), "$corr/T17_Δs0.T18_Δs15_19.pdf")

# Plots.savefig(fidelity_plot!(T17_Δs0_5, T20_Δs3_44), "$corr/T17_Δs0_5.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T20_Δs3_90), "$corr/T17_Δs0_5.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T18_Δs15_14), "$corr/T17_Δs0_5.T18_Δs15_14.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T18_Δs15_31), "$corr/T17_Δs0_5.T18_Δs15_31.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T18_Δs15_49), "$corr/T17_Δs0_5.T18_Δs15_49.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T18_Δs15_39), "$corr/T17_Δs0_5.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T17_Δs0_5, T18_Δs15_19), "$corr/T17_Δs0_5.T18_Δs15_19.pdf")

# Plots.savefig(fidelity_plot!(T19_Δs3_26, T17_Δs0), "$corr/T19_Δs3_26.T17_Δs0.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T17_Δs0_5), "$corr/T19_Δs3_26.T17_Δs0_5.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T20_Δs3_44), "$corr/T19_Δs3_26.T20_Δs3_44.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T20_Δs3_90), "$corr/T19_Δs3_26.T20_Δs3_90.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T18_Δs15_14), "$corr/T19_Δs3_26.T18_Δs15_14.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T18_Δs15_31), "$corr/T19_Δs3_26.T18_Δs15_31.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T18_Δs15_49), "$corr/T19_Δs3_26.T18_Δs15_49.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T18_Δs15_39), "$corr/T19_Δs3_26.T18_Δs15_39.pdf")
# Plots.savefig(fidelity_plot!(T19_Δs3_26, T18_Δs15_19), "$corr/T19_Δs3_26.T18_Δs15_19.pdf")


#= Make a specific fidelity plot, with proper axes and labels and stuff. =#
let
    Δt = 1/20

    Ψ1 = collect_statevector!(T18_Δs15_49)
    t1 = 0:Δt:18.0
    tick1 = 0:3:18
    label1 = L"Time (ns): $s = 1.5 \mathrm{ns}$"

    Ψ2 = collect_statevector!(T17_Δs0)
    t2 = 0:Δt:17.0
    tick2 = [0:3:15..., 17]
    label2 = L"Time (ns): $s = 0.05 \mathrm{ns}$"

    matrix = abs2.(Ψ2' * Ψ1)
    plt = Plots.heatmap(t1, t2, matrix;
        color=:hot,
        aspect_ratio=:equal, #size=(500,500),
        xlabel=label1,
        xlims=[minimum(t1), maximum(t1)],
        xticks=tick1,
        xminorticks=3,
        ylabel=label2,
        ylims=[minimum(t2), maximum(t2)],
        yticks=tick2,
        yminorticks=3,
        tickfontsize = 11,  # Match the default label size.
        colorbar_tickfontsize = 11,
    )

    Plots.savefig(plt, "figs/fidelity.pdf")
end