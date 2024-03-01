#=

We are interested in the following values as a function of d, from the "scan_LiH_$d" jobs:
1. HF
2. FCI
3. COR = HF - FCI
4. 1 - |⟨HF|FCI⟩|²
5. T_E: the time at which energy first gets below 1e-6
6. T_C: the time at which iterations are maxed

=#

ε0 = 1e-8       # ENERGY ERROR ON WHICH TO DEFINE eMET

struct Datum
    REF::JOB.Float      # Reference energy.
    FCI::JOB.Float      # Exactly diagonalized energy.
    COR::JOB.Float      # Correlation energy (just REF-FCI).
    iF::JOB.Float       # Infidelity: 1 - |⟨HF|FCI⟩|²
    TE::JOB.Float       # Effective MET: interpolated T when ε intersects with ε0
end

"""
    Datum(d)

Read in all the jobs corresponding to LiH at a particular bond distance.
"""
function Datum(d)
    # THE EASY PART: Load the system!
    code = "scan_LiH_$d"
    system = JOB.System(code)
    REF = system.REF
    FCI = system.FCI
    COR = REF - FCI
    iF = 1 - abs2( system.ψ0' * system.model.U[:,1] )

    # THE HARD PART: Load each job in the code directory, to identify TE and TC.

    T1 = 0      # Highest pulse duration for which ε > ε0.
    T2 = Inf    # Lowest pulse duration for which ε ≤ ε0.

    ε1 = 0.0    # Energy error at T1.
    ε2 = 0.0    # Energy error at T2.

    surveydir = "jobs/$(code)_resonant_ΔsMAX3.0"
    for jobdir in readdir(surveydir, join=true)
        # LOAD VARIABLES
        vars = JOB.load(jobdir)
        isnothing(vars.trace) && continue   # Skip if the job hasn't been run.

        # EXTRACT QUANTITIES OF INTEREST
        T = vars.setup.T
        ε = last(vars.trace.energy) - FCI

        # UPDATE eMET, IF APPROPRIATE
        (ε > ε0) && (T > T1) && (T1 = T; ε1 = ε)
        (ε ≤ ε0) && (T < T2) && (T2 = T; ε2 = ε)
    end

    # THE TRICKY PART: Interpolate T0 from Tl and Tr
    #= NOTE:
        In a couple of runs, there is an anomalous early success,
            meaning Tl (the longest failure) is to the right of Tr (the shortest success).
        But we'll just take our eMET as the interpolated value between the two either way.
    =#
    #= NOTE:
        I want a linear extrapolation on the semi-log plot, so we'll take log10 of each ε.
    =#
    m = (log10(ε2) - log10(ε1)) / (T2 - T1)
    T0 = T1 + (log10(ε0) - log10(ε1)) / m

    return Datum(REF, FCI, COR, iF, T0)
end

ds = 1.0:0.25:4.75

# LOAD DATA
data = [Datum(d) for d in ds]
REF = [datum.REF for datum in data]
FCI = [datum.FCI for datum in data]
COR = [datum.COR for datum in data]
iF = [datum.iF for datum in data]
TE = [datum.TE for datum in data]



# Try to plot COR, iF, and TE on the same plot. :/
import Plots


plt = Plots.plot(; xlabel="LiH Bond Distance (Å)", ylabel="Effective MET (ns)")
twin = Plots.twinx(plt)
Plots.plot!(plt; xlim=[0.75, 5.0])
Plots.plot!(plt, ds, TE; color=1, lw=2, label="Effective MET", shape=:circle)
Plots.plot!(twin, ds, iF; color=2, lw=2, ylabel="1-|⟨REF|FCI⟩|²", label=false)
Plots.plot!(plt, [0.0], [minimum(TE)]; color=2, lw=2, label="1-|⟨REF|FCI⟩|²")
Plots.savefig(plt, "figs/TEandIF.pdf")

plt = Plots.plot(; xlabel="LiH Bond Distance (Å)", ylabel="Effective MET (ns)")
twin = Plots.twinx(plt)
Plots.plot!(plt; xlim=[0.75, 5.0])
Plots.plot!(plt, ds, TE; color=1, lw=2, label="Effective MET", shape=:circle)
Plots.plot!(twin, ds, COR; color=2, lw=2, ylabel="Correlation Energy (Ha)", label=false)
Plots.plot!(plt, [0.0], [minimum(TE)]; color=2, lw=2, label="Correlation Energy")
Plots.savefig(plt, "figs/TEandCOR.pdf")

plt = Plots.plot(; xlabel="Correlation Energy (Ha)", ylabel="Effective MET (ns)", xlim=[0.0,maximum(COR)])
Plots.plot!(plt, COR, TE; color=1, lw=2, label=false, shape=:circle)
Plots.savefig(plt, "figs/TEvsCOR.pdf")

plt = Plots.plot(;
    xlabel="1-|⟨REF|FCI⟩|²",
    ylabel="Effective MET (ns)",
    xlim=[0.0,maximum(iF)],
)
Plots.plot!(plt, iF, TE; color=1, lw=2, label=false, shape=:circle)
Plots.savefig(plt, "figs/TEvsIF.pdf")