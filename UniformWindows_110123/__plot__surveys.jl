#= Plot data from many choices of `T` as a single curve. =#

import Plots

default_series = (
    label = false,
    color = :black,
    linestyle = :solid,
    linewidth = 2,
    shape = :circle,
    marksersize = 5,
)

""" kwargs for plot initialization which sets the T-axis from TMIN to TMAX,
        with minor gridlines every Δs ns.
    Major ticks are selected to have approximately no more than ten labels total.
"""
function __initplot__T_axis(TMIN, TMAX, Δs)
    labelskip = 1+ round(Int, (TMAX - TMIN) / Δs) ÷ 10
    return (
        xlims = [TMIN,TMAX],
        xticks = TMIN:labelskip*Δs:TMAX,
        xminorticks = labelskip,
        xminorgrid = true,
    )
end

__initplot__energy(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Energy Error",
    yscale = :log,
    ylims  = [1e-16, 1e2],
    yticks = 10.0 .^ (-16:2:0),
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    kwargs...
)

__initplot__iterations(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Iteration Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    kwargs...
)

function init_plots(TMIN, TMAX, Δs; kwargs...)
    T_axis = __initplot__T_axis(TMIN, TMAX, Δs)
    return (
        energy = __initplot__energy(; T_axis..., kwargs...),
        iterations = __initplot__iterations(; T_axis..., kwargs...),
    )
end

function plot_survey!(plots, surveydir; kwargs...)
    system = Ref{JOB.System}()

    # INITIALIZE RESULT VECTORS
    data = (
        T = JOB.Float[],
        iterations = Int[],
        energy = JOB.Float[],
    )

    for jobdir in readdir(surveydir, join=true)
        # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
        vars = JOB.load(jobdir)
        isnothing(vars.trace) && continue
        isassigned(system) || (system[] = JOB.System(vars.setup.code))

        # REGISTER THE DATA
        push!(data.T, vars.setup.T)
        push!(data.iterations, last(vars.trace.iterations))
        push!(data.energy, last(vars.trace.energy) - system[].FCI)  # Energy ERROR, really
    end

    # SORT OUTPUTS BY TIME
    σ = sortperm(data.T)
    for array in data; permute!(array, σ); end

    # PLOT DATA
    Plots.plot!(plots.energy, data.T, data.energy; default_series..., kwargs...)
    Plots.plot!(plots.iterations, data.T, data.iterations; default_series..., kwargs...)
end

function save_plots(plots, label)
    for (kind, plot) in pairs(plots)
        Plots.savefig(plot, "figs/$(label)__$(kind).pdf")
    end
end

# # HYDROGEN CASE STUDY
# plots = init_plots(0.0, 30.0, 3.0)
# plot_survey!(plots, "jobs/H2EPm_resonant_ΔsMAX3.0";
#     label="0.75Å", color=0, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/H215_resonant_ΔsMAX3.0";
#     label="1.50Å", color=0, linestyle=:dash, shape=:star)
# plot_survey!(plots, "jobs/H2SPm_resonant_ΔsMAX3.0";
#     label="3.00Å", color=0, linestyle=:solid, shape=:utriangle)
# save_plots(plots, "H2")

# # FULL SURVEY
# plots = init_plots(0.0, 100.0, 3.0)
# plot_survey!(plots, "jobs/H2EPm_resonant_ΔsMAX3.0";
#     label="H₂ 0.75Å", color=0, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/H215_resonant_ΔsMAX3.0";
#     label="H₂ 1.50Å", color=0, linestyle=:dash, shape=:star)
# plot_survey!(plots, "jobs/H2SPm_resonant_ΔsMAX3.0";
#     label="H₂ 3.00Å", color=0, linestyle=:solid, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label="LiH 1.50Å", color=1, linestyle=:dash, shape=:star)
# plot_survey!(plots, "jobs/H4EPm_resonant_ΔsMAX3.0";
#     label="H₄ 0.90Å", color=2, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/H4SPm_resonant_ΔsMAX3.0";
#     label="H₄ 3.00Å", color=2, linestyle=:solid, shape=:utriangle)
# save_plots(plots, "survey")

# # COMPARING DEVICE PARAMETERS
# plots = init_plots(0.0, 48.0, 3.0)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label="Baseline", color=1, linestyle=:solid, shape=:star)
# plot_survey!(plots, "jobs/lih30_resonant_largercoupling_ΔsMAX3.0";
#     label="Doubled Coupling", color=2, linestyle=:solid, shape=:star)
# plot_survey!(plots, "jobs/lih30_resonant_largerspacing_ΔsMAX3.0";
#     label="Doubled Spacing", color=3, linestyle=:solid, shape=:star)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX1.5";
#     label="Doubled Windows", color=4, linestyle=:solid, shape=:star)
# save_plots(plots, "device")


# # HUBBARD LATTICE REFERENCE/BASIS CHOICE
# plots = init_plots(42.0, 66.0, 3.0)#; legend=:bottomleft)
# plot_survey!(plots, "jobs/cL4tAPm_resonant_ΔsMAX3.0";
#     label="Small u, atomic orbitals", color=1, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/cL4tCPm_resonant_ΔsMAX3.0";
#     label="Small u, core orbitals", color=2, linestyle=:dot, shape=:square)
# plot_survey!(plots, "jobs/cL4UAPm_resonant_ΔsMAX3.0";
#     label="Large u, atomic orbitals", color=3, linestyle=:solid, shape=:circle)
# plot_survey!(plots, "jobs/cL4UCPm_resonant_ΔsMAX3.0";
#     label="Large u, core orbitals", color=4, linestyle=:dot, shape=:circle)
# save_plots(plots, "basischoice")

# # COMPARING `ΔsMAX` vs `perns` PARAMETERIZATIONS
# plots = init_plots(0.0, 48.0, 3.0)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label="Δs ≥ 3.0ns", color=1, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX1.5";
#     label="Δs ≥ 1.5ns", color=2, linestyle=:solid, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_resonant_perns20";
#     label="Continuous", color=:black, linestyle=:solid, shape=:circle)
# save_plots(plots, "windowspacing")

# COMPARING REAL vs COMPLEX, RESONANT vs DETUNED PARAMETERIZATIONS
plots = init_plots(0.0, 48.0, 3.0)
plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
    label="RI", color=1, linestyle=:solid, shape=:square)
plot_survey!(plots, "jobs/lih30_detuned_ΔsMAX3.0";
    label="RIF", color=2, linestyle=:solid, shape=:circle)
plot_survey!(plots, "jobs/lih30_resonant.polar_ΔsMAX3.0";
    label="MP", color=3, linestyle=:dot, shape=:square)
plot_survey!(plots, "jobs/lih30_detuned.polar_ΔsMAX3.0";
    label="MPF", color=4, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/lih30_resonant.real_ΔsMAX3.0";
#     label="A", color=3, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned.real_ΔsMAX3.0";
#     label="AF", color=4, linestyle=:dot, shape=:circle)
plot_survey!(plots, "jobs/lih30_resonant.real_ΔsMAX1.5";
    label="R", color=5, linestyle=:dash, shape=:square)
plot_survey!(plots, "jobs/lih30_detuned.real_ΔsMAX1.5";
    label="RF", color=6, linestyle=:dash, shape=:circle)
save_plots(plots, "realresonance")