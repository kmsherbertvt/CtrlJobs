#= Plot data from many choices of `T` as a single curve. =#

import Plots
import ColorSchemes
using LaTeXStrings

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
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__iterations(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Iteration Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

function init_plots(TMIN, TMAX, Δs; energy=Dict(), iter=Dict(), kwargs...)
    T_axis = __initplot__T_axis(TMIN, TMAX, Δs)
    return (
        energy = __initplot__energy(; T_axis..., kwargs..., energy...),
        iterations = __initplot__iterations(; T_axis..., kwargs..., iter...),
    )
end

function plot_survey!(plots, surveydir; energy=Dict(), iter=Dict(), kwargs...)
    system = Ref{JOB.System}()

    # INITIALIZE RESULT VECTORS
    data = (
        T = JOB.Float[],
        iterations = Int[],
        energy = JOB.Float[],
    )

    for jobdir in readdir(surveydir, join=true)
        endswith(jobdir, ".DS_Store") && continue
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
    Plots.plot!(plots.energy, data.T, data.energy;
            default_series..., kwargs..., energy...)
    Plots.plot!(plots.iterations, data.T, data.iterations;
            default_series..., kwargs..., iter...)
end

function save_plots(plots, label)
    for (kind, plot) in pairs(plots)
        Plots.savefig(plot, "figs/$(label)__$(kind).pdf")
    end
end

# # SINGLE LITHIUM HYDRIDE EXAMPLE
# plots = init_plots(0.0, 48.0, 3.0)
# #= Hack in some "pencil marks" indicating where the shotgun surveys landed. =#
# E0 = JOB.System("lih30").FCI
# Ts = [12.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 30.0, 48.0]
# Ws = [4, 6, 6, 6, 7, 7, 7, 8, 10, 16]
# xs = [[] for _ in Ts]   #           "       "
# εs = [[] for _ in Ts]   # Vectors to be filled with data.
# Cs = [[] for _ in Ts]   #           "       "
# for i in eachindex(Ts)
#     jobdir = "jobs/lih30_resonant_ΔsMAX3.0/T$(Ts[i])_W$(Ws[i])"
#     for pelletdir in readdir(jobdir)
#         startswith(pelletdir, "shotgun_") || continue
#         vars = load("$jobdir/$pelletdir")
#         push!(xs[i], Ts[i])
#         push!(εs[i], last(vars.trace.energy) - E0)
#         push!(Cs[i], last(vars.trace.iterations))
#     end
# end
# Plots.scatter!(plots.energy, vcat(xs...), vcat(εs...);
#     msw=2, shape=:hline, color=:lightgray, label="Random Initial Ω")
# Plots.scatter!(plots.iterations, vcat(xs...), vcat(Cs...);
#     msw=2, shape=:hline, color=:lightgray, label="Random Initial Ω")
# #= End of hack. =#
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label="Initial Ω=0", color=:black, linestyle=:solid, shape=:circle)
# save_plots(plots, "LiH3.0")


# # HYDROGEN CASE STUDY
# plots = init_plots(0.0, 30.0, 3.0)
# plot_survey!(plots, "jobs/H2EPm_resonant_ΔsMAX3.0";
#     label="0.75Å", color=0, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/H215_resonant_ΔsMAX3.0";
#     label="1.50Å", color=0, linestyle=:dash, shape=:star)
# plot_survey!(plots, "jobs/H2SPm_resonant_ΔsMAX3.0";
#     label="3.00Å", color=0, linestyle=:solid, shape=:utriangle)
# save_plots(plots, "H2")

# # FULL SURVEY UP TO H4
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

# # FULL SURVEY
# xticks = [3.0, 10.0, 20.0, 60.0, 400.0, 1500.0]
# plots = init_plots(3.0, 1536.0, 3.0;
#     xscale=:log10,
#     xticks=(xticks, string.(xticks)),
#     xlims=[3.0,2000.0],
#     xminorgrid=false,
#     # legend=:bottomleft,
#     legendfont=Plots.font(7),
#     iter=Dict(
#         :yscale => :log10,
#         :legend => :bottomright,
#     ),
# )
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
# plot_survey!(plots, "jobs/H6EPm_resonant_ΔsMAX3.0";
#     label="H₆ 0.95Å", color=3, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/H6SPm_resonant_ΔsMAX3.0";
#     label="H₆ 3.00Å", color=3, linestyle=:solid, shape=:utriangle)
# save_plots(plots, "longersurvey")

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
#     label=L"\{αβ\}", color=1, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX1.5";
#     label=L"\{αβ\}_2", color=2, linestyle=:solid, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_resonant_perns20";
#     label=L"\{αβ\}_∞", color=:black, linestyle=:solid, shape=:circle)
# save_plots(plots, "windowspacing")

# # COMPARING REAL vs COMPLEX, RESONANT vs DETUNED PARAMETERIZATIONS
# plots = init_plots(0.0, 48.0, 3.0)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label="RI", color=1, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned_ΔsMAX3.0";
#     label="RIF", color=2, linestyle=:solid, shape=:circle)
# plot_survey!(plots, "jobs/lih30_resonant.polar_ΔsMAX3.0";
#     label="MP", color=3, linestyle=:dot, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned.polar_ΔsMAX3.0";
#     label="MPF", color=4, linestyle=:dot, shape=:circle)
# # plot_survey!(plots, "jobs/lih30_resonant.real_ΔsMAX3.0";
# #     label="A", color=3, linestyle=:solid, shape=:square)
# # plot_survey!(plots, "jobs/lih30_detuned.real_ΔsMAX3.0";
# #     label="AF", color=4, linestyle=:dot, shape=:circle)
# plot_survey!(plots, "jobs/lih30_resonant.real_ΔsMAX1.5";
#     label="R", color=5, linestyle=:dash, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned.real_ΔsMAX1.5";
#     label="RF", color=6, linestyle=:dash, shape=:circle)
# save_plots(plots, "realresonance")

# # COMPARING POLAR vs CARTESIAN PARAMETERIZATIONS
# plots = init_plots(0.0, 48.0, 3.0)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label=L"\{αβ\}", color=1, linestyle=:solid, shape=:circle)
# plot_survey!(plots, "jobs/lih30_resonant.polar_ΔsMAX3.0";
#     label=L"\{Aφ\}", color=2, linestyle=:dash, shape=:utriangle)
# save_plots(plots, "polarcartesian")

# # COMPARING REAL vs COMPLEX, RESONANT vs DETUNED PARAMETERIZATIONS
# plots = init_plots(0.0, 48.0, 3.0)
# plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX3.0";
#     label=L"\{αβ\}", color=1, linestyle=:solid, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned_ΔsMAX3.0";
#     label=L"\{αβΔ\}", color=2, linestyle=:solid, shape=:circle)
# # NOTE: Halve window length so number of parameters is consistent.
# plot_survey!(plots, "jobs/lih30_resonant.real_ΔsMAX1.5";
#     label=L"\{α\}", color=3, linestyle=:dash, shape=:square)
# plot_survey!(plots, "jobs/lih30_detuned.real_ΔsMAX1.5";
#     label=L"\{αΔ\}", color=4, linestyle=:dash, shape=:circle)
# save_plots(plots, "resonantdetuned")

# # COMPARING BOND DISTANCES FOR LiH
# plots = init_plots(0.0, 48.0, 3.0)
# palette = ColorSchemes.roma
# ds = 1.0:0.25:4.75
# for (i, d) in enumerate(ds)
#     label = d % 1.0 == 0.0 ? "$d Å" : false
#     plot_survey!(plots, "jobs/scan_LiH_$(d)_resonant_ΔsMAX3.0";
#         label=label, color=palette[i/length(ds)], linealpha=0.5,
#         linewidth=2, linestyle=:solid, shape=:none)
# end
# save_plots(plots, "dissociation")

# # COMPARING BOND DISTANCES FOR H4
# plots = init_plots(0.0, 84.0, 3.0)
# palette = ColorSchemes.roma
# ds = 0.25:0.25:3.0
# for (i, d) in enumerate(ds)
#     label = d % 1.0 == 0.0 ? "$d Å" : false
#     plot_survey!(plots, "jobs/scan_H4_$(d)_resonant_ΔsMAX3.0";
#         label=label, color=palette[i/length(ds)], linealpha=0.5,
#         linewidth=2, linestyle=:solid, shape=:none)
# end
# save_plots(plots, "scan_H4")

# # COMPARING A FEW STRATEGIES OF INTRODUCING ANCILLAE - H2
# plots = init_plots(0.0, 21.0, 0.5)
# plot_survey!(plots, "jobs/H215_resonant_ΔsMAX0.5";
#     label="1xN", color=1, linestyle=:solid, shape=:circle)
# # plot_survey!(plots, "jobs/H215_ancilla.single_ΔsMAX0.5";
# #     label="1x(N+1)", color=2, linestyle=:dash, shape=:utriangle)
# # plot_survey!(plots, "jobs/H215_ancilla.linear_ΔsMAX0.5";
# #     label="1x(2N)", color=3, linestyle=:dash, shape=:utriangle)
# # plot_survey!(plots, "jobs/H215_ancilla_ΔsMAX0.5";
# #     label="2xN", color=4, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle_ΔsMAX0.5";
#     label="N>o", color=5, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle.projected_ΔsMAX0.5";
#     label="N>o Π", color=6, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle.projected.normalized_ΔsMAX0.5";
#     label="N>o Π/N", color=7, linestyle=:dash, shape=:utriangle)
# save_plots(plots, "ancillae")

# COMPARING A FEW STRATEGIES OF INTRODUCING ANCILLAE - LiH
plots = init_plots(0.0, 36.0, 0.5)
plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX0.5";
    label="1xN", color=1, linestyle=:solid, shape=:circle)
# plot_survey!(plots, "jobs/lih30_ancilla.single_ΔsMAX0.5";
#     label="1x(N+1)", color=2, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_ancilla.linear_ΔsMAX0.5";
#     label="1x(2N)", color=3, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_ancilla_ΔsMAX0.5";
#     label="2xN", color=4, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle_ΔsMAX0.5";
    label="N>o", color=5, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle.projected_ΔsMAX0.5";
    label="N>o Π", color=6, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle.projected.normalized_ΔsMAX0.5";
    label="N>o Π/N", color=7, linestyle=:dash, shape=:utriangle)
save_plots(plots, "ancillae_lih30")