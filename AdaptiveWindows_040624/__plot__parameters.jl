#= Show multiple trajectories side-by-side in a single plot vs parameter count.

Sensible things to plot:
- Energy error.
- Total iterations.

Be sure to allow a vertical line to mark "Hilbert dimension".

=#


import Plots
import ColorSchemes

default_series = (
    label = false,
    color = :black,
    linestyle = :solid,
    linewidth = 2,
    shape = :circle,
    marksersize = 5,
)

""" kwargs for plot initialization which sets the P-axis from PMIN to PMAX,
        with minor gridlines every ΔP units.
    Major ticks are selected to have approximately no more than ten labels total.
"""
function __initplot__P_axis(PMIN, PMAX, ΔP)
    labelskip = 1+ round(Int, (PMAX - PMIN) / ΔP) ÷ 10
    return (
        xlabel = "Parameter Count",
        xlims = [PMIN,PMAX],
        xticks = PMIN:labelskip*ΔP:PMAX,
        xminorticks = labelskip,
        xminorgrid = true,
    )
end

__initplot__energy(; kwargs...) = Plots.plot(;
    ylabel = "Energy Error",
    yscale = :log,
    ylims  = [1e-16, 1e2],
    yticks = 10.0 .^ (-16:2:0),
    yminorticks = 1,
    yminorgrid = true,
    legend = :bottomleft,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__iterations(; kwargs...) = Plots.plot(;
    ylabel = "Iteration Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topleft,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

function init_plots(
    PMIN, PMAX, ΔP;
    P0=[],
    energy=Dict(), iter=Dict(), kwargs...
)
    P_axis = __initplot__P_axis(PMIN, PMAX, ΔP)
    energy = __initplot__energy(; P_axis..., kwargs..., energy...)
    iterations = __initplot__iterations(; P_axis..., kwargs..., iter...)

    # VERTICAL MARKS AT EACH P0
    Plots.vline!(energy, P0; label=false, color=:black, lw=1, ls=:dot)
    Plots.vline!(iterations, P0; label=false, color=:black, lw=1, ls=:dot)

    return (
        energy = energy,
        iterations = iterations,
    )
end

function plot_job!(
    plots, jobdir;
    energy=Dict(), iter=Dict(), kwargs...
)
    # INITIALIZE RESULT VECTORS
    data = (
        parameters = Int[],
        energy = JOB.Float[],
        iterations = Int[],
    )

    vars = JOB.load(jobdir)
    system = JOB.System(vars.setup.code)
    nadapt = length(vars.trace.adaptations)

    for a in 1:nadapt
        unarchive!(vars, JOB.adaptid(a))
        push!(data.parameters, length(vars.state.x))

        i = vars.trace.adaptations[a]
        push!(data.iterations, i)
        ε = vars.trace.energy[i] - system.FCI
        push!(data.energy, ε < eps(JOB.Float) ? eps(JOB.Float) : ε)
    end

    # PLOT DATA
    Plots.plot!(plots.energy, data.parameters, data.energy;
            default_series..., kwargs..., energy...)
    Plots.plot!(plots.iterations, data.parameters, data.iterations;
            default_series..., kwargs..., iter...)
end

function save_plots(plots, label)
    for (kind, plot) in pairs(plots)
        Plots.savefig(plot, "figs/$(label)__param__$(kind).pdf")
    end
end

# CONTRASTING PULSE DURATIONS: OPTIMAL
plots = init_plots(0, 50, 2; P0=[30])
plot_job!(plots, "jobs/lih30_optimal/T18.0";
    label="18 ns", color=4, alpha=0.1, linestyle=:solid, shape=:ltriangle)
plot_job!(plots, "jobs/lih30_optimal/T19.0";
    label="19 ns", color=4, alpha=0.25, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_optimal/T20.0";
    label="20 ns", color=4, alpha=0.4, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_optimal/T21.0";
    label="21 ns", color=4, alpha=0.55, linestyle=:solid, shape=:rtriangle)
plot_job!(plots, "jobs/lih30_optimal/T22.0";
    label="22 ns", color=4, alpha=0.7, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_optimal/T23.0";
    label="23 ns", color=4, alpha=0.85, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_optimal/T24.0";
    label="24 ns", color=4, alpha=1.0, linestyle=:solid, shape=:circle)
save_plots(plots, "duration_optimal")

# CONTRASTING PULSE DURATIONS: ONENODE
plots = init_plots(0, 50, 2; P0=[30])
plot_job!(plots, "jobs/lih30_nodes.one/T18.0";
    label="18 ns", color=5, alpha=0.1, linestyle=:solid, shape=:ltriangle)
plot_job!(plots, "jobs/lih30_nodes.one/T19.0";
    label="19 ns", color=5, alpha=0.25, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_nodes.one/T20.0";
    label="20 ns", color=5, alpha=0.4, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_nodes.one/T21.0";
    label="21 ns", color=5, alpha=0.55, linestyle=:solid, shape=:rtriangle)
plot_job!(plots, "jobs/lih30_nodes.one/T22.0";
    label="22 ns", color=5, alpha=0.7, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_nodes.one/T23.0";
    label="23 ns", color=5, alpha=0.85, linestyle=:solid, shape=:none)
plot_job!(plots, "jobs/lih30_nodes.one/T24.0";
    label="24 ns", color=5, alpha=1.0, linestyle=:solid, shape=:circle)
save_plots(plots, "duration_nodes.one")

# CONTRASTING PULSE DURATIONS
function __sweep__(plots, surveydir, label, color)
    plot_job!(plots, "jobs/$surveydir/T18.0";
        label=false, color=color, alpha=0.4, linestyle=:dot, shape=:ltriangle)
    plot_job!(plots, "jobs/$surveydir/T22.0";
        label=label, color=color, alpha=0.8, linestyle=:solid, shape=:circle)
    plot_job!(plots, "jobs/$surveydir/T36.0";
        label=false, color=color, alpha=0.4, linestyle=:dash, shape=:rtriangle)
end
plots = init_plots(0, 50, 2; P0=[30])
__sweep__(plots, "lih30_bisection", "Bisectal", 1)
__sweep__(plots, "lih30_random", "Random", 2)
__sweep__(plots, "lih30_nodes", "Nodal", 3)
__sweep__(plots, "lih30_optimal", "Optimal", 4)
__sweep__(plots, "lih30_nodes.one", "One Node", 5)
# DUMMY PLOTS FOR LEGEND
for plot in plots
    Plots.plot!(plot, Int[], Int[]; default_series...,
        label="18 ns", color=:black, alpha=0.4, linestyle=:dot, shape=:ltriangle)
    Plots.plot!(plot, [], []; default_series...,
        label="22 ns", color=:black, alpha=0.8, linestyle=:solid, shape=:circle)
    Plots.plot!(plot, [], []; default_series...,
        label="36 ns", color=:black, alpha=0.4, linestyle=:dash, shape=:rtriangle)
end
save_plots(plots, "viability")

# CONTRASTING PARALLELISM: OPTIMAL
function __sweep__(plots, surveydir, label, color, ls)
    plot_job!(plots, "jobs/$surveydir/T18.0";
        label=false, color=color, alpha=0.4, linestyle=ls, shape=:ltriangle)
    plot_job!(plots, "jobs/$surveydir/T22.0";
        label=label, color=color, alpha=0.8, linestyle=ls, shape=:circle)
    plot_job!(plots, "jobs/$surveydir/T36.0";
        label=false, color=color, alpha=0.4, linestyle=ls, shape=:rtriangle)
end
plots = init_plots(0, 50, 2; P0=[30])
__sweep__(plots, "lih30_nodes.one", "Select One", 5, :solid)
__sweep__(plots, "lih30_nodes.oneeach", "Select Each", 5, :dash)
__sweep__(plots, "lih30_nodes.oneall", "Select All", 5, :dot)
__sweep__(plots, "lih30_uniform.each", "Uniform", :black, :dash)
# DUMMY PLOTS FOR LEGEND
for plot in plots
    Plots.plot!(plot, Int[], Int[]; default_series...,
        label="18 ns", color=:black, alpha=0.4, linestyle=:solid, shape=:ltriangle)
    Plots.plot!(plot, [], []; default_series...,
        label="22 ns", color=:black, alpha=0.8, linestyle=:solid, shape=:circle)
    Plots.plot!(plot, [], []; default_series...,
        label="36 ns", color=:black, alpha=0.4, linestyle=:solid, shape=:rtriangle)
end
save_plots(plots, "parallel")

# CONTRASTING THE MOST SENSIBLE STRATEGIES
function __sweep__(plots, surveydir, label, color)
    plot_job!(plots, "jobs/$surveydir/T18.0";
        label=false, color=color, alpha=0.4, linestyle=:dot, shape=:ltriangle)
    plot_job!(plots, "jobs/$surveydir/T22.0";
        label=label, color=color, alpha=0.8, linestyle=:solid, shape=:circle)
    plot_job!(plots, "jobs/$surveydir/T36.0";
        label=false, color=color, alpha=0.4, linestyle=:dash, shape=:rtriangle)
end
plots = init_plots(0, 50, 2; P0=[30])
__sweep__(plots, "lih30_uniform.each", "Uniform (2n/adapt)", :black)
__sweep__(plots, "lih30_bisection.one", "Bisectal (1/adapt)", 1)
__sweep__(plots, "lih30_nodes.one", "Optimal (1/adapt)", 5)
__sweep__(plots, "lih30_optimal", "Optimal (2/adapt)", 4)
# DUMMY PLOTS FOR LEGEND
for plot in plots
    Plots.plot!(plot, Int[], Int[]; default_series...,
        label="18 ns", color=:black, alpha=0.4, linestyle=:solid, shape=:ltriangle)
    Plots.plot!(plot, [], []; default_series...,
        label="22 ns", color=:black, alpha=0.8, linestyle=:solid, shape=:circle)
    Plots.plot!(plot, [], []; default_series...,
        label="36 ns", color=:black, alpha=0.4, linestyle=:solid, shape=:rtriangle)
end
save_plots(plots, "sensible")