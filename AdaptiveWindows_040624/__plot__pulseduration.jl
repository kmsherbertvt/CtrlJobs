#= Show multiple trajectories side-by-side in a single plot vs pulse duration.

Sensible things to plot:
- Energy error.
- Total parameters.
- Total iterations.
- Total adaptations.

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

__initplot__parameters(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Parameter Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topleft,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__iterations(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Iteration Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topleft,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__adaptations(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Adaptation Count",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topleft,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

function init_plots(
    TMIN, TMAX, Δs;
    energy=Dict(), param=Dict(), iter=Dict(), adapt=Dict(), kwargs...
)
    T_axis = __initplot__T_axis(TMIN, TMAX, Δs)
    return (
        energy = __initplot__energy(; T_axis..., kwargs..., energy...),
        parameters = __initplot__parameters(; T_axis..., kwargs..., param...),
        iterations = __initplot__iterations(; T_axis..., kwargs..., iter...),
        adaptations = __initplot__adaptations(; T_axis..., kwargs..., adapt...),
    )
end

function plot_survey!(
    plots, surveydir;
    energy=Dict(), param=Dict(), iter=Dict(), adapt=Dict(), kwargs...
)
    system = Ref{JOB.System}()

    # INITIALIZE RESULT VECTORS
    data = (
        T = JOB.Float[],
        energy = JOB.Float[],
        parameters = Int[],
        iterations = Int[],
        adaptations = Int[],
    )

    for jobdir in readdir(surveydir, join=true)
        endswith(jobdir, ".DS_Store") && continue
        # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
        vars = JOB.load(jobdir)
        isnothing(vars.trace) && continue
        isassigned(system) || (system[] = JOB.System(vars.setup.code))

        # REGISTER THE DATA
        push!(data.T, vars.setup.T)
        ε = last(vars.trace.energy) - system[].FCI
        push!(data.energy, ε < eps(JOB.Float) ? eps(JOB.Float) : ε)
            #= NOTE: Beware!
                ctrl-ADAPT-VQE does an even better job of finding eigenstates,
                    it seems, than "exact diagonalization".
                Meaning, we are actually getting energies LESS than that of Λ0.
                Like. 1e-14 less. But that is more than 1e-16 less. >_>
            =#
        push!(data.parameters, length(vars.state.x))
        push!(data.iterations, last(vars.trace.iterations))
        push!(data.adaptations, length(vars.trace.adaptations))

    end

    # SORT OUTPUTS BY TIME
    σ = sortperm(data.T)
    for array in data; permute!(array, σ); end

    # PLOT DATA
    Plots.plot!(plots.energy, data.T, data.energy;
            default_series..., kwargs..., energy...)
    Plots.plot!(plots.parameters, data.T, data.parameters;
            default_series..., kwargs..., param...)
    Plots.plot!(plots.iterations, data.T, data.iterations;
            default_series..., kwargs..., iter...)
    Plots.plot!(plots.adaptations, data.T, data.adaptations;
            default_series..., kwargs..., iter...)
end

function save_plots(plots, label)
    for (kind, plot) in pairs(plots)
        Plots.savefig(plot, "figs/$(label)__$(kind).pdf")
    end
end

# CONTRASTING VIABLILITY
plots = init_plots(0.0, 48.0, 3.0)
plot_survey!(plots, "jobs/lih30_bisection";
    label="Bisectal", color=1, linestyle=:solid, shape=:circle)
plot_survey!(plots, "jobs/lih30_random";
    label="Random", color=2, linestyle=:solid, shape=:hexagon)
plot_survey!(plots, "jobs/lih30_nodes";
    label="Nodal", color=3, linestyle=:solid, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_optimal";
    label="Optimal", color=4, linestyle=:solid, shape=:ltriangle)
save_plots(plots, "viability")

# CONTRASTING VIABLILITY - EACH
plots = init_plots(0.0, 48.0, 3.0)
plot_survey!(plots, "jobs/lih30_bisection.each";
    label="Bisectal", color=1, linestyle=:dash, shape=:circle)
plot_survey!(plots, "jobs/lih30_random.each";
    label="Random", color=2, linestyle=:dash, shape=:hexagon)
plot_survey!(plots, "jobs/lih30_nodes.each";
    label="Nodal", color=3, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_optimal.each";
    label="Optimal", color=4, linestyle=:dash, shape=:ltriangle)
plot_survey!(plots, "jobs/lih30_uniform.each";
    label="Uniform", color=:black, linestyle=:dash, shape=:square)
save_plots(plots, "viability.each")

# CONTRASTING VIABLILITY - ALL
plots = init_plots(0.0, 48.0, 3.0)
plot_survey!(plots, "jobs/lih30_bisection.all";
    label="Bisectal", color=1, linestyle=:dot, shape=:circle)
plot_survey!(plots, "jobs/lih30_random.all";
    label="Random", color=2, linestyle=:dot, shape=:hexagon)
plot_survey!(plots, "jobs/lih30_nodes.all";
    label="Nodal", color=3, linestyle=:dot, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_optimal.all";
    label="Optimal", color=4, linestyle=:dot, shape=:ltriangle)
save_plots(plots, "viability.all")

# CONTRASTING PARALLELISM
plots = init_plots(0.0, 48.0, 3.0)
plot_survey!(plots, "jobs/lih30_optimal";
    label="One", color=4, linestyle=:solid, shape=:ltriangle)
plot_survey!(plots, "jobs/lih30_optimal.each";
    label="One per Pulse", color=4, linestyle=:dash, shape=:ltriangle)
plot_survey!(plots, "jobs/lih30_optimal.all";
    label="One per Window", color=4, linestyle=:dot, shape=:ltriangle)
plot_survey!(plots, "jobs/lih30_uniform.each";
    label="Uniform", color=:black, linestyle=:dash, shape=:square)
save_plots(plots, "parallel")