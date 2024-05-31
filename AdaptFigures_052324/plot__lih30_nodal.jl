#= The idea here is to plot nparam and smallest window for each drive,
    at the two different times. =#

using AdaptFigures_052324

import Plots

prefix = "lih30_nodal"
JOB = SQUARE
SQUAREJOBS = "../AdaptiveWindows_040624/jobs"
DRIVE_STYLE(i,n) = (
    seriescolor = i % n == 0 ? n : i % n,
    linestyle = i ≤ n ? :solid : :dot,
)

function make_Pva(a, ε, P; n, T)
    plt = Plots.plot(); twin = Plots.twinx(plt)
    Plots.plot!(plt;
        PLOT_STYLE...,

        xlabel = "Adaptation",
        xlims = [0,26],
        xticks = 0:5:25,
        xminorticks = 5,

        ylabel = "Parameters",
        ylims  = [0, 6],
        yticks = 0:6,
        yminorticks = 1,

        legend = :topleft,
    )
    Plots.plot!(twin;
        ylabel = "Energy Error (Ha)",
        yscale = :log,
        ylims  = [1e-16, 1e2],
        yticks = 10.0 .^ (-16:2:0),
        yminorticks = 1,
        ygrid = false,
    )
    hline!_energy(twin; alpha=0.5)

    plot!_curve(twin, a, ε; alpha=0.5, shape=:none, label="Energy Error")


    for i in axes(P, 2)
        plot!_curve(plt, a, P[:,i]; DRIVE_STYLE(i,n)...)
    end

    savepdf(plt, prefix, "Pva_T$(T)")
end

function make_Wva(a, ε, W; n, T, fMAX)
    sMIN = 1.0 / (fMAX * T)
    sMAX = 1.0
    k2MIN = ceil(Int, log2(sMIN))

    plt = Plots.plot(); twin = Plots.twinx(plt)
    Plots.plot!(plt;
        PLOT_STYLE...,

        xlabel = "Adaptation",
        xlims = [0,26],
        xticks = 0:5:25,
        xminorticks = 5,

        ylabel = "Smallest Window Duration",
        yscale = :log2,
        ylims  = [sMIN, sMAX],
        yticks = (2.0 .^ (k2MIN:0), ("T/2^$(abs(k))" for k in k2MIN:0)),
        yminorticks = 1,

        legend = :topright,
    )
    Plots.plot!(twin;
        ylabel = "Energy Error (Ha)",
        yscale = :log,
        ylims  = [1e-16, 1e2],
        yticks = 10.0 .^ (-16:2:0),
        yminorticks = 1,
        ygrid = false,
    )
    hline!_energy(twin; alpha=0.5)

    plot!_curve(twin, a, ε; alpha=0.5, shape=:none, label="Energy Error")


    for i in axes(W, 2)
        plot!_curve(plt, a, W[:,i] ./ T; DRIVE_STYLE(i,n)...)
    end

    savepdf(plt, prefix, "Wva_T$(T)")
end

function make(jobdir)
    JOB.load!(jobdir)
    system = JOB.System(JOB._!.setup.code)

    n = system.n
    T = JOB._!.setup.T
    fMAX = JOB._!.setup.fMAX

    ε = JOB._!.trace.energy[JOB._!.trace.adaptations] .- system.FCI
    ε = [ε_ < eps(JOB.Float) ? eps(JOB.Float) : ε_ for ε_ in ε]

    a = collect(1:length(JOB._!.trace.adaptations))
    P = transpose(reduce(hcat, JOB._!.trace.parameters))
    W = transpose(reduce(hcat, JOB._!.trace.Δs_min))

    make_Pva(a, ε, P; n=n, T=T)
    make_Wva(a, ε, W; n=n, T=T, fMAX=fMAX)
end


make("$SQUAREJOBS/lih30_nodes.one/T22.0")
make("$SQUAREJOBS/lih30_nodes.one/T36.0")
