#= How much do parameters change over the course of an adaptive protocol?

Since ctrl-VQE doesn't involve non-commuting operators,
    we can sensibly assign a *single* number to each pool operator,
    and plot it over the course of all adaptations.
The question is, *once* we've added an operator, how stable is it?

This script will parse a job, check all its optimized states,
    fill out a matrix with (nadapt, ndrive*npool) dimensions,
    and plot them.

The ultimate goal is to contrast the "frozen" runs from normal runs.
I'm guessing parameters *aren't* stable,
    and the "frozen" funs start alternating operators in a pitiful attempt
    to emulate the full optimization that normal runs can do without thinking about it.

Alas, this script only works when each candidate mode has one parameter.
Too much effort to generalize...

=#

import Plots

""" Load the given jobdir to _! and construct a parameter matrix. """
function load_harmonics!(jobdir)
    load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl") catch end

    nD = CtrlVQE.ndrives(_!.work.device)
    nP = length(_!.work.pool)
    nA = length(_!.trace.adaptations)

    parameters = zeros(JOB.Float, nA, nD*nP)
    for a in 1:nA
        unarchive!(JOB.adaptid(a))

        cursor = 1
        for i in 1:nD
            for n in _!.state.n[i]
                c = (i-1)*nP + n        # WHICH COLUMN?
                parameters[a,c] = _!.state.x[cursor]
                cursor += 1
            end
        end
    end

    return parameters
end

""" Load the given jobdir to _! and construct a parameter matrix, for coupled runs. """
function load_coupled!(jobdir)
    load!(jobdir; run=false)
    try include("$(_!.outdir)/script.jl") catch end

    nA = length(_!.trace.adaptations)

    parameters = zeros(JOB.Float, nA, nA)
    for a in 1:nA
        unarchive!(JOB.adaptid(a))
        L = length(_!.state.x)
        parameters[a,1:L] .= _!.state.x
    end

    return parameters
end

""" Fill in the active jobdir with parameter trajectories, using a given name. """
function make_coupledplots!(x, name)
    JOB.require_work(_!)
    nD = CtrlVQE.ndrives(_!.work.device)
    nP = length(_!.work.pool)
    ΩMAX = _!.setup.ΩMAX

    plt = Plots.plot(;
        xlabel = "Adaptations",
        ylabel = "Frequency (GHz)",     # TODO: Normalize gradient signal max->1, so parameters are all on the same scale.
        ylims = [-ΩMAX, ΩMAX] ./ 2π,
        legend = :bottomleft,
    )

    for l in axes(x,2)
        Plots.plot!(plt, x[:,l] ./ 2π;
            color = l, ls = :solid, lw = 3, alpha=0.8, label="n=$l",
        )
    end

    Plots.savefig(plt, "$(_!.outdir)/$(name).pdf")
end

""" Fill in the active jobdir with parameter trajectories, using a given name. """
function make_plots!(x, name)
    JOB.require_work(_!)
    nD = CtrlVQE.ndrives(_!.work.device)
    nP = length(_!.work.pool)
    ΩMAX = _!.setup.ΩMAX

    for i in 1:nD
        plt = Plots.plot(;
            xlabel = "Adaptations",
            ylabel = "Frequency (GHz)",
            ylims = [-ΩMAX, ΩMAX] ./ 2π,
            legend = :bottomleft,
        )

        for n in 1:nP
            c = ((i-1)*nP) + n

            βmode = n > (nP >> 1)

            color = βmode ? n - (nP >> 1) : n
            ls = βmode ? :dot : :solid
            label = βmode ? false : "n=$n"

            Plots.plot!(plt, x[:,c] ./ 2π;
                color = color, ls = ls, lw = 3, alpha=0.8, label=label,
            )
        end

        Plots.savefig(plt, "$(_!.outdir)/$(name).$(i).pdf")
    end
end



#= Make parameter plots for a given job directory. =#

code = "H215"
T = 15.0

x  = load_harmonics!("jobs/$(code)_harmonics/T$(T)")
make_plots!(x, "x")

xG  = load_coupled!("jobs/$(code)_gradients.exact/T$(T)")
make_coupledplots!(xG, "x")

xF = load_harmonics!("jobs/$(code)_harmonics.frozen/T$(T)")
xFc = cumsum(xF; dims=1)
make_plots!(xF, "dx")
make_plots!(xFc, "x")

