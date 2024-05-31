module AdaptFigures_052324
    export SQUARE, MODAL
    export load_IvP, load_IvT, load_PvT, load_εvP, load_εvT
    export init_IvP, init_IvT, init_PvT, init_εvP, init_εvT
    export vline!_dimension, hline!_dimension, hline!_energy
    export label!, plot!_curve, savepdf
    export PLOT_STYLE, CURVE_STYLE

    import Plots, ColorSchemes
    import AdaptiveWindows_040624 as SQUARE
    import AdaptiveModes_041724 as MODAL

    const PLOT_STYLE = (
        color_palette = ColorSchemes.tab10,
        xminorgrid = true,
        yminorgrid = true,
        tickfontsize = 11,
        legendfontsize = 11,
    )

    const CURVE_STYLE = (
        seriescolor = :black,
        seriesalpha = 0.8,
        linestyle = :solid,
        linewidth = 3,
        markershape = :circle,
        markersize = 5,
        markerstrokewidth = 0,
        label = false,
    )

    ######################################################################################
    #= DATA LOADING =#

    function load_IvP(JOB, jobdir)
        # INITIALIZE DATA ARRAYS
        Ps = Int[]
        Is = Int[]

        # LOAD VARIABLES AND THE SYSTEM
        vars = JOB.load(jobdir)
        nadapt = length(vars.trace.adaptations)

        # LOAD THE STATE AT EACH ADAPTATION, AND REGISTER THE DATA
        for a in 1:nadapt
            JOB.unarchive!(vars, JOB.adaptid(a))
            push!(Ps, length(vars.state.x))
            push!(Is, vars.trace.adaptations[a])
        end

        # SORT OUTPUTS BY PARAMETERS
        σ = sortperm(Ps)
        permute!(Ps, σ)
        permute!(Is, σ)

        return Ps, Is
    end

    function load_IvT(JOB, surveydir; εtgt=0.00159)
        system = Ref{JOB.System}()

        # INITIALIZE DATA ARRAYS
        Ts = JOB.Float[]
        Is = Int[]

        for jobdir in readdir(surveydir, join=true)
            endswith(jobdir, ".DS_Store") && continue
            # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
            vars = JOB.load(jobdir)
            isnothing(vars.trace) && continue
            isassigned(system) || (system[] = JOB.System(vars.setup.code))

            # COMPUTE THE ENERGY ERRORS AT EACH ADAPTATION
            εs = vars.trace.energy[vars.trace.adaptations] .- system[].FCI
            a = findfirst(εs .≤ εtgt)       # FIRST ADAPTATION BELOW THRESHOLD
            isnothing(a) && continue        # SKIP IF THIS PULSE DURATION ALWAYS FAILED

            # REGISTER THE DATA
            push!(Ts, vars.setup.T)
            push!(Is, vars.trace.adaptations[a])
        end

        # SORT OUTPUTS BY TIME
        σ = sortperm(Ts)
        permute!(Ts, σ)
        permute!(Is, σ)

        return Ts, Is
    end

    function load_PvT(JOB, surveydir; εtgt=0.00159)
        system = Ref{JOB.System}()

        # INITIALIZE DATA ARRAYS
        Ts = JOB.Float[]
        Ps = Int[]

        for jobdir in readdir(surveydir, join=true)
            endswith(jobdir, ".DS_Store") && continue
            # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
            vars = JOB.load(jobdir)
            isnothing(vars.trace) && continue
            isassigned(system) || (system[] = JOB.System(vars.setup.code))

            # COMPUTE THE ENERGY ERRORS AT EACH ADAPTATION
            εs = vars.trace.energy[vars.trace.adaptations] .- system[].FCI
            a = findfirst(εs .≤ εtgt)       # FIRST ADAPTATION BELOW THRESHOLD
            isnothing(a) && continue        # SKIP IF THIS PULSE DURATION ALWAYS FAILED

            # REGISTER THE DATA
            push!(Ts, vars.setup.T)

            #= TODO: Decide between these.
                They should be the same except for "frozen" runs.
                The question is, do we care about cumulative parameter count
                    or instantaneous difficulty of optimiztion..?

                The latter seems "fairer" to the frozen methods.

                I should get Nick's opinion

                =#
            # push!(Ps, sum(vars.trace.parameters[a]))    # CUMULATIVE

            JOB.unarchive!(vars, JOB.adaptid(a))
            push!(Ps, length(vars.state.x))             # INSTANTANEOUS
        end

        # SORT OUTPUTS BY TIME
        σ = sortperm(Ts)
        permute!(Ts, σ)
        permute!(Ps, σ)

        return Ts, Ps
    end

    function load_εvP(JOB, jobdir)
        # INITIALIZE DATA ARRAYS
        Ps = Int[]
        εs = JOB.Float[]

        # LOAD VARIABLES AND THE SYSTEM
        vars = JOB.load(jobdir)
        system = JOB.System(vars.setup.code)
        nadapt = length(vars.trace.adaptations)

        # LOAD THE STATE AT EACH ADAPTATION, AND REGISTER THE DATA
        for a in 1:nadapt
            JOB.unarchive!(vars, JOB.adaptid(a))
            push!(Ps, length(vars.state.x))

            i = vars.trace.adaptations[a]
            ε = vars.trace.energy[i] - system.FCI
            push!(εs, ε < eps(JOB.Float) ? eps(JOB.Float) : ε)
        end

        # SORT OUTPUTS BY PARAMETERS
        σ = sortperm(Ps)
        permute!(Ps, σ)
        permute!(εs, σ)

        return Ps, εs
    end

    function load_εvT(JOB, surveydir)
        system = Ref{JOB.System}()

        # INITIALIZE DATA ARRAYS
        Ts = JOB.Float[]
        εs = JOB.Float[]

        for jobdir in readdir(surveydir, join=true)
            endswith(jobdir, ".DS_Store") && continue
            # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
            vars = JOB.load(jobdir)
            isnothing(vars.trace) && continue
            isassigned(system) || (system[] = JOB.System(vars.setup.code))

            # REGISTER THE DATA
            push!(Ts, vars.setup.T)
            ε = last(vars.trace.energy) - system[].FCI
            push!(εs, ε < eps(JOB.Float) ? eps(JOB.Float) : ε)
                #= NOTE: Beware!
                    ctrl-ADAPT-VQE does an even better job of finding eigenstates,
                        it seems, than "exact diagonalization".
                    Meaning, we are actually getting energies LESS than that of Λ0.
                    Like. 1e-14 less. But that is more than 1e-16 less. >_>
                =#
        end

        # SORT OUTPUTS BY TIME
        σ = sortperm(Ts)
        permute!(Ts, σ)
        permute!(εs, σ)

        return Ts, εs
    end

    ######################################################################################
    #= PLOT SETUP =#

    function init_IvP(; N=nothing, kwargs...)
        plt = Plots.plot(;
            PLOT_STYLE...,
            xlabel = "Parameter Count",
            ylabel = "Iteration Count",
            kwargs...
        )
        isnothing(N) || vline!_dimension(plt, N)
        return plt
    end

    function init_IvT(; kwargs...)
        plt = Plots.plot(;
            PLOT_STYLE...,
            xlabel = "Pulse Duration (ns)",
            ylabel = "BFGS Iterations to Chemical Accuracy",
            kwargs...
        )
        return plt
    end

    function init_PvT(; N=nothing, kwargs...)
        plt = Plots.plot(;
            PLOT_STYLE...,
            xlabel = "Pulse Duration (ns)",
            ylabel = "Parameters to Chemical Accuracy",
            kwargs...
        )
        isnothing(N) || hline!_dimension(plt, N)
        return plt
    end

    function init_εvP(; N=nothing, ε=0.00159, kwargs...)
        plt = Plots.plot(;
            PLOT_STYLE...,
            xlabel = "Parameter Count",
            ylabel = "Energy Error (Ha)",
            yscale = :log,
            ylims  = [1e-16, 1e2],
            yticks = 10.0 .^ (-16:2:0),
            yminorgrid = false,
            kwargs...
        )
        isnothing(N) || vline!_dimension(plt, N)
        isnothing(ε) || hline!_energy(plt, ε)
        return plt
    end

    function init_εvT(; ε=0.00159, kwargs...)
        plt = Plots.plot(;
            PLOT_STYLE...,
            xlabel = "Pulse Duration (ns)",
            ylabel = "Energy Error (Ha)",
            yscale = :log,
            ylims  = [1e-16, 1e2],
            yticks = 10.0 .^ (-16:2:0),
            yminorgrid = false,
            kwargs...
        )
        isnothing(ε) || hline!_energy(plt, ε)
        return plt
    end

    ######################################################################################
    #= COMMON RECIPES =#

    """ Shorthand to mark a key number of parameters on a _vP plot. """
    vline!_dimension(plt, N; kwargs...) = Plots.vline!(
        plt, [N];
        ls=:solid, lw=1, color=:black, label=false,
        kwargs...
    )

    """ Shorthand to mark a key number of parameters on a Pv_ plot. """
    hline!_dimension(plt, N; kwargs...) = Plots.hline!(
        plt, [N];
        ls=:solid, lw=1, color=:black, label=false,
        kwargs...
    )

    """ Shorthand to mark a key energy on a εv_ plot. """
    hline!_energy(plt, ε=0.00159; kwargs...) = Plots.hline!(
        plt, [ε];
        ls=:solid, lw=1, color=:black, label=false,
        kwargs...
    )

    """ Shorthand to hack in a label without any data. """
    label!(plt, label; kwargs...) = Plots.plot!(
        plt, [-1], [-1];
        CURVE_STYLE...,
        label=label,
        kwargs...
    )

    """ Shorthand to plot a curve with default style. """
    plot!_curve(plt, x, y; kwargs...) = Plots.plot!(
        plt, x, y;
        CURVE_STYLE...,
        kwargs...
    )

    """ Shorthand to save figure as a pdf to the figs folder. """
    savepdf(plt, prefix, plot) = Plots.savefig(plt, "figs/$(prefix)_$(plot).pdf")

end # module AdaptFigures_052324
