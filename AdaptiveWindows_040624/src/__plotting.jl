module Plotting
    import Plots

    import CtrlVQE

    import ..Vars, ..require_work, ..Calculations

    function plot_trace(vars::Vars, prefix::String, flags::Symbol...)
        require_work(vars)
        isnothing(vars.trace) && error("Trace is not initialized.")

        pathprefix = "$(vars.outdir)/$(prefix)"

        :optimization in flags && (
            plot = __trace__optimization(vars);
            Plots.savefig(plot, "$(pathprefix)_optimization.pdf")
        )
        :energies in flags && (
            plot = __trace__energies(vars);
            Plots.savefig(plot, "$(pathprefix)_energies.pdf")
        )
        :maxamplitudes in flags && (
            plot = __trace__maxamplitudes(vars);
            Plots.savefig(plot, "$(pathprefix)_maxamplitudes.pdf")
        )
        :frequencies in flags && (
            plot = __trace__frequencies(vars);
            Plots.savefig(plot, "$(pathprefix)_frequencies.pdf")
        )
        :adaptation in flags && (
            plot = __trace__adaptation(vars);
            Plots.savefig(plot, "$(pathprefix)_adaptation.pdf")
        )
        :parameters in flags && (
            plot = __trace__parameters(vars);
            Plots.savefig(plot, "$(pathprefix)_parameters.pdf")
        )
        :windows in flags && (
            plot = __trace__windows(vars);
            Plots.savefig(plot, "$(pathprefix)_windows.pdf")
        )

        return true
    end

    plot_trace(
        vars::Vars,
        prefix::String,
    ) = plot_trace(vars, prefix,
        :optimization, :energies, :maxamplitudes, :frequencies,
        :adaptation, :parameters, :windows,
    )

    function plot_pulse(
        vars::Vars, prefix::String, flags::Symbol...;
        only_α=false, only_E=false,
    )
        require_work(vars)

        pathprefix = "$(vars.outdir)/$(prefix)"

        :moduli in flags && (
            plot = __pulse__moduli(vars);
            Plots.savefig(plot, "$(pathprefix)_moduli.pdf")
        )
        :phases in flags && (
            plot = __pulse__phases(vars);
            Plots.savefig(plot, "$(pathprefix)_phases.pdf")
        )
        :amplitudes in flags && (
            plot = __pulse__amplitudes(vars; only_α=only_α);
            Plots.savefig(plot, "$(pathprefix)_amplitudes.pdf")
        )
        :gradients in flags && (
            plot = __pulse__gradients(vars; only_α=only_α);
            Plots.savefig(plot, "$(pathprefix)_gradients.pdf")
        )
        :trajectory in flags && (
            plot = __pulse__trajectory(vars; only_E=only_E);
            Plots.savefig(plot, "$(pathprefix)_trajectory.pdf")
        )

        return true
    end

    plot_pulse(
        vars::Vars,
        prefix::String;
        kwargs...
    ) = plot_pulse(vars, prefix, :moduli, :phases, :gradients, :trajectory; kwargs...)

    #= TODO: All the `plotargs` can get their own name and float to the top of each method.
    Or maybe even up here for full standardization? I like that alot. =#

    #= TODO: how to set custom color wheel? =#

    """ Plot physics-agnostic details of the optimization.

    The main feature is the gradient norm.
    It seems like its slope is the greatest litmus of optimization quality.

    Also of great interest is the degree to which the loss function decreases each step.
    Once it is consistently small, it is probably safe to terminate. Maybe...

    Each penalty function (modified by its weight) is plotted,
        to show if/when butting against the bounds slows the optimization down.

    Finally, the number of function and gradient calls are plotted.
    Sudden jumps in this plot diagnoses regions where the linesearch gets very difficult;
        often as not this diagnoses a bug in the code...

    """
    function __trace__optimization(vars; kwargs...)
        trace = vars.trace
        lossfn = vars.work.lossfn

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = 2e0
        Plots.plot!(plot;
            xlabel = "Iterations",
            xticks = integer_ticks(maximum(trace.iterations)),
            xminorticks = 5,
            xminorgrid = true,
            yscale = :log,
            ylims  = [1e-16, yMAX],
            yticks = 10.0 .^ (-16:2:0),
            yminorgrid = true,
            yminorticks = 1,
            legend = :left,
            kwargs...
        )

        twin = Plots.twinx(plot)
        twin_yMAX = ceil(Int, vars.setup.fnRATIO * maximum(trace.iterations))
        Plots.plot!(twin;
            ylabel = "Number of Calls",
            ylims = [0, twin_yMAX],
            yticks = integer_ticks(twin_yMAX),
            yminorticks = 5,
            legend = false,
        )

        # MARK WHERE EACH OPTIMIZATION TERMINATED
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.vline!(plot, trace.adaptations; label=false, plotargs...)

        # PLOT THE PRINCIPAL CURVES
        plotargs = (color=:black, lw=3, ls=:solid)
        Plots.plot!(plot, trace.iterations, trace.gd; label="|∇f|", plotargs...)

        plotargs = (color=:darkgray, lw=3, ls=:solid)
        Δf = abs.(diff(trace.fn))
        Plots.plot!(plot, trace.iterations[2:end], Δf; label="|Δf|", plotargs...)

        # PLOT (WEIGHTED) PENALTIES APPLIED THROUGHOUT ITERATION
        λ = lossfn.weights                                  # λ[i]  is weight of ith penalty
        Λ  = transpose(reduce(hcat, vars.trace.penalties))  # Λ[:,i] is trace of ith penalty
        for i in axes(Λ, 2)
            plotargs = (color=i, lw=2, ls=:solid)
            Plots.plot!(plot, trace.iterations, λ[i] .* Λ[:,i]; label="Λ$i", plotargs...)
        end

        # PLOT THE NUMBER OF CALLS
        plotargs = (color=:lightgray, lw=2, ls=:dash)
        Plots.plot!(plot, [0], [2yMAX]; label="# f", plotargs...)
        Plots.plot!(twin, trace.iterations, trace.f_calls; plotargs...)

        plotargs = (color=:lightgray, lw=2, ls=:dot)
        Plots.plot!(plot, [0], [2yMAX]; label="# ∇f", plotargs...)
        Plots.plot!(twin, trace.iterations, trace.g_calls; plotargs...)

        return plot
    end

    """ Plot how the energy changes over the optimization.

    The main conceit of this plot is mapping out the exact-diagonalized energy spectrum.
    If the energy trace hangs up anywhere,
        this tells us right away if it may be caught in an excited state.

    """
    function __trace__energies(vars; kwargs...)
        trace = vars.trace
        system = vars.work.system

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = 2e0
        Plots.plot!(plot;
            xlabel = "Iterations",
            xticks = integer_ticks(maximum(trace.iterations)),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Energy Error (Ha)",
            yscale = :log,
            ylims  = [1e-16, yMAX],
            yticks = 10.0 .^ (-16:2:0),
            yminorgrid = true,
            yminorticks = 1,
            legend = :bottomleft,
            kwargs...
        )

        # MARK WHERE EACH OPTIMIZATION TERMINATED
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.vline!(plot, trace.adaptations; label=false, plotargs...)

        # DOTTED LINES FOR EACH UNIQUE EXCITED STATE ENERGY IN RANGE
        Ex = system.model.Λ .- system.FCI   # EIGENSPECTRUM, SHIFTED
        filter!(x -> x < yMAX, Ex)          # DON'T TRY TO PLOT VALUES OFF THE PLOT
        unique!(Ex)                         # DON'T TRY TO PLOT THE SAME VALUE AGAIN
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.hline!(plot, Ex; label="Eigenspectrum", plotargs...)

        # THICKER MARK FOR REFERENCE STATE
        plotargs = (color=:black, lw=2, ls=:dash)
        Plots.hline!(plot, [system.REF-system.FCI]; label="Reference", plotargs...)

        # PLOT THE PRINCIPAL CURVE: ENERGY ERROR THROUGHOUT OPTIMIZATION
        ε = trace.energy .- system.FCI
        plotargs = (color=1, lw=3, ls=:solid)
        Plots.plot!(plot, trace.iterations, ε; label="Ansatz", plotargs...)

        return plot
    end

    """ Plot how the pulses grow throughout optimization.

    This is a more physics-inspired view
        into how the optimization handles the penalty terms.

    """
    function __trace__maxamplitudes(vars; kwargs...)
        trace = vars.trace
        ΩMAX = vars.setup.ΩMAX

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        Plots.plot!(plot;
            xlabel = "Iterations",
            xticks = integer_ticks(maximum(trace.iterations)),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Largest Pulse Amplitude (GHz)",
            ylims  = margined_lims(0, ΩMAX/2π),
            yminorgrid = true,
            yminorticks = 5,
            legend = false,
            kwargs...
        )

        # MARK WHERE EACH OPTIMIZATION TERMINATED
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.vline!(plot, trace.adaptations; label=false, plotargs...)

        # DOTTED LINE AT THE AMPLITUDE BOUND
        plotargs = (color=:black, lw=2, ls=:dot)
        Plots.hline!(plot, [ΩMAX/2π]; label="Amplitude Bound", plotargs...)

        # PLOT THE PRINCIPAL CURVES: Ωmax THROUGHOUT OPTIMIZATION, FOR EACH PULSE
        Ωmax = transpose(reduce(hcat, trace.Ωmax))
        for i in axes(Ωmax, 2)
            plotargs = (color=i, lw=3, ls=:solid)
            Plots.plot!(plot, trace.iterations, Ωmax[:,i] ./ 2π; label=false, plotargs...)
        end

        return plot
    end

    """ Plot how the drive frequencies change throughout optimization.

    The point of this plot is to visualize how drive frequencies relate,
        not just to the drive qubit, but to others.
    If state preparation is mostly entanglement between two particular qubits,
        we may be able to see a "cross-resonance" effect.

    """
    function __trace__frequencies(vars; kwargs...)
        trace = vars.trace
        ΔMAX = vars.setup.ΔMAX

        # FETCH ALL THE RESONANCE FREQUENCIES
        device = vars.work.device
        n = CtrlVQE.nqubits(device)
        ω = [CtrlVQE.resonancefrequency(device, q) for q in 1:n]

        # FETCH ALL THE DRIVE FREQUENCIES AND PICK YOUR BOUNDS
        ν = transpose(reduce(hcat, trace.ν))
        yMIN = min(minimum(ν), minimum(ω)) / 2π
        yMAX = max(maximum(ν), maximum(ω)) / 2π

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        Plots.plot!(plot;
            xlabel = "Iterations",
            xticks = integer_ticks(maximum(trace.iterations)),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Frequency (GHz)",
            ylims  = margined_lims(yMIN, yMAX),
            yminorgrid = true,
            yminorticks = 5,
            legend = :topleft,
            kwargs...
        )

        # MARK WHERE EACH OPTIMIZATION TERMINATED
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.vline!(plot, trace.adaptations; label=false, plotargs...)

        # DUMMY PLOTS FOR THE LEGEND
        ω_plotargs = (lw=:2, ls=:dash)
        ν_plotargs = (lw=:3, ls=:solid)
        Plots.plot!(plot, [0], [2yMAX]; label="Resonance", color=:black, ω_plotargs...)
        Plots.plot!(plot, [0], [2yMAX]; label="Drive", color=:black, ν_plotargs...)

        # PLOT THE ACTUAL FREQUENCIES THROUGHOUT OPTIMIZATION, FOR EACH PULSE
        for i in axes(ν, 2)
            q = CtrlVQE.drivequbit(device, i)
            Plots.hline!(plot, [ω[q]/2π]; label=false, color=i, ω_plotargs...)
            Plots.plot!(
                plot, trace.iterations, ν[:,i] ./ 2π;
                label=false, color=i, ν_plotargs...
            )
        end

        return plot
    end

    """ Plot how the energy and initial gradient converge with each adaptation.

    The main conceit of this plot is mapping out the exact-diagonalized energy spectrum.
    If the energy trace hangs up anywhere,
        this tells us right away if it may be caught in an excited state.

    """
    function __trace__adaptation(vars; kwargs...)
        trace = vars.trace
        system = vars.work.system

        nadapt = length(trace.adaptations)

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = 2e0
        Plots.plot!(plot;
            xlabel = "Adaptations",
            xticks = integer_ticks(nadapt),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Energy Error (Ha)",
            yscale = :log,
            ylims  = [1e-16, yMAX],
            yticks = 10.0 .^ (-16:2:0),
            yminorgrid = true,
            yminorticks = 1,
            legend = :bottomleft,
            kwargs...
        )

        # MARK OUT THE ADAPT CONVERGENCE CRITERION
        plotargs = (color=:gray, lw=2, ls=:dash, ms=2, msw=0, markershape=:circle)
        Plots.hline!(plot, [vars.meta.G_tol]; label=false, plotargs...)

        # PLOT THE LARGEST SCORE AT THE NEXT ADAPTATION
        G = trace.G_max
        plotargs = (color=:gray, lw=2, ls=:solid)
        Plots.plot!(plot, 1:nadapt, G; label="ADAPT Score", plotargs...)

        # DOTTED LINES FOR EACH UNIQUE EXCITED STATE ENERGY IN RANGE
        Ex = system.model.Λ .- system.FCI   # EIGENSPECTRUM, SHIFTED
        filter!(x -> x < yMAX, Ex)          # DON'T TRY TO PLOT VALUES OFF THE PLOT
        unique!(Ex)                         # DON'T TRY TO PLOT THE SAME VALUE AGAIN
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.hline!(plot, Ex; label="Eigenspectrum", plotargs...)

        # THICKER MARK FOR REFERENCE STATE
        plotargs = (color=:black, lw=2, ls=:dash)
        Plots.hline!(plot, [system.REF-system.FCI]; label="Reference", plotargs...)

        # PLOT THE PRINCIPAL CURVE: ENERGY ERROR FOR EACH OPTIMIZED ANSATZ
        ε = trace.energy[trace.adaptations] .- system.FCI
        plotargs = (color=1, lw=3, ls=:solid, ms=3, msw=0, markershape=:circle)
        Plots.plot!(plot, 1:nadapt, ε; label="Ansatz", plotargs...)

        return plot
    end

    """ Plot how the number of parameters in each pulse changes with each adaptation. """
    function __trace__parameters(vars; kwargs...)
        trace = vars.trace

        nadapt = length(trace.adaptations)
        parameters = transpose(reduce(hcat, vars.trace.parameters))
        poolsize = vars.trace.poolsize
        nP = vec(sum(parameters; dims=2))

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = max(maximum(nP), maximum(poolsize))
        Plots.plot!(plot;
            xlabel = "Adaptations",
            xticks = integer_ticks(nadapt),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Number of Parameters",
            yminorgrid = true,
            yticks = integer_ticks(yMAX),
            yminorticks = 5,
            legend = :topleft,
            kwargs...
        )

        # PLOT TOTAL POOL SIZES
        plotargs = (color=:gray, lw=2, ls=:solid, ms=3, msw=0, markershape=:circle)
        Plots.plot!(plot, 1:nadapt, poolsize; label="Pool Size", plotargs...)

        # PLOT TOTAL PARAMETER COUNTS
        plotargs = (color=:black, lw=2, ls=:solid, ms=3, msw=0, markershape=:circle)
        Plots.plot!(plot, 1:nadapt, nP; label="# Parameters", plotargs...)

        # PLOT PARAMETER COUNTS FOR EACH CURVE INDIVIDUALLY
        plotargs = (lw=2, ls=:solid, ms=3, msw=0, markershape=:circle)
        for i in axes(parameters, 2)
            Plots.plot!(plot, 1:nadapt, parameters[:,i]; color=i, label=false, plotargs...)
        end

        return plot
    end

    """ Plot how the number of parameters in each pulse changes with each adaptation. """
    function __trace__windows(vars; kwargs...)
        trace = vars.trace

        T = vars.setup.T
        nadapt = length(trace.adaptations)
        Δs_min = transpose(reduce(hcat, vars.trace.Δs_min))

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMIN = 1.0 / (vars.setup.fMAX * T)
        yMAX = 1.0

        k2MIN = ceil(Int, log2(yMIN))

        Plots.plot!(plot;
            xlabel = "Adaptations",
            xticks = integer_ticks(nadapt),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Smallest Window Duration",
            yscale = :log2,
            ylims  = [yMIN, yMAX],
            yticks = (2.0 .^ (k2MIN:0), ("T/2^$(abs(k))" for k in k2MIN:0)),
            yminorgrid = true,
            yminorticks = 1,
            legend = :bottomleft,
            kwargs...
        )

        plotargs = (lw=2, ls=:solid, ms=3, msw=0, markershape=:circle)
        for i in axes(Δs_min, 2)
            Plots.plot!(plot, 1:nadapt, Δs_min[:,i] ./ T; color=i, label=false, plotargs...)
        end

        return plot
    end




    """ Plot pulse amplitudes (absolute value) over pulse duration. """
    function __pulse__moduli(vars; kwargs...)
        ΩMAX = vars.setup.ΩMAX
        t = CtrlVQE.lattice(vars.work.grid)
        Ω, _, _ = Calculations.make_pulses(vars)

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = ΩMAX / 2π * 1.1
        Plots.plot!(plot;
            framestyle=:box,
            xlabel = "Time (ns)",
            xlims  = margined_lims(minimum(t), maximum(t); factor=0.00),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Pulse Modulus (GHz)",
            ylims  = margined_lims(0, yMAX; factor=0.00),
            yminorgrid = true,
            yminorticks = 5,
            legend = false,
            kwargs...
        )

        # DOTTED LINE AT THE AMPLITUDE BOUND
        plotargs = (color=:black, lw=2, ls=:dot)
        Plots.hline!(plot, [ΩMAX/2π]; label="Amplitude Bound", plotargs...)

        # PLOT |Ω| FOR EACH DRIVE
        plotargs = (ms=3, msw=0, markershape=:circle)
        for i in axes(Ω, 2)
            Plots.scatter!(plot, t, abs.(Ω[:,i]) ./ 2π; color=i, label=false, plotargs...)
        end

        return plot
    end

    """ Plot pulse phases over pulse duration. """
    function __pulse__phases(vars; kwargs...)
        ΩMAX = vars.setup.ΩMAX
        t = CtrlVQE.lattice(vars.work.grid)
        Ω, _, _ = Calculations.make_pulses(vars)

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        Plots.plot!(plot;
            framestyle = :box,
            xlabel = "Time (ns)",
            xlims  = margined_lims(minimum(t), maximum(t); factor=0.00),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Pulse Angle (π)",
            ylims  = margined_lims(-1.0, 1.0; factor=0.00),
            yticks = -1.0:0.5:1.0,
            yminorgrid = true,
            yminorticks = 5,
            legend = false,
            kwargs...
        )

        # STURDY LINES MARKING EACH QUADRANT
        plotargs = (color=:black, lw=0.75, ls=:solid)
        Plots.hline!(plot, -1.0:0.5:1.0; label=false, plotargs...)

        # PLOT ANGLE OF Ω FOR EACH DRIVE
        φ_args = (ms=3, msw=0, markershape=:circle)
        for i in axes(Ω, 2)
            Plots.scatter!(plot, t, angle.(Ω[:,i]) ./ π; color=i, label=false, φ_args...)
        end

        return plot
    end

    """ Plot real and imaginary parts of the pulse amplitudes over pulse duration. """
    function __pulse__amplitudes(vars; only_α=false, kwargs...)
        α_args = (lw=3, ls=:solid)
        β_args = (lw=3, ls=:dash)

        ΩMAX = vars.setup.ΩMAX
        t = CtrlVQE.lattice(vars.work.grid)
        _, α, β = Calculations.make_pulses(vars)

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = ΩMAX / 2π * 1.1
        Plots.plot!(plot;
            framestyle=:box,
            xlabel = "Time (ns)",
            xlims  = margined_lims(minimum(t), maximum(t); factor=0.00),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Pulse Amplitude (GHz)",
            ylims  =  margined_lims(-yMAX, yMAX),
            yminorgrid = true,
            yminorticks = 5,
            legend = only_α ? false : :topright,
            kwargs...
        )

        # DOTTED LINE AT THE AMPLITUDE BOUND
        plotargs = (color=:black, lw=2, ls=:dot)
        Plots.hline!(plot, [-ΩMAX/2π, ΩMAX/2π]; label="Amplitude Bound", plotargs...)

        # SOFT DOTTED LINE AT THE INSCRIBED AMPLITUDE BOUND
        if !only_α
            plotargs = (color=:black, lw=1, ls=:dot)
            ins = ΩMAX/2π/√2
            Plots.hline!(plot, [-ins, ins]; label="Inscribed Bound", plotargs...)

            # HACK IN THE LEGEND
            Plots.plot!(plot, [0], [2yMAX]; color=:black, label="α", α_args...)
            Plots.plot!(plot, [0], [2yMAX]; color=:black, label="β", β_args...)
        end

        # PLOT |Ω| FOR EACH DRIVE
        for i in axes(α, 2)
            Plots.plot!(plot, t, α[:,i] ./ 2π; color=i, label=false, α_args...)
            only_α || Plots.plot!(plot, t, β[:,i] ./ 2π; color=i, label=false, β_args...)
        end

        return plot
    end

    """ Plot gradient signals over pulse duration. """
    function __pulse__gradients(vars; only_α=false, kwargs...)
        α_args = (lw=3, ls=:solid)
        β_args = (lw=3, ls=:dash)

        t = CtrlVQE.lattice(vars.work.grid)
        _, ϕα, ϕβ = Calculations.make_gradientsignals(vars)

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMAX = max(maximum(abs.(ϕα)), maximum(abs.(ϕβ))) * 2π
        Plots.plot!(plot;
            framestyle = :box,
            xlabel = "Time (ns)",
            xlims  = margined_lims(minimum(t), maximum(t); factor=0.00),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Gradient Signal (Ha/GHz)",
            ylims  =  margined_lims(-yMAX, yMAX),
            yminorgrid = true,
            yminorticks = 5,
            legend = only_α ? false : :topright,
            kwargs...
        )

        # STURDY LINE MARKING ZERO
        plotargs = (color=:black, lw=0.75, ls=:solid)
        Plots.hline!(plot, [0.0]; label=false, plotargs...)

        # HACK IN THE LEGEND
        if !only_α
            Plots.plot!(plot, [0], [2yMAX]; color=:black, label="ϕα", α_args...)
            Plots.plot!(plot, [0], [2yMAX]; color=:black, label="ϕβ", β_args...)
        end

        # PLOT ANGLE OF Ω FOR EACH DRIVE
        for i in axes(ϕα, 2)
            Plots.plot!(plot, t, ϕα[:,i] .* 2π; color=i, label=false, α_args...)
            only_α || Plots.plot!(plot, t, ϕβ[:,i] .* 2π; color=i, label=false, β_args...)
        end

        return plot
    end

    """ Plot energy trajectory over pulse duration, perhaps with norm. """
    function __pulse__trajectory(vars; only_E=nothing, kwargs...)
        # BY DEFAULT, DO NOT BOTHER WITH NORM IF m=2
        isnothing(only_E) && (only_E = vars.setup.m == 2)

        t = CtrlVQE.lattice(vars.work.grid)
        E = Calculations.make_energytrajectory(vars)
        only_E || (F = Calculations.make_normtrajectory(vars))

        # INITIALIZE THE PLOT
        plot = Plots.plot()
        yMIN = vars.work.system.FCI
        yMAX = maximum(E)
        Plots.plot!(plot;
            framestyle = :box,
            xlabel = "Time (ns)",
            xlims  = margined_lims(minimum(t), maximum(t); factor=0.00),
            xminorticks = 5,
            xminorgrid = true,
            ylabel = "Energy (Ha)",
            ylims  =  margined_lims(yMIN, yMAX),
            yminorgrid = true,
            yminorticks = 5,
            legend = only_E ? false : :topright,
            kwargs...
        )

        if !only_E
            twin = Plots.twinx(plot)
            Plots.plot!(twin;
                ylabel = "Leakage",
                ylims = margined_lims(0.0, 1.0; factor=0.05),
                yminorticks = 5,
                legend = false,
            )
        end

        # DOTTED LINES FOR EACH UNIQUE EXCITED STATE ENERGY IN RANGE
        Ex = vars.work.system.model.Λ       # EIGENSPECTRUM, SHIFTED
        filter!(x -> x < yMAX, Ex)          # DON'T TRY TO PLOT VALUES OFF THE PLOT
        unique!(Ex)                         # DON'T TRY TO PLOT THE SAME VALUE AGAIN
        plotargs = (color=:black, lw=1, ls=:dot)
        Plots.hline!(plot, Ex; label="Eigenspectrum", plotargs...)

        # THICKER MARK FOR REFERENCE STATE
        plotargs = (color=:black, lw=2, ls=:dash)
        Plots.hline!(plot, [vars.work.system.REF]; label="Reference", plotargs...)

        # PLOT THE ENERGY TRAJECTORY
        plotargs = (color=1, lw=3, ls=:solid)
        Plots.plot!(plot, t, E; label="Trajectory", plotargs...)

        # PLOT THE LEAKAGE TRAJECTORY
        if !only_E
            plotargs = (color=2, lw=3, ls=:solid)
            Plots.plot!(plot, [0], [2yMAX]; label="Leakage", plotargs...)
            Plots.plot!(twin, t, 1 .- F; plotargs...)
        end

        return plot
    end







    """ Buffer the given min/max values so points at these values plot nicely. """
    function margined_lims(a, b; factor=0.1)
        rng = b - a
        margin = rng * factor
        return (a-margin, b+margin)
    end

    """ Calculate sensible ticks for a maximal integer.

    Ticks will range from 0 up to (not necessarily including) M.
    Intervals will be at either a power of 10, or half that.
    Therefore, if minor ticks are desired,
        selecting minorticks=5 will be a sensible choice.
    """
    function integer_ticks(M)
        expanded = digits(M)
        c = last(expanded)
        K = length(expanded) - 1
        K == 0 && return 0:M

        interval = c ≤ 5 ? (10^K)÷2 : (10^K)

        return 0:interval:M
    end

end