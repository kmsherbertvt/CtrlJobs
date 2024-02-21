module UniformWindows_110123
    export _!, unload!, load!, load, save
    export iterid, archive, unarchive!
    export inspect, report, plot_trace, plot_pulse

    const Float = Float64

    ######################################################################################
    #= Load matrices and reference states given a system code. =#
    include("__system.jl")
    import .Systems: MODEL_PATH, SECTORPATH, EIGEN_PATH
    import .Systems: Model, System

    ######################################################################################
    #= Define all the variables in the module, and how to display them. =#
    include("__variables.jl")
    import .Variables: Vars, SetupVars, MetaVars, TraceVars, StateVars, WorkVars
    import .Variables: require_work, initial_state, initial_trace
    import .Variables: update_trace!, is_converged, is_terminated
    import .Variables: indices

    # DECLARE AN ACTIVE VARIABLE, to be used by just about all functions.
    const _! = Vars()
    Variables.initial_state() = initial_state(_!)
    Variables.initial_trace() = initial_trace(_!)
    Variables.update_trace!(args...) = update_trace!(_!, args...)
    Variables.is_converged() = is_converged(_!)
    Variables.is_terminated() = is_terminated(_!)

    ######################################################################################
    #= Define physics-agnostic methods for manipulating variable objects. =#
    include("__bookkeeping.jl")
    import .BookKeeping: unload!, load!, load, save, archive, unarchive!, iterid
    BookKeeping.unload!() = unload!(_!)
    BookKeeping.load!(args...; kwargs...) = load!(_!, args...; kwargs...)
    BookKeeping.save() = save(_!)
    BookKeeping.archive(name) = archive(_!, name)
    BookKeeping.unarchive!(name) = unarchive!(_!, name)

    ######################################################################################
    #= Define some calculation-intensive methods needed for plotting and such. =#
    include("__calculations.jl")
    import .Calculations: make_energytrajectory, make_normtrajectory
    import .Calculations: make_pulses, make_gradientsignals
    Calculations.make_energytrajectory() = make_energytrajectory(_!)
    Calculations.make_normtrajectory() = make_normtrajectory(_!)
    Calculations.make_pulses() = make_pulses(_!)
    Calculations.make_gradientsignals() = make_gradientsignals(_!)

    ######################################################################################
    #= Standard ASCII displays for a given state. =#
    include("__printing.jl")
    import .Printing: inspect, report
    Printing.inspect(args...) = inspect(_!, args...)
    Printing.report(args...) = report(_!, args...)

    ######################################################################################
    #= Standard Plots plotting for a given state. =#
    include("__plotting.jl")
    import .Plotting: plot_trace, plot_pulse
    Plotting.plot_trace(args...) = plot_trace(_!, args...)
    Plotting.plot_pulse(args...; kwargs...) = plot_pulse(_!, args...; kwargs...)

    ######################################################################################
    #= Standard ssh workflow: Run `./juliarepl`, then any of the following functions. =#

    function init(outdir::String, script::String, code::String, T::Float, W::Int)
        ispath(outdir) && error("Data already exists in $outdir.")

        vars = Vars(outdir=outdir, setup = SetupVars(code=code, T=T, W=W))
        save(vars)
        cp("$script.jl", "$(vars.outdir)/script.jl")
    end

    function init_survey(mode, script::String, code::String, Ts::Float...)
        surveydir = "jobs/$(code)_$(script)_$(mode.first)$(mode.second)"
        !isdir(surveydir) && mkdir(surveydir)

        ARG = mode.second
        choose_W = (
            perns = T -> round(Int, T*ARG),
            ΔsMAX = T -> round(Int, T÷ARG),
        )

        for T in Ts
            W = choose_W[mode.first](T)
            outdir = "$surveydir/T$(T)_W$(W)"
            try init(outdir, script, code, T, W)
            catch ERROR println(ERROR.msg) end
        end
    end

    """ Construct a list of genuine "job" directories within a path. """
    function eachjob(path="jobs")
        jobs = String[]
        queue = [path]
        while !isempty(queue)
            candidate = pop!(queue)
            isdir(candidate) || continue
            append!(queue, readdir(candidate; join=true))

            # ATTEMPT TO LOAD VARIABLES - Success means this directory is a "job".
            try load(candidate) catch; continue; end
            push!(jobs, candidate)
        end
        return jobs
    end

    function status(path="jobs")
        println("""
        STATUS OF ALL JOBS IN $path:

        C    \t = Converged
        T *  \t = Terminated (but not converged)
        R ** \t = Running
        - *** \t = pry ought to be run ^_^

        """)
        for job in eachjob(path)
            vars = load(job)
            flag = is_converged(vars) ? "C   " :
                is_terminated(vars) ? "T *  " :
                isfile("$(vars.outdir)/running") ? "R ** " :
                "- *** "
            if !isnothing(vars.trace) && !isempty(vars.trace.iterations)
                I = last(vars.trace.iterations)
                E = last(vars.trace.energy)
                fn = last(vars.trace.fn)
                gd = last(vars.trace.gd)
                println("$flag\t$(vars.outdir)\t\tI=$I\tE=$E\tfn=$fn\tgd=$gd")
            else
                println("$flag\t$(vars.outdir)")
            end
        end
    end

    function unterminated_jobs(path="jobs")
        return filter(
            job -> !is_terminated(load(job)),
            eachjob(path),
        )
    end

    function inspect_all(path="jobs")
        for job in eachjob(path)
            inspect(load(job))
        end
    end

    function summarize_all(path="jobs")
        for job in eachjob(path)
            summarize(load(job))
        end
    end

    function start(path="jobs")
        for job in unterminated_jobs(path)
            try run(`./start $job`) catch end
        end
    end

end # module UniformWindows_110123
