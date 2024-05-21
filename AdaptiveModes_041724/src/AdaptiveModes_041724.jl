module AdaptiveModes_041724
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

    # DECLARE AN ACTIVE VARIABLE, to be used by just about all functions.
    const _! = Vars()
    Variables.initial_state() = initial_state(_!)
    Variables.initial_trace() = initial_trace(_!)

    ######################################################################################
    #= Define physics-agnostic methods for manipulating variable objects. =#
    include("__bookkeeping.jl")
    import .BookKeeping: unload!, load!, load, save, archive, unarchive!, iterid, adaptid
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
    #= Standard functions to make Optim tick. =#
    include("__optimizations.jl")
    import .Optimizations: make_objective, make_state, make_callback, make_options
    import .Optimizations: update_trace!
    import .Optimizations: optimization_is_converged, optimization_is_terminated
    Optimizations.make_objective(args...) = make_objective(_!, args...)
    Optimizations.make_state(args...) = make_state(_!, args...)
    Optimizations.make_callback(args...) = make_callback(_!, args...)
    Optimizations.make_options(args...) = make_options(_!, args...)
    Optimizations.update_trace!(args...) = update_trace!(_!, args...)
    Optimizations.optimization_is_converged() = optimization_is_converged(_!)
    Optimizations.optimization_is_terminated() = optimization_is_terminated(_!)

    ######################################################################################
    #= Standard functions to extend pulses adaptively. =#
    include("__adaptations.jl")
    import .Adaptations: makepool_harmonics, makepool_complexharmonics
    import .Adaptations: makepool_square, makepool_complexsquare
    import .Adaptations: upHk_naive, upHk_slate, upHk_recycle
    import .Adaptations: select_one, select_oneperpulse, select_iterative
    Adaptations.select_one(args...) = select_one(_!, args...)
    Adaptations.select_oneperpulse(args...) = select_oneperpulse(_!, args...)
    Adaptations.select_iterative(args...) = select_iterative(_!, args...)
    import .Adaptations: prepare_candidate, score_candidates
    import .Adaptations: adapt_trace!, adapt_state!, adapt_work!
    import .Adaptations: adapt_is_converged, adapt_is_terminated
    Adaptations.prepare_candidate(args...) = prepare_candidate(_!, args...)
    Adaptations.score_candidates(args...) = score_candidates(_!, args...)
    Adaptations.adapt_trace!(args...) = adapt_trace!(_!, args...)
    Adaptations.adapt_state!(args...) = adapt_state!(_!, args...)
    Adaptations.adapt_work!(args...) = adapt_work!(_!, args...)
    Adaptations.adapt_is_converged() = adapt_is_converged(_!)
    Adaptations.adapt_is_terminated() = adapt_is_terminated(_!)


    ######################################################################################
    #= Standard ssh workflow: Run `./juliarepl`, then any of the following functions. =#

    function init(outdir::String, script::String, code::String, T::Float)
        ispath(outdir) && error("Data already exists in $outdir.")

        vars = Vars(outdir=outdir, setup = SetupVars(code=code, T=T))
        save(vars)
        cp("$script.jl", "$(vars.outdir)/script.jl")
    end

    function init_survey(script::String, code::String, Ts::Float...)
        surveydir = "jobs/$(code)_$(script)"
        !isdir(surveydir) && mkdir(surveydir)

        for T in Ts
            outdir = "$surveydir/T$(T)"
            try init(outdir, script, code, T)
            catch ERROR println(ERROR.msg) end
        end
    end

    """ Duplicate an existing job, but modify the seed. """
    function init_shotgun(jobdir, seeds...)
        # LOAD TEMPLATE JOB
        if !ispath(jobdir) || !isfile("$jobdir/script.jl")
            # Also requires serialized setup and meta but it pry does if it has script.jl.
            error("Shotgunning a job requires a pre-initialized job.")
        end
        vars = load(jobdir; trace=false, state=false)

        for seed in seeds
            vars.outdir = "$jobdir/seed_$seed"
            if ispath(vars.outdir)
                println("Data already exists in $(vars.outdir).")
                continue
            end

            # ASSIGN THE SEED
            vars.setup.seed = seed

            # INITIALIZE THE NEW JOB DIR
            save(vars)
            cp("$jobdir/script.jl", "$(vars.outdir)/script.jl")
        end
    end

    """ Construct a list of genuine "job" directories within a path. """
    function eachjob(path="jobs")
        # Try to interpret `path` as an iterable of paths.
        pathlist = collect(path)
        eltype(pathlist) == String && return vcat((eachjob(item) for item in pathlist)...)

        # Interpret `path` as a single directory.
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
            flag = adapt_is_converged(vars) ? "C   " :
                adapt_is_terminated(vars) ? "T *  " :
                isfile("$(vars.outdir)/running") ? "R ** " :
                "- *** "

            entry = "$flag\t$(vars.outdir)"

            if !isnothing(vars.trace) && !isempty(vars.trace.iterations)
                I = last(vars.trace.iterations)
                E = last(vars.trace.energy)
                fn = last(vars.trace.fn)
                gd = last(vars.trace.gd)
                entry *= "\t\tI=$I\tE=$E\tfn=$fn\tgd=$gd"
            end

            if !isnothing(vars.trace) && !isempty(vars.trace.adaptations)
                a = length(vars.trace.adaptations)
                G = last(vars.trace.G_max)
                entry *= "\t\ta=$a\tG=$G"
            end

            println(entry)
        end
    end

    function unterminated_jobs(path="jobs")
        return filter(
            job -> !adapt_is_terminated(load(job)),
            eachjob(path),
        )
    end

    function inspect_all(path="jobs")
        for job in eachjob(path)
            inspect(load(job))
        end
    end

    function start(path="jobs")
        for job in unterminated_jobs(path)
            try run(`./start $job`) catch end
        end
    end

    function runlocal(vars=_!)
        include("$(vars.outdir)/script.jl")
    end

end # module AdaptiveModes_041724
