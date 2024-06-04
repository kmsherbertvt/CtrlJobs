#= In principle, `harmonics.complex.iterative` and `harmonics.complex.iterative.cold`,
    if the optimizations are "successful" (ie. find a global minimum)
    will produce energies which are identical.
The only differences between them should be the actual optimization trajectories,
    and the main question is, which is faster?
- Of course if you want the whole adaptive trajectory,
    generating a convergence plot as you add one parameter at a time,
    warm start will *obviously* be faster.
- If you only care about the result with a single set of parameters,
    then the single optimization from `cold` (independently runnable) may be faster.

I'd really like to know if the warm start ever *is* faster.
I suspect the picture may change as once changes g_tol (BFGS convergence),
    since a great deal of the cost of multiple optimizations is in a long tail.
A less strict gtol may favor warm start.
    Maybe. I mean, cold start gets a shorter tail too, so it might just be the same.
However, a less strict gtol will also destabilize ADAPT, and that's worth studying too.

This script provides functions for
- generating surveys of smaller g_tol within an existing job
- filling out the root job directory with plots that compare them.

It should be run from within a ./juliarepl.

=#

function init_shotgun_g_tol(jobdir, g_tols...)
    # LOAD TEMPLATE JOB
    if !ispath(jobdir) || !isfile("$jobdir/script.jl")
        # Also requires serialized setup and meta but it pry does if it has script.jl.
        error("Shotgunning a job requires a pre-initialized job.")
    end
    vars = load(jobdir; trace=false, state=false)

    for g_tol in g_tols
        vars.outdir = "$jobdir/g_tol_$g_tol"
        if ispath(vars.outdir)
            println("Data already exists in $(vars.outdir).")
            continue
        end

        # ASSIGN THE SEED
        vars.setup.g_tol = g_tol

        # INITIALIZE THE NEW JOB DIR
        save(vars)
        cp("$jobdir/script.jl", "$(vars.outdir)/script.jl")
    end
end

""" Compare iterations vs parameters for a warm-start and corresponding cold start job.

We'll use linestyle to designate "warm start", "cold start", and "diff cold start".
We'll extract all the g_tol jobs available, and distinguish them by color.

"""
function plot_IvP(warmdir, colddir; g_tols=nothing)
    color(g_tol) = round(Int, abs(log10(g_tol)))

    plt = Plots.plot(;
        xlabel = "Parameters",
        xlims = [0,35],
        ylabel = "BFGS Iterations",
        ylims = [0,400],
    )
    Plots.plot!(plt, [0],[0]; color=:black, lw=3, ls=:solid, label="Warm Start")
    Plots.plot!(plt, [0],[0]; color=:black, lw=3, ls=:dash, label="Cold Start")
    Plots.plot!(plt, [0],[0]; color=:black, lw=3, ls=:dot, label="Cold Start (Independent)")

    #################################
    #= Plot warm start data. =#
    for jobdir in [warmdir; readdir(warmdir; join=true)]
        jobdir == warmdir || startswith(basename(jobdir), "g_tol_") || continue
        vars = load(jobdir)
        isnothing(g_tols) || vars.setup.g_tol in g_tols || continue

        P = [sum(vars.trace.parameters[a]) for a in eachindex(vars.trace.adaptations)]
        I = vars.trace.adaptations

        Plots.plot!(plt, P, I;
            ls = :solid, lw=3, la=0.8,
            color = color(vars.setup.g_tol),
            label = "g_tol=1e-$(color(vars.setup.g_tol))",
        )

    end

    #################################
    #= Plot cold start data. =#
    for jobdir in [colddir; readdir(colddir; join=true)]
        jobdir == colddir || startswith(basename(jobdir), "g_tol_") || continue
        vars = load(jobdir)
        isnothing(g_tols) || vars.setup.g_tol in g_tols || continue

        P = [sum(vars.trace.parameters[a]) for a in eachindex(vars.trace.adaptations)]
        I = vars.trace.adaptations

        Plots.plot!(plt, P, I;
            ls = :dash, lw=3, la=0.8,
            color = color(vars.setup.g_tol),
            label = false,
        )

        Plots.plot!(plt, P, diff([0; I]);
            ls = :dot, lw=3, la=0.8,
            color = color(vars.setup.g_tol),
            label = false,
        )

    end

    return plt
end