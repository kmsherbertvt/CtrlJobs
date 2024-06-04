#= Casual observation suggests that `gradients.modalharmonics` and `gradients.exact`
    are really quite similar in performance.

Really, `gradients.exact` is the "continuous basis" limit of `gradients.modalharmonics`
    though even with infinite modes, the basis functions in `gradients.modalharmonics`
    are not in fact complete, so that's unfortunate.
But, that isn't so important, because the observation is that the default does WELL.
So the question is not "what's it take to converge?"
    (which mightn't be answerable, since the basis isn't complete)
    but, "what's it take to DIverge?"

That is, let's shotgun survey `gradients.modalharmonics` with smaller and smaller bases,
    until it *doesn't* match `gradients.exact`.

=#

function init_shotgun_nMAXs(jobdir, nMAXs...)
    # LOAD TEMPLATE JOB
    if !ispath(jobdir) || !isfile("$jobdir/script.jl")
        # Also requires serialized setup and meta but it pry does if it has script.jl.
        error("Shotgunning a job requires a pre-initialized job.")
    end
    vars = load(jobdir; trace=false, state=false)

    # Default to incrementing the number of modes from 1.
    if isempty(nMAXs)
        nMAX = floor(Int, vars.setup.fMAX * 2 * vars.setup.T)
        nMAXs = 1:nMAX
            # NOTE: For consistent interface, re-generate the original trajectory.
            #   (wasteful but convenient)
            # TODO: We could actually just...copy everything from the original in, here?
            #   for now just try to remember to do that manually...
    end

    for nMAX in nMAXs
        vars.outdir = "$jobdir/nMAX_$nMAX"
        if ispath(vars.outdir)
            println("Data already exists in $(vars.outdir).")
            continue
        end

        # ASSIGN THE SEED
        vars.setup.fMAX = nextfloat(nMAX / (2 * vars.setup.T))
            # NOTE: `nextfloat` evades round off errors when nMAX is calculated later on.

        # INITIALIZE THE NEW JOB DIR
        save(vars)
        cp("$jobdir/script.jl", "$(vars.outdir)/script.jl")
    end
end


""" Compare energy error vs parameters for a exact and incrementally complete modal.

We'll use black for "exact", and colors for each n.

"""
function plot_εvP(exactdir, modaldir; nMAXs=nothing)
    # Default to all nMAXs available.
    if isnothing(nMAXs)
        nMAXs = Int[]
        for jobdir in readdir(modaldir)
            re = match(r"nMAX_(\d+)", jobdir)
            isnothing(re) && continue
            push!(nMAXs, parse(Int, re[1]))
        end
        sort!(nMAXs)
    end

    plt = Plots.plot(;
        xlabel = "Parameter Count",
        ylabel = "Energy Error (Ha)",
        yscale = :log,
        ylims  = [1e-16, 1e2],
        yticks = 10.0 .^ (-16:2:0),
        yminorgrid = false,
        legend = :bottomleft,
    )
    Plots.hline!(plt, [0.00169]; color=:black, ls=:dot, label=false)
    Plots.vline!(plt, [30]; color=:black, ls=:dot, label=false)
    Plots.plot!(plt; xlims=[15,32], xticks=15:32)

    #################################
    #= Plot each limited modes data. =#
    vars = load(modaldir)
    FCI = JOB.System(vars.setup.code).FCI

    for nMAX in nMAXs
        jobdir = "$modaldir/nMAX_$nMAX"
        vars = load(jobdir)

        P = eachindex(vars.trace.adaptations) .- 1
        ε = vars.trace.energy[vars.trace.adaptations] .- FCI

        Plots.plot!(plt, P, ε;
            color = nMAX, lw=3, la=0.8,
            label = "n ≤ $nMAX",
        )
    end

    #################################
    #= Plot exact gradients data. =#
    vars = load(exactdir)
    FCI = JOB.System(vars.setup.code).FCI

    P = eachindex(vars.trace.adaptations) .- 1
    ε = vars.trace.energy[vars.trace.adaptations] .- FCI

    Plots.plot!(plt, P, ε;
        color = :black, lw=3, la=0.8,
        label = "Continuous",
    )

    return plt
end