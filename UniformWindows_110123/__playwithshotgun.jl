#= We have lotsa lotsa runs shotgunning the initial state.
    Let's try to plot some things.
=#

jobdirs = Dict(
    # :T17_s15 => "jobs/lih30_resonant_ΔsMAX1.5/T17.0_W11",
    # :T18_s15 => "jobs/lih30_resonant_ΔsMAX1.5/T18.0_W12",

    :T18_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T18.0_W6",
    :T19_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T19.0_W6",
    :T20_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T20.0_W6",

    # :T21_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T21.0_W7",
    # :T30_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T30.0_W10",
    # :T48_s3 => "jobs/lih30_resonant_ΔsMAX3.0/T48.0_W16",
)
ε0 = Dict()
ε = Dict()
Ni = Dict()
Nf = Dict()
dirs = Dict()

# LOAD THE ZERO-PULSE RUN AND FETCH ITS ENERGY
for (key, jobdir) in jobdirs
    load!(jobdir; run=false)
    try include("$jobdir/script.jl") catch end
    initstate = JOB.initial_state()
    xi0 = initstate.x
    xf0 = _!.state.x
    E0 = last(_!.trace.energy)
    FCI = _!.work.system.FCI


    # LOAD INITIAL PARAMETER VECTORS, INITIAL ENERGIES, AND FINAL ENERGIES FOR EACH PELLET
    xi = Vector{JOB.Float}[]
    xf = Vector{JOB.Float}[]
    E = JOB.Float[]
    pelletdirs = String[]
    for pelletdir in readdir(jobdir)
        startswith(pelletdir, "shotgun_") || continue
        push!(pelletdirs, pelletdir)
        vars = load("$jobdir/$pelletdir")
        vars.work = _!.work
        initstate = JOB.initial_state(vars)
        push!(xi, initstate.x)
        push!(xf, vars.state.x)
        push!(E, last(vars.trace.energy))
    end
    dirs[key] = pelletdirs

    # # CONVERT PARAMETER VECTORS INTO RATIOS OF ΩMAX (only makes sense for `resonant` job)
    # xn = reduce(hcat, v/_!.setup.ΩMAX for v in x)

    # COMPUTE NORMS FOR EACH TRIAL
    import LinearAlgebra: norm
    Ni[key] = [norm(xi[i]) for i in eachindex(xi)]
    Nf[key] = [norm(xf[i] .- xf0) for i in eachindex(xf)]

    # COMPUTE ENERGY ERRORS FOR EACH TRIAL
    ε0[key] = E0 - FCI
    ε[key] = E .- FCI
end

# # PLOT ENERGY AGAINST NORM OF INITIAL STATE ON A SCATTER PLOT
# import Plots
# plt = Plots.plot(;
#     xlabel = "Initial Distance from Zero-Pulse",
#     ylabel = "Final Energy Error",
#     yscale = :log10,
#     ylims = [1e-16,1e2],
#     yticks = 10.0 .^ (-16:2:0),
#     legend = :topleft,
# )
# for (i, key) in enumerate(keys(jobdirs))
#     Plots.scatter!(plt, Ni[key], ε[key]; color=i, label=string(key))
#     Plots.scatter!(plt, [0.0], [ε0[key]]; color=i, shape=:+, msw=2, label=false)
# end

# PLOT ENERGY AGAINST DISTANCE OF PULSE FROM OPTIMIZED ZERO-PULSE ON A SCATTER PLOT
import Plots
plt = Plots.plot(;
    xlabel = "Final Distance from Optimized Zero-Pulse Run",
    ylabel = "Final Energy Error",
    yscale = :log10,
    ylims = [1e-16,1e2],
    yticks = 10.0 .^ (-16:2:0),
    legend = :topleft,
)
for (i, key) in enumerate(keys(jobdirs))
    Plots.scatter!(plt, Nf[key], ε[key]; color=i, ma=0.5, label=string(key))
    Plots.scatter!(plt, [0.0], [ε0[key]]; color=i, ma=0.5, shape=:x, msw=4, label=false)
end

Plots.gui()

#= TODO: 20 trials with longer pulses: separate plot
    ...or maybe simply separate marker type? =#
#= TODO: Gotta do the 100 trials with continuous pulse also...
    but we'll need to not save Hessian in state!!!
    Maybe we hack in a "sparse save resonant" script that just shrinks the Hessian to a meaningless 0x0 matrix on job completion, and expands it into an identity matrix on job load. =#
