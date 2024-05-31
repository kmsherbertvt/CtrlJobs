using AdaptFigures_052324

include("lih30_surveys.jl");
surveys = [
    lih30_surveys.sub_naive,
    lih30_surveys.sub_adapt,
    lih30_surveys.sub_optim,
]

prefix = "lih30_subdivision"

T_axis = (
    xlims = [0.0,50.0],
    xticks = 0.0:6.0:48.0,
    xminorticks = 2,
)

P_axis = (
    N = 30,         # Dimension of Hilbert space
    xlims = [0,60],
    xticks = 0:6:60,
    xminorticks = 3,
)

bMET_STYLE = (
    seriesalpha=0.4,
    linestyle = :dot,
)


##########################################################################################
#= Energy vs pulse duration. =#


plt = init_εvT(; T_axis...,
    legend = :topright,
)

plot!(plt, survey; kwargs...) = plot!_curve(
    plt, load_εvT(survey.job, "$(survey.surveydir)")...;
    label = survey.label,
    survey.style...,
    kwargs...
)

foreach(survey -> plot!(plt, survey), surveys)

savepdf(plt, prefix, "εvT")



##########################################################################################
#= Energy vs parameters. =#

plt = init_εvP(; P_axis...,
    legend = :bottomleft,
)

label!(plt, "At eMET (22-24 ns)"; )
label!(plt, "Beyond eMET (36 ns)"; bMET_STYLE...)

plot!(plt, survey; kwargs...) = (
    plot!_curve(
        plt, load_εvP(survey.job, "$(survey.surveydir)/T$(survey.eMET)")...;
        label = survey.label,
        survey.style...,
        kwargs...
    );
    plot!_curve(
        plt, load_εvP(survey.job, "$(survey.surveydir)/T36.0")...;
        label = false,
        survey.style...,
        bMET_STYLE...,
        kwargs...
    );
)

foreach(survey -> plot!(plt, survey), surveys)

savepdf(plt, prefix, "εvP")





##########################################################################################
#= Parameters vs pulse duration. =#

plt = init_PvT(; N=P_axis.N, T_axis...,
    legend = :topright,
)

plot!(plt, survey; kwargs...) = plot!_curve(
    plt, load_PvT(survey.job, "$(survey.surveydir)")...;
    label = survey.label,
    survey.style...,
    kwargs...
)

foreach(survey -> plot!(plt, survey), surveys)

savepdf(plt, prefix, "PvT")

##########################################################################################
#= Iterations vs pulse duration. =#

plt = init_IvT(; T_axis...,
    legend = :topright,
)

plot!(plt, survey; kwargs...) = plot!_curve(
    plt, load_IvT(survey.job, "$(survey.surveydir)")...;
    label = survey.label,
    survey.style...,
    kwargs...
)

foreach(survey -> plot!(plt, survey), surveys)

savepdf(plt, prefix, "IvT")


##########################################################################################
#= Iterations vs parameters. =#

plt = init_IvP(; P_axis...,
    legend = :topleft,
)

label!(plt, "At eMET (22-24 ns)"; )
label!(plt, "Beyond eMET (36 ns)"; bMET_STYLE...)

plot!(plt, survey; kwargs...) = (
    plot!_curve(
        plt, load_IvP(survey.job, "$(survey.surveydir)/T$(survey.eMET)")...;
        label = survey.label,
        survey.style...,
        kwargs...
    );
    plot!_curve(
        plt, load_IvP(survey.job, "$(survey.surveydir)/T36.0")...;
        label = false,
        survey.style...,
        bMET_STYLE...,
        kwargs...
    );
)

foreach(survey -> plot!(plt, survey), surveys)

savepdf(plt, prefix, "IvP")