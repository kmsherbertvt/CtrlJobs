#= Try and plot out entanglement over final trajectories in ancilla runs. =#

##########################################################################################
#= PLOT ENTANGLEMENT ENTROPY AND PURITY OVER OPTIMIZED TRAJECTORY =#

#=

Stuff similar to RenyiADAPT.

=#

import TensorOperations: tensortrace
import LinearAlgebra: I, kron, tr, eigen

"""
    partial_trace(ρ::DensityMatrix, nH::Int)

--- Implementation courtesy of Jim Furches, in RenyiADAPT.jl. ---

Computes the partial trace of ρ, which is assumed to act on
subsystems V ⊗ H, where V is the visible system and H is the
hidden system of nH qubits.

# Parameters
- `ρ`: The full system whose subspace will be traced out
- `nH`: The number of hidden qubits to remove

# Returns
- The partial trace of ρ, removing the last `nH` qubits
"""
function partial_trace(ρ::Matrix, nH::Int)
    nH >= 0 || throw(ArgumentError("Hidden qubits must be positive"))

    # Treat the system as V ⊗ H, and trace out H
    qubits = trunc(Int, log2(size(ρ, 1)))
    nV = qubits - nH

    if nH == 0
        return ρ
    elseif nV == 0
        return tr(ρ)
    end

    # Fixme: I'm not sure why the subsystems are in this order. I thought
    # I should be able to reshape it to (V, H, V, H).
    ρ = reshape(ρ, (2^nH, 2^nV, 2^nH, 2^nV))
    return tensortrace(ρ, [1, -1, 1, -2])
end

realifclose(x::Complex) = (isapprox(imag(x), 0, atol=1e-9)) ? real(x) : x
purity(ρ::Matrix) = min(ρ^2 |> tr |> realifclose, 1)
von_neumann_entropy(ρ::Matrix) = (
    ΛU = eigen(ρ);
    Λ = ΛU.values;
    S = 0.0;
    for i in eachindex(Λ);
        (abs(Λ[i]) < eps(Float64) && continue);
        (real(Λ[i]) < eps(Float64) && continue);
            # Sometimes a slightly-larger-than-eps negative value leaks in.
        S += -Λ[i] * log(Λ[i]);
    end;
    S
)


#=

Stuff similar to src/__calculations.

=#

function calculate_projections(vars)
    JOB.require_work(vars)
    work = vars.work
    energyfn = work.lossfn.energyfn

    # CONSTRUCT THE PROJECTOR
    nH = CtrlVQE.QubitOperators.nqubits(energyfn.ψ0) - work.system.n
    In = Matrix{Bool}(I, 1<<work.system.n, 1<<work.system.n)
    # π̄H = zeros(Bool, 2, 2, nH); π̄H[1,1,:] .= 1
    # ΠH = CtrlVQE.LinearAlgebraTools.kron(π̄H)
    ΠH = zeros(Bool, 1<<nH, 1<<nH); ΠH[1,1] = 1
    Π = kron(In, ΠH)

    # π̄ = QubitOperators.localqubitprojectors(energyfn.device)
    # #= --- SPECIAL TO ANCILLAE --- =#
    # n = size(π̄,3); for q in 1+n-nH:n; π̄[2,2,q] = 0; end
    #     # Qubit projectors onto each ancilla will use |0⟩⟨0| instead of |0⟩⟨0|+|1⟩⟨1|
    # Π = CtrlVQE.LinearAlgebraTools.kron(π̄)

    # CONSTRUCT THE PROJECTED ENERGY OPERATOR
    O = copy(energyfn.O0)
    CtrlVQE.Devices.evolve!(
        energyfn.frame,
        energyfn.device,
        energyfn.basis,
        vars.setup.T,
        O,
    )
    ΠOΠ = CtrlVQE.LinearAlgebraTools.rotate!(Π, O)

    # RUN FULL TIME EVOLUTION - TODO: Generalize for other cost functions?
    CtrlVQE.Parameters.bind(energyfn.device, vars.state.x)
    ψ = CtrlVQE.Evolutions.evolve(
        energyfn.evolution,
        energyfn.device,
        energyfn.basis,
        energyfn.grid,
        energyfn.ψ0;
    )

    # MEASURE EXPECTATION VALUES
    N = real(CtrlVQE.LinearAlgebraTools.expectation(Π, ψ))
    EΠ = real(CtrlVQE.LinearAlgebraTools.expectation(ΠOΠ, ψ))

    return N, EΠ
end

#= TODO: We'll want another trajectory, tracking eigenpairs of ρS throughout. o_O =#

function trajectory_callback(vars, P, S)
    JOB.require_work(vars)
    work = vars.work
    energyfn = work.lossfn.energyfn

    workbasis = CtrlVQE.Evolutions.workbasis(energyfn.evolution)
    U = CtrlVQE.Devices.basisrotation(energyfn.basis, workbasis, work.device)
    ψ_ = similar(energyfn.ψ0)

    nH = CtrlVQE.QubitOperators.nqubits(energyfn.ψ0) - work.system.n

    return (i, t, ψ) -> (
        ψ_ .= ψ;
        CtrlVQE.LinearAlgebraTools.rotate!(U, ψ_);
            # ψ_ IS NOW IN MEASUREMENT BASIS
        CtrlVQE.Devices.evolve!(energyfn.frame, work.device, energyfn.basis, -t, ψ_);
            # ψ_ IS NOW IN ROTATING FRAME

        ρ = ψ_ * ψ_';               # TODO: Pre-allocate?
        ρS = partial_trace(ρ, nH);  # TODO: ...pre-allocate..? Seems unlikely...
        P[i] = purity(ρS);
        S[i] = von_neumann_entropy(ρS);
    )
end


function fill_entanglement_curves!(vars, P, S)
    callback = trajectory_callback(vars, P, S)
    CtrlVQE.cost_function(vars.work.lossfn.energyfn; callback=callback)(vars.state.x)
    return P, S
end


#=

Stuff similar to src/__plotting.

=#

import Plots

function plot_entanglement(vars, P, S; kwargs...)
    t = CtrlVQE.lattice(vars.work.grid)

    # INITIALIZE THE PLOT
    plot = Plots.plot()
    yMIN = 0.0
    yMAX = maximum(S)
    Plots.plot!(plot;
        framestyle = :box,
        xlabel = "Time (ns)",
        xlims  = JOB.Plotting.margined_lims(minimum(t), maximum(t); factor=0.00),
        xminorticks = 5,
        xminorgrid = true,
        ylabel = "Entanglement Entropy",
        ylims  =  JOB.Plotting.margined_lims(yMIN, yMAX),
        yminorgrid = true,
        yminorticks = 5,
        legend = :right,
        kwargs...
    )

    twin = Plots.twinx(plot)
    Plots.plot!(twin;
        ylabel = "Purity",
        ylims = JOB.Plotting.margined_lims(0.0, 1.0; factor=0.05),
        yminorticks = 5,
        legend = false,
    )

    # PLOT THE MAXIMUM ATTAINABLE ENTROPY
    # TODO: this. What is it, log(1<<nH) ?


    # PLOT THE ENERGY TRAJECTORY
    plotargs = (color=1, lw=3, ls=:solid)
    Plots.plot!(plot, t, S; label="Entropy", plotargs...)

    # PLOT THE LEAKAGE TRAJECTORY
    plotargs = (color=2, lw=3, ls=:solid)
    Plots.plot!(plot, [0], [2yMAX]; label="Purity", plotargs...)
    Plots.plot!(twin, t, P; plotargs...)

    return plot
end



#=

Stuff similar to __plot__surveys.

=#

import Plots
import ColorSchemes
using LaTeXStrings

default_series = (
    label = false,
    color = :black,
    linestyle = :solid,
    linewidth = 2,
    shape = :circle,
    marksersize = 5,
)

""" kwargs for plot initialization which sets the T-axis from TMIN to TMAX,
        with minor gridlines every Δs ns.
    Major ticks are selected to have approximately no more than ten labels total.
"""
function __initplot__T_axis(TMIN, TMAX, Δs)
    labelskip = 1+ round(Int, (TMAX - TMIN) / Δs) ÷ 10
    return (
        xlims = [TMIN,TMAX],
        xticks = TMIN:labelskip*Δs:TMAX,
        xminorticks = labelskip,
        xminorgrid = true,
    )
end

__initplot__ENORM(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Normalized Energy Error (Ha)",
    yscale = :log,
    ylims  = [1e-16, 1e2],
    yticks = 10.0 .^ (-16:2:0),
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__SEND(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Final Entropy",
    # yscale = :log,
    # ylims  = [1e-16, 1e2],
    # yticks = 10.0 .^ (-16:2:0),
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

__initplot__SMAX(; kwargs...) = Plots.plot(;
    xlabel = "Pulse Duration (ns)",
    ylabel = "Maximal Entropy",
    yminorticks = 1,
    yminorgrid = true,
    legend = :topright,
    tickfontsize = 11,  # Match the default label size.
    legendfontsize = 11,
    kwargs...
)

function init_plots(TMIN, TMAX, Δs; ENORM=Dict(), SEND=Dict(), SMAX=Dict(), kwargs...)
    T_axis = __initplot__T_axis(TMIN, TMAX, Δs)
    return (
        ENORM = __initplot__ENORM(; T_axis..., kwargs..., ENORM...),
        SEND = __initplot__SEND(; T_axis..., kwargs..., SEND...),
        SMAX = __initplot__SMAX(; T_axis..., kwargs..., SMAX...),
    )
end

function plot_survey!(plots, surveydir; energy=Dict(), iter=Dict(), kwargs...)
    # INITIALIZE RESULT VECTORS
    data = (
        T = JOB.Float[],        # Pulse duration
        N = JOB.Float[],        # Normalization ⟨ψ|Π|ψ⟩
        EΠ = JOB.Float[],       # Projected energy ⟨ψ|ΠHΠ|ψ⟩
        EN = JOB.Float[],       # Normalized energy error ⟨ψ|ΠHΠ|ψ⟩/⟨ψ|Π|ψ⟩ - FCI
        PEND = JOB.Float[],     # Purity at end of pulse.
        PMIN = JOB.Float[],     # Lowest purity attained throughout pulse.
        SEND = JOB.Float[],     # Entanglement entropy at end of pulse.
        SMAX = JOB.Float[],     # Highest entropy attained throughout pulse.
    )

    for jobdir in readdir(surveydir, join=true)
        endswith(jobdir, ".DS_Store") && continue
        # LOAD VARIABLES AND (IF NECESSARY) THE SYSTEM
        JOB.load!(jobdir; run=false)
        isnothing(_!.trace) && continue
        try include("$(_!.outdir)/script.jl") catch end

        # CALCULATE ENERGY TRAJECTORY
        N, EΠ = calculate_projections(_!)

        # CALCULATE ENTANGLEMENT TRAJECTORY
        P = Array{JOB.Float}(undef, _!.setup.r+1)
        S = Array{JOB.Float}(undef, _!.setup.r+1)
        fill_entanglement_curves!(_!, P, S)

        # GENERATE ENTANGLEMENT PLOT
        plt = plot_entanglement(_!, P, S)
        Plots.savefig(plt, "$(_!.outdir)/final_entanglement.pdf")
        # TODO: Don't actually need this usually..?

        # REGISTER THE DATA
        push!(data.T, _!.setup.T)
        push!(data.N, N)
        push!(data.EΠ, EΠ)
        push!(data.EN, EΠ/N - _!.work.system.FCI)
        push!(data.PEND, last(P))
        push!(data.PMIN, minimum(P))
        push!(data.SEND, last(S))
        push!(data.SMAX, maximum(S))
    end

    # SORT OUTPUTS BY TIME
    σ = sortperm(data.T)
    for array in data; permute!(array, σ); end

    # PLOT DATA
    Plots.plot!(plots.ENORM, data.T, data.EN;
            default_series..., kwargs..., energy...)
    Plots.plot!(plots.SEND, data.T, data.SEND;
            default_series..., kwargs..., energy...)
    Plots.plot!(plots.SMAX, data.T, data.SMAX;
            default_series..., kwargs..., iter...)
end

function save_plots(plots, label)
    for (kind, plot) in pairs(plots)
        Plots.savefig(plot, "figs/$(label)__$(kind).pdf")
    end
end








# #= Test functions with one trial. =#
# load!("jobs/H215_ancilla.single_ΔsMAX0.5/T10.0_W20"; run=false)
# try include("$(_!.outdir)/script.jl") catch end

# P = Array{JOB.Float}(undef, _!.setup.r+1)
# S = Array{JOB.Float}(undef, _!.setup.r+1)
# fill_entanglement_curves!(_!, P, S)

# plt = plot_entanglement(_!, P, S)
# Plots.savefig(plt, "$(_!.outdir)/final_entanglement.pdf")


# # COMPARING A FEW STRATEGIES OF INTRODUCING ANCILLAE
# plots = init_plots(0.0, 21.0, 0.5)
# plot_survey!(plots, "jobs/H215_resonant_ΔsMAX0.5";
#     label="1xN", color=1, linestyle=:solid, shape=:circle)
# # plot_survey!(plots, "jobs/H215_ancilla.single_ΔsMAX0.5";
# #     label="1x(N+1)", color=2, linestyle=:dash, shape=:utriangle)
# # plot_survey!(plots, "jobs/H215_ancilla.linear_ΔsMAX0.5";
# #     label="1x(2N)", color=3, linestyle=:dash, shape=:utriangle)
# # plot_survey!(plots, "jobs/H215_ancilla_ΔsMAX0.5";
# #     label="2xN", color=4, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle_ΔsMAX0.5";
#     label="N>o", color=5, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle.projected_ΔsMAX0.5";
#     label="N>o Π", color=6, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/H215_ancilla.fullsingle.projected.normalized_ΔsMAX0.5";
#     label="N>o Π/N", color=7, linestyle=:dash, shape=:utriangle)
# save_plots(plots, "ancillae")

# COMPARING A FEW STRATEGIES OF INTRODUCING ANCILLAE - LiH
plots = init_plots(0.0, 36.0, 0.5)
plot_survey!(plots, "jobs/lih30_resonant_ΔsMAX0.5";
    label="1xN", color=1, linestyle=:solid, shape=:circle)
# plot_survey!(plots, "jobs/lih30_ancilla.single_ΔsMAX0.5";
#     label="1x(N+1)", color=2, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_ancilla.linear_ΔsMAX0.5";
#     label="1x(2N)", color=3, linestyle=:dash, shape=:utriangle)
# plot_survey!(plots, "jobs/lih30_ancilla_ΔsMAX0.5";
#     label="2xN", color=4, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle_ΔsMAX0.5";
    label="N>o", color=5, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle.projected_ΔsMAX0.5";
    label="N>o Π", color=6, linestyle=:dash, shape=:utriangle)
plot_survey!(plots, "jobs/lih30_ancilla.fullsingle.projected.normalized_ΔsMAX0.5";
    label="N>o Π/N", color=7, linestyle=:dash, shape=:utriangle)
save_plots(plots, "ancillae_lih30")