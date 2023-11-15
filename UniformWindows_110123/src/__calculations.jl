module Calculations
    import LinearAlgebra

    import CtrlVQE
    import CtrlVQE: EvolutionType, CostFunctionType, DeviceType

    import ..Float, ..require_work

    ######################################################################################
    #= FUNCTIONS OF TIME =#

    """ Calculate E(t) = ⟨Ψ|H|Ψ⟩ throughout the pulse duration. """
    function make_energytrajectory(vars)
        require_work(vars)

        E = Array{Float}(undef, vars.setup.r+1) # Stores energy trajectory.

        callback = CtrlVQE.trajectory_callback(vars.work.lossfn.energyfn, E)
        CtrlVQE.cost_function(vars.work.lossfn.energyfn; callback=callback)(vars.state.x)

        return E
    end

    """ Calculate N(t) = ⟨Ψ|Π|Ψ⟩ throughout the pulse duration. """
    function make_normtrajectory(vars)
        require_work(vars)

        F = Array{Float}(undef, vars.setup.r+1) # Stores norm trajectory (for leakage).

        callback = CtrlVQE.trajectory_callback(vars.work.normfn, F)
        CtrlVQE.cost_function(vars.work.normfn; callback=callback)(vars.state.x)

        return F
    end

    """ Construct Ω(t) for each pulse, throughout pulse duration.

    TODO: Probably this guy only needs to make Ω, rather than α and β also.
        Only real reason to make them here is for the sake of minimizing memory,
            and we don't care about that in analysis generally.
        Then again, maybe it was for symmetry with gradient signal...so change that too. ;)

    """
    function make_pulses(vars)
        require_work(vars)

        CtrlVQE.Parameters.bind(vars.work.device, vars.state.x)
        t = CtrlVQE.lattice(vars.work.grid)

        nD = CtrlVQE.ndrives(vars.work.device)
        Ω = Array{Complex{Float}}(undef, vars.setup.r+1, nD)
        for i in 1:nD
            Ω[:,i] .= CtrlVQE.drivesignal(vars.work.device, i)(t)
        end

        return Ω, real.(Ω), imag.(Ω)
    end

    """ Construct ∂E/∂Ω(t) for each pulse, throughout pulse duration.

    TODO: Should we really be splitting ϕα and ϕβ here? Maybe better to have a separate function do that, yes?

    """
    function make_gradientsignals(vars)
        require_work(vars)

        nG = CtrlVQE.ngrades(vars.work.device)
        ϕ = Array{Float}(undef, vars.setup.r+1, nG)   # Stores gradient signals.
        CtrlVQE.grad_function(vars.work.lossfn.energyfn; ϕ=ϕ)(vars.state.x)
            # NOTE: Not compatible with NormalizedEnergy.

        nD = CtrlVQE.ndrives(vars.work.device)
        ϕα = Array{Float}(undef, vars.setup.r+1, nD)
        ϕβ = Array{Float}(undef, vars.setup.r+1, nD)
        for i in 1:nD
            j = 2(i-1) + 1
            ϕα[:,i] .= ϕ[:,j  ]
            ϕβ[:,i] .= ϕ[:,j+1]
        end

        return ϕ, ϕα, ϕβ
    end

    #= TODO: Might as well bring the FFT back, I think.
        Main motivation for removing it was clutter, which is now resolved. =#

end