
module Printing
    import Printf

    import CtrlVQE

    import ..Vars

    #= TODO (lo): Printing numbers with %s is dangerous, if float rep exceeds allocated length. But it's the only way to print numbers and headers with the same Format. Is there another way? =#

    ######################################################################################

    """ Describe variables as thoroughly as possible without doing any calculations. """
    function inspect(vars::Vars, flags::Symbol...)
        notepad = """

        Inspecting job variables
        ------------------------
        """
        notepad *= join((
            __setup__(vars.setup, flags),
             __meta__(vars.meta,  flags),
            __trace__(vars.trace, flags),
            __state__(vars.state, flags),
        ))

        println(notepad)
    end

    inspect(vars::Vars) = inspect(vars, :system, :optimization, :energies)

    function __fields__(object, header, fields)
        template = Printf.Format("%8s: %s\n")

        notepad = "\n\n***   $header   ***\n\n"
        for field in fields
            value = getfield(object, field)
            notepad *= Printf.format(template, field, value)
        end
        return notepad
    end

    function __setup__(setup, flags)
        notepad = ""
        :system in flags && (notepad *= __fields__(setup,
            "System Variables",
            (:code, :T, :W),
        ))
        :simulation in flags && (notepad *= __fields__(setup,
            "Simulation Variables",
            (:r, :m),
        ))
        :hardware in flags && (notepad *= __fields__(setup,
            "Hardware Bounds",
            (:ΩMAX, :ΔMAX),
        ))
        :init in flags && (notepad *= __fields__(setup,
            "Parameter Initialization",
            (:seed_Ω, :seed_φ, :seed_Δ, :kick_Ω, :kick_φ, :kick_Δ),
        ))
        :penalty in flags && (notepad *= __fields__(setup,
            "Penalty Variables",
            (:λΩ, :λΔ, :σΩ, :σΔ),
        ))
        return notepad
    end

    function __meta__(meta, flags)
        notepad = ""
        :meta in flags && (notepad *= __fields__(meta,
            "Meta-variables",
            (:g_tol, :maxiter, :fnRATIO, :update),
        ))
        return notepad
    end

    function __trace__(trace, flags; maxshow=5)
        isnothing(trace) && return "No trace is currently loaded."
        isempty(trace.iterations) && return "Trace is currently empty."

        nrecords = length(trace.iterations)
        nshow = min(maxshow, nrecords)

        notepad = ""
        :optimization in flags && (notepad *= __trace__optimization(trace, nrecords, nshow))
        :energies     in flags && (notepad *= __trace__energies(trace, nrecords, nshow))
        :amplitudes   in flags && (notepad *= __trace__amplitudes(trace, nrecords, nshow))
        :frequencies  in flags && (notepad *= __trace__frequencies(trace, nrecords, nshow))
        return notepad
    end

    function __trace__header(template, nrecords, nshow, args)
        notepad = Printf.format(template, args...)
        if nrecords > nshow
            # ADD ELLIPSES TO INDICATE THERE IS UN-REPORTED DATA
            notepad *= Printf.format(template, "⋮"^(length(args))...)
        end
        return notepad
    end

    function __trace__optimization(trace, nrecords, nshow)
        template = Printf.Format("%6s    "^3 * "%22.22s    "^2 * "\n")

        notepad = "\n\n***   Optimization Trace   ***\n\n"
        notepad *= __trace__header(
            template, nrecords, nshow,
            ("Iter", "# f", "# ∇f", "f(x)", "∇f(x)"),
        )
        for i in 1:nshow
            i += nrecords - nshow
            notepad *= Printf.format(template,
                trace.iterations[i],
                trace.f_calls[i],
                trace.g_calls[i],
                trace.fn[i],
                trace.gd[i],
            )
        end

        return notepad
    end

    function __trace__energies(trace, nrecords, nshow)
        nΛ = length(last(trace.penalties))
        template = Printf.Format("%22.22s    "^(1+nΛ) * "\n")

        notepad = "\n\n***   Energy and Penalty Functions (Ha)   ***\n\n"
        notepad *= __trace__header(
            template, nrecords, nshow,
            ("E", ("Λ$j" for j in 1:nΛ)...),
        )
        for i in 1:nshow
            i += nrecords - nshow
            notepad *= Printf.format(template, trace.energy[i], trace.penalties[i]...)
        end

        return notepad
    end

    function __trace__amplitudes(trace, nrecords, nshow)
        nD = length(first(trace.Ωmax))
        template = Printf.Format("%22.22s    "^nD * "\n")

        notepad = "\n\n***   Maximal Amplitudes (GHz)   ***\n\n"
        notepad *= __trace__header(
            template, nrecords, nshow,
            ("Pulse 1", ("$j" for j in 2:nD)...),
        )
        for i in 1:nshow
            i += nrecords - nshow
            notepad *= Printf.format(template, (trace.Ωmax[i] ./ 2π)...)
        end

        return notepad
    end

    function __trace__frequencies(trace, nrecords, nshow)
        nD = length(first(trace.Ωmax))
        template = Printf.Format("%22.22s    "^nD * "\n")

        notepad = "\n\n***   Drive Frequencies (GHz)   ***\n\n"
        notepad *= __trace__header(
            template, nrecords, nshow,
            ("Pulse 1", ("$j" for j in 2:nD)...),
        )
        for i in 1:nshow
            i += nrecords - nshow
            notepad *= Printf.format(template, (trace.ν[i] ./ 2π)...)
        end

        return notepad
    end

    function __state__(state, flags)
        isnothing(state) && return "No state is currently loaded."

        notepad = ""
        :counts  in flags && (notepad *= __state__counts(state))
        :arrays  in flags && (notepad *= __state__arrays(state))
        :hessian in flags && (notepad *= __state__hessian(state))
        return notepad
    end

    function __state__counts(state)
        template = Printf.Format("%8s: %s\n")
        notepad = "\n\n*** Parameter Counts ***\n\n"
        notepad *= Printf.format(template, "# Ω", length(state.Ω))
        notepad *= Printf.format(template, "# φ", length(state.φ))
        notepad *= Printf.format(template, "# ν", length(state.ν))
        return notepad
    end

    function __state__arrays(state)
        template = Printf.Format("%8s: %s\n")
        notepad = """

        All Parameters
        --------------
        """
        notepad *= Printf.format(template, "Ω", state.x[state.Ω])
        notepad *= Printf.format(template, "φ", state.x[state.φ])
        notepad *= Printf.format(template, "ν", state.x[state.ν])
        return notepad
    end

    function __state__hessian(state)
        notepad = """

        Inverse Hessian Approximation
        -----------------------------
        """
        io = IOBuffer()
        show(IOContext(io), "text/plain", state.Hk)
        notepad *= String(take!(io))
        return notepad
    end



    ######################################################################################

    """ Elaborate on a trace iteration using work variables or small calculation. """
    function report(vars::Vars, iter::Int, flags::Symbol...)
        isnothing(vars.work) && return  "Work variables are required to call `report`."
        isnothing(vars.trace) && return "No trace is currently loaded."
        length(vars.trace.iterations) < iter && return "Iteration $iter not in trace."

        notepad = """

        Reporting on iteration #$iter
        -----------------------------
        """
        :energy   in flags && (notepad *=   __report__energy(vars, iter))
        :cheating in flags && (notepad *= __report__cheating(vars, iter))
        :detuning in flags && (notepad *= __report__detuning(vars, iter))
        println(notepad)
    end

    report(vars::Vars, iter::Int) = report(vars, iter, :energy, :cheating, :detuning)

    report(vars::Vars) = isnothing(vars.trace) ?
            "No trace is currently loaded." :
            report(vars, lastindex(vars.trace.iterations))

    report(vars::Vars, flags::Symbol...) = isnothing(vars.trace) ?
            "No trace is currently loaded." :
            report(vars, lastindex(vars.trace.iterations), flags...)

    function __report__energy(vars, iter)
        E = vars.trace.energy[iter]
        system = vars.work.system
        model = system.model

        ERR = E          - system.FCI
        COR = system.REF - system.FCI
        GAP = model.Λ[2] - system.FCI

        template = Printf.Format("%13s: %22.22s\n")
        notepad = "\n\n***   Energy Report   ***\n\n"
        notepad *= Printf.format(template, "Energy (Ha)", E)
        notepad *= Printf.format(template, "Energy Error", ERR)
        notepad *= Printf.format(template, "% Corr Energy", 100*(1 - ERR/COR))
        notepad *= Printf.format(template, "%  Gap Energy", 100*(1 - ERR/GAP))
        return notepad
    end

    function __report__cheating(vars, iter)
        Ωmax = vars.trace.Ωmax[iter]
        χ = (Ωmax ./ vars.setup.ΩMAX) .- 1                  # CALCULATE CHEAT FACTOR
        for i in eachindex(χ); χ[i] ≤ 0 && (χ[i] = 0); end  # ZERO OUT "LEGAL" VALUES

        nD = length(Ωmax)
        template = Printf.Format("%22.22s    "^nD * "\n")
        notepad = "\n\n***   Cheat Factor (% over ΩMAX)   ***\n\n"
        notepad *= Printf.format(template, "Pulse 1", ("$j" for j in 2:nD)...)
        notepad *= Printf.format(template, (100 .* χ)...)   # SCALE TO PERCENT
        return notepad
    end

    function __report__detuning(vars, iter)
        ν = vars.trace.ν[iter]
        device = vars.work.device

        nD = length(ν)
        template = Printf.Format("%6.6s    " * "%22.22s    "^nD * "\n")
        notepad = "\n\n***   Detunings (GHz)   ***\n\n"
        notepad *= Printf.format(template, "Qubit", "Pulse 1", ("$j" for j in 2:nD)...)
        for q in 1:CtrlVQE.nqubits(device)
            ω = CtrlVQE.resonancefrequency(device, q)
            notepad *= Printf.format(template, q, ((ν .- ω) ./ 2π)...)
        end
        return notepad
    end

end