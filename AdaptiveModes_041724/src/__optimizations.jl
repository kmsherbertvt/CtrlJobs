module Optimizations
    import Optim, Printf

    import CtrlVQE

    import ..Float
    import ..require_work
    import ..save, ..archive, ..iterid
    import ..report
    import ..plot_trace, ..plot_pulse


    function make_objective(vars)
        require_work(vars)
        f = CtrlVQE.cost_function(vars.work.lossfn)
        g! = CtrlVQE.grad_function_inplace(vars.work.lossfn)
        return Optim.OnceDifferentiable(f, g!, vars.state.x)
    end

    function make_state(vars, optimizer, objective)
        return Optim.initial_state(
            optimizer,
            Optim.Options(),    # Not actually used, at least in BFGS interface.
            objective,
            vars.state.x,
        )
        bfgs_state.invH .= vars.state.Hk
    end

    function make_callback(vars, state)
        trace = vars.trace
        iteration = isempty(trace.iterations) ? Ref(0) : Ref(last(trace.iterations))

        return (cb_state) -> (
            iteration[] += 1;

            # PRINT ITERATION LINE
            println(Printf.format(
                Printf.Format("%8s    " * "%13.8g    "^3),
                cb_state.iteration,
                cb_state.value,
                cb_state.g_norm,
                cb_state.metadata["time"],
            ));

            # UPDATE TRACES
            update_trace!(vars, iteration[], cb_state.value, cb_state.g_norm);

            # UPDATE DATA
            vars.state.x  .= state.x;
            vars.state.Hk .= state.invH;
            save(vars);

            # MORE EXTENSIVE UPDATES EVERY SO OFTEN
            if iteration[] % vars.meta.update == 0;
                report(vars);
                plot_trace(vars, "");
                plot_pulse(vars, "", :amplitudes);
                archive(vars, iterid(iteration[]));
            end;

            optimization_is_terminated(vars)
        )
    end

    function make_options(vars, callback)
        return Optim.Options(
            f_tol = vars.setup.f_tol,
            g_tol = vars.setup.g_tol,
            iterations = vars.setup.maxiter,
            callback = callback,
        )
    end

    """ Update trace variables in the middle of an optimization.

    Scripts should call adapt_trace! when optimization completes.

    """
    function update_trace!(vars, iteration, fn, gd)
        require_work(vars)
        trace = vars.trace
        work  = vars.work

        push!(trace.iterations, iteration)
        push!(trace.f_calls, work.lossfn.f_counter[])
        push!(trace.g_calls, work.lossfn.g_counter[])

        push!(trace.fn, fn)
        push!(trace.gd, gd)

        push!(trace.energy, work.lossfn.energy[])
        push!(trace.penalties, copy(work.lossfn.penalties))

        nD = CtrlVQE.ndrives(work.device)
        Ωmax = Vector{Float}(undef, nD)
        ν    = Vector{Float}(undef, nD)
        for i in 1:nD
            # CALCULATE LARGEST Ω OVER ALL WINDOWS IN PULSE
            pulse = CtrlVQE.drivesignal(work.device, i)
            # TODO: Who signed off on this?! Needs to be type-agnostic.
            # Ωmax[i] = maximum(s -> abs(pulse(s)), pulse.starttimes)
            # TODO: This is a hack solution but we should use work variables more wisely.
            Ωmax[i] = maximum(abs.(CtrlVQE.valueat(pulse, CtrlVQE.lattice(work.grid))))

            # FETCH DRIVE FREQUENCY OF PULSE
            ν[i] = CtrlVQE.drivefrequency(vars.work.device, i)
        end
        push!(trace.Ωmax, Ωmax)
        push!(trace.ν, ν)
    end



    """ Inspect last record of the trace to see if it meets the convergence conditions. """
    function optimization_is_converged(vars)
        isnothing(vars.trace) && return false

        # SINGLE-ITERATION CHECKS
        isempty(vars.trace.iterations) && return false
        last(vars.trace.gd) ≤ vars.setup.g_tol && return true

        # TWO-ITERATION CHECKS
        length(vars.trace.iterations) < 2 && return false
        abs(vars.trace.fn[end-1] - vars.trace.fn[end]) ≤ vars.setup.f_tol && return true

        return false
    end

    """ Inspect last record of the trace to see if it meets the termination conditions. """
    function optimization_is_terminated(vars)
        isnothing(vars.trace) && return false
        isempty(vars.trace.iterations) && return false

        iterations = last(vars.trace.iterations)
        f_calls = last(vars.trace.f_calls)

        return any((
            optimization_is_converged(vars),
            # TOO MANY ITERATIONS - You can just run again from the final state.
            iterations >= vars.setup.maxiter,
            # TOO MANY FUNCTION EVALUATIONS - The linesearch seems to be unstable.
            iterations >= 10 && (f_calls > vars.setup.fnRATIO * iterations),
        ))
    end

end