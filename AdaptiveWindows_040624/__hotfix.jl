#=    -- Wizard's First Rule --

To use this file:
- You must be running in a ./juliarepl REPL, ie. the active project package is `JOB`.
- Include the file: `include("__hotfix.jl");`
- Check for errors: `Hotfix.scanfor_issues()`
- Fix fixable errors: `Hotfix.fixall_issues()`

=#

module Hotfix
    import ..JOB

    """ Sometimes the pulse gets kicked a bit too far into the forbidden region,
        and the cost function explodes.

    This is a real problem. I don't know how it can get kicked that much.
    For now, just terminate this trajectory; it got as far as it could.

    """
    function hasproblem_infinity(vars=_!)
        return last(vars.trace.fn) > 1e100
    end

    function hotfix_infinity!(vars=_!)
        vars.meta.maxadapt = 0
        JOB.save(vars)
        println("Modified $(vars.outdir)")
    end


    """ My `maxiter` attribute was meant to be a limit
        on iterations within a single optimization,
        but I used it for checks on the overall trajectory too. Oops.

    This is not a real problem. Just increase `_!.setup.maxiter`.

    Silly me, I'm just now realizing this attribute actually ought always
        have been in `meta`, not in `setup`. No big deal.

    """
    function hasproblem_iterations(vars=_!)
        maxiter = last(vars.trace.iterations)
        return maxiter ≥ vars.setup.maxiter
    end

    function hotfix_iterations!(vars=_!)
        maxiter = last(vars.trace.iterations)
        vars.setup.maxiter = 2 * maxiter
        JOB.save(vars)
        println("Modified $(vars.outdir)")
    end



    """ My `fnRATIO` attribute is meant to be a stopgap against
        not-quite-correct analytical gradients.
    But apparently sometimes the linesearch really does just take a lot.
    This actually is quite normal in gate-based ADAPT.

    This is not a real problem. Just increase `_!.setup.fnRATIO`.

    Silly me, I'm just now realizing this attribute actually ought always
        have been in `meta`, not in `setup`. No big deal.

    """
    function hasproblem_fnRATIO(vars=_!)
        fnRATIO = maximum(vars.trace.f_calls) / last(vars.trace.iterations)
        return fnRATIO ≥ vars.setup.fnRATIO
    end

    function hotfix_fnRATIO!(vars=_!)
        fnRATIO = maximum(vars.trace.f_calls) / last(vars.trace.iterations)
        vars.setup.fnRATIO = 2 * fnRATIO
        JOB.save(vars)
        println("Modified $(vars.outdir)")
    end


    """ Try to find a reason each job in the path is not terminated.

    If NONE reasons were found, it probably means your job just needs to be run.

    """
    function scanfor_issues(path="jobs")
        println("""

        I     \t = Optimization's just takin' awhile. Easy fix.
        fn    \t = Linesearches are just takin' awhile. Easy fix.
        ∞ *   \t = Penalty term exploded.   (the "fix" will just disable this run)
        - *** \t = NONE (so if there IS a problem, more investigation is needed)

        """)

        for job in JOB.unterminated_jobs(path)
            isfile("$job/running") && continue
            vars = JOB.load(job)

            flag = hasproblem_infinity(vars) ? "∞ *   " :
                hasproblem_iterations(vars) ? "I     " :
                hasproblem_fnRATIO(vars) ? "fn    " :
                "- *** "

            entry = "$flag\t$(vars.outdir)"
            println(entry)
        end
    end


    """ Modify config files in each job to fix a known problem.

    Be careful.

    """
    function fixall_issues!(path="jobs")
        for job in JOB.unterminated_jobs(path)
            isfile("$job/running") && continue
            vars = JOB.load(job)

            if hasproblem_infinity(vars)
                hotfix_infinity!(vars)
            elseif hasproblem_iterations(vars)
                hotfix_iterations!(vars)
            elseif hasproblem_fnRATIO(vars)
                hotfix_fnRATIO!(vars)
            else
                println("No fix available for $job (so maybe it's not broken?)")
            end
        end
    end

end