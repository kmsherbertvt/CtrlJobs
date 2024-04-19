
module BookKeeping
    import Serialization: serialize, deserialize

    import ..Variables: Vars, SetupVars, MetaVars, TraceVars, StateVars, WorkVars

    ######################################################################################
    #= Safe serialization of custom structs.

    (I could just serialize all my `Vars` directly,
        but then the data is lost if I ever change the struct structure.
    The "right" solution is to specify a version for the package,
        but...that is too much work.
    So instead I serialize them as dicts.
    Then, if something in the struct changes,
        I can still load previous data and (manually) compare or adapt it to the new form.)

    =#

    function as_dict(obj)
        return Dict(fld => getfield(obj, fld) for fld in fieldnames(typeof(obj)))
    end

    function deserialize_dict(file, type)
        kwargs = deserialize(file)
        isempty(kwargs) && return nothing
        return type(; kwargs...)
    end

    ######################################################################################
    #= Interface for serializing/deserializing variables to files. =#

    """ Replace each attribute of the master variable object with its default. """
    function unload!(vars)
        fresh = Vars()
        for field in fieldnames(Vars)
            setfield!(vars, field, getfield(fresh, field))
        end
        return vars
    end

    """ Load each serializable variable from its standard name in the given outdir.

    # Keyword Arguments
    - `newdir`: the master variable's new `outdir` value
    - `run`: the master variable's new `run` value
    - `meta`: whether to load the serialized `meta` variables
    - `trace`: whether to load the serialized `trace` variables
    - `state`: whether to load the serialized `state` variables

    """
    function load!(vars, outdir; newdir=outdir, run=true, meta=true, trace=true, state=true)
        vars.outdir = newdir
        vars.setup = deserialize_dict("$outdir/setup", SetupVars)
        meta  && (vars.meta  = deserialize_dict("$outdir/meta", MetaVars))
        trace && (vars.trace = deserialize_dict("$outdir/trace", TraceVars))
        state && (vars.state = deserialize_dict("$outdir/state", StateVars))
        vars.work = nothing
        vars.run = run
        return vars
    end

    load(outdir; kwargs...) = load!(Vars(), outdir; kwargs...)

    """ Save each serializable variable to a standard name in the assigned outdir. """
    function save(vars)
        !isdir(vars.outdir) && mkdir(vars.outdir)
        serialize("$(vars.outdir)/setup", as_dict(vars.setup))
        serialize("$(vars.outdir)/meta",  as_dict(vars.meta))
        serialize("$(vars.outdir)/trace", as_dict(vars.trace))
        serialize("$(vars.outdir)/state", as_dict(vars.state))
    end

    """ Save the state variable object to a given name in the assigned outdir. """
    function archive(vars, name)
        !isdir(vars.outdir) && mkdir(vars.outdir)
        serialize("$(vars.outdir)/$name", as_dict(vars.state))
    end

    """ Load the state variable object from a given name in the assigned outdir. """
    function unarchive!(vars, name)
        vars.state = deserialize_dict("$(vars.outdir)/$name", StateVars)
    end

    ######################################################################################
    #= Standardized file names for archiving states. =#

    """ Standard string for an iteration number.

    Width is fixed to 8, so that the iterations order lexicographically and ordinally,
        so long as iterations do not exceed a hundred-million...

    """
    iterid(iter) = lpad(iter, 8, "0")

    """ Standard string for an adaptation number.

    Width is fixed to 4, so that the iterations order lexicographically and ordinally,
        so long as iterations do not exceed ten thousand...

    """
    adaptid(iter) = "a"*lpad(iter, 4, "0")

end