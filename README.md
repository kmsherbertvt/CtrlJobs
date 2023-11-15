# CtrlJobs
Facilitates transfering files between desktop and remote HPC resources

Principles:
- Low-resource job initialization that puts everything a job needs to run in a job folder.
- Simple interface for starting a job, and almost-as-easy interface to do so with slurm.
- Jobs are NOT saved to git. Transfer data files with scp.
- Packages provide tools to _inspect_ job data _without_ any expensive calculations.
- Packages should be usable in top-level scripts, eg. to plot data from different packages.

Preferrred Workflow:
- Package defines state variables in an object, which can be archived throughout a job.
- Package defines an easy `load` function to load any archived state as the "active" one.
- The package defines a global variable `!_` to hold the active state.
- Plotting and other analysis methods act directly on the active state.

I expect this business with the `!_` variable will look very odd to any visitors,
    but I've found it works extraordinarily well for me.

# Project Naming Scheme

Something semantic, followed by something indicating the date it was first created.
This insures against having to come up with new semantic names when I want to overhaul something.
This also allows each project to set a particular "commit" of `CtrlVQE.jl` as its dependency,
    insuring against overhauls in that package without ruining my ability to quickly run new calculations.