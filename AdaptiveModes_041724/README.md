# README

This package, and its siblings, are designed to streamline my own use of ARC.
It was designed to make sense to *me*, and to be robust against stupid errors I am likely to make.
It may not make sense to you, so here's an attempt at a tutorial.

-Kyle Sherbert

P.S. For context, this is addressed to Hisham Amer, but it should be readable by anyone.

## Reference Guide

A short-list of step-by-steps to do something useful.

A line beginning with `$` just indicates typing the subsequent command into a shell terminal, and with `>` into a Julia REPL.

### Terminal Management
- Navigate to the project directory:

  ```
  $ cd <bla bla>/AdaptiveModes_041724
  ```

- Make sure it's all up to date:

  ```
  $ git pull
  ```

- Start a Julia REPL:

  ```
  $ ./juliarepl
  ```

### Job Management
- Initialize a set of jobs:

  ```
  > JOB.init_survey("<script>", "<code>", <T1>, <T2>, ...)
  ```

- Check on all the jobs for a given script/system:

  ```
  > JOB.status("jobs/<code>_<script>")
  ```

- Start all the incomplete jobs for a given script/system:

  ```
  > JOB.start("jobs/<code>_<script>")
  ```

### Data Management
- Load a specific job into the active memory:

  ```
  > load!("jobs/<code>_<script>/T<T>")
  ```

- Check out the parameters currently optimized for this job:

  ```
  > _!.trace.x
  ```

- Initialize work variables without actually running the job

  ```
  > _!.run = false
  > try include("$(_!.outdir)/script.jl"); catch end
  ```



## Tutorial Guide

Some explanation of what and why the above reference lists what it lists:

### Terminal Management
#### Navigate to the project directory

First thing's first: you want a terminal with the project directory (`AdaptiveModes_041724/`)
    as your present working directory.
A little awkwardly (and I think I kinda regret it now),
    the *git* project directory is `CtrlJobs`, one level up.
My personal preference is to open `CtrlJobs` in VS Code,
    so that the integrated terminal starts in `CtrlJobs`,
    and then I run `cd AdaptiveModes_041724`.

#### Make sure it's all up to date

Now we want to make sure your local repository is up to date with everything on the Internet.
That's just typing in `git pull` in the terminal.

The catches:
1. Your local repository may be cloned from a "fork"?
   So if I've made changes to my repo, you may need to first update your own fork.
   My personal preference is to handle fork management through the GitHub interface rather than the terminal.
2. If you have changed files that would be updated by `git pull`,
   you will need to first either *commit* those changes, discard them, or "stash" them.
   I have no idea how stashing works.
3. If you chose to *commit* those changes,
   and you changed any of the same *lines of code* that would be updated by `git pull`,
   you will need to deal with a *merge*. Good luck.
   (Git will step you through it, for the most part.)


#### Start a Julia REPL

Everything in this package is meant to be done through the Julia REPL.
Normally, you enter the Julia REPL by simply typing `julia` in the terminal.
But there are some standard bells and whistles that make this particular project easier to use,
    so I've written a shell script called `juliarepl` to use instead.

Normally, one calls shell scripts by typing something like `bash juliarepl`,
    where `bash` here is the shell language used to interpret the file.
However, I have this one indexed as an executable script,
    so `./juliarepl` is the preferred shorthand.
It'll use whatever shell language your terminal uses by default.
All the common shell languages use nearly identical syntax,
    so it *should* just work.

The main thing using the `juliarepl` script does for you
    is to import the package itself into the main namespace with the name `JOB`,
    and also to "use" the package so that all variables the package "exports"
    are also in the main namespace.
You can see the full list of exports at the top of `src/AdaptiveModes_041724.jl`.



### Job Management

#### Initialize a set of jobs

Using this repository, every run of ctrl-(ADAPT-)VQE gets its own "job directory".
Within this directory gets dumped all the variables used to define the run,
    a copy of the script used to actually *make* the run,
    a dangerously serialized states of the *current state* of the run
        which may or may not be completed,
    as well as any pdf plots generated during the run.

These directories are loosely organized within a `jobs` directory,
    which you may need to create if it doesn't already exist.
The organization is loosely enforced by always *creating* these directories
    with the `JOB.init_survey` function.

This function needs you to specify three things:
1. The script you want to run. This determines things like what pool you use, or whether or not to freeze parameters.
2. A "code" defining the molecular system you want to run.
3. A pulse duration.

Your choices of script are all the `<script>.jl` files in the top-level of the directory.
For example, `harmonics.jl` uses a trigonometric pool where real and imaginary components are considered separate modes,
    while `harmonics.complex.jl` combines real and imaginary components into a single mode.
You'll omit the `.jl` extension when passing the script to `JOB.init_survey`.
Note I've got some *other* scripts, which are not *job* scripts, that I always begin with a dunder (double underscore) `__`.
For example, `__study__frozen.jl` does some post-processing and generates new pdfs that go into a certain job directory.
Don't use it with `JOB.init_survey`.

Your choices of code are visible in `CtrlJobs/matrix/model`.
Every one of those files is a matrix, saved with numpy's serialization, with the name `<code>.npy`.
Like with the script, omit the `.npy` extension when passing the code to `JOB.init_survey`.
Adding in your own matrices to `CtrlJobs/matrix/model` is encouraged but non-trivial,
    and far beyond the scope of this tutorial, so talk to me before doing so.
Go-to proof of concept examples include "H215" for H2 at 1.5Å, or "lih30" for LiH at 3.0Å.

Finally, the pulse duration is just a float number, in nanoseconds.
I do not know what will happen if you use an integer for the pulse duration.
I do not recommend trying.
You can actually specify more than one at once, hence the name "survey",
    but you should use the `JOB.init_survey` function even when initializing just a single pulse duration,
    to ensure the job directory gets created where it belongs.

Here's an example, using the script `harmonics.jl` with the `H215.npy` matrix,
    for pulse durations of 24 and 30 nanoseconds:
```
> JOB.init_survey("harmonics", "H215", 24.0, 30.0)
```

This will initialize job directories in `jobs/H215_harmonics/T24.0/` and `jobs/H215_harmonics/T30.0/`.


#### Check on all the jobs for a given script/system

Job directories can be in several states:
1. Converged: This means the whole ctrl-ADAPT-VQE trajectory saved with this directory has attained some physically meaningful stopping condition.
2. Terminated: This means the trajectory has attained some computationally meaningful stopping condition.
3. Running: This means
4. None of the above

You can see the state of all jobs within a given directory via the `JOB.status` function.
Here's an example to get the status of the newly created jobs from the previous step,
    plus any others that happen to already exist with the same script and system:
```
> JOB.status("jobs/H215_harmonics/")
```

If a job is "converged", it means the whole ctrl-ADAPT-VQE trajectory saved within the job directory
    has attained some physically meaningful stopping condition.
If a job is "terminated", it means the trajectory has attained some, uh, *computationally* meaningful stopping condition.
Within the code, I consider any "converged" job to also be "terminated".
In this package, "convergence" is attained once either the pool is exhausted or all pool candidate scores fall below some threshold,
    while mere "termination" is additionally attained once some ridiculously large number of candidates have been added to the ansatz.

If a job is "running", it literally means that the job directory contains a *file* called `running`.
This has a *large correlation* with the job, um, running, in the background.
But the correlation factor is, uh, not necessarily 1.0 - more on that later...

The final "state" covers newly-created jobs,
    *as well as* any jobs that were interrupted before they could terminate.

Of course, if a job is in state (4), you probably want to *run* the job.
You *might* want to run the job if it is in state (2), but first you'd have to tweak its termination condition, which is a pain.
This really shouldn't be an issue with this particular package, so I won't bother to explain it here; ask if it becomes an issue.

Finally, you *might* want to run the job if it is in state (3) but you are absolutely 100% certain that it is not in fact already running.
In that case, you can simply delete the file called `running` and it will most certainly revert to state (4).
But I *do not know what happens,* if you do this to a job which *is* running.
It won't be good. It *might* crash your computer. It will certainly turn that job directory into junk. So I don't recommend testing it.

#### Start all the incomplete jobs for a given script/system

Since I do this on ARC all the time, I made starting jobs as painless as possible.
To start the newly created jobs from the previous step,
    plus any others that happen to already exist with the same script and system:
```
> JOB.start("jobs/H215_harmonics/")
```

It won't do a thing to converged and terminated jobs, and it displays a warning message for each "running" job.

But for all *other* jobs, it starts the job running *in the background*.

There's kinda a lot that goes behind the scenes here;
    I'll just let you trace the code and ask me questions about it at your leisure.
The main thing to know is that closing the terminal window where you ran `JOB.start` *might* interrupt the jobs (I'm not sure),
    closing VS Code (if you ran `JOB.start` from an integrated terminal) *should* interrupt the jobs (it had *better*),
    and there is some fancy-schmancy way of using the shell command `ps` to check up on it,
    but for unfortunate reasons partly my fault, its output is not very readable,
    and I've never really really gotten the hang of it.
(In my defence, the primary function of `JOB.start` is running things on ARC, where it works very smoothly,
    and the fact that it can also run things locally is merely a convenience.)

I guess the other thing you should know is that you absolutely can start multiple jobs running simultaneously,
    and I recommend doing so,
    but starting too many will slow your computer down or crash it.
Maybe keep it to within, I dunno, ten at a time?
Maybe less, depending on how many browser tabs you have open...



### Data Management

Now we come to some of the craziest, bad-practice things I've done with this `CtrlJobs` repo. ^_^

#### Load a specific job into the active memory:

Say you've got a job directory, it's run through ctrl-ADAPT-VQE, and it's got all these files in it.
Some of the pictures it generated are nice.
But all these other files, VS Code has no idea what to do with them.
You'd really like to be able to play around with the actual numbers and such.

Let's say you're interested in the 24ns run we created earlier.
"Load" it like this:
```
> load!("jobs/H215_harmonics/T24.0")
```
If you get a serialization error...yeah. That's because of my bad practices.
It means the job was generated with a different, presumably later, version of Julia.
If you start the REPL with the correct version, it should work.

*Assuming* it works, this dumps everything there is to see in a global variable called `_!`.

There are all sorts of things that `_!` gives you access to. Let's just look at a couple specific examples:

#### Check out the parameters currently optimized for this job

Let's say you want the actual numerical parameters that this job optimized to.
Here they are:
```
> _!.state.x
```

Well, this is just a big array. What mode does each element correspond to?
The attributes `_!.state.Ω` and `_!.state.n` help a little with that.
Try to find the inline comments where they're defined for details.
But to completely answer the question,
    you'll need more intimate knowledge of the job script and how signals are defined.
Trying to match `_!.state.x` against, say, the `_amplitudes.pdf` file may also help.
Just make good note of units - there's a factor of 2π in there somewhere.
Good luck!

One other note here: the file `state` in the job directory contains the *current* state.
But files `a0001`, `a0002`, etc. contain the optimized states at each adaptation.
You can switch between them using the `unarchive!` function.
See comments in code for details, if it becomes relevant.
(These files were honestly originally meant to be merely backup files
    - in case something ever corrupted `state`, I could manually restart from the last adaptation.
But I wound up using them as an integral part of the scripts which freeze parameters.
Very hacky!)

#### Initialize work variables without actually running the job

Simply calling `load!` won't *completely* fill out the `_!` variable.
All the heavy-weight objects, such as the device or the cost-function or so on,
    are what I call "work" variables.
I can't serialize them directly, but they *can* be constructed robustly.
Indeed, they *are* constructed when you run the script!

But we can't run the script in the background, like we did before.
We'll have to `include` the script itself.
But wait! We don't want the whole bloody job to run in the foreground,
    taking up time and space we could be using for other tasks.

(Actually, if the job is already converged, you *can* let the whole bloody job run -
    it will immediately see that it has already converged and just stop.
But maybe you are trying to debug something in the middle of a trajectory,
    and you really *don't* want it to run.)

One of the many things in the `_!` variable is a boolean flag, `_!.run`.
It defaults to true.
But every job script has a line of code right before the optimization begins,
    which says to throw an error if `_!.run` is false.

So, first thing's first: do `_!.run = false`.
(Alternatively, you could have included a keyword argument `run=false` when you called `load!`.)

Next, we're going to `include` the script.

Now, in our examples above, we used the script `harmonics.jl`.
So we *could* include *that* script.
But I don't recommend it, because *perhaps* that script has changed in the interim since the job was created.
This is why one of the tasks of `JOB.init_survey` was to *copy* the script, to `script.jl`.
This copy becomes the archive-worthy version of `harmonics.jl`.
Any time the job is run, it actually uses `script.jl`, not `harmonics.jl`.
Once the job is initialized, `harmonics` is just a sub-string in the survey directory.

One other thing to make this part more convenient:
    `_!.outdir` contains the path of the job directory itself.
So, here's how we include the script:
```
> include("$(_!.outdir)/script.jl");
```

This is going to generate an error message, saying "Setup Finished".
That error means it worked.
If you don't want to see the error message,
    you can wrap it in an inert try-catch block like this:
```
> try include("$(_!.outdir)/script.jl"); catch end
```
This is actually necessary if you ever write your own post-processing script
    that needs to load work variables,
    since otherwise the "Setup Finished" error will interrupt the post-processing script.

Before this point, `_!.work` was just `nothing`.
Now, you wll find it filled with a variety of very useful objects.
The details, I leave for you to discover.