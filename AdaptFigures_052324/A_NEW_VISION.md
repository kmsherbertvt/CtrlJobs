Nick and I put our heads together and came up with a nice organizational strategy for ctrl-ADAPT-VQE.


ctrl-VQE gives us maximal freedom to identify pulses for state preparation.
How can we greedily add parameters in to our ansatz to get the best energy for each parameter?

We can consider two broad paradigms:
- Subdivision - a top-down approach:
    Each pulse is piece-wise constant. Optimize it as one window. Then break it up into two, and so on.
- Superposition - a bottom-up approach:
    Each pulse is a linear combination of basis functions. Add in one function at a time.

Within each paradigm, we can formulate three strategies:
- Naive - Add in parameters sequentially, based on a pre-defined order.
- Adaptive - Select candidates from a pool based on the largest gradient, ala ADAPT-VQE.
- Optimal - Use the gradient particularly cleverly to make a particularly intelligent choice.

So there are six different protocols we'll discuss:
- Subdivision/Naive: "Uniform" (in code: `uniform.each`) (no warm start)
- Subdivision/Adaptive: "Bisectal" (in code: `bisection`)
- Subdivision/Optimal: "~~Optimal~~ Nodal" (in code: `nodes.one`) (requires gradient signal)
- Superposition/Naive: "Harmonics (Iterative)" (in code: `harmonics.complex.iterative`)
- Superposition/Adaptive: "Harmonics" (in code: `harmonics`)
- Superposition/Optimal: "Modal Harmonics" (in code: `gradients.modalharmonics`)

  This last is one that Nick just intuited,
    reasoning that the (locally) _optimal_ basis would be the gradient signals themselves.
  So the best use of whatever basis you actually have is the _linear combination_
    which best approximates the gradient signals.
  But lo! this is _precisely_ what the pool gradient gives you!

  So, at each adaptation, one records the entire pool gradient.
  And in each optimization, one sticks a coefficient in front of each pool gradient vector,
    and optimizes those linear combinations.
  Ah, that is technically doable with what I have,
    but it will require an even crazier hack than what I did to get the "frozen" runs to work.
  I'll need to load each adaptation, _calculate scores_, and construct a...
    ...sigh, a _weighted signal_ of _composite signals_.
  Each component in the composite signal will be a _constrained harmonic_.
  In particular, each component is a deepcopy of the pool signal,
    bound with the corresponding element of the pool gradient,
    and then (hack) adapted into a constrained harmonic with _all_ parameters frozen.

  Almost complicated enough to justify a whole new project directory, but not quite.
  But it is well worth the effort because lo! it's a perfectly viable strategy
    for coupling pulses together with a single parameter;
    Nick thinks it should work very well, and if it does,
    it might just obviate the need for the rather more daunting natural/normal modes project.
  Or at the very least teach us a great deal about how to best go about implementing it.


