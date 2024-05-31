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
- Subdivision/Adaptive: "Bisectal" (in code: `bisection.one`)
- Subdivision/Optimal: "~~Optimal~~ Nodal" (in code: `nodes.one`) (requires gradient signal)
- Superposition/Naive: "Harmonics (Iterative)" (in code: `harmonics.complex.iterative.cold`) (*could* be warm started, but let's not for the sake of consistency)
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

## Appendices

I'd like to consider the following subvariants, discussion of which may be relegated to appendices.
All these can just use LiH.

- Infinite basis limit for `sup/opt`: `gradients.exact`

  If the basis is complete, `sup/opt` is equivalent to taking the gradient signal itself as the next "direction".
  Like `div/opt`, this requires characterizing the complete gradient signal as a function of time.
  And, like `div/opt`, it doesn't actually appear to confer any particular advantage,
    at least for our choice of basis.

  TODO: The appendix should study how compactness converges with the size of the basis.
  We can fix pulse duration to one value (say 36ns), and manipulate fMAX, and run `gradients.modalharmonics`

  Mathematically, it's interesting to note that this protocol without freezing is entirely equivalent to GRAPE
    (gradient descent with continuously variable pulses),
    and with freezing, it is extremely evocative of a Krylov subspace method.
  Which suggests, perhaps, another variant akin to Lanczos.

  TODO: At any rate I think that the difference between Krylov and Lanczos is basically
    that Lanczos does some sort of SVD to get the "best" orthonormal basis for any given Krylov subspace.
  Which is a thing we could certainly do, though I'd definitely want to flesh it out in its own project,
    i.e. not right now.

- Cold start: `harmonics.complex.iterative`

  You can imagine doing `sup/naive` either iteratively adding on one mode at a time,
    (letting parameters relax from their previously obtained values, i.e. "warm start"),
    or, like in `div/naive`, just starting fresh for each choice of parameters.
  Intuitively, warm-start feels like it is more sophisticated, so it should be better, right?
  But if you know ahead of time how many parameters you can handle,
    maybe it's best to just do the one optimization (i.e. not adaptively at all).
  We should understand the tradeoffs.

  TODO: The tradeoffs here are probably a function of BFGS convergence,
    so we should run multiple surveys each at a different
  We can fix pulse duration to one value (say 36ns), and manipulate gtol,
    and run both `harmonics.complex.iterative` and `harmonics.complex.iterative.cold`.
  Contrast `warmtrace.adaptations[end]` with `diff(coldtrace.adaptations)[end]`.


## Even more Appendices

- Freezing: `harmonics.complex.iterative.frozen`

  Normally we let previously added parameters relax, after adding a new one.
  What if we just optimize the new parameter?

  (It's bad.)

  This is well-defined for all the methods except `div/naive`, though I don't think it makes sense for any `div`.
  I'll note that doing this on the infinite basis limit for `sup/opt` is utterly equivalent to plain-old GRAPE
    (ie. gradient descent on a continuously variable pulse).

  To gain intuition for the effects of freezing without throwing too much at once,
    I think I like the idea of doing this and cold start particularly for `sup/naive`.
  Since this is trivially bad, we can maybe omit this.


- TETRIS: `nodes.oneeach`, `nodes.oneall`

  Our non-naive methods all add exactly one parameter at each iteration.
  But really, all these choices in the superposition and subdivision paradigms always commute.
  We can very easily add one parameter _per qubit_, similar to TETRIS-ADAPT.
  We can even add one parameter _per window_, in subdivision.
  Any particular advantages when we do?

  In truth, the results are not particularly interesting and this one could pry be omitted entirely,
    if we think there is too much in the appendix.