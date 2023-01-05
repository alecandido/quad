# Criteria for simultaneous integration

Essentially, the core constituent of the `quad_vec` implementation.
Simultaneous criteria are defined on top of individual stopping criteria.

Possibly, define a general interface for criteria, in such a way that each
criterion can be applied to each integration algorithm, without the explicit
`M x N` support (i.e. keep it `M + N`).

## Simple

Basic simple criteria are:

- **and**: this criterion consists in keep integrating until all the individual
  components satisfy their individual criteria simultaneously
  - this is useful because all the points computed are always used to increase
    the accuracy of the integrals (assuming the expensive operation is to
    evaluate the array of functions, and not to transform)
  - it might be suboptimal, since adding points might deteriorate already
    converged integrals, and thus makes the stopping condition more difficult
    to satisfy; in the most extreme case, integrals that are separately
    convergent might not converge any longer, also for an arbitrary number of
    evaluations (oscillatory behavior might probe float precision at some point)
- **early stop**: this is the most conservative criterion, consisting in
  stopping the individual integration as soon as it converged with its
  individual criterion; new points are only used for dimensions not yet
  converged
  - though some evaluations are discarded, this algorithm is still a net
    improvement over the separate integration, since in the worst case one
    dimension is by far more expensive than any other, so you save a negligible
    amount of evaluations (corresponding to the very early-stopped dimensions),
    but it is only asymptotically equally expensive, and strictly better for any
    finite number of evaluations; in the common case in which an arbitrary
    rotation/transformation is applied to the initially evaluated vector,
    complications are spread over all the dimensions, and the evaluations
    required are comparable, in this case the gain factor would be equal to the
    number of dimensions of the final vector (the iterated would require as many
    integrations, and the simultaneous a single one, equivalent to all the
    others); this is also the maximally convenient case

These criteria are sufficient for _non-adaptive_ integrations, since only the
stopping has to be decided.
In the case of _adaptive_ algorithms also the positions of new evaluation points
have to be decided with a simultaneous criterion, so one has to be provided also
for that.

## Adaptive sampling criteria

## Advanced stopping criteria
