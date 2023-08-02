# Quadrature

Python bindings for [quad](../quad).

## Usage

The primary function is `qag`, which returns a `PyResult<QagsResult>`.
The interval extremes `a` and `b`, and the vector-valued function `f` must be provided to qag.
An example of use is the following:

```
a = 0.0
b = 1.0
f = lambda x: (math.sin(x), math.cos(x))
res = quad.qag(f, a, b)
```

Optional argument of `qag` are:

- epsabs: the absolute error,
- epsrel: the relative error,
- key: the Gauss-Kronrod rule that will be used,
- limit: the max number of subdivision allowed,
- points: a list of additional break-points,
- more_info: if set to `True` the `QagsResult` will contains also the number of function evaluation,
  the number of subdivision, and a list with all the sub-intervals and the respective error and result.
