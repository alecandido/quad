# References

Main theory:

- https://en.wikipedia.org/wiki/Numerical_integration
- https://en.wikipedia.org/wiki/Gaussian_quadrature

Main implementations:

- QUADPACK: https://nines.cs.kuleuven.be/software/QUADPACK/
  - structure and conventions: https://en.wikipedia.org/wiki/QUADPACK
- GSL: https://www.gnu.org/software/gsl/doc/html/integration.html
- scipy: https://scipy.github.io/devdocs/reference/integrate.html

## Papers & articles

- A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric integration
  over an N-dimensional rectangular region," J. Comput. Appl. Math., vol. 6 (no.
  4), 295-302 (1980)
- Gladwell, I. (1987). Vectorisation of One Dimensional Quadrature Codes. In:
  Keast, P., Fairweather, G. (eds) Numerical Integration. NATO ASI Series, vol 203.
  Springer, Dordrecht. https://doi.org/10.1007/978-94-009-3889-2_24
- J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm for the
  approximate calculation of multiple integrals," ACM Trans. Math. Soft., vol.
  17 (no. 4), 437-451 (1991)
- Jarle Berntsen, Terje O. Espelid, and Alan Genz. 1991. Algorithm 698: DCUHRE:
  an adaptive multidemensional integration routine for a vector of integrals.
  ACM Trans. Math. Softw. 17, 4 (Dec. 1991), 452–456.
  https://doi.org/10.1145/210232.210234
  - this is the reference implementation of Genz 1980 and Genz 1991
- J.M. Bull, T.L. Freeman, Parallel globally adaptive algorithms for
  multi-dimensional integration, Applied Numerical Mathematics, Volume 19,
  Issues 1–2, 1995, Pages 3-16, ISSN 0168-9274,
  https://doi.org/10.1016/0168-9274(95)00076-7

Cuba:

- Hahn, T., “CUBA—a library for multidimensional numerical integration”,
  Computer Physics Communications, vol. 168, no. 2, pp. 78–95, 2005.
  doi:10.1016/j.cpc.2005.01.010.
- Hahn, T., “Concurrent Cuba”, in Journal of Physics Conference Series, 2015,
  vol. 608, no. 1. doi:10.1088/1742-6596/608/1/012066.

Other Wikipedia relevant articles:

- https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas
- https://en.wikipedia.org/wiki/Category:Numerical_integration_(quadrature)
- https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula
- https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature
- https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula
- https://en.wikipedia.org/wiki/Truncation_error_(numerical_integration)
- https://en.wikipedia.org/wiki/Adaptive_quadrature
- https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
- https://en.wikipedia.org/wiki/Bayesian_quadrature

## More implementations

- QUADPACK2: https://github.com/jacobwilliams/quadpack
  - same QUADPACK (in theory), but refactored
  - modern Fortran (possibly more readable)
  - available on GitHub (instead of on Netlib, with a bit weird interface...)
  - inlined dependencies (QUADPACK depended on other libraries)
  - in practice: a single 9000 lines Fortran90 file
- Cubature: https://github.com/stevengj/cubature
  - implements vectorization following Gladwell
    https://github.com/stevengj/cubature#vectorized-interface
- Cuba library https://feynarts.de/cuba/

### Julia

There are many packages available in the Julia ecosystem, even though there is
no one as authoritative as SciPy or GSL:

One-dimensional:

- initially part of Julia Base: https://github.com/JuliaMath/QuadGK.jl
- https://github.com/JuliaApproximation/FastGaussQuadrature.jl
  - just provides weights
  - see [note in `QuadGK`](https://github.com/JuliaMath/QuadGK.jl#similar-packages)
  - inspired by [Chebfun](http://www.chebfun.org/)

Multi-dimensional:

- https://github.com/JuliaMath/HCubature.jl
  - implements Genz 1980
- https://github.com/pabloferz/NIntegration.jl
  - implements Genz 1991
- https://github.com/JuliaMath/Cubature.jl
  - this is also remarkable for the notes in the README.md, explaining something
    on the quadrature ecosystem
  - wraps the C library Cubature implementing Genz's papers and Gladwell for
    vectorization
- https://github.com/giordano/Cuba.jl
  - implements Cuba suite in Julia

## Relevant Rust crates

- Rust bindings to Cuba: https://github.com/benruijl/cuba
- Rust integrators traits: https://github.com/jhod0/integrators

## Further techniques

- Quasi-Monte Carlo in Python https://qmcpy.org/
