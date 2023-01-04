# References

Main theory:

- https://en.wikipedia.org/wiki/Gaussian_quadrature

Main implementations:

- QUADPACK: https://nines.cs.kuleuven.be/software/QUADPACK/
  - structure and conventions: https://en.wikipedia.org/wiki/QUADPACK
- GSL: https://www.gnu.org/software/gsl/doc/html/integration.html
- scipy: https://scipy.github.io/devdocs/reference/integrate.html

## More implementations

- QUADPACK2: https://github.com/jacobwilliams/quadpack
  - same QUADPACK (in theory), but refactored
  - modern Fortran (possibly more readable)
  - available on GitHub (instead of on Netlib, with a bit weird interface...)
  - inlined dependencies (QUADPACK depended on other libraries)
  - in practice: a single 9000 lines Fortran90 file
-

## More interesting software

- Rust bindings to Cuba: https://github.com/benruijl/cuba
- Rust integrators traits: https://github.com/jhod0/integrators
