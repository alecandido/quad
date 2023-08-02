# Quad

This library provides quadrature integration methods for vector-valued real function.

The adaptive quadrature are based on Gauss-Kronrod rules.
"A Gauss-Kronrod rule begins with a classical Gaussian quadrature rule of order m. This is extended with additional
points between each of the abscissae to give a higher order Kronrod rule of order 2m + 1.
The Kronrod rule is efficient because it reuses existing function evaluations from the Gaussian rule."

At every step of the algorithm up to 128 sub-interval could be bisected, thus allowing parallelization.
