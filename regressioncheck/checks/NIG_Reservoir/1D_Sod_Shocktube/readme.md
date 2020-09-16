# DSMC - Sod Shocktube
* Simulation of a Sod shock tube test case by Gary A Sod (1978), A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws, Journal of Computational Physics, Volume 27, Issue 1, 1-31, https://doi.org/10.1016/0021-9991(78)90023-2.
* Basic test case for Riemann solvers
* Regression check for 1D flow problems as well as the flow pattern
* Executed until 2.25E-5 seconds and sampled for the last 3% of execution time
* High amount of allowed differences in analyse.ini needed due to mostly zero velocity in y and z direction as well of an increased particle number in the reference file
