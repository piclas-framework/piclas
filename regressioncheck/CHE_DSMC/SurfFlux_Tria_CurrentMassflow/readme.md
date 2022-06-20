# DSMC - Emission Current/Mass flow through Surface Flux
* Particle emission using surface flux by defining an emission current of 2A or a mass flow of 1.1372E-11 kg/s
* Comparing the current and mass flow output (CalcSurfFluxInfo = T), which is based on the particle balance (in - out)
* Comparing the number of inserted particles per time step, testing with
  * Velocity: 0 and 1E6 m/s
  * Temperature: 5K and 500 000K
* Expected value is 6241/6242 simulation particles (corresponds to 1.25E+19 real particles per second)
  * w = 2e5
  * dt = 1e-10