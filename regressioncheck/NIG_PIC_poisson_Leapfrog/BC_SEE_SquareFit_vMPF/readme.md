# Secondary Electron Emission - Square-fit model
* Particle emission using surface flux by defining an emission current of 2A
* Testing the calculation of the yield (through the comparison of the SEE current at the surface) for a square-fit expression a*E + b*E^2 + c
* Testing two different approaches when using Part-vMPF = T
  * Part-Boundary1-SurfMod-vMPF = F: Using Poisson distribution to determine the number of secondaries from the yield
  * Part-Boundary1-SurfMod-vMPF = T: Emit 1 secondary and scale its weight by the yield
* Yield should be the same, however, the error in the energy conservation is slightly different at this energy level