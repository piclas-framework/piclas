# MCC - Species Time Step
* Particle emission using surface flux by defining an emission current of 2A
* Species specific time step for the electrons (factor 0.2)
  1. Part-VariableTimeStep-DisableForMCC = F: Collisions and chemistry are performed at the species-specific time step
  2. Part-VariableTimeStep-DisableForMCC = T: Collisions and chemistry are performed at the ManualTimeStep to accelerate convergence
* Comparing the number density with different reference files as
  1. Option is expected to ionize at a slower rate
  2. Reference file for the 2. option has been generated with a regular simulation at the lower electron time step (dt = 1e-11s)