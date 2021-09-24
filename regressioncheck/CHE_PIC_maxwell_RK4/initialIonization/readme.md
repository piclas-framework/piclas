# Initial Ionization
* Testing the initial ionization routine:
  * Read-in of neutral particles from restart file (1000 particles)
  * Ionization according to the user-defined charge average
* Checking bounds of particle positions to avoid particles outside the domain
* Compare the number of particles after a single time step with reference file
  * Charge average of 0.5 should result in particle numbers of 500 for each species