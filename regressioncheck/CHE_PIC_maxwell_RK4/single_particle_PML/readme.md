# Single Particle PML
- Maxwell solver with RK4 time integration
- PML surrounding the inner physical domain
- one particle that gyrates in the domain (only deposition is active, no interpolation)
- Test initial restart for load balance in combination with dt=dt_Analyze (tests that the time steps and analyze times are correctly set)
- Test "numberedmulti" variable to globally set the MPF via Part-Species$-MacroParticleFactor in combination with automatic restart
