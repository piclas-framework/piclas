# Flow through channel: Testing of surface output
* Constant inflow with a fixed emission current (2E-4A)
* Output in ParticleAnalyze.csv of emission current through surface flux should correspond to input value
* Output in SurfaceAnalyze.csv of current exiting the domain should be equal to the inflow
* Output in DSMCSurfState (CalcSurfaceImpact = T): absolute comparison of ImpactFlux