# BGK Reservoir - Relaxation of CH4
* Test the relaxation from thermal non-equilibrium with ESBGK (CollModel=1) and SBGK (CollModel=2)
* Continuous and quantized treatment of vibrational energy
* Restart from a given state file (additional test if VibQuantData is given)
* T_t = 10000K, T_r=7500K, T_v = 5000K
* Reference file was generated with a weighting factor of 100
* Relevant parameter
  * Particles-BGK-CollModel = 1,2
  * Particles-SBGK-EnergyConsMethod=1 (only relevant for CollModel=2)
  * Particles-BGK-DoVibRelaxation=T
  * Particles-BGK-UseQuantVibEn=T,F
  * Particles-BGK-DoCellAdaptation=F
  * Particles-BGK-MinPartsPerCell=20
  * PartDensity           = 1E23
  * MacroParticleFactor   = 500
  * RotRelaxProb          = 0.2
  * VibRelaxProb          = 0.05
* Theoretical equilibrium temperature is 5960K (see [Pfeiffer, M., Nizenkov, P., Mirza, A., and Fasoulas, S. (2016). Direct simulation Monte Carlo modeling of relaxation processes in polyatomic gases. Physics of Fluids 28(2), 027103.](http://dx.doi.org/10.1063/1.4940989))