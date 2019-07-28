# FP-Flow Reservoir - Relaxation of CH4
* Test the relaxation from thermal non-equilibrium with ESFP (CollModel=2)
* Continuous and quantized treatment of vibrational energy
* Restart from a given state file (additional test if VibQuantData is given)
* T_t = 10000K, T_r=7500K, T_v = 5000K
* Reference file was generated with a weighting factor of 100 (copied from BGK case)
* Relevant parameter
  * Particles-FP-CollModel = 2
  * Particles-ESFP-Model=1
  * Particles-FP-DoVibRelaxation=T
  * Particles-FP-UseQuantVibEn=T,F
  * Particles-FP-DoCellAdaptation=F
  * Particles-FP-MinPartsPerCell=20
  * PartDensity           = 1E23
  * MacroParticleFactor   = 500
  * RotRelaxProb          = 0.2
  * VibRelaxProb          = 0.05
* Theoretical equilibrium temperature is 5960K (see [Pfeiffer, M., Nizenkov, P., Mirza, A., and Fasoulas, S. (2016). Direct simulation Monte Carlo modeling of relaxation processes in polyatomic gases. Physics of Fluids 28(2), 027103.](http://dx.doi.org/10.1063/1.4940989))