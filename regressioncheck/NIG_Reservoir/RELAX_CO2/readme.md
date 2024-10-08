# Relaxation of CO2 (polyatomic)
* Test relaxation of CO2 from a thermal non-equilibrium to equilibrium (rotational, vibrational)
* Restart from a given state file (additional test if VibQuantData is given)
* T_t = 10000K, T_r=7500K, T_v = 5000K
* Settings
  * PartDensity           = 1E23
  * MacroParticleFactor   = 500
  * RotRelaxProb          = 0.2
  * VibRelaxProb          = 0.02
* Theoretical equilibrium temperature is 6549K (see [Pfeiffer, M., Nizenkov, P., Mirza, A., and Fasoulas, S. (2016). Direct simulation Monte Carlo modeling of relaxation processes in polyatomic gases. Physics of Fluids 28(2), 027103.](http://dx.doi.org/10.1063/1.4940989))
