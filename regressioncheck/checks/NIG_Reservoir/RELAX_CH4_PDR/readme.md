# Relaxation of CH4 - Prohibiting Double Relaxation
* Test relaxation of CH4 from a thermal non-equilibrium to equilibrium (rotational, vibrational) using an alternative relaxation procedure
* Restart from a given state file (additional test if VibQuantData is given)
* T_t = 10000K, T_r=7500K, T_v = 5000K
* Settings
  * PartDensity           = 1E23
  * MacroParticleFactor   = 500
  * RotRelaxProb          = 0.2
  * VibRelaxProb          = 0.006 (probability is not coupled to the rotational probability, value chosen to achieve equilibrium is time similar to the RELAX_CH4 case)
* Theoretical equilibrium temperature is 5960K (see [Pfeiffer, M., Nizenkov, P., Mirza, A., and Fasoulas, S. (2016). Direct simulation Monte Carlo modeling of relaxation processes in polyatomic gases. Physics of Fluids 28(2), 027103.](http://dx.doi.org/10.1063/1.4940989))