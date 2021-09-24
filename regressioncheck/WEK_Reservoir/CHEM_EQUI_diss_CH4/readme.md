# Chemistry
* Test dissociation and recombination (CH4 + M <-> CH3 + H + M) from a thermal and chemical non-equilibrium towards equilibrium
* Reference database compared with analytical solution (see PhD thesis, Nizenkov)
* Result with 5 x particles saved in PartAnalyze_T7000K_w100.csv
* Initial composition at 1.5e22 1/m3 with CH4 at T=7000K
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = false
  * DSMCReservoirStatistic = false
  * Particles-DSMC-SelectionProcedure=2
  * Particles-DSMC-PolyRelaxSingleMode=true
  * Particles-DSMC-RotRelaxProb=0.2
  * Particles-DSMC-VibRelaxProb=0.006