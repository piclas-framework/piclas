# Chemistry, 2D Axisymmetric with weighted particles
* Test case analogous to WEK_Reservoir/CHEM_EQUI_diss_CH4 but in a cylinder (L=4.64E-6, R=2.32E-6) with particle weights
* 
* Initial composition at 1.5e22 1/m3 with CH4 at T=7000K
* Settings
  * Symmetry
    * Particles-Symmetry2DAxisymmetric=T
  * Radial Weighting
    * Particles-DSMC-RadialWeighting=T
    * Particles-DSMC-RadialWeighting-PartScaleFactor=1
      * Weighting factor of 5E2 is increased to 1E2 towards the maximal y-dimension
    * Particles-DSMC-RadialWeighting-CloneMode=2
    * Particles-DSMC-RadialWeighting-CloneDelay=2
  * Chemistry
    * BackwardReacRate       = T
    * DSMCReservoirSim       = T
    * DSMCReservoirSimRate   = F
    * DSMCReservoirStatistic = F
    * Particles-DSMC-SelectionProcedure=1
    * Particles-DSMC-PolyRelaxSingleMode=F
    * Particles-DSMC-RotRelaxProb=0.2
    * Particles-DSMC-VibRelaxProb=0.02