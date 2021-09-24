# Chemistry
* Test the theoretical reaction rates for dissociation and recombination of N2
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * Forward (Deviation to reference database below 2%)
    * T=10000K: 8.9321949508785E-20
    * T=15000K: 2.03212315210182E-18
    * T=20000K: 8.46090790277447E-18
    * T=25000K: 1.83646349735145E-17
    * T=30000K: 2.91769489695039E-17
  * Backward (Deviation to reference database below 50%, largest deviation at T=30000K)
    * T=10000K: 5.9710E-46
    * T=15000K: 2.4908E-46
    * T=20000K: 1.2991E-46
    * T=25000K: 7.9302E-47
    * T=30000K: 5.3854E-47
* Initial composition: n=2E22 at equal parts N2 and N
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false
