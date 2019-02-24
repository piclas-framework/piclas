# Chemistry
* Test the theoretical reaction rates for recombination to CH4
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* Reference databases shows acceptable agreement with Arrhenius + partition function results (error below 6%)
  * T=500K:  5.48E-036
  * T=1000K: 6.96E-039
  * T=3000K: 5.40E-043
  * T=5000K: 1.05E-044
* Initial composition: n=2e22 at equal parts CH3 and H
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false
