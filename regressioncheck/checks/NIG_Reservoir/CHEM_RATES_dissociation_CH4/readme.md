# Chemistry
* Test the theoretical reaction rates for dissociation of CH4 with different collision partner (CH4, H2, H)
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * T=2000K: 9.30667349721537E-23
  * T=3000K: 6.45483346357596E-20
  * T=4000K: 8.47045706992641E-19
  * T=5000K: 2.6227470394402E-18
  * T=6000K: 4.23134040912595E-18
* Initial composition: n=3e22 at equal parts CH4, H2 and H
* Settings
  * BackwardReacRate       = false
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false
