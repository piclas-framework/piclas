# Chemistry
* Test the theoretical reaction rates for dissociation of CH4 with different collision partner (CH4, H2, H)
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * Forward (Deviation of reference database between 5% and 70% (greater error at lower temperatures))
    * T=2000K: 9.30667349721537E-23
    * T=3000K: 6.45483346357596E-20
    * T=4000K: 8.47045706992641E-19
    * T=5000K: 2.6227470394402E-18
    * T=6000K: 4.23134040912595E-18
  * Backward (Deviation of reference database below 6%)
    * T=2000K: 1.49066395794696E-41
    * T=3000K: 5.39890380048881E-43
    * T=4000K: 5.69483168886693E-44
    * T=5000K: 1.04555682001156E-44
    * T=6000K: 2.68952129104968E-45
* Initial composition: n=4e22 at equal parts CH4, CH3, H2 and H
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false
