# Chemistry
* Test the theoretical reaction rates for the exchange rate: CH4 + H <-> CH3 + H2
  * Forward reaction rate from: Baulch, D. L. (2005). Evaluated Kinetic Data for Combustion Modeling: Supplement II. Journal of Physical and Chemical Reference Data, 34(3), 757. https://doi.org/10.1063/1.1748524, p. 966 (p. 210 in PDF), "Preferred Values"
  * Backward reaction rate: Calculation through equilibrium constant and partition sums
* Verifying the correct TCE implementation and automatic backward reaction rates by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * Forward (Deviation to reference database below 7%)
    * T= 2500K: 4.62659882666551E-17
    * T= 5000K: 6.86959300337277E-16
    * T= 7500K: 2.61130267450647E-15
    * T=10000K: 6.29583132722603E-15
  * Backward (Deviation to reference database below 25%)
    * T= 2500K: 2.58049442614499E-18
    * T= 5000K: 5.11465974525123E-17
    * T= 7500K: 2.56848298297464E-16
    * T=10000K: 7.75498009039655E-16
* Initial composition: n=4E+22 at equal parts CH4, CH3, H2 and H
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false