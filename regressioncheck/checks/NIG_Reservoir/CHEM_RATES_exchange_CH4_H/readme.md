# Chemistry
* Test the theoretical reaction rates for the exchange rate: CH4 + H <-> CH3 + H2
  * Forward reaction rate from: Baulch, D. L. (2005). Evaluated Kinetic Data for Combustion Modeling: Supplement II. Journal of Physical and Chemical Reference Data, 34(3), 757. https://doi.org/10.1063/1.1748524, p. 966 (p. 210 in PDF), "Preferred Values"
  * Backward reaction rate: Calculation through equilibrium constant and partition sums
* Verifying the correct TCE implementation and automatic backward reaction rates by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * Forward: Arrhenius (Deviation to reference database)
    * T= 2500K: 4.62659882666551E-17 (4.39%)
    * T= 5000K: 6.86959300337277E-16 (9.49%)
    * T= 7500K: 2.61130267450647E-15 (41.61%)
    * T=10000K: 6.29583132722603E-15 (67.60%)
    * Large deviation at higher temperature due to reaction probabilities above 1 (which are limited by PICLas to 1)
  * Backward: Partition Function (Deviation to reference database)
    * T= 2500K: 2.58049442614499E-18 (26.58%)
    * T= 5000K: 5.11465974525123E-17 (17.19%)
    * T= 7500K: 2.56848298297464E-16 (8.69%)
    * T=10000K: 7.75498009039655E-16 (1.21%)
* Initial composition: n=4E+22 at equal parts CH4, CH3, H2 and H
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false