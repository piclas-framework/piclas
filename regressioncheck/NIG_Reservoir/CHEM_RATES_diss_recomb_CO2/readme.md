# Chemistry
* Test the theoretical reaction rates for dissociation/recombination of CO2
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* Reference databases show good agreement with Arrhenius results
  * Forward (Deviation to reference database below 2%)
    * T= 5000K: 1.03326632044077E-19
    * T= 7500K: 3.82136955663569E-18
    * T=10000K: 2.04589227106689E-17
    * T=12500K: 5.18992209119347E-17
    * T=15000K: 9.17947197464986E-17
  * Backward is compared against Arrhenius rate: 1.46E-46*EXP(2280.36/T)
    * NIST (1984WAR197C, Warnatz J. (1984) Rate Coefficients in the C/H/O System. In: Gardiner W.C. (eds) Combustion Chemistry. Springer, New York, NY)
    * Deviation is between 75% and 50%, acceptable qualitative agreement considering different approaches
    * T= 5000K: 2.30368195149908E-46
    * T= 7500K: 1.97878813874413E-46
    * T=10000K: 1.83395083063551E-46
    * T=12500K: 1.75218807938657E-46
    * T=15000K: 1.69971488272781E-46
* Initial composition: n=3e22 at equal parts CO2, CO and O
* Settings
  * BackwardReacRate       = true
  * DSMCReservoirSim       = true
  * DSMCReservoirSimRate   = true
  * DSMCReservoirStatistic = false