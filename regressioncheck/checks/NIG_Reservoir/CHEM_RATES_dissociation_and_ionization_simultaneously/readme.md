# Chemistry
* Test the theoretical reaction rates for dissociation and recombination as well as ionization and
    recombination of N2
* Verifying the correct Q-K implementation by determining the reaction rate through reaction probabilities
* Reference databases show acceptable agreement with Arrhenius results
  * Forward (Deviation to reference database below X%)
    * T=10000K: 
    * T=15000K: 
    * T=20000K: 
    * T=25000K: 
    * T=30000K: 
  * Backward (Deviation to reference database below X%, largest deviation at T=0K)
    * T=10000K: 
    * T=15000K: 
    * T=20000K: 
    * T=25000K: 
    * T=30000K: 
* Initial composition: n=4E23 at equal parts N2, e, N2+ and N
* Settings
  * BackwardReacRate       = T
  * DSMCReservoirSim       = T
  * DSMCReservoirSimRate   = T
  * DSMCReservoirStatistic = T
* Reactions

| Reaction   | Type     | Number    |
| :--------: | :------: | :-------: |
| N2+N2      | Diss.    | 1         |
| N2+e       | Ioni.    | 2         |
| N2+e       | Diss.    | 3         |
| N2+N2+     | Diss.    | 4         |
| N2+N       | Diss.    | 5         |
| N2         | R(D)     | 6         |
| e          | R(I)     | 7         |
| e          | R(D)     | 8         |
| N2+        | R(D)     | 9         |
| N          | R(D)     | 10        |
