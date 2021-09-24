# Chemistry
* Test the theoretical reaction rates for dissociation and electron-impact ionization of N2 + e
* Verifying the correct Q-K implementation for the possibility that both reactions can happen within a single collision
* Reference databases for dissociation show acceptable agreement with Arrhenius results (no Arrhenius rates for ionization)
  * Forward (Deviation to reference database below 33%)
    * T=10000K: 2.40523109598669E-17
    * T=15000K: 5.47203215244321E-16
    * T=20000K: 2.27832452157028E-15
    * T=25000K: 4.94516648457136E-15
    * T=30000K: 7.85666964653137E-15
* Initial composition: n=2E23 at equal parts N2 and e
* Settings
  * BackwardReacRate       = F
  * DSMCReservoirSim       = T
  * DSMCReservoirSimRate   = T
  * DSMCReservoirStatistic = T
* Reactions

| Reaction | Type  | Number |
| :------: | :---: | :----: |
|   N2+e   | Diss. |   1    |
|   N2+e   | Ioni. |   2    |