# Chemistry with Background Gas
* Test the theoretical reaction rates for dissociation of CO2 using a background gas
* Verifying the correct TCE implementation by determining the reaction rate through reaction probabilities
* CO2: 2E23 (background gas), CO and O: 1E21 (particles)
* Reference databases show good agreement with Arrhenius results without electronic excitation
  * M=CO: Deviation below 1%
  * M=O: Deviation below 5% for T=7500K-15000K and 12% at 5000K
* With electronic excitation deviation from reference is within 10%, worse agreement is probably due to the small sample size at these relatively low temperatures for electronic excitation

| T [K]  | Arrhenius: CO2 + CO [m3/s] | Arrhenius: CO2 + O [m3/s] |
| :----: | :-----------------------: | :----------------------: |
| 5000K  |         1.037E-19         |        2.065E-19         |
| 7500K  |         3.835E-18         |        7.638E-18         |
| 10000K |         2.053E-17         |        4.089E-17         |
| 12500K |         5.209E-17         |        1.037E-16         |
| 15000K |         9.213E-17         |        1.835E-16         |