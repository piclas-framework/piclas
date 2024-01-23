# DSMC - Hypersonic flow around a 70° Cone (Axisymmetric) - Radiation pipeline
* Simulation of a hypersonic N2-N flow around a 70° blunted cone
* Test case based on 70degCone-Reggie, however, since the radiation tool chain shall be tested, this one has the FIRE II inflow conditions at 76 km altitude. In addition, N was added as an inflow species
* Comparison of the heat flux with a reference surface state file
* The wake of the heat shield was not included due to strong fluctuations and long sampling duration prohibiting the use of the test case as a regression test
* Additionally, the "Particles-RadialWeighting-PartScaleFactor" is reduced to safe computational time
