# Chemistry
* Reservoir simulation of high-temperature air dissociating, testing the treatment of multiple reaction paths (e.g. N2 + O2 can result in two different reactions)
* Reference database compared with analytical solution (Haas1993_Temp.csv) of Haas, B. L., & McDonald, J. D. (1993). Validation of chemistry models employed in a particle simulation method. Journal of Thermophysics and Heat Transfer, 7(1), 42–48. https://doi.org/10.2514/3.11567
* Deviation towards t=1e-7s is due to the calculation of backward rates through equilibrium constant as opposed to the usage of Arrhenius rates (as was done by Haas1993), results in a different equilibrium states
* Initial composition at n=2.45E24 (79% N2, 21% O2) at T=30000K
* Relaxation probabilities set to 1 to emulate analytical solution
* Settings
  * BackwardReacRate       = true
  * Particles-DSMC-RotRelaxProb=1
  * Particles-DSMC-VibRelaxProb=1
