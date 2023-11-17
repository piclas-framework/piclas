# PIC - Thermionic Emission
* Particle emission using surface flux and thermionic emission
* Schottky effect reduces the work function (and thus the required wall temperature)
* Comparing the calculated current to the expected value of the Richardson Dushman equation for Tungsten
  * Input: W = 4.54 eV, A = 60 A/(cm^2 K^2), T_w = 2700 K (without Schottky), T_w = 2635.24 K
  * Output: j = 1.47 A / cm^2 -> I = 1.4675 A (A = 1 cm^2)
  * In the case with the Schottky effect, the current reduces slightly over time as the electrons in the domain reduce the potential difference
