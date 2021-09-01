# parallel plate testcase
- 5000 e- are accelerated against the right electrode where secondary electron emission (SEE-E) is performed
- SEE model from A.I. Morozov, "Structure of Steady-State Debye Layers in a Low-Density Plasma near a Dielectric Surface", 2004
  that is used for dielectric walls 
  - SEE-E model yields three possible outcomes upon electron impact (hence XXX electrons should be created in this setup)
  - (i) the incident electron disappears and the wall acquires a negative charge
  - (ii) one secondary electron with energy eps is knocked out of the wall, and 
  - (iii) two secondary electrons with energies eps1 and eps2 are knocked out of the wall
- Coupled Power
  - output is compared with a reference value that has been tuned to give ~ X J
  - parameter.ini for activating coupled Power output

        CalcCoupledPower = T
