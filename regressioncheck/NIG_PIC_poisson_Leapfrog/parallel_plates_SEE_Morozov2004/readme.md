# parallel plate testcase
- 10000 e- are accelerated against the right electrode where secondary electron emission (SEE-E) is performed
  - the 2nd case considers only 2500 e- with vMPF=T, they are instantly split with Part-Species1-vMPFSplitThreshold = 1000 into multiple particles with a smaller MPF
- SEE model from A.I. Morozov, "Structure of Steady-State Debye Layers in a Low-Density Plasma near a Dielectric Surface", 2004
  that is used for dielectric walls 
  - SEE-E model yields three possible outcomes upon electron impact (hence XXX electrons should be created in this setup)
  - (i) the incident electron disappears and the wall acquires a negative charge
  - (ii) one secondary electron with energy eps is knocked out of the wall, and 
  - (iii) two secondary electrons with energies eps1 and eps2 are knocked out of the wall
- Compare number of particles: PartAnalyze.csv
  - output is compared with a reference value for const. collision energy of primary electrons of 50, 100 and 200 eV
