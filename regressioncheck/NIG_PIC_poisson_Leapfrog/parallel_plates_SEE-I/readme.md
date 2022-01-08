# parallel plate testcase
- 100000 Ar+ ions (a lot are required to statistically meet the 1% SEE yield of one of the models) are accelerated against the right electrode where secondary electron emission (SEE-I) is performed
  - the 2nd case considers only 2500 Ar+ ions with vMPF=T, they are instantly split with Part-Species2-vMPFSplitThreshold = 10000 into multiple particles with a smaller MPF
- Two SEE-I models are tested
  - SEE model 7 from D. Depla "Magnetron sputter deposition: Linking discharge voltage with target properties", 2009
      - various metal surfaces with different yields
      - the SEE-I model yields an emission ratio of 13 % (hence 13000 electrons should be created in this setup)
  - SEE model 9
      - the SEE-I model yields an emission ratio of 1 % (hence 1000 electrons should be created in this setup)
- Coupled Power
    - parameter.ini for activating coupled Power output:

              CalcCoupledPower = T
- Count the number of emitted electrons from SEE surfaces
    - parameter.ini for activating this output:

              CalcElectronSEE = T
