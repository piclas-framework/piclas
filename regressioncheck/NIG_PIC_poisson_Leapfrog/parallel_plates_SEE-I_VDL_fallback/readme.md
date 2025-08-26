# parallel plate testcase

- 100000 Ar+ ions (a lot are required to statistically meet the 1% SEE yield of one of the models) are accelerated against the right electrode where secondary electron emission (SEE-I) is performed
  - the 2nd case considers only 2500 Ar+ ions with vMPF=T, they are instantly split with Part-Species2-vMPFSplitThreshold = 100 into multiple particles with a smaller MPF
  - when the electrons reach the left electrode, they accommodate with the wall temperature of 300 K, which effectively stops their movement

- Note that due to the applied electric field, the emitted electrons are accelerated after emission and therefore their velocity and
  energy distribution changes as compared with the distribution given by the emission model

- Two SEE-I models are tested
  - SEE model 9
      - the SEE-I model yields an emission ratio of 1 % (hence 1000 electrons should be created in this setup)

- Coupled Power
    - IMPORTANT: deactivated output of coupled Power and comparison in analysis because this value oscillates to much and is
      skipped during comparison therefore anyway
    - parameter.ini for activating coupled Power output:

              CalcCoupledPower = T

- Count the number of emitted electrons from SEE surfaces
    - parameter.ini for activating this output:

              CalcCurrentSEE = T

- Due to the VDL model and the slope of the right BC, some particles are shifted outside of the domain and are stored in
  parallel_plates_PartStateLost_000.00000050000000000.h5
