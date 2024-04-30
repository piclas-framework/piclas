# parallel plate testcase

- 100000 Ar+ ions (a lot are required to statistically meet the 1% SEE yield of one of the models) are accelerated against the right electrode where secondary electron emission (SEE-I) is performed
  - the 2nd case considers only 2500 Ar+ ions with vMPF=T, they are instantly split with Part-Species2-vMPFSplitThreshold = 100 into multiple particles with a smaller MPF
  - when the electrons reach the left electrode, they accommodate with the wall temperature of 300 K, which effectively stops their movement

- The following parameters are varied

        TrackingMethod                           = refmapping        , tracing        , triatracking
        Part-Boundary2-SurfaceModel              = 7                 , 7              , 9
        Part-Boundary2-SurfModEnergyDistribution = deltadistribution , uniform-energy
        Part-Boundary2-SurfModEmissionEnergy     = -1.0              , 7.0
        Part-Boundary2-SurfModEmissionYield      = 0.13              , 1.9
        Part-vMPF                                = F                 , T
        Part-SpeciesX-MacroParticleFactor        = 32e6              , 1.28e9
        Part-Species2-Init1-ParticleNumber       = 100000            , 2500

- Part-Boundary2-SurfModEnergyDistribution, Part-Boundary2-SurfModEmissionEnergy and Part-Boundary2-SurfModEmissionYield are only
  relevant for Part-Boundary2-SurfaceModel=7 (not for 9)

- Note that due to the applied electric field, the emitted electrons are accelerated after emission and therefore their velocity and
  energy distribution changes as compared with the distribution given by the emission model

- Two SEE-I models are tested
  - SEE model 7 from D. Depla "Magnetron sputter deposition: Linking discharge voltage with target properties", 2009
      - various metal surfaces with different yields
      - the SEE-I model yields an emission ratio of 13 % (hence 13000 electrons should be created in this setup)
  - SEE model 9
      - the SEE-I model yields an emission ratio of 1 % (hence 1000 electrons should be created in this setup)

- Coupled Power
    - IMPORTANT: deactivated output of coupled Power and comparison in analysis because this value oscillates to much and is
      skipped during comparison therefore anyway
    - parameter.ini for activating coupled Power output:

              CalcCoupledPower = T

- Count the number of emitted electrons from SEE surfaces
    - parameter.ini for activating this output:

              CalcElectronSEE = T
