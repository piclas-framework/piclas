# PIC (Poisson) Discharge
- Test the following functionalities in combination with p-adaption
- RF discharge between two electrodes with Helium
    - Ionization of background gas
- Secondary electron emission at x-direction boundaries
    - only happens rarely for this setup (but code coverage of SEE modules is still performed)
    - Emission model: SEE-I (bombarding electrons are removed, Argon ions on different materials is considered for secondary e- emission with 0.13 probability) by Depla2009
- The results (total number of particles) are shown in *n_parts_total.jpg* for *MPI=1,...,10*
- Count the number of electrons that impact on the left boundary and output to SurfaceAnalyze.csv (compare the last line with a reference)

        CalcBoundaryParticleOutput = T 
- using rounded values for all particle masses
    - otherwise fails when using CORE_SPLIT and MPI=5 due to
    "Message: CODE_ANALYZE: DSMC_Chemistry is not energy conserving for chemical reaction:       1"
    - not that rel. Energy difference : -1.18806002360921E-012
- Output of temporal derivative of the electric field

        CalcElectricTimeDerivative =T
- Output of min/amax/average electron energy (eV) in each cell

        CalcElectronEnergy = T

- Test both load balances schemes. Load balance via hdf5 I/O and via direct MPI communication without I/O

        UseH5IOLoadBalance = T,F

- Test record point analysis over line with defined in recordpoints.ini

        Line_nRP              = 50                  ! Number of RPs on line
        Line_xstart           = (/0.0    , 1.71E-05 , 1.71E-05/) ! Coordinates of start of line
        Line_xend             = (/6.7e-2 , 1.71E-05 , 1.71E-05/) ! Coordinates of end of line

  and parameter.ini

        RP_SamplingOffset        = 1 ! only each x time step
        RP_DefFile               = Line_RPSet.h5
        RP_inUse                 = T
