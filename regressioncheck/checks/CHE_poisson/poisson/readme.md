# PIC (Poisson) Discharge
- RF discharge between two electrodes with Helium
    - Ionization of background gas
- Secondary electron emission at x-direction boundaries
    - only happens rarely for this setup (but code coverage of SEE modules is still performed)
    - Emission model: SEE-I (bombarding electrons are removed, Argon ions on different materials is considered for secondary e- emission with 0.13 probability) by Depla2009
- The results (total number of particles) are shown in *n_parts_total.jpg* for *MPI=1,...,10*
- CalcBoundaryParticleOutput = T
    - count the number of electrons that impact on the left boundary and output to SurfaceAnalyze.csv (compare the last line with a
    reference)
- using rounded values for all particle masses
    - otherwise fails when using CORE_SPLIT and MPI=5 due to
    "Message: CODE_ANALYZE: DSMC_Chemistry is not energy conserving for chemical reaction:       1"
    - not that rel. Energy difference : -1.18806002360921E-012
