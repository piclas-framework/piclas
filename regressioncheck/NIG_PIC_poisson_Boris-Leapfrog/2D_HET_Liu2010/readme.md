# 2D HE-thruster Liu2010
- Boris-leapfrog time integration with cell volweight mean deposition on a 2D domain
- external magnetic field in radial direction (algebraic expression)

    PIC-AlgebraicExternalField      = 2 ! 2: 2D Liu 2010 magnetic + electric field
    PIC-AlgebraicExternalFieldDelta = 2 ! Integer

- MCC with variable background gas distribution from .h5 file

    Particles-BGGas-UseDistribution              = T
    Particles-MacroscopicRestart-Filename        = pre-BGGas/2D_HET_Liu2010_DSMCState_000.0007000000.h5

- SEE model with variable electron bulk temperature (which changes the behaviour of the SEE yield)
  - Morozov (2004) with variable SEE yield depending on electron impact properties

      Part-Boundary3-SurfaceModel         = 8 ! SEE-E (bombarding ions are reflected, e- on dielectric materials is considered for secondary e- emission with different probabilities for different outcomes) by Morozov2004
      Part-Species2-PartBound3-ResultSpec = 2 ! impacting e- (Part-Spec is 2) results in emission of e- (ResultSpec is 2)
      Part-SurfaceModel-SEE-Te            = 5.80226250308285e5 ! = 50 eV / Electron temperature in K: 5.80226250308285e5 K corresponds to 50 eV, 1.16045250061657e4 K corresponds to 1 eV

  - with variable electron bulk temperature (determined globally)

      Part-SurfaceModel-SEE-Te-automatic  = T ! Instead of using a fixed bulk electron temperature, determine the global temperature of the defined species (default is False). Note that Part-SurfaceModel-SEE-Te is used as initial value.
      Part-SurfaceModel-SEE-Te-Spec       = 2 ! For automatic bulk Te determination, state the species ID of the electrons

- Neutralization emission BC via keeping the exiting charge at the right BC zero over time (averaged)

    Part-Species2-Init2-SpaceIC = 2D_Liu2010_neutralization

  or by enforcing a neutral boundary layer at the right exit by emitting electrons if there is an ion surplus in the first row of
  elements

    Part-Species2-Init2-SpaceIC = 2D_Liu2010_neutralization_Szabo

- Output of particle flux, total electric current and emitted SEE over time into SurfaceAnalyze.csv for particle boundaries 1,2 and 3via

    CalcBoundaryParticleOutput = T
    BPO-NPartBoundaries        = 3         ! Nbr of boundaries
    BPO-PartBoundaries         = (/1,2,3/) ! Part-Boundary1 to Part-Boundary3
    BPO-NSpecies               = 2         ! Nbr of species
    BPO-Species                = (/2,3/)   ! electrons, Xe+
