# 3D HE-thruster Liu2010
- Boris-leapfrog time integration with cell volweight mean deposition on a 3D domain (90 deg segment of a cylinder)
- external magnetic field in radial direction (algebraic expression)

    PIC-AlgebraicExternalField      = 3 ! 3: 3D Liu 2010 magnetic + electric field
    PIC-AlgebraicExternalFieldDelta = 2 ! Integer

- MCC with variable background gas distribution from .h5 file

    Particles-BGGas-UseDistribution              = T
    Particles-MacroscopicRestart-Filename        = neutral-DSMC/3D_HET_Liu2010_DSMCState_000.00100000000000000.h5

- SEE model with variable electron bulk temperature (which changes the behaviour of the SEE yield)
  - Morozov (2004) with variable SEE yield depending on electron impact properties

      Part-Boundary6-SurfaceModel         = 8 ! SEE-E (bombarding ions are reflected, e- on dielectric materials is considered for secondary e- emission with different probabilities for different outcomes) by Morozov2004
      Part-Species2-PartBound6-ResultSpec = 2 ! impacting e- (Part-Spec is 2) results in emission of e- (ResultSpec is 2)
      Part-SurfaceModel-SEE-Te            = 5.80226250308285e5 ! = 50 eV / Electron temperature in K: 5.80226250308285e5 K corresponds to 50 eV, 1.16045250061657e4 K corresponds to 1 eV

  - with variable electron bulk temperature (determined globally)

      Part-SurfaceModel-SEE-Te-automatic  = T ! Instead of using a fixed bulk electron temperature, determine the global temperature of the defined species (default is False). Note that Part-SurfaceModel-SEE-Te is used as initial value.
      Part-SurfaceModel-SEE-Te-Spec       = 2 ! For automatic bulk Te determination, state the species ID of the electrons

  - surface charging of dielectric walls (inner BCs)

      Part-Boundary6-Dielectric           = T

- Dielectric: Outer and inner walls are dielectric via hollow circle type region

    DielectricTestCase     = HollowCircle ! Dielectric region is outside of outer radius and inside of inner radius

- Neutralization emission BC via keeping the exiting charge at the right BC zero over time (averaged)

    Part-Species2-Init2-SpaceIC = 3D_Liu2010_neutralization

  or by enforcing a neutral boundary layer at the right exit by emitting electrons if there is an ion surplus in the first row of
  elements

    Part-Species2-Init2-SpaceIC = 3D_Liu2010_neutralization_Szabo

- For 10 processors an asymmetric CVWM deposition occurs when some processors deposit charge on the dielectric interface, the nodes
  of which are not connected directly to their own mesh cells

- post-piclas-vMPF-restart
  - runs piclas and restarts from 3D_HET_Liu2010_State_000.00000000100000000.h5 (simulation without vMPF) with vMPF=T to test the restart functionality
  - uses Part-Species-vMPFMergeThreshold = 200 to merge electrons and ions to approx. 50 percent of their initial number (700k to  350k particles)
