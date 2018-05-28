MODULE MOD_Particle_Analyze_Vars
!===================================================================================================================================
! Contains global variables used by the Analyze modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: ParticleAnalyzeInitIsDone = .FALSE.
LOGICAL                       :: CalcNumSpec                         ! Calculate the number of simulated particles per species
LOGICAL                       :: CalcCollRates                       ! Calculate the collision rates per collision pair
LOGICAL                       :: CalcReacRates                       ! Calculate the reaction rate per reaction
LOGICAL                       :: CalcEkin                            ! Compute the kinetic energy of each species
LOGICAL                       :: CalcEint                            ! Compute the internal energy of each species
LOGICAL                       :: CalcTemp                            ! Computation of the temperature (trans, rot, vib, total)
LOGICAL                       :: CalcPartBalance                     ! Particle Power Balance - input and outflow energy of all
                                                                     ! particles
LOGICAL                       :: CalcSurfNumSpec                     ! Calculate the number of simulated particles per species 
                                                                     ! on surfaces
LOGICAL                       :: CalcEvaporation                     ! Calculate rate of evaporation [kg/s]
LOGICAL                       :: CalcSurfCoverage                    ! Calculate the surface coverages for each species
LOGICAL                       :: CalcAccomodation                    ! Calculate the surface accommodation coefficient
LOGICAL                       :: CalcAdsorbRates                     ! Calculate the adsorption probabilities of species
LOGICAL                       :: CalcSurfRates                       ! Calculate the surface reaction rate per reaction (k_r)
LOGICAL                       :: CalcVelos                           ! Computes the drift and thermal velocity of each species
LOGICAL                       :: VeloDirs(4)                         ! select the direction for velocity computation
LOGICAL                       :: TrackParticlePosition               ! track the particle movement
                                                                     ! stored in .csv format, debug only, no MPI 
INTEGER                       :: nSpecAnalyze                        ! number of analyzed species 1 or nSpecies+1
LOGICAL                       :: IsRestart                           ! check if restart, add data to Database
LOGICAL                       :: ChargeCalcDone                      ! check flag
LOGICAL                       :: CalcShapeEfficiency                 ! efficiency of shape function
CHARACTER(LEN=256)            :: CalcShapeEfficiencyMethod           ! Explanations in particle_analyze.f90
INTEGER                       :: ShapeEfficiencyNumber               ! Explanations in particle_analyze.f90
INTEGER                       :: PartAnalyzeStep                     ! Analyze is performed each Nth time step
INTEGER,ALLOCATABLE           :: nPartIn(:)                          ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartOut(:)                         ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartInTmp(:)                       ! Number of entry and leaving particles
REAL,ALLOCATABLE              :: PartEkinIn(:)                       ! energy and temperature of input particle
REAL,ALLOCATABLE              :: PartEkinOut(:)                      ! energy and temperature of input particle
REAL,ALLOCATABLE              :: PartEKinInTmp(:)                    ! energy and temperature of input particle

! get derived particle properties (for IMD/TTM initialization these values are calculated from the TTM grid values)
LOGICAL                       :: CalcDebyeLength                     ! Flag to compute the Debye length (min and max) in each cell
LOGICAL                       :: CalcPICTimeStep                     ! Flag to compute the PIC time step (min and max) in each cell
LOGICAL                       :: CalcElectronDensity                 ! Flag to compute the electron density in each cell
LOGICAL                       :: CalcElectronTemperature             ! Flag to compute the electron temperature in each cell
!LOGICAL                       :: ElectronTemperatureIsMaxwell        ! Assumption of Maxwell-Boltzmann or undistributed electrons 
LOGICAL                       :: CalcPlasmaFrequency                 ! Flag to compute the electron frequency in each cell
LOGICAL                       :: CalcPointsPerDebyeLength            ! Flag to compute the points per Debye length:
!                                                                    ! PPD=(p+1)lambda_D/L_cell
LOGICAL                       :: CalcIonizationDegree                ! Flag to compute the ionization degree per cell
REAL,ALLOCATABLE              :: IonizationCell(:)                   ! ionization degree cell value
REAL,ALLOCATABLE              :: PPDCell(:)                          ! Points per Debye length cell value
REAL,ALLOCATABLE              :: DebyeLengthCell(:)                  ! Debye length cell value
REAL,ALLOCATABLE              :: PICTimeStepCell(:)                  ! Approximated PIC Time Step per cell
REAL,ALLOCATABLE              :: ElectronDensityCell(:)              ! Electron density per cell
REAL,ALLOCATABLE              :: ElectronTemperatureCell(:)          ! Electron temperature per cell
REAL,ALLOCATABLE              :: PlasmaFrequencyCell(:)              ! plasma electron frequency per cell

LOGICAL                       :: CalcCharge                          ! Compute the whole deposited charge and abs and relative
                                                                     ! charge error
LOGICAL                       :: DoVerifyCharge                      ! validate the charge after each deposition and produces
                                                                     ! an output in std.out
REAL                          :: PartCharge(3)                       ! contains the whole deposited charge and its absolute
                                                                     ! and relative error
LOGICAL                       :: printDiff
REAL                          :: printDiffTime
REAL                          :: printDiffVec(6)
REAL                          :: ChemEnergySum
!===================================================================================================================================
END MODULE MOD_Particle_Analyze_Vars
