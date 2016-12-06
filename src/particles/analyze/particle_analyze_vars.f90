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
LOGICAL                       :: CalcNumSpec                           ! Calculate the number of simulated particles per species
LOGICAL                       :: CalcCollRates                         ! Calculate the collision rates per collision pair
LOGICAL                       :: CalcReacRates                         ! Calculate the reaction rate per reaction
LOGICAL                       :: CalcEpot                              ! Computation of the energy stored in the electric and
                                                                       ! magnetic field
LOGICAL                       :: CalcEkin                              ! Compute the kinetic energy of each species
LOGICAL                       :: CalcEint                              ! Compute the internal energy of each species
LOGICAL                       :: CalcTemp                              ! Computation of the temperature (trans, rot, vib, total)
LOGICAL                       :: CalcPartBalance                       ! Particle Power Balance - input and outflow energy of all
                                                                       ! particles
LOGICAL                       :: CalcVelos                             ! Computes the drift and thermal velocity of each species
LOGICAL                       :: VeloDirs(4)                           ! select the direction for velo computation
LOGICAL                       :: TrackParticlePosition                 ! track the particle movement
                                                                       ! stored in .csv format, debug only, no MPI 
INTEGER                       :: nEkin                                 ! number of kinetic energies 
LOGICAL                       :: IsRestart                             ! check if restart, add data to Database
LOGICAL                       :: ChargeCalcDone                        ! check flag
LOGICAL                       :: CalcShapeEfficiency                   ! efficiency of shape function
CHARACTER(LEN=256)            :: CalcShapeEfficiencyMethod             ! Explanations in particle_analyze.f90
INTEGER                       :: ShapeEfficiencyNumber                 ! Explanations in particle_analyze.f90
INTEGER                       :: PartAnalyzeStep                       ! Analyze is performed each Nth time step
INTEGER,ALLOCATABLE           :: nPartIn(:)                            ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartOut(:)                           ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartInTmp(:)                         ! Number of entry and leaving particles
REAL,ALLOCATABLE              :: PartEkinIn(:)                         ! energy and temperatur of input particle
REAL,ALLOCATABLE              :: PartEkinOut(:)                        ! energy and temperatur of input particle
REAL,ALLOCATABLE              :: PartEKinInTmp(:)                      ! energy and temperatur of input particle
LOGICAL                       :: CalcCharge                            ! Compute the whole deposited charge and abs and relative
                                                                       ! charge error
LOGICAL                       :: DoVerifyCharge                        ! validate the charge after each deposition and produces
                                                                       ! an output in std.out
REAL                          :: PartCharge(3)                         ! contains the whole deposited charge and its absolute
                                                                       ! and relative error
LOGICAL                       :: printDiff
REAL                          :: printDiffTime
REAL                          :: printDiffVec(6)
!===================================================================================================================================
END MODULE MOD_Particle_Analyze_Vars
