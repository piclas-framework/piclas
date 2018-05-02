MODULE MOD_LoadBalance_Vars
!===================================================================================================================================
! Variables needed for the evaluation of the record points
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: DoLoadBalance                              ! DoLoadBalance
INTEGER                             :: LoadBalanceSample                          ! Number of samples for loadbalance
LOGICAL                             :: PerformLBSample                            ! Flag for enabling time measurement in current
                                                                                  ! timestep (automatically set depending on LB
                                                                                  ! sampling method)
LOGICAL                             :: InitLoadBalanceIsDone                      ! switch for checking

! time measurement
REAL,ALLOCATABLE                    :: tTotal(:)                                  ! time measurement over whole dt_analyze 
                                                                                  ! measured elem-independent and later weighted
REAL,ALLOCATABLE                    :: tCurrent(:)                                ! time measurement over one step
                                                                                  ! measured elem-independent and later weighted
                                                                                  ! for indeces look into boltzplatz.h
! counter
INTEGER                             :: nLoadBalance                               ! number of load balances
INTEGER                             :: nLoadBalanceSteps                          ! number of performed  load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)                              ! Weighted load distribution of all procs
INTEGER,ALLOCATABLE                 :: PartDistri(:)                              ! Part distribution of all procs
REAL                                :: MaxWeight                                  ! Maximum Weight of proc on domain
REAL                                :: MinWeight                                  ! Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                               ! Imbalance after rebalance step

TYPE tData
  INTEGER, ALLOCATABLE :: offsetElemMPI(:)
  INTEGER              :: numOfCalls
  TYPE(tData), POINTER :: nextData => null()
END TYPE tData
TYPE(tData), POINTER :: firstData => null() !linked-list of old offsetElemMPI for WeightDistributionMethod 5 and 6

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: nSkipAnalyze                               ! Skip Analyze-Dt
REAL                                :: ParticleMPIWeight
REAL                                :: DeviationThreshold                         ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
REAL                                :: WeightSum                                  ! global sum of all weights
REAL                                :: targetWeight                               ! optimal weight for each proc

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nDeposPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nTracksPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacefluxPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerBCElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacePartsPerElem(:)


END MODULE MOD_LoadBalance_Vars
