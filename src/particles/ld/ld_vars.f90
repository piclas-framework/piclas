MODULE MOD_LD_Vars
!===================================================================================================================================
! Contains the LD variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
  LOGICAL                           :: UseLD
  REAL  , ALLOCATABLE               :: LD_RHS(:,:)  ! RHS of the LD Method/ deltaV (npartmax, direction)
  REAL                              :: LD_RepositionFak
  REAL                              :: LD_RelaxationFak
  LOGICAL  , ALLOCATABLE           :: IsDoneLagVelo(:)  ! (nSides) 
  REAL  , ALLOCATABLE              :: TempDens(:)
!!!  REAL  , ALLOCATABLE               :: NewNodePosIndx(:,:) ! (1:nDim,1:nNodes)  !!! nur für "Tetra-Methode"
  LOGICAL                           :: LD_CalcDelta_t
  LOGICAL                           :: LD_CalcResidual
  REAL, ALLOCATABLE                :: LD_Residual(:,:) ! Def. for LD Residual (number of Elements, 2nd index: 
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Velocity |u|
                                                                                                  ! 5.Mass Density
                                                                                                  ! 6.Temperature
TYPE tLD_SecantMethod
  REAL                              :: Guess        ! 2nd guess, plus user defined value, [m/s], (default 10 m/s)
  REAL                              :: MaxIter      ! Max. number of iterations for LAGRANGIAN vell calculation
  REAL                              :: Accuracy     ! accuracy for LAGRANGIAN vell calculation
END TYPE
TYPE (tLD_SecantMethod)             :: LD_SecantMeth

TYPE tBulkValues                                                  ! LD bulk values
  REAL                              :: CellV(3)                   ! Velocity (vx,vy,vz)
  REAL                              :: DegreeOfFreedom            ! Average number of internal degree of freedom
  REAL                              :: Beta                       ! Thermal speed scale
  REAL                              :: MassDens                   ! Mass Density
  REAL                              :: BulkTemperature            ! BulkTemperature
!  INTEGER                           :: CellType                   ! CellType (1:DSMC, 2: Bufferzone_A, 3: Bufferzone_B, 4: LD
END TYPE
TYPE(tBulkValues), ALLOCATABLE      :: BulkValues(:)              ! LD bulk values array (number of Elements)

TYPE tBulkValuesOpenBC                                                  ! LD bulk values
  REAL                              :: CellV(3)                   ! Velocity (vx,vy,vz)
  REAL                              :: DegreeOfFreedom            ! Average number of internal degree of freedom
  REAL                              :: Beta                       ! Thermal speed scale
  REAL                              :: MassDens                   ! Mass Density
END TYPE
TYPE(tBulkValuesOpenBC), ALLOCATABLE :: BulkValuesOpenBC(:)       ! LD bulk values array for open BCs (number of Elements)

TYPE tSurfLagValues                                               ! LD Lagrangian cellsurface values
  REAL                              :: LagVelo                    ! Lagrangian Velocity
  REAL                              :: Area                       ! area of cellsurface
  REAL                              :: DeltaM(3)                  ! momentumflux (Mx,My,Mz)
  REAL                              :: DeltaE                     ! energyflux
  REAL                              :: LagNormVec(3)              ! corresponding face normal vector, [-]
END TYPE
TYPE(tSurfLagValues), ALLOCATABLE  :: SurfLagValues(:,:,:)        ! LD Lagrangian cellsurface array (iLocSide, iElem, TriNum)

TYPE tMeanSurfValues                                              ! LD Lagrangian cellsurface values
  REAL                              :: MeanBaseD                  ! Mean Base (Base Vector in coordinat form of Mean Surface)
  REAL                              :: MeanBaseD2                 ! Mean Base for periodoc walls
  REAL                              :: MeanNormVec(3)             ! corresponding face normal vector, [-]
  REAL                              :: MeanLagVelo                ! Lagrangian Velocity for MeanSurf
END TYPE
TYPE(tMeanSurfValues), ALLOCATABLE  :: MeanSurfValues(:,:)          ! Mean Surface for LD-Particle push (iLocSide, iElem)
 
REAL    , ALLOCATABLE               :: PartStateBulkValues(:,:)   ! LD particle values (npartmax, with 2nd index: 
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Temperature
                                                                                                  ! 5.Degree of freedom
#ifdef MPI
  REAL  , ALLOCATABLE               :: MPINeighborBulkVal(:,:) ! LD values for cells on other procs (SideID, with 2nd index: 
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Temperature
                                                                                                  ! 5.Density
#endif

!===================================================================================================================================

END MODULE MOD_LD_Vars
