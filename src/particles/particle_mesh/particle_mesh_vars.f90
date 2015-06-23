#include "boltzplatz.h"

MODULE MOD_Particle_Mesh_Vars
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES


LOGICAL             :: ParticleMeshInitIsDone
!-----------------------------------------------------------------------------------------------------------------------------------
! Mesh info
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: SidePeriodicType(:)                                                ! periodic type of side
                                                                                          ! 0 - normal or BC side
                                                                                          ! >0 type of periodic displacement
REAL,ALLOCATABLE    :: SidePeriodicDisplacement(:,:)                                      ! dispacement vector
                                                                                          
INTEGER,ALLOCATABLE :: PartElemToSide(:,:,:)                                              ! containing the ElemToSide of my
                                                                                          ! geometry + halo information
                                                                                          
INTEGER,ALLOCATABLE :: PartSideToElem(:,:)
INTEGER,ALLOCATABLE :: PartNeighborElemID(:,:)
INTEGER,ALLOCATABLE :: PartNeighborlocSideID(:,:)
INTEGER             :: nTotalSides
INTEGER             :: nTotalElems

LOGICAL,ALLOCATABLE :: IsBCElem(:)
INTEGER             :: nTotalBCSides
INTEGER             :: nTotalBCElems
INTEGER,ALLOCATABLE :: PartBCSideList(:)
!-----------------------------------------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tFastInitBGM
  INTEGER                                :: nElem                             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)                        ! List of elements in BGM cell
#ifdef MPI     
  INTEGER, ALLOCATABLE                   :: ShapeProcs(:)                     ! first Entry: Number of Shapeprocs, 
                                                                              ! following: ShapeProcs
  INTEGER, ALLOCATABLE                   :: PaddingProcs(:)                   ! first Entry: Number of Paddingprocs, 
                                                                              ! following: PaddingProcs
  INTEGER, ALLOCATABLE                   :: SharedProcs(:)                    ! first Entry: Number of Sharedprocs, 
                                                                              ! following: SharedProcs
  INTEGER                                :: nBCSides                          ! number BC sides in BGM cell
#endif                     
END TYPE

INTEGER                                  :: FIBGMCellPadding(1:3)

TYPE tGeometry
  REAL                                   :: xminglob                          ! global minimum x coord of all nodes
  REAL                                   :: yminglob                          ! global minimum y coord of all nodes
  REAL                                   :: zminglob                          ! global minimum z coord of all nodes
  REAL                                   :: xmaxglob                          ! global max x coord of all nodes
  REAL                                   :: ymaxglob                          ! global max y coord of all nodes
  REAL                                   :: zmaxglob                          ! global max z coord of all nodes
  REAL                                   :: xmin                              ! minimum x coord of all nodes
  REAL                                   :: xmax                              ! maximum x coord of all nodes
  REAL                                   :: ymin                              ! minimum y coord of all nodes
  REAL                                   :: ymax                              ! maximum y coord of all nodes
  REAL                                   :: zmin                              ! minimum z coord of all nodes
  REAL                                   :: zmax                              ! maximum z coord of all nodes
  ! periodic
  INTEGER                                :: nPeriodicVectors                  ! Number of periodic Vectors
  REAL, ALLOCATABLE                      :: PeriodicVectors(:,:)              ! PeriodicVectors(1:3,1:nPeriodicVectors), 1:3=x,y,z
  INTEGER,ALLOCATABLE                    :: DirPeriodicVectors(:)             ! direction of periodic vectors
  LOGICAL                                :: directions(3)                     ! flag for direction
  ! FIBGM
  REAL                                   :: FIBGMdeltas(3)                    ! size of background mesh cell for particle init
  REAL                                   :: FactorFIBGM(3)                    ! scaling factor for FIBGM

  ! caution, possible pointer
  TYPE (tFastInitBGM),ALLOCATABLE        :: FIBGM(:,:,:)  !        =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: FIBGMimin                         ! smallest index of FastInitBGM (x)
  INTEGER                                :: FIBGMimax                         ! biggest index of FastInitBGM (x)
  INTEGER                                :: FIBGMjmin                         ! smallest index of FastInitBGM (y)
  INTEGER                                :: FIBGMjmax                         ! biggest index of FastInitBGM (y)
  INTEGER                                :: FIBGMkmin                         ! smallest index of FastInitBGM (z)
  INTEGER                                :: FIBGMkmax                         ! biggest index of FastInitBGM (z)
  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
  REAL, ALLOCATABLE                      :: DeltaEvMPF(:)                     ! Energy difference due to particle merge

!  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
END TYPE

TYPE (tGeometry)                         :: GEO


TYPE tBCElem
  INTEGER                                :: nInnerSides                       ! Number of BC-Sides of Element
  INTEGER                                :: lastSide                          ! total number of BC-Sides in eps-vicinity of element
  INTEGER, ALLOCATABLE                   :: BCSideID(:)                       ! List of elements in BGM cell
END TYPE

TYPE (tBCElem),ALLOCATABLE               :: BCElem(:)

TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleAnodeBC           = 4      ! = 4 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SimpleCathodeBC         = 5      ! = 5 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)!=>NULL() ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)!=>NULL() ! Link part 2 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: Map(:)            !=>NULL() ! Map from Boltzplatz BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)    !=>NULL() ! Map from Boltzplatz BCindex to Particle BC (NOT TO TYPE!)
  !INTEGER              , ALLOCATABLE     :: SideBCType(:)            ! list with boundary condition for each side
  REAL    , ALLOCATABLE                  :: MomentumACC(:)      
  REAL    , ALLOCATABLE                  :: WallTemp(:)     
  REAL    , ALLOCATABLE                  :: TransACC(:)     
  REAL    , ALLOCATABLE                  :: VibACC(:) 
  REAL    , ALLOCATABLE                  :: RotACC(:) 
  REAL    , ALLOCATABLE                  :: WallVelo(:,:) 
  REAL    , ALLOCATABLE                  :: Voltage(:)
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          !Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         !Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           !Species to be changed at wall (in, out), out=0: delete
  LOGICAL , ALLOCATABLE                  :: AmbientCondition(:)
  REAL    , ALLOCATABLE                  :: AmbientTemp(:)
  REAL    , ALLOCATABLE                  :: AmbientMeanPartMass(:)
  REAL    , ALLOCATABLE                  :: AmbientBeta(:)
  REAL    , ALLOCATABLE                  :: AmbientVelo(:,:)
  REAL    , ALLOCATABLE                  :: AmbientDens(:)
  REAL    , ALLOCATABLE                  :: AmbientDynamicVisc(:)               ! dynamic viscousity
  REAL    , ALLOCATABLE                  :: AmbientThermalCond(:)               ! thermal conuctivity
END TYPE

INTEGER                                  :: nPartBound                       ! number of particle boundaries
TYPE(tPartBoundary)                      :: PartBound                         ! Boundary Data for Particles




!===================================================================================================================================


END MODULE MOD_Particle_Mesh_Vars
