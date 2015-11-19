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

! periodic case
INTEGER, ALLOCATABLE                     :: casematrix(:,:)                  ! matrix to compute periodic cases
INTEGER                                  :: NbrOfCases                       ! Number of periodic cases


! general: periodic sides have to be Cartesian
INTEGER,ALLOCATABLE :: SidePeriodicType(:)                                                ! 1:nTotalSides, periodic type of side
                                                                                          ! 0 - normal or BC side
                                                                                          ! >0 type of periodic displacement
REAL,ALLOCATABLE    :: SidePeriodicDisplacement(:,:)                                      ! displacement vector
                                                                                          
INTEGER,ALLOCATABLE :: PartElemToSide(:,:,:)                                              ! extended list: 1:2,1:6,1:nTotalElems
                                                                                          ! ElemToSide: my geometry + halo
                                                                                          ! geometry + halo information
                                                                                          
INTEGER,ALLOCATABLE :: PartSideToElem(:,:)                                                ! extended list: 1:5,1:6,1:nTotalSides
                                                                                          ! SideToElem: my geometry + halo
                                                                                          ! geometry + halo information

INTEGER,ALLOCATABLE :: PartElemToElem(:,:,:)                                              ! Mapping from ElemToElem
                                                                                          ! 1:2,1:6,1:nTotalElems
                                                                                          ! 1 - E2E_NB_ELEM_ID 
                                                                                          ! 2 - E2E_NB_LOC_SIDE_ID 
INTEGER             :: nTotalSides                                                        ! total nb. of sides (my+halo)
INTEGER             :: nTotalElems                                                        ! total nb. of elems (my+halo)

LOGICAL,ALLOCATABLE :: IsBCElem(:)                                                        ! is a BC elem 
                                                                                          ! or BC in halo-eps distance to BC
INTEGER             :: nTotalBCSides                                                      ! total number of BC sides (my+halo)
INTEGER             :: nTotalBCElems                                                      ! total number of bc elems (my+halo)
INTEGER,ALLOCATABLE :: PartBCSideList(:)                                                  ! mapping from SideID to BCSideID

REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: XiEtaZetaBasis                                 ! element local basis vector (linear elem)
REAL,ALLOCATABLE,DIMENSION(:,:)         :: slenXiEtaZetaBasis                             ! inverse of length of basis vector
REAL,ALLOCATABLE,DIMENSION(:,:)         :: ElemBaryNGeo                                   ! element local basis: origin
REAL,ALLOCATABLE,DIMENSION(:)           :: ElemRadiusNGeo                                 ! radius of element 
REAL,ALLOCATABLE,DIMENSION(:)           :: ElemRadius2NGeo                                ! radius of element + 2% tolerance
INTEGER                                 :: MappingGuess                                   ! select guess for mapping into reference
                                                                                          ! element
REAL                                    :: epsMapping                                     ! tolerance for Netwton to get xi from X
REAL                                    :: epsInCell                                      ! tolerance for eps for particle 
                                                                                          ! inside of ref element
REAL                                    :: epsOneCell                                     ! tolerance for particle in 
                                                                                          ! inside ref element 1+epsinCell

!LOGICAL                                 :: DoRefMapping                  ! tracking by mapping particle into reference element
!-----------------------------------------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tFastInitBGM
  INTEGER                                :: nElem                             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)                        ! List of elements/physical cells in BGM cell
#ifdef MPI     
  INTEGER, ALLOCATABLE                   :: ShapeProcs(:)                     ! first Entry: Number of Shapeprocs, 
                                                                              ! following: ShapeProcs
  INTEGER, ALLOCATABLE                   :: PaddingProcs(:)                   ! first Entry: Number of Paddingprocs, 
                                                                              ! following: PaddingProcs
  INTEGER, ALLOCATABLE                   :: SharedProcs(:)                    ! first Entry: Number of Sharedprocs, 
                                                                              ! following: SharedProcs
  !INTEGER                                :: nBCSides                          ! number BC sides in BGM cell
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
  ! required for cartesian BGM for desposition
  INTEGER, ALLOCATABLE                   :: PeriodicBGMVectors(:,:)           ! = periodic vectors in backgroundmesh coords
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
  INTEGER,ALLOCATABLE                    :: ElemToFIBGM(:,:)                  ! range of FIGMB cells per element
                                                                              ! 1:6,1:nTotalElems, xmin,max,yminmax,...
  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
  REAL, ALLOCATABLE                      :: DeltaEvMPF(:)                     ! Energy difference due to particle merge
  INTEGER, ALLOCATABLE                   :: ElemToRegion(:)                   ! ElemToRegion(1:nElems)

  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
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
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)          ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)          ! Link part 2 for mapping Boltzplatz BCs to Particle BC
!  INTEGER              , ALLOCATABLE     :: Map(:)                      ! Map from Boltzplatz BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)              ! Map from Boltzplatz BCindex to Particle BC (NOT TO TYPE!)
  !!INTEGER              , ALLOCATABLE     :: SideBCType(:)            ! list with boundary condition for each side
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

INTEGER                                  :: NbrOfRegions      ! Nbr of regions to be mapped to Elems
REAL, ALLOCATABLE                        :: RegionBounds(:,:) ! RegionBounds ((xmin,xmax,ymin,...)|1:NbrOfRegions)


!===================================================================================================================================


END MODULE MOD_Particle_Mesh_Vars
