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
  ! FIBGM
  REAL                                   :: FIBGMdeltas(3)                    ! size of background mesh cell for particle init
  
  ! caution, possible pointer
  TYPE (tFastInitBGM),ALLOCATABLE        :: FIBGM(:,:,:)  !        =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: FIBGMimin                         ! smallest index of FastInitBGM (x)
  INTEGER                                :: FIBGMimax                         ! biggest index of FastInitBGM (x)
  INTEGER                                :: FIBGMjmin                         ! smallest index of FastInitBGM (y)
  INTEGER                                :: FIBGMjmax                         ! biggest index of FastInitBGM (y)
  INTEGER                                :: FIBGMkmin                         ! smallest index of FastInitBGM (z)
  INTEGER                                :: FIBGMkmax                         ! biggest index of FastInitBGM (z)

!  INTEGER, ALLOCATABLE                   :: ElemToNodeID(:,:)                 ! ElemToNodeID(1:nElemNodes,1:nElems)
!  INTEGER, ALLOCATABLE                   :: ElemSideNodeID(:,:,:)             ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
!                                                                              ! From element sides to node IDs
!  INTEGER, ALLOCATABLE                   :: PeriodicElemSide(:,:)             ! 0=not periodic side, others=PeriodicVectorsNum
!  INTEGER, ALLOCATABLE                   :: PeriodicBGMVectors(:,:)           ! = periodic vectors in backgroundmesh coords
!  LOGICAL, ALLOCATABLE                   :: ConcaveElemSide(:,:)              ! Whether LocalSide of Element is concave side
!  REAL, ALLOCATABLE                      :: NodeCoords(:,:)                   ! Node Coordinates (1:nDim,1:nNodes)
!  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
!  REAL, ALLOCATABLE                      :: DeltaEvMPF(:)                     ! Energy difference due to particle merge
!  REAL                                   :: FactorFIBGM(3)                    ! scaling factor for FIBGM
!  REAL                                   :: xminglob                          ! global minimum x coord of all nodes
!  REAL                                   :: yminglob                          ! global minimum y coord of all nodes
!  REAL                                   :: zminglob                          ! global minimum z coord of all nodes
!  REAL                                   :: xmaxglob                          ! global max x coord of all nodes
!  REAL                                   :: ymaxglob                          ! global max y coord of all nodes
!  REAL                                   :: zmaxglob                          ! global max z coord of all nodes
!  REAL                                   :: xmin                              ! minimum x coord of all nodes
!  REAL                                   :: xmax                              ! maximum x coord of all nodes
!  REAL                                   :: ymin                              ! minimum y coord of all nodes
!  REAL                                   :: ymax                              ! maximum y coord of all nodes
!  REAL                                   :: zmin                              ! minimum z coord of all nodes
!  REAL                                   :: zmax                              ! maximum z coord of all nodes
!  REAL                                   :: nnx(3),nny(3),nnz(3)              ! periodic vectors
!  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
END TYPE

TYPE (tGeometry)                         :: GEO
!===================================================================================================================================


END MODULE MOD_Particle_Mesh_Vars
