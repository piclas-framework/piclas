!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

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
INTEGER, ALLOCATABLE                     :: casematrix(:,:)   ! matrix to compute periodic cases
INTEGER                                  :: NbrOfCases        ! Number of periodic cases


! general: periodic sides have to be Cartesian
INTEGER,ALLOCATABLE :: SidePeriodicType(:)                    ! 1:nTotalSides, periodic type of side
                                                              ! 0 - normal or BC side
                                                              ! >0 type of periodic displacement
REAL,ALLOCATABLE    :: SidePeriodicDisplacement(:,:)          ! displacement vector

INTEGER,ALLOCATABLE :: PartElemToSide(:,:,:)                  ! extended list: 1:2,1:6,1:nTotalElems
                                                              ! ElemToSide: my geometry + halo
                                                              ! geometry + halo information


INTEGER,ALLOCATABLE :: PartSideToElem(:,:)                    ! extended list: 1:5,1:6,1:nTotalSides
                                                              ! SideToElem: my geometry + halo
                                                              ! geometry + halo information

INTEGER(KIND=8),ALLOCATABLE :: PartElemToElemGlob(:,:,:)          ! Mapping from ElemToElem
                                                              ! 1:4,1:6,1:nTotalElems
                                                              ! now in global-elem-ids !!!
INTEGER(KIND=4),ALLOCATABLE :: PartElemToElemAndSide(:,:,:)   ! Mapping from ElemToElem
                                                              ! 1:8,1:6,1:nTotalElems
                                                              ! [1]1:4 - MortarNeighborElemID
                                                              ! [1]5:8 -       Neighbor locSideID
                                                              ! [2]1:6 - locSideID
                                                              ! [3]    - nTotalElems
                                                              ! now in global-elem-ids !!!
INTEGER             :: nPartSides                             ! nPartSides - nSides+nPartPeriodicSides
INTEGER             :: nTotalSides                            ! total nb. of sides (my+halo)
INTEGER             :: nPartPeriodicSides                     ! total nb. of sides (my+halo)
INTEGER             :: nTotalElems                            ! total nb. of elems (my+halo)
INTEGER             :: nTotalNodes                            ! total nb. of nodes (my+halo)

INTEGER,ALLOCATABLE :: TracingBCInnerSides(:)                 ! number of local element boundary faces
                                                              ! used for tracing (connected to element)
INTEGER,ALLOCATABLE :: TracingBCTotalSides(:)                 ! total number of element boundary faces
                                                              ! used for tracing (loc faces + other
                                                              ! element faces that are possibly reached)
LOGICAL,ALLOCATABLE :: IsTracingBCElem(:)                     ! is an elem with BC sides for tracing
                                                              ! or BC in halo-eps distance to BC

LOGICAL,ALLOCATABLE :: IsLocalDepositionBCElem(:)             ! is an element where the deposition of a particle via a shape function
                                                              ! would result in the truncation of the shape function at the boundary.
TYPE tElemHaloInfo
  INTEGER,ALLOCATABLE            :: ElemHaloInfo(:)           !< Contains information regarding the halo region of each rank
                                                              !< ElemHaloInfo = -1            : element not in list
                                                              !<              = 0             : halo elements
                                                              !<              = 1 to PP_nElems: local elements
                                                              !< the default value of this variable is FALSE
END TYPE
TYPE(tElemHaloInfo),ALLOCATABLE      :: ElemHaloInfoProc(:)   ! ElemHaloInfo array for each rank

INTEGER,ALLOCATABLE :: ElemType(:)              !< Type of Element 1: only planar side, 2: one bilinear side 3. one curved side
LOGICAL,ALLOCATABLE :: ElemHasAuxBCs(:,:)
INTEGER             :: nTotalBCSides                          ! total number of BC sides (my+halo)
INTEGER             :: nTotalBCElems                          ! total number of bc elems (my+halo)
INTEGER,ALLOCATABLE :: PartBCSideList(:)                      ! mapping from SideID to BCSideID

REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: XiEtaZetaBasis     ! element local basis vector (linear elem)
REAL,ALLOCATABLE,DIMENSION(:,:)         :: slenXiEtaZetaBasis ! inverse of length of basis vector
REAL,ALLOCATABLE,DIMENSION(:)           :: ElemRadiusNGeo     ! radius of element
REAL,ALLOCATABLE,DIMENSION(:)           :: ElemRadius2NGeo    ! radius of element + 2% tolerance
INTEGER                                 :: RefMappingGuess    ! select guess for mapping into reference
                                                              ! element
                                                              ! 1 - Linear, cubical element
                                                              ! 2 - closest Gauss-Point
                                                              ! 3 - closest XCL-point
                                                              ! 4 - trivial guess - element origin
REAL                                    :: RefMappingEps      ! tolerance for Netwton to get xi from X
REAL                                    :: epsInCell          ! tolerance for eps for particle
                                                              ! inside of ref element
REAL,ALLOCATABLE                        :: epsOneCell(:)      ! tolerance for particle in
                                                              ! inside ref element 1+epsinCell

!LOGICAL                                 :: DoRefMapping      ! tracking by mapping particle into reference element
!-----------------------------------------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tFastInitBGM
  INTEGER                                :: nElem             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)        ! List of elements/physical cells in BGM cell
#ifdef MPI
  INTEGER, ALLOCATABLE                   :: ShapeProcs(:)     ! first Entry: Number of Shapeprocs,
                                                              ! following: ShapeProcs
  INTEGER, ALLOCATABLE                   :: PaddingProcs(:)   ! first Entry: Number of Paddingprocs,
                                                              ! following: PaddingProcs
  INTEGER, ALLOCATABLE                   :: SharedProcs(:)    ! first Entry: Number of Sharedprocs,
                                                              ! following: SharedProcs
  !INTEGER                                :: nBCSides         ! number BC sides in BGM cell
#endif
END TYPE

INTEGER                                  :: FIBGMCellPadding(1:3)

TYPE tNodeToElem
  INTEGER, ALLOCATABLE                   :: ElemID(:)
END TYPE

TYPE tGeometry
  REAL                                   :: xminglob                 ! global minimum x coord of all nodes
  REAL                                   :: yminglob                 ! global minimum y coord of all nodes
  REAL                                   :: zminglob                 ! global minimum z coord of all nodes
  REAL                                   :: xmaxglob                 ! global max x coord of all nodes
  REAL                                   :: ymaxglob                 ! global max y coord of all nodes
  REAL                                   :: zmaxglob                 ! global max z coord of all nodes
  REAL                                   :: xmin                     ! minimum x coord of all nodes
  REAL                                   :: xmax                     ! maximum x coord of all nodes
  REAL                                   :: ymin                     ! minimum y coord of all nodes
  REAL                                   :: ymax                     ! maximum y coord of all nodes
  REAL                                   :: zmin                     ! minimum z coord of all nodes
  REAL                                   :: zmax                     ! maximum z coord of all nodes
  ! periodic
  INTEGER                                :: nPeriodicVectors         ! Number of periodic Vectors
  REAL, ALLOCATABLE                      :: PeriodicVectors(:,:)     ! PeriodicVectors(1:3,1:nPeriodicVectors), 1:3=x,y,z
  REAL, ALLOCATABLE                      :: PeriodicVectorsLength(:) ! PeriodicVectorsLength(1:nPeriodicVectors)
  INTEGER,ALLOCATABLE                    :: DirPeriodicVectors(:)    ! direction of periodic vectors
  LOGICAL                                :: directions(3)            ! flag for direction
  ! required for cartesian BGM for desposition
  INTEGER, ALLOCATABLE                   :: PeriodicBGMVectors(:,:)  ! = periodic vectors in backgroundmesh coords
  ! FIBGM
  REAL                                   :: FIBGMdeltas(3)           ! size of background mesh cell for particle init
  REAL                                   :: FactorFIBGM(3)           ! scaling factor for FIBGM

  ! caution, possible pointer
  TYPE (tFastInitBGM),ALLOCATABLE        :: FIBGM(:,:,:)  !        =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: FIBGMimin                         ! smallest index of FastInitBGM (x)
  INTEGER                                :: FIBGMimax                         ! biggest index of FastInitBGM (x)
  INTEGER                                :: FIBGMjmin                         ! smallest index of FastInitBGM (y)
  INTEGER                                :: FIBGMjmax                         ! biggest index of FastInitBGM (y)
  INTEGER                                :: FIBGMkmin                         ! smallest index of FastInitBGM (z)
  INTEGER                                :: FIBGMkmax                         ! biggest index of FastInitBGM (z)

  TYPE (tFastInitBGM),ALLOCATABLE        :: TFIBGM(:,:,:)  !       =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: TFIBGMimin                        ! smallest index of FastInitBGM (x)
  INTEGER                                :: TFIBGMimax                        ! biggest index of FastInitBGM (x)
  INTEGER                                :: TFIBGMjmin                        ! smallest index of FastInitBGM (y)
  INTEGER                                :: TFIBGMjmax                        ! biggest index of FastInitBGM (y)
  INTEGER                                :: TFIBGMkmin                        ! smallest index of FastInitBGM (z)
  INTEGER                                :: TFIBGMkmax                        ! biggest index of FastInitBGM (z)

  INTEGER,ALLOCATABLE                    :: ElemToFIBGM(:,:)                  ! range of FIGMB cells per element
                                                                              ! 1:6,1:nTotalElems, xmin,max,yminmax,...
  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
  REAL, ALLOCATABLE                      :: CharLength(:)                     ! Characteristic length for each cell: L=V^(1/3)
  REAL                                   :: MeshVolume                        ! Total Volume of mesh
  REAL                                   :: LocalVolume                       ! Volume of proc
  REAL, ALLOCATABLE                      :: DeltaEvMPF(:)                     ! Energy difference due to particle merge
  INTEGER, ALLOCATABLE                   :: ElemToRegion(:)                   ! ElemToRegion(1:nElems)

  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
  INTEGER, ALLOCATABLE                   :: ElemToNodeID(:,:)                 ! ElemToNodeID(1:nElemNodes,1:nElems)
  !INTEGER, ALLOCATABLE                   :: ElemToNodeIDGlobal(:,:)           ! ElemToNodeID(1:nElemNodes,1:nElems)
  INTEGER, ALLOCATABLE                   :: ElemSideNodeID(:,:,:)             ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
                                                                              ! From element sides to node IDs
  INTEGER, ALLOCATABLE                   :: ElemsOnNode(:)                    ! number of elements on proc global node
  TYPE(tNodeToElem), ALLOCATABLE         :: NodeToElem(:)                     ! mapping of elements per node
  INTEGER, ALLOCATABLE                   :: NeighNodesOnNode(:)               ! number of neighbour nodes on proc global node
  TYPE(tNodeToElem), ALLOCATABLE         :: NodeToNeighNode(:)                ! mapping of neighbour nodes per proc global node
  INTEGER, ALLOCATABLE                   :: NumNeighborElems(:)               ! number of elements per element
  TYPE(tNodeToElem), ALLOCATABLE         :: ElemToNeighElems(:)               ! mapping of neighbour elements per element
                                                                              ! From element sides to node IDs
  INTEGER, ALLOCATABLE                   :: PeriodicElemSide(:,:)             ! 0=not periodic side, others=PeriodicVectorsNum
  LOGICAL, ALLOCATABLE                   :: ConcaveElemSide(:,:)              ! Whether LocalSide of Element is concave side
  REAL, ALLOCATABLE                      :: NodeCoords(:,:)                   ! Node Coordinates (1:nDim,1:nNodes)
END TYPE

TYPE (tGeometry)                         :: GEO

INTEGER                                  :: WeirdElems                        ! Number of Weird Elements (=Elements which are folded
                                                                              ! into themselves)
LOGICAL                                  :: FindNeighbourElems=.FALSE.        ! Flag defining if mapping for neighbour elements
                                                                              ! is build via nodes

TYPE tBCElem
  INTEGER                                :: nInnerSides                       ! Number of BC-Sides of Element
  INTEGER                                :: lastSide                          ! total number of BC-Sides in eps-vicinity of element
  INTEGER, ALLOCATABLE                   :: BCSideID(:)                       ! List of elements in BGM cell
  REAL,ALLOCATABLE                       :: ElemToSideDistance(:)             ! stores the distance between each element and the
                                                                              ! sides associated with this element
END TYPE

TYPE (tBCElem),ALLOCATABLE               :: BCElem(:)

INTEGER                                  :: NbrOfRegions      ! Nbr of regions to be mapped to Elems
REAL, ALLOCATABLE                        :: RegionBounds(:,:) ! RegionBounds ((xmin,xmax,ymin,...)|1:NbrOfRegions)
LOGICAL,ALLOCATABLE                      :: isTracingTrouble(:)
REAL,ALLOCATABLE                         :: ElemTolerance(:)
INTEGER, ALLOCATABLE                     :: ElemToGlobalElemID(:)  ! mapping form local-elemid to global-id
!===================================================================================================================================


END MODULE MOD_Particle_Mesh_Vars
