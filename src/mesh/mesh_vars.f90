!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Mesh_Vars
!===================================================================================================================================
!> Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! SwapMesh
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: DoSwapMesh                   !< flag for SwapMesh routines
CHARACTER(LEN=255):: SwapMeshExePath              !< path to swapmesh binary
INTEGER           :: SwapMeshLevel                !< 0: initial grid, 1: first swap mesh, 2: second swap mesh
!-----------------------------------------------------------------------------------------------------------------------------------
! basis
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: NGeo                        !< polynomial degree of geometric transformation
INTEGER          :: NGeoRef                     !< polynomial degree of geometric transformation
INTEGER          :: NGeoElevated                !< polynomial degree of elevated geometric transformation
REAL,ALLOCATABLE :: Xi_NGeo(:)                  !< 1D equidistant point positions for curved elements (during readin)
REAL             :: DeltaXi_NGeo
REAL,ALLOCATABLE :: Vdm_EQ_N(:,:)               !< Vandermonde mapping from equidistant (visu) to NodeType node set
REAL,ALLOCATABLE :: Vdm_N_EQ(:,:)               !< Vandermonde mapping from NodeType to equidistant (visu) node set
REAL,ALLOCATABLE :: Vdm_GL_N(:,:)               !< Vandermonde mapping from Gauss-Lobatto (analyze) to NodeType node set
REAL,ALLOCATABLE :: Vdm_N_GL(:,:)               !< Vandermonde mapping from NodeType to Gauss-Lobatto (analyze) node set
! check if these arrays are still used
REAL,ALLOCATABLE :: Vdm_CLN_GaussN(:,:)
REAL,ALLOCATABLE :: Vdm_CLNGeo_CLN(:,:)
REAL,ALLOCATABLE :: Vdm_CLNGeo_GaussN(:,:)
REAL,ALLOCATABLE :: Vdm_NGeo_CLNGeo(:,:)
REAL,ALLOCATABLE :: DCL_NGeo(:,:)
REAL,ALLOCATABLE :: DCL_N(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! will be used in the future
REAL,ALLOCATABLE,TARGET :: NodeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeo) of element interpolation points from meshfile
REAL,ALLOCATABLE,TARGET :: Elem_xGP(:,:,:,:,:)          !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: Face_xGP(:,:,:,:)            !< XYZ positions (first index 1:3) of the Boundary Face Gauss Point
REAL,DIMENSION(6):: xyzMinMax                    !< from Face_xGP points determined maximum domain extension (min/max of domain)
LOGICAL          :: GetMeshMinMaxBoundariesIsDone!< don't call twice the calculation of xyzMinMax
REAL,ALLOCATABLE,DIMENSION(:,:):: ElemBaryNGeo   !< element local basis: origin
!-----------------------------------------------------------------------------------------------------------------------------------
! Metrics on GaussPoints
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Metrics_fTilde(:,:,:,:,:) !< Metric Terms (first indices 3) on each GaussPoint
REAL,ALLOCATABLE :: Metrics_gTilde(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_hTilde(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET :: sJ(:,:,:,:)               !< 1/DetJac for each Gauss Point
!-----------------------------------------------------------------------------------------------------------------------------------
! PIC - for Newton localisation of particles in curved Elements
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE    :: wBaryCL_NGeo(:)
!< #ifdef PARTICLES
REAL,ALLOCATABLE    :: wBaryCL_NGeo1(:)
REAL,ALLOCATABLE    :: XiCL_NGeo1(:)
REAL,ALLOCATABLE    :: Vdm_CLNGeo1_CLNGeo(:,:)
!< #endif /*PARTICLES*/
REAL,ALLOCATABLE        :: XiCL_NGeo(:)
REAL,ALLOCATABLE,TARGET :: XCL_NGeo(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET :: dXCL_NGeo(:,:,:,:,:,:) !jacobi matrix of the mapping P\in NGeo
REAL,ALLOCATABLE        :: dXCL_N(:,:,:,:,:,:) !jacobi matrix of the mapping P\in NGeo
REAL,ALLOCATABLE        :: detJac_Ref(:,:,:,:,:)      !< determinant of the mesh Jacobian for each Gauss point at degree 3*NGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! surface vectors
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: NormVec(:,:,:,:)           !< normal vector for each side       (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec1(:,:,:,:)          !< tangential vector 1 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec2(:,:,:,:)          !< tangential vector 3 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: SurfElem(:,:,:)            !< surface area for each side        (    0:N,0:N,nSides)
REAL,ALLOCATABLE :: Ja_Face(:,:,:,:,:)         !< surface  metrics for each side
REAL,ALLOCATABLE :: nVecLoc(:,:,:,:,:)         !< element local normal vector       (1:3,0:N,0:N,1:6,1:nElems)
REAL,ALLOCATABLE :: SurfLoc(:,:,:,:)           !< element local surface element     (    0:N,0:N,1:6,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! mapping from GaussPoints to Side or Neighbor Volume
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: VolToSideA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: VolToSideIJKA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: VolToSide2A(:,:,:,:,:)
INTEGER,ALLOCATABLE :: CGNS_VolToSideA(:,:,:,:,:)
INTEGER,ALLOCATABLE :: CGNS_SideToVol2A(:,:,:,:)
INTEGER,ALLOCATABLE :: SideToVolA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: SideToVol2A(:,:,:,:,:)
INTEGER,ALLOCATABLE :: FS2M(:,:,:,:)     !< flip slave side to master and reverse
!-----------------------------------------------------------------------------------------------------------------------------------
! mapping from element to sides and sides to element
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: ElemInfo(:,:)           !< array containing the node and side connectivity of the elments as stored in the
                                               !< mesh file
INTEGER,ALLOCATABLE :: SideInfo(:,:)           !< array containing the connectivity, flip,... of the sides as stored in the
                                               !< mesh file
INTEGER,ALLOCATABLE :: SideToNonUniqueGlobalSide(:,:)     !< maps the local SideIDs to global SideIDs (for parallel HDG load balance currently)
INTEGER,ALLOCATABLE :: ElemToSide(:,:,:) !< SideID    = ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
                                         !< flip      = ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)
INTEGER,ALLOCATABLE :: SideToElem(:,:)   !< ElemID    = SideToElem(S2E_ELEM_ID,SideID)
                                         !< NB_ElemID = SideToElem(S2E_NB_ELEM_ID,SideID)
                                         !< locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
INTEGER,ALLOCATABLE :: BC(:)             !< BCIndex   = BC(SideID), 1:nBCSides
INTEGER,ALLOCATABLE :: GlobalUniqueSideID(:) !< SideInfo(SIDE_ID,iSide) = GlobalUniqueSideIDC(SideID), 1:nSides
INTEGER,ALLOCATABLE :: BoundaryType(:,:) !< BCType    = BoundaryType(BC(SideID),BC_TYPE)
                                         !< BCState   = BoundaryType(BC(SideID),BC_STATE)
INTEGER,ALLOCATABLE :: AnalyzeSide(:)    !< Marks, wheter a side belongs to a group of analyze sides (e.g. to a BC group)
                                         !< SurfIndex = AnalyzeSide(SideID), 1:nSides

INTEGER,PARAMETER :: NormalDirs(6) = (/ 3 , 2 , 1 , 2 , 1 , 3 /) !< normal vector direction for element local side
INTEGER,PARAMETER :: TangDirs(6)   = (/ 1 , 3 , 2 , 3 , 2 , 1 /) !< first tangential vector direction for element local side
REAL   ,PARAMETER :: NormalSigns(6)= (/-1.,-1., 1., 1.,-1., 1./) !< normal vector sign for element local side
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: nGlobalElems=0          !< number of elements in mesh
INTEGER          :: nElems=0                !< number of local elements
INTEGER          :: offsetElem=0            !< for MPI, until now=0 Elems pointer array range: [offsetElem+1:offsetElem+nElems]
INTEGER          :: nSides=0                !< =nInnerSides+nBCSides+nMPISides
INTEGER          :: nUniqueSides=0          !< =uniquesides for hdg output
INTEGER          :: nGlobalUniqueSides=0    !< =uniquesides for hdg output
INTEGER          :: nGlobalUniqueSidesFromMesh=0 !< =uniquesides read from mesh file
INTEGER          :: nGlobalMortarSides=0    !< global number of big mortar sides
INTEGER          :: offsetSide=0            !< for MPI, until now=0  Sides pointer array range
INTEGER          :: nSidesMaster=0          !< =sideIDMaster
INTEGER          :: nSidesSlave=0           !< =nInnerSides+nBCSides+nMPISides
INTEGER          :: nInnerSides=0           !< InnerSide index range: sideID [nBCSides+1:nBCSides+nInnerSides]
INTEGER          :: nBCSides=0              !< BCSide index range: sideID [1:nBCSides]
INTEGER          :: nAnalyzeSides=0         !< marker for each side (BC,analyze flag, periodic,...)
INTEGER          :: nMPISides=0             !< number of MPI sides in mesh
INTEGER          :: nMPISides_MINE=0        !< number of MINE MPI sides (on local processor)
INTEGER          :: nMPISides_YOUR=0        !< number of YOUR MPI sides (on neighbour processors)
INTEGER          :: nNodes=0                !< SIZE of Nodes pointer array, number of unique nodes
INTEGER          :: nBCs=0                  !< number of BCs in mesh
INTEGER          :: nUserBCs=0              !< number of BC in inifile
LOGICAL          :: ChangedPeriodicBC       !< is set true if BCs are changed from periodic to non-periodic (default is false)
!----------------------------------------------------------------------------------------------------------------------------------
! Define index ranges for all sides in consecutive order for easier access
INTEGER             :: firstBCSide             !< First SideID of BCs (in general 1)
INTEGER             :: firstMortarInnerSide    !< First SideID of Mortars (in general nBCSides+1)
INTEGER             :: firstInnerSide          !< First SideID of inner sides
INTEGER             :: firstMPISide_MINE       !< First SideID of MINE MPI sides (on local processor)
INTEGER             :: firstMPISide_YOUR       !< First SideID of YOUR MPI sides (on neighbour processor)
INTEGER             :: firstMortarMPISide      !< First SideID of Mortar MPI sides
INTEGER             :: lastBCSide              !< Last  SideID of BCs (in general nBCSides)
INTEGER             :: lastMortarInnerSide     !< Last  SideID of Mortars (in general nBCSides+nMortars)
INTEGER             :: lastInnerSide           !< Last  SideID of inner sides
INTEGER             :: lastMPISide_MINE        !< Last  SideID of MINE MPI sides (on local processor)
INTEGER             :: lastMPISide_YOUR        !< Last  SideID of YOUR MPI sides (on neighbour processor)
INTEGER             :: lastMortarMPISide       !< Last  SideID of Mortar MPI sides (in general nSides)
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nMortarSides=0          !< total number of mortar sides
INTEGER             :: nMortarInnerSides=0     !< number of inner mortar sides
INTEGER             :: nMortarMPISides=0       !< number of mortar MPI sides
INTEGER,ALLOCATABLE :: MortarType(:,:)         !< Side Info about mortars, [1:2,1:nSides], Type of mortar [1] :
                                               !< =-1: conforming side not belonging to mortar
                                               !< =0: small mortar side belonging to big side ,
                                               !< =-10: neighbor of small mortar side (exists only if its an MPI side, too)
                                               !< =1: bigside type 1-4, 2: bigside type 1-2 eta, 3: bigside type 1-2 xi
                                               !< [2] position index in mortarInfo list
INTEGER,ALLOCATABLE :: MortarInfo(:,:,:)       !< 1:2,1:4,1:nMortarSides: [1] nbSideID / flip, [2] max 4 mortar sides, [3] sides
INTEGER,ALLOCATABLE :: MortarSlave2MasterInfo(:) !< 1:nSides: map of slave mortar sides to belonging master mortar sides
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8),ALLOCATABLE     :: ElemToElemGlob(:,:,:)             !< mapping from element to neighbor element in global ids
                                                                     !< [1:4] (mortar) neighbors
                                                                     !< [1:6] local sides
                                                                     !< [OffSetElem+1:OffsetElem+PP_nElems]
INTEGER(KIND=8),ALLOCATABLE     :: ElemGlobalID(:)                   !< global element id of each element
INTEGER(KIND=8),ALLOCATABLE     :: myInvisibleRank(:)                !< global proc ID which the current proc cannot see (particle
                                                                     !< communication)
INTEGER(KIND=8),ALLOCATABLE     :: LostRotPeriodicSides(:)           !< Number of lost sides during rotational periodic search
!LOGICAL                         :: RotPeriodicReBuild                !< Force re-building of mapping (might already exist)
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE   :: BoundaryName(:)
CHARACTER(LEN=255)               :: MeshFile        !< name of hdf5 meshfile (write with ending .h5!)
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: useCurveds
LOGICAL          :: CrossProductMetrics=.FALSE.
!-----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES

TYPE tNodePtr
  TYPE(tNode),POINTER          :: np                     !< node pointer
END TYPE tNodePtr

TYPE tSidePtr
  TYPE(tSide),POINTER          :: sp              !< side pointer
END TYPE tSidePtr

TYPE tElemPtr
  TYPE(tElem),POINTER          :: ep              !< Local element pointer
END TYPE tElemPtr

TYPE tElem
  INTEGER                      :: ind             !< global element index
  INTEGER                      :: Type            !< element type (linear/bilinear/curved)
  INTEGER                      :: Zone
  TYPE(tNodePtr),POINTER       :: Node(:)
  TYPE(tSidePtr),POINTER       :: Side(:)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             !< global side ID
  INTEGER                      :: sideID          !< local side ID on Proc
  INTEGER                      :: tmp
  INTEGER                      :: NbProc
  INTEGER                      :: BCindex         !< index in BoundaryType array!
  INTEGER                      :: flip=-999
#ifdef PARTICLES
  LOGICAL                      :: InnerBCOutput   !< Logical if proc writes InnerBC information
  INTEGER                      :: BC_Alpha        !< inital value for periodic displacement before mapping in pos. bc-index range
#endif /*PARTICLES*/
  INTEGER                      :: nMortars        !< number of slave mortar sides associated with master mortar
  INTEGER                      :: MortarType      !< type of mortar from mesh file: =0: conforming side or small side of bigside
                                                  !< =1 : big side with 1-4 =2: big side with 1-2 in eta, =3: bigSide with 1-2 in xi
                                                  !< =-10: connected neighbor side of small mortar side
  TYPE(tSidePtr),POINTER       :: MortarSide(:)   !< array of side pointers to slave mortar sides
  TYPE(tNodePtr),POINTER       :: Node(:)
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
END TYPE tSide

TYPE tNode
  INTEGER                      :: NodeID=0        !< local proc specific node index
  INTEGER                      :: ind=0           !< global unique node index
  REAL                         :: x(3)=0.
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)
TYPE(tNodePtr),POINTER         :: Nodes(:)
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: MeshInitIsDone =.FALSE.
!===================================================================================================================================

INTERFACE getNewSide
  MODULE PROCEDURE getNewSide
END INTERFACE

INTERFACE getNewElem
  MODULE PROCEDURE getNewElem
END INTERFACE

INTERFACE deleteMeshPointer
  MODULE PROCEDURE deleteMeshPointer
END INTERFACE

CONTAINS

FUNCTION GETNEWSIDE()
!===================================================================================================================================
!<
!===================================================================================================================================
!< MODULES
!< IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!< INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!< OUTPUT VARIABLES
TYPE(tSide),POINTER :: getNewSide
!-----------------------------------------------------------------------------------------------------------------------------------
!< LOCAL VARIABLES
INTEGER             :: iNode
!===================================================================================================================================
ALLOCATE(getNewSide)
ALLOCATE(getNewSide%Node(4))
DO iNode=1,4
  NULLIFY(getNewSide%Node(iNode)%np)
END DO
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%MortarSide)
NULLIFY(getNewSide%connection)
getNewSide%sideID=0
getNewSide%ind=0
getNewSide%tmp=0
getNewSide%NbProc=-1
getNewSide%BCindex=0
getNewSide%flip=0
getNewSide%nMortars=0
getNewSide%MortarType=0
END FUNCTION GETNEWSIDE

FUNCTION GETNEWELEM()
!===================================================================================================================================
!<
!===================================================================================================================================
!< MODULES
!< IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!< INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!< OUTPUT VARIABLES
TYPE(tElem),POINTER :: getNewElem
!-----------------------------------------------------------------------------------------------------------------------------------
!< LOCAL VARIABLES
INTEGER             :: iNode,iLocSide
!===================================================================================================================================
ALLOCATE(getNewElem)
ALLOCATE(getNewElem%Node(8))
DO iNode=1,8
  NULLIFY(getNewElem%Node(iNode)%np)
END DO
ALLOCATE(getNewElem%Side(6))
DO iLocSide=1,6
  getNewElem%Side(iLocSide)%sp=>getNewSide()
END DO
getNewElem%ind=0
getNewElem%Zone=0
getNewElem%Type=0
END FUNCTION GETNEWELEM


SUBROUTINE deleteMeshPointer()
!===================================================================================================================================
!> Deallocates all pointers used for the mesh readin
!===================================================================================================================================
!< MODULES
!< IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: FirstElemInd,LastElemInd
INTEGER       :: iElem,iLocSide,iNbLocSide
INTEGER       :: iMortar,iNode
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!===================================================================================================================================
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems

DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iNode=1,8
    NULLIFY(aElem%Node(iNode)%np)
  END DO
  DEALLOCATE(aElem%Node)
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    ! Free nodes
    DO iNode=1,4
      NULLIFY(aSide%Node(iNode)%np)
    END DO
    DEALLOCATE(aSide%Node)
    ! Free mortar sides
    DO iMortar=1,aSide%nMortars
      ! Free MPI connection
      IF (ASSOCIATED(aSide%MortarSide(iMortar)%sp%connection) .AND. aSide%MortarSide(iMortar)%sp%NbProc.NE.-1) THEN
        ! Free the connected elem
        DO iNbLocSide=1,6
          DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Elem%Side(iNbLocSide)%sp%Node)
          DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Elem%Side(iNbLocSide)%sp)
        END DO
        DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Elem%Node)
        DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Elem%Side)
        DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Elem)
        ! Free the connected size
        DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection%Node)
        DEALLOCATE(aSide%MortarSide(iMortar)%sp%connection)
      END IF
      DEALLOCATE(aSide%MortarSide(iMortar)%sp%Node)
      DEALLOCATE(aSide%MortarSide(iMortar)%sp)
    END DO
    IF(ASSOCIATED(aSide%MortarSide)) DEALLOCATE(aSide%MortarSide)
    ! Free MPI connection
    IF (ASSOCIATED(aSide%connection) .AND. aSide%NbProc.NE.-1) THEN
      ! Free the connected elem
      DO iNbLocSide=1,6
        DEALLOCATE(aSide%connection%Elem%Side(iNbLocSide)%sp%Node)
        DEALLOCATE(aSide%connection%Elem%Side(iNbLocSide)%sp)
      END DO
      DEALLOCATE(aSide%connection%Elem%Node)
      DEALLOCATE(aSide%connection%Elem%Side)
      DEALLOCATE(aSide%connection%Elem)
      ! Free the connected size
      DEALLOCATE(aSide%connection%Node)
      DEALLOCATE(aSide%connection)
    END IF
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem%Side)
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)

! Free the node pointer
DO iNode=1,nNodes
  IF(ASSOCIATED(Nodes(iNode)%np))THEN
    DEALLOCATE(Nodes(iNode)%np)
  END IF
END DO
!DEALLOCATE(Nodes)

END SUBROUTINE deleteMeshPointer

END MODULE MOD_Mesh_Vars
