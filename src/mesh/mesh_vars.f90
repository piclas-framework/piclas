#include "boltzplatz.h"

MODULE MOD_Mesh_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! basis
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER           :: NGeo                        ! polynomial degree of geometric transformation
REAL,ALLOCATABLE  :: Xi_NGeo(:)                  ! 1D equidistant point positions for curved elements (during readin)
REAL,ALLOCATABLE  :: Vdm_CLN_GaussN(:,:)
REAL,ALLOCATABLE  :: Vdm_CLNGeo_CLN(:,:)
REAL,ALLOCATABLE  :: Vdm_CLNGeo_GaussN(:,:)  
REAL,ALLOCATABLE  :: Vdm_NGeo_CLNGeo(:,:)  
REAL,ALLOCATABLE  :: DCL_NGeo(:,:)  
REAL,ALLOCATABLE  :: DCL_N(:,:)  
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Elem_xGP(:,:,:,:,:)   ! XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: BCFace_xGP(:,:,:,:)   ! XYZ positions (first index 1:3) of the Boundary Face Gauss Point
!-----------------------------------------------------------------------------------------------------------------------------------
! Metrics on GaussPoints 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Metrics_fTilde(:,:,:,:,:) ! Metric Terms (first indices 3) on each GaussPoint
REAL,ALLOCATABLE :: Metrics_gTilde(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_hTilde(:,:,:,:,:)
REAL,ALLOCATABLE :: sJ(:,:,:,:)               ! 1/DetJac for each Gauss Point
!-----------------------------------------------------------------------------------------------------------------------------------
! surface vectors 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: NormVec(:,:,:,:)
REAL,ALLOCATABLE :: TangVec1(:,:,:,:)
REAL,ALLOCATABLE :: TangVec2(:,:,:,:)
REAL,ALLOCATABLE :: SurfElem(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: ElemToSide(:,:,:)
INTEGER,ALLOCATABLE :: SideToElem(:,:)
INTEGER,ALLOCATABLE :: SideToElem2(:,:)
INTEGER,ALLOCATABLE :: BC(:)
INTEGER,ALLOCATABLE :: BoundaryType(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: nGlobalElems=0      ! number of elements in mesh
INTEGER          :: nElems=0            ! number of local elements
INTEGER          :: offsetElem=0        ! for MPI, until now=0 Elems pointer array range: [offsetElem+1:offsetElem+nElems]
INTEGER          :: nSides=0            ! =nInnerSides+nBCSides+nMPISides
INTEGER          :: nInnerSides=0       ! InnerSide index range: sideID \in [nBCSides+1:nBCSides+nInnerSides]
INTEGER          :: nBCSides=0          ! BCSide index range: sideID \in [1:nBCSides]
INTEGER          :: nMPISides=0
INTEGER          :: nMPISides_MINE=0
INTEGER          :: nMPISides_YOUR=0
INTEGER          :: SideID_minus_lower  ! lower side ID of array U_minus/GradUx_minus...
INTEGER          :: SideID_minus_upper  ! upper side ID of array U_minus/GradUx_minus...
INTEGER          :: SideID_plus_lower   ! lower side ID of array U_plus/GradUx_plus...
INTEGER          :: SideID_plus_upper   ! upper side ID of array U_plus/GradUx_plus...
INTEGER          :: nNodes=0            ! SIZE of Nodes pointer array, number of unique nodes
INTEGER          :: nBCs=0              ! number of BCs in mesh
INTEGER          :: nUserBCs=0          ! number of BC in inifile
INTEGER          :: MeshType            !for DEBUGMESH
INTEGER          :: MaxBCState
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE   :: BoundaryName(:)
CHARACTER(LEN=255)               :: MeshFile        ! name of hdf5 meshfile (write with ending .h5!)
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: useCurveds
LOGICAL          :: CrossProductMetrics=.FALSE.
!-----------------------------------------------------------------------------------------------------------------------------------
! PoyntingVectorIntegral variables
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE    :: Face_xGP(:,:,:,:)     ! XYZ positions (first index 1:3) of the Face Gauss Point
INTEGER             :: nPoyntingIntSides=0   ! Sides for the calculation of the poynting vector
LOGICAL,ALLOCATABLE :: isPoyntingIntSide(:)  ! number of all PoyntingInt sides
INTEGER,ALLOCATABLE :: whichPoyntingPlane(:) ! number of plane used for calculation of poynting vector
!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
REAL            :: ParticleMPIWeight
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES 

TYPE tNodePtr
  TYPE(tNode),POINTER          :: np                     ! node pointer
END TYPE tNodePtr

TYPE tSidePtr
  TYPE(tSide),POINTER          :: sp                     ! side pointer
END TYPE tSidePtr

TYPE tElemPtr
  TYPE(tElem),POINTER          :: ep                     ! Local element pointer
END TYPE tElemPtr

TYPE tElem
  INTEGER                      :: ind             ! global element index
  INTEGER                      :: Type            ! element type (linear/bilinear/curved)
  INTEGER                      :: Zone
  TYPE(tNodePtr),POINTER       :: Node(:)
  TYPE(tSidePtr),POINTER       :: Side(:)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             ! global side ID 
  INTEGER                      :: sideID          ! local side ID on Proc 
  INTEGER                      :: tmp 
  INTEGER                      :: NbProc 
  INTEGER                      :: BCindex         ! index in BoundaryType array! 
  INTEGER                      :: flip 
  INTEGER                      :: nCurvedNodes 
  TYPE(tNodePtr),POINTER       :: Node(:)
  TYPE(tNodePtr),POINTER       :: CurvedNode(:)
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
END TYPE tSide

TYPE tNode
  INTEGER                      :: ind=0         ! global unique node index
  REAL                         :: x(3)=0.
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)
TYPE(tNodePtr),POINTER         :: Nodes(:)
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
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
!  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSide),POINTER :: getNewSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNode
!===================================================================================================================================
ALLOCATE(getNewSide)
ALLOCATE(getNewSide%Node(4))
DO iNode=1,4
  NULLIFY(getNewSide%Node(iNode)%np)
END DO
NULLIFY(getNewSide%CurvedNode)
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%connection)
getNewSide%nCurvedNodes=0
getNewSide%sideID=0
getNewSide%ind=0
getNewSide%tmp=0
getNewSide%NbProc=-1
getNewSide%BCindex=0
getNewSide%flip=0
END FUNCTION GETNEWSIDE

FUNCTION GETNEWELEM()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tElem),POINTER :: getNewElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
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


SUBROUTINE createSides(Elem)
!===================================================================================================================================
! if element nodes already assigned, create Sides using CGNS standard
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER :: Elem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!side 1
Elem%Side(1)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(1)%sp%Node(2)%np=>Elem%Node(4)%np
Elem%Side(1)%sp%Node(3)%np=>Elem%Node(3)%np
Elem%Side(1)%sp%Node(4)%np=>Elem%Node(2)%np
!side 2                                    
Elem%Side(2)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(2)%sp%Node(2)%np=>Elem%Node(2)%np
Elem%Side(2)%sp%Node(3)%np=>Elem%Node(6)%np
Elem%Side(2)%sp%Node(4)%np=>Elem%Node(5)%np
!side 3                                    
Elem%Side(3)%sp%Node(1)%np=>Elem%Node(2)%np
Elem%Side(3)%sp%Node(2)%np=>Elem%Node(3)%np
Elem%Side(3)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(3)%sp%Node(4)%np=>Elem%Node(6)%np
!side 4                                    
Elem%Side(4)%sp%Node(1)%np=>Elem%Node(3)%np
Elem%Side(4)%sp%Node(2)%np=>Elem%Node(4)%np
Elem%Side(4)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(4)%sp%Node(4)%np=>Elem%Node(7)%np
!side 5                                    
Elem%Side(5)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(5)%sp%Node(2)%np=>Elem%Node(5)%np
Elem%Side(5)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(5)%sp%Node(4)%np=>Elem%Node(4)%np
!side 6                                                
Elem%Side(6)%sp%Node(1)%np=>Elem%Node(5)%np
Elem%Side(6)%sp%Node(2)%np=>Elem%Node(6)%np
Elem%Side(6)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(6)%sp%Node(4)%np=>Elem%Node(8)%np
END SUBROUTINE createSides


SUBROUTINE deleteMeshPointer()
!===================================================================================================================================
! Deallocates all pointers used for the mesh readin
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: FirstElemInd,LastElemInd
INTEGER       :: iElem,iLocSide,iNode,nAssocNodes
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
    DO iNode=1,4
      NULLIFY(aSide%Node(iNode)%np)
    END DO
    DO iNode=1,aSide%nCurvedNodes
      NULLIFY(aSide%CurvedNode(iNode)%np)
    END DO
    DEALLOCATE(aSide%Node)
    IF(ASSOCIATED(aSide%CurvedNode))DEALLOCATE(aSide%CurvedNode)
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem%Side)
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)
nAssocNodes=0
DO iNode=1,nNodes
  IF(ASSOCIATED(Nodes(iNode)%np))THEN
    DEALLOCATE(Nodes(iNode)%np)
    nAssocNodes=nAssocNodes+1
  END IF
END DO
!WRITE(*,*)'DEBUG,nAssocNodes',nAssocNodes
DEALLOCATE(Nodes)
END SUBROUTINE deleteMeshPointer


END MODULE MOD_Mesh_Vars
