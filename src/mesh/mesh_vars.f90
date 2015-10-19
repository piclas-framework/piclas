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
REAL              :: DeltaXi_NGeo
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
! PIC - for Newton localisation of particles in curved Elements
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: wBaryCL_NGeo(:)
REAL,ALLOCATABLE :: XiCL_NGeo(:)
REAL,ALLOCATABLE :: XCL_NGeo(:,:,:,:,:)
REAL,ALLOCATABLE :: dXCL_NGeo(:,:,:,:,:,:) !jacobi matrix of the mapping P\in NGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! surface vectors 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: NormVec(:,:,:,:)
REAL,ALLOCATABLE :: TangVec1(:,:,:,:)
REAL,ALLOCATABLE :: TangVec2(:,:,:,:)
REAL,ALLOCATABLE :: SurfElem(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: ElemToSide(:,:,:) ! SideID    = ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
                                         ! flip      = ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)
INTEGER,ALLOCATABLE :: SideToElem(:,:)   ! ElemID    = SideToElem(S2E_ELEM_ID,SideID)
                                         ! NB_ElemID = SideToElem(S2E_NB_ELEM_ID,SideID)
                                         ! locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
INTEGER,ALLOCATABLE :: SideToElem2(:,:)  ! ElemID    = SideToElem2(S2E2_ELEM_ID,SideID)  
                                         ! SideID    = SideToElem2(S2E2_SIDE_ID,SideID)  
                                         ! locSideID = SideToElem2(S2E2_LOC_SIDE_ID,SideID)
                                         ! flip      = SideToElem2(S2E2_FLIP,SideID)
INTEGER,ALLOCATABLE :: BC(:)             ! BCIndex   = BC(SideID), 1:nBCSides
INTEGER,ALLOCATABLE :: BoundaryType(:,:) ! BCType    = BoundaryType(BC(SideID),BC_TYPE)
                                         ! BCState   = BoundaryType(BC(SideID),BC_STATE)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: nGlobalElems=0      ! number of elements in mesh
INTEGER          :: nElems=0            ! number of local elements
INTEGER          :: offsetElem=0        ! for MPI, until now=0 Elems pointer array range: [offsetElem+1:offsetElem+nElems]
INTEGER          :: offsetSurfElem=0        !
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
INTEGER          :: nNodes=0            ! total number of nodes in mesh nElems*(NGeo+1)**3
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
  TYPE(tSidePtr),POINTER       :: Side(:)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             ! global side ID 
  INTEGER                      :: sideID          ! local side ID on Proc 
  INTEGER                      :: tmp 
  INTEGER                      :: NbProc 
  INTEGER                      :: BCindex         ! index in BoundaryType array! 
  INTEGER                      :: flip 
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
END TYPE tSide

TYPE tNode
  INTEGER                      :: ind=0         ! global unique node index
  REAL                         :: x(3)=0.
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)
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
!===================================================================================================================================
ALLOCATE(getNewSide)
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%connection)
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
INTEGER             :: iLocSide
!===================================================================================================================================
ALLOCATE(getNewElem)
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
INTEGER       :: iElem,iLocSide
!===================================================================================================================================
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem%Side)
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)
END SUBROUTINE deleteMeshPointer


END MODULE MOD_Mesh_Vars
