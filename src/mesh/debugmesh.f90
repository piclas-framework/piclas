#include "boltzplatz.h"

MODULE MOD_DebugMesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDebugMesh
  MODULE PROCEDURE WriteDebugMesh
END INTERFACE

PUBLIC::WriteDebugMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDebugMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars, ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_Mesh_Vars,ONLY:nElems,Elem_xGP,ElemToSide
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_Tecplot,            ONLY:writeDataToTecplot
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,iLocSide
CHARACTER(LEN=32) :: VarNames(4)
REAL,ALLOCATABLE  :: debugVisu(:,:,:,:,:)
REAL,ALLOCATABLE  :: X_NVisu(:,:,:,:,:)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH...'
! WRITE Debugmesh.dat 
ALLOCATE(X_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))

DO iElem=1,nElems
  CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_Nvisu,Elem_xGP(:,:,:,:,iElem),X_NVisu(:,:,:,:,iElem))
END DO

VarNames(1)='ElemID'
VarNames(2)='SideID'
VarNames(3)='FLIP'
VarNames(4)='iLocSide'
ALLOCATE(debugVisu(4,0:NVisu,0:NVisu,0:NVisu,nElems))
debugVisu=-1.
DO iElem=1,nElems
  debugVisu(1,:,:,:,iElem)=REAL(iElem)
  DO iLocSide=1,6
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      debugVisu(2,0,:,:,iElem)=REAL(ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem))
      debugVisu(3,0,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_MINUS,iElem))
      debugVisu(4,0,:,:,iElem)=REAL(iLocSide)
    CASE(XI_PLUS)
      debugVisu(2,NVisu,:,:,iElem)=REAL(ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem))
      debugVisu(3,NVisu,:,:,iElem)=REAL(ElemToSide(E2S_FLIP,XI_PLUS,iElem))
      debugVisu(4,NVisu,:,:,iElem)=REAL(iLocSide)
    CASE(ETA_MINUS)
      debugVisu(2,:,0,:,iElem)=REAL(ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem))
      debugVisu(3,:,0,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_MINUS,iElem))
      debugVisu(4,:,0,:,iElem)=REAL(iLocSide)
    CASE(ETA_PLUS)
      debugVisu(2,:,NVisu,:,iElem)=REAL(ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem))
      debugVisu(3,:,NVisu,:,iElem)=REAL(ElemToSide(E2S_FLIP,ETA_PLUS,iElem))
      debugVisu(4,:,NVisu,:,iElem)=REAL(iLocSide)
    CASE(ZETA_MINUS)
      debugVisu(2,:,:,0,iElem)=REAL(ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))
      debugVisu(3,:,:,0,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem))
      debugVisu(4,:,:,0,iElem)=REAL(iLocSide)
    CASE(ZETA_PLUS)
      debugVisu(2,:,:,NVisu,iElem)=REAL(ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem))
      debugVisu(3,:,:,NVisu,iElem)=REAL(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem))
      debugVisu(4,:,:,NVisu,iElem)=REAL(iLocSide)
    END SELECT
  END DO
END DO
CALL  WriteDataToTecplot(NVisu,nElems,4,0,VarNames,X_NVisu,debugVisu(:,:,:,:,:),TRIM(INTSTAMP('Debugmesh',myRank))//'.dat')

SWRITE(UNIT_stdOut,'(A)')' WRITE DEBUGMESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WriteDebugMesh

END MODULE MOD_DebugMesh
