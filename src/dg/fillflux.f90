#include "boltzplatz.h"

MODULE MOD_FillFlux
!===================================================================================================================================
! Fills the inner, periodic and bc fluxes
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

INTERFACE FillFlux_BC
  MODULE PROCEDURE FillFlux_BC
END INTERFACE

PUBLIC::FillFlux,FillFlux_BC
!===================================================================================================================================

CONTAINS


SUBROUTINE FillFlux(Flux,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY: U_Minus,U_Plus
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Riemann,         ONLY: Riemann
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  firstSideID = nBCSides+nInnerSides+1
  lastSideID  = firstSideID-1+nMPISides_MINE 
ELSE
  ! fill only InnerSides
  firstSideID = nBCSides+1
  lastSideID  = firstSideID-1+nInnerSides 
END IF
!firstSideID=nBCSides+1
!lastSideID  =nBCSides+nInnerSides+nMPISides_MINE
DO SideID=firstSideID,lastSideID
  CALL Riemann(Flux(:,:,:,SideID),     U_Minus(:,:,:,SideID),     U_Plus(:,:,:,SideID), &
               NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID))
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO ! SideID

END SUBROUTINE FillFlux


SUBROUTINE FillFlux_BC(t,tDeriv,Flux)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY: U_Minus
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,SurfElem,BCFace_xGP,BC,BoundaryType
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
USE MOD_Equation_Vars,   ONLY: IniExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
INTEGER,INTENT(IN) :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            ::SideID,BCType,BCState,p,q
!===================================================================================================================================
! fill flux for boundary sides
DO SideID=1,nBCSides
  BCType=Boundarytype(BC(SideID),BC_TYPE)
  BCState=Boundarytype(BC(SideID),BC_STATE)
  IF (BCState.EQ.0)BCState=IniExactFunc
  CALL GetBoundaryFlux(Flux(:,:,:,SideID),BCType,BCState,BCFace_xGP(:,:,:,SideID),NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),  &
                       TangVec2(:,:,:,SideID),t,tDeriv,U_Minus(:,:,:,SideID))
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO! SideID
END SUBROUTINE FillFlux_BC

END MODULE MOD_FillFlux
