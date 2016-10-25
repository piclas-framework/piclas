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
#ifndef PP_HDG
INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!===================================================================================================================================

CONTAINS

SUBROUTINE FillFlux(t,tDeriv,Flux,U_Minus,U_Plus,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Riemann,         ONLY:Riemann
USE MOD_Mesh_Vars,       ONLY:SideID_plus_lower,SideID_plus_upper
USE MOD_Mesh_Vars,       ONLY:SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
REAL,INTENT(IN)    :: t           ! time
INTEGER,INTENT(IN) :: tDeriv      ! deriv
REAL,INTENT(IN)    :: U_Minus(PP_nVar,0:PP_N, 0:PP_N,SideID_Minus_lower:SideID_Minus_Upper)
REAL,INTENT(IN)    :: U_Plus( PP_nVar,0:PP_N, 0:PP_N,SideID_Plus_Lower :SideID_Plus_Upper )
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID,firstSideID2
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  firstSideID = nBCSides+nInnerSides+1
  firstSideID2= firstSideID
  lastSideID  = firstSideID-1+nMPISides_MINE 
ELSE
  ! fill only InnerSides
  firstSideID = nBCSides+1
  firstSideID2= 1
  lastSideID  = firstSideID-1+nInnerSides 
END IF
!firstSideID=nBCSides+1
!lastSideID  =nBCSides+nInnerSides+nMPISides_MINE

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  CALL Riemann(Flux(:,:,:,SideID),U_Minus( :,:,:,SideID),U_Plus(  :,:,:,SideID),NormVec(:,:,:,SideID))
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux           (1:PP_nVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_Minus        (1:PP_nVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! Apply surface element size
DO SideID=firstSideID2,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
END DO

END SUBROUTINE FillFlux
#endif

END MODULE MOD_FillFlux
