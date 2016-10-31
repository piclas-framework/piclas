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

SUBROUTINE FillFlux(t,tDeriv,Flux,U_master,U_slave,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Riemann,         ONLY:Riemann
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
USE MOD_Mesh_Vars,       ONLY:firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
REAL,INTENT(IN)    :: t           ! time
INTEGER,INTENT(IN) :: tDeriv      ! deriv
REAL,INTENT(IN)    :: U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: U_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  ! fill only flux for MINE MPISides (where the local proc is master) 
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
  lastSideID = lastInnerSide
END IF
!firstSideID=nBCSides+1
!lastSideID  =nBCSides+nInnerSides+nMPISides_MINE

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  CALL Riemann(Flux(:,:,:,SideID),U_master( :,:,:,SideID),U_slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux           (1:PP_nVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_master       (1:PP_nVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! Apply surface element size
DO SideID=firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
END DO

END SUBROUTINE FillFlux
#endif

END MODULE MOD_FillFlux
