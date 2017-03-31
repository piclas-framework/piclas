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
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides
USE MOD_Riemann,         ONLY:Riemann,RiemannPML
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
USE MOD_Mesh_Vars,       ONLY:firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_PML_vars,        ONLY:isPMLFace,FaceToPML,DoPML,isPMLFace,isPMLInterFace,PMLnVar
#ifdef maxwell
USE MOD_Equation,        ONLY:ExactFlux
#endif /*maxwell*/
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
REAL,INTENT(OUT)   :: Flux(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID
#ifdef maxwell
REAL               :: U_master_loc(PP_nVar,0:PP_N,0:PP_N)
REAL               :: U_slave_loc (PP_nVar,0:PP_N,0:PP_N)
#endif /*maxwell*/
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

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  IF(DoPML) THEN
    IF ( isPMLFace(SideID) )THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
                                 !     (flux-splitting!)
      IF(isPMLInterFace(SideID))THEN
#ifdef maxwell
        U_master_loc=U_Master(:,:,:,SideID)
        U_slave_loc =U_Slave(:,:,:,SideID)
        CALL ExactFlux(t,tDeriv,U_Master_loc,U_Slave_loc,NormVec(:,:,:,SideID),Face_xGP(:,:,:,SideID))
        CALL RiemannPML(Flux(1:32,:,:,SideID),U_Master_loc,U_Slave_loc, NormVec(:,:,:,SideID))
#endif /*maxwell*/
      ELSE
        CALL RiemannPML(Flux(1:32,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID), NormVec(:,:,:,SideID))
      END IF
    ELSE ! 2.) no PML, standard flux
      CALL Riemann(Flux(1:8,:,:,SideID), U_Master(:,:,:,SideID),  U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
    END IF
  ELSE ! no PML, standard flux
    CALL Riemann(Flux(:,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
  END IF ! DoPML
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux           (1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_master       (1:PP_nVar        ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) )
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
