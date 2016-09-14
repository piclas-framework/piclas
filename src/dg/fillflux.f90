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
USE MOD_Riemann,         ONLY:Riemann,RiemannPML
USE MOD_Mesh_Vars,       ONLY:SideID_plus_lower,SideID_plus_upper
USE MOD_Mesh_Vars,       ONLY:SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,BCFace_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
USE MOD_PML_vars,        ONLY:isPMLFace,FaceToPML,DoPML,isPMLFace,isPMLInterFace,PMLnVar
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
REAL,INTENT(OUT)   :: Flux(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
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
  IF(DoPML) THEN
    IF ( isPMLFace(SideID) )THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
                                 !     (flux-splitting!)
      CALL RiemannPML(Flux(1:32,:,:,SideID),U_Minus(:,:,:,SideID),U_Plus(:,:,:,SideID), NormVec(:,:,:,SideID))
      DO q=0,PP_N; DO p=0,PP_N
          Flux(1:PP_nVar+PMLnVar,p,q,SideID)=Flux(1:32,p,q,SideID)*SurfElem(p,q,SideID)
      END DO; END DO
    ELSE ! 2.) no PML, standard flux
      CALL Riemann(Flux(1:8,:,:,SideID), U_Minus(:,:,:,SideID),  U_Plus(:,:,:,SideID),NormVec(:,:,:,SideID))
      DO q=0,PP_N; DO p=0,PP_N
          Flux(1:PP_nVar,p,q,SideID)=Flux(1:PP_nVar,p,q,SideID)*SurfElem(p,q,SideID)
      END DO; END DO
      Flux(1+PP_nVar:PP_nVar+PMLnVar,:,:,SideID)=0. ! Flux for auxiliary variables at the interface is set to zero
    END IF
  ELSE ! no PML, standard flux
    CALL Riemann(Flux(:,:,:,SideID),U_Minus( :,:,:,SideID),U_Plus(  :,:,:,SideID),NormVec(:,:,:,SideID))
  END IF ! DoPML
END DO ! SideID
  
IF(doMPISides.EQV..FALSE.)THEN
  CALL GetBoundaryFlux(t,tDeriv, Flux ,U_Minus ,NormVec, TangVec1, TangVec2, BCFace_xGP)
END IF

! Apply surface element size
DO SideID=firstSideID2,lastSideID
  IF(DoPML) THEN
    IF ( isPMLFace(SideID) )THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
                                 !     (flux-splitting!)
      DO q=0,PP_N; DO p=0,PP_N
          Flux(1:PP_nVar+PMLnVar,p,q,SideID)=Flux(1:PP_nVar+PMLnVar,p,q,SideID)*SurfElem(p,q,SideID)
      END DO; END DO
    ELSE ! 2.) no PML, standard flux
      DO q=0,PP_N; DO p=0,PP_N
          Flux(1:PP_nVar,p,q,SideID)=Flux(1:PP_nVar,p,q,SideID)*SurfElem(p,q,SideID)
      END DO; END DO
      Flux(1+PP_nVar:PP_nVar+PMLnVar,:,:,SideID)=0. ! Flux for auxiliary variables at the interface is set to zero
    END IF
  ELSE ! no PML, standard flux
    DO q=0,PP_N; DO p=0,PP_N
        Flux(1:PP_nVar,p,q,SideID)=Flux(1:PP_nVar,p,q,SideID)*SurfElem(p,q,SideID)
    END DO; END DO
  END IF ! DoPML
END DO ! SideID

END SUBROUTINE FillFlux
#endif

END MODULE MOD_FillFlux
