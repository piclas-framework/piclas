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

SUBROUTINE FillFlux(t,tDeriv,Flux_Master,Flux_Slave,U_master,U_slave,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_GLobals
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides
USE MOD_Riemann,         ONLY:Riemann,RiemannPML
USE MOD_Riemann,         ONLY:RiemannDielectric,RiemannDielectricInterFace,RiemannDielectricInterFace2
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
USE MOD_Mesh_Vars,       ONLY:firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_PML_vars,        ONLY:PMLnVar
USE MOD_Dielectric_vars, ONLY:Dielectric_Master
USE MOD_Equation_Vars,   ONLY:DoExactFlux,isExactFluxInterFace
#ifdef maxwell
USE MOD_Riemann,         ONLY:ExactFlux
#endif /*maxwell*/
USE MOD_Interfaces_Vars, ONLY:InterfaceRiemann
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
REAL,INTENT(OUT)   :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(OUT)   :: Flux_Slave(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides (where the local proc is master) 
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
  lastSideID = lastInnerSide
END IF

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  SELECT CASE(InterfaceRiemann(SideID))
  ! Check every face and set the correct identifier for selecting the corresponding Riemann solver
  ! possible connections are (Master <-> Slave direction is important):
  !   - vaccuum    <-> vacuum       : RIEMANN_VACUUM            = 0
  !   - PML        <-> vacuum       : RIEMANN_PML               = 1
  !   - PML        <-> PML          : RIEMANN_PML               = 1
  !   - dielectric <-> dielectric   : RIEMANN_DIELECTRIC        = 2
  !   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC    = 3 ! for conservative fluxes (one flux)
  !   - vacuum      -> dielectri    : RIEMANN_VAC2DIELECTRIC    = 4 ! for conservative fluxes (one flux)
  !   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC_NC = 5 ! for non-conservative fluxes (two fluxes)
  !   - vacuum      -> dielectri    : RIEMANN_VAC2DIELECTRIC_NC = 6 ! for non-conservative fluxes (two fluxes)
  CASE(RIEMANN_VACUUM) 
    ! standard flux
    CALL Riemann(Flux_Master(1:8,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
  CASE(RIEMANN_PML) 
    ! RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations (flux-splitting!)
    CALL RiemannPML(Flux_Master(1:32,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
  CASE(RIEMANN_DIELECTRIC) 
    ! dielectric region <-> dielectric region
    CALL RiemannDielectric(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                           NormVec(:,:,:,SideID),Dielectric_Master(0:PP_N,0:PP_N,SideID))
  CASE(RIEMANN_DIELECTRIC2VAC)
    ! master is DIELECTRIC and slave PHYSICAL: A+(Eps0,Mu0) and A-(EpsR,MuR)
    CALL RiemannDielectricInterFace2(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                                     NormVec(:,:,:,SideID),Dielectric_Master(0:PP_N,0:PP_N,SideID))
  CASE(RIEMANN_VAC2DIELECTRIC) 
    ! master is PHYSICAL and slave DIELECTRIC: A+(EpsR,MuR) and A-(Eps0,Mu0)
    CALL RiemannDielectricInterFace(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                                    NormVec(:,:,:,SideID),Dielectric_Master(0:PP_N,0:PP_N,SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC)  ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! 1.) dielectric master side
    CALL RiemannDielectric(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                           NormVec(:,:,:,SideID),Dielectric_Master(0:PP_N,0:PP_N,SideID))
    ! 2.) vacuum slave side
    CALL Riemann(Flux_Slave(1:8,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
  CASE(RIEMANN_VAC2DIELECTRIC_NC) ! use non-conserving fluxes (two different fluxes for master and slave side) 
    ! 1.) dielectric slave side
    CALL RiemannDielectric(Flux_Slave(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                           NormVec(:,:,:,SideID),Dielectric_Master(0:PP_N,0:PP_N,SideID))
    ! 2.) vacuum master side
    CALL Riemann(Flux_Master(1:8,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
  CASE DEFAULT
    CALL abort(&
        __STAMP__&
        ,'Unknown interface type for Riemann solver (vacuum, dielectric, PML ...)')
  END SELECT
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux_Master    (1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_master       (1:PP_nVar        ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! Apply surface element size
DO SideID=firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    Flux_Master(:,p,q,SideID)=Flux_Master(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC,RIEMANN_VAC2DIELECTRIC_NC)
    ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! slaves sides have already been calculated
    DO q=0,PP_N; DO p=0,PP_N
      Flux_Slave(:,p,q,SideID)=Flux_Slave(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO; END DO
  CASE DEFAULT
    ! copy flux from Master side to slave side, DO not change sign
    Flux_slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)
  END SELECT
END DO


#ifdef maxwell
IF(DoExactFlux) THEN
  DO SideID=firstSideID,lastSideID
    IF (isExactFluxInterFace(SideID))THEN! CAUTION: Multiplication with SurfElem is done in ExactFlux
      CALL ExactFlux(t,tDeriv                                        &
                    , Flux_Master(1:PP_nVar+PMLnVar,:,:,SideID)      &
                    , Flux_Slave(1:PP_nVar+PMLnVar,:,:,SideID)       &
                    , U_Master(:,:,:,SideID)                         &
                    , U_Slave(:,:,:,SideID)                          &
                    , NormVec(:,:,:,SideID)                          &
                    , Face_xGP(1:3,:,:,SideID)                       &
                    , SurfElem(:,:,SideID)                           &
                    , SideID)
    END IF ! isExactFluxFace(SideID)
  END DO ! SideID
END IF                                           
#endif /*maxwell*/


END SUBROUTINE FillFlux
#endif

END MODULE MOD_FillFlux
