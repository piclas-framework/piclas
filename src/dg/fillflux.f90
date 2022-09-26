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
#if !(USE_HDG)
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
USE MOD_Mesh_Vars       ,ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars       ,ONLY: nSides,nBCSides
USE MOD_Riemann         ,ONLY: Riemann
USE MOD_Mesh_Vars       ,ONLY: NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux ,ONLY: GetBoundaryFlux
USE MOD_Mesh_Vars       ,ONLY: firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_PML_vars        ,ONLY: PMLnVar
USE MOD_Equation_Vars   ,ONLY: DoExactFlux,isExactFluxInterFace
#ifdef maxwell
USE MOD_Riemann         ,ONLY: ExactFlux
USE MOD_Interfaces_Vars ,ONLY: InterfaceRiemann
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
REAL,INTENT(OUT)   :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(OUT)   :: Flux_Slave(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
! Set the side range according to MPI or no MPI
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

! =============================
! Workflow:
!
!  1.  compute flux for non-BC sides
!  1.1) advective flux
!  ( 1.2) viscous flux )
!  ( 1.3) add up viscous flux to Flux_master )
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem
!  4.  copy flux from Flux_master to Flux_slave
!  ( 5.  convert FV flux to DG flux at mixed interfaces )
!  6. Exact flux determination (inner BC)
!==============================

! 1. compute flux for non-BC sides: Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
#ifdef maxwell
  CALL Riemann(Flux_Master(:,:,:,SideID),Flux_Slave(:,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID),SideID)
#else
  CALL Riemann(Flux_Master(:,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
#endif /* maxwell */
END DO ! SideID

! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux_Master    (1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_master       (1:PP_nVar        ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! 3. multiply by SurfElem: Apply surface element size
DO SideID=firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    Flux_Master(:,p,q,SideID)=Flux_Master(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
#ifdef maxwell
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC,RIEMANN_VAC2DIELECTRIC_NC)
    ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! slaves sides have already been calculated
    DO q=0,PP_N; DO p=0,PP_N
      Flux_Slave(:,p,q,SideID)=Flux_Slave(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO; END DO
  CASE DEFAULT
#endif /*maxwell*/
    ! 4. copy flux from master side to slave side: DO not change sign
    Flux_slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)
#ifdef maxwell
  END SELECT
#endif /*maxwell*/
END DO

#ifdef maxwell
!  6. Exact flux determination (inner BC)
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
