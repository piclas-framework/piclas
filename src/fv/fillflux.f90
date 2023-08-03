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

INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!===================================================================================================================================

CONTAINS

SUBROUTINE FillFlux(t,Flux_Master,Flux_Slave,U_master,U_slave,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_GLobals
USE MOD_PreProc
USE MOD_Mesh_Vars       ,ONLY: nSides,nBCSides
USE MOD_Riemann         ,ONLY: Riemann
USE MOD_Mesh_Vars       ,ONLY: NormVec_FV,TangVec1_FV, tangVec2_FV, SurfElem_FV,Face_xGP_FV,Elem_xGP_FV, SideToElem
USE MOD_GetBoundaryFlux ,ONLY: GetBoundaryFlux
USE MOD_Mesh_Vars       ,ONLY: firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
#ifdef drift_diffusion
USE MOD_Equation_Vars_FV,ONLY: EFluid_GradSide
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           ! time
REAL,INTENT(IN)    :: U_master(PP_nVar_FV,0:0,0:0,1:nSides)
REAL,INTENT(IN)    :: U_slave (PP_nVar_FV,0:0,0:0,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux_Master(1:PP_nVar_FV,0:0,0:0,nSides)
REAL,INTENT(OUT)   :: Flux_Slave(1:PP_nVar_FV,0:0,0:0,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,firstSideID_wo_BC,firstSideID ,lastSideID,ElemID
REAL               :: E(3)
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
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem
!  4.  copy flux from Flux_master to Flux_slave
!==============================

! 1. compute flux for non-BC sides: Compute fluxes on 0, no additional interpolation required
DO SideID=firstSideID_wo_BC,lastSideID
#ifdef drift_diffusion
  ElemID   = SideToElem(S2E_ELEM_ID,SideID)
  E=(/1.,0.,0./)
  ! print*, 'elemxgp', Elem_xGP_FV(:,:,:,:,ElemID)
  ! print*, 'face', Face_xGP_FV(:,:,:,SideID)
  ! print*, 'nv', NormVec_FV(:,:,:,SideID)
  ! print*, 'surfelem', SurfElem_FV(:,:,SideID)
  CALL Riemann(Flux_Master(:,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec_FV(:,:,:,SideID), &
               EFluid_GradSide(SideID),E)
#else
  CALL Riemann(Flux_Master(:,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec_FV(:,:,:,SideID))
#endif
END DO ! SideID

! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,0,Flux_Master    (1:PP_nVar_FV,0:0,0:0,1:nBCSides) &
                               ,U_master       (1:PP_nVar_FV     ,0:0,0:0,1:nBCSides) &
                               ,NormVec_FV        (1:3              ,0:0,0:0,1:nBCSides) &
                               ,TangVec1_FV       (1:3              ,0:0,0:0,1:nBCSides) &
                               ,TangVec2_FV       (1:3              ,0:0,0:0,1:nBCSides) &
                               ,Face_XGP_FV       (1:3              ,0:0,0:0,1:nBCSides) )
END IF

! 3. multiply by SurfElem: Apply surface element size
DO SideID=firstSideID,lastSideID
  Flux_Master(:,0,0,SideID)=Flux_Master(:,0,0,SideID)*SurfElem_FV(0,0,SideID)
  ! 4. copy flux from master side to slave side: DO not change sign
  Flux_slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)
END DO

END SUBROUTINE FillFlux

END MODULE MOD_FillFlux
