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

MODULE MOD_SurfInt
!===================================================================================================================================
! Contains the different Surface integral formulations
! Computes the Surface integral for all faces using U and updates Ut
! Computes only inner surface integrals!
! Surface integrals are separated for each direction
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
INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

PUBLIC::SurfInt
!===================================================================================================================================
CONTAINS


SUBROUTINE SurfInt(Flux_Master,Flux_Slave,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
!USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
REAL,INTENT(IN)    :: Flux_Master(1:PP_nVar_FV,0:0,0:0,nSides)
REAL,INTENT(IN)    :: Flux_Slave(1:PP_nVar_FV,0:0,0:0,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,Flip,SideID,locSideID
INTEGER            :: firstSideID,lastSideID
!INTEGER            :: offsetElemCNProc
!===================================================================================================================================
IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  ! ignore MPI-faces and boundary faces
  IF(ElemID.LT.0) CYCLE ! boundary side is BC or MPI side
  Ut(:,0,0,0,ElemID) = Ut(:,0,0,0,ElemID) - Flux_Slave(:,0,0,SideID)
END DO ! SideID=1,nSides


DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip      = 0
  IF(ElemID.LT.0) CYCLE ! if master is MPI side
  Ut(:,0,0,0,ElemID) = Ut(:,0,0,0,ElemID) + Flux_Master(:,0,0,SideID)
END DO ! SideID=1,nSides

END SUBROUTINE SurfInt

END MODULE MOD_SurfInt
