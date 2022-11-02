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

MODULE MOD_FV_Metrics
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE InitFV_Metrics
  MODULE PROCEDURE InitFV_Metrics
END INTERFACE

PUBLIC::InitFV_Metrics
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute the distance FV_dx from interface to cell center of the master/slave element
!==================================================================================================================================
SUBROUTINE InitFV_Metrics(dx_slave_temp,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_Mesh_Vars,          ONLY: SideToElem, Face_xGP, Elem_xGP
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,nSides,firstMortarMPISide,lastMortarMPISide
#if (PP_TimeDiscMethod==600)
USE MOD_Mesh_Vars,          ONLY: NormVec
USE MOD_Equation_Vars,      ONLY: DVMnVelos, DVMVelos, DVMSpeciesData
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL, INTENT(OUT)                      :: dx_slave_temp(1:PP_nVar+1,1:nSides)
LOGICAL, INTENT(IN)                    :: doMPISides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: SideID, ElemID, firstSideID, lastSideID
#if (PP_TimeDiscMethod==600)
INTEGER                                :: iVel, jVel, kVel, upos
#endif
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  Build FV Metrics ...'

IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = firstMPISide_YOUR
  lastSideID  = lastMPISide_YOUR
ELSE
  ! (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side !ElemID=-1 if not existing
  ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
  IF (ElemID.LT.0) CYCLE
  dx_slave_temp(PP_nVar+1,SideID)=VECNORM(Elem_xGP(:,0,0,0,ElemID)-Face_xGP(:,0,0,SideID))

#if (PP_TimeDiscMethod==600)
  DO kVel=1, DVMnVelos(3); DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    dx_slave_temp(upos,SideID) = NormVec(1,0,0,SideID)*DVMVelos(iVel,1) &
                               + NormVec(2,0,0,SideID)*DVMVelos(jVel,2) &
                               + NormVec(3,0,0,SideID)*DVMVelos(kVel,3)
    IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
      dx_slave_temp(NINT(PP_nVar/2.)+upos,SideID)=dx_slave_temp(upos,SideID)
    END IF
  END DO; END DO; END DO
#endif
END DO

! Second process Minus/Master sides, U_Minus is always MINE
! master side, flip=0
IF(doMPISides)THEN
  ! only MPI mortars
  firstSideID = firstMortarMPISide
   lastSideID =  lastMortarMPISide
ELSE
  ! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstBCSide
   lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  FV_dx_master(SideID)=VECNORM(Elem_xGP(:,0,0,0,ElemID)-Face_xGP(:,0,0,SideID))

#if (PP_TimeDiscMethod==600)
  DO kVel=1, DVMnVelos(3); DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    DVMtraj_master(upos,SideID) = NormVec(1,0,0,SideID)*DVMVelos(iVel,1) &
                                + NormVec(2,0,0,SideID)*DVMVelos(jVel,2) &
                                + NormVec(3,0,0,SideID)*DVMVelos(kVel,3)
    IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
      DVMtraj_master(NINT(PP_nVar/2.)+upos,SideID)=DVMtraj_master(upos,SideID)
    END IF
  END DO; END DO; END DO
#endif
END DO !SideID

SWRITE(UNIT_stdOut,'(A)')' Done !'

END SUBROUTINE InitFV_Metrics


END MODULE MOD_FV_Metrics
