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

MODULE MOD_Prolong_FV
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
PUBLIC::ProlongToFace_FV,ProlongToFace_ElemCopy
PUBLIC::ProlongToOutput
!===================================================================================================================================

CONTAINS

SUBROUTINE ProlongToFace_FV(Uvol,Uface_master,Uface_slave,doMPISides)
!===================================================================================================================================
! Performs finite volumes reconstruction by applying limited gradients from cell centers to faces
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Gradient_Vars,      ONLY: Grad_dx_master, Grad_dx_slave, Gradient_elem
USE MOD_Mesh_Vars,          ONLY: nSides, SideToElem
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide, lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,firstMortarMPISide,lastMortarMPISide
#ifdef discrete_velocity
USE MOD_TimeDisc_Vars,      ONLY: dt
USE MOD_Equation_Vars_FV,   ONLY: DVMSpecData, DVMnSpecies, DVMDim
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
#ifdef drift_diffusion
REAL,INTENT(IN)                 :: Uvol(PP_nVar_FV+3,0:0,0:0,0:0,1:PP_nElems)
#else
REAL,INTENT(IN)                 :: Uvol(PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#ifdef drift_diffusion
REAL,INTENT(INOUT)              :: Uface_master(PP_nVar_FV+3,0:0,0:0,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave(PP_nVar_FV+3,0:0,0:0,1:nSides)
#else
REAL,INTENT(INOUT)              :: Uface_master(PP_nVar_FV,0:0,0:0,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave(PP_nVar_FV,0:0,0:0,1:nSides)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID,SideID,firstSideID,lastSideID,iVel,jVel,kVel,upos
#ifdef discrete_velocity
INTEGER                         :: iSpec, vFirstID
#endif
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = firstMPISide_YOUR
  lastSideID  = lastMPISide_YOUR
ELSE
  ! (Mortar-)InnerSides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side !ElemID=-1 if not existing
  ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
  IF (ElemID.LT.0) CYCLE !mpi-mortar whatever
#ifdef discrete_velocity
  !DVM specific reconstruction
  vFirstID = 0
  DO iSpec=1,DVMnSpecies
    ASSOCIATE(Sp => DVMSpecData(iSpec))
    DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
      upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      Uface_slave(upos,0,0,SideID) = Uvol(upos,0,0,0,ElemID) &
                                    + Gradient_elem(1,upos,ElemID)*(Grad_dx_slave(1,SideID)-Sp%Velos(iVel,1)*dt/2.) &
                                    + Gradient_elem(2,upos,ElemID)*(Grad_dx_slave(2,SideID)-Sp%Velos(jVel,2)*dt/2.) &
                                    + Gradient_elem(3,upos,ElemID)*(Grad_dx_slave(3,SideID)-Sp%Velos(kVel,3)*dt/2.)
      IF (DVMDim.LT.3) THEN
        Uface_slave(Sp%nVar/2+upos,0,0,SideID) = Uvol(Sp%nVar/2+upos,0,0,0,ElemID) &
                                    + Gradient_elem(1,Sp%nVar/2+upos,ElemID)*(Grad_dx_slave(1,SideID)-Sp%Velos(iVel,1)*dt/2.) &
                                    + Gradient_elem(2,Sp%nVar/2+upos,ElemID)*(Grad_dx_slave(2,SideID)-Sp%Velos(jVel,2)*dt/2.) &
                                    + Gradient_elem(3,Sp%nVar/2+upos,ElemID)*(Grad_dx_slave(3,SideID)-Sp%Velos(kVel,3)*dt/2.)
      END IF
    END DO; END DO; END DO
    vFirstID = vFirstID + Sp%nVar
    END ASSOCIATE
  END DO
#else
  Uface_slave(1:PP_nVar_FV,0,0,SideID) = Uvol(1:PP_nVar_FV,0,0,0,ElemID) &
                            + Gradient_elem(1,1:PP_nVar_FV,ElemID)*Grad_dx_slave(1,SideID) &
                            + Gradient_elem(2,1:PP_nVar_FV,ElemID)*Grad_dx_slave(2,SideID) &
                            + Gradient_elem(3,1:PP_nVar_FV,ElemID)*Grad_dx_slave(3,SideID)
#endif
#ifdef drift_diffusion
  Uface_slave(PP_nVar_FV+1:PP_nVar_FV+3,0,0,SideID) = Uvol(PP_nVar_FV+1:PP_nVar_FV+3,0,0,0,ElemID)
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
  IF (ElemID.LT.0) CYCLE !for small mortar sides without info on big master element
#ifdef discrete_velocity
  !DVM specific reconstruction
  vFirstID = 0
  DO iSpec=1,DVMnSpecies
    ASSOCIATE(Sp => DVMSpecData(iSpec))
    DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
      upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      Uface_master(upos,0,0,SideID) = Uvol(upos,0,0,0,ElemID) &
                                    + Gradient_elem(1,upos,ElemID)*(Grad_dx_master(1,SideID)-Sp%Velos(iVel,1)*dt/2.) &
                                    + Gradient_elem(2,upos,ElemID)*(Grad_dx_master(2,SideID)-Sp%Velos(jVel,2)*dt/2.) &
                                    + Gradient_elem(3,upos,ElemID)*(Grad_dx_master(3,SideID)-Sp%Velos(kVel,3)*dt/2.)
      IF (DVMDim.LT.3) THEN
        Uface_master(Sp%nVar/2+upos,0,0,SideID) = Uvol(Sp%nVar/2+upos,0,0,0,ElemID) &
                                    + Gradient_elem(1,Sp%nVar/2+upos,ElemID)*(Grad_dx_master(1,SideID)-Sp%Velos(iVel,1)*dt/2.) &
                                    + Gradient_elem(2,Sp%nVar/2+upos,ElemID)*(Grad_dx_master(2,SideID)-Sp%Velos(jVel,2)*dt/2.) &
                                    + Gradient_elem(3,Sp%nVar/2+upos,ElemID)*(Grad_dx_master(3,SideID)-Sp%Velos(kVel,3)*dt/2.)
      END IF
    END DO; END DO; END DO
    vFirstID = vFirstID + Sp%nVar
    END ASSOCIATE
  END DO
#else
  Uface_master(1:PP_nVar_FV,0,0,SideID) = Uvol(1:PP_nVar_FV,0,0,0,ElemID) &
                            + Gradient_elem(1,1:PP_nVar_FV,ElemID)*Grad_dx_master(1,SideID) &
                            + Gradient_elem(2,1:PP_nVar_FV,ElemID)*Grad_dx_master(2,SideID) &
                            + Gradient_elem(3,1:PP_nVar_FV,ElemID)*Grad_dx_master(3,SideID)
#endif
#ifdef drift_diffusion
  Uface_master(PP_nVar_FV+1:PP_nVar_FV+3,0,0,SideID) = Uvol(PP_nVar_FV+1:PP_nVar_FV+3,0,0,0,ElemID)
#endif
END DO !SideID

END SUBROUTINE ProlongToFace_FV

SUBROUTINE ProlongToOutput(Uvol,Uout)
!===================================================================================================================================
! Performs finite volumes reconstruction by applying limited gradients from cell centers to output points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Gradient_Vars,      ONLY: Gradient_elem
USE MOD_Mesh_Vars_FV,       ONLY: Elem_xGP_PP_1,Elem_xGP_FV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Uvol(PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Uout(PP_nVar_FV,0:PP_1,0:PP_1,0:PP_1,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID,i,j,k
!===================================================================================================================================
DO ElemID=1,PP_nElems
  DO k=0,PP_1
    DO j=0,PP_1
      DO i=0,PP_1
        Uout(:,i,j,k,ElemID) = Uvol(:,0,0,0,ElemID) &
                            + Gradient_elem(1,:,ElemID)*(Elem_xGP_PP_1(1,i,j,k,ElemID)-Elem_xGP_FV(1,0,0,0,ElemID)) &
                            + Gradient_elem(2,:,ElemID)*(Elem_xGP_PP_1(2,i,j,k,ElemID)-Elem_xGP_FV(2,0,0,0,ElemID)) &
                            + Gradient_elem(3,:,ElemID)*(Elem_xGP_PP_1(3,i,j,k,ElemID)-Elem_xGP_FV(3,0,0,0,ElemID))
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ProlongToOutput

SUBROUTINE ProlongToFace_ElemCopy(VarDim,ElemVar,SideVar_master,SideVar_slave,doMPISides)
!===================================================================================================================================
! Simply copies element-based variable to sides
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,          ONLY: nSides, SideToElem, nElems
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide, lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
INTEGER,INTENT(IN)              :: VarDim
REAL,INTENT(IN)                 :: ElemVar(VarDim,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: SideVar_master(VarDim,1:nSides)
REAL,INTENT(INOUT)              :: SideVar_slave(VarDim,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID,SideID,firstSideID,lastSideID
!===================================================================================================================================
IF(doMPISides)THEN
! only YOUR MPI Sides are filled
firstSideID = firstMPISide_YOUR
lastSideID  = lastMPISide_YOUR
ELSE
! (Mortar-)InnerSides are filled
firstSideID = firstInnerSide
lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
! neighbor side !ElemID=-1 if not existing
ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
!Copy Bulk Velocity in U for Dense Fokker Planck instead of Finite Volumes Solution
IF (ElemID.LT.0) CYCLE !mpi-mortar whatever
SideVar_slave(:,SideID) = ElemVar(:,ElemID)
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
IF (ElemID.LT.0) CYCLE !for small mortar sides without info on big master element
SideVar_master(:,SideID)=ElemVar(:,ElemID)
END DO !SideID

END SUBROUTINE ProlongToFace_ElemCopy

END MODULE MOD_Prolong_FV

