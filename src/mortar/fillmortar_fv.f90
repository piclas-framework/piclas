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

MODULE MOD_FillMortar_FV
#if USE_FV
!===================================================================================================================================
! Add comments please!
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
PUBLIC:: U_mortar_FV, Flux_Mortar_FV
PUBLIC:: Dx_BigToSmall_Mortar, Dx_SmallToBig_Mortar
PUBLIC:: CalcDiff_Mortar
!===================================================================================================================================

CONTAINS

SUBROUTINE U_Mortar_FV(U_in_master,U_in_slave,doReco,doMPISides)
!===================================================================================================================================
!> fills small non-conforming sides with data for master side with data from the corresponding large side
!> if (doReco=T) reconstruction already done from big element center to face, continue it to small face centers
!> fluxes will be calculated on the (virtual) small master sides and then summed up on the large side in Flux_Mortar_FV
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars       ,ONLY: MortarType,MortarInfo, nSides, SideToElem
USE MOD_Mesh_Vars       ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars       ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars_FV    ,ONLY: Face_xGP_FV
USE MOD_Gradient_Vars   ,ONLY: Gradient_elem
#ifdef discrete_velocity
USE MOD_Equation_Vars_FV,ONLY: DVMSpecData,DVMnSpecies,DVMColl,DVMDim
USE MOD_TimeDisc_Vars   ,ONLY: dt
#endif /*discrete_velocity*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#ifdef drift_diffusion
REAL,INTENT(INOUT) :: U_in_master(1:PP_nVar_FV+3,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: U_in_slave( 1:PP_nVar_FV+3,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
#else
REAL,INTENT(INOUT) :: U_in_master(1:PP_nVar_FV,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: U_in_slave( 1:PP_nVar_FV,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
#endif /*drift_diffusion*/
LOGICAL,INTENT(IN) :: doReco                                     !< flag whether reconstruction is used
LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,locSide,flip,ElemID
REAL         :: FaceToFace(3)
#ifdef drift_diffusion
REAL         :: U_tmp( PP_nVar_FV+3)
#else
REAL         :: U_tmp( PP_nVar_FV)
#endif /*drift_diffusion*/
#ifdef discrete_velocity
INTEGER      :: iSpec,iVel,jVel,kVel,upos,vFirstID
#endif /*discrete_velocity*/
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    IF (doReco) THEN
      ElemID    = SideToElem(S2E_ELEM_ID,MortarSideID)
      ! Get distance from big face center to small face center
      FaceToFace(1:3) = Face_xGP_FV(1:3,SideID) - Face_xGP_FV(1:3,MortarSideID)
#ifdef discrete_velocity
      IF (DVMColl) THEN
        !DVM specific reconstruction
        vFirstID = 0
        DO iSpec=1,DVMnSpecies
          ASSOCIATE(Sp => DVMSpecData(iSpec))
          DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
            upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
            U_tmp(upos) = U_in_master(upos,MortarSideID) &
                                          + Gradient_elem(1,upos,ElemID)*(FaceToFace(1)-Sp%Velos(iVel,1)*dt/2.) &
                                          + Gradient_elem(2,upos,ElemID)*(FaceToFace(2)-Sp%Velos(jVel,2)*dt/2.) &
                                          + Gradient_elem(3,upos,ElemID)*(FaceToFace(3)-Sp%Velos(kVel,3)*dt/2.)
            IF (DVMDim.LT.3) THEN
              U_tmp(Sp%nVarReduced+upos) = U_in_master(Sp%nVarReduced+upos,MortarSideID) &
                                          + Gradient_elem(1,Sp%nVarReduced+upos,ElemID)*(FaceToFace(1)-Sp%Velos(iVel,1)*dt/2.) &
                                          + Gradient_elem(2,Sp%nVarReduced+upos,ElemID)*(FaceToFace(2)-Sp%Velos(jVel,2)*dt/2.) &
                                          + Gradient_elem(3,Sp%nVarReduced+upos,ElemID)*(FaceToFace(3)-Sp%Velos(kVel,3)*dt/2.)
            END IF
            IF (Sp%InterID.EQ.2.OR.Sp%InterID.EQ.20) THEN
            ! rotational energy reduced distribution
              U_tmp(Sp%nVarERotStart+upos) = U_in_master(Sp%nVarERotStart+upos,MortarSideID) &
                                          + Gradient_elem(1,Sp%nVarERotStart+upos,ElemID)*(FaceToFace(1)-Sp%Velos(iVel,1)*dt/2.) &
                                          + Gradient_elem(2,Sp%nVarERotStart+upos,ElemID)*(FaceToFace(2)-Sp%Velos(jVel,2)*dt/2.) &
                                          + Gradient_elem(3,Sp%nVarERotStart+upos,ElemID)*(FaceToFace(3)-Sp%Velos(kVel,3)*dt/2.)
            END IF
          END DO; END DO; END DO
          vFirstID = vFirstID + Sp%nVar
          END ASSOCIATE
        END DO
      ELSE
#endif /*discrete_velocity*/
        U_tmp(1:PP_nVar_FV) = U_in_master(1:PP_nVar_FV,MortarSideID) + Gradient_elem(1,1:PP_nVar_FV,ElemID)*FaceToFace(1) &
                                                                         + Gradient_elem(2,1:PP_nVar_FV,ElemID)*FaceToFace(2) &
                                                                         + Gradient_elem(3,1:PP_nVar_FV,ElemID)*FaceToFace(3)
#ifdef discrete_velocity
      END IF !DVMColl
#endif
#ifdef drift_diffusion
      ! Electric field is only copied, reconstruction only for e- density
      U_tmp(PP_nVar_FV+1:PP_nVar_FV) = U_in_master(PP_nVar_FV+1:PP_nVar_FV,MortarSideID)
#endif
    ELSE
      U_tmp(:) = U_in_master(:,MortarSideID) ! no reconstruction, just copy
    END IF
    SELECT CASE(flip)
      CASE(0) ! master side
        U_in_master(:,SideID)=U_tmp(:)
      CASE(1:4) ! slave side (only for MPI)
        U_in_slave(:,SideID) =U_tmp(:)
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar_FV


SUBROUTINE Flux_Mortar_FV(Flux_Master,Flux_Slave,doMPISides)
!===================================================================================================================================
! fills master side from big non-conforming sides, by simply summing up the contributions from the small sides
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef drift_diffusion
REAL,INTENT(INOUT) :: Flux_Master(1:PP_nVar_FV+3,1:nSides)
REAL,INTENT(INOUT) :: Flux_Slave(1:PP_nVar_FV+3,1:nSides)
#else
REAL,INTENT(INOUT) :: Flux_Master(1:PP_nVar_FV,1:nSides)
REAL,INTENT(INOUT) :: Flux_Slave(1:PP_nVar_FV,1:nSides)
#endif /*drift_diffusion*/
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,iSide,flip
#ifdef drift_diffusion
REAL         :: Flux_tmp( PP_nVar_FV+3,1:4)
#else
REAL         :: Flux_tmp( PP_nVar_FV,1:4)
#endif /*drift_diffusion*/
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  Flux_tmp=0. !nullify as the sum is always done assuming 4 small sides
  ! Loop Mortar faces
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      Flux_tmp(:,iMortar)=Flux_Master(:,SideID)
    CASE(1:4) ! slave sides (should only occur for MPI)
      Flux_tmp(:,iMortar)=-Flux_Slave(:,SideID)
    END SELECT
  END DO
  Flux_Master(:,MortarSideID)=SUM(Flux_tmp(:,:),2) ! sum up contributions from 2 or 4 small sides
END DO !MortarSideID

END SUBROUTINE Flux_Mortar_FV

SUBROUTINE Dx_BigToSmall_Mortar(Dx_master,Dx_slave,doMPISides)
!===================================================================================================================================
! fills master side from small non-conforming sides, by calculating the vector from big element center to small element face
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars_FV,ONLY: Face_xGP_FV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: Dx_Master(1:3,1:nSides)
REAL,INTENT(INOUT) :: Dx_Slave(1:3,1:nSides)
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iMortar,nMortars,locSide
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,flip
REAL     :: FaceToFace(3)
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    ! Get distance from big face center to small face center
    FaceToFace(1:3) = Face_xGP_FV(1:3,SideID) - Face_xGP_FV(1:3,MortarSideID)
    SELECT CASE(flip)
      CASE(0) ! master side
        Dx_master(:,SideID)=Dx_Master(1:3,MortarSideID)+FaceToFace(1:3)
      CASE(1:4) ! slave side (only for MPI)
        Dx_slave(:,SideID) =Dx_Master(1:3,MortarSideID)+FaceToFace(1:3)
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE Dx_BigToSmall_Mortar

SUBROUTINE Dx_SmallToBig_Mortar(Dx_Master,Dx_Slave,doMPISides)
!===================================================================================================================================
! fills master side from big non-conforming sides, by simply averaging the contributions from the small sides
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: Dx_Master(1:3,1:nSides)
REAL,INTENT(INOUT) :: Dx_Slave(1:3,1:nSides)
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,iSide,flip
REAL     :: Dx_tmp(3,1:4)
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  Dx_tmp=0. !nullify as the sum is always done assuming 4 small sides
  ! Loop Mortar faces
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      Dx_tmp(:,iMortar)=Dx_Master(:,SideID)
    CASE(1:4) ! slave sides (should only occur for MPI)
      Dx_tmp(:,iMortar)=Dx_Slave(:,SideID)
    END SELECT
  END DO
  Dx_Master(:,MortarSideID)=SUM(Dx_tmp(:,:),2)/FLOAT(nMortars) ! sum up contributions from 2 or 4 small sides
END DO !MortarSideID

END SUBROUTINE Dx_SmallToBig_Mortar

SUBROUTINE CalcDiff_Mortar(doMPISides)
!===================================================================================================================================
! fills master side from big non-conforming sides, by simply averaging the contributions from the small sides
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Gradient_Vars,ONLY: Diff_side
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,iSide,flip
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  Diff_side(:,MortarSideID)=0. !nullify as the sum is always done assuming 4 small sides
  ! Loop Mortar faces
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      Diff_side(:,MortarSideID)=Diff_side(:,MortarSideID) + Diff_side(:,SideID)/FLOAT(nMortars)
    CASE(1:4) ! slave sides (should only occur for MPI)
      Diff_side(:,MortarSideID)=Diff_side(:,MortarSideID) - Diff_side(:,SideID)/FLOAT(nMortars)
    END SELECT
  END DO
END DO !MortarSideID

END SUBROUTINE CalcDiff_Mortar

#endif /*USE_FV*/
END MODULE MOD_FillMortar_FV