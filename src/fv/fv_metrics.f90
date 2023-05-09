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

PUBLIC::InitFV_Metrics, Build_FVSideMatrix
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute the distance FV_dx from interface to cell center of the master/slave element
!==================================================================================================================================
SUBROUTINE InitFV_Metrics(dx_master_temp,dx_slave_temp,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem, Face_xGP, Elem_xGP
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide, lastBCSide, firstMPISide_MINE, lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,nSides,firstMortarMPISide,lastMortarMPISide
USE MOD_FV_Vars,            ONLY: FV_PerBoxMax,FV_PerBoxMin
#if (PP_TimeDiscMethod==600)
USE MOD_Mesh_Vars,          ONLY: NormVec
USE MOD_Equation_Vars,      ONLY: DVMnVelos, DVMVelos, DVMSpeciesData
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL, INTENT(OUT)                      :: dx_master_temp(3,1:PP_nVar+1,1:nSides)
REAL, INTENT(OUT)                      :: dx_slave_temp(3,1:PP_nVar+1,1:nSides)
LOGICAL, INTENT(IN)                    :: doMPISides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: SideID, ElemID, firstSideID, lastSideID, iCoord
REAL                                   :: Face_temp(3)
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
  ! (Mortar-)InnerSides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side !ElemID=-1 if not existing
  ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
  IF (ElemID.LT.0) print*, 'fvmetrics buuuuuuug'
  ! IF (SideID.EQ.9) print*, myRank, 'SIDE9slave', Face_xGP(:,0,0,SideID)
  ! IF (SideID.EQ.9) print*, myRank, 'SIDE9slave', ElemID, Elem_xGP(:,0,0,0,ElemID)
  DO iCoord=1,3
    IF (Face_xGP(iCoord,0,0,SideID).GE.FV_PerBoxMax(iCoord))  THEN
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)+FV_PerBoxMin(iCoord)-FV_PerBoxMax(iCoord)
    ELSE IF (Face_xGP(iCoord,0,0,SideID).LE.FV_PerBoxMin(iCoord)) THEN
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)+FV_PerBoxMax(iCoord)-FV_PerBoxMin(iCoord)
    ELSE
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)
    END IF
  END DO
  dx_slave_temp(:,PP_nVar+1,SideID)=Face_temp(:)-Elem_xGP(:,0,0,0,ElemID)

#if (PP_TimeDiscMethod==600)
  DO kVel=1, DVMnVelos(3); DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    dx_slave_temp(1,upos,SideID) = DVMVelos(iVel,1)
    dx_slave_temp(2,upos,SideID) = DVMVelos(jVel,2)                                             ! remove normvec communication
    dx_slave_temp(3,upos,SideID) = DVMVelos(kVel,3)                                             ! change the velo storage
    IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
      dx_slave_temp(:,PP_nVar/2+upos,SideID)=dx_slave_temp(:,upos,SideID)
    END IF
  END DO; END DO; END DO
#endif
END DO

! Second process Minus/Master sides, U_Minus is always MINE
! master side, flip=0
IF(doMPISides)THEN
  ! only MPI mortars
  firstSideID = firstMPISide_MINE
  lastSideID =  lastMortarMPISide
ELSE
  ! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstBCSide
  lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  ! IF (SideID.EQ.9) print*, myRank, 'SIDE9master', Face_xGP(:,0,0,SideID), firstBCSide, lastBCSide, firstMPISide_MINE
  ! IF (SideID.EQ.9) print*, myRank, 'SIDE9master', ElemID, Elem_xGP(:,0,0,0,ElemID)
  IF (ElemID.LT.0) CYCLE !small mortar sides don't have info of the big master element + skip MPI_YOUR sides
  dx_master_temp(:,PP_nVar+1,SideID)=Face_xGP(:,0,0,SideID)-Elem_xGP(:,0,0,0,ElemID)

#if (PP_TimeDiscMethod==600)
  DO kVel=1, DVMnVelos(3); DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    dx_master_temp(1,upos,SideID) = DVMVelos(iVel,1)
    dx_master_temp(2,upos,SideID) = DVMVelos(jVel,2)
    dx_master_temp(3,upos,SideID) = DVMVelos(kVel,3)
    IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
      dx_master_temp(:,PP_nVar/2+upos,SideID)=dx_master_temp(:,upos,SideID)
    END IF
  END DO; END DO; END DO
#endif
  ! mirror distance for ghost cells
  IF (SideID.LE.lastBCSide) dx_slave_temp(:,:,SideID) = -dx_master_temp(:,:,SideID)
END DO !SideID

SWRITE(UNIT_stdOut,'(A)')' Done !'

END SUBROUTINE InitFV_Metrics

!==================================================================================================================================
!> Build element-neighbour-distance matrix, invert it and multiply by slave->master vector
!==================================================================================================================================
SUBROUTINE Build_FVSideMatrix(dmaster,dslave)
  ! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: nElems, nSides, ElemToSide
USE MOD_FV_Vars            ,ONLY: FV_SysSol_slave, FV_SysSol_master
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL, INTENT(IN)                      :: dmaster(3,1:nSides)
REAL, INTENT(IN)                      :: dslave(3,1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: SideID, ElemID, locSideID, flip, info_dgesv, IPIV(PP_dim)
REAL                                   :: gradWeight, dElem(3), dMatrix(PP_dim,PP_dim)
!==================================================================================================================================

DO ElemID = 1, nElems
  dMatrix = 0.
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    ! print*, myRank, ElemID, SideID  
    ! print*, myRank, 'slave', SideID, dslave(:,SideID)
    ! print*, myRank, 'master', SideID, dmaster(:,SideID)
    dElem(:) = dslave(:,SideID) - dmaster(:,SideID) !doesnt matter if elem is slave or master, dElem is always "squared"
    gradWeight = 1/VECNORM(dElem)**2
    ! print*, 'elem', ElemID
    ! print*, 'side', SideID
    ! print*, 'dslave', dslave(:,SideID)
    ! print*, 'dmaster', dmaster(:,SideID)
    ! print*, 'dELem', dElem
    ! print*, 'gradweight', gradWeight
    ! print*, 'dMatrix', dMatrix(1,:)
    ! print*, '       ', dMatrix(2,:)
    ! ! print*, '       ', dMatrix(3,:)
    ! read*
    dMatrix(1:PP_dim,1) = dMatrix(1:PP_dim,1) + gradWeight*dElem(1)*dElem(1:PP_dim)
    dMatrix(1:PP_dim,2) = dMatrix(1:PP_dim,2) + gradWeight*dElem(2)*dElem(1:PP_dim)
#if PP_dim == 3
    dMatrix(1:PP_dim,3) = dMatrix(1:PP_dim,3) + gradWeight*dElem(3)*dElem(1:PP_dim)
#endif
  END DO
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    flip=ElemToSide(E2S_FLIP,locSideID,ElemID)
    !FV_SysSol needs to be calculated with elem -> neighbour direction
    IF (flip.EQ.0) THEN
      FV_SysSol_master(:,SideID) = dmaster(:,SideID) - dslave(:,SideID)
      ! print*, 'master', SideID, FV_SysSol_master(:,SideID)
      CALL DGESV(PP_dim,1,dMatrix,PP_dim,IPIV,FV_SysSol_master(1:PP_dim,SideID),PP_dim,info_dgesv)
      IF(info_dgesv.NE.0) CALL abort(__STAMP__,'FV metrics error')
    ELSE
      FV_SysSol_slave(:,SideID) = dslave(:,SideID) - dmaster(:,SideID)
      ! print*, 'slave', SideID, FV_SysSol_slave(:,SideID)
      CALL DGESV(PP_dim,1,dMatrix,PP_dim,IPIV,FV_SysSol_slave(1:PP_dim,SideID),PP_dim,info_dgesv)
      IF(info_dgesv.NE.0) CALL abort(__STAMP__,'FV metrics error')
    END IF
  END DO
END DO

END SUBROUTINE Build_FVSideMatrix

END MODULE MOD_FV_Metrics