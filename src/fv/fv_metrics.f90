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
SUBROUTINE InitFV_Metrics(FV_dx_master,FV_dx_slave,doMPISides)
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
REAL, INTENT(OUT)                      :: FV_dx_master(3,0:PP_N,0:PP_N,1:nSides)
REAL, INTENT(OUT)                      :: FV_dx_slave(3,0:PP_N,0:PP_N,1:nSides)
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
  ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
  IF (ElemID.LT.0) CALL abort(__STAMP__,'FV metrics error')
  DO iCoord=1,3
    IF (Face_xGP(iCoord,0,0,SideID).GE.FV_PerBoxMax(iCoord))  THEN
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)+FV_PerBoxMin(iCoord)-FV_PerBoxMax(iCoord)
    ELSE IF (Face_xGP(iCoord,0,0,SideID).LE.FV_PerBoxMin(iCoord)) THEN
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)+FV_PerBoxMax(iCoord)-FV_PerBoxMin(iCoord)
    ELSE
      Face_temp(iCoord) = Face_xGP(iCoord,0,0,SideID)
    END IF
  END DO
  FV_dx_slave(:,0,0,SideID)=Face_temp(:)-Elem_xGP(:,0,0,0,ElemID)
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
  IF (ElemID.LT.0) CYCLE !small mortar sides don't have info of the big master element + skip MPI_YOUR sides
  FV_dx_master(:,0,0,SideID)=Face_xGP(:,0,0,SideID)-Elem_xGP(:,0,0,0,ElemID)
  ! mirror distance for ghost cells
  IF (SideID.LE.lastBCSide) FV_dx_slave(:,0,0,SideID) = -FV_dx_master(:,0,0,SideID)
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
USE MOD_Mesh_Vars          ,ONLY: nElems, nSides, ElemToSide, lastBCSide
USE MOD_FV_Vars            ,ONLY: FV_SysSol_slave, FV_SysSol_master, FV_SysSol_BC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL, INTENT(IN)                      :: dmaster(3,1:nSides)
REAL, INTENT(IN)                      :: dslave(3,1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: SideID, ElemID, locSideID, flip, info_dgesv, IPIV(PP_dim), singleDir
REAL                                   :: gradWeight, dElem(3), dMatrix(PP_dim,PP_dim), dMatrixBC(PP_dim,PP_dim)
LOGICAL                                :: BCelem
!==================================================================================================================================

DO ElemID = 1, nElems

  BCelem = .FALSE.
  dMatrix = 0.
  dMatrixBC = 0.
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    dElem(:) = dslave(:,SideID) - dmaster(:,SideID) !doesnt matter if elem is slave or master, dElem is always "squared"
    gradWeight = 1/VECNORM(dElem)**2
    dMatrix(1:PP_dim,1) = dMatrix(1:PP_dim,1) + gradWeight*dElem(1)*dElem(1:PP_dim)
    dMatrix(1:PP_dim,2) = dMatrix(1:PP_dim,2) + gradWeight*dElem(2)*dElem(1:PP_dim)
#if PP_dim == 3
    dMatrix(1:PP_dim,3) = dMatrix(1:PP_dim,3) + gradWeight*dElem(3)*dElem(1:PP_dim)
#endif
    IF (SideID.LE.lastBCSide) THEN
      BCelem = .TRUE.
      dMatrixBC(1:PP_dim,1) = dMatrixBC(1:PP_dim,1) + gradWeight*dElem(1)*dElem(1:PP_dim)
      dMatrixBC(1:PP_dim,2) = dMatrixBC(1:PP_dim,2) + gradWeight*dElem(2)*dElem(1:PP_dim)
#if PP_dim == 3
      dMatrixBC(1:PP_dim,3) = dMatrixBC(1:PP_dim,3) + gradWeight*dElem(3)*dElem(1:PP_dim)
#endif
    END IF
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
      CALL DGESV(PP_dim,1,dMatrix,PP_dim,IPIV,FV_SysSol_master(1:PP_dim,SideID),PP_dim,info_dgesv)
      IF(info_dgesv.NE.0) CALL abort(__STAMP__,'FV metrics error')
    ELSE
      FV_SysSol_slave(:,SideID) = dslave(:,SideID) - dmaster(:,SideID)
      CALL DGESV(PP_dim,1,dMatrix,PP_dim,IPIV,FV_SysSol_slave(1:PP_dim,SideID),PP_dim,info_dgesv)
      IF(info_dgesv.NE.0) CALL abort(__STAMP__,'FV metrics error')
    END IF
  END DO


  IF (BCelem) THEN
    singleDir = 0
    dMatrixBC = dMatrix - dMatrixBC
#if PP_dim == 3
    IF (dMatrixBC(1,1).EQ.0.) THEN
      IF (dMatrixBC(2,2).EQ.0.) singleDir = 3
      IF (dMatrixBC(3,3).EQ.0.) singleDir = 2
    ELSE IF (dMatrixBC(2,2).EQ.0..AND.dMatrixBC(3,3).EQ.0.) THEN
      singleDir = 1
    END IF
    DO locSideID=1,6
#else
    IF (dMatrixBC(1,1).EQ.0.) THEN
      singleDir = 2
    ELSE IF (dMatrixBC(2,2).EQ.0.) THEN
      singleDir = 1
    END IF
    DO locSideID=2,5
#endif
      SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
      IF (SideID.LE.lastBCSide) CYCLE
      !FV_SysSol_BC calculated with slave -> master direction
      IF (singleDir.GT.0) THEN
        FV_SysSol_BC(:,SideID) = 0.
        IF (dMatrixBC(singleDir,singleDir).GT.0.) THEN
          FV_SysSol_BC(singleDir,SideID) = (dslave(singleDir,SideID) - dmaster(singleDir,SideID))/dMatrixBC(singleDir,singleDir)
        END IF
      ELSE
        FV_SysSol_BC(:,SideID) = dslave(:,SideID) - dmaster(:,SideID)
        CALL DGESV(PP_dim,1,dMatrixBC,PP_dim,IPIV,FV_SysSol_BC(1:PP_dim,SideID),PP_dim,info_dgesv)
        IF(info_dgesv.NE.0) CALL abort(__STAMP__,'FV metrics BC error')
      END IF
    END DO
  END IF

END DO
END SUBROUTINE Build_FVSideMatrix

END MODULE MOD_FV_Metrics