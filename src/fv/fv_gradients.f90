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

MODULE MOD_FV_Gradients
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


PUBLIC::CalcFVGradients, SumFVGradients
!==================================================================================================================================

CONTAINS

SUBROUTINE CalcFVGradients(doMPISides)
!===================================================================================================================================
! Compute gradients for fv reconstruction
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Vars                  ,ONLY: U_master, U_slave, FV_dx_master, FV_dx_slave, FV_gradU_side, FV_SysSol_BC
USE MOD_FV_Vars                  ,ONLY: LimiterTypeBC, FV_VktK
USE MOD_Mesh_Vars                ,ONLY: NormVec,Face_xGP
USE MOD_GetBoundaryGrad          ,ONLY: GetBoundaryGrad
USE MOD_Mesh_Vars                ,ONLY: firstBCSide,lastBCSide,firstInnerSide, lastInnerSide
USE MOD_Mesh_Vars                ,ONLY: firstMPISide_MINE,lastMPISide_MINE, SideToElem, ElemToSide, sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, SideID2, lastSideID, firstSideID_wo_BC, ElemID, locSideID2, flip, iVar
REAL                            :: gradUinside(3,PP_nVar), gradUinNorm(PP_nVar)
REAL                            :: gradWeight, minDiff(PP_nVar), maxDiff(PP_nVar), limiter(PP_nVar), w, a, b
!===================================================================================================================================

! Set the side range according to MPI or no MPI
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides (where the local proc is master)
  firstSideID_wo_BC = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  lastSideID = lastInnerSide
END IF

DO SideID=firstSideID_wo_BC,lastSideID
  FV_gradU_side(:,0,0,SideID) = U_master(:,0,0,SideID) - U_slave(:,0,0,SideID)
END DO

! 2. Compute the gradients at the boundary conditions: 1..nBCSides
IF (.NOT.doMPISides) THEN
  DO SideID=firstBCSide,lastBCSide
    gradUinside = 0.
    ElemID     = SideToElem(S2E_ELEM_ID,SideID)  !element is always master (slave=boundary ghost cell)
    maxDiff = 0.
    minDiff = 0.
#if PP_dim == 3
    DO locSideID2=1,6
#else
    DO locSideID2=2,5
#endif
      SideID2=ElemToSide(E2S_SIDE_ID,locSideID2,ElemID)
      flip=ElemToSide(E2S_FLIP,locSideID2,ElemID)
      IF (SideID2.LE.lastBCSide) CYCLE
      maxDiff = MAX(maxDiff,-FV_gradU_side(:,0,0,SideID2))
      minDiff = MIN(minDiff,-FV_gradU_side(:,0,0,SideID2))
      gradWeight = 1/VECNORM(FV_dx_master(:,0,0,SideID2)-FV_dx_slave(:,0,0,SideID2))**2
      gradUinside(1,:) = gradUinside(1,:) + gradWeight*FV_SysSol_BC(1,SideID2)*FV_gradU_side(:,0,0,SideID2)
      gradUinside(2,:) = gradUinside(2,:) + gradWeight*FV_SysSol_BC(2,SideID2)*FV_gradU_side(:,0,0,SideID2)
#if PP_dim == 3
      gradUinside(3,:) = gradUinside(3,:) + gradWeight*FV_SysSol_BC(3,SideID2)*FV_gradU_side(:,0,0,SideID2)
#endif
      IF (flip.EQ.0) THEN !master
        maxDiff = MAX(maxDiff,-FV_gradU_side(:,0,0,SideID2))
        minDiff = MIN(minDiff,-FV_gradU_side(:,0,0,SideID2))
      ELSE !slave
        maxDiff = MAX(maxDiff,FV_gradU_side(:,0,0,SideID2))
        minDiff = MIN(minDiff,FV_gradU_side(:,0,0,SideID2))
      END IF
    END DO

    !limiting
    limiter = 1.
    w = (FV_VktK**3)/sJ(0,0,0,ElemID)
#if PP_dim == 3
    DO locSideID2=1,6
#else
    DO locSideID2=2,5
#endif
      SideID2=ElemToSide(E2S_SIDE_ID,locSideID2,ElemID)
      flip=ElemToSide(E2S_FLIP,locSideID2,ElemID)
      IF (SideID2.LE.lastBCSide) CYCLE
      DO iVar=1,PP_nVar
        SELECT CASE(LimiterTypeBC)
        CASE(0)
          limiter(iVar) = 0.
        CASE(1) !minmax
          IF (flip.EQ.0) THEN !master
            b = DOT_PRODUCT(FV_dx_master(:,0,0,SideID2),gradUinside(:,iVar))
          ELSE !slave
            b = DOT_PRODUCT(FV_dx_slave(:,0,0,SideID2),gradUinside(:,iVar))
          END IF
          IF (b.GT.0.) THEN
            limiter(iVar) = MIN(limiter(iVar),MIN(1.,maxDiff(iVar)/b))
          ELSE iF (b.LT.0.) THEN
            limiter(iVar) = MIN(limiter(iVar),MIN(1.,minDiff(iVar)/b))
          ELSE
            limiter(iVar) = MIN(limiter(iVar),1.)
          END IF
        CASE(4) !venkatakrishnan
          IF (flip.EQ.0) THEN !master
            b = DOT_PRODUCT(FV_dx_master(:,0,0,SideID2),gradUinside(:,iVar))
          ELSE !slave
            b = DOT_PRODUCT(FV_dx_slave(:,0,0,SideID2),gradUinside(:,iVar))
          END IF
          IF (b.GT.0.) THEN
            a = maxDiff(iVar)
            limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
          ELSE iF (b.LT.0.) THEN
            a = minDiff(iVar)
            limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
          ELSE
            limiter(iVar) = MIN(limiter(iVar),1.)
          END IF
        CASE(9) ! no limiter (central)
          ! already equal to 1
        CASE DEFAULT
          CALL abort(__STAMP__,'Limiter type not implemented for boundary conditions.')
        END SELECT
      END DO
    END DO
    gradUinNorm(:) = limiter(:)*(-gradUinside(1,:)*FV_dx_master(1,0,0,SideID) &
                                 -gradUinside(2,:)*FV_dx_master(2,0,0,SideID) &
                                 -gradUinside(3,:)*FV_dx_master(3,0,0,SideID))

    CALL GetBoundaryGrad(SideID,FV_gradU_side(:,0,0,SideID),gradUinNorm,&
          U_master(:,0,0,SideID),&
          NormVec(:,0,0,SideID),&
          Face_xGP(:,0,0,SideID))
  END DO
END IF

END SUBROUTINE CalcFVGradients

SUBROUTINE SumFVGradients()
!===================================================================================================================================
! Sum up and limit gradients in one element for fv reconstruction
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars                  ,ONLY: FV_gradU_side, FV_gradU_elem, FV_dx_master, FV_dx_slave, FV_SysSol_slave, FV_SysSol_master
USE MOD_FV_Vars                  ,ONLY: FV_VktK, LimiterType
USE MOD_Mesh_Vars                ,ONLY: nElems, ElemToSide, sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID, SideID, locSideID, flip, iVar
REAL                            :: gradWeight, maxDiff(PP_nVar), minDiff(PP_nVar), limiter(PP_nVar), a, b, w
!===================================================================================================================================

DO ElemID = 1, nElems
  FV_gradU_elem(:,:,ElemID) = 0.
  maxDiff = 0.
  minDiff = 0.
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    flip=ElemToSide(E2S_FLIP,locSideID,ElemID)
    gradWeight = 1/VECNORM(FV_dx_master(:,0,0,SideID)-FV_dx_slave(:,0,0,SideID))**2
    IF (flip.EQ.0) THEN !master
      maxDiff = MAX(maxDiff,-FV_gradU_side(:,0,0,SideID))
      minDiff = MIN(minDiff,-FV_gradU_side(:,0,0,SideID))
      FV_gradU_elem(1,:,ElemID) = FV_gradU_elem(1,:,ElemID) &
                                          - gradWeight*FV_SysSol_master(1,SideID)*FV_gradU_side(:,0,0,SideID)
      FV_gradU_elem(2,:,ElemID) = FV_gradU_elem(2,:,ElemID) &
                                          - gradWeight*FV_SysSol_master(2,SideID)*FV_gradU_side(:,0,0,SideID)
#if PP_dim == 3
      FV_gradU_elem(3,:,ElemID) = FV_gradU_elem(3,:,ElemID) &
                                          - gradWeight*FV_SysSol_master(3,SideID)*FV_gradU_side(:,0,0,SideID)
#endif   
    ELSE !slave
      maxDiff = MAX(maxDiff,FV_gradU_side(:,0,0,SideID))
      minDiff = MIN(minDiff,FV_gradU_side(:,0,0,SideID))
      FV_gradU_elem(1,:,ElemID) = FV_gradU_elem(1,:,ElemID) &
                                          + gradWeight*FV_SysSol_slave(1,SideID)*FV_gradU_side(:,0,0,SideID)
      FV_gradU_elem(2,:,ElemID) = FV_gradU_elem(2,:,ElemID) &
                                          + gradWeight*FV_SysSol_slave(2,SideID)*FV_gradU_side(:,0,0,SideID)
#if PP_dim == 3
      FV_gradU_elem(3,:,ElemID) = FV_gradU_elem(3,:,ElemID) &
                                          + gradWeight*FV_SysSol_slave(3,SideID)*FV_gradU_side(:,0,0,SideID)
#endif
    END IF
  END DO

!limiting
  limiter = 1.
  w = (FV_VktK**3)/sJ(0,0,0,ElemID)
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    flip=ElemToSide(E2S_FLIP,locSideID,ElemID)
    DO iVar=1,PP_nVar
      SELECT CASE(LimiterType)
      CASE(0)
        limiter(iVar) = 0.
      CASE(1) !minmax
        IF (flip.EQ.0) THEN !master
          b = DOT_PRODUCT(FV_dx_master(:,0,0,SideID),FV_gradU_elem(:,iVar,ElemID))
        ELSE !slave
          b = DOT_PRODUCT(FV_dx_slave(:,0,0,SideID),FV_gradU_elem(:,iVar,ElemID))
        END IF
        IF (b.GT.0.) THEN
          limiter(iVar) = MIN(limiter(iVar),MIN(1.,maxDiff(iVar)/b))
        ELSE iF (b.LT.0.) THEN
          limiter(iVar) = MIN(limiter(iVar),MIN(1.,minDiff(iVar)/b))
        ELSE
          limiter(iVar) = MIN(limiter(iVar),1.)
        END IF
      CASE(4) !venkatakrishnan
        IF (flip.EQ.0) THEN !master
          b = DOT_PRODUCT(FV_dx_master(:,0,0,SideID),FV_gradU_elem(:,iVar,ElemID))
        ELSE !slave
          b = DOT_PRODUCT(FV_dx_slave(:,0,0,SideID),FV_gradU_elem(:,iVar,ElemID))
        END IF
        IF (b.GT.0.) THEN
          a = maxDiff(iVar)
          limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
        ELSE iF (b.LT.0.) THEN
          a = minDiff(iVar)
          limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
        ELSE
          limiter(iVar) = MIN(limiter(iVar),1.)
        END IF
      CASE(9) ! no limiter (central)
        ! already equal to 1
      CASE DEFAULT
        CALL abort(__STAMP__,'Limiter type not implemented for unstructured meshes.')
      END SELECT
    END DO
  END DO
  FV_gradU_elem(1,:,ElemID) = limiter(:)*FV_gradU_elem(1,:,ElemID)
  FV_gradU_elem(2,:,ElemID) = limiter(:)*FV_gradU_elem(2,:,ElemID)
  FV_gradU_elem(3,:,ElemID) = limiter(:)*FV_gradU_elem(3,:,ElemID)
  
END DO

END SUBROUTINE SumFVGradients

END MODULE MOD_FV_Gradients