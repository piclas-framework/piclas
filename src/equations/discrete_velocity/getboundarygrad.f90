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

#if (USE_FV)
MODULE MOD_GetBoundaryGrad
!===================================================================================================================================
! Contains FillBoundary (which depends on the considered equation)
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

INTERFACE GetBoundaryGrad
  MODULE PROCEDURE GetBoundaryGrad
END INTERFACE

PUBLIC:: GetBoundaryGrad
!===================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes the gradient at a boundary for Finite Volumes reconstruction (2nd order version).
!==================================================================================================================================
SUBROUTINE GetBoundaryGrad(SideID,gradU,gradUinside,UPrim_master,NormVec,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals          ,ONLY: Abort
USE MOD_Mesh_Vars        ,ONLY: BoundaryType,BC
USE MOD_Equation_FV      ,ONLY: ExactFunc_FV
USE MOD_Equation_Vars_FV ,ONLY: IniExactFunc_FV, RefState, DVMBGKModel
USE MOD_Equation_Vars_FV ,ONLY: DVMVelos, DVMnVelos, DVMSpeciesData, DVMVeloDisc, DVMVeloMax, DVMVeloMin, DVMDim, Pi, DVMWeights
USE MOD_TimeDisc_Vars    ,ONLY : dt, time
USE MOD_DistFunc         ,ONLY: MaxwellDistribution, MacroValuesFromDistribution, MaxwellScattering
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(OUT)  :: gradU       (PP_nVar_FV)
REAL,INTENT(IN)   :: gradUinside (PP_nVar_FV)
REAL,INTENT(IN)   :: UPrim_master(PP_nVar_FV)
REAL,INTENT(IN)   :: NormVec (3)
REAL,INTENT(IN)   :: Face_xGP(3)

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar_FV), fplus(1:PP_nVar_FV)
REAL    :: MacroVal(14), tau, vnormal, vwall, Sin, Sout, MovTerm, WallDensity, weight
INTEGER :: iVel, jVel, kVel, upos, upos_sp
!==================================================================================================================================
BCType  = BoundaryType(BC(SideID),BC_TYPE)
BCState = BoundaryType(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  IF(BCState.EQ.0) THEN ! Determine the exact BC state
    CALL ExactFunc_FV(IniExactFunc_FV,time,0,Face_xGP,UPrim_boundary)
  ELSE
    CALL MaxwellDistribution(RefState(:,BCState),UPrim_boundary)
  END IF
  gradU = UPrim_master - UPrim_boundary


CASE(3) ! specular reflection
  UPrim_boundary(:)=UPrim_master(:)-gradUinside(:)
  IF(DVMVeloDisc.EQ.2.AND.ANY((DVMVeloMin+DVMVeloMax).NE.0.)) THEN
    CALL abort(__STAMP__,'Specular reflection only implemented for zero-centered velocity grid')
  END IF
  IF (BCState.NE.0) CALL abort(__STAMP__,'DVM specular bc with moving wall not working (yet?)') !THEN
  !   IF (DVMVeloDisc.NE.3) CALL abort(__STAMP__,'DVM specular bc error: moving wall only with Gauss-Hermite disc')
  !   CALL MacroValuesFromDistribution(MacroVal,UPrim_boundary(:),dt/2.,tau,2)
  !   Sin=0.
  !   Sout=0.
  !   DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  !     upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  !     weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  !     vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
  !     vwall = DVMVelos(iVel,1)*RefState(2,BCState) + DVMVelos(jVel,2)*RefState(3,BCState) + DVMVelos(kVel,3)*RefState(4,BCState)
  !     IF (vnormal.LT.0.) THEN !inflow
  !       Sin = Sin + 2.*weight*vwall*EXP(-(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)/DVMSpeciesData%R_S/MacroVal(5)/2.) &
  !       /DVMSpeciesData%R_S/MacroVal(5)/(2.*Pi*DVMSpeciesData%R_S*MacroVal(5))**(DVMDim/2.)
  !     ELSE IF (vnormal.GT.0.) THEN
  !       Sout = Sout + weight*UPrim_boundary(upos)
  !     ELSE
  !       Sout = Sout + 2.*weight*UPrim_boundary(upos)
  !     END IF
  !   END DO; END DO; END DO
  !   WallDensity = Sout/(1-Sin)
  ! END IF
  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
    ! vwall = DVMVelos(iVel,1)*RefState(2,BCState) + DVMVelos(jVel,2)*RefState(3,BCState) + DVMVelos(kVel,3)*RefState(4,BCState)
    IF (ABS(ABS(NormVec(1)) - 1.).LE.1e-6) THEN !x-perpendicular boundary
      upos_sp=(DVMnVelos(1)+1-iVel)+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    ELSE IF (ABS(ABS(NormVec(2)) - 1.).LE.1e-6) THEN !y-perpendicular boundary
      upos_sp=iVel+(DVMnVelos(2)-jVel)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    ELSE IF (ABS(ABS(NormVec(3)) - 1.).LE.1e-6) THEN !z-perpendicular boundary
      upos_sp=iVel+(jVel-1)*DVMnVelos(1)+(DVMnVelos(3)-kVel)*DVMnVelos(1)*DVMnVelos(2)
    ELSE
      CALL abort(__STAMP__,'Specular reflection only implemented for boundaries perpendicular to velocity grid')
    END IF
    IF (vnormal.LT.0.) THEN !inflow
      ! IF (BCState.NE.0) THEN
      !   MovTerm = 2.*WallDensity*EXP(-(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)/DVMSpeciesData%R_S/MacroVal(5)/2.) &
      !   /DVMSpeciesData%R_S/MacroVal(5)/(2.*Pi*DVMSpeciesData%R_S*MacroVal(5))**(DVMDim/2.)
      ! ELSE
      !   MovTerm = 0.
      ! END IF
      gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos_sp))! - MovTerm)
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar_FV/2+upos) = 2.*(UPrim_master(PP_nVar_FV/2+upos) - UPrim_boundary(PP_nVar_FV/2+upos_sp))! - MovTerm)
      END IF
    ELSE
      gradU(upos) = 2.*gradUinside(upos)
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar_FV/2+upos) = 2.*gradUinside(PP_nVar_FV/2+upos)
      END IF
    END IF
  END DO; END DO; END DO


CASE(4) ! diffusive
  fplus(:)=UPrim_master(:)-gradUinside(:)
  MacroVal(:) = RefState(:,BCState)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  CALL MaxwellScattering(UPrim_boundary,fplus,NormVec,2,dt/2.)
  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
    IF (vnormal.LT.0.) THEN
      gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos))
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar_FV/2+upos) = 2.*(UPrim_master(PP_nVar_FV/2+upos) - UPrim_boundary(PP_nVar_FV/2+upos))
      END IF
    ELSE
      gradU(upos) = 2.*gradUinside(upos)
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar_FV/2+upos) = 2.*gradUinside(PP_nVar_FV/2+upos)
      END IF
    END IF
  END DO; END DO; END DO

CASE(14) ! diffusive order 1
  MacroVal(:) = RefState(:,BCState)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  CALL MaxwellScattering(UPrim_boundary,UPrim_master,NormVec,2,dt/2.)
  gradU = 2.*(UPrim_master-UPrim_boundary)

CASE(7) ! open outlet
  gradU = 0.

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
        'no BC defined in DVM/getboundarygrad.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad

END MODULE MOD_GetBoundaryGrad
#endif /*USE_FV*/
