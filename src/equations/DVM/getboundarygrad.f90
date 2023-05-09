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
SUBROUTINE GetBoundaryGrad(SideID,gradU,gradUinside,UPrim_master,NormVec,Face_xGP,dx_Face,dx_ElemIn)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc, RefState, DVMBGKModel, DVMVelos, DVMnVelos, DVMSpeciesData, DVMVeloDisc, DVMVeloMax, DVMVeloMin
USE MOD_TimeDisc_Vars, ONLY : dt, time
USE MOD_DistFunc      ,ONLY: MaxwellDistribution, ShakhovDistribution, MacroValuesFromDistribution, MaxwellScattering
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(IN)   :: UPrim_master(PP_nVar)
REAL,INTENT(OUT)  :: gradU       (PP_nVar)
REAL,INTENT(IN)   :: gradUinside (PP_nVar)
REAL,INTENT(IN)   :: NormVec (3)
REAL,INTENT(IN)   :: Face_xGP(3)
REAL,INTENT(IN)   :: dx_Face
REAL,INTENT(IN)   :: dx_ElemIn

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar), fplus(1:PP_nVar)!, f0incoming(1:PP_nVar), f0outgoing(1:PP_nVar)
REAL    :: MacroVal(8), tau, vnormal
INTEGER :: iVel, jVel, kVel, upos, upos_sp
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  fplus(:)=UPrim_master(:)-gradUinside(:)*dx_Face/dx_ElemIn
  ! f0outgoing(:)=UPrim_master(:)-2*dx_Face*gradUinside(:)
  MacroVal(:) = RefState(:,BCState)
  IF(BCState.EQ.0) THEN ! Determine the exact BC state
    CALL ExactFunc(IniExactFunc,time,0,Face_xGP,UPrim_boundary)
  ELSE
    CALL MaxwellDistribution(RefState(:,BCState),UPrim_boundary)
  END IF
  ! f0incoming=2*UPrim_boundary-UPrim_master
  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
    IF (vnormal.LT.0.) THEN
      ! gradU(upos) = (UPrim_master(upos) - f0incoming(upos))/dx_Face/2.
      gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos))
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*(UPrim_master(PP_nVar/2+upos) - UPrim_boundary(PP_nVar/2+upos))
      END IF
    ELSE
      ! gradU(upos) = (UPrim_master(upos) - f0outgoing(upos))/dx_Face/2.
      gradU(upos) = 2.*gradUinside(upos)*dx_Face/dx_ElemIn
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*gradUinside(PP_nVar/2+upos)*dx_Face/dx_ElemIn
      END IF
    END IF
  END DO; END DO; END DO


CASE(3) ! specular reflection
  UPrim_boundary(:)=UPrim_master(:)-gradUinside(:)*dx_Face/dx_ElemIn
  IF(DVMVeloDisc.EQ.2.AND.ANY((DVMVeloMin+DVMVeloMax).NE.0.)) THEN
    CALL abort(__STAMP__,'Specular reflection only implemented for zero-centered velocity grid')
  END IF
  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
    ! weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
    IF (ABS(NormVec(1)).EQ.1.) THEN !x-perpendicular boundary
      upos_sp=(DVMnVelos(1)+1-iVel)+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    ELSE IF (ABS(NormVec(2)).EQ.1.) THEN !y-perpendicular boundary
      upos_sp=iVel+(DVMnVelos(2)-jVel)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    ELSE IF (ABS(NormVec(3)).EQ.1.) THEN !z-perpendicular boundary
      upos_sp=iVel+(jVel-1)*DVMnVelos(1)+(DVMnVelos(3)-kVel)*DVMnVelos(1)*DVMnVelos(2)
    ELSE
      CALL abort(__STAMP__,'Specular reflection only implemented for boundaries perpendicular to velocity grid')
    END IF
    IF (vnormal.LT.0.) THEN !inflow
      gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos_sp))
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*(UPrim_master(PP_nVar/2+upos) - UPrim_boundary(PP_nVar/2+upos_sp))
      END IF
    ELSE
      gradU(upos) = 2.*gradUinside(upos)*dx_Face/dx_ElemIn
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*gradUinside(PP_nVar/2+upos)*dx_Face/dx_ElemIn
      END IF
    END IF
  END DO; END DO; END DO


CASE(4) ! diffusive
  fplus(:)=UPrim_master(:)-gradUinside(:)*dx_Face/dx_ElemIn
  ! f0outgoing(:)=UPrim_master(:)-2*dx_Face*gradUinside(:)
  MacroVal(:) = RefState(:,BCState)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  CALL MaxwellScattering(UPrim_boundary,fplus,NormVec,2,dt/2.)
  ! f0incoming=2*UPrim_boundary-UPrim_master
  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
    vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
    IF (vnormal.LT.0.) THEN
      ! gradU(upos) = (UPrim_master(upos) - f0incoming(upos))/dx_Face/2.
      gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos))
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*(UPrim_master(PP_nVar/2+upos) - UPrim_boundary(PP_nVar/2+upos))
      END IF
    ELSE
      ! gradU(upos) = (UPrim_master(upos) - f0outgoing(upos))/dx_Face/2.
      gradU(upos) = 2.*gradUinside(upos)*dx_Face/dx_ElemIn
      IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
        gradU(PP_nVar/2+upos) = 2.*gradUinside(PP_nVar/2+upos)*dx_Face/dx_ElemIn
      END IF
    END IF
  END DO; END DO; END DO


CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
        'no BC defined in DVM/getboundarygrad.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad
  
END MODULE MOD_GetBoundaryGrad
#endif /*USE_FV*/
