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

INTERFACE GetBoundaryGrad2
  MODULE PROCEDURE GetBoundaryGrad2
END INTERFACE

PUBLIC:: GetBoundaryGrad, GetBoundaryGrad2
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Computes the gradient at a boundary for Finite Volumes reconstruction.
!==================================================================================================================================
SUBROUTINE GetBoundaryGrad(SideID,gradU,UPrim_master,NormVec,Face_xGP,dx_Face)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc, RefState, DVMBGKModel
USE MOD_TimeDisc_Vars, ONLY : dt, time
USE MOD_DistFunc      ,ONLY: MaxwellDistribution, ShakhovDistribution, MacroValuesFromDistribution, MaxwellScattering
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(IN)   :: UPrim_master(PP_nVar)
REAL,INTENT(OUT)  :: gradU       (PP_nVar)
REAL,INTENT(IN)   :: NormVec (3)
REAL,INTENT(IN)   :: Face_xGP(3)
REAL,INTENT(IN)   :: dx_Face

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar), MacroVal(8), tau
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  IF(BCState.EQ.0) THEN ! Determine the exact BC state
    CALL ExactFunc(IniExactFunc,time,0,Face_xGP,UPrim_boundary)
  ELSE
    CALL MaxwellDistribution(RefState(:,BCState),UPrim_boundary)
  END IF
  gradU = (UPrim_master - UPrim_boundary) / dx_Face

CASE(3) ! specular reflection
  CALL MacroValuesFromDistribution(MacroVal,UPrim_master,dt/2.,tau,2)
  MacroVal(2:4) = MacroVal(2:4) - 2.*DOT_PRODUCT(NormVec(1:3),MacroVal(2:4))*NormVec(1:3)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  gradU = (UPrim_master - UPrim_boundary) /dx_Face

CASE(4) ! diffusive
  MacroVal(:) = RefState(:,BCState)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  CALL MaxwellScattering(UPrim_boundary,UPrim_master,NormVec,2,dt/2.)
  gradU = (UPrim_master - UPrim_boundary) /dx_Face

CASE(5) !constant static pressure+temperature inlet
  CALL MacroValuesFromDistribution(MacroVal,UPrim_master,dt/2.,tau,2)
  MacroVal(1)=RefState(1,BCState)
  MacroVal(5)=RefState(5,BCState)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  gradU = (UPrim_master - UPrim_boundary) /dx_Face

CASE(6) !constant static pressure outlet
  CALL MacroValuesFromDistribution(MacroVal,UPrim_master,dt/2.,tau,2)
  MacroVal(5)=RefState(5,BCState)*RefState(1,BCState)/MacroVal(1)
  CALL MaxwellDistribution(MacroVal,UPrim_boundary)
  gradU = (UPrim_master - UPrim_boundary) /dx_Face

CASE(7) !open outlet
  CALL MacroValuesFromDistribution(MacroVal,UPrim_master,dt/2.,tau,2)
  SELECT CASE (DVMBGKModel)
    CASE(1)
      CALL MaxwellDistribution(MacroVal,UPrim_boundary)
    CASE(2)
      CALL ShakhovDistribution(MacroVal,UPrim_boundary)
    CASE DEFAULT
      CALL abort(__STAMP__,'DVM BGK Model not implemented.',999,999.)
  END SELECT
  gradU = (UPrim_master - UPrim_boundary) /dx_Face

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad

!==================================================================================================================================
!> Computes the gradient at a boundary for Finite Volumes reconstruction (2nd order version).
!==================================================================================================================================
SUBROUTINE GetBoundaryGrad2(SideID,gradU,gradUinside,UPrim_master,NormVec,Face_xGP,dx_Face)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc, RefState, DVMBGKModel, DVMVelos, DVMnVelos
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

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar), fplus(1:PP_nVar)!, f0incoming(1:PP_nVar), f0outgoing(1:PP_nVar)
REAL    :: MacroVal(8), tau, vnormal
INTEGER :: iVel, jVel, kVel, upos
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!


CASE(4) ! diffusive
  fplus(:)=UPrim_master(:)-dx_Face*gradUinside(:)
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
      gradU(upos) = (UPrim_master(upos) - UPrim_boundary(upos))/dx_Face
    ELSE
      ! gradU(upos) = (UPrim_master(upos) - f0outgoing(upos))/dx_Face/2.
      gradU(upos) = gradUinside(upos)
    END IF
  END DO; END DO; END DO


CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
        'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad2
  
END MODULE MOD_GetBoundaryGrad
#endif /*USE_FV*/
