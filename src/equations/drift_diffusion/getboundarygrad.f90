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
SUBROUTINE GetBoundaryGrad(SideID,gradU,UPrim_master,NormVec,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BC
USE MOD_Mesh_Vars_FV  ,ONLY: BoundaryType_FV
USE MOD_Equation_FV   ,ONLY: ExactFunc_FV
USE MOD_Equation_Vars_FV ,ONLY: IniExactFunc_FV, RefState_FV
USE MOD_TimeDisc_Vars, ONLY : time
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(OUT)  :: gradU       (PP_nVar_FV)
REAL,INTENT(IN)   :: UPrim_master(PP_nVar_FV)
REAL,INTENT(IN)   :: NormVec (3)
REAL,INTENT(IN)   :: Face_xGP(3)

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar_FV)
!==================================================================================================================================
BCType  = BoundaryType_FV(BC(SideID),BC_TYPE)
BCState = BoundaryType_FV(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  IF(BCState.EQ.0) THEN ! Determine the exact BC state
    CALL ExactFunc_FV(IniExactFunc_FV,time,0,Face_xGP,UPrim_boundary)
  ELSE
    UPrim_boundary = RefState_FV(:,BCState)
  END IF
  gradU = UPrim_master - UPrim_boundary

CASE(3) ! von Neumann
  gradU = 0.

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
        'no BC defined in DVM/getboundarygrad.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad

END MODULE MOD_GetBoundaryGrad
#endif /*USE_FV*/
