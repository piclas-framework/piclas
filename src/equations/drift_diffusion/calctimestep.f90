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

MODULE MOD_CalcTimeStep
!===================================================================================================================================
! Low-Storage Runge-Kutta integration of degree 3 for one step.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!===================================================================================================================================

CONTAINS

FUNCTION CALCTIMESTEP()
!===================================================================================================================================
! Calculate the time step for the current update of U for the drift-diffusion equation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                         :: CalcTimeStep
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL abort(__STAMP__,'Use ManualTimeStep!')

END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
