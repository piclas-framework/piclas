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


#if (PP_TimeDiscMethod==700)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStep_DVM
!===================================================================================================================================

CONTAINS


SUBROUTINE TimeStep_DVM()
!===================================================================================================================================
! DUGKS/Exp diff DVM timestep with finite volumes
!===================================================================================================================================
! MODULES
USE MOD_FV_Vars               ,ONLY: U_FV,Ut_FV
USE MOD_TimeDisc_Vars         ,ONLY: dt,time,iter
USE MOD_FV                    ,ONLY: FV_main
USE MOD_DistFunc              ,ONLY: RescaleU, RescaleInit, ForceStep
USE MOD_Equation_Vars_FV      ,ONLY: IniExactFunc_FV, DVMAccel, DVMColl, DVMMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! initial rescaling if initialized with non-equilibrium flow
IF (DVMColl.AND.DVMMethod.GT.0..AND.iter.EQ.0) THEN
  IF (IniExactFunc_FV.EQ.4.OR.IniExactFunc_FV.EQ.5.OR.IniExactFunc_FV.EQ.6) CALL RescaleInit(dt)
END IF

IF (ANY(DVMAccel.NE.0.)) CALL ForceStep(dt)

IF (DVMColl) CALL RescaleU(1,dt)  ! ftilde -> fchapeau2
CALL FV_main(time,time,doSource=.FALSE.)  ! fchapeau2 -> ftilde2 -> Ut = flux of f
IF (DVMColl.AND.DVMMethod.GT.0) CALL RescaleU(2,dt/2.)  ! fchapeau2 -> fchapeau
U_FV = U_FV + Ut_FV*dt        ! fchapeau -> ftilde

IF (ANY(DVMAccel.NE.0.)) CALL ForceStep(dt) !two times for strang splitting -> second order

END SUBROUTINE TimeStep_DVM


END MODULE MOD_TimeStep
#endif /*PP_TimeDiscMethod==700*/
