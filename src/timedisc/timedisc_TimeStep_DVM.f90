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


#if (PP_TimeDiscMethod==600)
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
USE MOD_FV_Vars               ,ONLY: U,Ut
USE MOD_TimeDisc_Vars         ,ONLY: dt,time
USE MOD_FV                    ,ONLY: FV_main
USE MOD_DistFunc              ,ONLY: RescaleU, RescaleInit, ForceStep
USE MOD_Equation_Vars         ,ONLY: IniExactFunc, DVMForce
USE MOD_Mesh_Vars     ,ONLY: Elem_xGP

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
IF (IniExactFunc.EQ.4.AND.time.EQ.0.) CALL RescaleInit(dt) ! initial rescaling if simulation initialized with non-equilibrium flow

IF (ANY(DVMForce.NE.0.)) CALL ForceStep(dt)
! print*, Ut(95,0,0,0,31)
CALL RescaleU(1,dt)  ! ftilde -> fchapeau2
CALL FV_main(time,time,doSource=.FALSE.)  ! fchapeau2 -> ftilde2 -> Ut = flux de f
! print*, Ut(95,0,0,0,31)
! print*, Elem_xGP(1,0,0,0,31)
! ! read*
CALL RescaleU(2,dt/2.)  ! fchapeau2 -> fchapeau
U = U + Ut*dt        ! fchapeau -> ftilde

IF (ANY(DVMForce.NE.0.)) CALL ForceStep(dt) !two times for strang splitting -> second order

END SUBROUTINE TimeStep_DVM


END MODULE MOD_TimeStep
#endif /*PP_TimeDiscMethod==600*/
