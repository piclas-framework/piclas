!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

!===================================================================================================================================
!> Main Module for macroscopic bodies inside particle domain
!===================================================================================================================================
MODULE MOD_MacroBody
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC :: MacroBody_main
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> main routine for macroscopic bodies
!===================================================================================================================================
SUBROUTINE MacroBody_main()
! MODULES                                                                                                                          !
USE MOD_TimeDisc_Vars  ,ONLY: dt
USE MOD_MacroBody_Vars ,ONLY: nMacroBody, UseMacroBody, MacroSphere
USE MOD_MacroBody_Vars ,ONLY: MacroBodyFluxesEnabled, MacroBodyAccelerationEnabled
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                    :: iMB
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.UseMacroBody) RETURN

MacroSphere(:)%center(1) = MacroSphere(:)%center(1) + MacroSphere(:)%velocity(1)*dt
MacroSphere(:)%center(2) = MacroSphere(:)%center(2) + MacroSphere(:)%velocity(2)*dt
MacroSphere(:)%center(3) = MacroSphere(:)%center(3) + MacroSphere(:)%velocity(3)*dt
IF(MacroBodyAccelerationEnabled) THEN
  DO iMB=1,nMacroBody
    MacroSphere(iMB)%velocity(1:6) = MacroSphere(iMB)%velocity(1:6) + MacroSphere(iMB)%RHS(1:6)
    MacroSphere(iMB)%RHS(1:6)=0.
  END DO
END IF
IF(MacroBodyFluxesEnabled) THEN
  DO iMB=1,nMacroBody
    MacroSphere(iMB)%radius = MacroSphere(iMB)%radius + MacroSphere(iMB)%RHS(7)
    MacroSphere(iMB)%temp   = MacroSphere(iMB)%temp   + MacroSphere(iMB)%RHS(8)
    MacroSphere(iMB)%mass   = MacroSphere(iMB)%mass   + MacroSphere(iMB)%RHS(9)
    MacroSphere(iMB)%RHS(7:9)=0.
  END DO
END IF

END SUBROUTINE MacroBody_main


END MODULE MOD_MacroBody
