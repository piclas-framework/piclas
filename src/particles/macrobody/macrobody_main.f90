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
USE MOD_MacroBody_Vars ,ONLY: nMacroParticle, MacroPart, UseMacroPart
USE MOD_MacroBody_Vars ,ONLY: MacroPartFluxesEnabled, MacroPartAccelerationEnabled
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                    :: iMP
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.UseMacroPart) RETURN

MacroPart(:)%center(1) = MacroPart(:)%center(1) + MacroPart(:)%velocity(1)*dt
MacroPart(:)%center(2) = MacroPart(:)%center(2) + MacroPart(:)%velocity(2)*dt
MacroPart(:)%center(3) = MacroPart(:)%center(3) + MacroPart(:)%velocity(3)*dt
IF(MacroPartAccelerationEnabled) THEN
  DO iMP=1,nMacroParticle
    MacroPart(iMP)%velocity(1:6) = MacroPart(iMP)%velocity(1:6) + MacroPart(iMP)%RHS(1:6)
    MacroPart(iMP)%RHS(1:6)=0.
  END DO
END IF
IF(MacroPartFluxesEnabled) THEN
  DO iMP=1,nMacroParticle
    MacroPart(iMP)%radius = MacroPart(iMP)%radius + MacroPart(iMP)%RHS(7)
    MacroPart(iMP)%temp   = MacroPart(iMP)%temp   + MacroPart(iMP)%RHS(8)
    MacroPart(iMP)%mass   = MacroPart(iMP)%mass   + MacroPart(iMP)%RHS(9)
    MacroPart(iMP)%RHS(7:9)=0.
  END DO
END IF

END SUBROUTINE MacroBody_main


END MODULE MOD_MacroBody
