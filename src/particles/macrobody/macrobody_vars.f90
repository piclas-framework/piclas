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
!> Contains the variables necessary for MacroBodies inside particle Domain
!===================================================================================================================================
MODULE MOD_MacroBody_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tMacroParticle
  REAL    :: center(3)
  REAL    :: velocity(6)
  REAL    :: radius
  REAL    :: temp
  REAL    :: density
  REAL    :: mass
  REAL    :: RHS(1:9)
  REAL    :: momentumAcc
  REAL    :: transAcc
  REAL    :: vibAcc
  REAL    :: rotAcc
END TYPE

TYPE(tMacroParticle), ALLOCATABLE :: MacroPart(:)
INTEGER                           :: nMacroParticle
LOGICAL                           :: MacroPartFluxesEnabled
LOGICAL                           :: MacroPartAccelerationEnabled
LOGICAL                           :: MacroPartWriteElemData
LOGICAL                           :: UseMacroPart
LOGICAL,ALLOCATABLE               :: ElemHasMacroPart(:,:)
LOGICAL                           :: CalcMPVolumePortion

END MODULE MOD_MacroBody_Vars
