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
MODULE MOD_Equation_Vars
!===================================================================================================================================
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL              :: Pi
REAL              :: IniWavenumber(3) ! wavenumbers in 3 directions (sinus periodic with exactfunc=6)
INTEGER           :: IniExactFunc
REAL              :: IniCenter(3)
REAL              :: IniAmplitude
REAL              :: IniHalfwidth

REAL,ALLOCATABLE  :: chitens(:,:,:,:,:,:) !diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitensInv(:,:,:,:,:,:) ! inverse of diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitens_face(:,:,:,:,:) !diffusion 3x3 tensor on each face gausspoint

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames(3)=(/ CHARACTER(LEN=255) :: 'MagneticFieldX'           , &
                                                                                     'MagneticFieldY', &
                                                                                     'MagneticFieldZ'/)

LOGICAL           :: EquationInitIsDone=.FALSE.
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
REAL,ALLOCATABLE  :: B(:,:,:,:,:)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
