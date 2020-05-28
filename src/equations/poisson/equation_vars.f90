!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
REAL              :: IniWavenumber(3) ! wavenumbers in 3 directions (sinus periodic with exactfunc=6)
INTEGER           :: IniExactFunc
REAL              :: IniCenter(3)
REAL              :: IniAmplitude
REAL              :: IniHalfwidth
REAL              :: ACfrequency,ACamplitude ! used for alternate current BC

! needed for various stuff (compilation)
REAL              :: c_corr
REAL              :: c_corr2    !c_corr^2
REAL              :: c_corr_c   !c_corr*c
REAL              :: c_corr_c2  !c_corr*c^2
REAL              :: eta_c      !(c_corr -1 )*c
REAL              :: fDamping
LOGICAL           :: DoParabolicDamping

REAL,ALLOCATABLE  :: chitens(:,:,:,:,:,:)    ! diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitensInv(:,:,:,:,:,:) ! inverse of diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitens_face(:,:,:,:,:) ! diffusion 3x3 tensor on each face gausspoint

LOGICAL           :: EquationInitIsDone=.FALSE.
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
REAL,ALLOCATABLE  :: E(:,:,:,:,:)
! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE  :: nBCByType(:)
INTEGER,ALLOCATABLE  :: BCSideID(:,:)
! can specify BC state
CHARACTER(LEN=255):: BCStateFile

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames(4)=(/ CHARACTER(LEN=255) :: 'Phi'           , &
                                                                                     'ElectricFieldX', &
                                                                                     'ElectricFieldY', &
                                                                                     'ElectricFieldZ'/)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
