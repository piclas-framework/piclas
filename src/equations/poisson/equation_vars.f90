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

! needed for various stuff (compilation)
REAL              :: c_corr
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
REAL,ALLOCATABLE  :: Et(:,:,:,:,:) ! temporal derivative of E
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
INTEGER           :: nRefState     !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefState(:,:) !< refstates in primitive variables (as read from ini file)

! Special BC with linear potential ramp (constant in time)
REAL              :: LinPhiBasePoint(3)
REAL              :: LinPhiNormal(3)
REAL              :: LinPhiHeight
REAL              :: LinPhi

#if defined(PARTICLES)
! Special BC with floating potential that is defined by the absorbed power of the charged particles
REAL              :: CoupledPowerPotential(3)   ! (/min, start, max/) electric potential at all BoundaryType = (/2,2/)
REAL              :: CoupledPowerTarget         ! Target input power at all BoundaryType = (/2,2/)
LOGICAL           :: CalcPCouplElectricPotential! Switch calculation on/off
#endif /*defined(PARTICLES)*/
!===================================================================================================================================
END MODULE MOD_Equation_Vars
