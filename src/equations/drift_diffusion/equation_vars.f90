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
MODULE MOD_Equation_Vars_FV
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
LOGICAL           :: doCalcSource             !< Swith to calculate a source term or not, automatically set by calcsource itself

#if !(USE_HDG)
REAL                 :: c_corr
REAL                 :: fDamping
REAL                 :: WaveLength                             !> wave length
INTEGER,ALLOCATABLE  :: nBCByType(:)          !< Number of sides for each boundary
INTEGER,ALLOCATABLE  :: BCSideID(:,:)         !< SideIDs for BC types
#endif

REAL,ALLOCATABLE     :: EFluid_GradSide(:)

LOGICAL              :: DoExactFlux
LOGICAL,ALLOCATABLE  ::isExactFluxInterFace(:)

LOGICAL              :: EquationInitIsDone_FV=.FALSE.!< Init switch

INTEGER           :: IniExactFunc_FV             !< Number of exact function used for initialization
INTEGER           :: IniRefState_FV              !< RefState for initialization
INTEGER           :: nRefState_FV                !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefState_FV(:,:)        !< reference state

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames_FV(1)=(/ CHARACTER(LEN=255) :: 'ElectronDensity'/)

!===================================================================================================================================
END MODULE MOD_Equation_Vars_FV
