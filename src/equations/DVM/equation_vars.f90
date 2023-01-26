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
LOGICAL           :: doCalcSource             !< Swith to calculate a source term or not, automatically set by calcsource itself
REAL              :: AdvVel(3)                !< Advection velocity
REAL              :: DiffC                    !< Diffusion constant
INTEGER           :: IniExactFunc             !< Number of exact function used for initialization
INTEGER           :: IniRefState              !< RefState for initialization
INTEGER           :: nRefState                !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefState(:,:)        !< reference state

REAL              :: Pi

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)       !< Buffer array for BC data
INTEGER,ALLOCATABLE  :: nBCByType(:)          !< Number of sides for each boundary
INTEGER,ALLOCATABLE  :: BCSideID(:,:)         !< SideIDs for BC types

INTEGER              :: DVMnVelos(3)
INTEGER              :: DVMBGKModel
INTEGER              :: DVMMethod
INTEGER              :: DVMVeloDisc
REAL                 :: DVMGHTemp(3)
INTEGER              :: DVMNewtDeg(3)
REAL                 :: DVMVeloMin(3)
REAL                 :: DVMVeloMax(3)
INTEGER              :: DVMDim
REAL, ALLOCATABLE    :: DVMVelos(:,:)
REAL, ALLOCATABLE    :: DVMWeights(:,:)
REAL                 :: DVMForce(3)

REAL                 :: c_corr
REAL                 :: fDamping

REAL,ALLOCATABLE     :: DVMMomentSave(:,:)

TYPE tSpeciesData
  REAL            :: omegaVHS
  REAL            :: T_Ref
  REAL            :: d_Ref
  REAL            :: Internal_DOF
  REAL            :: mu_Ref
  REAL            :: Mass
  REAL            :: R_S
  REAL            :: Prandtl
END TYPE tSpeciesData

TYPE(tSpeciesData) :: DVMSpeciesData

CHARACTER(LEN=255),DIMENSION(9),PARAMETER :: StrVarNames = (/ CHARACTER(LEN=255) :: 'Density', &
                                                                                    'VeloX', &
                                                                                    'VeloY', &
                                                                                    'VeloZ', &
                                                                                    'Temp', &
                                                                                    'HeatfluxX', &
                                                                                    'HeatfluxY', &
                                                                                    'HeatfluxZ', &
                                                                                    'RelaxFact'/)

LOGICAL              :: EquationInitIsDone=.FALSE.!< Init switch
LOGICAL              :: DoExactFlux
LOGICAL,ALLOCATABLE  ::isExactFluxInterFace(:)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
