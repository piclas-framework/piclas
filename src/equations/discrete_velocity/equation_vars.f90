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
! Contains the DVM variables
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
INTEGER           :: IniExactFunc_FV             !< Number of exact function used for initialization
INTEGER           :: IniRefState_FV              !< RefState for initialization
INTEGER           :: nRefState_FV                !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefState_FV(:,:,:)        !< reference state

REAL              :: Pi

#if !(USE_HDG)
REAL                 :: c_corr
REAL                 :: fDamping
REAL                 :: WaveLength                             !> wave length
INTEGER,ALLOCATABLE  :: nBCByType(:)          !< Number of sides for each boundary
INTEGER,ALLOCATABLE  :: BCSideID(:,:)         !< SideIDs for BC types
#endif

LOGICAL              :: DVMColl
INTEGER              :: DVMnSpecies
INTEGER              :: DVMnMacro=14
INTEGER              :: DVMnInnerE
INTEGER              :: DVMBGKModel
INTEGER              :: DVMMethod
INTEGER              :: DVMDim
REAL                 :: DVMForce(3)
INTEGER,ALLOCATABLE  :: DVMVeloDisc(:)
REAL,ALLOCATABLE     :: DVMGHTemp(:,:)
INTEGER,ALLOCATABLE  :: DVMNewtDeg(:,:)

REAL,ALLOCATABLE     :: DVMMomentSave(:,:,:)
REAL,ALLOCATABLE     :: DVMInnerESave(:,:,:)

TYPE tSpeciesData
  CHARACTER(LEN=64):: Name
  LOGICAL          :: DoOverwriteParameters ! Flag to read in parameters manually
  REAL            :: omegaVHS
  REAL            :: T_Ref
  REAL            :: d_Ref
  REAL            :: mu_Ref
  REAL            :: Mass
  REAL            :: Charge
  REAL            :: R_S
  INTEGER         :: InterID
  INTEGER         :: Xi_Rot
  INTEGER         :: nVar
  INTEGER         :: nVarReduced
  INTEGER         :: nVarErotStart
  INTEGER         :: nVelos(3)
  REAL            :: VeloMin(3)
  REAL            :: VeloMax(3)
  REAL,ALLOCATABLE:: Velos(:,:)
  REAL,ALLOCATABLE:: Weights(:,:)
  REAL            :: GHTemp(3)
  INTEGER         :: NewtDeg(3)
END TYPE tSpeciesData

TYPE(tSpeciesData),ALLOCATABLE :: DVMSpecData(:)

CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames_FV(:)

LOGICAL              :: WriteDVMSurfaceValues
REAL,ALLOCATABLE     :: DVMSurfaceValues(:,:,:,:)
INTEGER              :: nVarDVMSurf=4

LOGICAL              :: EquationInitIsDone_FV=.FALSE.!< Init switch
LOGICAL              :: DoExactFlux
LOGICAL,ALLOCATABLE  ::isExactFluxInterFace(:)
!===================================================================================================================================
END MODULE MOD_Equation_Vars_FV
