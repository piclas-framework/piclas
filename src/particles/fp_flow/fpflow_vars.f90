!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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
MODULE MOD_FPFlow_Vars
!===================================================================================================================================
! Contains the FP Flow variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL               :: FPInitDone = .FALSE.

INTEGER               :: FPCollModel
INTEGER               :: ESFPModel
LOGICAL               :: FPDoVibRelaxation
LOGICAL               :: FPUseQuantVibEn
REAL                  :: FPDSMCSwitchDens
LOGICAL               :: CoupledFPDSMC

TYPE tSpecFP
  REAL, ALLOCATABLE          ::  CollFreqPreFactor(:)
END TYPE

TYPE(tSpecFP), ALLOCATABLE :: SpecFP(:)

REAL, ALLOCATABLE     :: FP_QualityFacSamp(:,:)
INTEGER               :: FP_MeanRelaxFactorCounter
REAL                  :: FP_MeanRelaxFactor
REAL                  :: FP_MaxRelaxFactor
REAL                  :: FP_MaxRotRelaxFactor
REAL                  :: FP_PrandtlNumber

!===================================================================================================================================
END MODULE MOD_FPFlow_Vars
