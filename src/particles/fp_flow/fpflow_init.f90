!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_FPFlow_Init
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitFPFlow
  MODULE PROCEDURE InitFPFlow
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitFPFlow, FP_BuildTransGaussNums, DefineParametersFPFlow
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FP-Flow
!==================================================================================================================================
SUBROUTINE DefineParametersFPFlow()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("FP-Flow")

CALL prms%CreateIntOption(    'Particles-FP-CollModel',       'TODO-DEFINE-PARAMETER. Define FP method used.\n'//&
                                                              '1: ...\n'//&
                                                              '2: ...\n'//&
                                                              '3: ...\n')
CALL prms%CreateIntOption(    'Particles-ESFP-Model',         'TODO-DEFINE-PARAMETER.\n'//&
                                                              '1: ...\n'//&
                                                              '2: ...\n', '1')
CALL prms%CreateLogicalOption('Particles-FP-DoVibRelaxation', 'Enable modelling of vibrational excitation','.TRUE.')
CALL prms%CreateLogicalOption('Particles-FP-UseQuantVibEn',   'Enable quantized modelling of vibrational energy levels','.TRUE.')
CALL prms%CreateLogicalOption('Particles-FP-DoCellAdaptation','Enables octree cell refinement until the given number of '//&
                                                              'particles is reached. Equal refinement in all three '//&
                                                              'directions (x,y,z)','.FALSE.')
CALL prms%CreateIntOption(    'Particles-FP-MinPartsPerCell', 'Define minimum number of particles per cell for octree '//&
                                                              'cell refinement')
CALL prms%CreateLogicalOption('Particles-CoupledFPDSMC',      'Perform a coupled DSMC-FP simulation with a given number density'//&
                                                              'as a switch parameter','.FALSE.')
CALL prms%CreateRealOption(   'Particles-FP-DSMC-SwitchDens', 'Number density [1/m3] above which the FP method is used, below'//&
                                                              'which DSMC is performed.','0.0')

END SUBROUTINE DefineParametersFPFlow

SUBROUTINE InitFPFlow()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc               ,ONLY: PP_N
USE MOD_Mesh_Vars             ,ONLY: NGeo
USE MOD_Globals_Vars          ,ONLY: PI, BoltzmannConst
USE MOD_ReadInTools
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_init_octree
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies, Species
USE MOD_FPFlow_Vars
USE MOD_BGK_Vars              ,ONLY: DoBGKCellAdaptation, BGKMinPartPerCell
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iSpec, iSpec2
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT FP Solver...'

ALLOCATE(SpecFP(nSpecies))

DO iSpec = 1, nSpecies
  ALLOCATE(SpecFP(iSpec)%CollFreqPreFactor(nSpecies))
  DO iSpec2=1, nSpecies
    SpecFP(iSpec)%CollFreqPreFactor(iSpec2)= 0.5*(SpecDSMC(iSpec)%DrefVHS + SpecDSMC(iSpec2)%DrefVHS)**2.0 & 
        * SQRT(2.*PI*BoltzmannConst*SpecDSMC(iSpec)%TrefVHS*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(Species(iSpec)%MassIC * Species(iSpec2)%MassIC))/SpecDSMC(iSpec)%TrefVHS**(-SpecDSMC(iSpec)%omegaVHS +0.5)
  END DO
END DO

FPCollModel = GETINT('Particles-FP-CollModel')
ESFPModel = GETINT('Particles-ESFP-Model')
DoBGKCellAdaptation = GETLOGICAL('Particles-FP-DoCellAdaptation')
IF(DoBGKCellAdaptation) THEN
  BGKMinPartPerCell = GETINT('Particles-FP-MinPartsPerCell')
  IF(.NOT.DSMC%UseOctree) THEN
    DSMC%UseOctree = .TRUE.
    IF(NGeo.GT.PP_N) CALL abort(&
__STAMP__&
,' Set PP_N to NGeo, else, the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
END IF
FPDoVibRelaxation = GETLOGICAL('Particles-FP-DoVibRelaxation')
FPUseQuantVibEn = GETLOGICAL('Particles-FP-UseQuantVibEn')
CoupledFPDSMC = GETLOGICAL('Particles-CoupledFPDSMC')
IF(CoupledFPDSMC) FPDSMCSwitchDens = GETREAL('Particles-FP-DSMC-SwitchDens')

END SUBROUTINE InitFPFlow

SUBROUTINE FP_BuildTransGaussNums(nPart, iRanPart)
!===================================================================================================================================
! Performs FP Momentum Evaluation
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: sumiRan(3), varianceiRan(3)
INTEGER                        :: iLoop
!===================================================================================================================================
sumiRan(1:3) = 0.0
varianceiRan(1:3) = 0.0
DO iLoop = 1, nPart
  iRanPart(1,iLoop) = rnor()
  iRanPart(2,iLoop) = rnor()
  iRanPart(3,iLoop) = rnor()
  sumiRan(1:3) = sumiRan(1:3) + iRanPart(1:3,iLoop)
END DO
sumiRan(1:3) = sumiRan(1:3)/nPart
DO iLoop = 1, nPart
  iRanPart(1:3,iLoop) = iRanPart(1:3,iLoop)-sumiRan(1:3)
  varianceiRan(1:3) = varianceiRan(1:3) + iRanPart(1:3,iLoop)*iRanPart(1:3,iLoop)
END DO
varianceiRan(1:3) = SQRT(varianceiRan(1:3)/nPart)

DO iLoop = 1, nPart
  iRanPart(1:3,iLoop) = iRanPart(1:3,iLoop)/varianceiRan(1:3)
END DO

END SUBROUTINE FP_BuildTransGaussNums

END MODULE MOD_FPFlow_Init
