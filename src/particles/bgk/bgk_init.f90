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

MODULE MOD_BGK_Init
!===================================================================================================================================
! Initialization of BGK
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitBGK
  MODULE PROCEDURE InitBGK
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitBGK, DefineParametersBGK
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for BGK
!==================================================================================================================================
SUBROUTINE DefineParametersBGK()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("BGK")

CALL prms%CreateIntOption(    'Particles-BGK-CollModel',            'TODO-DEFINE-PARAMETER. Define BGK method used.\n'//&
                                                                    '1: Ellipsoidal statistical\n'//&
                                                                    '2: Shakov\n'//&
                                                                    '3: Standard BGK\n'//&
                                                                    '4: Unified\n')
CALL prms%CreateIntOption(    'Particles-ESBGK-Model',              'TODO-DEFINE-PARAMETER.\n'//&
                                                                    '1: Approximative\n'//&
                                                                    '2: Exact\n'//&
                                                                    '3: MetropolisHastings', '1')
CALL prms%CreateIntOption(    'Particles-SBGK-EnergyConsMethod',    'SBGK energy conservation scheme', '1')
CALL prms%CreateRealOption(   'Particles-UnifiedBGK-Ces',           'TODO-DEFINE-PARAMETER', '1000.0')
CALL prms%CreateLogicalOption('Particles-BGK-DoCellAdaptation',     'Enables octree cell refinement until the given number of'//&
                                                                    'particles is reached. Equal refinement in all three'//&
                                                                    'directions (x,y,z)','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-MinPartsPerCell',      'Define minimum number of particles per cell for octree cell'//&
                                                                    'refinement', '10')
CALL prms%CreateLogicalOption('Particles-BGK-DoAveraging',          'Enable moving average of variables for the calculation of'//&
                                                                    'the cell temperature for the rotational and vibrational'//&
                                                                    'relaxation frequency for molecules','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-DoAveragingCorrection','Enable reset of the moving average after a certain number'//&
                                                                    'of iterations, defined by -AveragingLength','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-AveragingLength',      'Number of iterations after which the moving average is reset',&
                                                                    '5')
CALL prms%CreateRealOption(   'Particles-BGK-Acceleration',         'TODO-DEFINE-PARAMETER', '-9.81')
CALL prms%CreateRealOption(   'Particles-BGK-SplittingDens',        'TODO-DEFINE-PARAMETER', '0.0')
CALL prms%CreateLogicalOption('Particles-BGK-SampAdapFac',          'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-DoVibRelaxation',      'Enable modelling of vibrational excitation','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-UseQuantVibEn',        'Enable quantized treatment of vibrational energy levels',  &
                                                                    '.FALSE.')
CALL prms%CreateLogicalOption('Particles-CoupledBGKDSMC',           'Perform a coupled DSMC-BGK simulation with a given number'//&
                                                                    'density as a switch parameter','.FALSE.')
CALL prms%CreateRealOption(   'Particles-BGK-DSMC-SwitchDens',      'Number density [1/m3] above which the FP method is used'//&
                                                                    'below which DSMC is performed.','0.0')

END SUBROUTINE DefineParametersBGK

SUBROUTINE InitBGK()
!===================================================================================================================================
!> Init of BGK Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Vars      ,ONLY: nSpecies, Species
USE MOD_DSMC_Vars          ,ONLY: SpecDSMC, DSMC
USE MOD_Globals_Vars       ,ONLY: Pi, BoltzmannConst
USE MOD_ReadInTools
USE MOD_BGK_Vars
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix
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
SWRITE(UNIT_stdOut,'(A)') ' INIT BGK Solver...'
ALLOCATE(SpecBGK(nSpecies))
DO iSpec=1, nSpecies
  ALLOCATE(SpecBGK(iSpec)%CollFreqPreFactor(nSpecies))
  DO iSpec2=1, nSpecies
    SpecBGK(iSpec)%CollFreqPreFactor(iSpec2)= 0.5*(SpecDSMC(iSpec)%DrefVHS + SpecDSMC(iSpec2)%DrefVHS)**2.0 &
        * SQRT(2.*Pi*BoltzmannConst*SpecDSMC(iSpec)%TrefVHS*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(Species(iSpec)%MassIC * Species(iSpec2)%MassIC))/SpecDSMC(iSpec)%TrefVHS**(-SpecDSMC(iSpec)%omegaVHS +0.5)
  END DO
END DO

DoBGKCellAdaptation = GETLOGICAL('Particles-BGK-DoCellAdaptation')
BGKMinPartPerCell = GETINT('Particles-BGK-MinPartsPerCell')
BGKCollModel = GETINT('Particles-BGK-CollModel')
ESBGKModel = GETINT('Particles-ESBGK-Model')         ! 1: Approximative, 2: Exact, 3: MetropolisHastings
BGKUnifiedCes = GETREAL('Particles-UnifiedBGK-Ces')
BGKDSMCSwitchDens = GETREAL('Particles-BGK-DSMC-SwitchDens','0.')
CoupledBGKDSMC = GETLOGICAL('Particles-CoupledBGKDSMC','.FALSE.')
IF (BGKCollModel.EQ.2) THEN
  SBGKEnergyConsMethod = GETINT('Particles-SBGK-EnergyConsMethod')
ELSE
  SBGKEnergyConsMethod = 1
END IF
IF (BGKUnifiedCes.EQ.1000.) THEN
  BGKUnifiedCes = 1. - (6.-2.*SpecDSMC(1)%omegaVHS)*(4.- 2.*SpecDSMC(1)%omegaVHS)/30.
END IF
BGKAveragingLength = GETINT('Particles-BGK-AveragingLength')
BGKDoAveraging = GETLOGICAL('Particles-BGK-DoAveraging')
BGKDoAveragingCorrect = GETLOGICAL('Particles-BGK-DoAveragingCorrection')
BGKUseQuantVibEn = GETLOGICAL('Particles-BGK-UseQuantVibEn')
IF (BGKDoAveraging) CALL BGK_init_Averaging()
BGKAcceleration = GETREAL('Particles-BGK-Acceleration')
BGKDoVibRelaxation = GETLOGICAL('Particles-BGK-DoVibRelaxation')
BGKSplittingDens = GETREAL('Particles-BGK-SplittingDens')

IF(DSMC%CalcQualityFactors) THEN
  ALLOCATE(BGK_QualityFacSamp(1:5,nElems))
  BGK_QualityFacSamp(1:5,1:nElems) = 0.0
END IF

BGKInitDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)') ' INIT BGK DONE!'

END SUBROUTINE InitBGK


SUBROUTINE BGK_init_Averaging()
!===================================================================================================================================
!> Building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
USE MOD_BGK_Vars   ,ONLY: ElemNodeAveraging, BGKAveragingLength, BGKDoAveragingCorrect
USE MOD_Mesh_Vars  ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
!===================================================================================================================================
ALLOCATE(ElemNodeAveraging(nElems))
DO iElem = 1, nElems
  ALLOCATE(ElemNodeAveraging(iElem)%Root)
  IF (BGKDoAveragingCorrect) THEN
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(5,BGKAveragingLength))
     ElemNodeAveraging(iElem)%Root%AverageValues = 0.0
  ELSE
    BGKAveragingLength = 1
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(5,1))
    ElemNodeAveraging(iElem)%Root%AverageValues = 0.0
  END IF
  ElemNodeAveraging(iElem)%Root%CorrectStep = 0
END DO

END SUBROUTINE BGK_init_Averaging

END MODULE MOD_BGK_Init
