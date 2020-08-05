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
!> Initialization of the Bhatnagar-Gross-Krook method
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitBGK
  MODULE PROCEDURE InitBGK
END INTERFACE

INTERFACE FinalizeBGK
  MODULE PROCEDURE FinalizeBGK
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitBGK, DefineParametersBGK, FinalizeBGK
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

CALL prms%CreateIntOption(    'Particles-BGK-CollModel',            'Select the BGK method:\n'//&
                                                                    '1: Ellipsoidal statistical (ESBGK)\n'//&
                                                                    '2: Shakov (SBGK)\n'//&
                                                                    '3: Standard BGK\n'//&
                                                                    '4: Unified \n')
CALL prms%CreateIntOption(    'Particles-ESBGK-Model',              'Select sampling method for the ESBGK target distribution '//&
                                                                    'function:\n'//&
                                                                    '1: Approximative\n'//&
                                                                    '2: Exact\n'//&
                                                                    '3: Metropolis-Hastings', '1')
CALL prms%CreateIntOption(    'Particles-BGK-MixtureModel',         'Select model for mixture transport properties:\n'//&
                                                                    '1: Wilke\n'//&
                                                                    '2: Brokaw\n'//&
                                                                    '3: Collision Integrals Kestin', '1')
CALL prms%CreateIntOption(    'Particles-SBGK-EnergyConsMethod',    'Select the SBGK energy conservation scheme:\n'//&
                                                                    '1: Method includes all particles for energy conservation\n'//&
                                                                    '2: Number of particles included in the conservation scheme '//&
                                                                    'depends on the number of relaxing particles', '1')
CALL prms%CreateRealOption(   'Particles-UnifiedBGK-Ces',           'Parameter C_ES for the Unified BGK scheme. The default '//&
                                                                    'value 1000 enables the automatic calculation to reproduce '//&
                                                                    'the correct Prandtl number for equilibrium gas flows','1000.0')
CALL prms%CreateLogicalOption('Particles-BGK-DoCellAdaptation',     'Enables octree cell refinement until the given number of '//&
                                                                    'particles is reached. Equal refinement in all three '//&
                                                                    'directions (x,y,z)','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-MinPartsPerCell',      'Define minimum number of particles per cell for octree '//&
                                                                    'cell refinement')
CALL prms%CreateLogicalOption('Particles-BGK-MovingAverage',        'Enable a moving average of variables for the calculation '//&
                                                                    'of the cell temperature for relaxation frequencies','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-MovingAverageLength',  'Use the moving average with a fixed array length, where '//&
                                                                    'the first values are dismissed and the last values updated '//&
                                                                    'with current iteration to account for for unsteady flows','1')
CALL prms%CreateRealOption(   'Particles-BGK-SplittingDens',        'Octree-refinement will only be performed above this number '//&
                                                                    'density', '0.0')
CALL prms%CreateLogicalOption('Particles-BGK-DoVibRelaxation',      'Enable modelling of vibrational excitation','.TRUE.')
CALL prms%CreateLogicalOption('Particles-BGK-UseQuantVibEn',        'Enable quantized treatment of vibrational energy levels',  &
                                                                    '.TRUE.')
CALL prms%CreateLogicalOption('Particles-CoupledBGKDSMC',           'Perform a coupled DSMC-BGK simulation with a given number '//&
                                                                    'density as a switch parameter','.FALSE.')
CALL prms%CreateRealOption(   'Particles-BGK-DSMC-SwitchDens',      'Number density [1/m3], above which the BGK method is used, '//&
                                                                    'below which DSMC is performed','0.0')

END SUBROUTINE DefineParametersBGK

SUBROUTINE InitBGK()
!===================================================================================================================================
!> Initialization of the Bhatnagar-Gross-Krook method
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_BGK_Vars
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: nElems, NGeo
USE MOD_Particle_Vars         ,ONLY: nSpecies, Species, VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC, RadialWeighting
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_init_octree
USE MOD_Globals_Vars          ,ONLY: Pi, BoltzmannConst
USE MOD_Basis                 ,ONLY: PolynomialDerivativeMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iSpec2
REAL                  :: delta_ij
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT BGK Solver...'
ALLOCATE(SpecBGK(nSpecies))
DO iSpec=1, nSpecies
  ALLOCATE(SpecBGK(iSpec)%CollFreqPreFactor(nSpecies))
  DO iSpec2=1, nSpecies
    IF (iSpec.EQ.iSpec2) THEN
      delta_ij = 1.0
    ELSE
      delta_ij = 0.0
    END IF
    SpecBGK(iSpec)%CollFreqPreFactor(iSpec2)= 0.5*(2.-delta_ij)*(SpecDSMC(iSpec)%DrefVHS + SpecDSMC(iSpec2)%DrefVHS)**2.0 &
        * SQRT(2.*Pi*BoltzmannConst*SpecDSMC(iSpec)%TrefVHS*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(Species(iSpec)%MassIC * Species(iSpec2)%MassIC))/SpecDSMC(iSpec)%TrefVHS**(-SpecDSMC(iSpec)%omegaVHS +0.5)
  END DO
END DO

BGKCollModel = GETINT('Particles-BGK-CollModel')
BGKMixtureModel = GETINT('Particles-BGK-MixtureModel')
! ESBGK options
ESBGKModel = GETINT('Particles-ESBGK-Model')         ! 1: Approximative, 2: Exact, 3: MetropolisHastings
! Shakov BGK options
IF (BGKCollModel.EQ.2) THEN
  SBGKEnergyConsMethod = GETINT('Particles-SBGK-EnergyConsMethod')
  IF(SBGKEnergyConsMethod.EQ.2) THEN
    IF(ANY(SpecDSMC(:)%InterID.GT.1)) THEN
      CALL abort(&
__STAMP__&
,' ERROR SBGK: The chosen energy conservation method for SBGK was not tested with molecules!')
    END IF
  END IF
ELSE
  SBGKEnergyConsMethod = 1
END IF
! Unified BGK options
BGKUnifiedCes = GETREAL('Particles-UnifiedBGK-Ces')
IF (BGKUnifiedCes.EQ.1000.) THEN
  BGKUnifiedCes = 1. - (6.-2.*SpecDSMC(1)%omegaVHS)*(4.- 2.*SpecDSMC(1)%omegaVHS)/30.
ELSE IF((BGKUnifiedCes.LT.-0.5).OR.(BGKUnifiedCes.GE.1.0)) THEN
  CALL abort(&
__STAMP__&
,' ERROR Unified BGK: The parameter C_ES has to be in the range of -0.5 <= C_ES < 1 !')
END IF
! Coupled BGK with DSMC, use a number density as limit above which BGK is used, and below which DSMC is used
CoupledBGKDSMC = GETLOGICAL('Particles-CoupledBGKDSMC')
IF(CoupledBGKDSMC) THEN
  BGKDSMCSwitchDens = GETREAL('Particles-BGK-DSMC-SwitchDens')
ELSE
  IF(RadialWeighting%DoRadialWeighting) RadialWeighting%PerformCloning = .TRUE.
END IF
! Octree-based cell refinement, up to a certain number of particles
DoBGKCellAdaptation = GETLOGICAL('Particles-BGK-DoCellAdaptation')
IF(DoBGKCellAdaptation) THEN
  BGKMinPartPerCell = GETINT('Particles-BGK-MinPartsPerCell')
  IF(.NOT.DSMC%UseOctree) THEN
    DSMC%UseOctree = .TRUE.
    IF(NGeo.GT.PP_N) CALL abort(&
__STAMP__&
,' Set PP_N to NGeo, otherwise the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
END IF
BGKSplittingDens = GETREAL('Particles-BGK-SplittingDens')
! Moving Average
BGKMovingAverage = GETLOGICAL('Particles-BGK-MovingAverage')
IF(BGKMovingAverage) THEN
  BGKMovingAverageLength = GETINT('Particles-BGK-MovingAverageLength')
  CALL BGK_init_MovingAverage()
  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    CALL abort(&
__STAMP__&
,' ERROR BGK Init: Moving average is neither implemented with radial weighting nor variable time step!')
  END IF
END IF
! Vibrational modelling
BGKDoVibRelaxation = GETLOGICAL('Particles-BGK-DoVibRelaxation')
BGKUseQuantVibEn = GETLOGICAL('Particles-BGK-UseQuantVibEn')

IF(DSMC%CalcQualityFactors) THEN
  ALLOCATE(BGK_QualityFacSamp(1:7,nElems))
  BGK_QualityFacSamp(1:7,1:nElems) = 0.0
END IF

BGKInitDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)') ' INIT BGK DONE!'

END SUBROUTINE InitBGK


SUBROUTINE BGK_init_MovingAverage()
!===================================================================================================================================
!> Initialization of the arrays for the sampling of the moving average
!===================================================================================================================================
! MODULES
USE MOD_BGK_Vars   ,ONLY: ElemNodeAveraging, BGKMovingAverageLength
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
  IF (BGKMovingAverageLength.GT.1) THEN
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(5,BGKMovingAverageLength))
     ElemNodeAveraging(iElem)%Root%AverageValues = 0.0
  ELSE
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(5,1))
    ElemNodeAveraging(iElem)%Root%AverageValues = 0.0
  END IF
  ElemNodeAveraging(iElem)%Root%CorrectStep = 0
END DO

END SUBROUTINE BGK_init_MovingAverage


SUBROUTINE FinalizeBGK()
!----------------------------------------------------------------------------------------------------------------------------------!
!> Deallocating BGK-specific variables, for a complete finalization FinalizeDSMC is also required
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_BGK_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(SpecBGK)
SDEALLOCATE(ElemNodeAveraging)
SDEALLOCATE(BGK_QualityFacSamp)

END SUBROUTINE FinalizeBGK

END MODULE MOD_BGK_Init
