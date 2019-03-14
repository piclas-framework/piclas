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

MODULE MOD_ESBGK_Init
!===================================================================================================================================
! Initialization of ESBGK
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitESBGK
  MODULE PROCEDURE InitESBGK
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitESBGK, DefineParametersBGK
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
CALL prms%CreateLogicalOption('Particles-BGK-DoAveraging',          'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-DoAveragingCorrection','TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-AveragingLength',      'TODO-DEFINE-PARAMETER', '5')
CALL prms%CreateRealOption(   'Particles-BGK-Acceleration',         'TODO-DEFINE-PARAMETER', '-9.81')
CALL prms%CreateRealOption(   'Particles-BGK-SplittingDens',        'TODO-DEFINE-PARAMETER', '0.0')
CALL prms%CreateLogicalOption('Particles-BGK-DoBGKCellSplitting',   'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-SampAdapFac',          'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-DoVibRelaxation',      'Enable modelling of vibrational excitation','.FALSE.')
CALL prms%CreateLogicalOption('Particles-BGK-UseQuantVibEn',        'Enable quantized treatment of vibrational energy levels',  &
                                                                    '.FALSE.')
CALL prms%CreateLogicalOption('Particles-CoupledBGKDSMC',           'Perform a coupled DSMC-BGK simulation with a given number'//&
                                                                    'density as a switch parameter','.FALSE.')
CALL prms%CreateRealOption(   'Particles-BGK-DSMC-SwitchDens',      'Number density [1/m3] above which the FP method is used'//&
                                                                    'below which DSMC is performed.','0.0')

END SUBROUTINE DefineParametersBGK

SUBROUTINE InitESBGK()
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
USE MOD_ESBGK_Vars
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iSpec2, iElem
REAL                  :: b
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT BGK Solver...'
ALLOCATE(SpecESBGK(nSpecies))
DO iSpec=1, nSpecies
  ALLOCATE(SpecESBGK(iSpec)%CollFreqPreFactor(nSpecies))
  DO iSpec2=1, nSpecies
    SpecESBGK(iSpec)%CollFreqPreFactor(iSpec2)= 0.5*(SpecDSMC(iSpec)%DrefVHS + SpecDSMC(iSpec2)%DrefVHS)**2.0 &
        * SQRT(2.*Pi*BoltzmannConst*SpecDSMC(iSpec)%TrefVHS*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(Species(iSpec)%MassIC * Species(iSpec2)%MassIC))/SpecDSMC(iSpec)%TrefVHS**(-SpecDSMC(iSpec)%omegaVHS +0.5)
  END DO
END DO

b = (0.5 - SpecDSMC(1)%omegaVHS)
ESBGKTempCorrectFact = (2.-SpecDSMC(1)%omegaVHS)**b * GAMMA(2.-SpecDSMC(1)%omegaVHS) &
                        / GAMMA(2.-SpecDSMC(1)%omegaVHS+b)

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
DoBGKCellSplitting = GETLOGICAL('Particles-BGK-DoBGKCellSplitting')
BGKSplittingDens = GETREAL('Particles-BGK-SplittingDens')
IF (DoBGKCellSplitting) THEN
  DoBGKCellAdaptation = .FALSE.
  ALLOCATE(ElemSplitCells(nElems))
  CALL DefineElementOrientation()
  DO iElem = 1, nElems
    ElemSplitCells(iElem)%Splitnum(1:3) = (/0,4,0/) 
  END DO 
END IF

IF(DSMC%CalcQualityFactors) THEN
  ALLOCATE(BGK_QualityFacSamp(nElems,1:4))
  BGK_QualityFacSamp(1:nElems,1:4) = 0.0
END IF

BGKInitDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)') ' INIT BGK DONE!'

END SUBROUTINE InitESBGK


SUBROUTINE DefineElementOrientation()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars    ,ONLY: nElems,NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Globals_Vars ,ONLY: Pi
USE MOD_ESBGK_Vars   ,ONLY: ElemSplitCells
USE MOD_Eval_xyz     ,ONLY: TensorProductInterpolation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL      :: Vec(3), ResVec1(3,3), ResVec2(3,3), DiffVec(3,3), MinVec(3), angle(3,3), vecmag(3)
INTEGER   :: iElem, ii, jj, firstDir1, firstDir2
!===================================================================================================================================
DO iElem = 1, nElems
  Vec(1:3) = (/-0.99,0.,0./)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,1))
  Vec(1:3) = (/0.99,0.,0./)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,1))
  Vec(1:3) = (/0.,-0.99,0./)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,2))
  Vec(1:3) = (/0.,0.99,0./)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,2))
  Vec(1:3) = (/0.,0.,-0.99/)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,3))
  Vec(1:3) = (/0.,0.,0.99/)
  CALL TensorProductInterpolation(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,3))

  DO ii=1,3; DO jj=1,3
    DiffVec(ii,jj) = ResVec1(ii,jj)-ResVec2(ii,jj)
  END DO; END DO
  DO ii=1,3
    vecmag(ii) = SQRT(DiffVec(1,ii)*DiffVec(1,ii)+DiffVec(2,ii)*DiffVec(2,ii)+DiffVec(3,ii)*DiffVec(3,ii))
  END DO
  DO ii=1,3; DO jj=1,3
    angle(ii,jj) = ACOS(DiffVec(ii,jj)/ vecmag(jj))
    IF (angle(ii,jj).GT.Pi/2.) THEN
      angle(ii,jj) = Pi - angle(ii,jj)
    END IF
  END DO; END DO

  DO ii=1,3
    MinVec(ii) = MINVAL(angle(:,ii))
  END DO
  firstDir1 = MINLOC(MinVec,1)
  firstDir2 = MINLOC(angle(:,firstDir1),1)
  ElemSplitCells(iElem)%CellOrientation(firstDir2) = firstDir1
  DO ii=1,3
    angle(ii,firstDir1) = 1000
    angle(firstDir2,ii) = 1000
  END DO
  DO ii=1,3
    MinVec(ii) = MINVAL(angle(:,ii))
  END DO
  firstDir1 = MINLOC(MinVec,1)
  firstDir2 = MINLOC(angle(:,firstDir1),1)
  ElemSplitCells(iElem)%CellOrientation(firstDir2) = firstDir1
  DO ii=1,3
    angle(ii,firstDir1) = 1000
    angle(firstDir2,ii) = 1000
  END DO
  DO ii=1,3
    MinVec(ii) = MINVAL(angle(:,ii))
  END DO
  firstDir1 = MINLOC(MinVec,1)
  firstDir2 = MINLOC(angle(:,firstDir1),1)
  ElemSplitCells(iElem)%CellOrientation(firstDir2) = firstDir1
END DO
END SUBROUTINE DefineElementOrientation


SUBROUTINE BGK_init_Averaging()
!===================================================================================================================================
!> Building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
USE MOD_ESBGK_Vars ,ONLY: ElemNodeAveraging, BGKAveragingLength, BGKDoAveragingCorrect
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

END MODULE MOD_ESBGK_Init
