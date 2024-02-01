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
                                                                    '3: Standard BGK')
CALL prms%CreateIntOption(    'Particles-ESBGK-Model',              'Select sampling method for the ESBGK target distribution '//&
                                                                    'function:\n'//&
                                                                    '1: Approximative\n'//&
                                                                    '2: Exact\n'//&
                                                                    '3: Metropolis-Hastings', '1')
CALL prms%CreateIntOption(    'Particles-BGK-MixtureModel',         'Select model for mixture transport properties:\n'//&
                                                                    '1: Wilke\n'//&
                                                                    '2: Collision Integrals Kestin', '1')
CALL prms%CreateLogicalOption('Particles-BGK-DoCellAdaptation',     'Enables octree cell refinement until the given number of '//&
                                                                    'particles is reached. Equal refinement in all three '//&
                                                                    'directions (x,y,z)','.FALSE.')
CALL prms%CreateIntOption(    'Particles-BGK-MinPartsPerCell',      'Define minimum number of particles per cell for octree '//&
                                                                    'cell refinement')
CALL prms%CreateLogicalOption('Particles-BGK-MovingAverage',        'Enable a moving average of variables for the calculation '//&
                                                                    'of the cell temperature for relaxation frequencies','.FALSE.')
CALL prms%CreateRealOption(   'Particles-BGK-MovingAverageFac',     'Use the moving average of moments M with  '//&
                                                                    'M^n+1=AverageFac*M+(1-AverageFac)*M^n','0.01')
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
USE MOD_Particle_Vars         ,ONLY: nSpecies, Species, DoVirtualCellMerge
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC, RadialWeighting, CollInf
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_init_octree
USE MOD_Globals_Vars          ,ONLY: Pi, BoltzmannConst
USE MOD_Basis                 ,ONLY: PolynomialDerivativeMatrix
#if USE_MPI
USE MOD_Particle_MPI_Vars     ,ONLY: DoParticleLatencyHiding
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iSpec2
LOGICAL               :: MoleculePresent
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT BGK Solver...'
MoleculePresent = .FALSE.
ALLOCATE(SpecBGK(nSpecies))
DO iSpec=1, nSpecies
  IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) MoleculePresent = .TRUE.
  ALLOCATE(SpecBGK(iSpec)%CollFreqPreFactor(nSpecies))
  ! Calculation of the prefacor of the collision frequency per species
  ! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 87f
  DO iSpec2=1, nSpecies
    SpecBGK(iSpec)%CollFreqPreFactor(iSpec2)= 4.*CollInf%dref(iSpec,iSpec2)**2.0 &
        * SQRT(Pi*BoltzmannConst*CollInf%Tref(iSpec,iSpec2)*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(2.*(Species(iSpec)%MassIC * Species(iSpec2)%MassIC)))/CollInf%Tref(iSpec,iSpec2)**(-CollInf%omega(iSpec,iSpec2) +0.5)
  END DO
END DO

BGKCollModel = GETINT('Particles-BGK-CollModel')
IF ((nSpecies.GT.1).AND.(BGKCollModel.GT.1)) THEN
  CALL abort(__STAMP__,'ERROR Multispec only with ESBGK model!')
END IF
BGKMixtureModel = GETINT('Particles-BGK-MixtureModel')
! ESBGK options for sampling: 1: Approximative, 2: Exact, 3: MetropolisHastings
ESBGKModel = GETINT('Particles-ESBGK-Model')

! Coupled BGK with DSMC, use a number density as limit above which BGK is used, and below which DSMC is used
CoupledBGKDSMC = GETLOGICAL('Particles-CoupledBGKDSMC')
IF(CoupledBGKDSMC) THEN
  IF (DoVirtualCellMerge) THEN
    CALL abort(__STAMP__,'Virtual cell merge not implemented for coupled DSMC-BGK simulations!')
  END IF
#if USE_MPI
  IF (DoParticleLatencyHiding) THEN
    CALL abort(__STAMP__,'Particle latency hiding not implemented for coupled DSMC-BGK simulations!')
  END IF
#endif /*USE_MPI*/
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
    IF(NGeo.GT.PP_N) CALL abort(__STAMP__,'Set PP_N to NGeo, otherwise the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
END IF
BGKSplittingDens = GETREAL('Particles-BGK-SplittingDens')

! Moving Average
BGKMovingAverage = GETLOGICAL('Particles-BGK-MovingAverage')
IF(BGKMovingAverage) THEN
  BGKMovingAverageFac= GETREAL('Particles-BGK-MovingAverageFac')
  IF ((BGKMovingAverageFac.LE.0.0).OR.(BGKMovingAverageFac.GE.1.)) CALL abort(__STAMP__,'Particles-BGK-MovingAverageFac must be between 0 and 1!')
  CALL BGK_init_MovingAverage()
  IF(nSpecies.GT.1) CALL abort(__STAMP__,'nSpecies >1 and molecules not implemented for BGK averaging!')
END IF

IF(MoleculePresent) THEN
  ! Vibrational modelling
  BGKDoVibRelaxation = GETLOGICAL('Particles-BGK-DoVibRelaxation')
  BGKUseQuantVibEn = GETLOGICAL('Particles-BGK-UseQuantVibEn')
END IF

IF(DSMC%CalcQualityFactors) THEN
  ALLOCATE(BGK_QualityFacSamp(1:9,nElems))
  BGK_QualityFacSamp(1:9,1:nElems) = 0.0
END IF

BGKInitDone = .TRUE.

LBWRITE(UNIT_stdOut,'(A)') ' INIT BGK DONE!'

END SUBROUTINE InitBGK


SUBROUTINE BGK_init_MovingAverage()
!===================================================================================================================================
!> Initialization of the arrays for the sampling of the moving average
!===================================================================================================================================
! MODULES
USE MOD_BGK_Vars   ,ONLY: ElemNodeAveraging, BGKCollModel
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
  IF (BGKCollModel.EQ.1) THEN
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(8))
  ELSE IF (BGKCollModel.EQ.2) THEN
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(5))
  ELSE IF (BGKCollModel.EQ.3) THEN
    ALLOCATE(ElemNodeAveraging(iElem)%Root%AverageValues(2))
  END IF
  ElemNodeAveraging(iElem)%Root%AverageValues = 0.0
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
SDEALLOCATE(BGK_QualityFacSamp)
IF(BGKMovingAverage) CALL DeleteElemNodeAverage()

END SUBROUTINE FinalizeBGK


SUBROUTINE DeleteElemNodeAverage()
!----------------------------------------------------------------------------------------------------------------------------------!
! Delete the pointer tree ElemNodeVol
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_BGK_Vars
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Mesh_Vars              ,ONLY: nElems
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem
!===================================================================================================================================
IF(.NOT.DSMC%UseOctree) RETURN
DO iElem=1,nElems
  CALL DeleteNodeAverage(ElemNodeAveraging(iElem)%Root)
  DEALLOCATE(ElemNodeAveraging(iElem)%Root)
END DO
DEALLOCATE(ElemNodeAveraging)
END SUBROUTINE DeleteElemNodeAverage


RECURSIVE SUBROUTINE DeleteNodeAverage(NodeAverage)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if the Node has subnodes and delete them
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_BGK_Vars
USE MOD_Particle_Vars         ,ONLY: Symmetry
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
TYPE (tNodeAverage), INTENT(INOUT)  :: NodeAverage
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     ::  iLoop, nLoop
!===================================================================================================================================
nLoop = 2**Symmetry%Order
IF(ASSOCIATED(NodeAverage%SubNode)) THEN
  DO iLoop = 1, nLoop
    CALL DeleteNodeAverage(NodeAverage%SubNode(iLoop))
  END DO
  DEALLOCATE(NodeAverage%SubNode)
END IF

END SUBROUTINE DeleteNodeAverage

END MODULE MOD_BGK_Init
