#include "boltzplatz.h"

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
PUBLIC :: InitESBGK, DefineParametersBGK, ESBGK_BuildTransGaussNumsEnCon
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

CALL prms%CreateRealOption(     'Particles-BGKAdaptTimeStep'  ,     'TODO-DEFINE-PARAMETER.BGK adaptation timestep [sec]?', '0.0')
CALL prms%CreateLogicalOption(  'Particles-BGK-DoCellAdaptation',   'Enable BGK Cell Adaptation','.FALSE.')
CALL prms%CreateIntOption(      'Particles-BGK-MinPartsPerCell'  ,  'Define minimum number of particles per cell', '10')
CALL prms%CreateIntOption(      'Particles-BGKCollModel'  ,         'TODO-DEFINE-PARAMETER. Define BGK method used.\n'//&
                                '1: ...\n'//&
                                '2: ...\n'//&
                                '3: ...\n'//&
                                '4: ...\n', '1')
CALL prms%CreateIntOption(      'Particles-ESBGKModel'  ,         'TODO-DEFINE-PARAMETER.\n'//&
                                '1: Approximative\n'//&
                                '2: Exact\n'//&
                                '3: MetropolisHastings', '1')
CALL prms%CreateRealOption(     'Particles-UnifiedBGK-Ces'  ,     'TODO-DEFINE-PARAMETER', '1000.0')
CALL prms%CreateLogicalOption(  'Particles-BGKDoAveraging',       'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-BGKDoAveragingCorrection',       'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateIntOption(      'Particles-BGKAveragingLength'  ,  'TODO-DEFINE-PARAMETER', '5')
CALL prms%CreateLogicalOption(  'Particles-BGKUseQuantVibEn',       'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateRealOption(     'Particles-BGKAcceleration'  ,     'TODO-DEFINE-PARAMETER', '-9.81')
CALL prms%CreateLogicalOption(  'Particles-BGK-DoBGKCellSplitting',       'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-BGK-SampAdapFac',       'TODO-DEFINE-PARAMETER','.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-BGK-DoVibRelaxation',       'TODO-DEFINE-PARAMETER','.FALSE.')

END SUBROUTINE DefineParametersBGK

SUBROUTINE InitESBGK()
!===================================================================================================================================
!> Init of BGK Vars
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Vars      ,ONLY: nSpecies, Species, BoltzmannConst
USE MOD_DSMC_Vars          ,ONLY: SpecDSMC
USE MOD_Globals_Vars       ,ONLY: Pi
USE MOD_ReadInTools
USE MOD_ESBGK_Vars
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation_Vars ,ONLY: xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf
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

BGKAdaptTimeStep = GETINT('Particles-BGKAdaptTimeStep','0')
DoBGKCellAdaptation = GETLOGICAL('Particles-BGK-DoCellAdaptation','.FALSE.')
BGKMinPartPerCell = GETINT('Particles-BGK-MinPartsPerCell','10')
BGKCollModel = GETINT('Particles-BGKCollModel','1')
ESBGKModel = GETINT('Particles-ESBGKModel','1')         ! 1: Approximative, 2: Exact, 3: MetropolisHastings
BGKUnifiedCes = GETREAL('Particles-UnifiedBGK-Ces','1000.')
IF (BGKUnifiedCes.EQ.1000.) THEN
  BGKUnifiedCes = 1. - (6.-2.*SpecDSMC(1)%omegaVHS)*(4.- 2.*SpecDSMC(1)%omegaVHS)/30.
END IF
BGKAveragingLength = GETINT('Particles-BGKAveragingLength','5')
BGKDoAveraging = GETLOGICAL('Particles-BGKDoAveraging','.FALSE.')
BGKDoAveragingCorrect = GETLOGICAL('Particles-BGKDoAveragingCorrection','.FALSE.')
BGKUseQuantVibEn = GETLOGICAL('Particles-BGKUseQuantVibEn','.FALSE.')
IF (BGKDoAveraging) CALL BGK_init_Averaging()
BGKAcceleration = GETREAL('Particles-BGKAcceleration','-9.81')
BGKDoVibRelaxation = GETLOGICAL('Particles-BGK-DoVibRelaxation','.TRUE.')
DoBGKCellSplitting = GETLOGICAL('Particles-BGK-DoBGKCellSplitting','.FALSE.')
BGKSampAdapFac = GETLOGICAL('Particles-BGK-SampAdapFac','.FALSE.')
IF (DoBGKCellSplitting.OR.BGKSampAdapFac) ALLOCATE(ElemSplitCells(nElems))
IF (BGKSampAdapFac) THEN
  DO iElem = 1, nElems
    ALLOCATE(ElemSplitCells(iElem)%AdapFac(4))
    ElemSplitCells(iElem)%AdapFac = 0.0
  END DO
END IF
IF (DoBGKCellSplitting) THEN
  DoBGKCellAdaptation = .FALSE.
  CALL DefineElementOrientation()
  DO iElem = 1, nElems
    ElemSplitCells(iElem)%Splitnum(ElemSplitCells(iElem)%CellOrientation(1)) = 2
    ElemSplitCells(iElem)%Splitnum(ElemSplitCells(iElem)%CellOrientation(2)) = 0
    ElemSplitCells(iElem)%Splitnum(ElemSplitCells(iElem)%CellOrientation(3)) = 0
  END DO
END IF

SWRITE(UNIT_stdOut,'(A)') ' INIT BGK DONE!'

END SUBROUTINE InitESBGK


SUBROUTINE ESBGK_BuildTransGaussNumsEnCon(nPart, iRanPart)
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

END SUBROUTINE ESBGK_BuildTransGaussNumsEnCon



SUBROUTINE DefineElementOrientation()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems,NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Globals_Vars       ,ONLY: Pi
USE MOD_ESBGK_Vars         ,ONLY: ElemSplitCells
USE MOD_Eval_xyz,                    ONLY:Eval_XYZ_Poly
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL      :: Vec(3), P(3,8), ResVec1(3,3), ResVec2(3,3), DiffVec(3,3), MinVec(3), MaxDir(3), angle(3,3), vecmag(3)
INTEGER   :: iElem, iNode, ii, jj, firstDir1, firstDir2
!===================================================================================================================================
DO iElem = 1, nElems
  Vec(1:3) = (/-0.99,0.,0./)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,1))
  Vec(1:3) = (/0.99,0.,0./)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,1))
  Vec(1:3) = (/0.,-0.99,0./)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,2))
  Vec(1:3) = (/0.,0.99,0./)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,2))
  Vec(1:3) = (/0.,0.,-0.99/)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec1(1:3,3))
  Vec(1:3) = (/0.,0.,0.99/)
  CALL Eval_xyz_Poly(Vec(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),ResVec2(1:3,3))

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
INTEGER                       :: iElem, NodeDepth
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
