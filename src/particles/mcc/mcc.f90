!==================================================================================================================================
! Copyright (c) 2021 boltzplatz - numerical plasma dynamics GmbH
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

MODULE MOD_MCC
!===================================================================================================================================
! Module for MCC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: MCC
!===================================================================================================================================

CONTAINS

SUBROUTINE MCC(iElem)
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze            ,ONLY: CalcGammaVib, CalcMeanFreePath, DSMCMacroSampling, SummarizeQualityFactors
USE MOD_part_emission_tools     ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SpecXSec, DSMC_RHS
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, MCC_TotalPairNum, DSMCSumOfFormedParticles, XSec_NullCollision
USE MOD_DSMC_Vars               ,ONLY: PolyatomMolDSMC, VibQuantsPar, UseBGGasSplit, RadialWeighting
USE MOD_Part_Emission_Tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_Part_Pos_and_Velo       ,ONLY: SetParticleVelocity
USE MOD_Particle_Vars           ,ONLY: PEM, PDM, PartSpecies, nSpecies, PartState, Species, usevMPF, PartMPF, Species, PartPosRef
USE MOD_Particle_Vars           ,ONLY: VarTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Mesh_Vars               ,ONLY: nElems, offSetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars           ,ONLY: WriteMacroVolumeValues, WriteMacroSurfaceValues
USE MOD_DSMC_Collis             ,ONLY: DSMC_perform_collision
USE MOD_DSMC_Relaxation         ,ONLY: FinalizeCalcVibRelaxProb, SumVibRelaxProb, InitCalcVibRelaxProb
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, time
USE MOD_DSMC_CollisionProb      ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Relaxation         ,ONLY: CalcMeanVibQuaDiatomic
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_InsertParticles, AD_DeleteParticles
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_TimeDisc_Vars           ,ONLY: dt
USE MOD_DSMC_SpecXSec           ,ONLY: InterpolateCrossSection
USE MOD_DSMC_SpecXSec           ,ONLY: XSec_CalcReactionProb, XSec_CalcVibRelaxProb
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! LOCAL VARIABLES
INTEGER                       :: iPair, iPart, iLoop, nPart, iSpec, jSpec, bgSpec, PartIndex, bggPartIndex, PairCount, RandomPart
INTEGER                       :: cSpec1, cSpec2, iCase, SpecPairNumTemp, nPartAmbi, OldPairNum, CNElemID
INTEGER                       :: iPart_p1, iPart_p2, iPairNew, iPartSplit, SplitPartNum, SplitRestPart
INTEGER,ALLOCATABLE           :: iPartIndexSpec(:,:), SpecPartNum(:), SpecPairNum(:)
REAL                          :: iRan, ProbRest, SpecPairNumReal, MPF, Volume
INTEGER, ALLOCATABLE          :: iPartIndx_NodeTotalAmbiDel(:)
INTEGER, ALLOCATABLE, TARGET  :: iPartIndx_Node(:), iPartIndx_NodeTotalAmbi(:)
INTEGER, POINTER              :: iPartIndx_NodeTotal(:)
LOGICAL                       :: SplitInProgress
REAL                          :: CollProb, VeloBGGPart(1:3), CRela2, CollEnergy
REAL                          :: PartStateSplit(1:6), PartPosRefSplit(1:3), PartStateIntSplit(1:3), PartTimeStepSplit, PartMPFSplit
INTEGER, ALLOCATABLE          :: VibQuantsParSplit(:)
!===================================================================================================================================

CNElemID = GetCNElemID(iElem+offSetElem)
Volume = ElemVolume_Shared(CNElemID)
! Create particle index list for pairing
nPart = PEM%pNumber(iElem)
ALLOCATE(iPartIndx_Node(nPart))
iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  iPartIndx_Node(iLoop) = iPart
  iPart = PEM%pNext(iPart)
END DO

! Ambipolar Diffusion
IF (DSMC%DoAmbipolarDiff) THEN
  CALL AD_InsertParticles(iPartIndx_Node,nPart, iPartIndx_NodeTotalAmbi, nPartAmbi)
  ALLOCATE(iPartIndx_NodeTotalAmbiDel(1:nPartAmbi))
  iPartIndx_NodeTotalAmbiDel(1:nPartAmbi) = iPartIndx_NodeTotalAmbi(1:nPartAmbi)
  nPart = nPartAmbi
  iPartIndx_NodeTotal => iPartIndx_NodeTotalAmbi
ELSE
  iPartIndx_NodeTotal => iPartIndx_Node
END IF

MCC_TotalPairNum = 0

CollInf%Coll_SpecPartNum = 0.
CollInf%Coll_CaseNum = 0
CollInf%MeanMPF = 0.

ALLOCATE(iPartIndexSpec(nPart,nSpecies))
iPartIndexSpec = 0

ALLOCATE(SpecPartNum(nSpecies),SpecPairNum(CollInf%NumCase))
SpecPairNum = 0; SpecPairNumTemp = 0; SpecPairNumReal = 0.; SpecPartNum = 0
CALL InitCalcVibRelaxProb()

IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

! X.) Counting the number of particles per species and creating a species-specific particle index list
DO iLoop = 1, nPart
  iPart = iPartIndx_NodeTotal(iLoop)
  iSpec = PartSpecies(iPart)
  MPF = GetParticleWeight(iPart)
  SpecPartNum(iSpec) = SpecPartNum(iSpec) + 1
  ! Sum of the particle weights (in case the particle is split later, the sum of the weights remains constant and is equal to this greater weight added here)
  CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + MPF
  ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(1,iPart) * MPF
  ! Create species-specific particle index list for cross-section based pairing
  iPartIndexSpec(SpecPartNum(iSpec),iSpec) = iPart
END DO

! X.) Determine the particle number of the background species and calculate the cell temperature
DO bgSpec = 1, BGGas%NumberOfSpecies
  iSpec = BGGas%MapBGSpecToSpec(bgSpec)
  IF(usevMPF) THEN
    CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bgSpec)*Volume
  ELSE
    CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bgSpec)*Volume/Species(iSpec)%MacroParticleFactor
  END IF
END DO

IF(DSMC%CalcQualityFactors) THEN
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = 0.
  DO bgSpec = 1, BGGas%NumberOfSpecies
    iSpec = BGGas%MapBGSpecToSpec(bgSpec)
    DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) + BGGas%SpeciesFraction(bgSpec) &
                                                                            * Species(iSpec)%Init(1)%MWTemperatureIC
  END DO
END IF

! 2.) Determining the total number of pairs
DO iSpec = 1,nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) CYCLE    ! Loop over all non-background species
  DO bgSpec = 1, BGGas%NumberOfSpecies        ! Loop over all background species
    jSpec = BGGas%MapBGSpecToSpec(bgSpec)
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    IF(SpecXSec(iCase)%UseCollXSec.AND.XSec_NullCollision) THEN
      ! Collision cross-section: The maximum number of pairs to check is collision pair specific and depends on the null collision probability
      SpecPairNumReal = SpecPartNum(iSpec)*SpecXSec(iCase)%ProbNull
      SpecPairNumTemp = INT(SpecPartNum(iSpec)*SpecXSec(iCase)%ProbNull)
    ELSE
      ! Regular: The maximum number of pairs corresponds to the particle number
      SpecPairNumReal = BGGas%SpeciesFraction(bgSpec)*SpecPartNum(iSpec)
      SpecPairNumTemp = INT(BGGas%SpeciesFraction(bgSpec)*SpecPartNum(iSpec))
    END IF
    ! Avoid creating more pairs than currently particles in the simulation
    IF(SpecPairNum(iCase) + SpecPairNumTemp.LT.SpecPartNum(iSpec)) THEN
      ! Randomly deciding whether an additional pair is added based on the difference between the real and integer value
      ProbRest = SpecPairNumReal - REAL(SpecPairNumTemp)
      CALL RANDOM_NUMBER(iRan)
      IF (ProbRest.GT.iRan) SpecPairNumTemp = SpecPairNumTemp + 1
      ! Adding the number of pairs to the species-specific number and the cell total
      SpecPairNum(iCase) = SpecPairNum(iCase) + SpecPairNumTemp
      MCC_TotalPairNum = MCC_TotalPairNum + SpecPairNumTemp
    ELSE IF(SpecPairNum(iCase) + SpecPairNumTemp.EQ.SpecPartNum(iSpec)) THEN
      SpecPairNum(iCase) = SpecPairNum(iCase) + SpecPairNumTemp
      MCC_TotalPairNum = MCC_TotalPairNum + SpecPairNumTemp
    END IF
  END DO
END DO

ALLOCATE(Coll_pData(1))
Coll_pData%Ec = 0.
bggPartIndex = 0
SplitInProgress = .FALSE.
iPartSplit = 0
SplitPartNum = 0
SplitRestPart = 0
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) CYCLE    ! Loop over all non-background species
  DO bgSpec = 1, BGGas%NumberOfSpecies        ! Loop over all background species
    jSpec = BGGas%MapBGSpecToSpec(bgSpec)
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    iLoop = 1
    ! DO iLoop = 1, SpecPairNum(iCase)    ! Loop over all the number of pairs required for this species pairing
    DO WHILE(iLoop.LE.SpecPairNum(iCase))
      ! Choosing random particles from the available number of particles, getting the index of the simulation particle
      IF(.NOT.SplitInProgress) THEN
        CALL RANDOM_NUMBER(iRan)
        RandomPart = INT(SpecPartNum(iSpec)*iRan) + 1
        PartIndex = iPartIndexSpec(RandomPart,iSpec)
        iPartIndexSpec(RandomPart, iSpec) = iPartIndexSpec(SpecPartNum(iSpec),iSpec)
        SpecPartNum(iSpec) = SpecPartNum(iSpec) - 1
      END IF
      ! ==============================================================================================================================
      ! BGGasSplit
      ! ==============================================================================================================================
      IF(usevMPF.AND.UseBGGasSplit) THEN
        IF(SplitInProgress) THEN
          IF(iPartSplit.LT.SplitPartNum) THEN
            ! Clone the regular particle (re-using the index of the previous particle if it didn't collide)
            IF(PartIndex.EQ.0) THEN
              DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
              PartIndex = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
              IF (PartIndex.EQ.0) THEN
                CALL Abort(__STAMP__,'ERROR in MCC: MaxParticleNumber should be increased!')
              END IF
              PartSpecies(PartIndex) = iSpec
              PartState(1:6,PartIndex) = PartStateSplit(1:6)
              IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
                PartPosRef(1:3,PartIndex)=PartPosRefSplit(1:3)
              END IF
              IF(CollisMode.GT.1) THEN
                PartStateIntEn(1:2,PartIndex) = PartStateIntSplit(1:2)
                IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                  IF(ALLOCATED(VibQuantsPar(PartIndex)%Quants)) DEALLOCATE(VibQuantsPar(PartIndex)%Quants)
                  ALLOCATE(VibQuantsPar(PartIndex)%Quants(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%VibDOF))
                  VibQuantsPar(PartIndex)%Quants(:) = VibQuantsParSplit(:)
                END IF
                IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PartIndex) = PartStateIntSplit(3)
              END IF
              PEM%GlobalElemID(PartIndex)     = iElem + offSetElem
              PDM%ParticleInside(PartIndex)   = .TRUE.
              PEM%LastGlobalElemID(PartIndex) = PEM%GlobalElemID(PartIndex)
              PDM%IsNewPart(PartIndex)        = .TRUE.
              PDM%dtFracPush(PartIndex)       = .FALSE.
              PartMPF(PartIndex)              = PartMPFSplit
              IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(PartIndex) = PartTimeStepSplit
            END IF
            iPartSplit = iPartSplit + 1
          END IF
        ELSE IF(PartMPF(PartIndex).GT.Species(jSpec)%MacroParticleFactor) THEN
          ! Variante 1
          ! SplitPartNum = NINT(PartMPF(PartIndex) / Species(jSpec)%MacroParticleFactor) - 1
          ! Variante 2
          CALL RANDOM_NUMBER(iRan)
          SplitPartNum = INT(PartMPF(PartIndex) / Species(jSpec)%MacroParticleFactor+iRan) - 1
          IF(SplitPartNum.GT.0) THEN
            SplitInProgress  = .TRUE.
            SpecPairNum(iCase)      = SpecPairNum(iCase) + SplitPartNum
            PartStateSplit(1:6) = PartState(1:6,PartIndex)
            IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
              PartPosRefSplit(1:3) = PartPosRef(1:3,PartIndex)
            END IF
            IF(CollisMode.GT.1) THEN
              PartStateIntSplit(1:2) = PartStateIntEn(1:2,PartIndex)
              IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                ALLOCATE(VibQuantsParSplit(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%VibDOF))
                VibQuantsParSplit(:) = VibQuantsPar(PartIndex)%Quants(:)
              END IF
              IF(DSMC%ElectronicModel.GT.0) PartStateIntSplit(3) = PartStateIntEn(3,PartIndex)
            END IF
            IF(VarTimeStep%UseVariableTimeStep) PartTimeStepSplit = VarTimeStep%ParticleTimeStep(PartIndex)
            ! PartMPFSplit          = Species(jSpec)%MacroParticleFactor
            ! Set the new MPF based on the actual number of split particles
            PartMPFSplit          = PartMPF(PartIndex) / REAL(SplitPartNum+1)
            PartMPF(PartIndex)    = PartMPFSplit
          END IF
        END IF
      END IF
      ! ==============================================================================================================================
      ! Determine collision probability
      ! ==============================================================================================================================
      CALL CalcVelocity_maxwell_lpn(FractNbr=jSpec, Vec3D=VeloBGGPart(1:3), iInit=1)
      CRela2 = (PartState(4,PartIndex) - VeloBGGPart(1))**2 &
                          + (PartState(5,PartIndex) - VeloBGGPart(2))**2 &
                          + (PartState(6,PartIndex) - VeloBGGPart(3))**2
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! Using the relative kinetic energy of the particle pair (real energy value per particle pair, no weighting/scaling factors)
        CollEnergy = 0.5 * CollInf%MassRed(iCase) * CRela2
        ! Calculate the collision probability
        SpecXSec(iCase)%CrossSection = InterpolateCrossSection(iCase,CollEnergy)
        CollProb = (1. - EXP(-SQRT(CRela2) * SpecXSec(iCase)%CrossSection * BGGas%NumberDensity(bgSpec) * dt))
        ! Correct the collision probability in the case of the second species being a background species as the number of pairs
        ! is either determined based on the null collision probability or on the species fraction
        IF(XSec_NullCollision) THEN
          CollProb = CollProb / SpecXSec(iCase)%ProbNull
        ELSE
          CollProb = CollProb / BGGas%SpeciesFraction(bgSpec)
        END IF
      ELSE
        CollProb = 0.
      END IF
      ! ==============================================================================================================================
      ! Check whether a collision occurs
      ! ==============================================================================================================================
      CALL RANDOM_NUMBER(iRan)
      IF (CollProb.GE.iRan) THEN
        iPair = 1
        Coll_pData(iPair)%iPart_p1 = PartIndex
        ! Creating a new background gas particle
        IF(bggPartIndex.EQ.0) THEN
          DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
          bggPartIndex = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
          IF (bggPartIndex.EQ.0) THEN
            CALL Abort(__STAMP__,'ERROR in MCC: MaxParticleNumber should be increased!')
          END IF
        END IF
        ! Position the background particle at the simulation particle
        PartState(1:3,bggPartIndex) = PartState(1:3,PartIndex)
        PartState(4:6,bggPartIndex) = VeloBGGPart(1:3)
        IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
          PartPosRef(1:3,bggPartIndex)=PartPosRef(1:3,PartIndex)
        END IF
        ! Set the species of the background gas particle
        PartSpecies(bggPartIndex) = jSpec
        IF(CollisMode.GT.1) THEN
          IF(SpecDSMC(jSpec)%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(jSpec,1,bggPartIndex,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(jSpec,1,bggPartIndex,1)
          END IF
        END IF
        PEM%GlobalElemID(bggPartIndex) = iElem + offSetElem
        PDM%ParticleInside(bggPartIndex) = .TRUE.
        PEM%LastGlobalElemID(bggPartIndex) = PEM%GlobalElemID(PartIndex)
        PDM%IsNewPart(bggPartIndex)       = .TRUE.
        PDM%dtFracPush(bggPartIndex)      = .FALSE.
        ! Ambipolar diffusion: add the background particle to consider (as its index might be used for an ion after a reaction)
        IF(DSMC%DoAmbipolarDiff) THEN
          newAmbiParts = newAmbiParts + 1
          iPartIndx_NodeNewAmbi(newAmbiParts) = bggPartIndex
        END IF
        IF(usevMPF) PartMPF(bggPartIndex) = PartMPF(PartIndex)
        IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(bggPartIndex) = VarTimeStep%ParticleTimeStep(PartIndex)
        Coll_pData(iPair)%iPart_p2 = bggPartIndex
        ! Required in DSMC_Relax_Col_LauxTSHO
        Coll_pData(iPair)%cRela2 = CRela2
        ! Required in 
        Coll_pData(iPair)%Prob = CollProb
        ! Required in XSec_CalcReactionProb
        Coll_pData(iPair)%PairType = iCase
        ! Perform collision
        CALL DSMC_perform_collision(iPair,iElem)
        ! If the species of the particle has changed from the background species (e.g. during chemistry) then get a new index for the next particle
        IF(PartSpecies(bggPartIndex).NE.jSpec) THEN
          ! Add new particle to the index list
          PEM%pNext(PEM%pEnd(iElem)) = bggPartIndex
          PEM%pEnd(iElem) = bggPartIndex
          PEM%pNumber(iElem) = PEM%pNumber(iElem) + 1
          ! Set index to zero to get a new one for the next background gas particle
          bggPartIndex = 0
        END IF
        IF(usevMPF.AND.UseBGGasSplit) THEN
          ! Set index to zero to get a new one for the next split particle
          PartIndex = 0
        END IF
      ELSE  ! No collision
        IF(SplitInProgress) THEN
          ! Save the index of the first particle that did not collide
          IF(SplitRestPart.EQ.0) THEN 
            SplitRestPart = PartIndex
            ! Reset the PartIndex to use a new particle (unless it is the last particle, keep the index to check whether it can be deleted)
            IF(iPartSplit.NE.SplitPartNum) PartIndex = 0
          ELSE
            ! Add the particle to the others that did not collide (-> merging those particles)
            PartMPF(SplitRestPart) = PartMPF(SplitRestPart) + PartMPF(PartIndex)
          END IF
          ! Delete the last particle that did not collide (unless it is the first and last)
          IF(iPartSplit.EQ.SplitPartNum) THEN
            IF(SplitRestPart.NE.PartIndex) PDM%ParticleInside(PartIndex) = .FALSE.
          END IF
        END IF
      END IF  ! CollProb.GE.iRan
      IF(SplitInProgress) THEN
        ! Treatment at the end of the split
        IF(iPartSplit.EQ.SplitPartNum) THEN
          ! Reset the split counters
          iPartSplit = 0
          SplitPartNum = 0
          SplitRestPart = 0
          SplitInProgress = .FALSE.
          SDEALLOCATE(VibQuantsParSplit)
        END IF
      END IF
      iLoop = iLoop + 1
    END DO    ! DO WHILE(iLoop.LE.SpecPairNum(iCase))
  END DO      ! bgSpec = 1, BGGas%NumberOfSpecies
END DO        ! iSpec = 1, nSpecies
! Delete the dummy particle
IF(bggPartIndex.NE.0) THEN
  PDM%ParticleInside(bggPartIndex) = .FALSE.
END IF
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum), &
                          ElemVolume_Shared(GetCNElemID(iElem+offSetElem)), DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the MCS/MFP for the case without octree
    IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                    / DSMC%MeanFreePath
  END IF
  CALL SummarizeQualityFactors(iElem)
END IF
DEALLOCATE(iPartIndx_Node)
SDEALLOCATE(iPartIndx_NodeTotalAmbiDel)
DEALLOCATE(iPartIndexSpec)
DEALLOCATE(SpecPartNum)
DEALLOCATE(SpecPairNum)
DEALLOCATE(Coll_pData)

END SUBROUTINE MCC

END MODULE MOD_MCC
