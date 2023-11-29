!==================================================================================================================================
! Copyright (c) 2021 boltzplatz - numerical plasma dynamics GmbH
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

MODULE MOD_MCC
!===================================================================================================================================
! Module for Monte Carlo Collisions
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
PUBLIC :: MonteCarloCollision
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Monte Carlo Collision routine
!> 1. Counting the number of particles per species and creating a species-specific particle index list
!> 2. Determine the particle number of the background species and calculate the cell temperature
!> 3. Determining the total number of pairs of test
!> 4. Loop over the total number of pairs per species, determine collision probability and perform collision
!>    a. Draw random particle indices from the particle list
!>    b. (Trace species background) Clone particle to test against trace background gas species
!>    c. Get particle index (or re-use the previous) and get a velocity from the background gas
!>    d. Calculate the collision probability
!>    e. Perform collision: set particle information of the background gas particle, if the particle is not utilized in a reaction
!>       and did not changed its species, utilize it in the next pairing
!===================================================================================================================================
SUBROUTINE MonteCarloCollision(iElem)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
! VARIABLES
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, DSMCSumOfFormedParticles, PolyatomMolDSMC, VibQuantsPar
USE MOD_MCC_Vars                ,ONLY: SpecXSec, XSec_NullCollision
USE MOD_Particle_Vars           ,ONLY: PEM, PDM, PartSpecies, nSpecies, PartState, Species, usevMPF, PartMPF, Species, PartPosRef
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars           ,ONLY: WriteMacroVolumeValues
USE MOD_TimeDisc_Vars           ,ONLY: dt, TEnd, time
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
! ROUTINES
USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
USE MOD_DSMC_BGGas              ,ONLY: BGGas_AssignParticleProperties
USE MOD_part_tools              ,ONLY: GetParticleWeight, CalcVelocity_maxwell_particle, GetNextFreePosition
USE MOD_Part_Emission_Tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_Collis             ,ONLY: DSMC_perform_collision
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_InsertParticles, AD_DeleteParticles
USE MOD_MCC_XSec                ,ONLY: InterpolateCrossSection
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
INTEGER                       :: iPair, iPart, iLoop, nPart, iSpec, jSpec, bgSpec, PartIndex, bggPartIndex, RandomPart
INTEGER                       :: iCase, SpecPairNumTemp, nPartAmbi, CNElemID, GlobalElemID
INTEGER                       :: iLevel, nVib, iPartSplit, SplitPartNum, SplitRestPart
INTEGER,ALLOCATABLE           :: iPartIndexSpec(:,:), SpecPartNum(:), SpecPairNum(:), UseSpecPartNum(:)
REAL                          :: iRan, ProbRest, SpecPairNumReal, MPF, Volume, MPFRatio, BGGasNumDens, BGGasFraction
INTEGER, ALLOCATABLE          :: iPartIndx_NodeTotalAmbiDel(:)
INTEGER, ALLOCATABLE, TARGET  :: iPartIndx_Node(:), iPartIndx_NodeTotalAmbi(:)
INTEGER, POINTER              :: iPartIndx_NodeTotal(:)
LOGICAL                       :: SplitInProgress, GetInternalEnergy
REAL                          :: CollCaseNum, CollProb, VeloBGGPart(1:3), CRela2, CollEnergy, GammaFac, SumVibCrossSection
REAL                          :: PartStateSplit(1:6), PartPosRefSplit(1:3), PartStateIntSplit(1:3), PartTimeStepSplit, PartMPFSplit
INTEGER, ALLOCATABLE          :: VibQuantsParSplit(:), PartIndexCase(:)
REAL                          :: ProbNull, dtVar
!===================================================================================================================================

! Skip elements outside of any background gas regions
IF(BGGas%UseRegions) THEN
  IF(BGGas%RegionElemType(iElem).EQ.0) RETURN
END IF

GlobalElemID = iElem+offSetElem
CNElemID = GetCNElemID(GlobalElemID)
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

CollInf%Coll_SpecPartNum = 0.
CollInf%Coll_CaseNum = 0
CollInf%SumPairMPF = 0.

ALLOCATE(iPartIndexSpec(nPart,nSpecies))
iPartIndexSpec = 0

ALLOCATE(SpecPartNum(nSpecies),SpecPairNum(CollInf%NumCase))
ALLOCATE(UseSpecPartNum(nSpecies))
SpecPairNum = 0
SpecPairNumTemp = 0
SpecPairNumReal = 0.
SpecPartNum = 0
UseSpecPartNum = 0

IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

! 1.) Counting the number of particles per species and creating a species-specific particle index list
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

! 2.) Determine the particle number of the background species and calculate the cell temperature
DO bgSpec = 1, BGGas%NumberOfSpecies
  iSpec = BGGas%MapBGSpecToSpec(bgSpec)
  IF(BGGas%UseDistribution) THEN
    BGGasNumDens = BGGas%Distribution(bgSpec,7,iElem)
  ELSE
    BGGasNumDens = BGGas%NumberDensity(bgSpec)
  END IF
  IF(usevMPF) THEN
    CollInf%Coll_SpecPartNum(iSpec) = BGGasNumDens*Volume
  ELSE
    CollInf%Coll_SpecPartNum(iSpec) = BGGasNumDens*Volume/Species(iSpec)%MacroParticleFactor
  END IF
END DO

IF(DSMC%CalcQualityFactors) THEN
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = 0.
  DO bgSpec = 1, BGGas%NumberOfSpecies
    iSpec = BGGas%MapBGSpecToSpec(bgSpec)
    IF(BGGas%UseDistribution) THEN
      DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) + BGGas%SpeciesFractionElem(bgSpec,iElem) &
                                                                              * SUM(BGGas%Distribution(bgSpec,4:6,iElem)) / 3.
    ELSE
      DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) + BGGas%SpeciesFraction(bgSpec) &
                                                                              * Species(iSpec)%Init(1)%MWTemperatureIC
    END IF
  END DO
END IF

! 3.) Determining the total number of pairs
DO iSpec = 1,nSpecies
  IF(SpecPartNum(iSpec).EQ.0) CYCLE ! No particles for species iSpec present
  IF(BGGas%BackgroundSpecies(iSpec)) CYCLE    ! Loop over all non-background species
  DO bgSpec = 1, BGGas%NumberOfSpecies        ! Loop over all background species
    jSpec = BGGas%MapBGSpecToSpec(bgSpec)
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    IF(SpecXSec(iCase)%UseCollXSec.AND.XSec_NullCollision) THEN
      ! Collision cross-section: The maximum number of pairs to check is collision pair specific and depends on the null collision probability
      IF(BGGas%UseDistribution) THEN
        SpecPairNumReal = SpecPartNum(iSpec)*SpecXSec(iCase)%ProbNullElem(iElem)
      ELSE
        SpecPairNumReal = SpecPartNum(iSpec)*SpecXSec(iCase)%ProbNull
      END IF
    ELSE
      ! Regular: The maximum number of pairs corresponds to the particle number
      IF(BGGas%UseDistribution)THEN
        SpecPairNumReal = BGGas%SpeciesFractionElem(bgSpec,iElem)*SpecPartNum(iSpec)
      ELSE
        SpecPairNumReal = BGGas%SpeciesFraction(bgSpec)*SpecPartNum(iSpec)
      END IF ! BGGas%UseDistribution
    END IF
    SpecPairNumTemp = INT(SpecPairNumReal)
    ! Avoid creating more pairs than currently particles in the simulation
    IF(UseSpecPartNum(iSpec) + SpecPairNumTemp.LT.SpecPartNum(iSpec)) THEN
      ! Randomly deciding whether an additional pair is added based on the difference between the real and integer value
      ProbRest = SpecPairNumReal - REAL(SpecPairNumTemp)
      CALL RANDOM_NUMBER(iRan)
      IF (ProbRest.GT.iRan) SpecPairNumTemp = SpecPairNumTemp + 1
      ! Adding the number of pairs to the case-specific number of pairs
      SpecPairNum(iCase) = SpecPairNum(iCase) + SpecPairNumTemp
      ! Adding the number of pairs to the number of particles utilized per species to avoid defining more pairs than available partners
      UseSpecPartNum(iSpec) = UseSpecPartNum(iSpec) + SpecPairNumTemp
    ELSE IF(UseSpecPartNum(iSpec) + SpecPairNumTemp.EQ.SpecPartNum(iSpec)) THEN
      SpecPairNum(iCase) = SpecPairNum(iCase) + SpecPairNumTemp
      UseSpecPartNum(iSpec) = UseSpecPartNum(iSpec) + SpecPairNumTemp
    END IF
  END DO
END DO

! 4.) Loop over the total number of pairs per species, determine collision probability and perform collision
ALLOCATE(Coll_pData(1))
Coll_pData%Ec = 0.
bggPartIndex = 0
SplitInProgress = .FALSE.
GetInternalEnergy = .TRUE.
iPartSplit = 0
SplitPartNum = 0
SplitRestPart = 0

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) CYCLE    ! Loop over all non-background species
  DO bgSpec = 1, BGGas%NumberOfSpecies        ! Loop over all background species
    jSpec = BGGas%MapBGSpecToSpec(bgSpec)
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    MPFRatio = BGGas%MaxMPF/Species(jSpec)%MacroParticleFactor
    IF(SpecPairNum(iCase).EQ.0) CYCLE
    ! Randomly determining the particle index and storing it, determine the weighted number of pairs per case
    ! (only required for VHS collision modelling in combination with MCC/XSec and vMPF)
    CollCaseNum = 0.
    ALLOCATE(PartIndexCase(SpecPairNum(iCase)))
    DO iLoop = 1, SpecPairNum(iCase)
      CALL RANDOM_NUMBER(iRan)
      RandomPart = INT(SpecPartNum(iSpec)*iRan) + 1
      PartIndexCase(iLoop) = iPartIndexSpec(RandomPart,iSpec)
      CollCaseNum = CollCaseNum + GetParticleWeight(PartIndexCase(iLoop))
      iPartIndexSpec(RandomPart, iSpec) = iPartIndexSpec(SpecPartNum(iSpec),iSpec)
      SpecPartNum(iSpec) = SpecPartNum(iSpec) - 1
    END DO
    ! Loop over all the number of pairs required for this species pairing
    iLoop = 1
    iPart = 1
    DO WHILE(iLoop.LE.SpecPairNum(iCase))
      ! Getting the index of the simulation particle (previously randomly determined)
      IF(.NOT.SplitInProgress) THEN
        PartIndex = PartIndexCase(iPart)
        iPart = iPart + 1
      END IF
      ! ============================================================================================================================
      ! Treatment of background gas trace species:
      ! 1) Determining whether a collision with a trace species is imminent and enable split (clone particle to increase sample size)
      ! 2) If a split is in progress, get new index particle index and clone the simulation particle to test with the trace species
      ! ============================================================================================================================
      IF(usevMPF.AND.BGGas%TraceSpecies(jSpec)) THEN
        IF(SplitInProgress) THEN
          IF(iPartSplit.LT.SplitPartNum) THEN
            ! Clone the regular particle (re-using the index of the previous particle if it didn't collide)
            IF(PartIndex.EQ.0) THEN
              DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
              PartIndex = GetNextFreePosition()
              IF (PartIndex.EQ.0) THEN
                CALL Abort(__STAMP__,'ERROR in MCC: MaxParticleNumber should be increased!')
              END IF
              PartSpecies(PartIndex) = iSpec
              PartState(1:6,PartIndex) = PartStateSplit(1:6)
              IF(TrackingMethod.EQ.REFMAPPING) PartPosRef(1:3,PartIndex)=PartPosRefSplit(1:3)
              IF(CollisMode.GT.1) THEN
                PartStateIntEn(1:2,PartIndex) = PartStateIntSplit(1:2)
                IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                  IF(ALLOCATED(VibQuantsPar(PartIndex)%Quants)) DEALLOCATE(VibQuantsPar(PartIndex)%Quants)
                  ALLOCATE(VibQuantsPar(PartIndex)%Quants(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%VibDOF))
                  VibQuantsPar(PartIndex)%Quants(:) = VibQuantsParSplit(:)
                END IF
                IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PartIndex) = PartStateIntSplit(3)
              END IF
              ! Set global element indices
              PEM%GlobalElemID(PartIndex)     = GlobalElemID
              PEM%LastGlobalElemID(PartIndex) = PEM%GlobalElemID(PartIndex)
              ! Set simulation flags
              PDM%ParticleInside(PartIndex)   = .TRUE.
              PDM%IsNewPart(PartIndex)        = .TRUE.
              PDM%dtFracPush(PartIndex)       = .FALSE.
              ! Set particle weights
              PartMPF(PartIndex)              = PartMPFSplit
              IF(UseVarTimeStep) PartTimeStep(PartIndex) = PartTimeStepSplit
              ! Add new particle to the index list
              PEM%pNext(PEM%pEnd(iElem)) = PartIndex
              PEM%pEnd(iElem) = PartIndex
              PEM%pNumber(iElem) = PEM%pNumber(iElem) + 1
              ! Consider ambipolar diffusion
              IF (DSMC%DoAmbipolarDiff) THEN
                newAmbiParts = newAmbiParts + 1
                iPartIndx_NodeNewAmbi(newAmbiParts) = PartIndex
              END IF
            END IF
            iPartSplit = iPartSplit + 1
          END IF
        ! Enable split if the weighting of the particle species is greater than the trace species weighting factor
        ELSE IF(PartMPF(PartIndex).GT.Species(jSpec)%MacroParticleFactor) THEN
          ! Determine the number of additional split particles, avoid splitting the particle more often than the ratio between the
          ! largest background gas weighting factor and the trace species weighting factor
          CALL RANDOM_NUMBER(iRan)
          IF(PartMPF(PartIndex) / Species(jSpec)%MacroParticleFactor.GT.MPFRatio) THEN
            SplitPartNum = INT(MPFRatio+iRan) - 1
          ELSE
            SplitPartNum = INT(PartMPF(PartIndex) / Species(jSpec)%MacroParticleFactor+iRan) - 1
          END IF
          IF(SplitPartNum.GT.0) THEN
            SplitInProgress  = .TRUE.
            SpecPairNum(iCase)  = SpecPairNum(iCase) + SplitPartNum
            PartStateSplit(1:6) = PartState(1:6,PartIndex)
            IF(TrackingMethod.EQ.REFMAPPING) PartPosRefSplit(1:3) = PartPosRef(1:3,PartIndex)
            IF(CollisMode.GT.1) THEN
              PartStateIntSplit(1:2) = PartStateIntEn(1:2,PartIndex)
              IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                ALLOCATE(VibQuantsParSplit(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%VibDOF))
                VibQuantsParSplit(:) = VibQuantsPar(PartIndex)%Quants(:)
              END IF ! SpecDSMC(iSpec)%PolyatomicMol
              IF(DSMC%ElectronicModel.GT.0) PartStateIntSplit(3) = PartStateIntEn(3,PartIndex)
            END IF ! CollisMode.GT.1
            IF(UseVarTimeStep) PartTimeStepSplit = PartTimeStep(PartIndex)
            ! Set the new MPF based on the actual number of split particles
            PartMPFSplit          = PartMPF(PartIndex) / REAL(SplitPartNum+1)
            PartMPF(PartIndex)    = PartMPFSplit
          END IF ! SplitPartNum.GT.0
        END IF ! SplitInProgress
      END IF ! usevMPF.AND.BGGas%TraceSpecies(jSpec)

      ! ==============================================================================================================================
      ! Determine collision probability
      ! ==============================================================================================================================
      ! Determine the particle velocity
      IF(BGGas%UseDistribution) THEN
        VeloBGGPart(1:3) = CalcVelocity_maxwell_particle(jSpec,BGGas%Distribution(bgSpec,4:6,iElem)) &
                           + BGGas%Distribution(bgSpec,1:3,iElem)
      ELSE
        CALL CalcVelocity_maxwell_lpn(FractNbr=jSpec, Vec3D=VeloBGGPart(1:3), iInit=1)
      END IF
      CRela2 = (PartState(4,PartIndex) - VeloBGGPart(1))**2 &
             + (PartState(5,PartIndex) - VeloBGGPart(2))**2 &
             + (PartState(6,PartIndex) - VeloBGGPart(3))**2

      IF(BGGas%UseDistribution) THEN
        BGGasNumDens  = BGGas%Distribution(bgSpec,7,iElem)
        BGGasFraction = BGGas%SpeciesFractionElem(bgSpec,iElem)
        IF(XSec_NullCollision) ProbNull = SpecXSec(iCase)%ProbNullElem(iElem)
      ELSE
        BGGasNumDens  = BGGas%NumberDensity(bgSpec)
        BGGasFraction = BGGas%SpeciesFraction(bgSpec)
        IF(XSec_NullCollision) ProbNull = SpecXSec(iCase)%ProbNull
      END IF

      ! Relative kinetic energy of the particle pair (real energy value per particle pair, no weighting/scaling factors)
      IF(CRela2 .LT. RelativisticLimit) THEN
        CollEnergy = 0.5 * CollInf%MassRed(iCase) * CRela2
      ELSE
        ! Relativistic treatment under the assumption that the velocity of the background species is zero or negligible
        GammaFac = CRela2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        CollEnergy = (GammaFac-1.) * CollInf%MassRed(iCase) * c2
      END IF

      ! Set the time step in case of species-specific time stepping
      IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
        dtVar = dt * Species(iSpec)%TimeStepFactor
      ELSE
        dtVar = dt
      END IF

      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! ==========================================================================================
        ! XSec
        ! ==========================================================================================
        ! Interpolate cross-section at the collision energy
        SpecXSec(iCase)%CrossSection = InterpolateCrossSection(SpecXSec(iCase)%CollXSecData,CollEnergy)
        ! Calculate the collision probability
        CollProb = (1. - EXP(-SQRT(CRela2) * SpecXSec(iCase)%CrossSection * BGGasNumDens * dtVar))
        ! Correct the collision probability in the case of the second species being a background species as the number of pairs
        ! is either determined based on the null collision probability or on the species fraction
        IF(XSec_NullCollision) THEN
          CollProb = CollProb / ProbNull
        ELSE
          CollProb = CollProb / BGGasFraction
        END IF
      ELSE
        ! ==========================================================================================
        ! DSMC
        ! ==========================================================================================
        CollProb = CollInf%Coll_SpecPartNum(iSpec)*BGGasNumDens/(1+CollInf%KronDelta(iCase))*CollInf%Cab(iCase) &
                  / CollCaseNum * CRela2 ** (0.5-CollInf%omega(iSpec,jSpec)) * dtVar
        IF(CollisMode.EQ.3) THEN
          ! Chemical reaction with cross-section based probability
          IF(ChemReac%CollCaseInfo(iCase)%HasXSecReaction) THEN
            IF(bggPartIndex.EQ.0) THEN
              DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
              bggPartIndex = GetNextFreePosition()
              IF (bggPartIndex.EQ.0) THEN
                CALL Abort(__STAMP__,'ERROR in MCC: MaxParticleNumber should be increased!')
              END IF
            END IF
            ! If standard collision modelling is used, the reaction probability is added to the collision probability
            CALL MCC_CalcReactionProb(iCase,bgSpec,CRela2,CollEnergy,PartIndex,bggPartIndex,iElem)
            CollProb = CollProb + SUM(ChemReac%CollCaseInfo(iCase)%ReactionProb(:))
            ! If a collision occurs, re-use the energy values set in MCC_CalcReactionProb
            GetInternalEnergy = .FALSE.
          END IF
        END IF

        ! Vibrational excitation
        IF(SpecXSec(iCase)%UseVibXSec) THEN
          ! Calculate the total vibrational cross-section
          nVib = SIZE(SpecXSec(iCase)%VibMode)
          SumVibCrossSection = 0.
          DO iLevel = 1, nVib
            SumVibCrossSection = SumVibCrossSection + InterpolateCrossSection(SpecXSec(iCase)%VibMode(iLevel)%XSecData,CollEnergy)
          END DO
          ! Calculate the total vibrational relaxation probability
          SpecXSec(iCase)%VibProb = 1. - EXP(-SQRT(CRela2) * SumVibCrossSection * BGGasNumDens * dtVar)
          ! Correct the collision probability in the case of the second species being a background species as the number of pairs
          ! is determined based on the species fraction
          SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb / BGGasFraction
          CollProb = CollProb + SpecXSec(iCase)%VibProb
        END IF

        ! Electronic excitation
        IF(SpecXSec(iCase)%UseElecXSec) THEN
          DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
            IF(CollEnergy.GT.SpecXSec(iCase)%ElecLevel(iLevel)%Threshold) THEN
              ! Interpolate the electronic cross-section
              SpecXSec(iCase)%ElecLevel(iLevel)%Prob = InterpolateCrossSection(SpecXSec(iCase)%ElecLevel(iLevel)%XSecData,CollEnergy)
              ! Calculate the electronic excitation probability
              SpecXSec(iCase)%ElecLevel(iLevel)%Prob = 1. - EXP(-SQRT(CRela2) * SpecXSec(iCase)%ElecLevel(iLevel)%Prob &
                                                            * BGGasNumDens * dtVar)
              ! Correct the collision probability in the case of the second species being a background species as the number of pairs
              ! is determined based on the species fraction
              SpecXSec(iCase)%ElecLevel(iLevel)%Prob = SpecXSec(iCase)%ElecLevel(iLevel)%Prob / BGGasFraction
              CollProb = CollProb + SpecXSec(iCase)%ElecLevel(iLevel)%Prob
            END IF
          END DO
        END IF
      END IF ! SpecXSec(iCase)%UseCollXSec
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
          bggPartIndex = GetNextFreePosition()
          IF (bggPartIndex.EQ.0) THEN
            CALL Abort(__STAMP__,'ERROR in MCC: MaxParticleNumber should be increased!')
          END IF
        END IF
        ! Set collision pair index
        Coll_pData(iPair)%iPart_p2 = bggPartIndex
        ! Assign properties
        CALL BGGas_AssignParticleProperties(jSpec,PartIndex,bggPartIndex,GetVelocity_opt=.FALSE.,GetInternalEnergy_opt=GetInternalEnergy)
        PartState(4:6,bggPartIndex) = VeloBGGPart(1:3)
        ! Required in DSMC_Relax_Col_LauxTSHO
        Coll_pData(iPair)%cRela2 = CRela2
        ! Required in ReactionDecision
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
          ! Ambipolar diffusion: add the background particle to consider (as its index might be used for an ion after a reaction)
          IF(DSMC%DoAmbipolarDiff) THEN
            newAmbiParts = newAmbiParts + 1
            iPartIndx_NodeNewAmbi(newAmbiParts) = bggPartIndex
          END IF
          ! Set index to zero to get a new one for the next background gas particle
          bggPartIndex = 0
          GetInternalEnergy = .TRUE.
        END IF
        ! Set index to zero to get a new one for the next split particle
        IF(usevMPF.AND.BGGas%TraceSpecies(jSpec)) PartIndex = 0
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

      ! ==============================================================================================
      ! Determine collision probabilities
      IF(DSMC%CalcQualityFactors) THEN
        DSMC%CollProbMax = MAX(CollProb, DSMC%CollProbMax)
        ! Remove the correction factor for the mean collision probability
        IF(SpecXSec(iSpec)%UseCollXSec) THEN
          IF(XSec_NullCollision) THEN
            CollProb = CollProb * ProbNull
          ELSE
            CollProb = CollProb * BGGasFraction
          END IF
        END IF
        DSMC%CollProbMean = DSMC%CollProbMean + CollProb
        DSMC%CollProbMeanCount = DSMC%CollProbMeanCount + 1
      END IF ! DSMC%CalcQualityFactors
      ! Reservoir simulation: determination of the reaction probabilities
      IF (DSMC%ReservoirSimu) THEN
        ! Sum of collision probabilities for the collision pair, required for the correct reaction rate
        IF(ChemReac%NumOfReact.GT.0) THEN
          IF (ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths.GT.0) THEN
            IF(SpecXSec(iSpec)%UseCollXSec) THEN
              ! Calculate the collision probability for the null collision probability case
              IF(XSec_NullCollision) THEN
                CollProb = CollProb * ProbNull
              ELSE
                CollProb = CollProb * BGGasFraction
              END IF
            END IF
            ChemReac%ReacCollMean(iCase) = ChemReac%ReacCollMean(iCase) + CollProb
          END IF
        END IF  ! ChemReac%NumOfReact.GT.0
      END IF    ! DSMC%ReservoirSimu
    END DO    ! DO WHILE(iLoop.LE.SpecPairNum(iCase))
    SDEALLOCATE(PartIndexCase)
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
                          ElemVolume_Shared(CNElemID), DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the MCS/MFP for the case without octree
    IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                    / DSMC%MeanFreePath
  END IF
END IF

IF (DSMC%DoAmbipolarDiff) THEN
  CALL AD_DeleteParticles(iPartIndx_NodeTotalAmbiDel,nPart)
END IF

DEALLOCATE(iPartIndx_Node)
SDEALLOCATE(iPartIndx_NodeTotalAmbiDel)
DEALLOCATE(iPartIndexSpec)
DEALLOCATE(SpecPartNum)
DEALLOCATE(UseSpecPartNum)
DEALLOCATE(SpecPairNum)
DEALLOCATE(Coll_pData)

END SUBROUTINE MonteCarloCollision


!===================================================================================================================================
!> Calculate the reaction probability if collision cross-section data is used (only with a background gas from the MCC routine)
!===================================================================================================================================
SUBROUTINE MCC_CalcReactionProb(iCase,bgSpec,CRela2,CollEnergy_in,PartIndex,bggPartIndex,iElem)
! MODULES
USE MOD_Particle_Vars         ,ONLY: Species, PartSpecies, VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, BGGas, ChemReac, DSMC, PartStateIntEn
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_part_tools            ,ONLY: CalcERot_particle, CalcEVib_particle, CalcEElec_particle
USE MOD_MCC_XSec              ,ONLY: InterpolateCrossSection
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase,bgSpec,PartIndex,bggPartIndex,iElem
REAL,INTENT(IN)               :: CRela2, CollEnergy_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: jSpec, iPath, ReacTest, EductReac(1:3), ProductReac(1:4), iProd
INTEGER                       :: NumWeightProd
REAL                          :: EZeroPoint_Educt, EZeroPoint_Prod, CollEnergy
REAL                          :: CrossSection, dtVar
REAL                          :: Temp_Rot, Temp_Vib, Temp_Elec, BGGasNumDens, BGGasFraction
!===================================================================================================================================
NumWeightProd = 2

jSpec = BGGas%MapBGSpecToSpec(bgSpec)

! Set the time step in case of species-specific time stepping
IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
  dtVar = dt * Species(PartSpecies(PartIndex))%TimeStepFactor
ELSE
  dtVar = dt
END IF

DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
  IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'XSec') THEN
    EductReac(1:3) = ChemReac%Reactants(ReacTest,1:3); ProductReac(1:4) = ChemReac%Products(ReacTest,1:4)

    ! Sum of the zero-point energies of the reactants
    EZeroPoint_Educt = 0.0; EZeroPoint_Prod = 0.0
    IF((SpecDSMC(EductReac(1))%InterID.EQ.2).OR.(SpecDSMC(EductReac(1))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(1))%EZeroPoint
    END IF
    IF((SpecDSMC(EductReac(2))%InterID.EQ.2).OR.(SpecDSMC(EductReac(2))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(2))%EZeroPoint
    END IF
    ! Sum of the zero-point energies of the products
    IF(ProductReac(4).NE.0) THEN
      ! 4 Products
      NumWeightProd = 4
    ELSE IF(ProductReac(3).NE.0) THEN
      ! 3 Products
      NumWeightProd = 3
    END IF
    DO iProd = 1, NumWeightProd
      IF((SpecDSMC(ProductReac(iProd))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(iProd))%InterID.EQ.20)) THEN
        EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(iProd))%EZeroPoint
      END IF
    END DO
    ! Adding the internal energy of particle species
    CollEnergy = CollEnergy_in + PartStateIntEn(1,PartIndex) + PartStateIntEn(2,PartIndex)
    ! Internal energy of background species
    IF((SpecDSMC(jSpec)%InterID.EQ.2).OR.(SpecDSMC(jSpec)%InterID.EQ.20)) THEN
      IF(BGGas%UseDistribution) THEN
        Temp_Vib   = BGGas%Distribution(bgSpec,8,iElem)
        Temp_Rot   = BGGas%Distribution(bgSpec,9,iElem)
      ELSE
        Temp_Vib   = SpecDSMC(jSpec)%Init(1)%TVib
        Temp_Rot   = SpecDSMC(jSpec)%Init(1)%TRot
      END IF
      PartStateIntEn(1,bggPartIndex) = CalcEVib_particle(jSpec,Temp_Vib,bggPartIndex)
      PartStateIntEn(2,bggPartIndex) = CalcERot_particle(jSpec,Temp_Rot)
      CollEnergy = CollEnergy + PartStateIntEn(1,bggPartIndex) + PartStateIntEn(2,bggPartIndex)
    END IF
    IF ((DSMC%ElectronicModel.GT.0).AND.(.NOT.SpecDSMC(jSpec)%FullyIonized)) THEN
      IF(BGGas%UseDistribution) THEN
        Temp_Elec = BGGas%Distribution(bgSpec,10,iElem)
      ELSE
        Temp_Elec = SpecDSMC(jSpec)%Init(1)%TElec
      END IF
      PartStateIntEn(3,bggPartIndex) = CalcEElec_particle(jSpec,Temp_Elec,bggPartIndex)
      CollEnergy = CollEnergy + PartStateIntEn(3,PartIndex) + PartStateIntEn(3,bggPartIndex)
    END IF
    ! Check first if sufficient energy is available for the products after the reaction
    IF(((CollEnergy-EZeroPoint_Prod).GE.-ChemReac%EForm(ReacTest))) THEN
      CollEnergy = CollEnergy - EZeroPoint_Educt
      CrossSection = InterpolateCrossSection(SpecXSec(iCase)%ReactionPath(iPath)%XSecData,CollEnergy)
      ASSOCIATE( ReactionProb => ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) )
        IF(SpecXSec(iCase)%UseCollXSec) THEN
          ! Interpolate the reaction cross-section
          ReactionProb = CrossSection
        ELSE
          ! Calculate the reaction probability
          IF(BGGas%UseDistribution) THEN
            BGGasNumDens  = BGGas%Distribution(bgSpec,7,iElem)
            BGGasFraction = BGGas%SpeciesFractionElem(bgSpec,iElem)
          ELSE
            BGGasNumDens  = BGGas%NumberDensity(bgSpec)
            BGGasFraction = BGGas%SpeciesFraction(bgSpec)
          END IF
          ReactionProb = 1. - EXP(-SQRT(CRela2) * dtVar * BGGasNumDens * CrossSection)
          ! Correct the reaction probability in the case of the second species being a background species as the number of pairs
          ! is based on the species fraction
          ReactionProb = ReactionProb / BGGasFraction
        END IF
      END ASSOCIATE
    ELSE
      ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = 0.
    END IF
    ! Reservoir simulation: Calculation of reaction rate coefficient
    IF (DSMC%ReservoirSimu.AND..NOT.DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(ReacTest) = ChemReac%NumReac(ReacTest) + ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath)
      ChemReac%ReacCount(ReacTest) = ChemReac%ReacCount(ReacTest) + 1
    END IF
  END IF
END DO

END SUBROUTINE MCC_CalcReactionProb

END MODULE MOD_MCC
