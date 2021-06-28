!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_DSMC_BGGas
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGGas_InsertParticles
  MODULE PROCEDURE BGGas_InsertParticles
END INTERFACE

INTERFACE DSMC_pairing_bggas
  MODULE PROCEDURE DSMC_pairing_bggas
END INTERFACE

INTERFACE BGGas_ControlParticles
  MODULE PROCEDURE BGGas_ControlParticles
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGGas_Initialize, BGGas_InsertParticles, DSMC_pairing_bggas, MCC_pairing_bggas, BGGas_ControlParticles
PUBLIC :: BGGas_PhotoIonization
!===================================================================================================================================

CONTAINS

SUBROUTINE BGGas_Initialize()
!===================================================================================================================================
!> Initialization of the background gas: compatibility check, array allocation, background species to species mapping and
!> calculation of the molar fraction
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_DSMC_Vars              ,ONLY: BGGas
USE MOD_Particle_Vars          ,ONLY: PDM, Symmetry, Species, nSpecies, VarTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec, counterSpec
REAL              :: SpeciesDensTemp(1:nSpecies)
!===================================================================================================================================

! 1.) Check compatibility with other features and whether required parameters have been read-in
IF((Symmetry%Order.EQ.2).OR.VarTimeStep%UseVariableTimeStep) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: 2D/Axisymmetric and variable timestep are not implemented with a background gas yet!')
END IF

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    IF (BGGas%NumberDensity(iSpec).EQ.0.) CALL abort(__STAMP__&
                                          ,'ERROR: NumberDensity is zero but must be defined for a background gas!')
    IF (Species(iSpec)%NumberOfInits.NE.1) &
      CALL abort(&
        __STAMP__&
        ,'ERROR: BGG species can be used ONLY for BGG!')
  END IF
END DO

! 2.) Allocation
ALLOCATE(BGGas%PairingPartner(PDM%maxParticleNumber))
BGGas%PairingPartner = 0
ALLOCATE(BGGas%MapSpecToBGSpec(nSpecies))
BGGas%MapSpecToBGSpec = 0
SpeciesDensTemp(1:nSpecies) = BGGas%NumberDensity(1:nSpecies)
DEALLOCATE(BGGas%NumberDensity)
ALLOCATE(BGGas%NumberDensity(BGGas%NumberOfSpecies))
ALLOCATE(BGGas%SpeciesFraction(BGGas%NumberOfSpecies))

! 3.) Create a mapping of background species to regular species and calculate the molar fraction
counterSpec = 0
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    counterSpec = counterSpec + 1
    BGGas%MapSpecToBGSpec(iSpec) = counterSpec
    BGGas%NumberDensity(counterSpec) = SpeciesDensTemp(iSpec)
    BGGas%SpeciesFraction(counterSpec) = BGGas%NumberDensity(counterSpec) / SUM(SpeciesDensTemp)
    IF(counterSpec.GT.BGGas%NumberOfSpecies) THEN
      CALL Abort(&
        __STAMP__&
        ,'ERROR in BGGas: More background species detected than previously defined!')
    END IF
  END IF
END DO

END SUBROUTINE BGGas_Initialize


INTEGER FUNCTION BGGas_GetSpecies()
!===================================================================================================================================
!> Get a species index of the background gas by randomly choosing a species based on the molar fraction
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: BGGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: iRan
INTEGER           :: iSpec
!===================================================================================================================================

CALL RANDOM_NUMBER(iRan)
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    IF(BGGas%NumberOfSpecies.GT.1) THEN
      IF(SUM(BGGas%SpeciesFraction(1:BGGas%MapSpecToBGSpec(iSpec))).GT.iRan) THEN
        BGGas_GetSpecies = iSpec
        RETURN
      END IF
    ELSE
      BGGas_GetSpecies = iSpec
    END IF
  END IF
END DO

END FUNCTION BGGas_GetSpecies

SUBROUTINE BGGas_InsertParticles()
!===================================================================================================================================
!> Creating particles of the background species for each actual simulation particle
!> 1. Initialize particles (loop over the current ParticleVecLength):
!>    a) Get the new index from the nextFreePosition array
!>    b) Same position as the non-BGG particle, set species and initialize the internal energy (if required)
!>    c) Include BGG-particles in PEM%p-Lists (counting towards the amount of particles per cell stored in PEM%pNumber)
!>    d) Map BGG-particle to the non-BGG particle for particle pairing
!> 2. Call SetParticleVelocity: loop over the newly created particles to set the thermal velocity, using nextFreePosition to get
!>    the same indices, possible since the CurrentNextFreePosition was not updated yet
!> 3. Adjust ParticleVecLength and currentNextFreePosition
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_Vars              ,ONLY: BGGas, SpecDSMC, CollisMode
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_PARTICLE_Vars          ,ONLY: PDM, PartSpecies, PartState, PEM, PartPosRef, usevMPF, PartMPF, VarTimeStep
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr, iSpec, LocalElemID
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

iNewPart=0
PositionNbr = 0
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Skip background particles that have been created within this loop
    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) CYCLE
    iNewPart = iNewPart + 1
    PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
    IF (PositionNbr.EQ.0) THEN
      CALL Abort(&
        __STAMP__&
        ,'ERROR in BGGas: MaxParticleNumber should be twice the expected number of particles, to account for the BGG particles!')
    END IF
    PartState(1:3,PositionNbr) = PartState(1:3,iPart)
    IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
      PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
    END IF
    iSpec = BGGas_GetSpecies()
    PartSpecies(PositionNbr) = iSpec
    IF(CollisMode.GT.1) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        CALL DSMC_SetInternalEnr_Poly(iSpec,1,PositionNbr,1)
      ELSE
        CALL DSMC_SetInternalEnr_LauxVFD(iSpec,1,PositionNbr,1)
      END IF
    END IF        
    IF(usevMPF) PartMPF(PositionNbr) = PartMPF(iPart)
    IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(PositionNbr) = VarTimeStep%ParticleTimeStep(iPart)
    PEM%GlobalElemID(PositionNbr) = PEM%GlobalElemID(iPart)
    LocalElemID = PEM%LocalElemID(PositionNbr)
    PDM%ParticleInside(PositionNbr) = .true.
    PEM%pNext(PEM%pEnd(LocalElemID)) = PositionNbr     ! Next Particle of same Elem (Linked List)
    PEM%pEnd(LocalElemID) = PositionNbr
    PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
    BGGas%PairingPartner(iPart) = PositionNbr
    CALL CalcVelocity_maxwell_lpn(FractNbr=iSpec, Vec3D=PartState(4:6,PositionNbr), iInit=1)
  END IF
END DO
PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE BGGas_InsertParticles


SUBROUTINE DSMC_pairing_bggas(iElem)
!===================================================================================================================================
! Building of pairs for the background gas
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze          ,ONLY: CalcGammaVib, CalcMeanFreePath
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SelectionProc
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, RadialWeighting
USE MOD_DSMC_PolyAtomicModel  ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_emission_tools   ,ONLY: DSMC_SetInternalEnr_LauxVFD, CalcVelocity_maxwell_lpn
USE MOD_Particle_Vars         ,ONLY: PEM,PartSpecies,nSpecies,PartState,Species,usevMPF,PartMPF,Species, WriteMacroVolumeValues
USE MOD_Particle_Vars         ,ONLY: PartPosRef, PDM, VarTimeStep
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_DSMC_Collis           ,ONLY: DSMC_perform_collision
USE MOD_DSMC_Relaxation       ,ONLY: FinalizeCalcVibRelaxProb, SumVibRelaxProb, InitCalcVibRelaxProb
USE MOD_TimeDisc_Vars         ,ONLY: TEnd, time
USE MOD_DSMC_CollisionProb    ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Relaxation       ,ONLY: CalcMeanVibQuaDiatomic
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nPair, iPair, iPart, iLoop, nPart, iLoop2, iNewPart, iSplit, LocalElemID, CNElemID
INTEGER                       :: cSpec1, cSpec2, iCase, iSpec, iSpec2, iSpec3, iSplitPart
REAL                          :: iRan, MPF
LOGICAL                       :: NeedToSplit
INTEGER                       :: NewPairNum, SplitPartNum, PositionNbr
INTEGER, ALLOCATABLE          :: TempIndxArray(:)
!===================================================================================================================================
nPart = PEM%pNumber(iElem)
NeedToSplit = .FALSE.
IF(usevMPF) THEN
  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    iSpec  = PartSpecies(iPart)! Skip background particles that have been created within this loop
    IF(BGGas%BackgroundSpecies(iSpec)) CYCLE
    iSpec2 = PartSpecies(BGGas%PairingPartner(iPart))
    IF(PartMPF(iPart).GT.Species(iSpec2)%MacroParticleFactor) THEN 
      NeedToSplit = .TRUE.
      EXIT
    END IF
    iPart = PEM%pNext(iPart)
  END DO
END IF

IF(NeedToSplit) THEN
  iPart = PEM%pStart(iElem)
  NewPairNum = INT(nPart/2)
  CollInf%Coll_CaseNum = 0
  iNewPart = 0

  DO iLoop = 1, nPart
    iSpec  = PartSpecies(iPart)! Skip background particles that have been created within this loop
    IF(BGGas%BackgroundSpecies(iSpec)) CYCLE
    iSpec2 = PartSpecies(BGGas%PairingPartner(iPart))
    IF(PartMPF(iPart).GT.Species(iSpec2)%MacroParticleFactor) THEN
      SplitPartNum = NINT(PartMPF(iPart) / Species(iSpec2)%MacroParticleFactor)
      NewPairNum   = NewPairNum + (SplitPartNum - 1)          ! New Number of Pairs  
      PositionNbr = 0
      PartMPF(iPart) = Species(iSpec2)%MacroParticleFactor ! set new MPF for test particle
      PartMPF(BGGas%PairingPartner(iPart)) = PartMPF(iPart)
      ALLOCATE(TempIndxArray(SplitPartNum))
      iSplitPart = 0
      DO iLoop2 = 1, 2    ! 1: split test particle, 2: Creating particles of the background species
        DO iSplit = 2, SplitPartNum
          iNewPart = iNewPart + 1
          iSplitPart = iSplitPart + 1
          PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr.EQ.0) THEN
            CALL Abort(&
              __STAMP__&
              ,'ERROR in BGGas: MaxParticleNumber should be twice the expected number of particles, to account for the BGG particles!')
          END IF
          PartState(1:3,PositionNbr) = PartState(1:3,iPart)
          IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
            PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
          END IF
          IF (iLoop2.EQ.1) THEN
            iSpec3 = iSpec
            PartState(4:6,PositionNbr) = PartState(4:6,iPart)
            IF(CollisMode.GT.1) THEN
              PartStateIntEn(1:2,PositionNbr) = PartStateIntEn(1:2,iPart)
              IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PositionNbr) = PartStateIntEn(3,iPart)
            END IF
          ELSE
            iSpec3 = iSpec2
            CALL CalcVelocity_maxwell_lpn(FractNbr=iSpec3, Vec3D=PartState(4:6,PositionNbr), iInit=1)
            IF(CollisMode.GT.1) THEN
              IF(SpecDSMC(iSpec3)%PolyatomicMol) THEN
                CALL DSMC_SetInternalEnr_Poly(iSpec3,1,PositionNbr,1)
              ELSE
                CALL DSMC_SetInternalEnr_LauxVFD(iSpec3,1,PositionNbr,1)
              END IF
            END IF
          END IF
          PartSpecies(PositionNbr) = iSpec3
          PartMPF(PositionNbr) = PartMPF(iPart)
          IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(PositionNbr) = VarTimeStep%ParticleTimeStep(iPart)
          PEM%GlobalElemID(PositionNbr) = PEM%GlobalElemID(iPart)
          PEM%LastGlobalElemID(PositionNbr) = PEM%GlobalElemID(iPart)
          LocalElemID = PEM%LocalElemID(PositionNbr)
          PDM%ParticleInside(PositionNbr) = .TRUE.
          PDM%IsNewPart(PositionNbr)       = .TRUE.
          PDM%dtFracPush(PositionNbr)      = .FALSE.
          PEM%pNext(PEM%pEnd(LocalElemID)) = PositionNbr     ! Next Particle of same Elem (Linked List)
          PEM%pEnd(LocalElemID) = PositionNbr
          PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
          IF (iLoop2.EQ.1) THEN
            TempIndxArray(iSplitPart) = PositionNbr
          ELSE
            BGGas%PairingPartner(TempIndxArray(iSplitPart-(SplitPartNum - 1))) = PositionNbr
          END IF  
        END DO   
      END DO 
      DEALLOCATE(TempIndxArray)
    END IF
    iPart = PEM%pNext(iPart)
  END DO
  nPair = NewPairNum
  PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart  
  nPart = PEM%pNumber(iElem)
  CALL InitCalcVibRelaxProb()
  ALLOCATE(Coll_pData(NewPairNum))
  Coll_pData%Ec=0
  CollInf%Coll_SpecPartNum = 0

  iPair = 1

  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iSpec = PartSpecies(iPart)
    ! Counting the number of particles per species
    MPF = GetParticleWeight(iPart)
    CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + MPF
    ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(1,iPart) * MPF
    ! Creating pairs for species, which are not the background species
    IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN 
      Coll_pData(iPair)%iPart_p1 = iPart
      Coll_pData(iPair)%iPart_p2 = BGGas%PairingPartner(iPart)
      iPair = iPair + 1
    END IF
    iPart = PEM%pNext(iPart)
  END DO
  
    
ELSE IF(usevMPF) THEN    ! No need for additional split but variable weighting
  nPair = INT(nPart/2)
  CALL InitCalcVibRelaxProb()

  CollInf%Coll_SpecPartNum = 0.
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  Coll_pData%Ec=0
  iPair = 1

  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iSpec = PartSpecies(iPart)
    ! Counting the number of particles per species
    MPF = GetParticleWeight(iPart)
    CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + MPF
    ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(1,iPart) * MPF  ! siehe dsmc_relaxation.f90 line 112
    ! Creating pairs for species, which are not the background species
    IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN 
      Coll_pData(iPair)%iPart_p1 = iPart
      Coll_pData(iPair)%iPart_p2 = BGGas%PairingPartner(iPart)
      iPair = iPair + 1
    END IF
    iPart = PEM%pNext(iPart)
  END DO
  
  
ELSE    ! standard 
  
  
  nPart = PEM%pNumber(iElem)
  nPair = INT(nPart/2)
  CALL InitCalcVibRelaxProb()

  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  Coll_pData%Ec=0
  iPair = 1

  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iSpec = PartSpecies(iPart)
    ! Counting the number of particles per species
    MPF = GetParticleWeight(iPart)
    CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + MPF
    ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(1,iPart) * MPF
    ! Creating pairs for species, which are not the background species
    IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN
      Coll_pData(iPair)%iPart_p1 = iPart
      Coll_pData(iPair)%iPart_p2 = BGGas%PairingPartner(iPart)
      iPair = iPair + 1
    END IF
    iPart = PEM%pNext(iPart)
  END DO
END IF



IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors) THEN
  ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
  !           relaxation probability (CalcGammaVib)
  ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
  ! 3. Case: Temperature required for the mean free path with the VHS model
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = 0.
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) &
                                    + BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(iSpec)) * Species(iSpec)%Init(1)%MWTemperatureIC
    END IF
  END DO
  IF(SelectionProc.EQ.2) CALL CalcGammaVib()
END IF

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    CNElemID = GetCNElemID(iElem+offSetElem)
    iSpec2   = BGGas%MapSpecToBGSpec(iSpec)
    IF(usevMPF) THEN
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(iSpec2)*ElemVolume_Shared(CNElemID)
    ELSE
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(iSpec2)*ElemVolume_Shared(CNElemID)/Species(iSpec)%MacroParticleFactor
    END IF
  END IF
END DO

CollInf%MeanMPF = 0.

DO iPair = 1, nPair
  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + GetParticleWeight(Coll_pData(iPair)%iPart_p1)
  Coll_pData(iPair)%CRela2 = (PartState(4,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(4,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(5,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(5,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(6,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(6,Coll_pData(iPair)%iPart_p2))**2
  Coll_pData(iPair)%PairType = iCase
  Coll_pData(iPair)%NeedForRec = .FALSE.
END DO

IF(CollisMode.EQ.3) CALL CalcMeanVibQuaDiatomic()

! 6.) Calculate the collision probability and perform the collision if necessary

DO iPair = 1, nPair
  CALL SumVibRelaxProb(iPair)
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    CALL DSMC_prob_calc(iElem, iPair)
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
      CALL DSMC_perform_collision(iPair,iElem)
    END IF
  END IF
END DO
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum), &
                          ElemVolume_Shared(GetCNElemID(iElem+offSetElem)), DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the MCS/MFP for the case without octree
    IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                    / DSMC%MeanFreePath
  END IF
END IF
DEALLOCATE(Coll_pData)
CALL FinalizeCalcVibRelaxProb(iElem)

END SUBROUTINE DSMC_pairing_bggas


SUBROUTINE MCC_pairing_bggas(iElem)
!===================================================================================================================================
!> Routine to create pairs with particles from the background gas using read-in collision cross sections (null collision method) and
!> conventional VHS model. If the null collision method is used, the number of pairs is determined by a constant collision
!> probability at the maximal collision frequency. For the regular background gas case, for every particle a pair is created.
!> 1.) Counting the number of particles per species and creating a species-specific particle index list
!> 2.) Determining the total number of pairs
!> 3a.) Creating the background particles as required by the determined numbers of collision pairs
!> 3b.) Pairing the newly created background particles with the actual simulation particles
!> 4.) Determine the particle number of the background species and calculate the cell temperature
!> 5.) Calculate the square of the relative collision velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze            ,ONLY: CalcGammaVib, CalcMeanFreePath
USE MOD_part_emission_tools     ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SpecXSec
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, MCC_TotalPairNum, DSMCSumOfFormedParticles, XSec_NullCollision
USE MOD_Part_Emission_Tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_Part_Pos_and_Velo       ,ONLY: SetParticleVelocity
USE MOD_Particle_Vars           ,ONLY: PEM, PDM, PartSpecies, nSpecies, PartState, Species, usevMPF, PartMPF, Species, PartPosRef
USE MOD_Particle_Vars           ,ONLY: VarTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars           ,ONLY: WriteMacroVolumeValues
USE MOD_DSMC_Collis             ,ONLY: DSMC_perform_collision
USE MOD_DSMC_Relaxation         ,ONLY: FinalizeCalcVibRelaxProb, SumVibRelaxProb, InitCalcVibRelaxProb
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, time
USE MOD_DSMC_CollisionProb      ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Relaxation         ,ONLY: CalcMeanVibQuaDiatomic
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_InsertParticles, AD_DeleteParticles
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPair, iPart, iLoop, nPart, iSpec, jSpec, bgSpec, PartIndex, bggPartIndex, PairCount, RandomPart
INTEGER                       :: cSpec1, cSpec2, iCase, SpecPairNumTemp, nPartAmbi, OldPairNum, CNElemID
INTEGER                       :: iPart_p1, iPart_p2, iPairNew
INTEGER,ALLOCATABLE           :: iPartIndexSpec(:,:), SpecPartNum(:), SpecPairNum(:)
REAL                          :: iRan, ProbRest, SpecPairNumReal, MPF
INTEGER, ALLOCATABLE          :: iPartIndx_NodeTotalAmbiDel(:)
INTEGER, ALLOCATABLE, TARGET  :: iPartIndx_Node(:), iPartIndx_NodeTotalAmbi(:)
INTEGER, POINTER              :: iPartIndx_NodeTotal(:)
LOGICAL, ALLOCATABLE          :: NeedToSplit(:)
REAL, ALLOCATABLE             :: CollpDataMPF(:)
INTEGER, ALLOCATABLE          :: CollpDataPairIndex(:,:), SplitPartNum(:)
!===================================================================================================================================
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

! 2.) Determining the total number of pairs
DO iSpec = 1,nSpecies
  IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN    ! Loop over all non-background species
    DO jSpec = 1, nSpecies
      IF(BGGas%BackgroundSpecies(jSpec)) THEN     ! Loop over all background species
        iCase = CollInf%Coll_Case(iSpec,jSpec)
        bgSpec = BGGas%MapSpecToBGSpec(jSpec)
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
      END IF
    END DO
  END IF
END DO

ALLOCATE(CollpDataPairIndex(1:2,1:MCC_TotalPairNum))
ALLOCATE(CollpDataMPF(1:MCC_TotalPairNum))
! ALLOCATE(PairingPartner(MCC_TotalPairNum))
! PairingPartner = 0
PartIndex = 0
PairCount = 0
ALLOCATE(NeedToSplit(MCC_TotalPairNum))
NeedToSplit = .FALSE.
ALLOCATE(SplitPartNum(MCC_TotalPairNum))
SplitPartNum = 0

OldPairNum = MCC_TotalPairNum

! 3a.) Creating the background particles as required by the determined numbers of collision pairs
! 3b.) Pairing the newly created background particles with the actual simulation particles
DO iSpec = 1,nSpecies                             ! Loop over all non-background species
  IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN
    DO jSpec = 1, nSpecies                        ! Loop over all background species
      IF(BGGas%BackgroundSpecies(jSpec)) THEN
        iCase = CollInf%Coll_Case(iSpec,jSpec)
        DO iLoop = 1, SpecPairNum(iCase)    ! Loop over all the number of pairs required for this species pairing
          ! Choosing random particles from the available number of particles, getting the index of the simulation particle
          IF(SpecPartNum(iSpec).GT.0) THEN
            CALL RANDOM_NUMBER(iRan)
            RandomPart = INT(SpecPartNum(iSpec)*iRan) + 1
            PartIndex = iPartIndexSpec(RandomPart,iSpec)
            iPartIndexSpec(RandomPart, iSpec) = iPartIndexSpec(SpecPartNum(iSpec),iSpec)
            SpecPartNum(iSpec) = SpecPartNum(iSpec) - 1
          END IF
          ! Creating a new background gas particle
          DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
          bggPartIndex = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
          IF (bggPartIndex.EQ.0) THEN
            CALL Abort(__STAMP__,&
              'ERROR in MCC: MaxParticleNumber should be twice the expected number of particles, to account for the BGG/MCC particles!')
          END IF
          ! Position the background particle at the simulation particle
          PartState(1:3,bggPartIndex) = PartState(1:3,PartIndex)
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
          ! Determine the particle velocity
          CALL CalcVelocity_maxwell_lpn(FractNbr=jSpec, Vec3D=PartState(4:6,bggPartIndex), iInit=1)
          ! Advance the total count
          PairCount = PairCount + 1
          ! Pairing
          CollpDataPairIndex(1,PairCount) = PartIndex
          CollpDataPairIndex(2,PairCount) = bggPartIndex
          ! Ambipolar diffusion: add the background particle to consider (as its index might be used for an ion after a reaction)
          IF(DSMC%DoAmbipolarDiff) THEN
            newAmbiParts = newAmbiParts + 1
            iPartIndx_NodeNewAmbi(newAmbiParts) = bggPartIndex
          END IF
          IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(bggPartIndex) = VarTimeStep%ParticleTimeStep(PartIndex)
          ! Species-specific vMPF
          IF(usevMPF) THEN
            ! Set the background gas species and mark it for split
            IF(PartMPF(PartIndex).GT.Species(jSpec)%MacroParticleFactor) THEN
              PartMPF(bggPartIndex)   = Species(jSpec)%MacroParticleFactor
              NeedToSplit(PairCount)  = .TRUE.
              SplitPartNum(PairCount) = NINT(PartMPF(PartIndex) / PartMPF(bggPartIndex)) - 1
              MCC_TotalPairNum        = MCC_TotalPairNum + SplitPartNum(PairCount)
              PartMPF(PartIndex)      = PartMPF(bggPartIndex)
            ELSE
              PartMPF(bggPartIndex) = PartMPF(PartIndex)
            END IF
            CollpDataMPF(PairCount) = PartMPF(bggPartIndex)
          END IF
        END DO
      END IF
    END DO
  END IF
END DO

ALLOCATE(Coll_pData(MCC_TotalPairNum))
Coll_pData%Ec = 0.

! Storing the temporary variables in the global Coll_pData
Coll_pData(1:OldPairNum)%iPart_p1 = CollpDataPairIndex(1,1:OldPairNum)
Coll_pData(1:OldPairNum)%iPart_p2 = CollpDataPairIndex(2,1:OldPairNum)

DEALLOCATE(CollpDataPairIndex)
DEALLOCATE(CollpDataMPF)

PairCount = OldPairNum

IF(ANY(NeedToSplit)) THEN
  DO iPair = 1, OldPairNum
    IF(NeedToSplit(iPair)) THEN
      iPart_p1 = Coll_pData(iPair)%iPart_p1
      iPart_p2 = Coll_pData(iPair)%iPart_p2
      DO iPairNew = 1, SplitPartNum(iPair)
        ! Increase number of total pair counts
        PairCount = PairCount + 1
        ! Creating a new split gas particle (cloning)
        DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
        PartIndex = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
        IF (PartIndex.EQ.0) THEN
          CALL Abort(__STAMP__,&
            'ERROR in MCC: MaxParticleNumber should be twice the expected number of particles, to account for the BGG/MCC particles!')
        END IF
        !
        iSpec = PartSpecies(iPart_p1)
        PartSpecies(PartIndex) = iSpec
        PartState(1:6,PartIndex) = PartState(1:6,iPart_p1)
        IF(CollisMode.GT.1) THEN
          PartStateIntEn(1:2,PartIndex) = PartStateIntEn(1:2,iPart_p1)
          IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PartIndex) = PartStateIntEn(3,iPart_p1)
        END IF
        !
        PartMPF(PartIndex) = PartMPF(iPart_p1)
        IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(PartIndex) = VarTimeStep%ParticleTimeStep(iPart_p1)
        PEM%GlobalElemID(PartIndex) = PEM%GlobalElemID(iPart_p1)
        PEM%LastGlobalElemID(PartIndex) = PEM%GlobalElemID(iPart_p1)
        PDM%ParticleInside(PartIndex)  = .TRUE.
        PDM%IsNewPart(PartIndex)       = .TRUE.
        PDM%dtFracPush(PartIndex)      = .FALSE.
        ! Creating a new background gas particle
        DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
        bggPartIndex = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
        IF (bggPartIndex.EQ.0) THEN
          CALL Abort(__STAMP__,&
            'ERROR in MCC: MaxParticleNumber should be twice the expected number of particles, to account for the BGG/MCC particles!')
        END IF
        ! 
        bgSpec = PartSpecies(iPart_p2)
        PartSpecies(bggPartIndex) = bgSpec
        PartState(1:3,bggPartIndex) = PartState(1:3,iPart_p1)
        CALL CalcVelocity_maxwell_lpn(FractNbr=bgSpec, Vec3D=PartState(4:6,bggPartIndex), iInit=1)
        IF(CollisMode.GT.1) THEN
          IF(SpecDSMC(bgSpec)%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(bgSpec,1,bggPartIndex,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(bgSpec,1,bggPartIndex,1)
          END IF
        END IF
        ! 
        PartMPF(bggPartIndex) = PartMPF(PartIndex)
        IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(bggPartIndex) = VarTimeStep%ParticleTimeStep(iPart_p1)
        PEM%GlobalElemID(bggPartIndex) = PEM%GlobalElemID(iPart_p1)
        PEM%LastGlobalElemID(bggPartIndex) = PEM%GlobalElemID(iPart_p1)
        PDM%ParticleInside(bggPartIndex)  = .TRUE.
        PDM%IsNewPart(bggPartIndex)       = .TRUE.
        PDM%dtFracPush(bggPartIndex)      = .FALSE.
        !
        Coll_pData(PairCount)%iPart_p1 = PartIndex
        Coll_pData(PairCount)%iPart_p2 = bggPartIndex
      END DO
    END IF
  END DO
END IF

! 4.) Determine the particle number of the background species and calculate the cell temperature
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    CNElemID = GetCNElemID(iElem+offSetElem)
    bgSpec   = BGGas%MapSpecToBGSpec(iSpec)
    IF(usevMPF) THEN
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bgSpec)*ElemVolume_Shared(CNElemID)
    ELSE
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bgSpec)*ElemVolume_Shared(CNElemID)/Species(iSpec)%MacroParticleFactor
    END IF
  END IF
END DO

IF(DSMC%CalcQualityFactors) THEN
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = 0.
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) &
                                    + BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(iSpec)) * Species(iSpec)%Init(1)%MWTemperatureIC
    END IF
  END DO
END IF

! 5.) Calculate the square of the relative collision velocity
DO iPair = 1, MCC_TotalPairNum
  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + GetParticleWeight(Coll_pData(iPair)%iPart_p1)
  CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
  Coll_pData(iPair)%CRela2 = (PartState(4,Coll_pData(iPair)%iPart_p1) - PartState(4,Coll_pData(iPair)%iPart_p2))**2 &
                           + (PartState(5,Coll_pData(iPair)%iPart_p1) - PartState(5,Coll_pData(iPair)%iPart_p2))**2 &
                           + (PartState(6,Coll_pData(iPair)%iPart_p1) - PartState(6,Coll_pData(iPair)%iPart_p2))**2
  Coll_pData(iPair)%PairType = iCase
  Coll_pData(iPair)%NeedForRec = .FALSE.
END DO

IF(CollisMode.EQ.3) CALL CalcMeanVibQuaDiatomic()

! 6.) Calculate the collision probability and perform the collision if necessary
DO iPair = 1, MCC_TotalPairNum
  CALL SumVibRelaxProb(iPair)
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    CALL DSMC_prob_calc(iElem, iPair)
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.GE.iRan) THEN
      CALL DSMC_perform_collision(iPair,iElem)
    END IF
  END IF
END DO
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum), &
                          ElemVolume_Shared(GetCNElemID(iElem+offSetElem)), DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the MCS/MFP for the case without octree
    IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                    / DSMC%MeanFreePath
  END IF
END IF
DEALLOCATE(Coll_pData)
CALL FinalizeCalcVibRelaxProb(iElem)

IF (DSMC%DoAmbipolarDiff) THEN
  CALL AD_DeleteParticles(iPartIndx_NodeTotalAmbiDel,nPart)
END IF

END SUBROUTINE MCC_pairing_bggas


SUBROUTINE BGGas_ControlParticles()
!===================================================================================================================================
! Deletes all background gas particles and updates the particle index list
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies,vMPFNewPartNum
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
USE MOD_vMPF,               ONLY : SplitMerge_main
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart
!===================================================================================================================================
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) PDM%ParticleInside(iPart) = .FALSE.
  END IF
END DO
BGGas%PairingPartner = 0

IF(vMPFNewPartNum.NE.0) CALL SplitMerge_main

CALL UpdateNextFreePosition()

END SUBROUTINE BGGas_ControlParticles


SUBROUTINE BGGas_PhotoIonization(iSpec,iInit,TotalNbrOfReactions)
!===================================================================================================================================
!> Particle emission through photo ionization, defined by chemical reactions of the type "phIon" and an emission with SpaceIC
!> "photon_cylinder" for the 
!> 0.) Determine the total ionization cross-section
!> 1.) Compute the number of photoionization events in the local domain of each proc
!> 2.) Delete left-over inserted particles
!> 3.) Insert the products of the photoionization rections
!> 4.) Perform the reaction, distribute the collision energy (including photon energy) and emit electrons perpendicular
!>     to the photon's path
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze           ,ONLY: CalcGammaVib, CalcMeanFreePath
USE MOD_DSMC_Vars              ,ONLY: Coll_pData, CollisMode, ChemReac, PartStateIntEn, DSMC
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars              ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_Particle_Vars          ,ONLY: PEM, PDM, PartSpecies, PartState, Species, usevMPF, PartMPF, Species, PartPosRef
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_pos_and_velo      ,ONLY: SetParticleVelocity
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_ChemReact         ,ONLY: CalcPhotoIonizationNumber
USE MOD_DSMC_ChemReact         ,ONLY: PhotoIonization_InsertProducts
USE MOD_DSMC_AmbipolarDiffusion,ONLY: AD_DeleteParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec,iInit,TotalNbrOfReactions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iPair, iNewPart, iReac, ParticleIndex, NewParticleIndex, bgSpec, NbrOfParticle
REAL                          :: RandVal,NumTmp,ProbRest
INTEGER                       :: TotalNbrOfReactionsTmp,iCrossSection,NbrCrossSections
INTEGER                       :: NumPhotoIonization(ChemReac%NumOfReact)
REAL                          :: SumCrossSections
!===================================================================================================================================
NumPhotoIonization = 0

IF(TotalNbrOfReactions.LE.0) RETURN
TotalNbrOfReactionsTmp = TotalNbrOfReactions

! Ambipolar diffusion (up to four new particles can be produced during a single photo-ionization event)
newAmbiParts = 0
IF (ALLOCATED(iPartIndx_NodeNewAmbi)) DEALLOCATE(iPartIndx_NodeNewAmbi)
ALLOCATE(iPartIndx_NodeNewAmbi(4*TotalNbrOfReactions))

!> 0.) Calculate sum of cross-sections for photoionization
SumCrossSections = 0.
NbrCrossSections = 0
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  SumCrossSections = SumCrossSections + ChemReac%CrossSection(iReac)
  NbrCrossSections = NbrCrossSections + 1
END DO ! iReac = 1, ChemReac%NumOfReact

!> 1.) Compute the number of photoionization events in the local domain of each proc
iCrossSection = 0
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  iCrossSection  = iCrossSection + 1
  IF(iCrossSection.EQ.NbrCrossSections)THEN
    NumPhotoIonization(iReac) = TotalNbrOfReactionsTmp
    EXIT
  END IF ! iCrossSection.EQ.NbrCrossSections
  NumTmp = TotalNbrOfReactionsTmp*ChemReac%CrossSection(iReac)/SumCrossSections
  SumCrossSections = SumCrossSections - ChemReac%CrossSection(iReac)
  NumPhotoIonization(iReac) = INT(NumTmp)
  ProbRest = NumTmp - REAL(NumPhotoIonization(iReac))
  CALL RANDOM_NUMBER(RandVal)
  IF (ProbRest.GT.RandVal) NumPhotoIonization(iReac) = NumPhotoIonization(iReac) + 1
  TotalNbrOfReactionsTmp = TotalNbrOfReactionsTmp - NumPhotoIonization(iReac)
END DO

!> 2.) Delete left-over inserted particles
IF(TotalNbrOfReactions.GT.SUM(NumPhotoIonization)) THEN
  DO iPart = SUM(NumPhotoIonization)+1,TotalNbrOfReactions
    PDM%ParticleInside(PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)) = .FALSE.
  END DO
ELSE IF(TotalNbrOfReactions.LT.SUM(NumPhotoIonization)) THEN
  CALL Abort(&
    __STAMP__&
    ,'ERROR in PhotoIonization: Something is wrong, trying to perform more reactions than anticipated!')
END IF

IF(SUM(NumPhotoIonization).EQ.0) RETURN

!> 3.) Insert the products of the photoionization rections
NbrOfParticle = SUM(NumPhotoIonization)

ALLOCATE(Coll_pData(NbrOfParticle))
Coll_pData%Ec=0.
DSMCSumOfFormedParticles = 0

iNewPart = 0; iPair = 0

DO iPart = 1, NbrOfParticle
  ! Loop over the particles with a set position (from SetParticlePosition)
  ParticleIndex = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = ParticleIndex
  END IF
  iNewPart = iNewPart + 1
  ! Get a new index for the second product
  NewParticleIndex = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition+NbrOfParticle)
  IF (NewParticleIndex.EQ.0) THEN
    CALL Abort(&
      __STAMP__&
      ,'ERROR in PhotoIonization: MaxParticleNumber should be increased!')
  END IF
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = NewParticleIndex
  END IF
  PartState(1:3,NewParticleIndex) = PartState(1:3,ParticleIndex)
  IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
    PartPosRef(1:3,NewParticleIndex)=PartPosRef(1:3,ParticleIndex)
  END IF
  ! Species index given from the initialization
  PartSpecies(ParticleIndex) = iSpec
  ! Get the species index of the background gas
  bgSpec = BGGas_GetSpecies()
  PartSpecies(NewParticleIndex) = bgSpec
  IF(CollisMode.GT.1) THEN
    IF(SpecDSMC(bgSpec)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(bgSpec,1,NewParticleIndex,1)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(bgSpec,1,NewParticleIndex,1)
    END IF
  END IF
  CALL CalcVelocity_maxwell_lpn(FractNbr=bgSpec, Vec3D=PartState(4:6,NewParticleIndex), iInit=1)
  ! Particle flags
  PDM%ParticleInside(NewParticleIndex)  = .TRUE.
  PDM%IsNewPart(NewParticleIndex)       = .TRUE.
  PDM%dtFracPush(NewParticleIndex)      = .FALSE.
  ! Particle element
  PEM%GlobalElemID(NewParticleIndex) = PEM%GlobalElemID(ParticleIndex)
  ! Last element ID
  PEM%LastGlobalElemID(NewParticleIndex) = PEM%GlobalElemID(NewParticleIndex)
  PEM%LastGlobalElemID(ParticleIndex) = PEM%GlobalElemID(ParticleIndex)
  ! Pairing (first particle is the background gas species)
  Coll_pData(iNewPart)%iPart_p1 = NewParticleIndex
  Coll_pData(iNewPart)%iPart_p2 = ParticleIndex
  ! Relative velocity is not required as the relative translational energy will not be considered
  Coll_pData(iNewPart)%CRela2 = 0.
  ! Weighting factor
  IF(usevMPF) THEN
    PartMPF(NewParticleIndex) = Species(bgSpec)%MacroParticleFactor
    PartMPF(ParticleIndex) = Species(iSpec)%MacroParticleFactor
  END IF
  ! Velocity (set it to zero, as it will be substracted in the chemistry module)
  PartState(4:6,ParticleIndex) = 0.
  ! Internal energies (set it to zero)
  PartStateIntEn(1:2,ParticleIndex) = 0.
  IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,ParticleIndex) = 0.
END DO

! Add the particles initialized through the emission and the background particles
PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle + iNewPart
! Update the current next free position
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle + iNewPart

IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) THEN
  CALL Abort(&
    __STAMP__&
    ,'ERROR in PhotoIonization: ParticleVecLength greater than MaxParticleNumber! Increase the MaxParticleNumber to at least: ' &
    , IntInfoOpt=PDM%ParticleVecLength)
END IF

!> 4.) Perform the reaction, distribute the collision energy (including photon energy) and emit electrons perpendicular
!>     to the photon's path
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  DO iPart = 1, NumPhotoIonization(iReac)
    iPair = iPair + 1
    CALL PhotoIonization_InsertProducts(iPair, iReac, iInit, iSpec)
  END DO
END DO

! Advance particle vector length and the current next free position with newly created particles
PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles

DSMCSumOfFormedParticles = 0

DEALLOCATE(Coll_pData)

IF (DSMC%DoAmbipolarDiff) THEN
  ! Every particle created during photo-ionization is stored in the iPartIndx_NodeNewAmbi array, which is treated in the routine
  CALL AD_DeleteParticles()
END IF

END SUBROUTINE BGGas_PhotoIonization

END MODULE MOD_DSMC_BGGas
