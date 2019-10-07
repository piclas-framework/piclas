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

INTERFACE DSMC_InitBGGas
  MODULE PROCEDURE DSMC_InitBGGas
END INTERFACE

INTERFACE DSMC_pairing_bggas
  MODULE PROCEDURE DSMC_pairing_bggas
END INTERFACE

INTERFACE DSMC_FinalizeBGGas
  MODULE PROCEDURE DSMC_FinalizeBGGas
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_InitBGGas, DSMC_pairing_bggas, MCC_pairing_bggas, DSMC_FinalizeBGGas
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_InitBGGas()
!===================================================================================================================================
!> Initialization of background gas (BGG): Creating particles of the background species for each actual simulation particle
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
USE MOD_DSMC_Init              ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_Vars              ,ONLY: BGGas, SpecDSMC
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_PARTICLE_Vars          ,ONLY: PDM, PartSpecies, PartState, PEM, PartPosRef
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF
USE MOD_part_pos_and_velo      ,ONLY: SetParticleVelocity
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefmapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr
!===================================================================================================================================
iNewPart=0
PositionNbr = 0
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart).AND.(PartSpecies(iPart).NE.BGGas%BGGasSpecies)) THEN
    iNewPart = iNewPart + 1
    PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
    IF (PositionNbr.EQ.0) THEN
      CALL Abort(&
__STAMP__&
,'ERROR in BGGas: MaxParticleNumber should be twice the expected number of particles, to account for the BGG particles!')
    END IF
    PartState(PositionNbr,1:3) = PartState(iPart,1:3)
    IF(DoRefMapping)THEN ! here Nearst-GP is missing
      PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
    END IF
    PartSpecies(PositionNbr) = BGGas%BGGasSpecies
    IF(SpecDSMC(BGGas%BGGasSpecies)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(BGGas%BGGasSpecies,0,PositionNbr,1)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(BGGas%BGGasSpecies,0,PositionNbr,1)
    END IF
    PEM%Element(PositionNbr) = PEM%Element(iPart)
    PDM%ParticleInside(PositionNbr) = .true.
    PEM%pNext(PEM%pEnd(PEM%Element(PositionNbr))) = PositionNbr     ! Next Particle of same Elem (Linked List)
    PEM%pEnd(PEM%Element(PositionNbr)) = PositionNbr
    PEM%pNumber(PEM%Element(PositionNbr)) = &                       ! Number of Particles in Element
    PEM%pNumber(PEM%Element(PositionNbr)) + 1
    BGGas%PairingPartner(iPart) = PositionNbr
  END IF
END DO
CALL SetParticleVelocity(BGGas%BGGasSpecies,0,iNewPart,1) ! Properties of BG gas are stored in iInit=0
PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart

END SUBROUTINE DSMC_InitBGGas


SUBROUTINE DSMC_pairing_bggas(iElem)
!===================================================================================================================================
! Building of pairs for the background gas
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze        ,ONLY: CalcGammaVib
USE MOD_DSMC_Vars           ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SelectionProc
USE MOD_DSMC_Vars           ,ONLY: VarVibRelaxProb
USE MOD_Particle_Vars       ,ONLY: PEM,PartSpecies,nSpecies,PartState,Species,usevMPF,PartMPF,Species
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: nPair, iPair, iPart, iLoop, nPart, iSpec
  INTEGER                       :: cSpec1, cSpec2, iCase
!===================================================================================================================================
  nPart = PEM%pNumber(iElem)
  nPair = INT(nPart/2)
  IF(DSMC%VibRelaxProb.EQ.2.0) THEN ! Set summs for variable vibrational relaxation to zero
    DO iSpec=1,nSpecies
      VarVibRelaxProb%ProbVibAvNew(iSpec) = 0
      VarVibRelaxProb%nCollis(iSpec) = 0
    END DO
  END IF

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
    CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + 1
    ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(iPart,1)
    ! Creating pairs for species, which are not the background species nor treated with a MCC model
    IF(iSpec.NE.BGGas%BGGasSpecies) THEN
      Coll_pData(iPair)%iPart_p1 = iPart
      Coll_pData(iPair)%iPart_p2 = BGGas%PairingPartner(iPart)
      iPair = iPair + 1
    END IF
    iPart = PEM%pNext(iPart)
  END DO

  ! Setting Number of BGGas Particles per Cell
  IF (Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensityFileID.GT.0) THEN
    BGGas%BGColl_SpecPartNum = Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensity(iElem) * GEO%Volume(iElem)      &
                                               / Species(BGGas%BGGasSpecies)%MacroParticleFactor
  ELSE
    BGGas%BGColl_SpecPartNum = BGGas%BGGasDensity * GEO%Volume(iElem)      &
                                               / Species(BGGas%BGGasSpecies)%MacroParticleFactor
  END IF

  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
    ! 3. Case: Temperature required for the mean free path with the VHS model
    ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
    ! an impact on the temperature, a background gas should not be utilized in the first place.
    DSMC%InstantTransTemp(nSpecies+1) = Species(BGGas%BGGasSpecies)%Init(0)%MWTemperatureIC
    IF(SelectionProc.EQ.2) CALL CalcGammaVib()
  END IF

  CollInf%Coll_SpecPartNum(BGGas%BGGasSpecies) = BGGas%BGColl_SpecPartNum

  DO iPair = 1, nPair
    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
    IF (usevMPF) PartMPF(Coll_pData(iPair)%iPart_p2) = PartMPF(Coll_pData(iPair)%iPart_p1)
    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
    Coll_pData(iPair)%PairType = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
  END DO

END SUBROUTINE DSMC_pairing_bggas


SUBROUTINE MCC_pairing_bggas(iElem)
!===================================================================================================================================
! Building of pairs for the background gas
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze            ,ONLY: CalcGammaVib
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SelectionProc
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, SpecMCC, MCC
USE MOD_Particle_Vars           ,ONLY: PEM,PDM,PartSpecies,nSpecies,PartState,Species,usevMPF,PartMPF,Species,PartPosRef
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_DSMC_Init               ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_emission_tools     ,ONLY: SetParticleChargeAndMass,SetParticleMPF
USE MOD_part_pos_and_velo       ,ONLY: SetParticleVelocity
USE MOD_part_tools              ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefmapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPair, iPart, iLoop, nPart, iSpec, iNewPart, PositionNbr
INTEGER                       :: cSpec1, cSpec2, iCase, SpecPairNum(nSpecies), SpecPairNumCounter(nSpecies)
INTEGER,ALLOCATABLE           :: iPartIndex(:), PairingPartner(:)
!===================================================================================================================================
nPart = PEM%pNumber(iElem)
MCC%TotalPairNum = 0.

ALLOCATE(iPartIndex(nPart))
CollInf%Coll_SpecPartNum = 0
CollInf%Coll_CaseNum = 0

SpecPairNum = 0
SpecPairNumCounter = 0

IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  iSpec = PartSpecies(iPart)
  ! Counting the number of particles per species
  CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + 1
  ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(iPart,1)
  ! Create particle index list for pairing
  iPartIndex(iLoop) = iPart
  ! Counting the number of pairs for regular DSMC background gas
  IF(.NOT.SpecMCC(iSpec)%UseCollXSec) THEN
    MCC%TotalPairNum = MCC%TotalPairNum + 1
  END IF
  iPart = PEM%pNext(iPart)
END DO

! Determining the number of pairs for MCC
DO iSpec = 1,nSpecies
  IF(SpecMCC(iSpec)%UseCollXSec) THEN
    SpecPairNum(iSpec) = NINT(CollInf%Coll_SpecPartNum(iSpec)*SpecMCC(iSpec)%ProbNull)
    MCC%TotalPairNum = MCC%TotalPairNum + SpecPairNum(iSpec)
  END IF
END DO

ALLOCATE(Coll_pData(MCC%TotalPairNum))
Coll_pData%Ec = 0.
ALLOCATE(PairingPartner(MCC%TotalPairNum))
PairingPartner = 0
! 
iNewPart=0
PositionNbr = 0
DO iLoop = 1, MCC%TotalPairNum
  ! Taking a random particle from the cell to get the position of the new particle
  iPart = iPartIndex(iLoop)
  ! Creating a new background gas particle
  iNewPart = iNewPart + 1
  PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
  IF (PositionNbr.EQ.0) THEN
    CALL Abort(&
__STAMP__&
,'ERROR in MCC: MaxParticleNumber should be twice the expected number of particles, to account for the BGG/MCC particles!')
  END IF
  PartState(PositionNbr,1:3) = PartState(iPart,1:3)
  IF(DoRefMapping)THEN ! here Nearst-GP is missing
    PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
  END IF
  PartSpecies(PositionNbr) = BGGas%BGGasSpecies
  IF(SpecDSMC(BGGas%BGGasSpecies)%PolyatomicMol) THEN
    CALL DSMC_SetInternalEnr_Poly(BGGas%BGGasSpecies,0,PositionNbr,1)
  ELSE
    CALL DSMC_SetInternalEnr_LauxVFD(BGGas%BGGasSpecies,0,PositionNbr,1)
  END IF
  PEM%Element(PositionNbr) = iElem
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PEM%pNext(PEM%pEnd(iElem)) = PositionNbr     ! Next Particle of same Elem (Linked List)
  PEM%pEnd(iElem) = PositionNbr
  PEM%pNumber(iElem) = PEM%pNumber(iElem) + 1
  ! Saving the particle index to later
  PairingPartner(iLoop) = PositionNbr
END DO

CALL SetParticleVelocity(BGGas%BGGasSpecies,0,iNewPart,1) ! Properties of BG gas are stored in iInit=0
PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart

iPair = 0
DO iLoop = 1, nPart
  iPart = iPartIndex(iLoop)
  iSpec = PartSpecies(iPart)
  IF(SpecMCC(iSpec)%UseCollXSec) THEN
    IF(SpecPairNumCounter(iSpec).GE.SpecPairNum(iSpec)) CYCLE
    SpecPairNumCounter(iSpec) = SpecPairNumCounter(iSpec) + 1
  END IF
  iPair = iPair + 1
  Coll_pData(iPair)%iPart_p1 = iPart
  Coll_pData(iPair)%iPart_p2 = PairingPartner(iPair)
END DO  

! Setting Number of BGGas Particles per Cell
IF (Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensityFileID.GT.0) THEN
  BGGas%BGColl_SpecPartNum = Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensity(iElem) * GEO%Volume(iElem)      &
                                              / Species(BGGas%BGGasSpecies)%MacroParticleFactor
ELSE
  BGGas%BGColl_SpecPartNum = BGGas%BGGasDensity * GEO%Volume(iElem)      &
                                              / Species(BGGas%BGGasSpecies)%MacroParticleFactor
END IF

IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
  ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
  !           relaxation probability (CalcGammaVib)
  ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
  ! 3. Case: Temperature required for the mean free path with the VHS model
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = Species(BGGas%BGGasSpecies)%Init(0)%MWTemperatureIC
  IF(SelectionProc.EQ.2) CALL CalcGammaVib()
END IF

CollInf%Coll_SpecPartNum(BGGas%BGGasSpecies) = BGGas%BGColl_SpecPartNum

DO iPair = 1, MCC%TotalPairNum
  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
  IF (usevMPF) PartMPF(Coll_pData(iPair)%iPart_p2) = PartMPF(Coll_pData(iPair)%iPart_p1)
  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
  Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                            -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                            + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                            -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                            + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                            -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
  Coll_pData(iPair)%PairType = iCase
  Coll_pData(iPair)%NeedForRec = .FALSE.
END DO

END SUBROUTINE MCC_pairing_bggas


SUBROUTINE DSMC_FinalizeBGGas()
!===================================================================================================================================
! Kills all BG Particles
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  WHERE (PartSpecies(1:PDM%ParticleVecLength).EQ.BGGas%BGGasSpecies) &
         PDM%ParticleInside(1:PDM%ParticleVecLength) = .FALSE.
  BGGas%PairingPartner = 0
  CALL UpdateNextFreePosition()

END SUBROUTINE DSMC_FinalizeBGGas

END MODULE MOD_DSMC_BGGas
