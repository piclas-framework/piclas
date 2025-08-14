!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_DSMC_Symmetry
!===================================================================================================================================
!> Routines for 2D (planar/axisymmetric) simulations
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
PUBLIC :: AdjustParticleWeight, SetInClones, DSMC_TreatIdenticalParticles
!===================================================================================================================================

CONTAINS

SUBROUTINE AdjustParticleWeight(iPart,iElem)
!===================================================================================================================================
!> Routine for the treatment of particles with enabled radial/linear weighting (weighting factor is increasing linearly with
!> increasing y/along user-defined vector)
!> 1.) Determine the new particle weight and decide whether to clone or to delete the particle
!> 2a.) Particle cloning, if the local weighting factor is smaller than the previous
!> 2b.) Particle deletion, if the local weighting factor is greater than the previous
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_DSMC_Vars               ,ONLY: ParticleWeighting, DSMC, PartStateIntEn, useDSMC, CollisMode, AmbipolElecVelo
USE MOD_DSMC_Vars               ,ONLY: ClonedParticles, VibQuantsPar, SpecDSMC, PolyatomMolDSMC, ElectronicDistriPart
USE MOD_DSMC_Vars               ,ONLY: DoRadialWeighting
USE MOD_Particle_Vars           ,ONLY: PartMPF, PartSpecies, PartState, Species, LastPartPos
USE MOD_TimeDisc_Vars           ,ONLY: iter
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF, CalcVarWeightMPF
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SpecID, iPolyatMole, cloneIndex, DelayCounter
REAL                            :: DeleteProb, iRan, NewMPF, CloneProb, OldMPF
LOGICAL                         :: DoCloning
!===================================================================================================================================
DoCloning = .FALSE.
DeleteProb = 0.
SpecID = PartSpecies(iPart)

! 1.) Determine the new particle weight and decide whether to clone or to delete the particle
IF(DoRadialWeighting) THEN
  IF (.NOT.(PartMPF(iPart).GT.Species(SpecID)%MacroParticleFactor)) RETURN
  NewMPF = CalcRadWeightMPF(PartState(2,iPart),SpecID,iPart)
ELSE
  NewMPF = CalcVarWeightMPF(PartState(:,iPart),(iElem-offSetElem),iPart)
END IF
OldMPF = PartMPF(iPart)
CloneProb = (OldMPF/NewMPF)-INT(OldMPF/NewMPF)
CALL RANDOM_NUMBER(iRan)
IF((CloneProb.GT.iRan).AND.(NewMPF.LT.OldMPF)) THEN
  DoCloning = .TRUE.
  IF(INT(OldMPF/NewMPF).GT.1) THEN
    IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
    CALL Abort(__STAMP__,&
      'ERROR in particle weighting: More than one clone per particle is not allowed! Reduce the time step or'//&
        ' the radial/linear weighting factor! Cloning probability is:',RealInfoOpt=CloneProb)
  END IF
END IF
PartMPF(iPart) = NewMPF

IF(DoCloning) THEN
  ! 2a.) Particle cloning, if the local weighting factor is smaller than the previous (particle travelling downwards)
  ! Get the list number to store the clones, depending on the chosen clone mode
  SELECT CASE(ParticleWeighting%CloneMode)
  CASE(1)
  ! ######## Clone Delay ###################################################################################################
  ! Insertion of the clones after a defined delay, all clones are collected in a single list and inserted before boundary
  ! treatment in the next time step at their original positions
    DelayCounter = MOD((INT(iter,4)+ParticleWeighting%CloneDelayDiff-1),ParticleWeighting%CloneInputDelay)
  CASE(2)
  ! ######## Clone Random Delay #############################################################################################
  ! A list, which is ParticleWeighting%CloneInputDelay + 1 long, is filled with clones to be inserted. After the list
  ! is full, NextClone gives the empty particle list, whose clones were inserted during the last SetInClones step
    IF((INT(iter,4)+ParticleWeighting%CloneDelayDiff).LE.ParticleWeighting%CloneInputDelay) THEN
      DelayCounter = INT(iter,4)+ParticleWeighting%CloneDelayDiff
    ELSE
      DelayCounter = ParticleWeighting%NextClone
    END IF
  END SELECT
  ! Storing the particle information
  IF(ParticleWeighting%ClonePartNum(DelayCounter)+1.GT.ParticleWeighting%CloneVecLength) CALL IncreaseClonedParticlesType()
  ParticleWeighting%ClonePartNum(DelayCounter) = ParticleWeighting%ClonePartNum(DelayCounter) + 1
  cloneIndex = ParticleWeighting%ClonePartNum(DelayCounter)
  ClonedParticles(cloneIndex,DelayCounter)%PartState(1:6)= PartState(1:6,iPart)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(1:2) = PartStateIntEn(1:2,iPart)
    IF(DSMC%ElectronicModel.GT.0) THEN
      ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(3) =   PartStateIntEn(3,iPart)
      IF ((DSMC%ElectronicModel.EQ.2).AND.(.NOT.((Species(SpecID)%InterID.EQ.4).OR.SpecDSMC(SpecID)%FullyIonized))) THEN
        IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc)) &
          DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc)
        ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc(1:SpecDSMC(SpecID)%MaxElecQuant))
        ClonedParticles(cloneIndex,DelayCounter)%DistriFunc(:) = ElectronicDistriPart(iPart)%DistriFunc(:)
      END IF
    END IF
    IF ((DSMC%DoAmbipolarDiff).AND.(Species(SpecID)%ChargeIC.GT.0.0)) THEN
      IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo)) &
        DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo)
      ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo(1:3))
      ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo(1:3) = AmbipolElecVelo(iPart)%ElecVelo(1:3)
    END IF
    IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
      IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)) &
        DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)
      ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      ClonedParticles(cloneIndex,DelayCounter)%VibQuants(:) = VibQuantsPar(iPart)%Quants(:)
    END IF
  END IF
  ClonedParticles(cloneIndex,DelayCounter)%Species = SpecID
  ClonedParticles(cloneIndex,DelayCounter)%Element = iElem
  ClonedParticles(cloneIndex,DelayCounter)%LastPartPos(1:3) = LastPartPos(1:3,iPart)
  ClonedParticles(cloneIndex,DelayCounter)%WeightingFactor = PartMPF(iPart)
ELSE
! ######## Particle Delete #######################################################################################################
! 2b.) Particle deletion, if the local weighting factor is greater than the previous (particle travelling upwards)
  IF(NewMPF.GT.OldMPF) THEN
    ! Start deleting particles after the clone delay has passed and particles are also inserted
    IF((INT(iter,4)+ParticleWeighting%CloneDelayDiff).LE.ParticleWeighting%CloneInputDelay) RETURN
    DeleteProb = 1. - CloneProb
    IF (DeleteProb.GT.0.5) THEN
      IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
      CALL abort(__STAMP__,&
        'ERROR in Radial Weighting of 2D/Axisymmetric: The deletion probability is higher than 0.5! Reduce the time step or'//&
        ' the radial weighting factor! Deletion probability is:',RealInfoOpt=DeleteProb)
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF(DeleteProb.GT.iRan) THEN
      CALL RemoveParticle(iPart)
    END IF
  END IF
END IF

END SUBROUTINE AdjustParticleWeight


SUBROUTINE SetInClones()
!===================================================================================================================================
!> Insertion of cloned particles during the previous time steps. Clones insertion is delayed by at least one time step to avoid the
!> avalanche phenomenon (identical particles travelling on the same path, not colliding due to zero relative velocity).
!> 1.) Chose which list to insert depending on the clone mode
!> 2.) Insert the clones at the position they were created
!> 3.) Reset the list
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: ClonedParticles, PartStateIntEn, useDSMC, CollisMode, DSMC, ParticleWeighting
USE MOD_DSMC_Vars               ,ONLY: AmbipolElecVelo, DoRadialWeighting
USE MOD_DSMC_Vars               ,ONLY: VibQuantsPar, SpecDSMC, PolyatomMolDSMC, SamplingActive, ElectronicDistriPart
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartSpecies, PartState, LastPartPos, PartMPF, WriteMacroVolumeValues, Species
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: iter
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart, PositionNbr, iPolyatMole, DelayCounter, locElemID
REAL                            :: iRan
!===================================================================================================================================

! 1.) Chose which list to insert depending on the clone mode
SELECT CASE(ParticleWeighting%CloneMode)
CASE(1)
  ! During the first iterations the delay counter refers to the empty clone array (which is filled during the following tracking)
  ! Afterwards, the MODULUS counts up from zero to CloneInputDelay-1
  DelayCounter = MOD((INT(iter,4)+ParticleWeighting%CloneDelayDiff-1),ParticleWeighting%CloneInputDelay)
CASE(2)
  ! During the first iterations, check if number of iterations is less than the input delay and leave routine. Afterwards, a
  ! random clone list from the previous time steps is chosen.
  IF((INT(iter,4)+ParticleWeighting%CloneDelayDiff).GT.ParticleWeighting%CloneInputDelay) THEN
    CALL RANDOM_NUMBER(iRan)
    ! Choosing random clone between 0 and CloneInputDelay
    DelayCounter = INT((ParticleWeighting%CloneInputDelay+1)*iRan)
    DO WHILE (DelayCounter.EQ.ParticleWeighting%NextClone)
      CALL RANDOM_NUMBER(iRan)
      DelayCounter = INT((ParticleWeighting%CloneInputDelay+1)*iRan)
    END DO
    ! Save the chosen list as the next available list to store clones in the next time step
    ParticleWeighting%NextClone = DelayCounter
  ELSE
    RETURN
  END IF
END SELECT

IF(ParticleWeighting%ClonePartNum(DelayCounter).EQ.0) RETURN

! 2.) Insert the clones at the position they were created
DO iPart = 1, ParticleWeighting%ClonePartNum(DelayCounter)
  PositionNbr = GetNextFreePosition()
  ! Copy particle parameters
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PDM%IsNewPart(PositionNbr) = .TRUE.
  PDM%dtFracPush(PositionNbr) = .FALSE.
  PartState(1:5,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartState(1:5)
  IF (DoRadialWeighting) THEN
    ! Creating a relative velocity in the z-direction
    PartState(6,PositionNbr) = -ClonedParticles(iPart,DelayCounter)%PartState(6)
  ELSE
    PartState(6,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartState(6)
  END IF
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    PartStateIntEn(1:2,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(1:2)
    IF(DSMC%ElectronicModel.GT.0) THEN
      PartStateIntEn(3,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(3)
      IF ((DSMC%ElectronicModel.EQ.2).AND.(.NOT.((Species(ClonedParticles(iPart,DelayCounter)%Species)%InterID.EQ.4) &
          .OR.SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%FullyIonized))) THEN
        IF(ALLOCATED(ElectronicDistriPart(PositionNbr)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(PositionNbr)%DistriFunc)
        ALLOCATE(ElectronicDistriPart(PositionNbr)%DistriFunc(1:SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%MaxElecQuant))
        ElectronicDistriPart(PositionNbr)%DistriFunc(:) = ClonedParticles(iPart,DelayCounter)%DistriFunc(:)
      END IF
    END IF
    IF ((DSMC%DoAmbipolarDiff).AND.(Species(ClonedParticles(iPart,DelayCounter)%Species)%ChargeIC.GT.0.0)) THEN
      IF(ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(1:3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:2) = ClonedParticles(iPart,DelayCounter)%AmbiPolVelo(1:2)
      IF(DoRadialWeighting) THEN
        ! Creating a relative velocity in the z-direction
        AmbipolElecVelo(PositionNbr)%ElecVelo(3) = -ClonedParticles(iPart,DelayCounter)%AmbiPolVelo(3)
      ELSE
        AmbipolElecVelo(PositionNbr)%ElecVelo(3) = ClonedParticles(iPart,DelayCounter)%AmbiPolVelo(3)
      END IF
    END IF
    IF(SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%SpecToPolyArray
      IF(ALLOCATED(VibQuantsPar(PositionNbr)%Quants)) DEALLOCATE(VibQuantsPar(PositionNbr)%Quants)
      ALLOCATE(VibQuantsPar(PositionNbr)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      VibQuantsPar(PositionNbr)%Quants(:) = ClonedParticles(iPart,DelayCounter)%VibQuants(:)
    END IF
  END IF
  PartSpecies(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Species
  ! Set the global element number with the offset
  PEM%GlobalElemID(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Element
  PEM%LastGlobalElemID(PositionNbr) = 0
  locElemID = PEM%LocalElemID(PositionNbr)
  LastPartPos(1:3,PositionNbr) = PartState(1:3,PositionNbr)
  PartMPF(PositionNbr) =  ClonedParticles(iPart,DelayCounter)%WeightingFactor
  IF (UseVarTimeStep) THEN
    PartTimeStep(PositionNbr) = GetParticleTimeStep(PartState(1,PositionNbr),PartState(2,PositionNbr),locElemID)
  END IF
  ! Counting the number of clones per cell
  IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
    IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(locElemID,5) = DSMC%QualityFacSamp(locElemID,5) + 1
  END IF
  IF(CalcPartBalance) THEN
    nPartIn(PartSpecies(PositionNbr))=nPartIn(PartSpecies(PositionNbr)) + 1
  END IF ! CalcPartBalance
END DO

! 3.) Reset the list
ParticleWeighting%ClonePartNum(DelayCounter) = 0

! 3.1) Reduce ClonedParticles if necessary
CALL ReduceClonedParticlesType()

END SUBROUTINE SetInClones


SUBROUTINE DSMC_TreatIdenticalParticles(iPair, nPair, nPart, iElem, iPartIndx_Node)
!===================================================================================================================================
!> Check if particle pairs have a zero relative velocity (and thus a collision probability of zero), if they do, break up the pair
!> and use either a left-over particle (uneven number of particles in a cell) or swap the collision partners with the next pair in
!> the list.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, DSMC, ChemReac, SamplingActive, CollInf, CollisMode
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartState, WriteMacroVolumeValues
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(IN)           :: nPair
INTEGER, INTENT(IN)           :: nPart
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)        :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart_p1, iPart_p2, tempPart, cSpec1, cSpec2, iCase
!===================================================================================================================================
! Two particles with the exact same velocities at the same positions -> clones that did not interact with other particles/walls
IF (Coll_pData(iPair)%CRela2.EQ.0.0) THEN
  IF ((CollisMode.LT.3).AND.(nPart.EQ.1)) THEN
    ! Uneven number of particles in the cell, a single particle is left without a pair
    ! Removing the pairs from the weighting factor and the case num sums
    CollInf%SumPairMPF(Coll_pData(iPair)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair)%PairType) &
      -(GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
    CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
    ! Swapping particle without a pair with the first particle of the current pair
    tempPart = Coll_pData(iPair)%iPart_p1
    Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(1)
    iPartIndx_Node(1) = tempPart
    IF (CollisMode.EQ.3) ChemReac%RecombParticle = iPartIndx_Node(1)
    IF (CollInf%ProhibitDoubleColl)  CollInf%OldCollPartner(iPartIndx_Node(1)) = 0
    ! Increase the appropriate case number and set the right pair type
    iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
    cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
    Coll_pData(iPair)%PairType = iCase
    ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
    CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))*0.5
    ! Calculation of the relative velocity for the new first pair
    Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                             + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                             + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
  ELSE IF (iPair.LT.nPair) THEN
    IF (.NOT.Coll_pData(iPair+1)%NeedForRec) THEN
    ! "Partner-Tausch": if there are pairs ahead in the pairing list, the next is pair is broken up and collision partners
    ! are swapped
      CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
      CollInf%SumPairMPF(Coll_pData(iPair)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair)%PairType) &
        - 0.5 * (GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))
      CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) - 1
      CollInf%SumPairMPF(Coll_pData(iPair+1)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair+1)%PairType) &
        - 0.5 * (GetParticleWeight(Coll_pData(iPair+1)%iPart_p1) + GetParticleWeight(Coll_pData(iPair+1)%iPart_p2))
      ! Breaking up the next pair and swapping partners
      tempPart = Coll_pData(iPair)%iPart_p1
      Coll_pData(iPair)%iPart_p1 = Coll_pData(iPair + 1)%iPart_p1
      Coll_pData(iPair + 1)%iPart_p1 = tempPart
      ! Increase the appropriate case number and set the right pair type
      iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
      cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
      iCase = CollInf%Coll_Case(cSpec1, cSpec2)
      CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      Coll_pData(iPair)%PairType = iCase
      ! Calculation of the relative velocity for the new first pair
      Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                               + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                               + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      IF(Coll_pData(iPair)%CRela2.EQ.0.0) THEN
        ! If the relative velocity is still zero, add the pair to the identical particles count
        IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
          IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
        END IF
      END IF
      ! Increase the appropriate case number and set the right pair type
      iPart_p1 = Coll_pData(iPair+1)%iPart_p1; iPart_p2 = Coll_pData(iPair+1)%iPart_p2
      cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
      iCase = CollInf%Coll_Case(cSpec1, cSpec2)
      CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      ! Calculation of the relative velocity for the new follow-up pair
      Coll_pData(iPair+1)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                                 + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                                 + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      Coll_pData(iPair+1)%PairType = iCase
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
    ELSE
      ! For the last pair, only invert the velocity in z and calculate new relative velocity
      iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
      PartState(6,iPart_p1) = - PartState(6,iPart_p1)
      Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                               + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                               + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      ! Add the pair to the number of identical particles output
      IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
        IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
      END IF
    END IF  ! nPart.EQ.1/iPair.LT.nPair
  ELSE
    ! For the last pair, only invert the velocity in z and calculate new relative velocity
    iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
    PartState(6,iPart_p1) = - PartState(6,iPart_p1)
    Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                             + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                             + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
    ! Add the pair to the number of identical particles output
    IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
      IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
    END IF
  END IF  ! nPart.EQ.1/iPair.LT.nPair
END IF    ! Coll_pData(iPair)%CRela2.EQ.0.0

END SUBROUTINE DSMC_TreatIdenticalParticles


SUBROUTINE IncreaseClonedParticlesType()
!===================================================================================================================================
!> Increases CloneVecLength and the ClonedParticles(iPart,iDelay) type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: ParticleWeighting, ClonedParticles, tClonedParticles
USE MOD_Particle_Vars         ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: NewSize,i,ii,ALLOCSTAT
TYPE (tClonedParticles), ALLOCATABLE :: ClonedParticles_new(:,:)
!===================================================================================================================================
NewSize = MAX(CEILING(ParticleWeighting%CloneVecLength * (1+PDM%MaxPartNumIncrease)),ParticleWeighting%CloneVecLength+ParticleWeighting%CloneVecLengthDelta)

SELECT CASE(ParticleWeighting%CloneMode)
CASE(1)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:(ParticleWeighting%CloneInputDelay-1)),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,ParticleWeighting%CloneInputDelay-1
    DO i=1,ParticleWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE(2)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:ParticleWeighting%CloneInputDelay),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,ParticleWeighting%CloneInputDelay
    DO i=1,ParticleWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE DEFAULT
  CALL Abort(__STAMP__,'ERROR in IncreaseClonedParticlesType: The selected cloning mode is not available! Choose between 1 and 2.'//&
                      ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

ParticleWeighting%CloneVecLength = NewSize

END SUBROUTINE IncreaseClonedParticlesType


SUBROUTINE ReduceClonedParticlesType()
!===================================================================================================================================
!> Reduces ParticleWeighting%CloneVecLength and the ClonedParticles(iPart,iDelay) type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: ParticleWeighting, ClonedParticles, tClonedParticles
USE MOD_Particle_Vars         ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: NewSize,i,ii,ALLOCSTAT
TYPE (tClonedParticles), ALLOCATABLE :: ClonedParticles_new(:,:)
!===================================================================================================================================
IF (MAXVAL(ParticleWeighting%ClonePartNum(:)).GE.PDM%maxParticleNumber/(1.+PDM%MaxPartNumIncrease)**2) RETURN

NewSize = MAX(CEILING(MAXVAL(ParticleWeighting%ClonePartNum(:))*(1.+PDM%MaxPartNumIncrease)),1)

IF (NewSize.GT.ParticleWeighting%CloneVecLength-ParticleWeighting%CloneVecLengthDelta) RETURN

SELECT CASE(ParticleWeighting%CloneMode)
CASE(1)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:(ParticleWeighting%CloneInputDelay-1)),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,ParticleWeighting%CloneInputDelay-1
    DO i=1,ParticleWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE(2)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:ParticleWeighting%CloneInputDelay),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,ParticleWeighting%CloneInputDelay
    DO i=1,ParticleWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE DEFAULT
  CALL Abort(__STAMP__,'ERROR in ReduceClonedParticlesType: The selected cloning mode is not available! Choose between 1 and 2.'//&
                      ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

ParticleWeighting%CloneVecLength = NewSize

END SUBROUTINE ReduceClonedParticlesType

END MODULE MOD_DSMC_Symmetry
