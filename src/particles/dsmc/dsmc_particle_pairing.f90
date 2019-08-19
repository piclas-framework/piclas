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

MODULE MOD_DSMC_ParticlePairing
!===================================================================================================================================
! Module including different particle pairing algorithms
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_pairing_octree
  MODULE PROCEDURE DSMC_pairing_octree
END INTERFACE

INTERFACE DSMC_pairing_statistical
  MODULE PROCEDURE DSMC_pairing_statistical
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_pairing_statistical, DSMC_pairing_octree, DSMC_init_octree, DSMC_pairing_quadtree, DSMC_CalcSubNodeVolumes
PUBLIC :: DSMC_CalcSubNodeVolumes2D, GeoCoordToMap2D
!===================================================================================================================================

CONTAINS


SUBROUTINE FindNearestNeigh(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : CollInf, tTreeNode, CollisMode, ChemReac, PartStateIntEn, Coll_pData, SelectionProc
  USE MOD_DSMC_Vars,              ONLY : DSMC, PairE_vMPF, SpecDSMC, VarVibRelaxProb
  USE MOD_Particle_Vars,          ONLY : PartState, nSpecies, PartSpecies, usevMPF, PartMPF, WriteMacroVolumeValues, VarTimeStep
  USE MOD_DSMC_Relaxation,        ONLY : SetMeanVibQua
  USE MOD_DSMC_Analyze,           ONLY : CalcGammaVib, CalcInstantTransTemp, CalcMeanFreePath
  USE MOD_Particle_Analyze_Vars,  ONLY : CalcEkin
  USE MOD_DSMC_CollisionProb,     ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,            ONLY : DSMC_perform_collision, DSMC_calc_var_P_vib
  USE MOD_vmpf_collision,         ONLY : DSMC_vmpf_prob
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
  USE MOD_part_tools,             ONLY : GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                        :: NodeVolume
  INTEGER, INTENT(IN)                     :: PartNum
  INTEGER, INTENT(IN)                     :: iElem
  INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPair, iPart1, iPart2, iLoop, iPart, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase , PairNum_Node
  REAL                          :: Dist1, Dist2, iRan, VibProb
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

  PairNum_Node = INT(PartNum/2)
  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  IF(VarTimeStep%UseVariableTimeStep) CollInf%MeanMPF = 0.

  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
    ChemReac%nPairForRec = 0
    ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
    DO iPart = 1, PartNum
      ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart)))=ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) &
                                              + PartStateIntEn(iPartIndx_Node(iPart),1)*GetParticleWeight(iPartIndx_Node(iPart))
    END DO
  END IF

  DO iPart = 1, PartNum
    CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
              CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + GetParticleWeight(iPartIndx_Node(iPart))
  END DO

  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
    ! 3. Case: Temperature required for the mean free path with the VHS model
    CALL CalcInstantTransTemp(iPartIndx_Node,PartNum)
    IF((SelectionProc.EQ.2).OR.(DSMC%VibRelaxProb.EQ.2.0)) CALL CalcGammaVib()
  END IF

  ALLOCATE(Coll_pData(PairNum_Node))
  nPart = PartNum

  IF (usevMPF) MPFFac = 1
  DO iPair = 1, PairNum_Node
    CALL RANDOM_NUMBER(iRan)
    iPart1 = 1 + INT(nPart * iRan)
    Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
    iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
    nPart = nPart - 1
    iPart2 = 1
    Dist1 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
           - PartState(iPartIndx_Node(iPart2),1))**2 &
           +(PartState(Coll_pData(iPair)%iPart_p1,2) &
           - PartState(iPartIndx_Node(iPart2),2))**2 &
           +(PartState(Coll_pData(iPair)%iPart_p1,3) &
           - PartState(iPartIndx_Node(iPart2),3))**2
    DO iLoop = 2, nPart
      Dist2 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
             - PartState(iPartIndx_Node(iLoop),1))**2 &
             +(PartState(Coll_pData(iPair)%iPart_p1,2) &
             - PartState(iPartIndx_Node(iLoop),2))**2 &
             +(PartState(Coll_pData(iPair)%iPart_p1,3) &
             - PartState(iPartIndx_Node(iLoop),3))**2
      IF (Dist2.LT.Dist1) THEN
        iPart2 = iLoop
        Dist1 = Dist2
      END IF
    END DO
    Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
    iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

    IF (usevMPF) THEN
      TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
      IF (TempMPFFac .GE. MPFFac) THEN
          MPFFac = TempMPFFac
          PairE_vMPF(1) = iPair
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
        ELSE
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
        END IF
      END IF
    END IF

    iCase = CollInf%Coll_Case(cSpec1, cSpec2)

    IF(VarTimeStep%UseVariableTimeStep) THEN
      CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(Coll_pData(iPair)%iPart_p1) &
                                                           + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
    END IF

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
  IF ((PairNum_Node.NE.0).AND.(CollisMode.EQ.3).and.(MOD(PartNum, PairNum_Node).NE.0)) THEN
    ChemReac%RecombParticle = iPartIndx_Node(1)
  END IF

  IF (CollisMode.EQ.3) THEN
    CALL SetMeanVibQua()
  END IF

  DO iPair = 1,  PairNum_Node
    IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
      IF (usevMPF) THEN            ! calculation of collision prob
        CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
      ELSE
        CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
      END IF
      ! variable vibrational relaxation probability has to average of all collisions
      IF(DSMC%VibRelaxProb.EQ.2.0) THEN
        cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
        cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
        IF((SpecDSMC(cSpec1)%InterID.EQ.2).OR.(SpecDSMC(cSpec1)%InterID.EQ.20)) THEN
          CALL DSMC_calc_var_P_vib(cSpec1,cSpec2,iPair,VibProb)
          VarVibRelaxProb%ProbVibAvNew(cSpec1) = VarVibRelaxProb%ProbVibAvNew(cSpec1) + VibProb
          VarVibRelaxProb%nCollis(cSpec1) = VarVibRelaxProb%nCollis(cSpec1) + 1
        END IF
        IF((SpecDSMC(cSpec2)%InterID.EQ.2).OR.(SpecDSMC(cSpec2)%InterID.EQ.20)) THEN
          CALL DSMC_calc_var_P_vib(cSpec2,cSpec1,iPair,VibProb)
          VarVibRelaxProb%ProbVibAvNew(cSpec2) = VarVibRelaxProb%ProbVibAvNew(cSpec2) + VibProb
          VarVibRelaxProb%nCollis(cSpec2) = VarVibRelaxProb%nCollis(cSpec2) + 1
        END IF
      END IF
      CALL RANDOM_NUMBER(iRan)
      IF (Coll_pData(iPair)%Prob.GE.iRan) THEN
#if (PP_TimeDiscMethod==42)
        IF(CalcEkin.OR.DSMC%ReservoirSimu) THEN
#else
        IF(CalcEkin) THEN
#endif
          DSMC%NumColl(Coll_pData(iPair)%PairType) = DSMC%NumColl(Coll_pData(iPair)%PairType) + 1
          DSMC%NumColl(CollInf%NumCase + 1) = DSMC%NumColl(CollInf%NumCase + 1) + 1
        END IF
        CALL DSMC_perform_collision(iPair, iElem, NodeVolume, PartNum)  ! call collision from octree
      END IF
    END IF
  END DO

  IF(DSMC%CalcQualityFactors) THEN
    IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
      ! Calculation of the mean free path with VHS model and the current translational temperature in the cell
      DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                                            SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
      ! Determination of the maximum MCS/MFP for the cell
    IF((DSMC%CollSepCount.GT.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = &
                                                    MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
    END IF
  END IF

  DEALLOCATE(Coll_pData)

END SUBROUTINE FindNearestNeigh


SUBROUTINE DSMC_pairing_statistical(iElem)
!===================================================================================================================================
! Classic statistical pairing method
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, CollisMode, PartStateIntEn, ChemReac, PairE_vMPF, CRelaMax, CRelaAv
  USE MOD_DSMC_Vars,              ONLY : DSMC, SelectionProc, RadialWeighting
  USE MOD_DSMC_Analyze,           ONLY : CalcGammaVib, CalcInstantTransTemp
  USE MOD_Particle_Vars,          ONLY : PEM, PartSpecies, nSpecies, PartState, usevMPF, PartMPF, VarTimeStep
  USE MOD_Particle_Vars,          ONLY : KeepWallParticles, PDM
  USE MOD_part_tools,             ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: nPair, iPair, iPart, iLoop, cPart1, cPart2, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase
  INTEGER, ALLOCATABLE          :: iPartIndx(:) ! List of particles in the cell nec for stat pairing
  REAL                          :: iRan
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================
  IF (KeepWallParticles) THEN
    nPart = PEM%pNumber(iElem)-PEM%wNumber(iElem)
  ELSE
    nPart = PEM%pNumber(iElem)
  END IF
  nPair = INT(nPart/2)
  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
  END IF

  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) CollInf%MeanMPF = 0.

  CollInf%Coll_SpecPartNum = 0.
  CollInf%Coll_CaseNum = 0
  CRelaMax = 0
  CRelaAv = 0
  ALLOCATE(Coll_pData(nPair))
  ALLOCATE(iPartIndx(nPart))
  Coll_pData%Ec=0
  iPartIndx = 0
  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    ! check if particle is on wall and chose next particle until particle is not at wall
    IF (KeepWallParticles) THEN
      DO WHILE (PDM%ParticleAtWall(iPart))
        iPart = PEM%pNext(iPart)
      END DO
    END IF
    iPartIndx(iLoop) = iPart
    ! Counter for part num of spec per cell
    CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    ! Calculation of mean evib per cell and iter, necessary for disso prob
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) = ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) &
                                                                  + PartStateIntEn(iPart,1) * GetParticleWeight(iPart)
    ! Choose next particle in Element
    iPart = PEM%pNext(iPart)
  END DO

  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
    ! 3. Case: Temperature required for the mean free path with the VHS model
    CALL CalcInstantTransTemp(iPartIndx,nPart)
    IF(SelectionProc.EQ.2) CALL CalcGammaVib()
  END IF

  IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) MPFFac = 1

  DO iPair = 1, nPair                               ! statistical pairing
    CALL RANDOM_NUMBER(iRan)
    cPart1 = 1 + INT(nPart * iRan)                       ! first pair particle
    Coll_pData(iPair)%iPart_p1 = iPartIndx(cPart1)
    iPartIndx(cPart1) = iPartIndx(nPart)
    nPart = nPart - 1
    CALL RANDOM_NUMBER(iRan)
    cPart2 = 1 + INT(nPart * iRan)                       ! second pair particle
    Coll_pData(iPair)%iPart_p2 = iPartIndx(cPart2)
    iPartIndx(cPart2) = iPartIndx(nPart)
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

    IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
      TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
      IF (TempMPFFac .GE. MPFFac) THEN
          MPFFac = TempMPFFac
          PairE_vMPF(1) = iPair
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
        ELSE
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
        END IF
      END IF
    END IF

    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(Coll_pData(iPair)%iPart_p1) &
                                                          + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
    END IF
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
    Coll_pData(iPair)%PairType = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
    ! get maximum and average relative velocity
    CRelaMax = MAX(CRelaMax, SQRT(Coll_pData(iPair)%CRela2))
    CRelaAv = CRelaAv + SQRT(Coll_pData(iPair)%CRela2)
  END DO
  IF(nPair.NE.0)  CRelaAv = CRelaAv / nPair
  IF ((nPair.NE.0).AND.(CollisMode.EQ.3).AND.(MOD(nPart, nPair).NE.0)) THEN
    ChemReac%RecombParticle = iPartIndx(1)
  END IF

  DEALLOCATE(iPartIndx)
END SUBROUTINE DSMC_pairing_statistical


SUBROUTINE FindNearestNeigh2D(iPartIndx_Node, PartNum, iElem, NodeVolume, MidPoint, Depth)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY: CollInf, tTreeNode, CollisMode, ChemReac, PartStateIntEn, Coll_pData, SelectionProc
USE MOD_DSMC_Vars,              ONLY: DSMC, PairE_vMPF, RadialWeighting, SamplingActive, SpecDSMC, VarVibRelaxProb
USE MOD_DSMC_Symmetry2D,        ONLY: CalcRadWeightMPF
USE MOD_Particle_Vars,          ONLY: PartState, nSpecies, PartSpecies, usevMPF, PartMPF, WriteMacroVolumeValues, VarTimeStep
USE MOD_DSMC_Relaxation,        ONLY: SetMeanVibQua
USE MOD_DSMC_Analyze,           ONLY: CalcGammaVib, CalcInstantTransTemp, CalcMeanFreePath
USE MOD_Particle_Analyze_Vars,  ONLY: CalcEkin
USE MOD_DSMC_CollisionProb,     ONLY: DSMC_prob_calc
USE MOD_DSMC_Collis,            ONLY: DSMC_perform_collision, DSMC_calc_var_P_vib
USE MOD_vmpf_collision,         ONLY: DSMC_vmpf_prob
USE MOD_TimeDisc_Vars,          ONLY: TEnd, Time
USE MOD_part_tools,             ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: NodeVolume, MidPoint(3)
INTEGER, INTENT(IN)             :: PartNum
INTEGER, INTENT(IN)             :: iElem, Depth
INTEGER, INTENT(INOUT)          :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPair, iPart1, iPart2, iLoop, iPart, nPart, loopStart
INTEGER                         :: cSpec1, cSpec2, iCase , PairNum_Node, tempPart
REAL                            :: Dist1, Dist2, iRan, VibProb
REAL                            :: TempMPFFac, MPFFac, iRanVec(2), NodeLength
!===================================================================================================================================

PairNum_Node = INT(PartNum/2)
CollInf%Coll_SpecPartNum = 0
CollInf%Coll_CaseNum = 0
IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) CollInf%MeanMPF = 0.

IF (CollisMode.EQ.3) THEN
  ChemReac%RecombParticle = 0
  ChemReac%nPairForRec = 0
END IF

IF (CollisMode.EQ.3) THEN
  ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
  DO iPart = 1, PartNum
    ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart)))=ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) &
      + PartStateIntEn(iPartIndx_Node(iPart),1) * GetParticleWeight(iPartIndx_Node(iPart))
  END DO
END IF


DO iPart = 1, PartNum
  CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
            CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + GetParticleWeight(iPartIndx_Node(iPart))
END DO

IF((CollisMode.GT.1).AND.(SelectionProc.EQ.2).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
  CALL CalcInstantTransTemp(iPartIndx_Node,PartNum)
  IF(SelectionProc.EQ.2) CALL CalcGammaVib()
END IF

ALLOCATE(Coll_pData(PairNum_Node))
nPart = PartNum

IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) MPFFac = 1
DO iPair = 1, PairNum_Node
  loopStart = 0
  CALL RANDOM_NUMBER(iRan)
  iPart1 = 1 + INT(nPart * iRan)
  Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
  iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
  nPart = nPart - 1
  iPart2 = 1
  IF (CollInf%ProhibitDoubleColl) THEN
    IF (nPart.GT.1) THEN
      IF (iPartIndx_Node(ipart2).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
        iPart2 = 2
        loopStart = 1
      END IF
    END IF
  END IF
  Dist1 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
          - PartState(iPartIndx_Node(iPart2),1))**2 &
          +(PartState(Coll_pData(iPair)%iPart_p1,2) &
          - PartState(iPartIndx_Node(iPart2),2))**2
  DO iLoop = 2 + loopStart, nPart
    IF (CollInf%ProhibitDoubleColl) THEN
        IF (iPartIndx_Node(iLoop).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
          CYCLE
        END IF
    END IF
    Dist2 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
            - PartState(iPartIndx_Node(iLoop),1))**2 &
            +(PartState(Coll_pData(iPair)%iPart_p1,2) &
            - PartState(iPartIndx_Node(iLoop),2))**2
    IF (Dist2.LT.Dist1) THEN
      iPart2 = iLoop
      Dist1 = Dist2
    END IF
  END DO
  Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
  iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
  nPart = nPart - 1

  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

  IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
    TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
    IF (TempMPFFac .GE. MPFFac) THEN
        MPFFac = TempMPFFac
        PairE_vMPF(1) = iPair
      IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
        PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
      ELSE
        PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
      END IF
    END IF
  END IF

  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(Coll_pData(iPair)%iPart_p1) &
                                                        + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
  END IF

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

IF (CollInf%ProhibitDoubleColl.AND.(nPart.EQ.1)) CollInf%OldCollPartner(iPartIndx_Node(1)) = 0

IF ((PairNum_Node.NE.0).AND.(CollisMode.EQ.3).and.(MOD(PartNum, PairNum_Node).NE.0)) THEN
  ChemReac%RecombParticle = iPartIndx_Node(1)
END IF

IF (CollisMode.EQ.3) THEN
  CALL SetMeanVibQua()
END IF
DO iPair = 1,  PairNum_Node
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    IF(RadialWeighting%DoRadialWeighting) THEN
      ! Two particles with the exact same velocities at the same positions -> sick clones
      IF (Coll_pData(iPair)%CRela2.EQ.0.0) THEN
        IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
          IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
        END IF
        IF (iPair.LT.PairNum_Node) THEN
          ! "Partner-Tausch": if there are pairs ahead in the pairing list, the next is pair is broken up and collision partners
          ! are swapped but first, changing z-direction
          PartState(Coll_pData(iPair)%iPart_p1,6) = - PartState(Coll_pData(iPair)%iPart_p1,6)
          ! Removing the pairs from the weighting factor and the case num sums
          CollInf%MeanMPF(Coll_pData(iPair)%PairType) = CollInf%MeanMPF(Coll_pData(iPair)%PairType) &
            -(GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
          CollInf%MeanMPF(Coll_pData(iPair+1)%PairType) = CollInf%MeanMPF(Coll_pData(iPair+1)%PairType) &
            - (GetParticleWeight(Coll_pData(iPair+1)%iPart_p1) + GetParticleWeight(Coll_pData(iPair+1)%iPart_p2))*0.5
          CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
          CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) - 1
          ! Breaking up the next pair and swapping partners
          tempPart = Coll_pData(iPair)%iPart_p1
          Coll_pData(iPair)%iPart_p1 = Coll_pData(iPair + 1)%iPart_p1
          Coll_pData(iPair + 1)%iPart_p1 = tempPart
          ! Calculation of the relative velocity for the new first pair
          cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
          cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
          iCase = CollInf%Coll_Case(cSpec1, cSpec2)
          ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
          CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(cSpec1) + GetParticleWeight(cSpec2))*0.5
          CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
          Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                                    -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                                    + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                                    -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                                    + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                                    -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
          Coll_pData(iPair)%PairType = iCase
          ! Calculation of the relative velocity for the new follow-up pair
          cSpec1 = PartSpecies(Coll_pData(iPair+1)%iPart_p1)
          cSpec2 = PartSpecies(Coll_pData(iPair+1)%iPart_p2)
          iCase = CollInf%Coll_Case(cSpec1, cSpec2)
          ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
          CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(cSpec1) + GetParticleWeight(cSpec2))*0.5
          CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
          Coll_pData(iPair+1)%CRela2 = (PartState(Coll_pData(iPair+1)%iPart_p1,4) &
                                      -  PartState(Coll_pData(iPair+1)%iPart_p2,4))**2 &
                                      + (PartState(Coll_pData(iPair+1)%iPart_p1,5) &
                                      -  PartState(Coll_pData(iPair+1)%iPart_p2,5))**2 &
                                      + (PartState(Coll_pData(iPair+1)%iPart_p1,6) &
                                      -  PartState(Coll_pData(iPair+1)%iPart_p2,6))**2
          Coll_pData(iPair+1)%PairType = iCase
        ELSE
          ! For the last pair, the z-direction is changed and a random position for the first particle is chosen (last resort)
          CALL RANDOM_NUMBER(iRanVec)
          NodeLength = 2./2.**(Depth)
          PartState(Coll_pData(iPair)%iPart_p1,1) = MidPoint(1) - NodeLength + 2.*NodeLength*iRanVec(1)
          PartState(Coll_pData(iPair)%iPart_p1,2) = MidPoint(2) - NodeLength + 2.*NodeLength*iRanVec(2)
          PartState(Coll_pData(iPair)%iPart_p1,1:2) = MapToGeo2D(PartState(Coll_pData(iPair)%iPart_p1,1:2),iElem)
          PartState(Coll_pData(iPair)%iPart_p1,3) = 0.0
          ! creating relative velocity
          PartState(Coll_pData(iPair)%iPart_p1,6) = - PartState(Coll_pData(iPair)%iPart_p1,6)
          Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,6) - PartState(Coll_pData(iPair)%iPart_p2,6))**2
          PartMPF(Coll_pData(iPair)%iPart_p1) = &
                              CalcRadWeightMPF(PartState(Coll_pData(iPair)%iPart_p1,2),cSpec1,Coll_pData(iPair)%iPart_p1)
        END IF
      END IF
    END IF
    IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN            ! calculation of collision prob
      CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
    ELSE
      CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
    END IF
    ! variable vibrational relaxation probability has to average of all collisions
    IF(DSMC%VibRelaxProb.EQ.2.0) THEN
      cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
      cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
      IF((SpecDSMC(cSpec1)%InterID.EQ.2).OR.(SpecDSMC(cSpec1)%InterID.EQ.20)) THEN
        CALL DSMC_calc_var_P_vib(cSpec1,cSpec2,iPair,VibProb)
        VarVibRelaxProb%ProbVibAvNew(cSpec1) = VarVibRelaxProb%ProbVibAvNew(cSpec1) + VibProb
        VarVibRelaxProb%nCollis(cSpec1) = VarVibRelaxProb%nCollis(cSpec1) + 1
      END IF
      IF((SpecDSMC(cSpec2)%InterID.EQ.2).OR.(SpecDSMC(cSpec2)%InterID.EQ.20)) THEN
        CALL DSMC_calc_var_P_vib(cSpec2,cSpec1,iPair,VibProb)
        VarVibRelaxProb%ProbVibAvNew(cSpec2) = VarVibRelaxProb%ProbVibAvNew(cSpec2) + VibProb
        VarVibRelaxProb%nCollis(cSpec2) = VarVibRelaxProb%nCollis(cSpec2) + 1
      END IF
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
#if (PP_TimeDiscMethod==42)
      IF(CalcEkin.OR.DSMC%ReservoirSimu) THEN
#else
      IF(CalcEkin) THEN
#endif
        DSMC%NumColl(Coll_pData(iPair)%PairType) = DSMC%NumColl(Coll_pData(iPair)%PairType) + 1
        DSMC%NumColl(CollInf%NumCase + 1) = DSMC%NumColl(CollInf%NumCase + 1) + 1
      END IF
      IF (CollInf%ProhibitDoubleColl) THEN
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1) = Coll_pData(iPair)%iPart_p2
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p2) = Coll_pData(iPair)%iPart_p1
      END IF
      CALL DSMC_perform_collision(iPair, iElem, NodeVolume, PartNum)  ! call collision from octree
    ELSE
      IF (CollInf%ProhibitDoubleColl) THEN
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1) = 0
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p2) = 0
      END IF
    END IF
  END IF
END DO

IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path with VHS model and the current translational temperature in the cell
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                                          SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the maximum MCS/MFP for the cell
    IF((DSMC%CollSepCount.GT.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = &
                                                    MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
  END IF
END IF

DEALLOCATE(Coll_pData)

END SUBROUTINE FindNearestNeigh2D


SUBROUTINE DSMC_pairing_octree(iElem)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
  USE MOD_DSMC_Vars               ,ONLY: tTreeNode, DSMC, ElemNodeVol, VarVibRelaxProb
  USE MOD_Particle_Vars           ,ONLY: PEM, PartState, nSpecies, PartSpecies,PartPosRef
  USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
  USE MOD_Particle_Tracking_vars  ,ONLY: DoRefMapping
  USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
  USE MOD_part_tools,             ONLY : GetParticleWeight
  USE MOD_TimeDisc_Vars,          ONLY : iter
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, nPart, iSpec
  REAL                          :: SpecPartNum(nSpecies)
  TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

SpecPartNum = 0.
nPart = PEM%pNumber(iElem)
IF(DSMC%VibRelaxProb.EQ.2.0) THEN
  DO iSpec=1,nSpecies
    VarVibRelaxProb%ProbVibAvNew(iSpec) = 0
    VarVibRelaxProb%nCollis(iSpec) = 0
  END DO
END IF

IF (nPart.GT.1) THEN
  NULLIFY(TreeNode)

  ALLOCATE(TreeNode)
  ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing
  TreeNode%iPartIndx_Node(1:nPart) = 0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    TreeNode%iPartIndx_Node(iLoop) = iPart
    iPart = PEM%pNext(iPart)
    ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
    SpecPartNum(PartSpecies(TreeNode%iPartIndx_Node(iLoop))) = &
              SpecPartNum(PartSpecies(TreeNode%iPartIndx_Node(iLoop))) + GetParticleWeight(TreeNode%iPartIndx_Node(iLoop))
  END DO

  DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum, SUM(SpecPartNum), GEO%Volume(iElem))
  ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
  IF(nPart.GE.DSMC%PartNumOctreeNodeMin) THEN
    ! Additional check afterwards if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
    ! the side length of a cube (approximation) with same volume as the actual cell -> octree
    IF((DSMC%MeanFreePath.LT.(GEO%CharLength(iElem))) .OR.(nPart.GT.DSMC%PartNumOctreeNode)) THEN
      ALLOCATE(TreeNode%MappedPartStates(1:nPart, 1:3))
      TreeNode%PNum_Node = nPart
      iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
      IF (DoRefMapping) THEN
        DO iLoop = 1, nPart
          TreeNode%MappedPartStates(iLoop,1:3)=PartPosRef(1:3,iPart)
          iPart = PEM%pNext(iPart)
        END DO
      ELSE ! position in reference space [-1,1] has to be computed
        DO iLoop = 1, nPart
          CALL GetPositionInRefElem(PartState(iPart,1:3),TreeNode%MappedPartStates(iLoop,1:3),iElem)
          iPart = PEM%pNext(iPart)
        END DO
      END IF ! DoRefMapping
      TreeNode%NodeDepth = 1
      TreeNode%MidPoint(1:3) = (/0.0,0.0,0.0/)
      CALL AddOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
      DEALLOCATE(TreeNode%MappedPartStates)
    ELSE
      IF(DSMC%UseNearestNeighbour) THEN
        CALL FindNearestNeigh(TreeNode%iPartIndx_Node, nPart, iElem, GEO%Volume(iElem))
      ELSE
        CALL FindStatisticalNeigh(TreeNode%iPartIndx_Node,nPart,iElem, GEO%Volume(iElem))
      END IF
    END IF
  ELSE  IF (nPart.GT.1) THEN
    IF(DSMC%UseNearestNeighbour) THEN
      CALL FindNearestNeigh(TreeNode%iPartIndx_Node, nPart, iElem, GEO%Volume(iElem))
    ELSE
      CALL FindStatisticalNeigh(TreeNode%iPartIndx_Node,nPart,iElem, GEO%Volume(iElem))
    END IF
  END IF

  DEALLOCATE(TreeNode%iPartIndx_Node)
  DEALLOCATE(TreeNode)
END IF !nPart > 0

IF(DSMC%VibRelaxProb.EQ.2.0) THEN
  DO iSpec=1,nSpecies
    IF(VarVibRelaxProb%nCollis(iSpec).NE.0) THEN
      VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VarVibRelaxProb%ProbVibAv(iElem,iSpec) &
                                             * VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec)) &
                                             + (1.-VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec))) &
                                             / (VarVibRelaxProb%nCollis(iSpec)) * VarVibRelaxProb%ProbVibAvNew(iSpec) 
      ! open(unit=226,file='ProbVibAv.csv',action='write',position='append')
      !   WRITE(226,*) Iter,',',ProbVibAv(iElem),',',2*nCollis
      ! CLOSE(Unit=226)
    END IF
  END DO
END IF

END SUBROUTINE DSMC_pairing_octree

RECURSIVE SUBROUTINE AddOctreeNode(TreeNode, iElem, NodeVol)
!===================================================================================================================================
! Adds additional octree node/branch (fancy)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Analyze,           ONLY : CalcMeanFreePath
  USE MOD_DSMC_Vars,              ONLY : tTreeNode, DSMC, tNodeVolume
  USE MOD_Particle_Vars,          ONLY : nSpecies, PartSpecies
  USE MOD_DSMC_Vars,              ONLY : ElemNodeVol
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                     :: iElem
  TYPE(tTreeNode),INTENT(INOUT), POINTER     :: TreeNode
  TYPE(tNodeVolume),INTENT(IN), POINTER   :: NodeVol
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, iPartIndx, localDepth, SpecPartNum(nSpecies)
  INTEGER, ALLOCATABLE          :: iPartIndx_ChildNode(:,:)
  REAL, ALLOCATABLE             :: MappedPart_ChildNode(:,:,:)
  INTEGER                       :: PartNumChildNode(8)
  REAL                          :: NodeVolumeTemp(8)
!===================================================================================================================================

  ALLOCATE(iPartIndx_ChildNode(8,TreeNode%PNum_Node))
  ALLOCATE(MappedPart_ChildNode(8,TreeNode%PNum_Node,3))
  PartNumChildNode(:) = 0
  IF (ABS(TreeNode%MidPoint(1)) .EQ. 1.0) THEN
    CALL Abort(&
__STAMP__&
,'ERROR in Octree Pairing: Too many branches, machine precision reached')
  END IF
  !         Numbering of the 8 ChildNodes of the Octree
  !          __________
  !         /    /    /|   |z
  !        /  7 / 6  / |   |
  !       /----/----/|6|   |
  !      /_3__/_2__/ |/|   |______Y
  !     |    |    |2 |5|   /
  !     | 3  | 2  |/ |/   /
  !     |----|----|1 /   /x
  !     | 4  | 1  | /
  !     |____|____|/

  ! particle to Octree ChildNode sorting
  DO iPart=1,TreeNode%PNum_Node
    iPartIndx = TreeNode%iPartIndx_Node(iPart)
    IF ((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
        .AND.(TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) &
        .AND.(TreeNode%MappedPartStates(iPart,3).LE.TreeNode%MidPoint(3))) THEN
      PartNumChildNode(1) = PartNumChildNode(1) + 1
      iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
      MappedPart_ChildNode(1,PartNumChildNode(1),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
        .AND.(TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2))) THEN
      PartNumChildNode(2) = PartNumChildNode(2) + 1
      iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
      MappedPart_ChildNode(2,PartNumChildNode(2),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
        .AND.(TreeNode%MappedPartStates(iPart,3).GE.TreeNode%MidPoint(3))) THEN
      PartNumChildNode(3) = PartNumChildNode(3) + 1
      iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
      MappedPart_ChildNode(3,PartNumChildNode(3),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF (TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) THEN
      PartNumChildNode(4) = PartNumChildNode(4) + 1
      iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
      MappedPart_ChildNode(4,PartNumChildNode(4),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF((TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) &
        .AND.(TreeNode%MappedPartStates(iPart,3).LE.TreeNode%MidPoint(3))) THEN
      PartNumChildNode(5) = PartNumChildNode(5) + 1
      iPartIndx_ChildNode(5,PartNumChildNode(5)) = iPartIndx
      MappedPart_ChildNode(5,PartNumChildNode(5),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF (TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) THEN
      PartNumChildNode(6) = PartNumChildNode(6) + 1
      iPartIndx_ChildNode(6,PartNumChildNode(6)) = iPartIndx
      MappedPart_ChildNode(6,PartNumChildNode(6),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE IF (TreeNode%MappedPartStates(iPart,3).GE.TreeNode%MidPoint(3)) THEN
      PartNumChildNode(7) = PartNumChildNode(7) + 1
      iPartIndx_ChildNode(7,PartNumChildNode(7)) = iPartIndx
      MappedPart_ChildNode(7,PartNumChildNode(7),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    ELSE
      PartNumChildNode(8) = PartNumChildNode(8) + 1
      iPartIndx_ChildNode(8,PartNumChildNode(8)) = iPartIndx
      MappedPart_ChildNode(8,PartNumChildNode(8),1:3) = TreeNode%MappedPartStates(iPart,1:3)
    END IF
  END DO

  IF((.NOT.ASSOCIATED(NodeVol)).OR.(.NOT.ASSOCIATED(NodeVol%SubNode1))) THEN
    localDepth = TreeNode%NodeDepth
    CALL DSMC_CalcSubNodeVolumes(iElem, localDepth, ElemNodeVol(iElem)%Root)
  END IF

  NodeVolumeTemp(1) = NodeVol%SubNode1%Volume
  NodeVolumeTemp(2) = NodeVol%SubNode2%Volume
  NodeVolumeTemp(3) = NodeVol%SubNode3%Volume
  NodeVolumeTemp(4) = NodeVol%SubNode4%Volume
  NodeVolumeTemp(5) = NodeVol%SubNode5%Volume
  NodeVolumeTemp(6) = NodeVol%SubNode6%Volume
  NodeVolumeTemp(7) = NodeVol%SubNode7%Volume
  NodeVolumeTemp(8) = NodeVol%SubNode8%Volume

  DO iLoop = 1, 8
    IF (PartNumChildNode(iLoop).GT.1) THEN
      ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
      SpecPartNum = 0
      DO iPart = 1, PartNumChildNode(iLoop)
        SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) = &
          SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) + 1
      END DO
      DSMC%MeanFreePath = CalcMeanFreePath(REAL(SpecPartNum),REAL(PartNumChildNode(iLoop)),NodeVolumeTemp(iLoop))
    END IF
    ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
    IF(PartNumChildNode(iLoop).GE.DSMC%PartNumOctreeNodeMin) THEN
      ! Additional check if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
      ! the side length of a cube (approximation) with same volume as the actual cell -> octree
      IF((DSMC%MeanFreePath.LT.(NodeVolumeTemp(iLoop)**(1./3.))) &
                                                               .OR.(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode)) THEN
        NULLIFY(TreeNode%ChildNode)
        ALLOCATE(TreeNode%ChildNode)
        ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
        ALLOCATE(TreeNode%ChildNode%MappedPartStates(PartNumChildNode(iLoop),1:3))
        TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
        TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
        TreeNode%ChildNode%MappedPartStates(1:PartNumChildNode(iLoop),1:3)= &
              MappedPart_ChildNode(iLoop,1:PartNumChildNode(iLoop),1:3)
        IF (iLoop.LT.5) THEN
          TreeNode%ChildNode%MidPoint(1) = 1.0
          IF (iLoop.LT.3) THEN
            TreeNode%ChildNode%MidPoint(2) = 1.0
          ELSE
            TreeNode%ChildNode%MidPoint(2) = -1.0
          END IF
        ELSE
          TreeNode%ChildNode%MidPoint(1) = -1.0
          IF (iLoop.LT.7) THEN
            TreeNode%ChildNode%MidPoint(2) = 1.0
          ELSE
            TreeNode%ChildNode%MidPoint(2) = -1.0
          END IF
        END IF
        IF ((iLoop.EQ.1).OR.(iLoop.EQ.4).OR.(iLoop.EQ.5).OR.(iLoop.EQ.8)) THEN
          TreeNode%ChildNode%MidPoint(3) = -1.0
        ELSE
          TreeNode%ChildNode%MidPoint(3) = 1.0
        END IF
        TreeNode%ChildNode%MidPoint(1:3) = TreeNode%MidPoint(1:3) &
                                         + TreeNode%ChildNode%MidPoint(1:3)*2.0/REAL(2.0**(TreeNode%NodeDepth+1))
        TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
        ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further octree division)
        IF (iLoop.EQ.1) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1)
        IF (iLoop.EQ.2) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2)
        IF (iLoop.EQ.3) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3)
        IF (iLoop.EQ.4) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4)
        IF (iLoop.EQ.5) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode5)
        IF (iLoop.EQ.6) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode6)
        IF (iLoop.EQ.7) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode7)
        IF (iLoop.EQ.8) CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode8)
        DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
        DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
        DEALLOCATE(TreeNode%ChildNode)
      ELSE
        IF(DSMC%UseNearestNeighbour) THEN
        CALL FindNearestNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
        ELSE
        CALL FindStatisticalNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
        END IF
      END IF
    ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
      IF(DSMC%UseNearestNeighbour) THEN
      CALL FindNearestNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
      ELSE
      CALL FindStatisticalNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
      END IF
    END IF
  END DO

END SUBROUTINE AddOctreeNode


SUBROUTINE DSMC_pairing_quadtree(iElem)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, DSMC, ElemNodeVol, CollInf, VarVibRelaxProb
USE MOD_Particle_Vars           ,ONLY: PEM, PartState, nSpecies, PartSpecies
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_part_tools              ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, nPart, iSpec
REAL                          :: SpecPartNum(nSpecies), Volume
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

  Volume = GEO%Volume(iElem)
  IF(DSMC%VibRelaxProb.EQ.2.0) THEN
    DO iSpec=1,nSpecies
      VarVibRelaxProb%ProbVibAvNew(iSpec) = 0
      VarVibRelaxProb%nCollis(iSpec) = 0
    END DO
  END IF
  SpecPartNum = 0.

  NULLIFY(TreeNode)
  nPart = PEM%pNumber(iElem)
  IF(nPart.GT.0) THEN
    ALLOCATE(TreeNode)
    ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing
    TreeNode%iPartIndx_Node(1:nPart) = 0

    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing

    DO iLoop = 1, nPart
      TreeNode%iPartIndx_Node(iLoop) = iPart
      ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
      SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
      iPart = PEM%pNext(iPart)
    END DO

    DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum, SUM(SpecPartNum), Volume)

    ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
    IF(nPart.GE.DSMC%PartNumOctreeNodeMin) THEN
      ! Additional check afterwards if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
      ! the side length of a cube (approximation) with same volume as the actual cell -> octree
      IF((DSMC%MeanFreePath.LT.GEO%CharLength(iElem)).OR.(nPart.GT.DSMC%PartNumOctreeNode)) THEN
        ALLOCATE(TreeNode%MappedPartStates(1:nPart, 1:2))
        TreeNode%PNum_Node = nPart
        iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
        DO iLoop = 1, nPart
          CALL GeoCoordToMap2D(PartState(iPart,1:2), TreeNode%MappedPartStates(iLoop,1:2), iElem)
          iPart = PEM%pNext(iPart)
        END DO
        TreeNode%NodeDepth = 1
        TreeNode%MidPoint(1:3) = (/0.0,0.0,0.0/)
        CALL AddQuadTreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
        DEALLOCATE(TreeNode%MappedPartStates)
      ELSE
        IF(DSMC%UseNearestNeighbour) THEN
          CALL FindNearestNeigh2D(TreeNode%iPartIndx_Node, nPart, iElem, GEO%Volume(iElem),  (/0.0,0.0,0.0/), 1)
        ELSE
          CALL FindStatisticalNeigh(TreeNode%iPartIndx_Node,nPart,iElem, GEO%Volume(iElem))
        END IF
      END IF
    ELSE  IF (nPart.GT.1) THEN
      IF(DSMC%UseNearestNeighbour) THEN
        CALL FindNearestNeigh2D(TreeNode%iPartIndx_Node, nPart, iElem, GEO%Volume(iElem),  (/0.0,0.0,0.0/), 1)
      ELSE
        CALL FindStatisticalNeigh(TreeNode%iPartIndx_Node,nPart,iElem, GEO%Volume(iElem))
      END IF
    ELSE IF (CollInf%ProhibitDoubleColl.AND.(nPart.EQ.1)) THEN
      CollInf%OldCollPartner(TreeNode%iPartIndx_Node(1)) = 0
    END IF

    DEALLOCATE(TreeNode%iPartIndx_Node)
    DEALLOCATE(TreeNode)
  END IF

  IF(DSMC%VibRelaxProb.EQ.2.0) THEN
    DO iSpec=1,nSpecies
      IF(VarVibRelaxProb%nCollis(iSpec).NE.0) THEN
        VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VarVibRelaxProb%ProbVibAv(iElem,iSpec) &
                                               * VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec)) &
                                               + (1.-VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec))) &
                                               / (VarVibRelaxProb%nCollis(iSpec)) * VarVibRelaxProb%ProbVibAvNew(iSpec) 
        ! open(unit=226,file='ProbVibAv.csv',action='write',position='append')
        !   WRITE(226,*) Iter,',',ProbVibAv(iElem),',',2*nCollis
        ! CLOSE(Unit=226)
      END IF
    END DO
  END IF

END SUBROUTINE DSMC_pairing_quadtree


RECURSIVE SUBROUTINE AddQuadTreeNode(TreeNode, iElem, NodeVol)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze      ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars         ,ONLY: tTreeNode, DSMC, tNodeVolume, RadialWeighting, CollInf
USE MOD_Particle_Vars     ,ONLY: nSpecies, PartSpecies, VarTimeStep
USE MOD_DSMC_Vars         ,ONLY: ElemNodeVol
USE MOD_part_tools        ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                     :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER     :: TreeNode
TYPE(tNodeVolume),INTENT(IN), POINTER   :: NodeVol
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iPartIndx, localDepth
INTEGER, ALLOCATABLE          :: iPartIndx_ChildNode(:,:)
REAL, ALLOCATABLE             :: MappedPart_ChildNode(:,:,:)
INTEGER                       :: PartNumChildNode(4)
REAL                          :: NodeVolumeTemp(4), FaceVolumeTemp(4), SpecPartNum(nSpecies,4), RealParts(4), Volume(4)
LOGICAL                       :: ForceNearestNeigh
!===================================================================================================================================
ForceNearestNeigh = .FALSE.
ALLOCATE(iPartIndx_ChildNode(4,TreeNode%PNum_Node))
ALLOCATE(MappedPart_ChildNode(4,TreeNode%PNum_Node,2))
PartNumChildNode(:) = 0
IF (ABS(TreeNode%MidPoint(1)) .EQ. 1.0) THEN
  CALL Abort(&
    __STAMP__,&
    'ERROR in QuadTree Pairing: Too many branches, machine precision reached')
END IF

!         Numbering of the 4 ChildNodes of the QuadTree
!      _________
!     |    |    |   |y
!     | 3  | 2  |   |
!     |----|----|   |
!     | 4  | 1  |   |______x
!     |____|____|

DO iPart=1,TreeNode%PNum_Node
  iPartIndx = TreeNode%iPartIndx_Node(iPart)
  IF ((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(iPart,2).LE.TreeNode%MidPoint(2))) THEN
    PartNumChildNode(1) = PartNumChildNode(1) + 1
    iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    MappedPart_ChildNode(1,PartNumChildNode(1),1:2) = TreeNode%MappedPartStates(iPart,1:2)
  ELSE IF(TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) THEN
    PartNumChildNode(2) = PartNumChildNode(2) + 1
    iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    MappedPart_ChildNode(2,PartNumChildNode(2),1:2) = TreeNode%MappedPartStates(iPart,1:2)
  ELSE IF(TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) THEN
    PartNumChildNode(3) = PartNumChildNode(3) + 1
    iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
    MappedPart_ChildNode(3,PartNumChildNode(3),1:2) = TreeNode%MappedPartStates(iPart,1:2)
  ELSE
    PartNumChildNode(4) = PartNumChildNode(4) + 1
    iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
    MappedPart_ChildNode(4,PartNumChildNode(4),1:2) = TreeNode%MappedPartStates(iPart,1:2)
  END IF
END DO

IF((.NOT.ASSOCIATED(NodeVol)).OR.(.NOT.ASSOCIATED(NodeVol%SubNode1))) THEN
  localDepth = TreeNode%NodeDepth
  CALL DSMC_CalcSubNodeVolumes2D(iElem, localDepth, ElemNodeVol(iElem)%Root)
END IF

NodeVolumeTemp(1) = NodeVol%SubNode1%Volume
NodeVolumeTemp(2) = NodeVol%SubNode2%Volume
NodeVolumeTemp(3) = NodeVol%SubNode3%Volume
NodeVolumeTemp(4) = NodeVol%SubNode4%Volume
FaceVolumeTemp(1) = NodeVol%SubNode1%Area
FaceVolumeTemp(2) = NodeVol%SubNode2%Area
FaceVolumeTemp(3) = NodeVol%SubNode3%Area
FaceVolumeTemp(4) = NodeVol%SubNode4%Area

Volume(1:4) = NodeVolumeTemp(1:4)

IF(DSMC%MergeSubcells) THEN
  IF (ALL(PartNumChildNode.LT.DSMC%PartNumOctreeNodeMin)) THEN
    ForceNearestNeigh =.TRUE.
    DO iLoop = 1, 3
      IF (PartNumChildNode(iLoop).LT.7) THEN
        DO iPart=1, PartNumChildNode(iLoop)
          iPartIndx_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart) = iPartIndx_ChildNode(iLoop,iPart)
          MappedPart_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart,1:2) = MappedPart_ChildNode(iLoop,iPart,1:2)
        END DO
        PartNumChildNode(iLoop+1) = PartNumChildNode(iLoop+1) + PartNumChildNode(iLoop)
        PartNumChildNode(iLoop) = 0
        NodeVolumeTemp(iLoop+1) = NodeVolumeTemp(iLoop+1) + NodeVolumeTemp(iLoop)
        NodeVolumeTemp(iLoop) = 0.0
      END IF
    END DO
    IF (PartNumChildNode(4).LT.7) THEN
      DO iPart=1, PartNumChildNode(4)
        iPartIndx_ChildNode(1,PartNumChildNode(1)+iPart) = iPartIndx_ChildNode(4,iPart)
        MappedPart_ChildNode(1,PartNumChildNode(1)+iPart,1:2) = MappedPart_ChildNode(4,iPart,1:2)
      END DO
      PartNumChildNode(1) = PartNumChildNode(1) + PartNumChildNode(4)
      PartNumChildNode(4) = 0
      NodeVolumeTemp(1) = NodeVolumeTemp(1) + NodeVolumeTemp(4)
      NodeVolumeTemp(4) = 0.0
    END IF
  END IF
END IF

DO iLoop = 1, 4
! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
  ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
  SpecPartNum(:,iLoop) = 0.
  RealParts(iLoop) = 0.
  DO iPart = 1, PartNumChildNode(iLoop)
    RealParts(iLoop) = RealParts(iLoop) + GetParticleWeight(iPartIndx_ChildNode(iLoop,iPart))
    SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart)),iLoop) = &
        SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart)),iLoop) + GetParticleWeight(iPartIndx_ChildNode(iLoop,iPart))
  END DO
END DO

DO iLoop = 1, 4
  ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
  IF((PartNumChildNode(iLoop).GE.DSMC%PartNumOctreeNodeMin).AND.(.NOT.ForceNearestNeigh)) THEN
    ! Additional check if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
    ! the side length of a cube (approximation) with same volume as the actual cell -> octree
    IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop), RealParts(iLoop), Volume(iLoop))
    ELSE
      DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop),REAL(PartNumChildNode(iLoop)), Volume(iLoop))
    END IF
    IF((DSMC%MeanFreePath.LT.(FaceVolumeTemp(iLoop)**(1./2.))).OR.(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode)) THEN
      NULLIFY(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
      ALLOCATE(TreeNode%ChildNode%MappedPartStates(PartNumChildNode(iLoop),1:2))
      TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
      TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
      TreeNode%ChildNode%MappedPartStates(1:PartNumChildNode(iLoop),1:2)= &
            MappedPart_ChildNode(iLoop,1:PartNumChildNode(iLoop),1:2)
      IF (iLoop.LT.3) THEN
        TreeNode%ChildNode%MidPoint(1) = 1.0
        IF (iLoop.EQ.1) THEN
          TreeNode%ChildNode%MidPoint(2) = -1.0
        ELSE
          TreeNode%ChildNode%MidPoint(2) = 1.0
        END IF
      ELSE
        TreeNode%ChildNode%MidPoint(1) = -1.0
        IF (iLoop.EQ.3) THEN
          TreeNode%ChildNode%MidPoint(2) = 1.0
        ELSE
          TreeNode%ChildNode%MidPoint(2) = -1.0
        END IF
      END IF
      TreeNode%ChildNode%MidPoint(3) = 0.0
      TreeNode%ChildNode%MidPoint(1:3) = TreeNode%MidPoint(1:3) &
                                        + TreeNode%ChildNode%MidPoint(1:3)*2.0/(2.0**(TreeNode%NodeDepth+1.0))
      TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
      ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further quadtree division)
      IF (iLoop.EQ.1) CALL AddQuadTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1)
      IF (iLoop.EQ.2) CALL AddQuadTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2)
      IF (iLoop.EQ.3) CALL AddQuadTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3)
      IF (iLoop.EQ.4) CALL AddQuadTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4)
      DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
      DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
      DEALLOCATE(TreeNode%ChildNode)
    ELSE
      IF(DSMC%UseNearestNeighbour) THEN
        CALL FindNearestNeigh2D(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop), TreeNode%MidPoint(1:3), TreeNode%NodeDepth)
      ELSE
        CALL FindStatisticalNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)),PartNumChildNode(iLoop), &
                                  iElem, NodeVolumeTemp(iLoop))
      END IF
    END IF
  ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
    IF(DSMC%UseNearestNeighbour) THEN
      CALL FindNearestNeigh2D(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
            PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop), TreeNode%MidPoint(1:3), TreeNode%NodeDepth)
    ELSE
      CALL FindStatisticalNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)),PartNumChildNode(iLoop), &
                                  iElem, NodeVolumeTemp(iLoop))
    END IF
  ELSE IF (CollInf%ProhibitDoubleColl.AND.(PartNumChildNode(iLoop).EQ.1)) THEN
    CollInf%OldCollPartner(iPartIndx_ChildNode(iLoop, 1)) = 0
  END IF
END DO

END SUBROUTINE AddQuadTreeNode


SUBROUTINE DSMC_init_octree()
!===================================================================================================================================
! Building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : ElemNodeVol
  USE MOD_Mesh_Vars,              ONLY : nElems
  USE MOD_Particle_Vars,          ONLY : Symmetry2D
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
  ALLOCATE(ElemNodeVol(nElems))

  !Calculate recursive Volumes
  DO iElem = 1, nElems
    ALLOCATE(ElemNodeVol(iElem)%Root)
    DO NodeDepth = 1, 2
      IF (Symmetry2D) THEN
        CALL DSMC_CalcSubNodeVolumes2D(iElem, NodeDepth, ElemNodeVol(iElem)%Root)
      ELSE
        CALL DSMC_CalcSubNodeVolumes(iElem, NodeDepth, ElemNodeVol(iElem)%Root)
      END IF
    END DO
  END DO

END SUBROUTINE DSMC_init_octree

SUBROUTINE DSMC_CalcSubNodeVolumes2D(iElem, NodeDepth, Node)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : OctreeVdm, tNodeVolume, SymmetrySide
USE MOD_Mesh_Vars,              ONLY : SurfElem, Face_xGP
USE MOD_PreProc,                ONLY : PP_N
USE MOD_ChangeBasis,            ONLY : ChangeBasis2D
USE MOD_PreProc,                ONLY : PP_N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                       :: iElem
INTEGER, INTENT(INOUT)                    :: NodeDepth
TYPE (tNodeVolume), INTENT(OUT), POINTER  :: Node
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: j, k , NumOfPoints,SideID
REAL                                    :: DetLocal(1,0:PP_N,0:PP_N)
REAL                                    :: FaceLocal(2,0:PP_N,0:PP_N)
REAL, ALLOCATABLE                       :: DetJac(:,:,:)
REAL, ALLOCATABLE                       :: FacexGP(:,:,:)
REAL, ALLOCATABLE                       :: LocalVdm(:,:)
INTEGER                                 :: LocalDepth
!===================================================================================================================================
  LocalDepth = 1
  NumOfPoints = 2**NodeDepth
  SideID = SymmetrySide(iElem,1)

  DO j=0, PP_N; DO k=0, PP_N
    DetLocal(1,j,k)=SurfElem(j,k,SideID)
  END DO; END DO

  DO j=0, PP_N; DO k=0, PP_N
    FaceLocal(1:2,j,k) = Face_xGP(1:2,j,k,SideID)
  END DO; END DO

  ALLOCATE( DetJac(1,0:NumOfPoints - 1,0:NumOfPoints - 1))
  ALLOCATE(LocalVdm(0:NumOfPoints - 1,0:PP_N))
  ALLOCATE(FacexGP(2,0:NumOfPoints - 1,0:NumOfPoints - 1))
  CALL InitVanderOct(LocalVdm, NodeDepth, LocalDepth, OctreeVdm)
  CALL ChangeBasis2D(1, PP_N, NumOfPoints - 1, LocalVdm ,DetLocal(:,:,:),DetJac(:,:,:))
  CALL ChangeBasis2D(2, PP_N, NumOfPoints - 1, LocalVdm ,FaceLocal(:,:,:),FacexGP(:,:,:))
  CALL AddNodeVolumes2D(NodeDepth, Node, DetJac, OctreeVdm, iElem, FacexGP)

END SUBROUTINE DSMC_CalcSubNodeVolumes2D


SUBROUTINE DSMC_CalcSubNodeVolumes(iElem, NodeDepth, Node)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : OctreeVdm, tNodeVolume
  USE MOD_Mesh_Vars,              ONLY : sJ
  USE MOD_PreProc,                ONLY : PP_N
  USE MOD_ChangeBasis,            ONLY : ChangeBasis3D
  USE MOD_PreProc,                ONLY : PP_N
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                       :: iElem
  INTEGER, INTENT(INOUT)                    :: NodeDepth
  TYPE (tNodeVolume), INTENT(OUT), POINTER  :: Node
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                 :: j, k ,l, NumOfPoints
  REAL                                    :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N)
  REAL, ALLOCATABLE                       :: DetJac(:,:,:,:)
  REAL, ALLOCATABLE                       :: LocalVdm(:,:)
  INTEGER                                 :: LocalDepth
!===================================================================================================================================
  LocalDepth = 1
  NumOfPoints = 2**NodeDepth

  DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
    DetLocal(1,j,k,l)=1./sJ(j,k,l,iElem)
  END DO; END DO; END DO

  ALLOCATE( DetJac(1,0:NumOfPoints - 1,0:NumOfPoints - 1,0:NumOfPoints - 1))
  ALLOCATE(LocalVdm(0:NumOfPoints - 1,0:PP_N))
  CALL InitVanderOct(LocalVdm, NodeDepth, LocalDepth, OctreeVdm)
  CALL ChangeBasis3D(1,PP_N, NumOfPoints - 1, LocalVdm, DetLocal(:,:,:,:),DetJac(:,:,:,:))
  CALL AddNodeVolumes(NodeDepth, Node, DetJac, OctreeVdm)

END SUBROUTINE DSMC_CalcSubNodeVolumes

RECURSIVE SUBROUTINE InitVanderOct(LocalVdm, NodeDepth, LocalDepth, OctreeVdmLoc)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : tOctreeVdm
  USE MOD_PreProc,                ONLY : PP_N
  USE MOD_Interpolation_Vars,     ONLY : xGP, wBary
  USE MOD_Basis,                  ONLY : InitializeVandermonde
  USE MOD_PreProc,                ONLY : PP_N
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(INOUT)                    :: NodeDepth, LocalDepth
  REAL, INTENT(OUT)                         :: LocalVdm(0:2**NodeDepth-1,0:PP_N)
  TYPE (tOctreeVdm), POINTER, INTENT(INOUT) :: OctreeVdmLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                   :: i, NodePointNum
!===================================================================================================================================
  IF (.NOT.ASSOCIATED(OctreeVdmLoc)) THEN
    NULLIFY(OctreeVdmLoc)
    ALLOCATE(OctreeVdmLoc)
    NodePointNum = 2**LocalDepth - 1
    ALLOCATE(OctreeVdmLoc%Vdm(0:NodePointNum, 0:PP_N), OctreeVdmLoc%xGP(0:NodePointNum))
    DO i=0,NodePointNum
      OctreeVdmLoc%xGP(i) = -1.0 + 2./(1.+NodePointNum) * ((REAL(i)+1.) - 0.5)
    END DO
    OctreeVdmLoc%wGP = 2./REAL(1.0+NodePointNum)
    CALL InitializeVandermonde(PP_N,NodePointNum,wBary,xGP,OctreeVdmLoc%xGP,OctreeVdmLoc%Vdm)
  END IF
  IF (LocalDepth.EQ.NodeDepth) THEN
    LocalVdm = OctreeVdmLoc%Vdm
  ELSE
    LocalDepth = LocalDepth + 1
    CALL InitVanderOct(LocalVdm, NodeDepth, LocalDepth, OctreeVdmLoc%SubVdm)
  END IF
END SUBROUTINE InitVanderOct

RECURSIVE SUBROUTINE AddNodeVolumes(NodeDepth, Node, DetJac, VdmLocal, SubNodesIn)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : tOctreeVdm, tNodeVolume
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(INOUT)                    :: NodeDepth
  INTEGER, INTENT(IN), OPTIONAL             :: SubNodesIn(:)
  TYPE (tNodeVolume), INTENT(OUT), POINTER  :: Node
  REAL, INTENT(INOUT)                       :: DetJac(1,0:2**NodeDepth-1,0:2**NodeDepth-1,0:2**NodeDepth-1)
  TYPE (tOctreeVdm), INTENT(OUT), POINTER   :: VdmLocal
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER, ALLOCATABLE                       :: SubNodesOut(:)
  INTEGER                                    :: OldNodeNum, NewNodeNum, j, VolPos(3)
!===================================================================================================================================
IF (PRESENT(SubNodesIn)) THEN
  OldNodeNum = SIZE(SubNodesIn)
  NewNodeNum = OldNodeNum + 1
  ALLOCATE(SubNodesOut(NewNodeNum))
  SubNodesOut(1:OldNodeNum) = SubNodesIn(1:OldNodeNum)
ELSE
  OldNodeNum = 0
  NewNodeNum = OldNodeNum + 1
  ALLOCATE(SubNodesOut(NewNodeNum))
END IF

IF (OldNodeNum.NE.NodeDepth) THEN
  IF (OldNodeNum.EQ.0) THEN
    SubNodesOut(NewNodeNum) = 1
    IF(.NOT.ASSOCIATED(Node%SubNode1)) ALLOCATE(Node%SubNode1)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode1, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 2
    IF(.NOT.ASSOCIATED(Node%SubNode2)) ALLOCATE(Node%SubNode2)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode2, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 3
    IF(.NOT.ASSOCIATED(Node%SubNode3)) ALLOCATE(Node%SubNode3)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode3, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 4
    IF(.NOT.ASSOCIATED(Node%SubNode4)) ALLOCATE(Node%SubNode4)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode4, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 5
    IF(.NOT.ASSOCIATED(Node%SubNode5)) ALLOCATE(Node%SubNode5)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode5, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 6
    IF(.NOT.ASSOCIATED(Node%SubNode6)) ALLOCATE(Node%SubNode6)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode6, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 7
    IF(.NOT.ASSOCIATED(Node%SubNode7)) ALLOCATE(Node%SubNode7)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode7, DetJac, VdmLocal, SubNodesOut)
    SubNodesOut(NewNodeNum) = 8
    IF(.NOT.ASSOCIATED(Node%SubNode8)) ALLOCATE(Node%SubNode8)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode8, DetJac, VdmLocal, SubNodesOut)
  ELSE
    SubNodesOut(NewNodeNum) = 1
    IF(.NOT.ASSOCIATED(Node%SubNode1)) ALLOCATE(Node%SubNode1)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode1, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 2
    IF(.NOT.ASSOCIATED(Node%SubNode2)) ALLOCATE(Node%SubNode2)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode2, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 3
    IF(.NOT.ASSOCIATED(Node%SubNode3)) ALLOCATE(Node%SubNode3)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode3, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 4
    IF(.NOT.ASSOCIATED(Node%SubNode4)) ALLOCATE(Node%SubNode4)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode4, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 5
    IF(.NOT.ASSOCIATED(Node%SubNode5)) ALLOCATE(Node%SubNode5)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode5, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 6
    IF(.NOT.ASSOCIATED(Node%SubNode6)) ALLOCATE(Node%SubNode6)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode6, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 7
    IF(.NOT.ASSOCIATED(Node%SubNode7)) ALLOCATE(Node%SubNode7)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode7, DetJac, VdmLocal%SubVdm, SubNodesOut)
    SubNodesOut(NewNodeNum) = 8
    IF(.NOT.ASSOCIATED(Node%SubNode8)) ALLOCATE(Node%SubNode8)
    CALL AddNodeVolumes(NodeDepth, Node%SubNode8, DetJac, VdmLocal%SubVdm, SubNodesOut)
   END IF
ELSE
  VolPos(:) = 0
  DO j=1, NodeDepth
    ! x direction
    IF (SubNodesOut(j).LT.5) THEN
      VolPos(1)= VolPos(1) + 2**NodeDepth/2**j
    END IF
    ! y direction
    IF ((SubNodesOut(j).LT.3).OR.((SubNodesOut(j).LT.7).AND.(SubNodesOut(j).GT.4))) THEN
      VolPos(2)= VolPos(2) + 2**NodeDepth/2**j
    END IF
    ! z direction
    IF ((SubNodesOut(j).EQ.2).OR.(SubNodesOut(j).EQ.3).OR.(SubNodesOut(j).EQ.6).OR.(SubNodesOut(j).EQ.7)) THEN
      VolPos(3)= VolPos(3) + 2**NodeDepth/2**j
    END IF
  END DO
  Node%Volume = DetJac(1, VolPos(1), VolPos(2), VolPos(3)) * VdmLocal%wGP**3
END IF

END SUBROUTINE AddNodeVolumes


RECURSIVE SUBROUTINE AddNodeVolumes2D(NodeDepth, Node, DetJac, VdmLocal, iElem, FacexGP, SubNodesIn)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars    ,ONLY: Pi
USE MOD_DSMC_Vars       ,ONLY: tOctreeVdm, tNodeVolume
USE MOD_Particle_Vars   ,ONLY: Symmetry2DAxisymmetric
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)                    :: NodeDepth
INTEGER, INTENT(IN)                       :: iElem
INTEGER, INTENT(IN), OPTIONAL             :: SubNodesIn(:)
TYPE (tNodeVolume), INTENT(OUT), POINTER  :: Node
REAL, INTENT(INOUT)                       :: DetJac(1,0:2**NodeDepth-1,0:2**NodeDepth-1)
REAL, INTENT(INOUT)                       :: FacexGP(2,0:2**NodeDepth-1,0:2**NodeDepth-1)
TYPE (tOctreeVdm), INTENT(OUT), POINTER   :: VdmLocal
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE                      :: SubNodesOut(:)
INTEGER                                   :: OldNodeNum, NewNodeNum, j, VolPos(2), dirX, dirY, realDirX, realDirY
REAL                                      :: Pos(2)
REAL                                      :: Distance, TempDistance
!===================================================================================================================================
IF (PRESENT(SubNodesIn)) THEN
  OldNodeNum = SIZE(SubNodesIn)
  NewNodeNum = OldNodeNum + 1
  ALLOCATE(SubNodesOut(NewNodeNum))
  SubNodesOut(1:OldNodeNum) = SubNodesIn(1:OldNodeNum)
ELSE
  OldNodeNum = 0
  NewNodeNum = OldNodeNum + 1
  ALLOCATE(SubNodesOut(NewNodeNum))
END IF

IF (OldNodeNum.NE.NodeDepth) THEN
  IF (OldNodeNum.EQ.0) THEN
    SubNodesOut(NewNodeNum) = 1
    IF(.NOT.ASSOCIATED(Node%SubNode1)) ALLOCATE(Node%SubNode1)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode1, DetJac, VdmLocal, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 2
    IF(.NOT.ASSOCIATED(Node%SubNode2)) ALLOCATE(Node%SubNode2)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode2, DetJac, VdmLocal, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 3
    IF(.NOT.ASSOCIATED(Node%SubNode3)) ALLOCATE(Node%SubNode3)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode3, DetJac, VdmLocal, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 4
    IF(.NOT.ASSOCIATED(Node%SubNode4)) ALLOCATE(Node%SubNode4)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode4, DetJac, VdmLocal, iElem, FacexGP, SubNodesOut)
  ELSE
    SubNodesOut(NewNodeNum) = 1
    IF(.NOT.ASSOCIATED(Node%SubNode1)) ALLOCATE(Node%SubNode1)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode1, DetJac, VdmLocal%SubVdm, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 2
    IF(.NOT.ASSOCIATED(Node%SubNode2)) ALLOCATE(Node%SubNode2)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode2, DetJac, VdmLocal%SubVdm, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 3
    IF(.NOT.ASSOCIATED(Node%SubNode3)) ALLOCATE(Node%SubNode3)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode3, DetJac, VdmLocal%SubVdm, iElem, FacexGP, SubNodesOut)
    SubNodesOut(NewNodeNum) = 4
    IF(.NOT.ASSOCIATED(Node%SubNode4)) ALLOCATE(Node%SubNode4)
    CALL AddNodeVolumes2D(NodeDepth, Node%SubNode4, DetJac, VdmLocal%SubVdm, iElem, FacexGP, SubNodesOut)
   END IF
ELSE
  VolPos(:) = 0
  DO j=1, NodeDepth
    ! x direction
    IF (SubNodesOut(j).LT.3) THEN
      VolPos(1)= VolPos(1) + 2**NodeDepth/2**j
    END IF
    ! y direction
    IF ((SubNodesOut(j).LT.4).AND.(SubNodesOut(j).GT.1)) THEN
      VolPos(2)= VolPos(2) + 2**NodeDepth/2**j
    END IF
  END DO
  Pos(1) = VdmLocal%xGP(VolPos(1))
  Pos(2) = VdmLocal%xGP(VolPos(2))
  Pos(1:2) = MapToGeo2D(Pos, iElem)
  Distance = HUGE(Distance)
  TempDistance = 0.0
  realDirX = 0
  realDirY = 0
  DO dirY = 0,2**NodeDepth-1; DO dirX = 0,2**NodeDepth-1
    TempDistance = (FacexGP(1,dirX,dirY)-Pos(1))**2 + (FacexGP(2,dirX,dirY)-Pos(2))**2
    IF(Distance.GT.TempDistance) THEN
      Distance = TempDistance
      realDirX = dirX
      realDirY = dirY
    END IF
  END DO; END DO
  IF (Symmetry2DAxisymmetric) THEN
    Node%Volume = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2 * 2. * Pi * Pos(2)
  ELSE
    Node%Volume = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2
  END IF
  Node%Area = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2
END IF

END SUBROUTINE AddNodeVolumes2D


SUBROUTINE FindStatisticalNeigh(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
! Classic statistical pairing method for the use in the octree routines
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Relaxation       ,ONLY: SetMeanVibQua
USE MOD_DSMC_CollisionProb    ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Collis           ,ONLY: DSMC_perform_collision
USE MOD_vmpf_collision        ,ONLY: DSMC_vmpf_prob
USE MOD_DSMC_Vars             ,ONLY: Coll_pData,CollInf,CollisMode,PartStateIntEn,ChemReac,PairE_vMPF,BGGas,DSMC,RadialWeighting
USE MOD_DSMC_Vars             ,ONLY: SamplingActive, SelectionProc, SpecDSMC
USE MOD_Particle_Vars         ,ONLY: PartSpecies, nSpecies, PartState, usevMPF, PartMPF, WriteMacroVolumeValues, VarTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: TEnd, time
USE MOD_DSMC_Analyze          ,ONLY: CalcGammaVib, CalcInstantTransTemp, CalcMeanFreePath
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(IN)                     :: PartNum
INTEGER, INTENT(IN)                     :: iElem
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nPair, iPair, iPart, cPart1, cPart2, nPart, tempPart
INTEGER                       :: cSpec1, cSpec2, iCase
REAL                          :: iRan
REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

nPart = PartNum
nPair = INT(nPart/2)
IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) CollInf%MeanMPF = 0.

IF (CollisMode.EQ.3) THEN
  ChemReac%RecombParticle = 0
  ChemReac%nPairForRec = 0
END IF

CollInf%Coll_SpecPartNum = 0
CollInf%Coll_CaseNum = 0

ALLOCATE(Coll_pData(nPair))
Coll_pData%Ec=0

IF (CollisMode.EQ.3) THEN
! Determination of the mean vibrational energy for the cell, only needed for chemical reactions
  ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
  DO iPart = 1, PartNum
    ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) = ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) &
      + PartStateIntEn(iPartIndx_Node(iPart),1) * GetParticleWeight(iPartIndx_Node(iPart))
  END DO
END IF

DO iPart = 1, PartNum
  CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) &
                                                                  + GetParticleWeight(iPartIndx_Node(iPart))
END DO

IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
  ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
  !           relaxation probability (CalcGammaVib)
  ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
  ! 3. Case: Temperature required for the mean free path with the VHS model
  CALL CalcInstantTransTemp(iPartIndx_Node,PartNum)
  IF(SelectionProc.EQ.2) CALL CalcGammaVib()
END IF

IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) MPFFac = 1

IF (CollInf%ProhibitDoubleColl.AND.(nPair.EQ.1)) THEN
! Do not get stuck in an endless loop if only two particles/one pair are present in the cell
  CollInf%OldCollPartner(iPartIndx_Node(1)) = 0
  CollInf%OldCollPartner(iPartIndx_Node(2)) = 0
END IF

DO iPair = 1, nPair                               ! statistical pairing
  CALL RANDOM_NUMBER(iRan)
  cPart1 = 1 + INT(nPart * iRan)                       ! first pair particle
  Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(cPart1)
  iPartIndx_Node(cPart1) = iPartIndx_Node(nPart)
  nPart = nPart - 1

  CALL RANDOM_NUMBER(iRan)
  cPart2 = 1 + INT(nPart * iRan)                       ! second pair particle

  IF(CollInf%ProhibitDoubleColl) THEN
    IF(iPartIndx_Node(cPart2).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
      IF(iPair.LT.nPair) THEN
        DO WHILE (iPartIndx_Node(cPart2).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1))
          CALL RANDOM_NUMBER(iRan)
          cPart2 = 1 + INT(nPart * iRan)
        END DO
!        ELSE ! Last pair did collide last iteration, might require additional treatment, split up another pair
      END IF
    END IF
  END IF

  Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(cPart2)
  iPartIndx_Node(cPart2) = iPartIndx_Node(nPart)
  nPart = nPart - 1

  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

  IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
    TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
    IF (TempMPFFac .GE. MPFFac) THEN
        MPFFac = TempMPFFac
        PairE_vMPF(1) = iPair
      IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
        PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
      ELSE
        PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
      END IF
    END IF
  END IF

  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  ! Summation of the average weighting factor of the collision pairs for each case (AA, AB, BB)
  IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(Coll_pData(iPair)%iPart_p1) &
                                                        + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
  END IF

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
IF ((nPair.NE.0).AND.(CollisMode.EQ.3).AND.(MOD(nPart, nPair).NE.0)) THEN
  ChemReac%RecombParticle = iPartIndx_Node(1)
END IF

! Resetting the previous collision partner of the remaining particle due to uneven PartNum
IF (CollInf%ProhibitDoubleColl.AND.(nPart.EQ.1)) CollInf%OldCollPartner(iPartIndx_Node(1)) = 0

IF (CollisMode.EQ.3) THEN
  CALL SetMeanVibQua()
END IF

IF(RadialWeighting%DoRadialWeighting) THEN
  DO iPair = 1, nPair
    ! Two particles with the exact same velocities at the same positions -> clones that did not interact with other particles/walls
    IF (Coll_pData(iPair)%CRela2.EQ.0.0) THEN
      IF(DSMC%CalcQualityFactors.AND.(SamplingActive.OR.WriteMacroVolumeValues)) THEN
        DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
      END IF
      ! "Partner-Tausch": if there are pairs ahead in the pairing list, the next is pair is broken up and collision partners
      ! are swapped but first, changing z-direction
      IF (nPart.EQ.1) THEN
        ! Uneven number of particles in the cell, a single particle is left without a pair
        PartState(Coll_pData(iPair)%iPart_p1,6) = - PartState(Coll_pData(iPair)%iPart_p1,6)
        ! Removing the pairs from the weighting factor and the case num sums
        CollInf%MeanMPF(Coll_pData(iPair)%PairType) = CollInf%MeanMPF(Coll_pData(iPair)%PairType) &
          -(GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
        CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
        ! Swapping particle without a pair with the first particle of the current pair
        tempPart = Coll_pData(iPair)%iPart_p1
        Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(1)
        iPartIndx_Node(1) = tempPart
        IF (CollisMode.EQ.3) ChemReac%RecombParticle = iPartIndx_Node(1)
        IF (CollInf%ProhibitDoubleColl)  CollInf%OldCollPartner(iPartIndx_Node(1)) = 0
        ! Calculation of the relative velocity for the new first pair
        cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
        cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
        iCase = CollInf%Coll_Case(cSpec1, cSpec2)
        ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
        CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(cSpec1) + GetParticleWeight(cSpec2))*0.5
        CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
        Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                                  + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                                  + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
        Coll_pData(iPair)%PairType = iCase
      ELSE IF (iPair.LT.nPair) THEN
        PartState(Coll_pData(iPair)%iPart_p1,6) = - PartState(Coll_pData(iPair)%iPart_p1,6)
        ! Removing the pairs from the weighting factor and the case num sums
        CollInf%MeanMPF(Coll_pData(iPair)%PairType) = CollInf%MeanMPF(Coll_pData(iPair)%PairType) &
          -(GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
        CollInf%MeanMPF(Coll_pData(iPair+1)%PairType) = CollInf%MeanMPF(Coll_pData(iPair+1)%PairType) &
          - (GetParticleWeight(Coll_pData(iPair+1)%iPart_p1) + GetParticleWeight(Coll_pData(iPair+1)%iPart_p2))*0.5
        CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
        CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) - 1
        ! Breaking up the next pair and swapping partners
        tempPart = Coll_pData(iPair)%iPart_p1
        Coll_pData(iPair)%iPart_p1 = Coll_pData(iPair + 1)%iPart_p1
        Coll_pData(iPair + 1)%iPart_p1 = tempPart
        ! Calculation of the relative velocity for the new first pair
        cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
        cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
        iCase = CollInf%Coll_Case(cSpec1, cSpec2)
        ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
        CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(cSpec1) + GetParticleWeight(cSpec2))*0.5
        CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
        Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                                  + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                                  + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                                  -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
        Coll_pData(iPair)%PairType = iCase
        ! Calculation of the relative velocity for the new follow-up pair
        cSpec1 = PartSpecies(Coll_pData(iPair+1)%iPart_p1)
        cSpec2 = PartSpecies(Coll_pData(iPair+1)%iPart_p2)
        iCase = CollInf%Coll_Case(cSpec1, cSpec2)
        ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
        CollInf%MeanMPF(iCase) = CollInf%MeanMPF(iCase) + (GetParticleWeight(cSpec1) + GetParticleWeight(cSpec2))*0.5
        CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
        Coll_pData(iPair+1)%CRela2 = (PartState(Coll_pData(iPair+1)%iPart_p1,4) &
                                    -  PartState(Coll_pData(iPair+1)%iPart_p2,4))**2 &
                                    + (PartState(Coll_pData(iPair+1)%iPart_p1,5) &
                                    -  PartState(Coll_pData(iPair+1)%iPart_p2,5))**2 &
                                    + (PartState(Coll_pData(iPair+1)%iPart_p1,6) &
                                    -  PartState(Coll_pData(iPair+1)%iPart_p2,6))**2
        Coll_pData(iPair+1)%PairType = iCase
      END IF
    END IF    ! Coll_pData(iPair)%CRela2.EQ.0.0
  END DO      ! DO iPair=1, nPair
END IF        ! DoRadialWeighting

DO iPair = 1, nPair
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    IF (usevMPF.AND.(BGGas%BGGasSpecies.EQ.0).AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
      CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
    ELSE
      CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
      CALL DSMC_perform_collision(iPair,iElem, NodeVolume, PartNum)
      IF (CollInf%ProhibitDoubleColl) THEN
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1) = Coll_pData(iPair)%iPart_p2
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p2) = Coll_pData(iPair)%iPart_p1
      END IF
    ELSE
      IF (CollInf%ProhibitDoubleColl) THEN
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1) = 0
        CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p2) = 0
      END IF
    END IF
  END IF
END DO

IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path with VHS model and the current translational temperature in the cell
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                                          SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the maximum MCS/MFP for the cell
    IF((DSMC%CollSepCount.GT.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = &
                                                    MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
  END IF
END IF

DEALLOCATE(Coll_pData)

END SUBROUTINE FindStatisticalNeigh


FUNCTION Calc_F2D(xi,x,P)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: xi(2)
REAL,INTENT(IN)          :: x(2)
REAL,INTENT(IN)          :: P(2,4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_F2D(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Calc_F2D = 0.25 *(P(:,1)*(1-xi(1)) * (1-xi(2)) &
              + P(:,2)*(1+xi(1)) * (1-xi(2)) &
              + P(:,3)*(1+xi(1)) * (1+xi(2)) &
              + P(:,4)*(1-xi(1)) * (1+xi(2))) &
              - x;
END FUNCTION Calc_F2D


FUNCTION Calc_dF_inv2D(xi,P)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
  USE MOD_Globals,        ONLY : Abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: xi(2)
REAL,INTENT(IN)          :: P(2,4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_dF_inv2D(2,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: dF(2,2)
REAL                     :: dF_inv(2,2)
REAL                     :: detjb
!===================================================================================================================================
dF = 0.
dF(:,1)= 0.25 * ((P(:,2)-P(:,1))*(1-xi(2))+(P(:,3)-P(:,4))*(1+xi(2)))
dF(:,2)= 0.25 * ((P(:,4)-P(:,1))*(1-xi(1))+(P(:,3)-P(:,2))*(1+xi(1)))

! Determines the determinant of xj and checks for zero values
!
detjb = dF(1, 1) * dF(2, 2) - dF(1, 2) * dF(2, 1)

IF ( detjb <= 0.d0 ) then
  WRITE(*,*)"Negative determinant of Jacobian in calc_df_inv:"
  WRITE(*,*)"dF",dF
  WRITE(*,*)"Determinant is:",detjb
  CALL Abort(&
       __STAMP__,&
      'Negative determinant of Jacobian in calc_df_inv')
END IF
!
! Determines the inverse of xj
!
dF_inv(1, 1) = dF(2, 2) / detjb
dF_inv(1, 2) = - dF(1, 2) / detjb
dF_inv(2, 1) = - dF(2, 1) / detjb
dF_inv(2, 2) = dF(1, 1) / detjb

Calc_dF_inv2D = dF_inv

END FUNCTION Calc_dF_inv2D


FUNCTION Calc_inv2D(M)
!===================================================================================================================================
!> calc inverse of M
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: M(2,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_inv2D(2,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: M_inv(2,2)
REAL                     :: detjb
!===================================================================================================================================
M_inv = 0.
! Determines the determinant of xj and checks for zero values
detjb = M (1, 1) * M (2, 2) - M (1, 2) * M (2, 1)
IF ( detjb == 0.d0 ) then
  IPWRITE(UNIT_errOut,*)"Determinant is:",detjb
  IPWRITE(UNIT_errOut,*)"KM:",M_inv
  CALL abort(__STAMP__, &
        "Zero determinant of Jacobian in M_inv")
END IF
! Determines the inverse of xj
M_inv (1, 1) = M (2, 2)/ detjb
M_inv (1, 2) = - M (1, 2) / detjb
M_inv (2, 1) = - M (2, 1) / detjb
M_inv (2, 2) = M (1, 1) / detjb
Calc_inv2D = M_inv
END FUNCTION Calc_inv2D


SUBROUTINE GeoCoordToMap2D(x_in,xi_Out,iElem)
!===================================================================================================================================
!> interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
!> first get xi,eta,zeta from x,y,z...then do tenso product interpolation
!> xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_DSMC_Vars             ,ONLY: SymmetrySide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem                                 !< elem index
REAL,INTENT(IN)               :: x_in(2)                               !< physical position of particle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: xi_Out(2)  ! Interpolated Pos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,j,k, iNode, SideID
REAL                          :: xi(2)
REAL                          :: P(2,4), F(2), dF_inv(2,2), s(2)
REAL, PARAMETER               :: EPS=1E-10
REAL                          :: T_inv(2,2), DP(2), T(2,2)
!===================================================================================================================================
! --------------------------------------------------
! 1.) Mapping: get xi,eta,zeta value from x,y,z
! --------------------------------------------------
! 1.1.) initial guess from linear part:
SideID = SymmetrySide(iElem,2)
DO iNode = 1,4
  P(1:2,iNode) = GEO%NodeCoords(1:2,GEO%ElemSideNodeID(iNode,SideID,iElem))
END DO
T(:,1) = 0.5 * (P(:,2)-P(:,1))
T(:,2) = 0.5 * (P(:,4)-P(:,1))
T_inv = Calc_inv2D(T)

! transform also the physical coordinate of the point into the unit element (this is the solution of the linear problem already)
xi = 0.
DP = x_in - P(:,1)
DO i=1,2
  DO j=1,2
    xi(i)= xi(i) + T_inv(i,j) * DP(j)
  END DO
END DO

IF ((xi(1).GE.0.0.AND.xi(1).LE.2.0).AND.(xi(2).GE.0.0.AND.xi(2).LE.2.0)) THEN
  xi = xi - (/1.,1./)
ELSE
  xi = (/0.,0./)
END IF

! 1.2.) Newton-Method to solve non-linear part (if linear elements then F should become 0 and no Newton step is required).

F = Calc_F2D(xi,x_in,P)
DO WHILE(SUM(ABS(F)).GE.EPS)
  dF_inv = Calc_dF_inv2D(xi,P)
  s=0.
  DO j = 1,2
    DO k = 1,2
      s(j) = s(j) + dF_inv(j,k) * F(k)
    END DO ! k
  END DO ! j
  xi = xi - s
  F = Calc_F2D(xi,x_in,P)
END DO ! i
IF ((xi(1).GE.-1.0.AND.xi(1).LE.1.0).AND.(xi(2).GE.-1.0.AND.xi(2).LE.1.0)) THEN
  xi_Out = xi
ELSE IF ((xi(1).LE.-1.0)) THEN
  xi_Out = xi
  xi_Out(1) = -0.9999999999999
ELSE IF ((xi(1).GE.1.0)) THEN
  xi_Out = xi
  xi_Out(1) = 0.9999999999999
ELSE IF ((xi(2).LE.-1.0)) THEN
  xi_Out = xi
  xi_Out(2) = -0.9999999999999
ELSE IF ((xi(2).GE.1.0)) THEN
  xi_Out = xi
  xi_Out(2) = 0.9999999999999
END IF

END SUBROUTINE GeoCoordToMap2D


FUNCTION MapToGeo2D(xi,iElem)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_DSMC_Vars               ,ONLY: SymmetrySide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: xi(2)
INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, iNode
REAL                            :: MapToGeo2D(2),P(2,4)
!===================================================================================================================================
SideID = SymmetrySide(iElem,2)
DO iNode = 1,4
  P(1:2,iNode) = GEO%NodeCoords(1:2,GEO%ElemSideNodeID(iNode,SideID,iElem))
END DO

MapToGeo2D =0.25*(P(:,1)*(1-xi(1)) * (1-xi(2)) &
              + P(:,2)*(1+xi(1)) * (1-xi(2)) &
              + P(:,3)*(1+xi(1)) * (1+xi(2)) &
              + P(:,4)*(1-xi(1)) * (1+xi(2)))

END FUNCTION MapToGeo2D

END MODULE MOD_DSMC_ParticlePairing
