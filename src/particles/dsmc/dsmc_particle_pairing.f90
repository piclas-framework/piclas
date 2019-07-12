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
PUBLIC :: DSMC_pairing_statistical, DSMC_pairing_octree, DSMC_init_octree, DSMC_CalcSubNodeVolumes
!===================================================================================================================================

CONTAINS


SUBROUTINE FindNearestNeigh(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : CollInf, tTreeNode, CollisMode, ChemReac, PartStateIntEn, Coll_pData, SelectionProc
  USE MOD_DSMC_Vars,              ONLY : DSMC, PairE_vMPF, SpecDSMC
  USE MOD_Particle_Vars,          ONLY : PartState, nSpecies, PartSpecies, usevMPF, PartMPF, WriteMacroVolumeValues
  USE MOD_DSMC_Relaxation,        ONLY : SetMeanVibQua
  USE MOD_DSMC_Analyze,           ONLY : CalcGammaVib, CalcInstantTransTemp, CalcMeanFreePath
  USE MOD_Particle_Analyze_Vars,  ONLY : CalcEkin
  USE MOD_DSMC_CollisionProb,     ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,            ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,         ONLY : DSMC_vmpf_prob
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
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
  REAL                          :: Dist1, Dist2, iRan
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

  PairNum_Node = INT(PartNum/2)
  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
    ChemReac%nPairForRec = 0
  END IF

  IF (CollisMode.EQ.3) THEN
    ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
    DO iPart = 1, PartNum
      ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) = &
        ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) + PartStateIntEn(iPartIndx_Node(iPart),1)
    END DO
  END IF

  DO iPart = 1, PartNum
    IF(usevMPF) THEN
      CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
                CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + PartMPF(iPartIndx_Node(iPart))
    ELSE
      CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
                CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + 1
    END IF
  END DO

  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
    ! 3. Case: Temperature required for the mean free path with the VHS model
    CALL CalcInstantTransTemp(iPartIndx_Node,PartNum)
    IF(SelectionProc.EQ.2) CALL CalcGammaVib()
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
    CollInf%Coll_CaseNum(iCase)  = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2     = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                                 -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                                 + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                                 -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                                 + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                                 -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
    Coll_pData(iPair)%PairType   = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
  END DO
  IF ((PairNum_Node.NE.0).AND.(CollisMode.EQ.3).AND.(MOD(PartNum, PairNum_Node).NE.0)) THEN
    ChemReac%RecombParticle = iPartIndx_Node(1)
  END IF

  IF (CollisMode.EQ.3) THEN
    CALL SetMeanVibQua()
  END IF

  DO iPair = 1, PairNum_Node
    IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
      IF (usevMPF) THEN            ! calculation of collision prob
        CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
      ELSE
        CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
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
      IF(CollInf%collModel.EQ.0) THEN
        ! Calculation of the mean free path with VHS model and the current translational temperature in the cell
        DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                            SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
      ELSE ! Calculation of mean free path with VSS model accordingly
        DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                            opt_temp=DSMC%InstantTransTemp(nSpecies+1))
      END IF
      ! Determination of the maximum MCS/MFP for the cell
      IF (DSMC%CollSepCount.GT.0 .AND. DSMC%MeanFreePath.GT.0.0) THEN
        DSMC%MCSoverMFP = MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
      END IF
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
  USE MOD_DSMC_Vars,              ONLY : DSMC, SelectionProc
  USE MOD_DSMC_Analyze,           ONLY : CalcGammaVib, CalcInstantTransTemp
  USE MOD_Particle_Vars,          ONLY : PEM, PartSpecies, nSpecies, PartState, usevMPF, PartMPF
  USE MOD_Particle_Vars,          ONLY : KeepWallParticles, PDM
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

  CollInf%Coll_SpecPartNum = 0
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
    IF(usevMPF) THEN
      CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + PartMPF(iPart)
    ELSE
      CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + 1
    END IF
    ! counter for part num of spec per cell
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) = &
                          ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) &
                        + PartStateIntEn(iPart,1) !Calculation of mean evib per cell and iter, necessary for disso prob
    ! choose next particle in Element
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

  IF (usevMPF) MPFFac = 1

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


SUBROUTINE DSMC_pairing_octree(iElem)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Analyze,           ONLY : CalcMeanFreePath
  USE MOD_DSMC_Vars,              ONLY : tTreeNode, DSMC, ElemNodeVol
  USE MOD_Particle_Vars,          ONLY : PEM, PartState, nSpecies, PartSpecies,PartPosRef
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_Particle_Tracking_vars, ONLY : DoRefMapping
  USE MOD_Eval_xyz,               ONLY : GetPositionInRefElem
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, nPart, SpecPartNum(nSpecies)
  TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

SpecPartNum = 0
nPart = PEM%pNumber(iElem)

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
              SpecPartNum(PartSpecies(TreeNode%iPartIndx_Node(iLoop))) + 1
  END DO

  DSMC%MeanFreePath = CalcMeanFreePath(REAL(SpecPartNum), REAL(nPart), GEO%Volume(iElem))
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
      CALL FindNearestNeigh(TreeNode%iPartIndx_Node, nPart &
                                , iElem, GEO%Volume(iElem))
    END IF
  ELSE  IF (nPart.GT.1) THEN
    CALL FindNearestNeigh(TreeNode%iPartIndx_Node, nPart &
                              , iElem, GEO%Volume(iElem))
  END IF

  DEALLOCATE(TreeNode%iPartIndx_Node)
  DEALLOCATE(TreeNode)
END IF !nPart > 0

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
      SpecPartNum = 0
      DO iPart = 1, PartNumChildNode(iLoop)
        SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) = &
          SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) + 1
      END DO
      DSMC%MeanFreePath = CalcMeanFreePath(REAL(SpecPartNum),REAL(PartNumChildNode(iLoop)),NodeVolumeTemp(iLoop))
    END IF
    ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
    IF(PartNumChildNode(iLoop).GE.DSMC%PartNumOctreeNodeMin) THEN
      ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
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
        CALL FindNearestNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
      END IF
    ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
      CALL FindNearestNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
    END IF
  END DO

END SUBROUTINE AddOctreeNode

SUBROUTINE DSMC_init_octree()
!===================================================================================================================================
! Building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : ElemNodeVol
  USE MOD_Mesh_Vars,              ONLY : nElems
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
      CALL DSMC_CalcSubNodeVolumes(iElem, NodeDepth, ElemNodeVol(iElem)%Root)
    END DO
  END DO

END SUBROUTINE DSMC_init_octree

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


#ifdef DONTCOMPILETHIS
SUBROUTINE FindStatisticalNeigh(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
! Classic statistical pairing method for the use in the octree routines
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Relaxation,       ONLY : SetMeanVibQua
  USE MOD_DSMC_CollisionProb,    ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,           ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,        ONLY : DSMC_vmpf_prob
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, CollInf, CollisMode, PartStateIntEn, ChemReac, PairE_vMPF, BGGas, DSMC
  USE MOD_Particle_Vars,         ONLY : PartSpecies, nSpecies, PartState, usevMPF, PartMPF, WriteMacroVolumeValues
  USE MOD_TimeDisc_Vars,         ONLY : TEnd, time
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
  INTEGER                       :: nPair, iPair, iPart, cPart1, cPart2, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase
  REAL                          :: iRan
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

  nPart = PartNum
  nPair = INT(nPart/2)

  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
    ChemReac%nPairForRec = 0
  END IF

  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  Coll_pData%Ec=0

  IF (CollisMode.EQ.3) THEN
    ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
    DO iPart = 1, PartNum
      ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) = &
        ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_Node(iPart))) + PartStateIntEn(iPartIndx_Node(iPart),1)
    END DO
  END IF

  DO iPart = 1, PartNum
    CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
              CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + 1
  END DO

  IF (usevMPF) MPFFac = 1

  DO iPair = 1, nPair                               ! statistical pairing
    CALL RANDOM_NUMBER(iRan)
    cPart1 = 1 + INT(nPart * iRan)                       ! first pair particle
    Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(cPart1)
    iPartIndx_Node(cPart1) = iPartIndx_Node(nPart)
    nPart = nPart - 1
    CALL RANDOM_NUMBER(iRan)
    cPart2 = 1 + INT(nPart * iRan)                       ! second pair particle
    Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(cPart2)
    iPartIndx_Node(cPart2) = iPartIndx_Node(nPart)
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

  IF (CollisMode.EQ.3) THEN
    CALL SetMeanVibQua()
  END IF

  DO iPair = 1, nPair
    IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
      IF (usevMPF.AND.(BGGas%BGGasSpecies.EQ.0)) THEN            ! calculation of collision prob
        CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
      ELSE
        CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
      END IF
      CALL RANDOM_NUMBER(iRan)
      IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
        CALL DSMC_perform_collision(iPair,iElem, NodeVolume, PartNum)
      END IF
    END IF
  END DO

  IF(DSMC%CalcQualityFactors) THEN
    IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
      ! Determination of the maximum MCS/MFP for the cell
      IF (DSMC%CollSepCount.GT.0 .AND. DSMC%MeanFreePath.GT.0.0) THEN
        DSMC%MCSoverMFP = MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
      END IF
    END IF
  END IF

  DEALLOCATE(Coll_pData)

END SUBROUTINE FindStatisticalNeigh
#endif /*DONTCOMPILETHIS,NOOOW*/


END MODULE MOD_DSMC_ParticlePairing
