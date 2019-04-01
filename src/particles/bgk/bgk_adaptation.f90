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

MODULE MOD_BGK_Adaptation
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGK_octree_adapt
  MODULE PROCEDURE BGK_octree_adapt
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGK_octree_adapt
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_octree_adapt(iElem)
!===================================================================================================================================
!> Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node 
!> or start nearest neighbour search
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars,          ONLY: TEnd, Time
USE MOD_DSMC_Vars              ,ONLY: tTreeNode, ElemNodeVol, DSMC
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Vars          ,ONLY: PEM, PartState, PartPosRef,Species,WriteMacroVolumeValues
USE MOD_Particle_Tracking_vars ,ONLY: DoRefMapping
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_BGK_Vars               ,ONLY: BGKMinPartPerCell, BGKDoAveraging, ElemNodeAveraging, BGKAveragingLength,BGKSplittingDens
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperatorOctree
USE MOD_BGK_Vars               ,ONLY: BGKInitDone,BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor
USE MOD_BGK_Vars               ,ONLY: BGK_QualityFacSamp, BGK_MaxRotRelaxFactor
USE MOD_FPFlow_Vars            ,ONLY: FPInitDone, FP_PrandtlNumber, FP_QualityFacSamp
USE MOD_FPFlow_Vars            ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, nPart
REAL                          :: vBulk(3), Dens
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF

nPart = PEM%pNumber(iElem)
IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

NULLIFY(TreeNode)

ALLOCATE(TreeNode)
ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing
TreeNode%iPartIndx_Node(1:nPart) = 0

vBulk(1:3) = 0.0
iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
DO iLoop = 1, nPart
  TreeNode%iPartIndx_Node(iLoop) = iPart
  vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)
  iPart = PEM%pNext(iPart)
END DO
vBulk = vBulk / nPart

Dens = nPart * Species(1)%MacroParticleFactor / GEO%Volume(iElem)
! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
IF(nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)) THEN
  ! Additional check afterwards if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
  ! the side length of a cube (approximation) with same volume as the actual cell -> octree
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
      !CALL GeoCoordToMap(PartState(iPart,1:3), TreeNode%MappedPartStates(iLoop,1:3), iElem)
      iPart = PEM%pNext(iPart)
    END DO
  END IF ! DoRefMapping
  TreeNode%NodeDepth = 1
  TreeNode%MidPoint(1:3) = (/0.0,0.0,0.0/)
  IF (BGKDoAveraging) THEN
    CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root, ElemNodeAveraging(iElem)%Root)
  ELSE
    CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
  END IF
  DEALLOCATE(TreeNode%MappedPartStates)
ELSE
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperatorOctree(TreeNode%iPartIndx_Node, nPart &
                              , GEO%Volume(iElem), vBulk)
#else
  IF (BGKDoAveraging) THEN
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, &
              GEO%Volume(iElem), vBulk, &
             ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKAveragingLength), &
             CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
  ELSE
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart &
                            , GEO%Volume(iElem), vBulk)
  END IF
#endif
END IF

IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    IF(BGKInitDone) THEN
      BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
      BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
      BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
      BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
      BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
    END IF
    IF(FPInitDone) THEN
      FP_QualityFacSamp(1,iElem) = FP_QualityFacSamp(1,iElem) + FP_MeanRelaxFactor
      FP_QualityFacSamp(2,iElem) = FP_QualityFacSamp(2,iElem) + REAL(FP_MeanRelaxFactorCounter)
      FP_QualityFacSamp(3,iElem) = FP_QualityFacSamp(3,iElem) + FP_MaxRelaxFactor
      FP_QualityFacSamp(4,iElem) = FP_QualityFacSamp(4,iElem) + 1.
      FP_QualityFacSamp(5,iElem) = FP_QualityFacSamp(5,iElem) + FP_MaxRotRelaxFactor
      FP_QualityFacSamp(6,iElem) = FP_QualityFacSamp(6,iElem) + FP_PrandtlNumber
    END IF
  END IF
END IF

DEALLOCATE(TreeNode%iPartIndx_Node)
DEALLOCATE(TreeNode)

END SUBROUTINE BGK_octree_adapt


RECURSIVE SUBROUTINE AddBGKOctreeNode(TreeNode, iElem, NodeVol, Averaging)
!===================================================================================================================================
!> Adds additional octree node/branch (fancy)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: tTreeNode, tNodeVolume, ElemNodeVol
USE MOD_Particle_Vars         ,ONLY: PartState
USE MOD_BGK_CollOperator      ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_CalcSubNodeVolumes
USE MOD_BGK_Vars              ,ONLY: BGKMinPartPerCell,tNodeAverage, BGKDoAveraging
USE MOD_BGK_Vars              ,ONLY: BGKAveragingLength
USE MOD_FP_CollOperator       ,ONLY: FP_CollisionOperatorOctree
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                                 :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER                 :: TreeNode
TYPE(tNodeVolume),INTENT(IN), POINTER               :: NodeVol
TYPE(tNodeAverage),INTENT(INOUT), POINTER, OPTIONAL :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, iLoop, iPartIndx, localDepth, iLoop2
INTEGER, ALLOCATABLE         :: iPartIndx_ChildNode(:,:)
REAL, ALLOCATABLE            :: MappedPart_ChildNode(:,:,:)
INTEGER                      :: PartNumChildNode(8)
REAL                         :: NodeVolumeTemp(8), vBulk(3,8)
LOGICAL                      :: ForceBGK
!===================================================================================================================================
IF (TreeNode%PNum_Node.LE.8.*BGKMinPartPerCell) THEN
  ForceBGK = .TRUE.
ELSE
  ForceBGK = .FALSE.
END IF
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
vBulk = 0.0
DO iPart=1,TreeNode%PNum_Node
  iPartIndx = TreeNode%iPartIndx_Node(iPart)
  IF ((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) &
      .AND.(TreeNode%MappedPartStates(iPart,3).LE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(1) = PartNumChildNode(1) + 1
    iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    vBulk(1:3,1) = vBulk(1:3,1) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(1,PartNumChildNode(1),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2))) THEN
    PartNumChildNode(2) = PartNumChildNode(2) + 1
    iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    vBulk(1:3,2) = vBulk(1:3,2) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(2,PartNumChildNode(2),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF((TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(iPart,3).GE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(3) = PartNumChildNode(3) + 1
    iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
    vBulk(1:3,3) = vBulk(1:3,3) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(3,PartNumChildNode(3),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF (TreeNode%MappedPartStates(iPart,1).GE.TreeNode%MidPoint(1)) THEN
    PartNumChildNode(4) = PartNumChildNode(4) + 1
    iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
    vBulk(1:3,4) = vBulk(1:3,4) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(4,PartNumChildNode(4),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF((TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) &
      .AND.(TreeNode%MappedPartStates(iPart,3).LE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(5) = PartNumChildNode(5) + 1
    iPartIndx_ChildNode(5,PartNumChildNode(5)) = iPartIndx
    vBulk(1:3,5) = vBulk(1:3,5) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(5,PartNumChildNode(5),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF (TreeNode%MappedPartStates(iPart,2).GE.TreeNode%MidPoint(2)) THEN
    PartNumChildNode(6) = PartNumChildNode(6) + 1
    iPartIndx_ChildNode(6,PartNumChildNode(6)) = iPartIndx
    vBulk(1:3,6) = vBulk(1:3,6) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(6,PartNumChildNode(6),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE IF (TreeNode%MappedPartStates(iPart,3).GE.TreeNode%MidPoint(3)) THEN
    PartNumChildNode(7) = PartNumChildNode(7) + 1
    iPartIndx_ChildNode(7,PartNumChildNode(7)) = iPartIndx
    vBulk(1:3,7) = vBulk(1:3,7) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(7,PartNumChildNode(7),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  ELSE
    PartNumChildNode(8) = PartNumChildNode(8) + 1
    iPartIndx_ChildNode(8,PartNumChildNode(8)) = iPartIndx
    vBulk(1:3,8) = vBulk(1:3,8) + PartState(iPartIndx,4:6)
    MappedPart_ChildNode(8,PartNumChildNode(8),1:3) = TreeNode%MappedPartStates(iPart,1:3)
  END IF
END DO

IF(ANY(PartNumChildNode.LT.BGKMinPartPerCell)) ForceBGK = .TRUE.

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

IF (BGKDoAveraging) THEN
  IF (.NOT.ASSOCIATED(Averaging%SubNode1)) THEN
    CALL BGK_AllocateAveragingNode(Averaging)
  END IF
END IF

IF(ForceBGK) THEN
  DO iLoop = 1, 7
    IF (PartNumChildNode(iLoop).LT.BGKMinPartPerCell) THEN
      DO iPart=1, PartNumChildNode(iLoop)
        iPartIndx_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart) = iPartIndx_ChildNode(iLoop,iPart)
        MappedPart_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart,1:3) = MappedPart_ChildNode(iLoop,iPart,1:3)
        vBulk(1:3,iLoop+1) = vBulk(1:3,iLoop+1) + PartState(iPartIndx_ChildNode(iLoop,iPart),4:6)
      END DO
      PartNumChildNode(iLoop+1) = PartNumChildNode(iLoop+1) + PartNumChildNode(iLoop)
      PartNumChildNode(iLoop) = 0
      NodeVolumeTemp(iLoop+1) = NodeVolumeTemp(iLoop+1) + NodeVolumeTemp(iLoop)
      NodeVolumeTemp(iLoop) = 0.0
    END IF
  END DO
  IF (PartNumChildNode(8).LT.BGKMinPartPerCell) THEN
    DO iLoop = 1, 7
     iLoop2 = iLoop
     IF (PartNumChildNode(iLoop).GT.0) EXIT
    END DO
    DO iPart=1, PartNumChildNode(8)
      iPartIndx_ChildNode(iLoop2,PartNumChildNode(iLoop2)+iPart) = iPartIndx_ChildNode(8,iPart)
      MappedPart_ChildNode(iLoop2,PartNumChildNode(iLoop2)+iPart,1:3) = MappedPart_ChildNode(8,iPart,1:3)
      vBulk(1:3,iLoop2) = vBulk(1:3,iLoop2) + PartState(iPartIndx_ChildNode(8,iPart),4:6)
    END DO
    PartNumChildNode(iLoop2) = PartNumChildNode(iLoop2) + PartNumChildNode(8)
    PartNumChildNode(8) = 0
    NodeVolumeTemp(iLoop2) = NodeVolumeTemp(iLoop2) + NodeVolumeTemp(8)
    NodeVolumeTemp(8) = 0.0
  END IF
END IF


DO iLoop = 1, 8
  ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
  IF ((PartNumChildNode(iLoop).GE.(2.*BGKMinPartPerCell)).AND.(.NOT.ForceBGK)) THEN
    ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
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
                                     + TreeNode%ChildNode%MidPoint(1:3)*2.0/(2.0**(TreeNode%NodeDepth+1.0))
    TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
    ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further octree division)
    IF (BGKDoAveraging) THEN
      IF (iLoop.EQ.1) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1, Averaging%SubNode1)
      IF (iLoop.EQ.2) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2, Averaging%SubNode2)
      IF (iLoop.EQ.3) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3, Averaging%SubNode3)
      IF (iLoop.EQ.4) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4, Averaging%SubNode4)
      IF (iLoop.EQ.5) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode5, Averaging%SubNode5)
      IF (iLoop.EQ.6) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode6, Averaging%SubNode6)
      IF (iLoop.EQ.7) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode7, Averaging%SubNode7)
      IF (iLoop.EQ.8) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode8, Averaging%SubNode8)
    ELSE
      IF (iLoop.EQ.1) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1)
      IF (iLoop.EQ.2) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2)
      IF (iLoop.EQ.3) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3)
      IF (iLoop.EQ.4) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4)
      IF (iLoop.EQ.5) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode5)
      IF (iLoop.EQ.6) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode6)
      IF (iLoop.EQ.7) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode7)
      IF (iLoop.EQ.8) CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode8)
    END IF
    DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
    DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
    DEALLOCATE(TreeNode%ChildNode)
  ELSE IF (PartNumChildNode(iLoop).GE.2) THEN
    vBulk(1:3,iLoop) = vBulk(1:3,iLoop) / PartNumChildNode(iLoop)
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperatorOctree(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop))
#else
    IF (BGKDoAveraging) THEN
      SELECT CASE(iLoop)
        CASE(1)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode1%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode1%CorrectStep)
        CASE(2)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode2%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode2%CorrectStep)
        CASE(3)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode3%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode3%CorrectStep)
        CASE(4)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode4%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode4%CorrectStep)
        CASE(5)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode5%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode5%CorrectStep)
        CASE(6)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode6%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode6%CorrectStep)
        CASE(7)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode7%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode7%CorrectStep)
        CASE(8)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop), &
                Averaging%SubNode8%AverageValues(1:5,1:BGKAveragingLength), &
                CorrectStep = Averaging%SubNode8%CorrectStep)
      END SELECT
    ELSE
      CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), vBulk(1:3,iLoop))
    END IF
#endif
  END IF
END DO

END SUBROUTINE AddBGKOctreeNode


integer function lcm(a,b)
integer:: a,b
    lcm = a*b / gcd(a,b)
end function lcm

integer function gcd(a,b)
integer :: a,b,t
    do while (b/=0)
        t = b
        b = mod(a,b)
        a = t
    end do
    gcd = abs(a)
end function gcd


SUBROUTINE BGK_AllocateAveragingNode(Averaging)
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES
USE MOD_BGK_Vars,               ONLY :tNodeAverage, BGKAveragingLength
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodeAverage),INTENT(INOUT), POINTER, OPTIONAL :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(Averaging%SubNode1)
ALLOCATE(Averaging%SubNode2)
ALLOCATE(Averaging%SubNode3)
ALLOCATE(Averaging%SubNode4)
ALLOCATE(Averaging%SubNode5)
ALLOCATE(Averaging%SubNode6)
ALLOCATE(Averaging%SubNode7)
ALLOCATE(Averaging%SubNode8)

ALLOCATE(Averaging%SubNode1%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode2%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode3%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode4%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode5%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode6%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode7%AverageValues(5,BGKAveragingLength))
ALLOCATE(Averaging%SubNode8%AverageValues(5,BGKAveragingLength))

Averaging%SubNode1%AverageValues = 0.0
Averaging%SubNode2%AverageValues = 0.0
Averaging%SubNode3%AverageValues = 0.0
Averaging%SubNode4%AverageValues = 0.0
Averaging%SubNode5%AverageValues = 0.0
Averaging%SubNode6%AverageValues = 0.0
Averaging%SubNode7%AverageValues = 0.0
Averaging%SubNode8%AverageValues = 0.0

Averaging%SubNode1%CorrectStep = 0
Averaging%SubNode2%CorrectStep = 0
Averaging%SubNode3%CorrectStep = 0
Averaging%SubNode4%CorrectStep = 0
Averaging%SubNode5%CorrectStep = 0
Averaging%SubNode6%CorrectStep = 0
Averaging%SubNode7%CorrectStep = 0
Averaging%SubNode8%CorrectStep = 0

END SUBROUTINE BGK_AllocateAveragingNode

END MODULE MOD_BGK_Adaptation
