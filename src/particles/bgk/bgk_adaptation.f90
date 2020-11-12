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
!> Module containing routines for the recursive octree cell refinement algorithm for the BGK and FP-Flow particle methods.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGK_octree_adapt
  MODULE PROCEDURE BGK_octree_adapt
END INTERFACE

INTERFACE BGK_quadtree_adapt
  MODULE PROCEDURE BGK_quadtree_adapt
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGK_octree_adapt, BGK_quadtree_adapt
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_octree_adapt(iElem)
!===================================================================================================================================
!> Main octree routine for BGK and FP-Flow: Checks whether particle number (or density) is above the given limit and performs
!> a recursive octree algorithm to subdivide the cell until the limit is reached.
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, Time
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, ElemNodeVol, DSMC, RadialWeighting
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Vars           ,ONLY: PEM, PartState, PartPosRef,Species,WriteMacroVolumeValues, usevMPF, PartSpecies
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping
USE MOD_BGK_CollOperator        ,ONLY: BGK_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKMinPartPerCell,BGKMovingAverage,ElemNodeAveraging,BGKSplittingDens,BGKMovingAverageLength
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_FP_CollOperator         ,ONLY: FP_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKInitDone,BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor
USE MOD_BGK_Vars                ,ONLY: BGK_QualityFacSamp, BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_FPFlow_Vars             ,ONLY: FPInitDone, FP_PrandtlNumber, FP_QualityFacSamp
USE MOD_FPFlow_Vars             ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
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

INTEGER                       :: iPart, iLoop, nPart
REAL                          :: Dens, partWeight, totalWeight
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF

! Skip cell if number of particles is less than 2, create particle list (iPartIndx_Node) and sum-up bulk velocity
nPart = PEM%pNumber(iElem)
IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

NULLIFY(TreeNode)
ALLOCATE(TreeNode)
ALLOCATE(TreeNode%iPartIndx_Node(nPart))
TreeNode%iPartIndx_Node(1:nPart) = 0

totalWeight = 0.0
iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  TreeNode%iPartIndx_Node(iLoop) = iPart
  partWeight = GetParticleWeight(iPart)
  totalWeight = totalWeight + partWeight
  iPart = PEM%pNext(iPart)
END DO

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  Dens = totalWeight / GEO%Volume(iElem)
ELSE
  Dens = totalWeight * Species(1)%MacroParticleFactor / GEO%Volume(iElem)
END IF

! The octree refinement is performed if either the particle number or number density is above a user-given limit
IF(nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)) THEN
  ALLOCATE(TreeNode%MappedPartStates(1:3,1:nPart))
  TreeNode%PNum_Node = nPart
  IF (DoRefMapping) THEN
    DO iLoop = 1, nPart
      TreeNode%MappedPartStates(1:3,iLoop)=PartPosRef(1:3,TreeNode%iPartIndx_Node(iLoop))
    END DO
  ELSE ! position in reference space [-1,1] has to be computed
    DO iLoop = 1, nPart
      CALL GetPositionInRefElem(PartState(1:3,TreeNode%iPartIndx_Node(iLoop)),TreeNode%MappedPartStates(1:3,iLoop),iElem)
    END DO
  END IF ! DoRefMapping
  TreeNode%NodeDepth = 1
  TreeNode%MidPoint(1:3) = (/0.0,0.0,0.0/)
  ! Start of the recursive routine, which will descend further down the octree until the aforementioned criteria are fulfilled
  IF (BGKMovingAverage) THEN
    CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root, ElemNodeAveraging(iElem)%Root)
  ELSE
    CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
  END IF
  DEALLOCATE(TreeNode%MappedPartStates)
ELSE ! No octree refinement: Call of the respective collision operator
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPart &
                              , GEO%Volume(iElem))
#else
  IF (BGKMovingAverage) THEN
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, &
              GEO%Volume(iElem), &
             ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKMovingAverageLength), &
             CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
  ELSE
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart &
                            , GEO%Volume(iElem))
  END IF
#endif
END IF ! nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)

! Sampling of quality factors for BGK and FP-Flow methods
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    IF(BGKInitDone) THEN
      BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
      BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
      BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
      BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
      BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
      BGK_QualityFacSamp(6,iElem) = BGK_QualityFacSamp(6,iElem) + BGK_PrandtlNumber
      BGK_QualityFacSamp(7,iElem) = BGK_QualityFacSamp(7,iElem) + BGK_ExpectedPrandtlNumber
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
!> Adds an additional octree node/branch until either the particle number or number density is above a user-given limit
!> 1.) Sorting the particles into the subcells (octree child nodes)
!> 2.) Calculate the volumes of the subcells
!> 3.) Combines subcells, if the particle number within the subcell is below the limit
!> 4.) Check each child node if a further refinement is required
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: tTreeNode, tNodeVolume, ElemNodeVol
USE MOD_Particle_Vars         ,ONLY: PartState
USE MOD_BGK_CollOperator      ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_CalcSubNodeVolumes
USE MOD_BGK_Vars              ,ONLY: BGKMinPartPerCell,tNodeAverage, BGKMovingAverage, BGKMovingAverageLength
USE MOD_FP_CollOperator       ,ONLY: FP_CollisionOperator
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
REAL                         :: NodeVolumeTemp(8)
LOGICAL                      :: CombineChildNodes
!===================================================================================================================================
! 0. Determine if subcells with less particles than the given limit might occur (MARCEL FRAGEN)
IF (TreeNode%PNum_Node.LE.8.*BGKMinPartPerCell) THEN
  CombineChildNodes = .TRUE.
ELSE
  CombineChildNodes = .FALSE.
END IF
ALLOCATE(iPartIndx_ChildNode(8,TreeNode%PNum_Node))
ALLOCATE(MappedPart_ChildNode(1:3,TreeNode%PNum_Node,1:8))
PartNumChildNode(:) = 0
IF (ABS(TreeNode%MidPoint(1)) .EQ. 1.0) THEN
  CALL Abort(&
__STAMP__&
,'ERROR in Octree Pairing: Too many branches, machine precision reached')
END IF

! 1.) Sorting the received particles to the respective octree child node
!
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
!

DO iPart=1,TreeNode%PNum_Node
  iPartIndx = TreeNode%iPartIndx_Node(iPart)
  IF ((TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(2,iPart).GE.TreeNode%MidPoint(2)) &
      .AND.(TreeNode%MappedPartStates(3,iPart).LE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(1) = PartNumChildNode(1) + 1
    iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(1),1) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF((TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(2,iPart).GE.TreeNode%MidPoint(2))) THEN
    PartNumChildNode(2) = PartNumChildNode(2) + 1
    iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(2),2) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF((TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(3,iPart).GE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(3) = PartNumChildNode(3) + 1
    iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(3),3) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF (TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) THEN
    PartNumChildNode(4) = PartNumChildNode(4) + 1
    iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(4),4) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF((TreeNode%MappedPartStates(2,iPart).GE.TreeNode%MidPoint(2)) &
      .AND.(TreeNode%MappedPartStates(3,iPart).LE.TreeNode%MidPoint(3))) THEN
    PartNumChildNode(5) = PartNumChildNode(5) + 1
    iPartIndx_ChildNode(5,PartNumChildNode(5)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(5),5) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF (TreeNode%MappedPartStates(2,iPart).GE.TreeNode%MidPoint(2)) THEN
    PartNumChildNode(6) = PartNumChildNode(6) + 1
    iPartIndx_ChildNode(6,PartNumChildNode(6)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(6),6) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE IF (TreeNode%MappedPartStates(3,iPart).GE.TreeNode%MidPoint(3)) THEN
    PartNumChildNode(7) = PartNumChildNode(7) + 1
    iPartIndx_ChildNode(7,PartNumChildNode(7)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(7),7) = TreeNode%MappedPartStates(1:3,iPart)
  ELSE
    PartNumChildNode(8) = PartNumChildNode(8) + 1
    iPartIndx_ChildNode(8,PartNumChildNode(8)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(8),8) = TreeNode%MappedPartStates(1:3,iPart)
  END IF
END DO

! Check if any of the subcells has less particles than the limit, if so perform a recombination of cells (3.)
IF(ANY(PartNumChildNode.LT.BGKMinPartPerCell)) CombineChildNodes = .TRUE.

! 2.) Calculate the subcell volume (if necessary)
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

IF (BGKMovingAverage) THEN
  IF (.NOT.ASSOCIATED(Averaging%SubNode1)) THEN
    CALL BGK_AllocateAveragingNode(Averaging)
  END IF
END IF

! 3.) Combine subcells together if the particle number is less than the limit (BGKMinPartPerCell). Go through the first 7 subcells
!    and if the subcell is below the limit, add the particles and the volume to the next subcell and delete them from the original.
!    For the last subcell: if it has less than the limit, find a cell, which still has particles and add them to it.
IF(CombineChildNodes) THEN
  DO iLoop = 1, 7
    IF (PartNumChildNode(iLoop).LT.BGKMinPartPerCell) THEN
      DO iPart=1, PartNumChildNode(iLoop)
        iPartIndx_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart) = iPartIndx_ChildNode(iLoop,iPart)
        MappedPart_ChildNode(1:3,PartNumChildNode(iLoop+1)+iPart,iLoop+1) = MappedPart_ChildNode(1:3,iPart,iLoop)
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
      MappedPart_ChildNode(1:3,PartNumChildNode(iLoop2)+iPart,iLoop2) = MappedPart_ChildNode(1:3,iPart,8)
    END DO
    PartNumChildNode(iLoop2) = PartNumChildNode(iLoop2) + PartNumChildNode(8)
    PartNumChildNode(8) = 0
    NodeVolumeTemp(iLoop2) = NodeVolumeTemp(iLoop2) + NodeVolumeTemp(8)
    NodeVolumeTemp(8) = 0.0
  END IF
END IF

! 4.) Check each child node if a further refinement is required. If no further refinement is necessary or if cells were combined
!    -> Perform the collision operator
DO iLoop = 1, 8
  IF ((PartNumChildNode(iLoop).GE.(2.*BGKMinPartPerCell)).AND.(.NOT.CombineChildNodes)) THEN
    NULLIFY(TreeNode%ChildNode)
    ALLOCATE(TreeNode%ChildNode)
    ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
    ALLOCATE(TreeNode%ChildNode%MappedPartStates(1:3,PartNumChildNode(iLoop)))
    TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
    TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
    TreeNode%ChildNode%MappedPartStates(1:3,1:PartNumChildNode(iLoop))= &
                   MappedPart_ChildNode(1:3,1:PartNumChildNode(iLoop),iLoop)
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
    IF (BGKMovingAverage) THEN
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
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop))
#else
    IF (BGKMovingAverage) THEN
      SELECT CASE(iLoop)
        CASE(1)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode1%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode1%CorrectStep)
        CASE(2)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode2%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode2%CorrectStep)
        CASE(3)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode3%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode3%CorrectStep)
        CASE(4)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode4%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode4%CorrectStep)
        CASE(5)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode5%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode5%CorrectStep)
        CASE(6)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode6%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode6%CorrectStep)
        CASE(7)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode7%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode7%CorrectStep)
        CASE(8)
          CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
                Averaging%SubNode8%AverageValues(1:5,1:BGKMovingAverageLength), &
                CorrectStep = Averaging%SubNode8%CorrectStep)
      END SELECT
    ELSE
      CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop))
    END IF
#endif
  END IF
END DO

END SUBROUTINE AddBGKOctreeNode

SUBROUTINE BGK_AllocateAveragingNode(Averaging)
!===================================================================================================================================
!> Allocation of the arrays and iteration counter required for the sampling of the moving average in the octree subnodes
!===================================================================================================================================
! MODULES
USE MOD_BGK_Vars,               ONLY :tNodeAverage, BGKMovingAverageLength
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

ALLOCATE(Averaging%SubNode1%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode2%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode3%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode4%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode5%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode6%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode7%AverageValues(5,BGKMovingAverageLength))
ALLOCATE(Averaging%SubNode8%AverageValues(5,BGKMovingAverageLength))

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

SUBROUTINE BGK_quadtree_adapt(iElem)
!===================================================================================================================================
!> Main quadtree routine for BGK and FP-Flow: Checks whether particle number (or density) is above the given limit and performs
!> a recursive quadtree algorithm to subdivide the cell until the limit is reached.
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, Time
USE MOD_DSMC_ParticlePairing    ,ONLY: GeoCoordToMap2D
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, ElemNodeVol, DSMC, RadialWeighting
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Vars           ,ONLY: PEM, PartState, Species,WriteMacroVolumeValues, usevMPF
USE MOD_BGK_CollOperator        ,ONLY: BGK_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKMinPartPerCell,BGKSplittingDens!,BGKMovingAverage,ElemNodeAveraging,BGKMovingAverageLength
USE MOD_FP_CollOperator         ,ONLY: FP_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKInitDone,BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor
USE MOD_BGK_Vars                ,ONLY: BGK_QualityFacSamp, BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_FPFlow_Vars             ,ONLY: FPInitDone, FP_PrandtlNumber, FP_QualityFacSamp
USE MOD_FPFlow_Vars             ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
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
INTEGER                       :: iPart, iLoop, nPart
REAL                          :: Dens, partWeight, totalWeight
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF

! Skip cell if number of particles is less than 2, create particle list (iPartIndx_Node) and sum-up bulk velocity
nPart = PEM%pNumber(iElem)
IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

NULLIFY(TreeNode)
ALLOCATE(TreeNode)
ALLOCATE(TreeNode%iPartIndx_Node(nPart))
TreeNode%iPartIndx_Node(1:nPart) = 0

totalWeight = 0.0
iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  TreeNode%iPartIndx_Node(iLoop) = iPart
  partWeight = GetParticleWeight(iPart)
  totalWeight = totalWeight + partWeight
  iPart = PEM%pNext(iPart)
END DO

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  Dens = totalWeight / GEO%Volume(iElem)
ELSE
  Dens = totalWeight * Species(1)%MacroParticleFactor / GEO%Volume(iElem)
END IF

! The quadtree refinement is performed if either the particle number or number density is above a user-given limit
IF(nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)) THEN
  ALLOCATE(TreeNode%MappedPartStates(1:2,1:nPart))
  TreeNode%PNum_Node = nPart
  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    CALL GeoCoordToMap2D(PartState(1:2,iPart), TreeNode%MappedPartStates(1:2,iLoop), iElem)
    iPart = PEM%pNext(iPart)
  END DO
  TreeNode%NodeDepth = 1
  TreeNode%MidPoint(1:3) = (/0.0,0.0,0.0/)
  ! Start of the recursive routine, which will descend further down the quadtree until the aforementioned criteria are fulfilled
  ! IF (BGKMovingAverage) THEN
  !   CALL AddBGKQuadtreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root, ElemNodeAveraging(iElem)%Root)
  ! ELSE
    CALL AddBGKQuadtreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
  ! END IF
  DEALLOCATE(TreeNode%MappedPartStates)
ELSE ! No quadtree refinement: Call of the respective collision operator
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPart &
                              , GEO%Volume(iElem))
#else
  ! IF (BGKMovingAverage) THEN
  !   CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, &
  !             GEO%Volume(iElem), &
  !            ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKMovingAverageLength), &
  !            CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
  ! ELSE
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart &
                            , GEO%Volume(iElem))
  ! END IF
#endif
END IF ! nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)

! Sampling of quality factors for BGK and FP-Flow methods
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    IF(BGKInitDone) THEN
      BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
      BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
      BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
      BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
      BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
      BGK_QualityFacSamp(6,iElem) = BGK_QualityFacSamp(6,iElem) + BGK_PrandtlNumber
      BGK_QualityFacSamp(7,iElem) = BGK_QualityFacSamp(7,iElem) + BGK_ExpectedPrandtlNumber
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

END SUBROUTINE BGK_quadtree_adapt


! RECURSIVE SUBROUTINE AddBGKQuadtreeNode(TreeNode, iElem, NodeVol, Averaging)
RECURSIVE SUBROUTINE AddBGKQuadtreeNode(TreeNode, iElem, NodeVol)
!===================================================================================================================================
!> Adds an additional quadtree node/branch until either the particle number or number density is above a user-given limit
!> 1.) Sorting the particles into the subcells (quadtree child nodes)
!> 2.) Calculate the volumes of the subcells
!> 3.) Combines subcells, if the particle number within the subcell is below the limit
!> 4.) Check each child node if a further refinement is required
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: tTreeNode, tNodeVolume, ElemNodeVol
USE MOD_Particle_Vars         ,ONLY: PartState
USE MOD_BGK_CollOperator      ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_CalcSubNodeVolumes2D
USE MOD_BGK_Vars              ,ONLY: BGKMinPartPerCell    !,tNodeAverage, BGKMovingAverage, BGKMovingAverageLength
USE MOD_FP_CollOperator       ,ONLY: FP_CollisionOperator
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                                 :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER                 :: TreeNode
TYPE(tNodeVolume),INTENT(IN), POINTER               :: NodeVol
! TYPE(tNodeAverage),INTENT(INOUT), POINTER, OPTIONAL :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, iLoop, iPartIndx, localDepth, iLoop2
INTEGER, ALLOCATABLE         :: iPartIndx_ChildNode(:,:)
REAL, ALLOCATABLE            :: MappedPart_ChildNode(:,:,:)
INTEGER                      :: PartNumChildNode(1:4)
REAL                         :: NodeVolumeTemp(1:4)
LOGICAL                      :: CombineChildNodes
!===================================================================================================================================
! 0. Determine if subcells with less particles than the given limit might occur (MARCEL FRAGEN)
IF (TreeNode%PNum_Node.LE.4.*BGKMinPartPerCell) THEN
  CombineChildNodes = .TRUE.
ELSE
  CombineChildNodes = .FALSE.
END IF
ALLOCATE(iPartIndx_ChildNode(4,TreeNode%PNum_Node))
ALLOCATE(MappedPart_ChildNode(1:3,TreeNode%PNum_Node,1:4))
PartNumChildNode(:) = 0
IF (ABS(TreeNode%MidPoint(1)) .EQ. 1.0) THEN
  CALL Abort(&
__STAMP__&
,'ERROR in BGK/FP Quadtree Refinement: Too many branches, machine precision reached')
END IF

! 1.) Sorting the received particles to the respective QuadTree child node
!
!      Numbering of the 4 ChildNodes of the QuadTree
!      _________
!     |    |    |   |y
!     | 3  | 2  |   |
!     |----|----|   |
!     | 4  | 1  |   |______x
!     |____|____|

DO iPart=1,TreeNode%PNum_Node
  iPartIndx = TreeNode%iPartIndx_Node(iPart)
  IF ((TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(2,iPart).LE.TreeNode%MidPoint(2))) THEN
    PartNumChildNode(1) = PartNumChildNode(1) + 1
    iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(1),1) = TreeNode%MappedPartStates(1:2,iPart)
  ELSE IF(TreeNode%MappedPartStates(1,iPart).GE.TreeNode%MidPoint(1)) THEN
    PartNumChildNode(2) = PartNumChildNode(2) + 1
    iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(2),2) = TreeNode%MappedPartStates(1:2,iPart)
  ELSE IF(TreeNode%MappedPartStates(2,iPart).GE.TreeNode%MidPoint(2)) THEN
    PartNumChildNode(3) = PartNumChildNode(3) + 1
    iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(3),3) = TreeNode%MappedPartStates(1:2,iPart)
  ELSE
    PartNumChildNode(4) = PartNumChildNode(4) + 1
    iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(4),4) = TreeNode%MappedPartStates(1:2,iPart)
  END IF
END DO

! Check if any of the subcells has less particles than the limit, if so perform a recombination of cells (3.)
IF(ANY(PartNumChildNode.LT.BGKMinPartPerCell)) CombineChildNodes = .TRUE.

! 2.) Calculate the subcell volume (if necessary)
IF((.NOT.ASSOCIATED(NodeVol)).OR.(.NOT.ASSOCIATED(NodeVol%SubNode1))) THEN
  localDepth = TreeNode%NodeDepth
  CALL DSMC_CalcSubNodeVolumes2D(iElem, localDepth, ElemNodeVol(iElem)%Root)
END IF

NodeVolumeTemp(1) = NodeVol%SubNode1%Volume
NodeVolumeTemp(2) = NodeVol%SubNode2%Volume
NodeVolumeTemp(3) = NodeVol%SubNode3%Volume
NodeVolumeTemp(4) = NodeVol%SubNode4%Volume

!---- NOT IMPLEMENTED YET ---
! IF (BGKMovingAverage) THEN
!   IF (.NOT.ASSOCIATED(Averaging%SubNode1)) THEN
!     CALL BGK_AllocateAveragingNode(Averaging)
!   END IF
! END IF

! 3.) Combine subcells together if the particle number is less than the limit (BGKMinPartPerCell). Go through the first 7 subcells
!    and if the subcell is below the limit, add the particles and the volume to the next subcell and delete them from the original.
!    For the last subcell: if it has less than the limit, find a cell, which still has particles and add them to it.
IF(CombineChildNodes) THEN
  DO iLoop = 1, 3
    IF (PartNumChildNode(iLoop).LT.BGKMinPartPerCell) THEN
      DO iPart=1, PartNumChildNode(iLoop)
        iPartIndx_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart) = iPartIndx_ChildNode(iLoop,iPart)
        MappedPart_ChildNode(1:3,PartNumChildNode(iLoop+1)+iPart,iLoop+1) = MappedPart_ChildNode(1:3,iPart,iLoop)
      END DO
      PartNumChildNode(iLoop+1) = PartNumChildNode(iLoop+1) + PartNumChildNode(iLoop)
      PartNumChildNode(iLoop) = 0
      NodeVolumeTemp(iLoop+1) = NodeVolumeTemp(iLoop+1) + NodeVolumeTemp(iLoop)
      NodeVolumeTemp(iLoop) = 0.0
    END IF
  END DO
  IF (PartNumChildNode(4).LT.BGKMinPartPerCell) THEN
    DO iLoop = 1, 3
     iLoop2 = iLoop
     IF (PartNumChildNode(iLoop).GT.0) EXIT
    END DO
    DO iPart=1, PartNumChildNode(4)
      iPartIndx_ChildNode(iLoop2,PartNumChildNode(iLoop2)+iPart) = iPartIndx_ChildNode(4,iPart)
      MappedPart_ChildNode(1:3,PartNumChildNode(iLoop2)+iPart,iLoop2) = MappedPart_ChildNode(1:3,iPart,4)
    END DO
    PartNumChildNode(iLoop2) = PartNumChildNode(iLoop2) + PartNumChildNode(4)
    PartNumChildNode(4) = 0
    NodeVolumeTemp(iLoop2) = NodeVolumeTemp(iLoop2) + NodeVolumeTemp(4)
    NodeVolumeTemp(4) = 0.0
  END IF
END IF

! 4.) Check each child node if a further refinement is required. If no further refinement is necessary or if cells were combined
!    -> Perform the collision operator
DO iLoop = 1, 4
  IF ((PartNumChildNode(iLoop).GE.(2.*BGKMinPartPerCell)).AND.(.NOT.CombineChildNodes)) THEN
    NULLIFY(TreeNode%ChildNode)
    ALLOCATE(TreeNode%ChildNode)
    ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
    ALLOCATE(TreeNode%ChildNode%MappedPartStates(1:3,PartNumChildNode(iLoop)))
    TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
    TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
    TreeNode%ChildNode%MappedPartStates(1:3,1:PartNumChildNode(iLoop))= &
          MappedPart_ChildNode(1:3,1:PartNumChildNode(iLoop),iLoop)
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
    ! IF (BGKMovingAverage) THEN
    !   IF (iLoop.EQ.1) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1, Averaging%SubNode1)
    !   IF (iLoop.EQ.2) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2, Averaging%SubNode2)
    !   IF (iLoop.EQ.3) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3, Averaging%SubNode3)
    !   IF (iLoop.EQ.4) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4, Averaging%SubNode4)
    ! ELSE
    IF (iLoop.EQ.1) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode1)
    IF (iLoop.EQ.2) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode2)
    IF (iLoop.EQ.3) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode3)
    IF (iLoop.EQ.4) CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode4)
    ! END IF
    DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
    DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
    DEALLOCATE(TreeNode%ChildNode)
  ELSE IF (PartNumChildNode(iLoop).GE.2) THEN
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop))
#else
    ! IF (BGKMovingAverage) THEN
    !   SELECT CASE(iLoop)
    !     CASE(1)
    !       CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
    !             PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
    !             Averaging%SubNode1%AverageValues(1:5,1:BGKMovingAverageLength), &
    !             CorrectStep = Averaging%SubNode1%CorrectStep)
    !     CASE(2)
    !       CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
    !             PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
    !             Averaging%SubNode2%AverageValues(1:5,1:BGKMovingAverageLength), &
    !             CorrectStep = Averaging%SubNode2%CorrectStep)
    !     CASE(3)
    !       CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
    !             PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
    !             Averaging%SubNode3%AverageValues(1:5,1:BGKMovingAverageLength), &
    !             CorrectStep = Averaging%SubNode3%CorrectStep)
    !     CASE(4)
    !       CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
    !             PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), &
    !             Averaging%SubNode4%AverageValues(1:5,1:BGKMovingAverageLength), &
    !             CorrectStep = Averaging%SubNode4%CorrectStep)
    !   END SELECT
    ! ELSE
      CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
              PartNumChildNode(iLoop), NodeVolumeTemp(iLoop))
    ! END IF
#endif
  END IF
END DO

END SUBROUTINE AddBGKQuadtreeNode

END MODULE MOD_BGK_Adaptation
