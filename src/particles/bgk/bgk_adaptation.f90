!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_BGK_Adaptation
!===================================================================================================================================
!> Module containing routines for the recursive octree cell refinement algorithm for the BGK and FP-Flow particle methods.
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
USE MOD_Particle_Vars           ,ONLY: PEM, PartPosRef,Species,WriteMacroVolumeValues, usevMPF, LastPartPos, VirtMergedCells
USE MOD_Particle_Vars           ,ONLY: DoVirtualCellMerge
#if PP_TimeDiscMethod==300
USE MOD_Particle_Vars           ,ONLY: PartState
#endif /*PP_TimeDiscMethod==300*/
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_BGK_CollOperator        ,ONLY: BGK_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKMinPartPerCell,BGKSplittingDens
USE MOD_BGK_Vars                ,ONLY: BGKMovingAverage,ElemNodeAveraging
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_FP_CollOperator         ,ONLY: FP_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKInitDone,BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor
USE MOD_BGK_Vars                ,ONLY: BGK_QualityFacSamp, BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_BGK_Vars                ,ONLY: BGK_Viscosity, BGK_ThermalConductivity
USE MOD_FPFlow_Vars             ,ONLY: FPInitDone, FP_PrandtlNumber, FP_QualityFacSamp
USE MOD_FPFlow_Vars             ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars               ,ONLY: offsetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iLoop,nPart,CNElemID,GlobalElemID, nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem
REAL                          :: Dens,partWeight,totalWeight
LOGICAL                       :: DoMergedCell
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.; BGK_Viscosity=0.; BGK_ThermalConductivity=0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF
DoMergedCell = .FALSE.
! Skip cell if number of particles is less than 2, create particle list (iPartIndx_Node) and sum-up bulk velocity
nPart = PEM%pNumber(iElem)
IF (DoVirtualCellMerge) THEN
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    nPartMerged = nPart
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
    END DO
    IF (nPartMerged.LE.1) RETURN
    NULLIFY(TreeNode)
    ALLOCATE(TreeNode)
    ALLOCATE(TreeNode%iPartIndx_Node(nPartMerged))
    iPart = PEM%pStart(iElem)
    iLoopLoc = 0
    DO iLoop = 1, nPart
      iLoopLoc = iLoopLoc + 1
      TreeNode%iPartIndx_Node(iLoopLoc) = iPart
      iPart = PEM%pNext(iPart)
    END DO
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      locElem = VirtMergedCells(iElem)%MergedCellID(iMergeElem)
      nPartLoc = PEM%pNumber(locElem)
      iPart = PEM%pStart(locElem)
      DO iLoop = 1, nPartLoc
        iLoopLoc = iLoopLoc + 1
        TreeNode%iPartIndx_Node(iLoopLoc) = iPart
        iPart = PEM%pNext(iPart)
      END DO
    END DO
    DoMergedCell = .TRUE.
  END IF  
ELSE IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

IF (DoMergedCell) THEN
#if (PP_TimeDiscMethod==300)  
  CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged, VirtMergedCells(iELem)%MergedVolume)
#else
  IF (BGKMovingAverage) THEN
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged,VirtMergedCells(iELem)%MergedVolume, &
            ElemNodeAveraging(iElem)%Root%AverageValues(:))
  ELSE
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged, VirtMergedCells(iELem)%MergedVolume)
  END IF
#endif
ELSE
  GlobalElemID = iElem+offSetElem
  CNElemID     = GetCNElemID(GlobalElemID)

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
    Dens = totalWeight / ElemVolume_Shared(CNElemID)
  ELSE
    Dens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
  END IF

! The octree refinement is performed if either the particle number or number density is above a user-given limit
  IF(nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)) THEN
    ALLOCATE(TreeNode%MappedPartStates(1:3,1:nPart))
    TreeNode%PNum_Node = nPart
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      DO iLoop = 1, nPart
        TreeNode%MappedPartStates(1:3,iLoop)=PartPosRef(1:3,TreeNode%iPartIndx_Node(iLoop))
      END DO
    ELSE ! position in reference space [-1,1] has to be computed
      DO iLoop = 1, nPart
#if PP_TimeDiscMethod==300
        CALL GetPositionInRefElem(PartState(1:3,TreeNode%iPartIndx_Node(iLoop)),TreeNode%MappedPartStates(1:3,iLoop),GlobalElemID)
#else
      ! Attention: LastPartPos is the reference position here
       TreeNode%MappedPartStates(1:3,iLoop)=LastPartPos(1:3,TreeNode%iPartIndx_Node(iLoop))
#endif /*PP_TimeDiscMethod==300*/
      END DO
    END IF ! TrackingMethod.EQ.REFMAPPING
    TreeNode%NodeDepth = 1
    ElemNodeVol(iElem)%Root%NodeDepth = 1
    ElemNodeVol(iElem)%Root%MidPoint(1:3) = (/0.0,0.0,0.0/)
    ! Start of the recursive routine, which will descend further down the octree until the aforementioned criteria are fulfilled
      IF (BGKMovingAverage) THEN
        CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root, ElemNodeAveraging(iElem)%Root)
      ELSE
        CALL AddBGKOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
      END IF
    DEALLOCATE(TreeNode%MappedPartStates)
  ELSE ! No octree refinement: Call of the respective collision operator
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
#else
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart,ElemVolume_Shared(CNElemID), &
                ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
    END IF
#endif
  END IF ! nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)
END IF

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
      BGK_QualityFacSamp(8,iElem) = BGK_QualityFacSamp(8,iElem) + BGK_Viscosity
      BGK_QualityFacSamp(9,iElem) = BGK_QualityFacSamp(9,iElem) + BGK_ThermalConductivity
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
USE MOD_BGK_CollOperator      ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_CalcSubNodeVolumes3D, OCTANTCUBEID, OCTANTCUBEMIDPOINT
USE MOD_BGK_Vars              ,ONLY: BGKMinPartPerCell,tNodeAverage, BGKMovingAverage
USE MOD_FP_CollOperator       ,ONLY: FP_CollisionOperator
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                                 :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER                 :: TreeNode
CLASS(tNodeVolume),INTENT(INOUT)                    :: NodeVol
CLASS(tNodeAverage),INTENT(INOUT), OPTIONAL         :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, iLoop, iPartIndx, iLoop2, ChildNodeID
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
IF (ABS(NodeVol%MidPoint(1)) .EQ. 1.0) THEN
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
    ChildNodeID = OCTANTCUBEID(NodeVol%MidPoint(:),TreeNode%MappedPartStates(:,iPart))
    PartNumChildNode(ChildNodeID) = PartNumChildNode(ChildNodeID) + 1
    iPartIndx_ChildNode(ChildNodeID,PartNumChildNode(ChildNodeID)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(ChildNodeID),ChildNodeID) = TreeNode%MappedPartStates(1:3,iPart)
  END DO

! Check if any of the subcells has less particles than the limit, if so perform a recombination of cells (3.)
IF(ANY(PartNumChildNode.LT.BGKMinPartPerCell)) CombineChildNodes = .TRUE.

! 2.) Calculate the subcell volume (if necessary)
IF(.NOT.ASSOCIATED(NodeVol%SubNode)) THEN
  CALL DSMC_CalcSubNodeVolumes3D(iElem, TreeNode%NodeDepth, ElemNodeVol(iElem)%Root)
END IF

DO iLoop = 1, 8
  NodeVolumeTemp(iLoop) = NodeVol%SubNode(iLoop)%Volume
END DO
!IF (BGKMovingAverage) THEN
! IF (.NOT.ASSOCIATED(Averaging%SubNode)) THEN
!   CALL BGK_AllocateAveragingNode(Averaging)
! END IF
!END IF

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
    TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
    ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further octree division)
    IF (BGKMovingAverage) THEN
      CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop), Averaging%SubNode(iLoop))
    ELSE
      CALL AddBGKOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop))
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
     CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
           PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), Averaging%SubNode(iLoop)%AverageValues(:))

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
USE MOD_BGK_Vars,               ONLY: tNodeAverage,BGKCollModel
USE MOD_Particle_Vars,          ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CLASS(tNodeAverage),INTENT(INOUT)    :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: numSubNodes, iNode
!===================================================================================================================================
IF (Symmetry%Order.EQ.3) THEN
  numSubNodes = 8
ELSE IF (Symmetry%Order.EQ.2) THEN
  numSubNodes = 4
END IF
ALLOCATE(Averaging%SubNode(numSubNodes))
DO iNode =1, numSubNodes
  IF (BGKCollModel.EQ.1) THEN
    ALLOCATE(Averaging%SubNode(iNode)%AverageValues(10))
  ELSE IF (BGKCollModel.EQ.2) THEN
    ALLOCATE(Averaging%SubNode(iNode)%AverageValues(8))
  ELSE IF (BGKCollModel.EQ.3) THEN
    ALLOCATE(Averaging%SubNode(iNode)%AverageValues(5))
  END IF
  Averaging%SubNode(iNode)%CorrectStep = Averaging%CorrectStep
  Averaging%SubNode(iNode)%AverageValues = Averaging%AverageValues
END DO

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
USE MOD_Particle_Vars           ,ONLY: PEM, Species,WriteMacroVolumeValues, usevMPF, VirtMergedCells, DoVirtualCellMerge
USE MOD_BGK_CollOperator        ,ONLY: BGK_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKMinPartPerCell,BGKSplittingDens,BGKMovingAverage,ElemNodeAveraging
USE MOD_FP_CollOperator         ,ONLY: FP_CollisionOperator
USE MOD_BGK_Vars                ,ONLY: BGKInitDone,BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor
USE MOD_BGK_Vars                ,ONLY: BGK_QualityFacSamp, BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_BGK_Vars                ,ONLY: BGK_Viscosity, BGK_ThermalConductivity
USE MOD_FPFlow_Vars             ,ONLY: FPInitDone, FP_PrandtlNumber, FP_QualityFacSamp
USE MOD_FPFlow_Vars             ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars               ,ONLY: offsetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
#if PP_TimeDiscMethod==300
USE MOD_Particle_Vars           ,ONLY: PartState
#else
USE MOD_Particle_Vars           ,ONLY: LastPartPos
#endif /*PP_TimeDiscMethod==300*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, nPart, CNElemID, nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem
REAL                          :: Dens, partWeight, totalWeight
LOGICAL                       :: DoMergedCell
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.; BGK_Viscosity=0.; BGK_ThermalConductivity=0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF
DoMergedCell = .FALSE.
! Skip cell if number of particles is less than 2, create particle list (iPartIndx_Node) and sum-up bulk velocity
nPart = PEM%pNumber(iElem)
IF (DoVirtualCellMerge) THEN
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    nPartMerged = nPart
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
    END DO
    IF (nPartMerged.LE.1) RETURN
    NULLIFY(TreeNode)
    ALLOCATE(TreeNode)
    ALLOCATE(TreeNode%iPartIndx_Node(nPartMerged))
    iPart = PEM%pStart(iElem)
    iLoopLoc = 0
    DO iLoop = 1, nPart
      iLoopLoc = iLoopLoc + 1
      TreeNode%iPartIndx_Node(iLoopLoc) = iPart
      iPart = PEM%pNext(iPart)
    END DO
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      locElem = VirtMergedCells(iElem)%MergedCellID(iMergeElem)
      nPartLoc = PEM%pNumber(locElem)
      iPart = PEM%pStart(locElem)
      DO iLoop = 1, nPartLoc
        iLoopLoc = iLoopLoc + 1
        TreeNode%iPartIndx_Node(iLoopLoc) = iPart
        iPart = PEM%pNext(iPart)
      END DO
    END DO
    DoMergedCell = .TRUE.
  END IF  
ELSE IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

IF (DoMergedCell) THEN
#if (PP_TimeDiscMethod==300)
  CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged, VirtMergedCells(iELem)%MergedVolume)
#else
  IF (BGKMovingAverage) THEN
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged,VirtMergedCells(iELem)%MergedVolume, &
            ElemNodeAveraging(iElem)%Root%AverageValues(:))
  ELSE
    CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPartMerged, VirtMergedCells(iELem)%MergedVolume)
  END IF
#endif
ELSE
  CNElemID = GetCNElemID(iElem+offSetElem)

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
    Dens = totalWeight / ElemVolume_Shared(CNElemID)
  ELSE
    Dens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
  END IF

  ! The quadtree refinement is performed if either the particle number or number density is above a user-given limit
  IF(nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)) THEN
    ALLOCATE(TreeNode%MappedPartStates(1:2,1:nPart))
    TreeNode%PNum_Node = nPart
    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
    DO iLoop = 1, nPart
#if PP_TimeDiscMethod==300
      CALL GeoCoordToMap2D(PartState(1:2,iPart), TreeNode%MappedPartStates(1:2,iLoop), iElem)
#else
      ! Attention: LastPartPos is the reference position here
      TreeNode%MappedPartStates(1:2,iLoop)=LastPartPos(1:2,iPart)
#endif /*PP_TimeDiscMethod==300*/
      iPart = PEM%pNext(iPart)
    END DO
    TreeNode%NodeDepth = 1
    ElemNodeVol(iElem)%Root%NodeDepth = 1
    ElemNodeVol(iElem)%Root%MidPoint(1:3) = (/0.0,0.0,0.0/)
    ! Start of the recursive routine, which will descend further down the quadtree until the aforementioned criteria are fulfilled
    IF (BGKMovingAverage) THEN
      CALL AddBGKQuadtreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root, ElemNodeAveraging(iElem)%Root)
    ELSE
      CALL AddBGKQuadtreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
    END IF
    DEALLOCATE(TreeNode%MappedPartStates)
  ELSE ! No quadtree refinement: Call of the respective collision operator
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(TreeNode%iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
#else
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, &
               ElemVolume_Shared(CNElemID), ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(TreeNode%iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
    END IF
#endif
  END IF ! nPart.GE.(2.*BGKMinPartPerCell).AND.(Dens.GT.BGKSplittingDens)
END IF

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
      BGK_QualityFacSamp(8,iElem) = BGK_QualityFacSamp(8,iElem) + BGK_Viscosity
      BGK_QualityFacSamp(9,iElem) = BGK_QualityFacSamp(9,iElem) + BGK_ThermalConductivity
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


RECURSIVE SUBROUTINE AddBGKQuadtreeNode(TreeNode, iElem, NodeVol, Averaging)
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
USE MOD_BGK_CollOperator      ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_CalcSubNodeVolumes2D, QUADCUBEMIDPOINT
USE MOD_BGK_Vars              ,ONLY: BGKMinPartPerCell,tNodeAverage, BGKMovingAverage
USE MOD_FP_CollOperator       ,ONLY: FP_CollisionOperator
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                                 :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER                 :: TreeNode
CLASS(tNodeVolume),INTENT(INOUT)                    :: NodeVol
CLASS(tNodeAverage),INTENT(INOUT), OPTIONAL         :: Averaging
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, iLoop, iPartIndx, iLoop2
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
IF (ABS(NodeVol%MidPoint(1)) .EQ. 1.0) THEN
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
  IF ((TreeNode%MappedPartStates(1,iPart).GE.NodeVol%MidPoint(1)) &
      .AND.(TreeNode%MappedPartStates(2,iPart).LE.NodeVol%MidPoint(2))) THEN
    PartNumChildNode(1) = PartNumChildNode(1) + 1
    iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(1),1) = TreeNode%MappedPartStates(1:2,iPart)
  ELSE IF(TreeNode%MappedPartStates(1,iPart).GE.NodeVol%MidPoint(1)) THEN
    PartNumChildNode(2) = PartNumChildNode(2) + 1
    iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    MappedPart_ChildNode(1:2,PartNumChildNode(2),2) = TreeNode%MappedPartStates(1:2,iPart)
  ELSE IF(TreeNode%MappedPartStates(2,iPart).GE.NodeVol%MidPoint(2)) THEN
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
IF (.NOT.ASSOCIATED(NodeVol%SubNode)) THEN
  CALL DSMC_CalcSubNodeVolumes2D(iElem, TreeNode%NodeDepth, ElemNodeVol(iElem)%Root)
END IF

DO iLoop = 1, 4
  NodeVolumeTemp(iLoop) = NodeVol%SubNode(iLoop)%Volume
END DO

IF (BGKMovingAverage) THEN
  IF (.NOT.ASSOCIATED(Averaging%SubNode)) THEN
    CALL BGK_AllocateAveragingNode(Averaging)
  END IF
END IF

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
    TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
    ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further quadtree division)
    IF (BGKMovingAverage) THEN
      CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop), Averaging%SubNode(iLoop))
    ELSE
      CALL AddBGKQuadtreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop))
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
    CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
       PartNumChildNode(iLoop), NodeVolumeTemp(iLoop), Averaging%SubNode(iLoop)%AverageValues(:))
  ELSE
    CALL BGK_CollisionOperator(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
          PartNumChildNode(iLoop), NodeVolumeTemp(iLoop))
  END IF
#endif
  END IF
END DO

END SUBROUTINE AddBGKQuadtreeNode

END MODULE MOD_BGK_Adaptation
