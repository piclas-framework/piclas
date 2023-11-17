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

INTERFACE DSMC_pairing_dotree
  MODULE PROCEDURE DSMC_pairing_dotree
END INTERFACE

INTERFACE DSMC_pairing_standard
  MODULE PROCEDURE DSMC_pairing_standard
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_pairing_standard, DSMC_pairing_octree, DSMC_init_octree, DSMC_pairing_quadtree, DSMC_CalcSubNodeVolumes3D
PUBLIC :: DSMC_CalcSubNodeVolumes2D, GeoCoordToMap2D, DSMC_pairing_dotree, OCTANTCUBEID, QUADCUBEMIDPOINT, OCTANTCUBEMIDPOINT
!===================================================================================================================================

CONTAINS


SUBROUTINE DSMC_pairing_standard(iElem)
!===================================================================================================================================
!> Collisions within a single cell: creates mapping between particle number within the cell to the particle index on the processor.
!> Calls the pairing and collision routine, where the actual pairing with a random partner or nearest neighbour as well as the
!> collision procedure is performed.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PEM, VirtMergedCells, DoVirtualCellMerge
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, nPart, nPartMerged, iMergeElem, iLoopLoc, locElem, nPartLoc
INTEGER, ALLOCATABLE          :: iPartIndx(:)                 !< List of particles in the cell required for pairing
REAL                          :: elemVolume
!===================================================================================================================================

nPart = PEM%pNumber(iElem)
IF (DoVirtualCellMerge) THEN
  ! 1.) Create particle index list for pairing in the case of virtually merged cells. So, all particles from the merged cells are
  !   used for the pairing and the collisions.
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  nPartMerged = nPart
  DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
    nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
  END DO
  ALLOCATE(iPartIndx(nPartMerged))
  iPart = PEM%pStart(iElem)
  iLoopLoc = 0
  DO iLoop = 1, nPart
    iLoopLoc = iLoopLoc + 1
    iPartIndx(iLoopLoc) = iPart
    iPart = PEM%pNext(iPart)
  END DO
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      locElem = VirtMergedCells(iElem)%MergedCellID(iMergeElem)
      nPartLoc = PEM%pNumber(locElem)
      iPart = PEM%pStart(locElem)
      DO iLoop = 1, nPartLoc
        iLoopLoc = iLoopLoc + 1
        iPartIndx(iLoopLoc) = iPart
        iPart = PEM%pNext(iPart)
      END DO
    END DO
    elemVolume = VirtMergedCells(iELem)%MergedVolume
  ELSE
    elemVolume = ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
  END IF
ELSE
  nPartMerged = nPart
  ALLOCATE(iPartIndx(nPart))
  iPartIndx = 0
  ! 1.) Create particle index list for pairing
  !     Using PEM%pStart to get the first particle index based on the element index and PEM%pNext to get the next particle index based
  !     on the previous particle index. This mapping is done in the UpdateNextFreePosition routine.
  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    iPartIndx(iLoop) = iPart
    ! Choose next particle in the element
    iPart = PEM%pNext(iPart)
  END DO
  elemVolume = ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
END IF

! 2.) Perform pairing (random pairing or nearest neighbour pairing) and collision (including the decision for a reaction/relaxation)
CALL PerformPairingAndCollision(iPartIndx, nPartMerged, iElem , elemVolume)
DEALLOCATE(iPartIndx)

END SUBROUTINE DSMC_pairing_standard


SUBROUTINE DSMC_pairing_octree(iElem)
!===================================================================================================================================
!> Collisions within a cell/subcell using a recursive octree mesh refinement based on the mean free path criterion for DSMC and/or
!> the number of particles per cell (to accelerate a potential nearest neighbour search). Creates mapping between particle number
!> within the cell/subcell to the particle index on the processor. Calls the pairing and collision routine, where the actual pairing
!> with a random partner or nearest neighbour as well as the collision procedure is performed.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, DSMC, ElemNodeVol
USE MOD_Particle_Vars           ,ONLY: PEM, nSpecies, PartSpecies,PartPosRef,LastPartPos, VirtMergedCells, DoVirtualCellMerge
USE MOD_Particle_Tracking_vars  ,ONLY: TrackingMethod
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
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
INTEGER                       :: iPart, iLoop, nPart, CNElemID, nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem
REAL                          :: SpecPartNum(nSpecies)
TYPE(tTreeNode), POINTER      :: TreeNode
REAL                          :: elemVolume
LOGICAL                       :: DoMergedCell
!===================================================================================================================================

SpecPartNum = 0.
nPart = PEM%pNumber(iElem)
DoMergedCell = .FALSE.
IF (DoVirtualCellMerge) THEN
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    nPartMerged = nPart
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
    END DO
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
END IF

IF (DoMergedCell) THEN
  CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPartMerged, iElem, VirtMergedCells(iELem)%MergedVolume)
ELSE
  CNElemID = GetCNElemID(iElem+offSetElem)

  NULLIFY(TreeNode)

  ALLOCATE(TreeNode)
  ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing
  TreeNode%iPartIndx_Node(1:nPart) = 0

  ! 1.) Create particle index list for pairing
  !     Using PEM%pStart to get the first particle index based on the element index and PEM%pNext to get the next particle index based
  !     on the previous particle index. This mapping is done in the UpdateNextFreePosition routine.
  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    TreeNode%iPartIndx_Node(iLoop) = iPart
    iPart = PEM%pNext(iPart)
    ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
    SpecPartNum(PartSpecies(TreeNode%iPartIndx_Node(iLoop))) = &
              SpecPartNum(PartSpecies(TreeNode%iPartIndx_Node(iLoop))) + GetParticleWeight(TreeNode%iPartIndx_Node(iLoop))
  END DO

  elemVolume=ElemVolume_Shared(CNElemID)

  ! 2.) Octree cell refinement algorithm
  DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum, SUM(SpecPartNum), elemVolume)
  ! Octree can only performed if nPart is greater than the defined value (default = 28 for 2D/axisymmetric or = 50 for 3D)
  IF(nPart.GE.DSMC%PartNumOctreeNodeMin) THEN
    ! Additional check afterwards if nPart is greater than PartNumOctreeNode (default = 40  for 2D/axisymmetric or = 80 for 3D)
    ! or the mean free path is less than the side length of a cube (approximation) with same volume as the actual cell
    IF((DSMC%MeanFreePath.LT.ElemCharLength_Shared(CNElemID)) .OR.(nPart.GT.DSMC%PartNumOctreeNode)) THEN
      ALLOCATE(TreeNode%MappedPartStates(1:3,1:nPart))
      TreeNode%PNum_Node = nPart
      iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
      IF (TrackingMethod.EQ.REFMAPPING) THEN
        DO iLoop = 1, nPart
          TreeNode%MappedPartStates(1:3,iLoop)=PartPosRef(1:3,iPart)
          iPart = PEM%pNext(iPart)
        END DO
      ELSE ! position in reference space [-1,1] has to be computed
        DO iLoop = 1, nPart
          ! Attention: LastPartPos is the reference position here
          TreeNode%MappedPartStates(1:3,iLoop) = LastPartPos(1:3,iPart)
          iPart = PEM%pNext(iPart)
        END DO
      END IF ! TrackingMethod.EQ.REFMAPPING
      TreeNode%NodeDepth = 1
      ElemNodeVol(iElem)%Root%NodeDepth = 1
      ElemNodeVol(iElem)%Root%MidPoint(1:3) = (/0.0,0.0,0.0/)
      CALL AddOctreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
      DEALLOCATE(TreeNode%MappedPartStates)
    ELSE
      CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, ElemVolume_Shared(CNElemID))
    END IF
  ELSE
    CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, ElemVolume_Shared(CNElemID))
  END IF
END IF

DEALLOCATE(TreeNode%iPartIndx_Node)
DEALLOCATE(TreeNode)

END SUBROUTINE DSMC_pairing_octree


SUBROUTINE FindNearestNeigh3D(iPartIndx_Node, nPart, iPair)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars              ,ONLY: Coll_pData
  USE MOD_Particle_Vars          ,ONLY: PartState
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(INOUT)                  :: nPart
  INTEGER, INTENT(IN)                     :: iPair
  INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart1, iPart2, iLoop
REAL                          :: Dist1, Dist2, iRan
!===================================================================================================================================

CALL RANDOM_NUMBER(iRan)
iPart1 = 1 + INT(nPart * iRan)
Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
nPart = nPart - 1
iPart2 = 1
Dist1 = (PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iPart2)))**2 &
      + (PartState(2,Coll_pData(iPair)%iPart_p1) - PartState(2,iPartIndx_Node(iPart2)))**2 &
      + (PartState(3,Coll_pData(iPair)%iPart_p1) - PartState(3,iPartIndx_Node(iPart2)))**2
DO iLoop = 2, nPart
  Dist2 = (PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iLoop)))**2 &
        + (PartState(2,Coll_pData(iPair)%iPart_p1) - PartState(2,iPartIndx_Node(iLoop)))**2 &
        + (PartState(3,Coll_pData(iPair)%iPart_p1) - PartState(3,iPartIndx_Node(iLoop)))**2
  IF (Dist2.LT.Dist1) THEN
    iPart2 = iLoop
    Dist1 = Dist2
  END IF
END DO
Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
nPart = nPart - 1

END SUBROUTINE FindNearestNeigh3D


SUBROUTINE FindNearestNeigh2D(iPartIndx_Node, nPart, iPair)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf
USE MOD_Particle_Vars,          ONLY: PartState
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)          :: nPart
INTEGER, INTENT(IN)             :: iPair
INTEGER, INTENT(INOUT)          :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart1, iPart2, iLoop, loopStart
REAL                            :: Dist1, Dist2, iRan
!===================================================================================================================================

loopStart = 0
CALL RANDOM_NUMBER(iRan)
iPart1 = 1 + INT(nPart * iRan)
Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
nPart = nPart - 1
iPart2 = 1
IF (CollInf%ProhibitDoubleColl) THEN
  IF (nPart.GT.1) THEN
    IF (iPartIndx_Node(iPart2).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
      iPart2 = 2
      loopStart = 1
    END IF
  END IF
END IF
Dist1 = (PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iPart2)))**2 &
      + (PartState(2,Coll_pData(iPair)%iPart_p1) - PartState(2,iPartIndx_Node(iPart2)))**2
DO iLoop = 2 + loopStart, nPart
  IF (CollInf%ProhibitDoubleColl) THEN
      IF (iPartIndx_Node(iLoop).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
        CYCLE
      END IF
  END IF
  Dist2 = (PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iLoop)))**2 &
        + (PartState(2,Coll_pData(iPair)%iPart_p1) - PartState(2,iPartIndx_Node(iLoop)))**2
  IF (Dist2.LT.Dist1) THEN
    iPart2 = iLoop
    Dist1 = Dist2
  END IF
END DO
Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
nPart = nPart - 1

END SUBROUTINE FindNearestNeigh2D


SUBROUTINE FindNearestNeigh1D(iPartIndx_Node, nPart, iPair)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf
USE MOD_Particle_Vars,          ONLY: PartState
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)          :: nPart
INTEGER, INTENT(IN)             :: iPair
INTEGER, INTENT(INOUT)          :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart1, iPart2, iLoop, loopStart
REAL                            :: Dist1, Dist2, iRan
!===================================================================================================================================

loopStart = 0
CALL RANDOM_NUMBER(iRan)
iPart1 = 1 + INT(nPart * iRan)
Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
nPart = nPart - 1
iPart2 = 1
IF (CollInf%ProhibitDoubleColl) THEN
  IF (nPart.GT.1) THEN
    IF (iPartIndx_Node(iPart2).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
      iPart2 = 2
      loopStart = 1
    END IF
  END IF
END IF
Dist1 = ABS(PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iPart2)))
DO iLoop = 2 + loopStart, nPart
  IF (CollInf%ProhibitDoubleColl) THEN
      IF (iPartIndx_Node(iLoop).EQ.CollInf%OldCollPartner(Coll_pData(iPair)%iPart_p1)) THEN
        CYCLE
      END IF
  END IF
  Dist2 = ABS(PartState(1,Coll_pData(iPair)%iPart_p1) - PartState(1,iPartIndx_Node(iLoop)))
  IF (Dist2.LT.Dist1) THEN
    iPart2 = iLoop
    Dist1 = Dist2
  END IF
END DO
Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
nPart = nPart - 1

END SUBROUTINE FindNearestNeigh1D


SUBROUTINE PerformPairingAndCollision(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
!> Main pairing and collision routine performed in a cell/subcell: calls the statistical and nearest neighbour pairing routines
!> as well as the collision probability calculation and actual collision execution
!> 1). Reset collision and pair specific variables
!> 2.) Calculate cell/subcell local variables and count the number of particles per species
!> 3.) Perform the particle pairing (statistical or nearest neighbour) and determine the relative velocity
!> 4.) Perform additional operations for radial weighting and chemistry (AFTER pairing)
!> 5). Calculate the collision probability and perform the collision (if required)
!> 6.) Calculate the mean free path and the mean collision separation distance within a cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_CollisionProb      ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Collis             ,ONLY: DSMC_perform_collision
USE MOD_DSMC_Vars               ,ONLY: Coll_pData,CollInf,CollisMode,PartStateIntEn,ChemReac,DSMC,RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: SelectionProc, useRelaxProbCorrFactor, iPartIndx_NodeNewElecRelax, newElecRelaxParts
USE MOD_DSMC_Vars               ,ONLY: iPartIndx_NodeElecRelaxChem,nElecRelaxChemParts
USE MOD_Particle_Vars           ,ONLY: PartSpecies, nSpecies, PartState, WriteMacroVolumeValues, UseVarTimeStep, Symmetry, usevMPF
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, time
USE MOD_DSMC_Analyze            ,ONLY: CalcGammaVib, CalcInstantTransTemp, CalcMeanFreePath, CalcInstantElecTempXi
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_DSMC_Relaxation         ,ONLY: CalcMeanVibQuaDiatomic,SumVibRelaxProb
USE MOD_DSMC_Symmetry           ,ONLY: DSMC_2D_TreatIdenticalParticles
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_InsertParticles, AD_DeleteParticles
USE MOD_DSMC_ElectronicModel    ,ONLY: LT_ElectronicEnergyExchange, LT_ElectronicExc_ConstructPartList
USE MOD_DSMC_ElectronicModel    ,ONLY: CalcProbCorrFactorElec, LT_ElectronicEnergyExchangeChem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)              :: NodeVolume
INTEGER, INTENT(IN)           :: iElem
INTEGER, INTENT(INOUT)        :: PartNum
INTEGER, INTENT(INOUT), TARGET:: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nPair, iPair, iPart, nPart, TotalPartNum, nPartElecRelac
INTEGER                       :: cSpec1, cSpec2, iCase
REAL                          :: iRan
INTEGER, ALLOCATABLE, TARGET  :: iPartIndx_NodeTotalAmbi(:)
INTEGER, POINTER              :: iPartIndx_NodeTotal(:)
INTEGER, ALLOCATABLE          :: iPartIndx_NodeTotalAmbiDel(:)
INTEGER, ALLOCATABLE          :: iPartIndx_NodeTotalElecExc(:), iPartIndx_NodeTotalElecRelax(:)
!===================================================================================================================================

! 0). Ambipolar Diffusion
IF (DSMC%DoAmbipolarDiff) THEN
  nPart = PartNum
  CALL AD_InsertParticles(iPartIndx_Node,nPart, iPartIndx_NodeTotalAmbi, TotalPartNum)
  ALLOCATE(iPartIndx_NodeTotalAmbiDel(1:TotalPartNum))
  iPartIndx_NodeTotalAmbiDel(1:TotalPartNum) = iPartIndx_NodeTotalAmbi(1:TotalPartNum)
  nPart = TotalPartNum
  iPartIndx_NodeTotal => iPartIndx_NodeTotalAmbi
ELSE
  nPart = PartNum
  TotalPartNum = PartNum
  iPartIndx_NodeTotal => iPartIndx_Node
END IF

IF (DSMC%ElectronicModel.EQ.4) THEN
  ALLOCATE(iPartIndx_NodeTotalElecRelax(TotalPartNum))
  iPartIndx_NodeTotalElecRelax = iPartIndx_NodeTotal
  IF (CollisMode.EQ.3) THEN
    newElecRelaxParts = 0; nElecRelaxChemParts = 0
    ALLOCATE(iPartIndx_NodeNewElecRelax(2*PartNum), iPartIndx_NodeElecRelaxChem(2*PartNum))
  END IF
END IF

! 1). Reset collision and pair specific variables
nPair = INT(TotalPartNum/2)
CollInf%Coll_SpecPartNum = 0.
CollInf%Coll_CaseNum = 0
ALLOCATE(Coll_pData(nPair))
Coll_pData%Ec=0

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) CollInf%SumPairMPF = 0.

! 2.) Calculate cell/subcell local variables and count the number of particles per species
DO iPart = 1, TotalPartNum
  CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_NodeTotal(iPart))) = CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_NodeTotal(iPart))) &
                                                                  + GetParticleWeight(iPartIndx_NodeTotal(iPart))
END DO

IF (CollisMode.EQ.3) THEN
  ChemReac%RecombParticle = 0
  ChemReac%nPairForRec = 0
  ChemReac%LastPairForRec = 0
! Determination of the mean vibrational energy for the cell
  ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
  DO iPart = 1, TotalPartNum
    ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_NodeTotal(iPart)))=ChemReac%MeanEVib_PerIter(PartSpecies(iPartIndx_NodeTotal(iPart))) &
      + PartStateIntEn(1,iPartIndx_NodeTotal(iPart)) * GetParticleWeight(iPartIndx_NodeTotal(iPart))
  END DO
  CALL CalcMeanVibQuaDiatomic()
END IF

IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors &
.OR.(useRelaxProbCorrFactor.AND.(CollisMode.GT.1)).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.4)) THEN
  ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
  !           relaxation probability (CalcGammaVib)
  ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
  ! 3. Case: Temperature required for the mean free path with the VHS model
  ! 4. Case: Needed to calculate the correction factor
  CALL CalcInstantTransTemp(iPartIndx_NodeTotal,TotalPartNum)
  IF ((DSMC%ElectronicModel.EQ.2).OR.useRelaxProbCorrFactor) CALL CalcInstantElecTempXi(iPartIndx_NodeTotal,TotalPartNum)
  IF((SelectionProc.EQ.2).OR.(useRelaxProbCorrFactor)) CALL CalcGammaVib()
  IF (useRelaxProbCorrFactor.AND.(DSMC%ElectronicModel.EQ.1)) CALL CalcProbCorrFactorElec()
END IF

IF (CollInf%ProhibitDoubleColl.AND.(nPair.EQ.1)) THEN
! Do not get stuck in an endless loop if only two particles/one pair are present in the cell
  CollInf%OldCollPartner(iPartIndx_NodeTotal(1)) = 0
  CollInf%OldCollPartner(iPartIndx_NodeTotal(2)) = 0
END IF

! 3.) Perform the particle pairing (statistical or nearest neighbour) and determine the relative velocity
DO iPair = 1, nPair
  IF(DSMC%UseNearestNeighbour) THEN
    IF(Symmetry%Order.EQ.3) THEN
      CALL FindNearestNeigh3D(iPartIndx_NodeTotal, nPart, iPair)
    ELSE IF (Symmetry%Order.EQ.2) THEN
      CALL FindNearestNeigh2D(iPartIndx_NodeTotal, nPart, iPair)
    ELSE
      CALL FindNearestNeigh1D(iPartIndx_NodeTotal, nPart, iPair)
    END IF
  ELSE
    CALL FindRandomPartner(iPartIndx_NodeTotal, nPart, iPair, nPair)
  END IF

  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  ! Summation of the average weighting factor of the collision pairs for each case (AA, AB, BB)
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
    CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + (GetParticleWeight(Coll_pData(iPair)%iPart_p1) &
                                                        + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
  END IF

  CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
  Coll_pData(iPair)%CRela2 =  (PartState(4,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(4,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(5,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(5,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(6,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(6,Coll_pData(iPair)%iPart_p2))**2
  Coll_pData(iPair)%PairType = iCase
  Coll_pData(iPair)%NeedForRec = .FALSE.
END DO

! 4.) Perform additional operations for chemistry (AFTER pairing)
! If a third particle is required of a recombination, the last particle due to uneven nPart is used
IF(CollisMode.EQ.3) THEN
  IF(nPart.EQ.1) ChemReac%RecombParticle = iPartIndx_NodeTotal(1)
END IF
! Resetting the previous collision partner of the remaining particle due to uneven nPart
IF (CollInf%ProhibitDoubleColl.AND.(nPart.EQ.1)) CollInf%OldCollPartner(iPartIndx_NodeTotal(1)) = 0

! 5). Calculate the collision probability and perform the collision (if required)
DO iPair = 1, nPair
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    CALL SumVibRelaxProb(iPair)
    ! 2D axisymmetric with radial weighting: split up pairs of identical particles
    IF(RadialWeighting%DoRadialWeighting) CALL DSMC_2D_TreatIdenticalParticles(iPair, nPair, nPart, iElem, iPartIndx_NodeTotal)
    ! Calculate the collision probability and test it against a random number
    CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.GE.iRan) THEN
      CALL DSMC_perform_collision(iPair,iElem, NodeVolume, TotalPartNum)
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

! 6.) Calculate the mean free path and the mean collision separation distance within a cell
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path with VHS model and the current translational temperature in the cell
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum), REAL(SUM(CollInf%Coll_SpecPartNum)), NodeVolume, &
                                          DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the maximum MCS/MFP for the cell
    IF((DSMC%CollSepCount.GT.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = &
                                                    MAX(DSMC%MCSoverMFP,(DSMC%CollSepDist/DSMC%CollSepCount)/DSMC%MeanFreePath)
  END IF
END IF

DEALLOCATE(Coll_pData)

IF (DSMC%ElectronicModel.EQ.4) THEN
  IF (CollisMode.EQ.3) THEN
    CALL LT_ElectronicEnergyExchangeChem(iPartIndx_NodeElecRelaxChem, nElecRelaxChemParts)
    CALL LT_ElectronicExc_ConstructPartList(iPartIndx_NodeTotalElecRelax, iPartIndx_NodeTotalElecExc,  TotalPartNum, nPartElecRelac)
    CALL LT_ElectronicEnergyExchange(iPartIndx_NodeTotalElecExc, nPartElecRelac, NodeVolume)
    DEALLOCATE(iPartIndx_NodeTotalElecExc, iPartIndx_NodeNewElecRelax, iPartIndx_NodeElecRelaxChem)
  ELSE
    CALL LT_ElectronicEnergyExchange(iPartIndx_NodeTotalElecRelax, TotalPartNum, NodeVolume)
  END IF
END IF

IF (DSMC%DoAmbipolarDiff) THEN
  CALL AD_DeleteParticles(iPartIndx_NodeTotalAmbiDel,TotalPartNum)
END IF

END SUBROUTINE PerformPairingAndCollision


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
  TYPE(tTreeNode),INTENT(INOUT), POINTER  :: TreeNode
  CLASS(tNodeVolume),INTENT(INOUT)        :: NodeVol
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, iPartIndx, SpecPartNum(nSpecies), childNodeID
  INTEGER, ALLOCATABLE          :: iPartIndx_ChildNode(:,:)
  REAL, ALLOCATABLE             :: MappedPart_ChildNode(:,:,:)
  INTEGER                       :: PartNumChildNode(1:8)
!===================================================================================================================================

  ALLOCATE(iPartIndx_ChildNode(1:8,TreeNode%PNum_Node))
  ALLOCATE(MappedPart_ChildNode(1:3,TreeNode%PNum_Node,1:8))
  PartNumChildNode(:) = 0
  IF (ABS(NodeVol%MidPoint(1)) .EQ. 1.0) THEN
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
    ChildNodeID = OCTANTCUBEID(NodeVol%MidPoint(:),TreeNode%MappedPartStates(:,iPart))
    PartNumChildNode(ChildNodeID) = PartNumChildNode(ChildNodeID) + 1
    iPartIndx_ChildNode(ChildNodeID,PartNumChildNode(ChildNodeID)) = iPartIndx
    MappedPart_ChildNode(1:3,PartNumChildNode(ChildNodeID),ChildNodeID) = TreeNode%MappedPartStates(1:3,iPart)
  END DO

  IF(.NOT.ASSOCIATED(NodeVol%SubNode)) THEN
    CALL DSMC_CalcSubNodeVolumes3D(iElem, TreeNode%NodeDepth, ElemNodeVol(iElem)%Root)
  END IF

  DO iLoop = 1, 8
    IF (PartNumChildNode(iLoop).GT.1) THEN
      ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
      SpecPartNum = 0
      DO iPart = 1, PartNumChildNode(iLoop)
        SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) = &
          SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart))) + 1
      END DO
      DSMC%MeanFreePath = CalcMeanFreePath(REAL(SpecPartNum),REAL(PartNumChildNode(iLoop)),NodeVol%SubNode(iLoop)%Volume)
    END IF
    ! Octree can only performed if nPart is greater than the defined value (default=20), otherwise nearest neighbour pairing
    IF(PartNumChildNode(iLoop).GE.DSMC%PartNumOctreeNodeMin) THEN
      ! Additional check if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
      ! the side length of a cube (approximation) with same volume as the actual cell -> octree
      IF((DSMC%MeanFreePath.LT.(NodeVol%SubNode(iLoop)%Volume**(1./3.))) &
                                                               .OR.(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode)) THEN
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
        CALL AddOctreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop))
        DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
        DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
        DEALLOCATE(TreeNode%ChildNode)
      ELSE
          CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                  PartNumChildNode(iLoop), iElem, NodeVol%SubNode(iLoop)%Volume)
      END IF
    ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
        CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
                PartNumChildNode(iLoop), iElem, NodeVol%SubNode(iLoop)%Volume)
    END IF
  END DO

END SUBROUTINE AddOctreeNode


SUBROUTINE DSMC_pairing_quadtree(iElem)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, DSMC, ElemNodeVol
USE MOD_Particle_Vars           ,ONLY: PEM, nSpecies, PartSpecies, LastPartPos,VirtMergedCells, DoVirtualCellMerge
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
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
INTEGER                       :: iPart, iLoop, nPart, CNElemID, nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem
REAL                          :: SpecPartNum(nSpecies), Volume
TYPE(tTreeNode), POINTER      :: TreeNode
LOGICAL                       :: DoMergedCell
!===================================================================================================================================

CNElemID = GetCNElemID(iElem+offSetElem)
Volume = ElemVolume_Shared(CNElemID)
SpecPartNum = 0.
nPart = PEM%pNumber(iElem)
DoMergedCell = .FALSE.
IF (DoVirtualCellMerge) THEN
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    nPartMerged = nPart
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
    END DO
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
END IF

IF (DoMergedCell) THEN
  CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPartMerged, iElem, VirtMergedCells(iELem)%MergedVolume)
ELSE

  NULLIFY(TreeNode)
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
    IF((DSMC%MeanFreePath.LT.ElemCharLength_Shared(CNElemID)).OR.(nPart.GT.DSMC%PartNumOctreeNode)) THEN
      ALLOCATE(TreeNode%MappedPartStates(1:2,1:nPart))
      TreeNode%PNum_Node = nPart
      iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
      DO iLoop = 1, nPart
        ! Attention: LastPartPos is the reference position here, set in timedisc_TimeStep_DSMC.f90 / timedisc_TimeStep_BGK.f90
        TreeNode%MappedPartStates(1:2,iLoop) = LastPartPos(1:2,iPart)
        iPart = PEM%pNext(iPart)
      END DO
      TreeNode%NodeDepth = 1
      ElemNodeVol(iElem)%Root%MidPoint(1:3) = (/0.0,0.0,0.0/)
      ElemNodeVol(iElem)%Root%NodeDepth = 1
      CALL AddQuadTreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
      DEALLOCATE(TreeNode%MappedPartStates)
    ELSE
      CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, Volume)
    END IF
  ELSE
    CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, Volume)
  END IF
END IF

DEALLOCATE(TreeNode%iPartIndx_Node)
DEALLOCATE(TreeNode)

END SUBROUTINE DSMC_pairing_quadtree


RECURSIVE SUBROUTINE AddQuadTreeNode(TreeNode, iElem, NodeVol)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze      ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars         ,ONLY: tTreeNode, DSMC, tNodeVolume, RadialWeighting, CollInf
USE MOD_Particle_Vars     ,ONLY: nSpecies, PartSpecies, UseVarTimeStep, usevMPF
USE MOD_DSMC_Vars         ,ONLY: ElemNodeVol
USE MOD_part_tools        ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                     :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER     :: TreeNode
CLASS(tNodeVolume),INTENT(INOUT)        :: NodeVol
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iPartIndx
INTEGER, ALLOCATABLE          :: iPartIndx_ChildNode(:,:)
REAL, ALLOCATABLE             :: MappedPart_ChildNode(:,:,:)
INTEGER                       :: PartNumChildNode(1:4)
REAL                          :: NodeVolumeTemp(1:4), FaceVolumeTemp(1:4), SpecPartNum(nSpecies,1:4), RealParts(1:4)
LOGICAL                       :: ForceNearestNeigh
!===================================================================================================================================
ForceNearestNeigh = .FALSE.
ALLOCATE(iPartIndx_ChildNode(1:4,TreeNode%PNum_Node))
ALLOCATE(MappedPart_ChildNode(1:2,TreeNode%PNum_Node,1:4))
PartNumChildNode(:) = 0
IF (ABS(NodeVol%MidPoint(1)) .EQ. 1.0) THEN
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

IF(.NOT.ASSOCIATED(NodeVol%SubNode)) THEN
  CALL DSMC_CalcSubNodeVolumes2D(iElem, TreeNode%NodeDepth, ElemNodeVol(iElem)%Root)
END IF

DO iLoop = 1, 4
  NodeVolumeTemp(iLoop) = NodeVol%SubNode(iLoop)%Volume
  FaceVolumeTemp(iLoop) = NodeVol%SubNode(iLoop)%Area
END DO

IF(DSMC%MergeSubcells) THEN
  IF (ALL(PartNumChildNode.LT.DSMC%PartNumOctreeNodeMin)) THEN
    ForceNearestNeigh =.TRUE.
!    DO iLoop = 1, 3
!      IF (PartNumChildNode(iLoop).LT.7) THEN
!        DO iPart=1, PartNumChildNode(iLoop)
!          iPartIndx_ChildNode(iLoop+1,PartNumChildNode(iLoop+1)+iPart) = iPartIndx_ChildNode(iLoop,iPart)
!          MappedPart_ChildNode(1:2,PartNumChildNode(iLoop+1)+iPart,iLoop+1) = MappedPart_ChildNode(1:2,iPart,iLoop)
!        END DO
!        PartNumChildNode(iLoop+1) = PartNumChildNode(iLoop+1) + PartNumChildNode(iLoop)
!        PartNumChildNode(iLoop) = 0
!        NodeVolumeTemp(iLoop+1) = NodeVolumeTemp(iLoop+1) + NodeVolumeTemp(iLoop)
!        NodeVolumeTemp(iLoop) = 0.0
!      END IF
!    END DO
!    IF (PartNumChildNode(4).LT.7) THEN
!      DO iPart=1, PartNumChildNode(4)
!        iPartIndx_ChildNode(1,PartNumChildNode(1)+iPart) = iPartIndx_ChildNode(4,iPart)
!        MappedPart_ChildNode(1:2,PartNumChildNode(1)+iPart,1) = MappedPart_ChildNode(1:2,iPart,4)
!      END DO
!      PartNumChildNode(1) = PartNumChildNode(1) + PartNumChildNode(4)
!      PartNumChildNode(4) = 0
!      NodeVolumeTemp(1) = NodeVolumeTemp(1) + NodeVolumeTemp(4)
!      NodeVolumeTemp(4) = 0.0
!    END IF
    IF (ANY(PartNumChildNode.LT.5)) THEN
      DO iPart=1, PartNumChildNode(1)
        iPartIndx_ChildNode(2,PartNumChildNode(2)+iPart) = iPartIndx_ChildNode(1,iPart)
        MappedPart_ChildNode(1:2,PartNumChildNode(2)+iPart,2) = MappedPart_ChildNode(1:2,iPart,1)
      END DO
      PartNumChildNode(2) = PartNumChildNode(2) + PartNumChildNode(1)
      PartNumChildNode(1) = 0
      NodeVolumeTemp(2) = NodeVolumeTemp(2) + NodeVolumeTemp(1)
      NodeVolumeTemp(1) = 0.0
      FaceVolumeTemp(2) = FaceVolumeTemp(2) + FaceVolumeTemp(1)
      FaceVolumeTemp(1) = 0.0
      DO iPart=1, PartNumChildNode(3)
        iPartIndx_ChildNode(4,PartNumChildNode(4)+iPart) = iPartIndx_ChildNode(3,iPart)
        MappedPart_ChildNode(1:2,PartNumChildNode(4)+iPart,4) = MappedPart_ChildNode(1:2,iPart,3)
      END DO
      PartNumChildNode(4) = PartNumChildNode(4) + PartNumChildNode(3)
      PartNumChildNode(3) = 0
      NodeVolumeTemp(4) = NodeVolumeTemp(4) + NodeVolumeTemp(3)
      NodeVolumeTemp(3) = 0.0
      FaceVolumeTemp(4) = FaceVolumeTemp(4) + FaceVolumeTemp(3)
      FaceVolumeTemp(3) = 0.0
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
    IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
      DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop), RealParts(iLoop), NodeVolumeTemp(iLoop))
    ELSE
      DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop),REAL(PartNumChildNode(iLoop)), NodeVolumeTemp(iLoop))
    END IF
    IF((DSMC%MeanFreePath.LT.(FaceVolumeTemp(iLoop)**(1./2.))).OR.(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode)) THEN
      NULLIFY(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
      ALLOCATE(TreeNode%ChildNode%MappedPartStates(1:2,PartNumChildNode(iLoop)))
      TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
      TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
      TreeNode%ChildNode%MappedPartStates(1:2,1:PartNumChildNode(iLoop))= &
                     MappedPart_ChildNode(1:2,1:PartNumChildNode(iLoop),iLoop)
      TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
      ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further quadtree division)
      CALL AddQuadTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop))
      DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
      DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
      DEALLOCATE(TreeNode%ChildNode)
    ELSE
      CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
            PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
    END IF
  ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
    CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
          PartNumChildNode(iLoop), iElem, NodeVolumeTemp(iLoop))
  ELSE IF (CollInf%ProhibitDoubleColl.AND.(PartNumChildNode(iLoop).EQ.1)) THEN
    CollInf%OldCollPartner(iPartIndx_ChildNode(iLoop, 1)) = 0
  END IF
END DO

END SUBROUTINE AddQuadTreeNode


SUBROUTINE DSMC_pairing_dotree(iElem)
!===================================================================================================================================
! Pairing subroutine for dotree, decides whether to create a new dotree node or start particle pairing
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Analyze            ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars               ,ONLY: tTreeNode, DSMC, ElemNodeVol
USE MOD_Particle_Vars           ,ONLY: PEM, PartState, nSpecies, PartSpecies
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
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
INTEGER                       :: iPart, iLoop, nPart, CNElemID
REAL                          :: SpecPartNum(nSpecies), Volume
TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

CNElemID = GetCNElemID(iElem+offSetElem)
Volume = ElemVolume_Shared(CNElemID)
SpecPartNum = 0.

NULLIFY(TreeNode)
nPart = PEM%pNumber(iElem)

ALLOCATE(TreeNode)
ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell necessary for stat pairing
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
  IF((DSMC%MeanFreePath.LT.ElemCharLength_Shared(CNElemID)).OR.(nPart.GT.DSMC%PartNumOctreeNode)) THEN
    ALLOCATE(TreeNode%MappedPartStates(1,1:nPart))
    TreeNode%PNum_Node = nPart
    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
    DO iLoop = 1, nPart
      CALL GeoCoordToMap1D(PartState(1,iPart), TreeNode%MappedPartStates(1,iLoop), iElem)
      iPart = PEM%pNext(iPart)
    END DO
    TreeNode%NodeDepth = 1
    ElemNodeVol(iElem)%Root%NodeDepth = 1
    ElemNodeVol(iElem)%Root%MidPoint(1:3) = (/0.0,0.0,0.0/)
    CALL AddDoTreeNode(TreeNode, iElem, ElemNodeVol(iElem)%Root)
    DEALLOCATE(TreeNode%MappedPartStates)
  ELSE
    CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, Volume)
  END IF
ELSE
  CALL PerformPairingAndCollision(TreeNode%iPartIndx_Node, nPart, iElem, Volume)
END IF

DEALLOCATE(TreeNode%iPartIndx_Node)
DEALLOCATE(TreeNode)

END SUBROUTINE DSMC_pairing_dotree


RECURSIVE SUBROUTINE AddDoTreeNode(TreeNode, iElem, NodeVol)
!===================================================================================================================================
!> Adds additional dotree node/branch and decide if this node has to be divided again or start particle pairing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze      ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars         ,ONLY: tTreeNode, DSMC, tNodeVolume, CollInf
USE MOD_Particle_Vars     ,ONLY: nSpecies, PartSpecies
USE MOD_DSMC_Vars         ,ONLY: ElemNodeVol
USE MOD_part_tools        ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                     :: iElem
TYPE(tTreeNode),INTENT(IN), POINTER     :: TreeNode
CLASS(tNodeVolume),INTENT(INOUT)        :: NodeVol
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iPartIndx, ChildNodeID
INTEGER, ALLOCATABLE          :: iPartIndx_ChildNode(:,:)
REAL, ALLOCATABLE             :: MappedPart_ChildNode(:,:,:)
INTEGER                       :: PartNumChildNode(1:2)
REAL                          :: SpecPartNum(nSpecies,1:2), RealParts(1:2)
LOGICAL                       :: ForceNearestNeigh
!===================================================================================================================================
ForceNearestNeigh = .FALSE.
ALLOCATE(iPartIndx_ChildNode(1:2,TreeNode%PNum_Node))
ALLOCATE(MappedPart_ChildNode(1,TreeNode%PNum_Node,1:2))
PartNumChildNode(:) = 0
IF (ABS(NodeVol%MidPoint(1)) .EQ. 1.0) THEN
  CALL Abort(&
    __STAMP__,&
    'ERROR in DoTree Pairing: Too many branches, machine precision reached')
END IF

!         Numbering of the 2 ChildNodes of the DoTree
!      _________
!     |    |    |
!     | 2  | 1  |   ______x
!     |____|____|

DO iPart=1,TreeNode%PNum_Node
  iPartIndx = TreeNode%iPartIndx_Node(iPart)
  ChildNodeID = DOTANTCUBEID(NodeVol%MidPoint(:),TreeNode%MappedPartStates(:,iPart))
  PartNumChildNode(ChildNodeID) = PartNumChildNode(ChildNodeID) + 1
  iPartIndx_ChildNode(ChildNodeID,PartNumChildNode(ChildNodeID)) = iPartIndx
  MappedPart_ChildNode(1,PartNumChildNode(ChildNodeID),ChildNodeID) = TreeNode%MappedPartStates(1,iPart)
END DO

IF(.NOT.ASSOCIATED(NodeVol%SubNode)) THEN
  CALL AddNodeVolumes1D(TreeNode%NodeDepth, ElemNodeVol(iElem)%Root, iElem)
END IF

DO iLoop = 1, 2
! DoTree can only performed if nPart is greater than the defined value, otherwise nearest neighbour pairing
  ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
  SpecPartNum(:,iLoop) = 0.
  RealParts(iLoop) = 0.
  DO iPart = 1, PartNumChildNode(iLoop)
    RealParts(iLoop) = RealParts(iLoop) + GetParticleWeight(iPartIndx_ChildNode(iLoop,iPart))
    SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart)),iLoop) = &
        SpecPartNum(PartSpecies(iPartIndx_ChildNode(iLoop,iPart)),iLoop) + GetParticleWeight(iPartIndx_ChildNode(iLoop,iPart))
  END DO
  ! DoTree can only performed if nPart is greater than the defined value, otherwise nearest neighbour pairing
  IF((PartNumChildNode(iLoop).GE.DSMC%PartNumOctreeNodeMin).AND.(.NOT.ForceNearestNeigh)) THEN
    ! Additional check if nPart is greater than PartNumOctreeNode (default=80) or the mean free path is less than
    ! the side length of a cube (approximation) with same volume as the actual cell -> DoTree
    ! IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
    !   DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop), RealParts(iLoop), Volume(iLoop))
    ! ELSE
      DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum(:,iLoop),REAL(PartNumChildNode(iLoop)), NodeVol%SubNode(iLoop)%Volume)
    ! END IF
    IF((DSMC%MeanFreePath.LT.(NodeVol%SubNode(iLoop)%Length)).OR.(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode)) THEN
      NULLIFY(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
      ALLOCATE(TreeNode%ChildNode%MappedPartStates(1:2,PartNumChildNode(iLoop)))
      TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
      TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
      TreeNode%ChildNode%MappedPartStates(1,1:PartNumChildNode(iLoop))= &
                     MappedPart_ChildNode(1,1:PartNumChildNode(iLoop),iLoop)
      TreeNode%ChildNode%NodeDepth = TreeNode%NodeDepth + 1
      ! Determination of the sub node number for the correct pointer handover (pointer acts as root for further quadtree division)
      CALL AddDoTreeNode(TreeNode%ChildNode, iElem, NodeVol%SubNode(iLoop))
      DEALLOCATE(TreeNode%ChildNode%MappedPartStates)
      DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
      DEALLOCATE(TreeNode%ChildNode)
    ELSE
      CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
            PartNumChildNode(iLoop), iElem, NodeVol%SubNode(iLoop)%Volume)
    END IF
  ELSE IF (PartNumChildNode(iLoop).GT.1) THEN
    CALL PerformPairingAndCollision(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), &
         PartNumChildNode(iLoop), iElem, NodeVol%SubNode(iLoop)%Volume)
  ELSE IF (CollInf%ProhibitDoubleColl.AND.(PartNumChildNode(iLoop).EQ.1)) THEN
    CollInf%OldCollPartner(iPartIndx_ChildNode(iLoop, 1)) = 0
  END IF
END DO

END SUBROUTINE AddDoTreeNode


SUBROUTINE DSMC_init_octree()
!===================================================================================================================================
! Read-in of octree variables and building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_DSMC_Vars             ,ONLY: DSMC, ElemNodeVol
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Vars         ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem!, NodeDepth
!===================================================================================================================================

! If number of particles is greater than OctreePartNumNode, cell is going to be divided for performance of nearest neighbour
IF(Symmetry%Order.EQ.3) THEN
  DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','80')
ELSE IF(Symmetry%Order.EQ.2) THEN
  DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','40')
ELSE
  DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','20')
END IF
! If number of particles is less than OctreePartNumNodeMin, cell is NOT going to be split even if mean free path is not resolved
! 3D: 50/8; 2D: 28/4 -> ca. 6-7 particles per cell
IF(Symmetry%Order.EQ.3) THEN
  DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','50')
  IF (DSMC%PartNumOctreeNodeMin.LT.20) CALL abort(__STAMP__,'ERROR: Given Particles-OctreePartNumNodeMin is less than 20!')
ELSE IF(Symmetry%Order.EQ.2) THEN
  DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','28')
  IF (DSMC%PartNumOctreeNodeMin.LT.10) CALL abort(__STAMP__,'ERROR: Given Particles-OctreePartNumNodeMin is less than 10!')
ELSE
  DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','14')
  IF (DSMC%PartNumOctreeNodeMin.LT.5) CALL abort(__STAMP__,'ERROR: Given Particles-OctreePartNumNodeMin is less than 5!')
END IF

ALLOCATE(ElemNodeVol(nElems))

!Calculate recursive Volumes
DO iElem = 1, nElems
  ALLOCATE(ElemNodeVol(iElem)%Root)
!  DO NodeDepth = 1, 2
!    IF (Symmetry%Order.EQ.3) THEN
!      CALL DSMC_CalcSubNodeVolumes3D(iElem, NodeDepth, ElemNodeVol(iElem)%Root)
!    ELSE IF(Symmetry%Order.EQ.2) THEN
!      CALL DSMC_CalcSubNodeVolumes2D(iElem, NodeDepth, ElemNodeVol(iElem)%Root)
!    ELSE
!      CALL AddNodeVolumes1D(NodeDepth, ElemNodeVol(iElem)%Root, iElem)
!    END IF
!  END DO
END DO

END SUBROUTINE DSMC_init_octree

SUBROUTINE DSMC_CalcSubNodeVolumes2D(iElem, NodeDepth, Node)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : OctreeVdm, tNodeVolume, SymmetrySide
USE MOD_Mesh_Vars,              ONLY : SurfElem, Face_xGP
USE MOD_Preproc
USE MOD_ChangeBasis,            ONLY : ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                       :: iElem
INTEGER, INTENT(IN)                       :: NodeDepth
CLASS(tNodeVolume), INTENT(INOUT)         :: Node
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


SUBROUTINE DSMC_CalcSubNodeVolumes3D(iElem, NodeDepth, Node)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : OctreeVdm, tNodeVolume
  USE MOD_Mesh_Vars,              ONLY : sJ
  USE MOD_Preproc
  USE MOD_ChangeBasis,            ONLY : ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                       :: iElem
  INTEGER, INTENT(IN)                       :: NodeDepth
  CLASS(tNodeVolume), INTENT(INOUT)         :: Node
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

END SUBROUTINE DSMC_CalcSubNodeVolumes3D

RECURSIVE SUBROUTINE InitVanderOct(LocalVdm, NodeDepth, LocalDepth, OctreeVdmLoc)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : tOctreeVdm
  USE MOD_Preproc
  USE MOD_Interpolation_Vars,     ONLY : xGP, wBary
  USE MOD_Basis,                  ONLY : InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                       :: NodeDepth
  INTEGER, INTENT(INOUT)                    :: LocalDepth
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
  INTEGER, INTENT(IN)                       :: NodeDepth
  INTEGER, INTENT(IN), OPTIONAL             :: SubNodesIn(:)
  CLASS(tNodeVolume), INTENT(INOUT)         :: Node
  REAL, INTENT(INOUT)                       :: DetJac(1,0:2**NodeDepth-1,0:2**NodeDepth-1,0:2**NodeDepth-1)
  TYPE (tOctreeVdm), INTENT(OUT), POINTER   :: VdmLocal
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER, ALLOCATABLE                       :: SubNodesOut(:)
  INTEGER                                    :: OldNodeNum, NewNodeNum, j, VolPos(3), iNode
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
    IF(.NOT.ASSOCIATED(Node%SubNode)) THEN
      ALLOCATE(Node%SubNode(8))
      DO iNode = 1, 8
        Node%SubNode(iNode)%NodeDepth = Node%NodeDepth + 1
        Node%SubNode(iNode)%MidPoint(1:3) = OCTANTCUBEMIDPOINT(iNode,Node%NodeDepth,Node%MidPoint(1:3))
      END DO
    END IF
    DO iNode = 1, 8
      SubNodesOut(NewNodeNum) = iNode
      CALL AddNodeVolumes(NodeDepth, Node%SubNode(iNode), DetJac, VdmLocal, SubNodesOut)
    END DO
  ELSE
    IF(.NOT.ASSOCIATED(Node%SubNode)) THEN
      ALLOCATE(Node%SubNode(8))
      DO iNode = 1, 8
        Node%SubNode(iNode)%NodeDepth = Node%NodeDepth + 1
        Node%SubNode(iNode)%MidPoint(1:3) = OCTANTCUBEMIDPOINT(iNode,Node%NodeDepth,Node%MidPoint(1:3))
      END DO
    END IF
    DO iNode = 1, 8
      SubNodesOut(NewNodeNum) = iNode
      CALL AddNodeVolumes(NodeDepth, Node%SubNode(iNode), DetJac, VdmLocal%SubVdm, SubNodesOut)
    END DO
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
USE MOD_Particle_Vars   ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                       :: NodeDepth
INTEGER, INTENT(IN)                       :: iElem
INTEGER, INTENT(IN), OPTIONAL             :: SubNodesIn(:)
CLASS(tNodeVolume), INTENT(INOUT)         :: Node
REAL, INTENT(INOUT)                       :: DetJac(1,0:2**NodeDepth-1,0:2**NodeDepth-1)
REAL, INTENT(INOUT)                       :: FacexGP(2,0:2**NodeDepth-1,0:2**NodeDepth-1)
TYPE (tOctreeVdm), INTENT(OUT), POINTER   :: VdmLocal
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE                      :: SubNodesOut(:)
INTEGER                                   :: OldNodeNum, NewNodeNum, j, VolPos(2), dirX, dirY, realDirX, realDirY, iNode
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
    IF(.NOT.ASSOCIATED(Node%SubNode)) THEN
      ALLOCATE(Node%SubNode(4))
      DO iNode = 1, 4
        Node%SubNode(iNode)%NodeDepth = Node%NodeDepth + 1
        Node%SubNode(iNode)%MidPoint(1:3) = QUADCUBEMIDPOINT(iNode,Node%NodeDepth,Node%MidPoint(1:3))
      END DO
    END IF
    DO iNode = 1, 4
      SubNodesOut(NewNodeNum) = iNode
      CALL AddNodeVolumes2D(NodeDepth, Node%SubNode(iNode), DetJac, VdmLocal, iElem, FacexGP, SubNodesOut)
    END DO
  ELSE
    IF(.NOT.ASSOCIATED(Node%SubNode)) THEN
      ALLOCATE(Node%SubNode(4))
      DO iNode = 1, 4
        Node%SubNode(iNode)%NodeDepth = Node%NodeDepth + 1
        Node%SubNode(iNode)%MidPoint(1:3) = QUADCUBEMIDPOINT(iNode,Node%NodeDepth,Node%MidPoint(1:3))
      END DO
    END IF
    DO iNode = 1, 4
      SubNodesOut(NewNodeNum) = iNode
      CALL AddNodeVolumes2D(NodeDepth, Node%SubNode(iNode), DetJac, VdmLocal%SubVdm, iElem, FacexGP, SubNodesOut)
    END DO
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
  IF (Symmetry%Axisymmetric) THEN
    Node%Volume = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2 * 2. * Pi * Pos(2)
  ELSE
    Node%Volume = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2
  END IF
  Node%Area = DetJac(1, realDirX, realDirY) * VdmLocal%wGP**2
END IF

END SUBROUTINE AddNodeVolumes2D


SUBROUTINE FindRandomPartner(iPartIndx_Node, nPart, iPair, nPair)
!===================================================================================================================================
! Classic statistical pairing method for the use in the octree routines
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: Coll_pData,CollInf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)        :: nPart
INTEGER, INTENT(IN)           :: iPair, nPair
INTEGER, INTENT(INOUT)        :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: cPart1, cPart2
REAL                          :: iRan
!===================================================================================================================================

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

END SUBROUTINE FindRandomPartner


PPURE INTEGER FUNCTION OCTANTCUBEID(centerPoint,coord)
!===================================================================================================================================
!> determine position of Coord in a cube Octant in reference to a 3D centerpoint
!>         Numbering of the 8 Octant IDs (octree)
!>          __________
!>         /    /    /|   |z
!>        /  7 / 6  / |   |
!>       /----/----/|6|   |
!>      /_3__/_2__/ |/|   |______Y
!>     |    |    |2 |5|   /
!>     | 3  | 2  |/ |/   /
!>     |----|----|1 /   /x
!>     | 4  | 1  | /
!>     |____|____|/
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN) :: centerPoint(1:3)
REAL,INTENT(IN) :: coord(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
IF ((coord(1).GE.centerPoint(1)) .AND. (coord(2).GE.centerPoint(2)) .AND.(coord(3).LE.centerPoint(3))) THEN
  OCTANTCUBEID=1
ELSE IF((coord(1).GE.centerPoint(1)) .AND.(coord(2).GE.centerPoint(2))) THEN
  OCTANTCUBEID=2
ELSE IF((coord(1).GE.centerPoint(1)) .AND.(coord(3).GE.centerPoint(3))) THEN
  OCTANTCUBEID=3
ELSE IF (coord(1).GE.centerPoint(1)) THEN
  OCTANTCUBEID=4
ELSE IF((coord(2).GE.centerPoint(2)) .AND.(coord(3).LE.centerPoint(3))) THEN
  OCTANTCUBEID=5
ELSE IF (coord(2).GE.centerPoint(2)) THEN
  OCTANTCUBEID=6
ELSE IF (coord(3).GE.centerPoint(3)) THEN
  OCTANTCUBEID=7
ELSE
  OCTANTCUBEID=8
END IF

END FUNCTION OCTANTCUBEID


PPURE FUNCTION OCTANTCUBEMIDPOINT(CubeID,octantDepth,octantCenter)
!===================================================================================================================================
!> determines the position of the center of the given cubeID of an Octant for a given 3D centerpoint and depth
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: CubeID
INTEGER,INTENT(IN) :: octantDepth
REAL,INTENT(IN)    :: octantCenter(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, DIMENSION(3) :: OCTANTCUBEMIDPOINT
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
IF (CubeID.LT.5) THEN
  OCTANTCUBEMIDPOINT(1) = 1.0
  IF (CubeID.LT.3) THEN
    OCTANTCUBEMIDPOINT(2) = 1.0
  ELSE
    OCTANTCUBEMIDPOINT(2) = -1.0
  END IF
ELSE
  OCTANTCUBEMIDPOINT(1) = -1.0
  IF (CubeID.LT.7) THEN
    OCTANTCUBEMIDPOINT(2) = 1.0
  ELSE
    OCTANTCUBEMIDPOINT(2) = -1.0
  END IF
END IF
IF ((CubeID.EQ.1).OR.(CubeID.EQ.4).OR.(CubeID.EQ.5).OR.(CubeID.EQ.8)) THEN
  OCTANTCUBEMIDPOINT(3) = -1.0
ELSE
  OCTANTCUBEMIDPOINT(3) = 1.0
END IF
OCTANTCUBEMIDPOINT(1:3) = octantCenter(1:3) + OCTANTCUBEMIDPOINT(1:3)*2.0/(2.0**(REAL(octantDepth)+1.))
END FUNCTION OCTANTCUBEMIDPOINT


PPURE FUNCTION QUADCUBEMIDPOINT(CubeID,octantDepth,octantCenter)
!===================================================================================================================================
!> determines the position of the center of the given cubeID of an Octant for a given 3D centerpoint and depth
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: CubeID
INTEGER,INTENT(IN) :: octantDepth
REAL,INTENT(IN)    :: octantCenter(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, DIMENSION(3) :: QUADCUBEMIDPOINT
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
IF (CubeID.LT.3) THEN
  QUADCUBEMIDPOINT(1) = 1.0
  IF (CubeID.EQ.1) THEN
    QUADCUBEMIDPOINT(2) = -1.0
  ELSE
    QUADCUBEMIDPOINT(2) = 1.0
  END IF
ELSE
  QUADCUBEMIDPOINT(1) = -1.0
  IF (CubeID.EQ.3) THEN
    QUADCUBEMIDPOINT(2) = 1.0
  ELSE
    QUADCUBEMIDPOINT(2) = -1.0
  END IF
END IF
QUADCUBEMIDPOINT(1:2) = octantCenter(1:2) + QUADCUBEMIDPOINT(1:2)*2.0/(2.0**(REAL(octantDepth)+1.0))
QUADCUBEMIDPOINT(3) = 0.0
END FUNCTION QUADCUBEMIDPOINT


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
USE MOD_DSMC_Vars             ,ONLY: SymmetrySide
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
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
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,SideID,GetCNElemID(iElem+offSetElem))+1)
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
USE MOD_DSMC_Vars               ,ONLY: SymmetrySide
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared
USE MOD_Mesh_Vars               ,ONLY: offsetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
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
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,SideID,GetCNElemID(iElem+offSetElem))+1)
END DO

MapToGeo2D =0.25*(P(:,1)*(1-xi(1)) * (1-xi(2)) &
              + P(:,2)*(1+xi(1)) * (1-xi(2)) &
              + P(:,3)*(1+xi(1)) * (1+xi(2)) &
              + P(:,4)*(1-xi(1)) * (1+xi(2)))

END FUNCTION MapToGeo2D


INTEGER FUNCTION DOTANTCUBEID(centerPoint,coord)
!===================================================================================================================================
!> determine position of Coord in a cube Dotant in reference to a 3D centerpoint
!>         Numbering of the 2 Dotant IDs (dotree)
!>      _________
!>     |    |    |
!>     | 2  | 1  |   ______x
!>     |____|____|

!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN) :: centerPoint(1:3)
REAL,INTENT(IN) :: coord(1)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================

IF(coord(1).GE.centerpoint(1)) THEN
  DOTANTCUBEID = 1
ELSE
  DOTANTCUBEID = 2
END IF

END FUNCTION DOTANTCUBEID

FUNCTION DOTANTCUBEMIDPOINT(CubeID,octantDepth,octantCenter)
!===================================================================================================================================
!> determines the position of the center of the given cubeID of an Dotant for a given 3D centerpoint and depth
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: CubeID
INTEGER,INTENT(IN) :: octantDepth
REAL,INTENT(IN)    :: octantCenter(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, DIMENSION(3) :: DOTANTCUBEMIDPOINT
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
IF (CubeID.EQ.1) THEN
  DOTANTCUBEMIDPOINT(1) = 1.0
ELSE
  DOTANTCUBEMIDPOINT(1) = -1.0
END IF
DOTANTCUBEMIDPOINT(2:3) = 0.0
DOTANTCUBEMIDPOINT(1) = octantCenter(1) + DOTANTCUBEMIDPOINT(1)*2.0/REAL(2.0**(octantDepth+1))
END FUNCTION DOTANTCUBEMIDPOINT


SUBROUTINE GeoCoordToMap1D(x_in,xi_Out,iElem)
!===================================================================================================================================
!> get the x-Coordinate in the Ref Element by linear interpolation
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem                                 !< elem index
REAL,INTENT(IN)               :: x_in                                  !< physical position of particle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: xi_Out                                ! Position in Ref Space
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

xi_Out = -1. + (2.)/(GEO%XMinMax(2,iElem)-GEO%XMinMax(1,iElem)) * (x_in-GEO%XMinMax(1,iElem))

END SUBROUTINE GeoCoordToMap1D


RECURSIVE SUBROUTINE AddNodeVolumes1D(NodeDepth, Node, iElem, SubNodesIn)
!===================================================================================================================================
!> calculate volume and x-length of an DotreeNode
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: Pi
USE MOD_DSMC_Vars           ,ONLY: tOctreeVdm, tNodeVolume
USE MOD_Particle_Vars       ,ONLY: Symmetry
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                       :: NodeDepth
INTEGER, INTENT(IN)                       :: iElem
INTEGER, INTENT(IN), OPTIONAL             :: SubNodesIn(:)
CLASS(tNodeVolume), INTENT(INOUT)         :: Node
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE                      :: SubNodesOut(:)
INTEGER                                   :: OldNodeNum, NewNodeNum, iDepth, iNode
REAL                                      :: Xmin, Xmax
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
  IF(.NOT.ASSOCIATED(Node%SubNode)) THEN
    ALLOCATE(Node%SubNode(2))
    DO iNode = 1, 2
      Node%SubNode(iNode)%NodeDepth = Node%NodeDepth + 1
      Node%SubNode(iNode)%MidPoint(1:3) = DOTANTCUBEMIDPOINT(iNode,Node%NodeDepth,Node%MidPoint(1:3))
    END DO
  END IF
  DO iNode = 1, 2
    SubNodesOut(NewNodeNum) = iNode
    CALL AddNodeVolumes1D(NodeDepth, Node%SubNode(iNode), iElem, SubNodesOut)
  END DO
ELSE
  IF(Symmetry%Axisymmetric) THEN
    Xmin=GEO%XMinMax(1,iElem)
    Xmax=GEO%XMinMax(2,iElem)
    DO iDepth=1,NodeDepth
      SELECT CASE(SubNodesIn(NodeDepth))
      CASE(1)
        Xmax = 0.5*(Xmax-Xmin) + Xmin
      CASE(2)
        XMin = 0.5*(Xmax-Xmin) + Xmin
      CASE DEFAULT
        CALL ABORT(__STAMP__&
                  ,'ERROR: Wrong SubNode Index',IntInfoOpt=SubNodesIn(NodeDepth))
      END SELECT
    END DO
    Node%Length = Xmax-Xmin
    IF(Symmetry%Axisymmetric) THEN
      Node%Volume = (Xmax-Xmin) * Pi * (Xmax+Xmin)
    ELSE
      Node%Volume = 4./3. * PI * ( Xmax**3 - Xmin**3 )
    END IF
  ELSE
    Node%Volume = (GEO%XMinMax(2,iElem) - GEO%XMinMax(1,iElem)) / 2.**REAL(NodeDepth)
    Node%Length = Node%Volume
  END IF
END IF

END SUBROUTINE AddNodeVolumes1D

END MODULE MOD_DSMC_ParticlePairing
