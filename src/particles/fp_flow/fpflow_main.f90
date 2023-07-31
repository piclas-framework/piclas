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

MODULE MOD_FPFlow
!===================================================================================================================================
! Module for the main Fokker-Planck routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE FPFlow_main
  MODULE PROCEDURE FPFlow_main
END INTERFACE

INTERFACE FP_DSMC_main
  MODULE PROCEDURE FP_DSMC_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: FPFlow_main, FP_DSMC_main
!===================================================================================================================================

CONTAINS

SUBROUTINE FP_DSMC_main()
!===================================================================================================================================
!> Coupled FP and DSMC routine: Cell-local decision with FPDSMCSwitchDens
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Mesh_Vars           ,ONLY: nElems, offsetElem
USE MOD_Particle_Vars       ,ONLY: PEM, Species, WriteMacroVolumeValues, Symmetry, usevMPF
USE MOD_FP_CollOperator     ,ONLY: FP_CollisionOperator
USE MOD_FPFlow_Vars         ,ONLY: FPDSMCSwitchDens, FP_QualityFacSamp, FP_PrandtlNumber
USE MOD_FPFlow_Vars         ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_DSMC_Vars           ,ONLY: DSMC, RadialWeighting
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt, BGK_quadtree_adapt
USE MOD_DSMC                ,ONLY: DSMC_main
USE MOD_Part_Tools          ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart, CNElemID
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
LOGICAL               :: DoElement(nElems)
REAL                  :: dens, partWeight, totalWeight
!===================================================================================================================================
DoElement = .FALSE.

DO iElem = 1, nElems
  CNElemID = GetCNElemID(iElem + offsetElem)
  nPart = PEM%pNumber(iElem)
  totalWeight = 0.0
  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    partWeight = GetParticleWeight(iPart)
    totalWeight = totalWeight + partWeight
    iPart = PEM%pNext(iPart)
  END DO

  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
    dens = totalWeight / ElemVolume_Shared(CNElemID)
  ELSE
    dens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
  END IF
  IF (dens.LT.FPDSMCSwitchDens) THEN
    DoElement(iElem) = .TRUE.
    CYCLE
  END IF
  IF (nPart.LT.3) CYCLE

  IF (DoBGKCellAdaptation) THEN
    IF(Symmetry%Order.EQ.2) THEN
      CALL BGK_quadtree_adapt(iElem)
    ELSE
      CALL BGK_octree_adapt(iElem)
    END IF
  ELSE
    ALLOCATE(iPartIndx_Node(nPart))
    iPart = PEM%pStart(iElem)
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      iPart = PEM%pNext(iPart)
    END DO

    IF(DSMC%CalcQualityFactors) THEN
      FP_MeanRelaxFactorCounter=0; FP_MeanRelaxFactor=0.; FP_MaxRelaxFactor=0.; FP_MaxRotRelaxFactor=0.; FP_PrandtlNumber=0.
    END IF

    CALL FP_CollisionOperator(iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
    DEALLOCATE(iPartIndx_Node)
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
        FP_QualityFacSamp(1,iElem) = FP_QualityFacSamp(1,iElem) + FP_MeanRelaxFactor
        FP_QualityFacSamp(2,iElem) = FP_QualityFacSamp(2,iElem) + REAL(FP_MeanRelaxFactorCounter)
        FP_QualityFacSamp(3,iElem) = FP_QualityFacSamp(3,iElem) + FP_MaxRelaxFactor
        FP_QualityFacSamp(4,iElem) = FP_QualityFacSamp(4,iElem) + 1.
        FP_QualityFacSamp(5,iElem) = FP_QualityFacSamp(5,iElem) + FP_MaxRotRelaxFactor
        FP_QualityFacSamp(6,iElem) = FP_QualityFacSamp(6,iElem) + FP_PrandtlNumber
      END IF
    END IF
  END IF
END DO

CALL DSMC_main(DoElement)

END SUBROUTINE FP_DSMC_main


SUBROUTINE FPFlow_main()
!===================================================================================================================================
!> Main routine for the FP model
!> 1.) Loop over all elements, call of octree refinement or directly of the collision operator
!> 2.) Sampling of macroscopic variables with DSMC routines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Mesh_Vars           ,ONLY: nElems, offsetElem
USE MOD_Particle_Vars       ,ONLY: PEM, WriteMacroVolumeValues, WriteMacroSurfaceValues, Symmetry, DoVirtualCellMerge, VirtMergedCells
USE MOD_FP_CollOperator     ,ONLY: FP_CollisionOperator
USE MOD_DSMC_Vars           ,ONLY: DSMC
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt, BGK_quadtree_adapt
USE MOD_FPFlow_Vars         ,ONLY: FP_QualityFacSamp, FP_PrandtlNumber
USE MOD_FPFlow_Vars         ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
USE MOD_DSMC_Analyze        ,ONLY: DSMCMacroSampling
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iElem, nPart, iPart, iLoop, CNElemID, nPartMerged, iMergeElem, iLoopLoc, locElem, nPartLoc
INTEGER, ALLOCATABLE        :: iPartIndx_Node(:)
REAL                        :: elemVolume
!===================================================================================================================================

IF (DoBGKCellAdaptation) THEN
  DO iElem = 1, nElems
    IF(Symmetry%Order.EQ.2) THEN
      CALL BGK_quadtree_adapt(iElem)
    ELSE
      CALL BGK_octree_adapt(iElem)
    END IF
  END DO
ELSE
  DO iElem = 1, nElems
    CNElemID = GetCNElemID(iElem + offsetElem)
    nPart = PEM%pNumber(iElem)
    IF (DoVirtualCellMerge) THEN
      IF(VirtMergedCells(iElem)%isMerged) CYCLE      
      nPartMerged = nPart
      DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
        nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
      END DO
      IF ((nPartMerged.EQ.0).OR.(nPartMerged.EQ.1)) CYCLE
      ALLOCATE(iPartIndx_Node(nPartMerged))
      iPart = PEM%pStart(iElem)
      iLoopLoc = 0
      DO iLoop = 1, nPart
        iLoopLoc = iLoopLoc + 1
        iPartIndx_Node(iLoopLoc) = iPart
        iPart = PEM%pNext(iPart)
      END DO
      IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
        DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
          locElem = VirtMergedCells(iElem)%MergedCellID(iMergeElem)
          nPartLoc = PEM%pNumber(locElem)
          iPart = PEM%pStart(locElem)
          DO iLoop = 1, nPartLoc
            iLoopLoc = iLoopLoc + 1
            iPartIndx_Node(iLoopLoc) = iPart
            iPart = PEM%pNext(iPart)
          END DO
        END DO
        elemVolume = VirtMergedCells(iELem)%MergedVolume
      ELSE
        elemVolume = ElemVolume_Shared(CNElemID)
      END IF        
    ELSE      
      nPartMerged = nPart   
      IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE
      ALLOCATE(iPartIndx_Node(nPart))
      iPart = PEM%pStart(iElem)
      DO iLoop = 1, nPart
        iPartIndx_Node(iLoop) = iPart
        iPart = PEM%pNext(iPart)
      END DO
      elemVolume = ElemVolume_Shared(CNElemID)
    END IF

    IF(DSMC%CalcQualityFactors) THEN
      FP_MeanRelaxFactorCounter=0; FP_MeanRelaxFactor=0.; FP_MaxRelaxFactor=0.; FP_MaxRotRelaxFactor=0.; FP_PrandtlNumber=0.
    END IF

    CALL FP_CollisionOperator(iPartIndx_Node, nPartMerged, elemVolume)
    DEALLOCATE(iPartIndx_Node)
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
        FP_QualityFacSamp(1,iElem) = FP_QualityFacSamp(1,iElem) + FP_MeanRelaxFactor
        FP_QualityFacSamp(2,iElem) = FP_QualityFacSamp(2,iElem) + REAL(FP_MeanRelaxFactorCounter)
        FP_QualityFacSamp(3,iElem) = FP_QualityFacSamp(3,iElem) + FP_MaxRelaxFactor
        FP_QualityFacSamp(4,iElem) = FP_QualityFacSamp(4,iElem) + 1.
        FP_QualityFacSamp(5,iElem) = FP_QualityFacSamp(5,iElem) + FP_MaxRotRelaxFactor
        FP_QualityFacSamp(6,iElem) = FP_QualityFacSamp(6,iElem) + FP_PrandtlNumber
      END IF
    END IF
  END DO
END IF

! Sampling of macroscopic values
! (here for a continuous average; average over N iterations is performed in src/analyze/analyze.f90)
IF (.NOT.WriteMacroVolumeValues .AND. .NOT.WriteMacroSurfaceValues) THEN
  CALL DSMCMacroSampling()
END IF

END SUBROUTINE FPFlow_main

END MODULE MOD_FPFLOW
