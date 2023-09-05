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

MODULE MOD_BGK
#if (PP_TimeDiscMethod==400)
!===================================================================================================================================
!> Main module for the the Bhatnagar-Gross-Krook method
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGK_main
  MODULE PROCEDURE BGK_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGK_main, BGK_DSMC_main
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_DSMC_main(stage_opt)
!===================================================================================================================================
!> Coupled BGK and DSMC routine: Cell-local decision with BGKDSMCSwitchDens
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt, BGK_quadtree_adapt
USE MOD_Particle_Vars       ,ONLY: PEM, Species, WriteMacroVolumeValues, Symmetry, usevMPF
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation,BGKDSMCSwitchDens
USE MOD_BGK_Vars            ,ONLY: BGKMovingAverage,ElemNodeAveraging
USE MOD_BGK_Vars            ,ONLY: BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor,BGK_QualityFacSamp
USE MOD_BGK_Vars            ,ONLY: BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_BGK_Vars            ,ONLY: BGK_Viscosity, BGK_ThermalConductivity
USE MOD_BGK_CollOperator    ,ONLY: BGK_CollisionOperator
USE MOD_DSMC                ,ONLY: DSMC_main
USE MOD_DSMC_Vars           ,ONLY: DSMC, RadialWeighting
USE MOD_Mesh_Vars           ,ONLY: nElems, offsetElem
USE MOD_Part_Tools          ,ONLY: GetParticleWeight
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_Particle_Mesh_Vars  ,ONLY: IsExchangeElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: stage_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart, CNElemID, stage
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
LOGICAL               :: DoElement(nElems)
REAL                  :: dens, partWeight, totalWeight
!===================================================================================================================================
IF (PRESENT(stage_opt)) THEN
  stage = stage_opt
ELSE
  stage = 0
END IF

DoElement = .FALSE.

DO iElem = 1, nElems
#if USE_MPI
  IF (stage.EQ.1) THEN
    IF (IsExchangeElem(iElem)) CYCLE
  ELSE IF (stage.EQ.2) THEN
    IF (.NOT.IsExchangeElem(iELem)) CYCLE
  END IF
#endif
  nPart = PEM%pNumber(iElem)
  CNElemID = GetCNElemID(iElem + offsetElem)
  IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE

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

  IF (dens.LT.BGKDSMCSwitchDens) THEN
    DoElement(iElem) = .TRUE.
    CYCLE
  END IF

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
      BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
      BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.; BGK_Viscosity=0.; BGK_ThermalConductivity=0.
    END IF
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID), ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, ElemVolume_Shared(CNElemID))
    END IF
    DEALLOCATE(iPartIndx_Node)
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
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
    END IF
  END IF
END DO

CALL DSMC_main(DoElement)

END SUBROUTINE BGK_DSMC_main


SUBROUTINE BGK_main(stage_opt)
!===================================================================================================================================
!> Main routine for the BGK model
!> 1.) Loop over all elements, call of octree refinement or directly of the collision operator
!> 2.) Sampling of macroscopic variables with DSMC routines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Mesh_Vars           ,ONLY: nElems, offsetElem
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt, BGK_quadtree_adapt
USE MOD_Particle_Vars       ,ONLY: PEM, WriteMacroVolumeValues, WriteMacroSurfaceValues, Symmetry, DoVirtualCellMerge, VirtMergedCells
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation, BGKMovingAverage, ElemNodeAveraging
USE MOD_BGK_Vars            ,ONLY: BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor,BGK_QualityFacSamp
USE MOD_BGK_Vars            ,ONLY: BGK_MaxRotRelaxFactor, BGK_PrandtlNumber, BGK_ExpectedPrandtlNumber
USE MOD_BGK_Vars            ,ONLY: BGK_Viscosity, BGK_ThermalConductivity
USE MOD_BGK_CollOperator    ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_Analyze        ,ONLY: DSMCMacroSampling
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_DSMC_Vars           ,ONLY: DSMC
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_Particle_Mesh_Vars  ,ONLY: IsExchangeElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: stage_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iElem, nPart, iLoop, iPart, CNElemID, stage, nPartMerged, iMergeElem, iLoopLoc, locElem, nPartLoc
INTEGER, ALLOCATABLE        :: iPartIndx_Node(:)
REAL                        :: elemVolume
!===================================================================================================================================
IF (PRESENT(stage_opt)) THEN
  stage = stage_opt
ELSE
  stage = 0
END IF

IF (DoBGKCellAdaptation) THEN
  DO iElem = 1, nElems
#if USE_MPI
    IF (stage.EQ.1) THEN
      IF (IsExchangeElem(iElem)) CYCLE
    ELSE IF (stage.EQ.2) THEN
      IF (.NOT.IsExchangeElem(iELem)) CYCLE
    END IF
#endif
    IF(Symmetry%Order.EQ.2) THEN
      CALL BGK_quadtree_adapt(iElem)
    ELSE
      CALL BGK_octree_adapt(iElem)
    END IF
  END DO
ELSE ! No octree cell refinement
  DO iElem = 1, nElems
#if USE_MPI
    IF (stage.EQ.1) THEN
      IF (IsExchangeElem(iElem)) CYCLE
    ELSE IF (stage.EQ.2) THEN
      IF (.NOT.IsExchangeElem(iELem)) CYCLE
    END IF
#endif
    CNElemID = GetCNElemID(iElem + offsetElem)
    nPart = PEM%pNumber(iElem)
    IF (DoVirtualCellMerge) THEN
      IF(VirtMergedCells(iElem)%isMerged) CYCLE      
      nPartMerged = nPart
      DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
        nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
      END DO
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
      BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
      BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.; BGK_Viscosity=0.; BGK_ThermalConductivity=0.
    END IF

    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(iPartIndx_Node, nPartMerged, elemVolume,ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(iPartIndx_Node, nPartMerged, elemVolume)
    END IF
    DEALLOCATE(iPartIndx_Node)
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
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
    END IF
  END DO
END IF ! DoBGKCellAdaptation

! Sampling of macroscopic values
! (here for a continuous average; average over N iterations is performed in src/analyze/analyze.f90)
IF (.NOT.WriteMacroVolumeValues .AND. .NOT.WriteMacroSurfaceValues) THEN
  IF ((stage.EQ.0).OR.(stage.EQ.2)) CALL DSMCMacroSampling()
END IF

END SUBROUTINE BGK_main

#endif /*(PP_TimeDiscMethod==400)*/
END MODULE MOD_BGK
