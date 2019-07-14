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

MODULE MOD_BGK
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

SUBROUTINE BGK_DSMC_main()
!===================================================================================================================================
!> Coupled BGK and DSMC routine: Cell-local decision with BGKDSMCSwitchDens
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_DSMC_Vars           ,ONLY: DSMC_RHS, DSMC
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_Particle_Vars       ,ONLY: PEM, PartState, PartSpecies, Species, WriteMacroVolumeValues
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation,BGKMovingAverage,ElemNodeAveraging,BGKMovingAverageLength,BGKDSMCSwitchDens
USE MOD_BGK_Vars            ,ONLY: BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor,BGK_QualityFacSamp
USE MOD_BGK_Vars            ,ONLY: BGK_MaxRotRelaxFactor
USE MOD_BGK_CollOperator    ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_Analyze        ,ONLY: DSMCHO_data_sampling
USE MOD_DSMC                ,ONLY: DSMC_main
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart, iSpec
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
LOGICAL               :: DoElement(nElems)
REAL                  :: vBulk(3), TotalMass, dens
!===================================================================================================================================
DSMC_RHS = 0.0
DoElement = .FALSE.

DO iElem = 1, nElems
  nPart = PEM%pNumber(iElem)
  IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE
  dens = nPart * Species(1)%MacroParticleFactor / GEO%Volume(iElem) 
  IF (dens.LT.BGKDSMCSwitchDens) THEN
    DoElement(iElem) = .TRUE.
    CYCLE
  END IF

  IF (DoBGKCellAdaptation) THEN
    CALL BGK_octree_adapt(iElem)
  ELSE  
    ALLOCATE(iPartIndx_Node(nPart))
    TotalMass = 0.0
    vBulk(1:3) = 0.0
    iPart = PEM%pStart(iElem)
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      iSpec = PartSpecies(iPart)
      vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)*Species(iSpec)%MassIC
      TotalMass = TotalMass + Species(iSpec)%MassIC
      iPart = PEM%pNext(iPart)
    END DO
    vBulk = vBulk / TotalMass

    IF(DSMC%CalcQualityFactors) THEN
      BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    END IF
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk, &
          ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKMovingAverageLength), &
               CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
    ELSE 
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk)
    END IF
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
        BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
        BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
        BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
        BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
        BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
      END IF
    END IF
    DEALLOCATE(iPartIndx_Node)
  END IF
END DO

CALL DSMC_main(DoElement)

END SUBROUTINE BGK_DSMC_main


SUBROUTINE BGK_main()
!===================================================================================================================================
!> Main routine for the BGK model
!> 1.) Loop over all elements, call of octree refinement or directly of the collision operator
!> 2.) Sampling of macroscopic variables with DSMC routines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: TEnd, Time
USE MOD_Mesh_Vars           ,ONLY: nElems, MeshFile
USE MOD_DSMC_Vars           ,ONLY: DSMC_RHS, DSMC, SamplingActive
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_Particle_Vars       ,ONLY: PEM, PartState, WriteMacroVolumeValues, WriteMacroSurfaceValues, Species, PartSpecies
USE MOD_Restart_Vars        ,ONLY: RestartTime
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation, BGKMovingAverage, ElemNodeAveraging, BGKMovingAverageLength
USE MOD_BGK_Vars            ,ONLY: BGK_MeanRelaxFactor,BGK_MeanRelaxFactorCounter,BGK_MaxRelaxFactor,BGK_QualityFacSamp
USE MOD_BGK_Vars            ,ONLY: BGK_MaxRotRelaxFactor
USE MOD_BGK_CollOperator    ,ONLY: BGK_CollisionOperator
USE MOD_DSMC_Analyze        ,ONLY: DSMCHO_data_sampling,CalcSurfaceValues,WriteDSMCHOToHDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart, nOutput, iSpec
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
REAL                  :: vBulk(3), TotalMass
!===================================================================================================================================
DSMC_RHS = 0.0

IF (DoBGKCellAdaptation) THEN
  DO iElem = 1, nElems
    CALL BGK_octree_adapt(iElem)
  END DO
ELSE ! No octree cell refinement
  DO iElem = 1, nElems
    nPart = PEM%pNumber(iElem)
    IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE
    ALLOCATE(iPartIndx_Node(nPart))
    vBulk(1:3) = 0.0
    TotalMass = 0.0
    iPart = PEM%pStart(iElem)
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      iSpec = PartSpecies(iPart)
      vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)*Species(iSpec)%MassIC
      TotalMass = TotalMass + Species(iSpec)%MassIC
      iPart = PEM%pNext(iPart)
    END DO
    vBulk = vBulk / TotalMass

    IF(DSMC%CalcQualityFactors) THEN
      BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    END IF

    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk, &
          ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKMovingAverageLength), &
               CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
    ELSE
      CALL BGK_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk)
    END IF
    IF(DSMC%CalcQualityFactors) THEN
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
        BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
        BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
        BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
        BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
        BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
      END IF
    END IF
    DEALLOCATE(iPartIndx_Node)
  END DO
END IF ! DoBGKCellAdaptation

IF((.NOT.WriteMacroVolumeValues) .AND. (.NOT.WriteMacroSurfaceValues)) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).AND.(.NOT.SamplingActive))  THEN
    SamplingActive=.TRUE.
    SWRITE(*,*)'Sampling active'
  END IF
END IF

IF(SamplingActive) THEN
  CALL DSMCHO_data_sampling()
  IF(DSMC%NumOutput.NE.0) THEN
    nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
    IF(Time.GE.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
      DSMC%NumOutput = DSMC%NumOutput - 1
      ! Skipping outputs immediately after the first few iterations
      IF(RestartTime.LT.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * REAL(nOutput))) THEN 
        CALL WriteDSMCHOToHDF5(TRIM(MeshFile),time)
        IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues(during_dt_opt=.TRUE.)
      END IF
    END IF
  END IF
END IF

END SUBROUTINE BGK_main

END MODULE MOD_BGK
