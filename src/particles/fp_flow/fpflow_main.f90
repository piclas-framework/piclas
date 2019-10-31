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
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Particle_Vars       ,ONLY: PEM, PartState, Species, WriteMacroVolumeValues
USE MOD_FP_CollOperator     ,ONLY: FP_CollisionOperator
USE MOD_FPFlow_Vars         ,ONLY: FPDSMCSwitchDens, FP_QualityFacSamp, FP_PrandtlNumber
USE MOD_FPFlow_Vars         ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_DSMC_Vars           ,ONLY: DSMC_RHS, DSMC
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt
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
INTEGER               :: iElem, nPart, iLoop, iPart
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
LOGICAL               :: DoElement(nElems)
REAL                  :: vBulk(3), dens
!===================================================================================================================================
DSMC_RHS = 0.0
DoElement = .FALSE.

DO iElem = 1, nElems
  nPart = PEM%pNumber(iElem)
  dens = nPart * Species(1)%MacroParticleFactor / GEO%Volume(iElem)
  IF (dens.LT.FPDSMCSwitchDens) THEN
    DoElement(iElem) = .TRUE.
    CYCLE
  END IF
  IF (nPart.LT.3) CYCLE

  IF (DoBGKCellAdaptation) THEN
    CALL BGK_octree_adapt(iElem)
  ELSE
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 1.
    END IF

    ALLOCATE(iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing

    vBulk(1:3) = 0.0
    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      vBulk(1:3)  =  vBulk(1:3) + PartState(4:6,iPart)
      iPart = PEM%pNext(iPart)
    END DO
    vBulk = vBulk / nPart

    IF(DSMC%CalcQualityFactors) THEN
      FP_MeanRelaxFactorCounter=0; FP_MeanRelaxFactor=0.; FP_MaxRelaxFactor=0.; FP_MaxRotRelaxFactor=0.; FP_PrandtlNumber=0.
    END IF

    CALL FP_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk)
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
USE MOD_Mesh_Vars           ,ONLY: nElems, MeshFile
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_Particle_Vars       ,ONLY: PEM, PartState, WriteMacroVolumeValues, WriteMacroSurfaceValues
USE MOD_FP_CollOperator     ,ONLY: FP_CollisionOperator
USE MOD_DSMC_Vars           ,ONLY: DSMC_RHS, DSMC, SamplingActive
USE MOD_BGK_Vars            ,ONLY: DoBGKCellAdaptation
USE MOD_BGK_Adaptation      ,ONLY: BGK_octree_adapt
USE MOD_DSMC_Analyze        ,ONLY: DSMCHO_data_sampling,WriteDSMCHOToHDF5,CalcSurfaceValues
USE MOD_Restart_Vars        ,ONLY: RestartTime
USE MOD_FPFlow_Vars         ,ONLY: FP_QualityFacSamp, FP_PrandtlNumber
USE MOD_FPFlow_Vars         ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, nPart, iPart, iLoop, nOutput
REAL                :: vBulk(1:3)
INTEGER, ALLOCATABLE   :: iPartIndx_Node(:)
!===================================================================================================================================
DSMC_RHS = 0.0

IF (DoBGKCellAdaptation) THEN
  DO iElem = 1, nElems
    CALL BGK_octree_adapt(iElem)
  END DO
ELSE
  DO iElem = 1, nElems
    nPart = PEM%pNumber(iElem)
    IF (nPart.LT.3) CYCLE

    ALLOCATE(iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing

    vBulk(1:3) = 0.0
    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      vBulk(1:3)  =  vBulk(1:3) + PartState(4:6,iPart)
      iPart = PEM%pNext(iPart)
    END DO
    vBulk = vBulk / nPart

    IF(DSMC%CalcQualityFactors) THEN
      FP_MeanRelaxFactorCounter=0; FP_MeanRelaxFactor=0.; FP_MaxRelaxFactor=0.; FP_MaxRotRelaxFactor=0.; FP_PrandtlNumber=0.
    END IF

    CALL FP_CollisionOperator(iPartIndx_Node, nPart, GEO%Volume(iElem), vBulk)
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

END SUBROUTINE FPFlow_main

END MODULE MOD_FPFLOW
