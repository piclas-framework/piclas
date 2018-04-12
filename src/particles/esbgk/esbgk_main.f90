#include "boltzplatz.h"

MODULE MOD_ESBGK
!===================================================================================================================================
! Module for ESBGK Flow
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ESBGK_main
  MODULE PROCEDURE ESBGK_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ESBGK_main
!===================================================================================================================================

CONTAINS


SUBROUTINE ESBGK_main()
!===================================================================================================================================
!> Performs ESBGK Momentum Evaluation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars      ,ONLY: TEnd, Time
USE MOD_Mesh_Vars          ,ONLY: nElems, MeshFile
USE MOD_DSMC_Vars          ,ONLY: DSMC_RHS, DSMC, SamplingActive
USE MOD_ESBGK_Adaptation   ,ONLY: ESBGK_octree_adapt!, ESBGKSplitCells
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_Particle_Vars      ,ONLY: PEM, PartState, WriteMacroVolumeValues, WriteMacroSurfaceValues
USE MOD_Restart_Vars       ,ONLY: RestartTime
USE MOD_ESBGK_Vars         ,ONLY: DoBGKCellAdaptation, BGKDoAveraging, ElemNodeAveraging, BGKAveragingLength
USE MOD_ESBGK_Vars         ,ONLY: DoBGKCellSplitting
USE MOD_ESBGK_CollOperator ,ONLY: ESBGK_CollisionOperatorOctree
USE MOD_DSMC_Analyze       ,ONLY: DSMCHO_data_sampling,CalcSurfaceValues,WriteDSMCHOToHDF5
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart, nOutput
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
REAL                  :: vBulk(3)
!===================================================================================================================================
DSMC_RHS = 0.0

IF (DoBGKCellAdaptation) THEN
  DO iElem = 1, nElems
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 1.
    END IF
    CALL ESBGK_octree_adapt(iElem)
    IF(DSMC%CalcQualityFactors) THEN
      IF(Time.GE.(1-DSMC%TimeFracSamp)*TEnd) THEN
        DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
      END IF
    END IF
  END DO
!ELSE IF (DoBGKCellSplitting) THEN
!  DO iElem = 1, nElems
!    CALL ESBGKSplitCells(iElem)
!  END DO
ELSE
  DO iElem = 1, nElems
    nPart = PEM%pNumber(iElem)
    IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE

    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 1.
    END IF

    ALLOCATE(iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing

    vBulk(1:3) = 0.0
    iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
    DO iLoop = 1, nPart
      iPartIndx_Node(iLoop) = iPart
      vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)
      iPart = PEM%pNext(iPart)
    END DO
    vBulk = vBulk / nPart

    IF (BGKDoAveraging) THEN
      CALL ESBGK_CollisionOperatorOctree(iPartIndx_Node, nPart, iElem, GEO%Volume(iElem), vBulk, &
          ElemNodeAveraging(iElem)%Root%AverageValues(1:5,1:BGKAveragingLength), &
               CorrectStep = ElemNodeAveraging(iElem)%Root%CorrectStep)
    ELSE
      CALL ESBGK_CollisionOperatorOctree(iPartIndx_Node, nPart, iElem, GEO%Volume(iElem), vBulk)
    END IF
    IF(DSMC%CalcQualityFactors) THEN
      IF(Time.GE.(1-DSMC%TimeFracSamp)*TEnd) THEN
        DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
      END IF
    END IF
    DEALLOCATE(iPartIndx_Node)
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

END SUBROUTINE ESBGK_main

SUBROUTINE BGKEuler_main()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_Particle_Vars      ,ONLY: PEM, PartState
USE MOD_DSMC_Vars          ,ONLY: DSMC_RHS
USE MOD_ESBGK_Vars         ,ONLY: DoBGKCellAdaptation
USE MOD_ESBGK_CollOperator ,ONLY: ESBGK_Euler
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
REAL                  :: vBulk(3)
!===================================================================================================================================
DSMC_RHS = 0.0

DO iElem = 1, nElems
  nPart = PEM%pNumber(iElem)
  IF ((nPart.EQ.0).OR.(nPart.EQ.1)) CYCLE

  ALLOCATE(iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing

  vBulk(1:3) = 0.0
  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iPartIndx_Node(iLoop) = iPart
    vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)
    iPart = PEM%pNext(iPart)
  END DO
  vBulk = vBulk / nPart

  CALL ESBGK_Euler(iPartIndx_Node(1:nPart), nPart, iElem, GEO%Volume(iElem), vBulk)
  DEALLOCATE(iPartIndx_Node)
END DO

END SUBROUTINE BGKEuler_main


END MODULE MOD_ESBGK
