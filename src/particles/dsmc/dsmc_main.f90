#include "boltzplatz.h"

MODULE MOD_DSMC
!===================================================================================================================================
! Module for DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_main
  MODULE PROCEDURE DSMC_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_main
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_main()
!===================================================================================================================================
! Performs DSMC routines (containing loop over all cells)
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_BGGas,            ONLY : DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, DSMCSumOfFormedParticles, BGGas, CollisMode
  USE MOD_DSMC_Vars,             ONLY : ChemReac
  USE MOD_Particle_Vars,         ONLY : PEM, Time, PDM, usevMPF
  USE MOD_Particle_Analyze_Vars, ONLY : CalcEkin
  !USE MOD_DSMC_Analyze,          ONLY : DSMC_data_sampling,DSMC_output_calc,CalcSurfaceValues,WriteOutputMeshSamp
  USE MOD_DSMC_Analyze,          ONLY : DSMC_data_sampling,DSMC_output_calc,WriteOutputMeshSamp
  USE MOD_TimeDisc_Vars,         ONLY : TEnd
  USE MOD_DSMC_ChemReact,        ONLY : SetMeanVibQua
  USE MOD_DSMC_ParticlePairing,  ONLY : DSMC_pairing_octree, DSMC_pairing_statistical
  USE MOD_DSMC_CollisionProb,    ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,           ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,        ONLY : DSMC_vmpf_prob
#if (PP_TimeDiscMethod==1001)
  USE MOD_LD_Vars,               ONLY : BulkValues, LD_DSMC_RHS
#endif
#if (PP_TimeDiscMethod!=1001) /* --- LD-DSMC Output in timedisc */
  USE MOD_Particle_Vars,         ONLY : WriteMacroValues
  USE MOD_Restart_Vars,          ONLY : RestartTime
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: iElem, nPart, nPair, iPair
  REAL              :: iRan
#if (PP_TimeDiscMethod!=1001)
  INTEGER           :: nOutput
#endif
!===================================================================================================================================

  DSMC_RHS(1:PDM%ParticleVecLength,1) = 0
  DSMC_RHS(1:PDM%ParticleVecLength,2) = 0
  DSMC_RHS(1:PDM%ParticleVecLength,3) = 0
  DSMCSumOfFormedParticles =0

  IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_InitBGGas 
  DO iElem = 1, nElems ! element/cell main loop
#if (PP_TimeDiscMethod==1001)
    IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN  ! --- DSMC Cell ?
#endif
    DSMC%CollMean = 0.0
    DSMC%CollMeanCount = 0
    IF(BGGas%BGGasSpecies.NE.0) THEN
      CALL DSMC_pairing_bggas(iElem)
    ELSE IF (DSMC%UseOctree) THEN
      CALL DSMC_pairing_octree(iElem)
    ELSE
      CALL DSMC_pairing_statistical(iElem)  ! pairing of particles per cell
    END IF

    IF (.NOT.DSMC%UseOctree) THEN                                                               ! no octree
      !Calc the mean evib per cell and iter, necessary for dissociation probability
      IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) THEN
        CALL SetMeanVibQua()
      END IF
      
      nPart = PEM%pNumber(iElem)
      nPair = int(nPart/2)

      DO iPair = 1, nPair
        IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
          IF (usevMPF.AND.(BGGas%BGGasSpecies.EQ.0)) THEN            ! calculation of collision prob
            CALL DSMC_vmpf_prob(iElem, iPair)
          ELSE
            CALL DSMC_prob_calc(iElem, iPair)
          END IF
          CALL RANDOM_NUMBER(iRan)
          IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
#if (PP_TimeDiscMethod==42)
            IF(CalcEkin.OR.DSMC%ReservoirSimu) THEN
#else
            IF(CalcEkin) THEN
#endif
              DSMC%NumColl(Coll_pData(iPair)%PairType) = DSMC%NumColl(Coll_pData(iPair)%PairType) + 1
              DSMC%NumColl(CollInf%NumCase + 1) = DSMC%NumColl(CollInf%NumCase + 1) + 1
            END IF
            CALL DSMC_perform_collision(iPair,iElem)
          END IF
        END IF
      END DO
      DEALLOCATE(Coll_pData)
    END IF                                                                                     ! no end octree
    IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN ! mean collision probability of all collision pairs
      IF (DSMC%CollMeanCount.GT.0) DSMC%CollProbOut(iElem,2) = DSMC%CollProbOut(iElem,2) + DSMC%CollMean / REAL(DSMC%CollMeanCount)
    END IF
#if (PP_TimeDiscMethod==1001)
    END IF  ! --- END DSMC Cell?
#endif
  END DO ! iElem Loop
! Output!
#if (PP_TimeDiscMethod!=1001) /* --- LD-DSMC Output in timedisc */
  PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles 
  IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_FinalizeBGGas 
#if (PP_TimeDiscMethod==42)
  IF ((.NOT.DSMC%ReservoirSimu).AND.(.NOT.WriteMacroValues)) THEN
#else
  IF (.NOT.WriteMacroValues) THEN
#endif
    IF(Time.GE.(1-DSMC%TimeFracSamp)*TEnd) THEN
      CALL DSMC_data_sampling()  ! Data sampling for output
      IF(DSMC%NumOutput.NE.0) THEN
        nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
        IF(Time.GE.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
          DSMC%NumOutput = DSMC%NumOutput - 1
          ! Skipping outputs immediately after the first few iterations
          IF(RestartTime.LT.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * REAL(nOutput))) THEN 
            CALL DSMC_output_calc
            IF (DSMC%OutputMeshSamp) CALL WriteOutputMeshSamp() !EmType6
            !IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues
          END IF
        END IF
      END IF
    END IF
  END IF
#else /* --- LD-DSMC? */
  LD_DSMC_RHS(1:PDM%ParticleVecLength,1) = DSMC_RHS(1:PDM%ParticleVecLength,1)
  LD_DSMC_RHS(1:PDM%ParticleVecLength,2) = DSMC_RHS(1:PDM%ParticleVecLength,2)
  LD_DSMC_RHS(1:PDM%ParticleVecLength,3) = DSMC_RHS(1:PDM%ParticleVecLength,3)
#endif /* --- END LD-DSMC? */
END SUBROUTINE DSMC_main

END MODULE MOD_DSMC
