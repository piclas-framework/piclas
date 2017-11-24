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
  USE MOD_TimeDisc_Vars,         ONLY : time, iter, TEnd
  USE MOD_Globals
  USE MOD_DSMC_BGGas,            ONLY : DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
  USE MOD_Mesh_Vars,             ONLY : nElems, MeshFile
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, DSMCSumOfFormedParticles, BGGas, CollisMode
  USE MOD_DSMC_Vars,             ONLY : ChemReac
  USE MOD_DSMC_Vars,             ONLY : UseQCrit, SamplingActive, QCritTestStep, QCritLastTest, UseSSD
  USE MOD_DSMC_SteadyState,      ONLY : QCrit_evaluation, SteadyStateDetection_main
  USE MOD_Particle_Vars,         ONLY : PEM, PDM, usevMPF, BoltzmannConst, WriteMacroVolumeValues
  USE MOD_Particle_Analyze_Vars, ONLY : CalcEkin
  USE MOD_DSMC_Analyze,          ONLY : DSMCHO_data_sampling,CalcSurfaceValues, WriteDSMCHOToHDF5, CalcGammaVib
  USE MOD_DSMC_Relaxation,       ONLY : SetMeanVibQua
  USE MOD_DSMC_ParticlePairing,  ONLY : DSMC_pairing_octree, DSMC_pairing_statistical
  USE MOD_DSMC_CollisionProb,    ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,           ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,        ONLY : DSMC_vmpf_prob
  USE MOD_Particle_Vars,         ONLY : KeepWallParticles
#if (PP_TimeDiscMethod==1001)
  USE MOD_LD_Vars,               ONLY : BulkValues, LD_DSMC_RHS
#endif
#if (PP_TimeDiscMethod!=1001) /* --- LD-DSMC Output in timedisc */
  USE MOD_Particle_Vars,         ONLY : WriteMacroVolumeValues
  USE MOD_Restart_Vars,          ONLY : RestartTime
#endif
#ifdef MPI
  USE MOD_LoadBalance_Vars,      ONLY: ElemTime
#endif /*MPI*/
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
#ifdef MPI
! load balance
  REAL                             :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

  DSMC_RHS(1:PDM%ParticleVecLength,1) = 0
  DSMC_RHS(1:PDM%ParticleVecLength,2) = 0
  DSMC_RHS(1:PDM%ParticleVecLength,3) = 0
  DSMCSumOfFormedParticles = 0

  IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_InitBGGas 
  DO iElem = 1, nElems ! element/cell main loop
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
#if (PP_TimeDiscMethod==1001)
    IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN  ! --- DSMC Cell ?
#endif
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 0.0 
      DSMC%CollProbMean = 0.0
      DSMC%CollProbMeanCount = 0
      DSMC%CollSepDist = 0.0
      DSMC%CollSepCount = 0
    END IF
#if (PP_TimeDiscMethod==42)
    IF (ChemReac%NumOfReact.GT.0) THEN
      ChemReac%ReacCount = 0
      ChemReac%ReacCollMean = 0.0
      ChemReac%ReacCollMeanCount = 0
    END IF
#endif
    IF (CollisMode.NE.0) THEN
      ChemReac%nPairForRec = 0
      IF(BGGas%BGGasSpecies.NE.0) THEN
        CALL DSMC_pairing_bggas(iElem)
      ELSE IF (DSMC%UseOctree) THEN
        CALL DSMC_pairing_octree(iElem)
      ELSE
        CALL DSMC_pairing_statistical(iElem)  ! pairing of particles per cell
      END IF

      IF (.NOT.DSMC%UseOctree) THEN                                                               ! no octree
        !Calc the mean evib per cell and iter, necessary for dissociation probability
        IF (CollisMode.EQ.3) THEN
          CALL SetMeanVibQua()
        END IF

        IF (KeepWallParticles) THEN
          nPart = PEM%pNumber(iElem)-PEM%wNumber(iElem)
        ELSE
          nPart = PEM%pNumber(iElem)
        END IF
        nPair = INT(nPart/2)

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
      END IF                                                                                     ! end no octree
      IF(DSMC%CalcQualityFactors) THEN
        IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
            ! mean collision probability of all collision pairs
            IF(DSMC%CollProbMeanCount.GT.0) THEN
              DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
              DSMC%QualityFacSamp(iElem,2) = DSMC%QualityFacSamp(iElem,2) + DSMC%CollProbMean / REAL(DSMC%CollProbMeanCount)
            END IF
            ! mean collision separation distance of actual collisions
            IF(DSMC%CollSepCount.GT.0) DSMC%QualityFacSamp(iElem,3) = DSMC%QualityFacSamp(iElem,3) &
                                                                            + DSMC%CollSepDist / REAL(DSMC%CollSepCount)
        END IF
      END IF
    END IF  ! --- CollisMode.NE.0
#if (PP_TimeDiscMethod==1001)
    END IF  ! --- END DSMC Cell?
#endif
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
  END DO ! iElem Loop
! Output!
#if (PP_TimeDiscMethod!=1001) /* --- LD-DSMC Output in timedisc */
  PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles 
  IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_FinalizeBGGas 
#if (PP_TimeDiscMethod==42)
  IF ((.NOT.DSMC%ReservoirSimu).AND.(.NOT.WriteMacroVolumeValues)) THEN
#else
  IF (.NOT.WriteMacroVolumeValues) THEN
#endif
    IF(UseQCrit) THEN
      ! Use QCriterion (Burt,Boyd) for steady - state detection
      IF((.NOT.SamplingActive).AND.(iter-QCritLastTest.EQ.QCritTestStep)) THEN
        CALL QCrit_evaluation()
        QCritLastTest=iter
        IF(SamplingActive) THEN
          SWRITE(*,*)'Sampling active'
          ! Set TimeFracSamp and DeltaTimeOutput -> correct number of outputs
          DSMC%TimeFracSamp = (TEnd-Time)/TEnd
          DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
        ENDIF
      ENDIF
    ELSEIF(UseSSD) THEN
      ! Use SSD for steady - state detection
      IF((.NOT.SamplingActive)) THEN
        CALL SteadyStateDetection_main()
        IF(SamplingActive) THEN
          SWRITE(*,*)'Sampling active'
          ! Set TimeFracSamp and DeltaTimeOutput -> correct number of outputs
          DSMC%TimeFracSamp = (TEnd-Time)/TEnd
          DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
        ENDIF
      ENDIF
    ELSE
      ! Use user given TimeFracSamp
      IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).AND.(.NOT.SamplingActive))  THEN
        SamplingActive=.TRUE.
        SWRITE(*,*)'Sampling active'
      ENDIF
    ENDIF
    !
    ! Calculate Entropy using Theorem of Boltzmann
    !CALL EntropyCalculation()
    !

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
  END IF
#else /* --- LD-DSMC? */
  LD_DSMC_RHS(1:PDM%ParticleVecLength,1) = DSMC_RHS(1:PDM%ParticleVecLength,1)
  LD_DSMC_RHS(1:PDM%ParticleVecLength,2) = DSMC_RHS(1:PDM%ParticleVecLength,2)
  LD_DSMC_RHS(1:PDM%ParticleVecLength,3) = DSMC_RHS(1:PDM%ParticleVecLength,3)
#endif /* --- END LD-DSMC? */
END SUBROUTINE DSMC_main

END MODULE MOD_DSMC
