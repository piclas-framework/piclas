MODULE MOD_DSMC
!===================================================================================================================================
! module including collisions
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

PUBLIC :: DSMC_main
!===================================================================================================================================

CONTAINS


SUBROUTINE DSMC_main()

  USE MOD_DSMC_BGGas,            ONLY : DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, DSMCSumOfFormedParticles, BGGas, CollisMode
  USE MOD_DSMC_Vars,             ONLY : ChemReac, CollMean
  USE MOD_Particle_Vars,         ONLY : PEM, Time, PDM, nSpecies, WriteMacroValues, usevMPF
  USE MOD_Particle_Analyze_Vars, ONLY : CalcEkin
  USE MOD_DSMC_Analyze,          ONLY : DSMC_data_sampling, DSMC_output_calc, CalcSurfaceValues, OutputMaxCollProb
  USE MOD_TimeDisc_Vars,         ONLY : TEnd
  USE MOD_DSMC_ChemReact,        ONLY : SetMeanVibQua
  USE MOD_DSMC_ParticlePairing,  ONLY : DSMC_pairing_octree, DSMC_pairing_statistical
  USE MOD_DSMC_CollisionProb,    ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,           ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,        ONLY : DSMC_vmpf_prob

!--------------------------------------------------------------------------------------------------!
! main DSMC routine
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER           :: iElem, nPart, nPair, iPair, ProbMeanCellCount
  REAL              :: iRan, CollMeanCell, ProbMeanCell
  INTEGER           :: nOutput
!--------------------------------------------------------------------------------------------------!

DSMC_RHS(1:PDM%ParticleVecLength,1) = 0
DSMC_RHS(1:PDM%ParticleVecLength,2) = 0
DSMC_RHS(1:PDM%ParticleVecLength,3) = 0
DSMCSumOfFormedParticles =0

IF (DSMC%CollProbMaxOut) THEN
  DSMC%CollProbMax = 0.0  
  DSMC%CollMean = 0.0
  DSMC%CollMeanCount = 0.0
END IF
IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_InitBGGas 
DO iElem = 1, nElems ! element/cell main loop 
  CollMeanCell=0
  ProbMeanCell=0
  ProbMeanCellCount=0
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
          IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN
            ProbMeanCell = ProbMeanCell + Coll_pData(iPair)%Prob    ! collision probabilities of actual collisions
            ProbMeanCellCount = ProbMeanCellCount + 1               ! counting collisions
          END IF
          CALL DSMC_perform_collision(iPair,iElem)
        END IF
      END IF
      
      IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN                   ! all collision probabilities
        CollMeanCell = CollMeanCell + Coll_pData(iPair)%Prob
      END IF

    END DO
    IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN ! mean collision probability of cell
      IF (ProbMeanCellCount.eq.0) THEN
        ELSE                                              
        CollMean(iElem,1)=CollMean(iElem,1)+ProbMeanCell/ProbMeanCellCount
      END IF
        CollMean(iElem,2)=CollMean(iElem,2)+CollMeanCell/nPair
    END IF
    DEALLOCATE(Coll_pData)
  END IF                                                                                     ! no end octree
END DO
IF (DSMC%CollProbMaxOut) CALL OutputMaxCollProb(Time)
PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles 
IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_FinalizeBGGas 
! Output!
#if (PP_TimeDiscMethod==42)
IF ((.NOT.DSMC%ReservoirSimu).AND.(.NOT.WriteMacroValues)) THEN
#else
IF (.NOT.WriteMacroValues) THEN
#endif
  IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN
    CALL DSMC_data_sampling()  ! Data sampling for output
    IF(DSMC%NumOutput.NE.0) THEN
      nOutput = (DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput-DSMC%NumOutput + 1
      IF(Time.ge.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
        DSMC%NumOutput = DSMC%NumOutput - 1
        CALL DSMC_output_calc(nOutput)
        IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues(nOutput)
      END IF
    END IF
  END IF
END IF
END SUBROUTINE DSMC_main

!--------------------------------------------------------------------------------------------------!


END MODULE MOD_DSMC
