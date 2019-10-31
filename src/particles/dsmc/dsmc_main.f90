!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

SUBROUTINE DSMC_main(DoElement)
!===================================================================================================================================
!> Performs DSMC routines (containing loop over all cells)
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: time, TEnd
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_BGGas            ,ONLY: DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, DSMC_RHS, DSMC, CollInf, DSMCSumOfFormedParticles, BGGas, CollisMode
USE MOD_DSMC_Vars             ,ONLY: ChemReac, SpecDSMC, VarVibRelaxProb, ConsiderVolumePortions
USE MOD_DSMC_Analyze          ,ONLY: CalcMeanFreePath
USE MOD_DSMC_SteadyState      ,ONLY: QCrit_evaluation, SteadyStateDetection_main
USE MOD_Particle_Vars         ,ONLY: PEM, PDM, WriteMacroVolumeValues, nSpecies, Symmetry2D, PartSpecies
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_Particle_Analyze_Vars ,ONLY: CalcEkin
USE MOD_DSMC_Analyze          ,ONLY: DSMCHO_data_sampling,CalcSurfaceValues, WriteDSMCHOToHDF5, CalcGammaVib, &
                                     SamplingRotVibRelaxProb
USE MOD_DSMC_Relaxation       ,ONLY: SetMeanVibQua
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_pairing_octree, DSMC_pairing_statistical, DSMC_pairing_quadtree
USE MOD_DSMC_CollisionProb    ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Collis           ,ONLY: DSMC_perform_collision, DSMC_calc_var_P_vib
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_Mesh_Vars             ,ONLY: MeshFile
USE MOD_TimeDisc_Vars         ,ONLY: iter
USE MOD_DSMC_Vars             ,ONLY: UseQCrit, SamplingActive, QCritTestStep, QCritLastTest, UseSSD
USE MOD_Particle_Vars         ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars,         ONLY: KeepWallParticles
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,OPTIONAL  :: DoElement(nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: iElem, nPart, nPair, iPair, cSpec1, cSpec2, iSpec
  REAL              :: iRan, VibProb
INTEGER           :: nOutput
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

IF(.NOT.PRESENT(DoElement)) THEN
  DSMC_RHS(1,1:PDM%ParticleVecLength) = 0
  DSMC_RHS(2,1:PDM%ParticleVecLength) = 0
  DSMC_RHS(3,1:PDM%ParticleVecLength) = 0
END IF
DSMCSumOfFormedParticles = 0

IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_InitBGGas
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO iElem = 1, nElems ! element/cell main loop
  IF(PRESENT(DoElement)) THEN
    IF (.NOT.DoElement(iElem)) CYCLE
  END IF
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 0.0; DSMC%CollProbMean = 0.0; DSMC%CollProbMeanCount = 0; DSMC%CollSepDist = 0.0; DSMC%CollSepCount = 0
      DSMC%MeanFreePath = 0.0; DSMC%MCSoverMFP = 0.0
      IF(DSMC%RotRelaxProb.GT.2) THEN
        DSMC%CalcRotProb = 0.
      END IF
      IF(DSMC%VibRelaxProb.EQ.2) THEN
        DSMC%CalcVibProb = 0.
      END IF
    END IF
    IF (CollisMode.NE.0) THEN
      ChemReac%nPairForRec = 0
      IF(BGGas%BGGasSpecies.NE.0) THEN
        CALL DSMC_pairing_bggas(iElem)
      ELSE IF (DSMC%UseOctree) THEN
        IF(Symmetry2D) THEN
          CALL DSMC_pairing_quadtree(iElem)
        ELSE
          CALL DSMC_pairing_octree(iElem)
        END IF
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
          ! variable vibrational relaxation probability has to average of all collisions
          IF(DSMC%VibRelaxProb.EQ.2.0) THEN
            cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1)
            cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2)
            IF((SpecDSMC(cSpec1)%InterID.EQ.2).OR.(SpecDSMC(cSpec1)%InterID.EQ.20)) THEN
              CALL DSMC_calc_var_P_vib(cSpec1,cSpec2,iPair,VibProb)
              VarVibRelaxProb%ProbVibAvNew(cSpec1) = VarVibRelaxProb%ProbVibAvNew(cSpec1) + VibProb
              VarVibRelaxProb%nCollis(cSpec1) = VarVibRelaxProb%nCollis(cSpec1) + 1
              IF(DSMC%CalcQualityFactors) THEN
                DSMC%CalcVibProb(cSpec1,2) = MAX(DSMC%CalcVibProb(cSpec1,2),VibProb)
              END IF
            END IF
            IF((SpecDSMC(cSpec2)%InterID.EQ.2).OR.(SpecDSMC(cSpec2)%InterID.EQ.20)) THEN
              CALL DSMC_calc_var_P_vib(cSpec2,cSpec1,iPair,VibProb)
              VarVibRelaxProb%ProbVibAvNew(cSpec2) = VarVibRelaxProb%ProbVibAvNew(cSpec2) + VibProb
              VarVibRelaxProb%nCollis(cSpec2) = VarVibRelaxProb%nCollis(cSpec2) + 1
              IF(DSMC%CalcQualityFactors) THEN
                DSMC%CalcVibProb(cSpec2,2) = MAX(DSMC%CalcVibProb(cSpec2,2),VibProb)
              END IF
            END IF
          END IF
          IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
            CALL DSMC_prob_calc(iElem, iPair)
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
        IF(DSMC%CalcQualityFactors) THEN
          IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
            ! Calculation of the mean free path
            IF (ConsiderVolumePortions) THEN
              DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum)&
                          ,GEO%Volume(iElem)*(1.-GEO%MPVolumePortion(iElem)),SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
            ELSE
              DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum),GEO%Volume(iElem), &
                                                    SpecDSMC(1)%omegaVHS,DSMC%InstantTransTemp(nSpecies+1))
            END IF
            ! Determination of the MCS/MFP for the case without octree
            IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                            / DSMC%MeanFreePath
          END IF
        END IF
        DEALLOCATE(Coll_pData)
        IF(DSMC%VibRelaxProb.EQ.2.0) THEN
          DO iSpec=1,nSpecies
            IF(VarVibRelaxProb%nCollis(iSpec).NE.0) THEN ! Calc new vibrational relaxation probability
              VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VarVibRelaxProb%ProbVibAv(iElem,iSpec) &
                                                     * VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec)) &
                                                     + (1.-VarVibRelaxProb%alpha**(VarVibRelaxProb%nCollis(iSpec))) &
                                                     / (VarVibRelaxProb%nCollis(iSpec)) * VarVibRelaxProb%ProbVibAvNew(iSpec)
            END IF
          END DO
        END IF
      END IF ! end no octree
      IF(DSMC%CalcQualityFactors) THEN
        IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
            ! mean collision probability of all collision pairs
            IF(DSMC%CollProbMeanCount.GT.0) THEN
              DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
              DSMC%QualityFacSamp(iElem,2) = DSMC%QualityFacSamp(iElem,2) + DSMC%CollProbMean / REAL(DSMC%CollProbMeanCount)
            END IF
            ! mean collision separation distance of actual collisions
            IF(DSMC%CollSepCount.GT.0) DSMC%QualityFacSamp(iElem,3) = DSMC%QualityFacSamp(iElem,3) + DSMC%MCSoverMFP
            ! Counting sample size
            DSMC%QualityFacSamp(iElem,4) = DSMC%QualityFacSamp(iElem,4) + 1.
            ! Sample rotation relaxation probability
            IF((DSMC%RotRelaxProb.EQ.2).OR.(DSMC%VibRelaxProb.EQ.2)) CALL SamplingRotVibRelaxProb(iElem)
          END IF
        END IF
      END IF  ! --- CollisMode.NE.0
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem Loop
  ! Output!
  PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles
  IF(BGGas%BGGasSpecies.NE.0) CALL DSMC_FinalizeBGGas
#if (PP_TimeDiscMethod==42)
    IF ((.NOT.DSMC%ReservoirSimu).AND.(.NOT.WriteMacroVolumeValues).AND.(.NOT.WriteMacroSurfaceValues)) THEN
#else
    IF (.NOT.WriteMacroVolumeValues .AND. .NOT.WriteMacroSurfaceValues) THEN
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
END SUBROUTINE DSMC_main

END MODULE MOD_DSMC
