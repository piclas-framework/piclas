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
USE MOD_Globals
USE MOD_DSMC_BGGas            ,ONLY: BGGas_InsertParticles, DSMC_pairing_bggas, BGGas_DeleteParticles
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_DSMC_Vars             ,ONLY: DSMC, CollInf, DSMCSumOfFormedParticles, BGGas, CollisMode, ElecRelaxPart
USE MOD_DSMC_Analyze          ,ONLY: SummarizeQualityFactors, DSMCMacroSampling
USE MOD_DSMC_Relaxation       ,ONLY: FinalizeCalcVibRelaxProb, InitCalcVibRelaxProb
USE MOD_Particle_Vars         ,ONLY: PEM, PDM, WriteMacroVolumeValues, Symmetry
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_pairing_standard, DSMC_pairing_octree, DSMC_pairing_quadtree, DSMC_pairing_dotree
USE MOD_Particle_Vars         ,ONLY: WriteMacroSurfaceValues
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
USE MOD_MCC_Vars              ,ONLY: UseMCC
USE MOD_MCC                   ,ONLY: MonteCarloCollision
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,OPTIONAL  :: DoElement(nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, nPart
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! Reset the number of particles created during the DSMC loop
DSMCSumOfFormedParticles = 0
! Insert background gas particles for every simulation particle (except when using MCC and free-molecular flow)
IF((BGGas%NumberOfSpecies.GT.0).AND.(.NOT.UseMCC).AND.(CollisMode.NE.0)) CALL BGGas_InsertParticles()

IF ((DSMC%ElectronicModel.EQ.4).AND.(CollisMode.EQ.3)) THEN
  ElecRelaxPart(1:PDM%ParticleVecLength) = .TRUE.
END IF

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF (CollisMode.NE.0) THEN
  DO iElem = 1, nElems ! element/cell main loop
    IF(PRESENT(DoElement)) THEN
      IF (.NOT.DoElement(iElem)) CYCLE
    END IF
    nPart = PEM%pNumber(iElem)
    IF (nPart.LT.1) CYCLE
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CollProbMax = 0.0; DSMC%CollProbMean = 0.0; DSMC%CollProbMeanCount = 0; DSMC%CollSepDist = 0.0; DSMC%CollSepCount = 0
      DSMC%MeanFreePath = 0.0; DSMC%MCSoverMFP = 0.0
      IF(DSMC%RotRelaxProb.GT.2) DSMC%CalcRotProb = 0.
      DSMC%CalcVibProb = 0.
    END IF
    CALL InitCalcVibRelaxProb()
    IF(BGGas%NumberOfSpecies.GT.0) THEN
      ! Decide between MCC and DSMC-based background gas
      IF(UseMCC) THEN
        CALL MonteCarloCollision(iElem)
      ELSE
        CALL DSMC_pairing_bggas(iElem)
      END IF
    ELSE IF (nPart.GT.1) THEN
      IF (DSMC%UseOctree) THEN
        ! On-the-fly cell refinement and pairing within subcells
        IF(Symmetry%Order.EQ.3) THEN
          CALL DSMC_pairing_octree(iElem)
        ELSE IF(Symmetry%Order.EQ.2) THEN
          CALL DSMC_pairing_quadtree(iElem)
        ELSE
          CALL DSMC_pairing_dotree(iElem)
        END IF
      ELSE ! NOT DSMC%UseOctree
        ! Standard pairing of particles within a cell
        CALL DSMC_pairing_standard(iElem)
      END IF ! DSMC%UseOctree
    ELSE ! less than 2 particles
      IF (CollInf%ProhibitDoubleColl.AND.(nPart.EQ.1)) CollInf%OldCollPartner(PEM%pStart(iElem)) = 0
      CYCLE ! next element
    END IF
    CALL FinalizeCalcVibRelaxProb(iElem)
    IF(DSMC%CalcQualityFactors) CALL SummarizeQualityFactors(iElem)
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem Loop
END IF ! CollisMode.NE.0

IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) THEN
  CALL Abort(__STAMP__&
    ,'ERROR in DSMC: ParticleVecLength greater than MaxParticleNumber! Increase the MaxParticleNumber to at least: ' &
    , IntInfoOpt=PDM%ParticleVecLength)
END IF

! Delete background gas particles
IF((BGGas%NumberOfSpecies.GT.0).AND.(.NOT.UseMCC)) CALL BGGas_DeleteParticles()

! Sampling of macroscopic values
! (here for a continuous average; average over N iterations is performed in src/analyze/analyze.f90)
IF (.NOT.WriteMacroVolumeValues .AND. .NOT.WriteMacroSurfaceValues) THEN
  CALL DSMCMacroSampling()
END IF

END SUBROUTINE DSMC_main

END MODULE MOD_DSMC
