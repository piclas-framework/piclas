!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Macro_Restart
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC         :: MacroRestart_InsertParticles
!===================================================================================================================================
CONTAINS

SUBROUTINE MacroRestart_InsertParticles()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: Pi
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting, DSMC
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Particle_Vars           ,ONLY: Species, PDM, nSpecies, PartState, VarTimeStep
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,BoundsOfElem_Shared
USE MOD_Particle_Tracking       ,ONLY: ParticleInsideCheck
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iElem,iSpec,iPart,nPart,locnPart,iHeight,yPartitions,GlobalElemID
REAL                                :: iRan, RandomPos(3), PartDens, TempMPF, MaxPosTemp, MinPosTemp
REAL                                :: TempVol, Volume
LOGICAL                             :: InsideFlag
!===================================================================================================================================

SWRITE(UNIT_stdOut,*) 'PERFORMING MACROSCOPIC RESTART...'

locnPart = 1

DO iElem = 1, nElems
  GlobalElemID = iElem + offsetElem
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GlobalElemID) ) ! 1-2: Min, Max value; 1-3: x,y,z
! #################### 2D ##########################################################################################################
    IF (Symmetry%Axisymmetric) THEN
      IF (RadialWeighting%DoRadialWeighting) THEN
        DO iSpec = 1, nSpecies
          IF (DSMC%DoAmbipolarDiff) THEN
            IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
          END IF
          yPartitions = 6
          PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS)
          ! Particle weighting
          DO iHeight = 1, yPartitions
            MinPosTemp = Bounds(1,2) + (Bounds(2,2) - Bounds(1,2))/ yPartitions *(iHeight-1.)
            MaxPosTemp = Bounds(1,2) + (Bounds(2,2) - Bounds(1,2))/ yPartitions *iHeight
            TempVol =  (MaxPosTemp-MinPosTemp)*(Bounds(2,1)-Bounds(1,1)) * Pi * (MaxPosTemp+MinPosTemp)
            TempMPF = CalcRadWeightMPF((MaxPosTemp+MinPosTemp)*0.5,iSpec)
            IF(VarTimeStep%UseVariableTimeStep) THEN
              TempMPF = TempMPF * CalcVarTimeStep((Bounds(2,1)+Bounds(1,1))*0.5, (MaxPosTemp+MinPosTemp)*0.5, iElem)
            END IF
            CALL RANDOM_NUMBER(iRan)
            nPart = INT(PartDens / TempMPF  * TempVol + iRan)
            DO iPart = 1, nPart
              InsideFlag=.FALSE.
              CALL RANDOM_NUMBER(RandomPos)
              RandomPos(1) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
              RandomPos(2) = MinPosTemp + RandomPos(2)*(MaxPosTemp-MinPosTemp)
              RandomPos(3) = 0.0
              InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
              IF (InsideFlag) THEN
                PartState(1:3,locnPart) = RandomPos(1:3)
                CALL MacroRestart_InitializeParticle_Maxwell(locnPart,iSpec,iElem)
                locnPart = locnPart + 1
              END IF
            END DO ! nPart
          END DO ! yPartitions
        END DO ! nSpecies
      ELSE ! No RadialWeighting
        DO iSpec = 1, nSpecies
          IF (DSMC%DoAmbipolarDiff) THEN
            IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
          END IF
          CALL RANDOM_NUMBER(iRan)
          TempMPF = Species(iSpec)%MacroParticleFactor
          IF(VarTimeStep%UseVariableTimeStep) THEN
            TempMPF = TempMPF * CalcVarTimeStep((Bounds(2,1)+Bounds(1,1))*0.5, (Bounds(2,2)+Bounds(1,2))*0.5, iElem)
          END IF
          nPart = INT(MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / TempMPF * ElemVolume_Shared(GlobalElemID) + iRan)
          DO iPart = 1, nPart
            InsideFlag=.FALSE.
            DO WHILE (.NOT.InsideFlag)
              CALL RANDOM_NUMBER(RandomPos)
              RandomPos(1) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
              RandomPos(2) = SQRT(RandomPos(2)*(Bounds(2,2)**2-Bounds(1,2)**2)+Bounds(1,2)**2)
              RandomPos(3) = 0.0
              InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
            END DO
            PartState(1:3,locnPart) = RandomPos(1:3)
            CALL MacroRestart_InitializeParticle_Maxwell(locnPart,iSpec,iElem)
            locnPart = locnPart + 1
          END DO ! nPart
        END DO ! nSpecies
      END IF ! RadialWeighting: YES/NO
    ELSE IF(Symmetry%Order.EQ.2) THEN
      Volume = (Bounds(2,2) - Bounds(1,2))*(Bounds(2,1) - Bounds(1,1))
      DO iSpec = 1, nSpecies
        IF (DSMC%DoAmbipolarDiff) THEN
          IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
        END IF
        CALL RANDOM_NUMBER(iRan)
        TempMPF = Species(iSpec)%MacroParticleFactor
        IF(VarTimeStep%UseVariableTimeStep) THEN
          TempMPF = TempMPF * CalcVarTimeStep((Bounds(2,1)+Bounds(1,1))*0.5, (Bounds(2,2)+Bounds(1,2))*0.5, iElem)
        END IF
        nPart = INT(MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / TempMPF * Volume + iRan)
        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          CALL RANDOM_NUMBER(RandomPos(1:2))
          RandomPos(1:2) = Bounds(1,1:2) + RandomPos(1:2)*(Bounds(2,1:2)-Bounds(1,1:2))
          RandomPos(3) = 0.0
          InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          IF (InsideFlag) THEN
            PartState(1:3,locnPart) = RandomPos(1:3)
            CALL MacroRestart_InitializeParticle_Maxwell(locnPart,iSpec,iElem)
            locnPart = locnPart + 1
          END IF
        END DO ! nPart
      END DO ! nSpecies
    ELSE IF(Symmetry%Order.EQ.1) THEN
      Volume = (Bounds(2,1) - Bounds(1,1))
      DO iSpec = 1, nSpecies
        IF (DSMC%DoAmbipolarDiff) THEN
          IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
        END IF
        CALL RANDOM_NUMBER(iRan)
        TempMPF = Species(iSpec)%MacroParticleFactor
        IF(VarTimeStep%UseVariableTimeStep) THEN
          TempMPF = TempMPF * CalcVarTimeStep((Bounds(2,1)+Bounds(1,1))*0.5, (Bounds(2,2)+Bounds(1,2))*0.5, iElem)
        END IF
        nPart = INT(MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / TempMPF * Volume + iRan)
        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          CALL RANDOM_NUMBER(RandomPos(1))
          RandomPos(1:2) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
          RandomPos(2) = 0.0
          RandomPos(3) = 0.0
          InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          IF (InsideFlag) THEN
            PartState(1:3,locnPart) = RandomPos(1:3)
            CALL MacroRestart_InitializeParticle_Maxwell(locnPart,iSpec,iElem)
            locnPart = locnPart + 1
          END IF
        END DO ! nPart
      END DO ! nSpecies
    ELSE
! #################### 3D ##########################################################################################################
      Volume = (Bounds(2,3) - Bounds(1,3))*(Bounds(2,2) - Bounds(1,2))*(Bounds(2,1) - Bounds(1,1))
      DO iSpec = 1, nSpecies
        IF (DSMC%DoAmbipolarDiff) THEN
          IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
        END IF
        CALL RANDOM_NUMBER(iRan)
        TempMPF = Species(iSpec)%MacroParticleFactor
        IF(VarTimeStep%UseVariableTimeStep) THEN
          TempMPF = TempMPF * CalcVarTimeStep(iElem=iElem)
        END IF
        nPart = INT(MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / TempMPF * Volume + iRan)
        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          CALL RANDOM_NUMBER(RandomPos)
          RandomPos(1:3) = Bounds(1,1:3) + RandomPos(1:3)*(Bounds(2,1:3)-Bounds(1,1:3))
          InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          IF (InsideFlag) THEN
            PartState(1:3,locnPart) = RandomPos(1:3)
            CALL MacroRestart_InitializeParticle_Maxwell(locnPart,iSpec,iElem)
            locnPart = locnPart + 1
          END IF
        END DO ! nPart
      END DO ! nSpecies
    END IF ! 1D/2D/Axisymmetric/3D
  END ASSOCIATE
END DO ! nElems

IF(locnPart.GE.PDM%maxParticleNumber) THEN
  CALL abort(__STAMP__,&
    'ERROR in MacroRestart: Increase maxParticleNumber!', locnPart)
END IF

PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart

END SUBROUTINE MacroRestart_InsertParticles


SUBROUTINE MacroRestart_InitializeParticle_Maxwell(iPart,iSpec,iElem)
!===================================================================================================================================
!> Initialize a particle from a given macroscopic result, requires the macroscopic velocity, translational and internal temperatures
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_Particle_Vars           ,ONLY: PDM, PartSpecies, PartState, PEM, VarTimeStep, PartMPF, Species
USE MOD_DSMC_Vars               ,ONLY: DSMC, PartStateIntEn, CollisMode, SpecDSMC, RadialWeighting, AmbipolElecVelo
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF, CalcEElec_particle, CalcEVib_particle, CalcERot_particle
USE MOD_part_tools              ,ONLY: CalcVelocity_maxwell_particle
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart, iSpec, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! 1) Set particle velocity from macroscopic bulk velocity and translational temperature in the cell
PartState(4:6,iPart) = CalcVelocity_maxwell_particle(iSpec,MacroRestartValues(iElem,iSpec,4:6)) &
                          + MacroRestartValues(iElem,iSpec,1:3)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(iSpec)%ChargeIC.GT.0.0) THEN
    IF (ALLOCATED(AmbipolElecVelo(iPart)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(iPart)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(iPart)%ElecVelo(3))
    AmbipolElecVelo(iPart)%ElecVelo(1:3) = CalcVelocity_maxwell_particle(DSMC%AmbiDiffElecSpec, &
          MacroRestartValues(iElem,DSMC%AmbiDiffElecSpec,4:6)) + MacroRestartValues(iElem,DSMC%AmbiDiffElecSpec,1:3)
  END IF
END IF
! 2) Set internal energies (rotational, vibrational, electronic)
IF(CollisMode.GT.1) THEN
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    PartStateIntEn(1,iPart) = CalcEVib_particle(iSpec,MacroRestartValues(iElem,iSpec,DSMC_TVIB),iPart)
    PartStateIntEn(2,iPart) = CalcERot_particle(iSpec,MacroRestartValues(iElem,iSpec,DSMC_TROT))
  ELSE
    PartStateIntEn(1:2,iPart) = 0.0
  END IF
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      PartStateIntEn(3,iPart) = CalcEElec_particle(iSpec,MacroRestartValues(iElem,iSpec,DSMC_TELEC),iPart)
    ELSE
      PartStateIntEn(3,iPart) = 0.0
    END IF
  END IF
END IF

! 3) Set the species and element number
PartSpecies(iPart) = iSpec
PEM%GlobalElemID(iPart) = iElem+offSetElem
PEM%LastGlobalElemID(iPart) = iElem+offSetElem
PDM%ParticleInside(iPart) = .TRUE.

! 4) Set particle weights (if required)
IF (VarTimeStep%UseVariableTimeStep) THEN
  VarTimeStep%ParticleTimeStep(iPart) = CalcVarTimeStep(PartState(1,iPart),PartState(2,iPart),iElem)
END IF
IF (RadialWeighting%DoRadialWeighting) THEN
  PartMPF(iPart) = CalcRadWeightMPF(PartState(2,iPart),iSpec,iPart)
END IF

END SUBROUTINE MacroRestart_InitializeParticle_Maxwell

END MODULE MOD_Macro_Restart
