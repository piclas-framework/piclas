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
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF,InitializeParticleMaxwell, IncreaseMaxParticleNumber
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_Particle_Vars           ,ONLY: Species, PDM, nSpecies, PartState, Symmetry, UseVarTimeStep
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,BoundsOfElem_Shared, ElemMidPoint_Shared
USE MOD_Particle_Tracking       ,ONLY: ParticleInsideCheck
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
INTEGER                             :: iElem,iSpec,iPart,nPart,locnPart,iHeight,yPartitions,CNElemID,GlobalElemID
REAL                                :: iRan, RandomPos(3), PartDens, MaxPosTemp, MinPosTemp
REAL                                :: TempVol, Volume
LOGICAL                             :: InsideFlag
!===================================================================================================================================

SWRITE(UNIT_stdOut,*) 'PERFORMING MACROSCOPIC RESTART...'

locnPart = 1

DO iElem = 1, nElems
  GlobalElemID = iElem + offsetElem
  CNElemID = GetCNElemID(GlobalElemID)
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GlobalElemID) ) ! 1-2: Min, Max value; 1-3: x,y,z
! #################### 2D ##########################################################################################################
    IF (Symmetry%Axisymmetric) THEN
      IF (RadialWeighting%DoRadialWeighting) THEN
        DO iSpec = 1, nSpecies
          IF (DSMC%DoAmbipolarDiff) THEN
            IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
          END IF
          yPartitions = 6
          ! Particle weighting
          DO iHeight = 1, yPartitions
            MinPosTemp = Bounds(1,2) + (Bounds(2,2) - Bounds(1,2))/ yPartitions *(iHeight-1.)
            MaxPosTemp = Bounds(1,2) + (Bounds(2,2) - Bounds(1,2))/ yPartitions *iHeight
            TempVol =  (MaxPosTemp-MinPosTemp)*(Bounds(2,1)-Bounds(1,1)) * Pi * (MaxPosTemp+MinPosTemp)
            PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / CalcRadWeightMPF((MaxPosTemp+MinPosTemp)*0.5,iSpec)
            IF(UseVarTimeStep) THEN
              PartDens = PartDens * GetParticleTimeStep((Bounds(2,1)+Bounds(1,1))*0.5, (MaxPosTemp+MinPosTemp)*0.5, iElem)
            END IF
            CALL RANDOM_NUMBER(iRan)
            nPart = INT(PartDens  * TempVol + iRan)
            DO iPart = 1, nPart
              InsideFlag=.FALSE.
              CALL RANDOM_NUMBER(RandomPos)
              RandomPos(1) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
              RandomPos(2) = MinPosTemp + RandomPos(2)*(MaxPosTemp-MinPosTemp)
              RandomPos(3) = 0.0
              InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
              IF (InsideFlag) THEN
                IF (locnPart.GE.PDM%maxParticleNumber) CALL IncreaseMaxParticleNumber()
                PartState(1:3,locnPart) = RandomPos(1:3)
                CALL InitializeParticleMaxwell(locnPart,iSpec,iElem,Mode=1)
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
          PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / Species(iSpec)%MacroParticleFactor
          CALL RANDOM_NUMBER(iRan)
          IF(UseVarTimeStep) THEN
            PartDens = PartDens / GetParticleTimeStep(ElemMidPoint_Shared(1,CNElemID), ElemMidPoint_Shared(2,CNElemID), iElem)
          END IF
          nPart = INT(PartDens * ElemVolume_Shared(CNElemID) + iRan)
          CALL IncreaseMaxParticleNumber(nPart)
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
            CALL InitializeParticleMaxwell(locnPart,iSpec,iElem,Mode=1)
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
        PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / Species(iSpec)%MacroParticleFactor
        CALL RANDOM_NUMBER(iRan)
        IF(UseVarTimeStep) THEN
          PartDens = PartDens / GetParticleTimeStep(ElemMidPoint_Shared(1,CNElemID), ElemMidPoint_Shared(2,CNElemID), iElem)
        END IF
        nPart = INT(PartDens * ElemVolume_Shared(CNElemID) + iRan)
        CALL IncreaseMaxParticleNumber(nPart)

        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          DO WHILE(.NOT.InsideFlag)
            CALL RANDOM_NUMBER(RandomPos(1:2))
            RandomPos(1:2) = Bounds(1,1:2) + RandomPos(1:2)*(Bounds(2,1:2)-Bounds(1,1:2))
            RandomPos(3) = 0.0
            InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          END DO
          PartState(1:3,locnPart) = RandomPos(1:3)
          CALL InitializeParticleMaxwell(locnPart,iSpec,iElem,Mode=1)
          locnPart = locnPart + 1
        END DO ! nPart
      END DO ! nSpecies
    ELSE IF(Symmetry%Order.EQ.1) THEN
      Volume = (Bounds(2,1) - Bounds(1,1))
      DO iSpec = 1, nSpecies
        IF (DSMC%DoAmbipolarDiff) THEN
          IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
        END IF
        PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / Species(iSpec)%MacroParticleFactor
        CALL RANDOM_NUMBER(iRan)
        IF(UseVarTimeStep) THEN
          PartDens = PartDens / GetParticleTimeStep(ElemMidPoint_Shared(1,CNElemID), ElemMidPoint_Shared(2,CNElemID), iElem)
        END IF
        nPart = INT(PartDens * ElemVolume_Shared(CNElemID) + iRan)
        CALL IncreaseMaxParticleNumber(nPart)
        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          DO WHILE(.NOT.InsideFlag)
            CALL RANDOM_NUMBER(RandomPos(1))
            RandomPos(1:2) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
            RandomPos(2) = 0.0
            RandomPos(3) = 0.0
            InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          END DO
          PartState(1:3,locnPart) = RandomPos(1:3)
          CALL InitializeParticleMaxwell(locnPart,iSpec,iElem,Mode=1)
          locnPart = locnPart + 1
        END DO ! nPart
      END DO ! nSpecies
    ELSE
! #################### 3D ##########################################################################################################
      Volume = (Bounds(2,3) - Bounds(1,3))*(Bounds(2,2) - Bounds(1,2))*(Bounds(2,1) - Bounds(1,1))
      DO iSpec = 1, nSpecies
        IF (DSMC%DoAmbipolarDiff) THEN
          IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
        END IF
        PartDens = MacroRestartValues(iElem,iSpec,DSMC_NUMDENS) / Species(iSpec)%MacroParticleFactor
        CALL RANDOM_NUMBER(iRan)
        IF(UseVarTimeStep) THEN
          PartDens = PartDens / GetParticleTimeStep(ElemMidPoint_Shared(1,CNElemID), ElemMidPoint_Shared(2,CNElemID), iElem)
        END IF
        nPart = INT(PartDens * ElemVolume_Shared(CNElemID) + iRan)
        CALL IncreaseMaxParticleNumber(nPart)
        DO iPart = 1, nPart
          InsideFlag=.FALSE.
          DO WHILE(.NOT.InsideFlag)
            CALL RANDOM_NUMBER(RandomPos)
            RandomPos(1:3) = Bounds(1,1:3) + RandomPos(1:3)*(Bounds(2,1:3)-Bounds(1,1:3))
            InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
          END DO
          PartState(1:3,locnPart) = RandomPos(1:3)
          CALL InitializeParticleMaxwell(locnPart,iSpec,iElem,Mode=1)
          locnPart = locnPart + 1
        END DO ! nPart
      END DO ! nSpecies
    END IF ! 1D/2D/Axisymmetric/3D
  END ASSOCIATE
END DO ! nElems

PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
#ifdef CODE_ANALYZE
IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
  IF (PDM%ParticleInside(iPart)) THEN
    IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
    CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
  END IF
END DO
#endif

END SUBROUTINE MacroRestart_InsertParticles


END MODULE MOD_Macro_Restart
