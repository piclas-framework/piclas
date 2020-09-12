!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_macro_restart
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
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry           ,ONLY: CalcRadWeightMPF
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping, TriaTracking
USE MOD_Particle_Localization   ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Tools     ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemEpsOneCell
USE MOD_Particle_Vars           ,ONLY: Species, PDM, nSpecies, PartState, Symmetry, VarTimeStep
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,BoundsOfElem_Shared
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
REAL                                :: iRan, RandomPos(3), PartDens, TempMPF, MaxPosTemp, MinPosTemp
REAL                                :: TempVol, Volume, Det(6,2), RefPos(1:3)
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
              IF (DoRefMapping) THEN
                CALL GetPositionInRefElem(RandomPos,RefPos,GlobalElemID)
                IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
              ELSE
                IF (TriaTracking) THEN
                  CALL ParticleInsideQuad3D(RandomPos,GlobalElemID,InsideFlag,Det)
                ELSE
                  CALL PartInElemCheck(RandomPos,iPart,GlobalElemID,InsideFlag)
                END IF
              END IF
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
              IF (DoRefMapping) THEN
                CALL GetPositionInRefElem(RandomPos,RefPos,GlobalElemID)
                IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
              ELSE
                IF (TriaTracking) THEN
                  CALL ParticleInsideQuad3D(RandomPos,GlobalElemID,InsideFlag,Det)
                ELSE
                  CALL PartInElemCheck(RandomPos,iPart,GlobalElemID,InsideFlag)
                END IF
              END IF
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
          IF (DoRefMapping) THEN
            CALL GetPositionInRefElem(RandomPos,RefPos,GlobalElemID)
            IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
          ELSE
            IF (TriaTracking) THEN
              CALL ParticleInsideQuad3D(RandomPos,GlobalElemID,InsideFlag,Det)
            ELSE
              CALL PartInElemCheck(RandomPos,iPart,GlobalElemID,InsideFlag)
            END IF
          END IF
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
          IF (DoRefMapping) THEN
            CALL GetPositionInRefElem(RandomPos,RefPos,iElem)
            IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
          ELSE
            IF (TriaTracking) THEN
              CALL ParticleInsideQuad3D(RandomPos,iElem,InsideFlag,Det)
            ELSE
              CALL PartInElemCheck(RandomPos,iPart,iElem,InsideFlag)
            END IF
          END IF
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
          IF (DoRefMapping) THEN
            CALL GetPositionInRefElem(RandomPos,RefPos,GlobalElemID)
            IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
          ELSE
            IF (TriaTracking) THEN
              CALL ParticleInsideQuad3D(RandomPos,GlobalElemID,InsideFlag,Det)
            ELSE
              CALL PartInElemCheck(RandomPos,iPart,GlobalElemID,InsideFlag)
            END IF
          END IF
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
USE MOD_Particle_Vars           ,ONLY: PDM, PartSpecies, PartState, PEM, VarTimeStep, PartMPF
USE MOD_DSMC_Vars               ,ONLY: DSMC, PartStateIntEn, CollisMode, SpecDSMC, RadialWeighting
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_DSMC_Symmetry           ,ONLY: CalcRadWeightMPF
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

! 2) Set internal energies (rotational, vibrational, electronic)
IF(CollisMode.GT.1) THEN
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    PartStateIntEn(1,iPart) = CalcEVib_particle(iSpec,iPart,MacroRestartValues(iElem,iSpec,DSMC_TVIB))
    PartStateIntEn(2,iPart) = CalcERot_particle(iSpec,MacroRestartValues(iElem,iSpec,DSMC_TROT))
  ELSE
    PartStateIntEn(1:2,iPart) = 0.0
  END IF
  IF(DSMC%ElectronicModel) THEN
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      PartStateIntEn(3,iPart) = CalcEElec_particle(iSpec,MacroRestartValues(iElem,iSpec,DSMC_TELEC))
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


FUNCTION CalcVelocity_maxwell_particle(iSpec,Temp)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species
USE Ziggurat,                   ONLY : rnor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec
REAL, INTENT(IN)                :: Temp(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: CalcVelocity_maxwell_particle(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CalcVelocity_maxwell_particle(1:3) = 0.0

IF(Temp(1).GT.0.0) CalcVelocity_maxwell_particle(1) = rnor()*SQRT(BoltzmannConst*Temp(1)/Species(iSpec)%MassIC)
IF(Temp(2).GT.0.0) CalcVelocity_maxwell_particle(2) = rnor()*SQRT(BoltzmannConst*Temp(2)/Species(iSpec)%MassIC)
IF(Temp(3).GT.0.0) CalcVelocity_maxwell_particle(3) = rnor()*SQRT(BoltzmannConst*Temp(3)/Species(iSpec)%MassIC)

END FUNCTION CalcVelocity_maxwell_particle


REAL FUNCTION CalcEVib_particle(iSpec,iPart,TempVib)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC, PolyatomMolDSMC, VibQuantsPar, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec, iPart
REAL, INTENT(IN)          :: TempVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: iRan
INTEGER                   :: iQuant, iDOF, iPolyatMole
!===================================================================================================================================

IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
  ! set vibrational energy
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  CalcEVib_particle = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TempVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
    DO WHILE (iQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*TempVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
    END DO
    CalcEVib_particle = CalcEVib_particle &
                                + (iQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
    VibQuantsPar(iPart)%Quants(iDOF)=iQuant
  END DO
ELSE
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT(-LOG(iRan)*TempVib/SpecDSMC(iSpec)%CharaTVib)
  DO WHILE (iQuant.GE.SpecDSMC(iSpec)%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TempVib/SpecDSMC(iSpec)%CharaTVib)
  END DO
  CalcEVib_particle = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
END IF

RETURN

END FUNCTION CalcEVib_particle


REAL FUNCTION CalcERot_particle(iSpec,TempRot)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec
REAL, INTENT(IN)          :: TempRot
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: PartStateTempVar, NormProb, iRan2
!===================================================================================================================================

IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
  CALL RANDOM_NUMBER(iRan2)
  CalcERot_particle = -BoltzmannConst*TempRot*LOG(iRan2)
ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateTempVar = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
  NormProb = SQRT(PartStateTempVar)*EXP(-PartStateTempVar)/(SQRT(0.5)*EXP(-0.5))
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan2)
    PartStateTempVar = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateTempVar)*EXP(-PartStateTempVar)/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
  END DO
  CalcERot_particle = PartStateTempVar*BoltzmannConst*TempRot
END IF

RETURN

END FUNCTION CalcERot_particle


REAL FUNCTION CalcEElec_particle(iSpec,TempElec)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec
REAL, INTENT(IN)          :: TempElec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iQua
REAL                      :: iRan, ElectronicPartition, ElectronicPartitionTemp
!===================================================================================================================================

ElectronicPartition  = 0.
ElectronicPartitionTemp = 0.
IF(TempElec.GT.0.0) THEN
  ! calculate sum over all energy levels == partition function for temperature Telec
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TempElec)
    IF ( ElectronicPartitionTemp .GT. ElectronicPartition ) THEN
      ElectronicPartition = ElectronicPartitionTemp
    END IF
  END DO
  ElectronicPartitionTemp = 0.
  ! select level
  CALL RANDOM_NUMBER(iRan)
  DO WHILE ( iRan .GE. ElectronicPartitionTemp / ElectronicPartition )
    CALL RANDOM_NUMBER(iRan)
    iQua = int( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
    ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TempElec)
    CALL RANDOM_NUMBER(iRan)
  END DO
ELSE
  iQua = 0
END IF
CalcEElec_particle = BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)

RETURN

END FUNCTION CalcEElec_particle


END MODULE MOD_macro_restart
