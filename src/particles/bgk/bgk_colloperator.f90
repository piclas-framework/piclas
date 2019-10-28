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

MODULE MOD_BGK_CollOperator
!===================================================================================================================================
! Module approximating the collision operator using the Bhatnagar-Gross-Krook model
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGK_CollisionOperator
  MODULE PROCEDURE BGK_CollisionOperator
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGK_CollisionOperator, ARShakhov, CalcTEquiPoly, CalcTEqui
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_CollisionOperator(iPartIndx_Node, nPart, NodeVolume, vBulkAll, AveragingPara, CorrectStep)
!===================================================================================================================================
!> Subroutine for the cell-local BGK collision operator:
!> 1.) Moment calculation: Summing up the relative velocities and their squares
!> 2.) Calculation of the relaxation frequency of the distribution function towards the target distribution function
!> 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency
!> 4.) Determine the number of particles undergoing a relaxation (including vibration and rotation)
!> 5.) Sample new particle velocities from the target distribution function, depending on the chosen model
!> 6.) Determine the new bulk velocity and the new relative velocity of the particles
!> 7.) Treatment of the vibrational energy of molecules
!> 8.) Determine the new DSMC_RHS (for molecules, including rotational energy)
!> 9.) Scaling of the rotational energy of molecules
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, VarTimeStep, usevMPF
USE MOD_DSMC_Vars             ,ONLY: DSMC_RHS, SpecDSMC, DSMC, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar, RadialWeighting
USE MOD_DSMC_Vars             ,ONLY: CollInf
USE MOD_DSMC_Analyze          ,ONLY: CalcTVibPoly
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Globals_Vars          ,ONLY: Pi, BoltzmannConst
USE MOD_BGK_Vars              ,ONLY: SpecBGK, ESBGKModel, BGKCollModel, BGKUnifiedCes, BGKMovingAverageLength, BGKMovingAverage
USE MOD_BGK_Vars              ,ONLY: BGKUseQuantVibEn, BGKDoVibRelaxation, SBGKEnergyConsMethod
USE MOD_BGK_Vars              ,ONLY: BGK_MeanRelaxFactor, BGK_MeanRelaxFactorCounter, BGK_MaxRelaxFactor, BGK_MaxRotRelaxFactor
USE MOD_part_tools            ,ONLY: GetParticleWeight
#ifdef CODE_ANALYZE
USE MOD_Globals               ,ONLY: abort,unit_stdout,myrank
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(IN)                        :: vBulkAll(3)
REAL, INTENT(INOUT), OPTIONAL           :: AveragingPara(5,BGKMovingAverageLength)
INTEGER, INTENT(INOUT), OPTIONAL        :: CorrectStep
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: KronDelta, tempVelo(3), vBulk(3), u0ij(3,3), SMat(3,3), u2, V_rel(3), vmag2, u0i(3), u2i(3)
REAL                  :: alpha, CellTemp, dens, InnerDOF, dynamicvis, iRan, NewEn, OldEn, Prandtl, relaxfreq
REAL                  :: rotrelaxfreq, vibrelaxfreq, collisionfreq, ProbAddPart, Evib, Tvib, Xi_vib, TEqui, Xi_Vib_old, Xi_rot
REAL                  :: MaxColQua, ERot    !, TEquiV, TEquiR
INTEGER               :: iLoop, nRelax, fillMa1, fillMa2, iQuant, iQuaMax, iDOF, iPolyatMole
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelax(:),iPartIndx_NodeRelaxTemp(:),iPartIndx_NodeRelaxRot(:),iPartIndx_NodeRelaxVib(:)
REAL, ALLOCATABLE     :: iRanPart(:,:), Xi_vib_DOF(:), VibEnergyDOF(:,:)
REAL                  :: A(3,3), Work(1000), W(3), trace, CShak
INTEGER               :: INFO, nNotRelax, nRotRelax, nVibRelax
REAL                  :: TRot, betaV, OldEnRot, RotExp, VibExp, NewEnRot, NewEnVib, vBulkRelaxOld(3),vBulkRelax(3)
REAL                  :: CellTempRelax, vBulkAver(3), u2Aver, nPartAver, dtCell
REAL                  :: partWeight, totalWeight, totalWeightRelax
#ifdef CODE_ANALYZE
REAL                  :: Energy_old,Energy_new,Momentum_old(3),Momentum_new(3)
INTEGER               :: iMom
#endif /* CODE_ANALYZE */
!===================================================================================================================================
#ifdef CODE_ANALYZE
! Momentum and energy conservation check: summing up old values
Momentum_new = 0.0; Momentum_old = 0.0; Energy_new = 0.0; Energy_old = 0.0
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  Momentum_old(1:3) = Momentum_old(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
  Energy_old = Energy_old + (PartState(iPartIndx_Node(iLoop),4)**2. + PartState(iPartIndx_Node(iLoop),5)**2. &
           + PartState(iPartIndx_Node(iLoop),6)**2.)*0.5*Species(1)%MassIC * partWeight
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    Energy_old = Energy_old + (PartStateIntEn(iPartIndx_Node(iLoop),1) + PartStateIntEn(iPartIndx_Node(iLoop),2)) * partWeight
  END IF
END DO
#endif

NewEn = 0.; OldEn = 0.
OldEnRot = 0.; NewEnRot = 0.; NewEnVib = 0.
u2 = 0.0; u0ij = 0.0; u0i = 0.0; u2i = 0.0
Evib = 0.0; ERot = 0.0
u2Aver = 0.0; vBulkRelax = 0.0; vBulkRelaxOld = 0.0
totalWeight = 0.0; dtCell = 0.0; totalWeightRelax = 0.0

! 1.) Summing up the relative velocities and their square to calculate the moments
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulkAll(1:3)
  IF (BGKMovingAverage) u2Aver = u2Aver + PartState(iPartIndx_Node(iLoop),4)**2. &
      +PartState(iPartIndx_Node(iLoop),5)**2. + PartState(iPartIndx_Node(iLoop),6)**2.
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  u2= u2 + vmag2 * partWeight
  IF ((BGKCollModel.EQ.1).OR.(BGKCollModel.EQ.4)) THEN
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ij(fillMa1, fillMa2)= u0ij(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2) * partWeight
      END DO
    END DO
    u0i(1:3) = u0i(1:3) + V_rel(1:3) * partWeight
  END IF
  IF ((BGKCollModel.EQ.2).OR.(BGKCollModel.EQ.4)) u2i(1:3) = u2i(1:3) + V_rel(1:3)*vmag2 * partWeight
  IF (SBGKEnergyConsMethod.NE.2) OldEn = OldEn + 0.5*Species(1)%MassIC * vmag2 * partWeight
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    IF(BGKDoVibRelaxation) Evib = Evib + (PartStateIntEn(iPartIndx_Node(iLoop),1) - SpecDSMC(1)%EZeroPoint) * partWeight
    ERot = ERot + PartStateIntEn(iPartIndx_Node(iLoop),2) * partWeight
  END IF
  IF(VarTimeStep%UseVariableTimeStep) THEN
    dtCell = dtCell + VarTimeStep%ParticleTimeStep(iPartIndx_Node(iLoop))
  END IF
END DO

IF(VarTimeStep%UseVariableTimeStep) THEN
  dtCell = dt * dtCell / nPart
ELSE
  dtCell = dt
END IF

IF ((BGKCollModel.EQ.2).OR.(BGKCollModel.EQ.4)) THEN
  IF(nPart.LE.2) RETURN
  u2i = u2i*nPart*nPart/((nPart-1.)*(nPart-2.)*totalWeight)
END IF

u2 = u2*nPart/((nPart-1.)*totalWeight)
u0ij = u0ij/totalWeight
u0i(1:3) = u0i(1:3)/totalWeight
trace = u0ij(1,1)+u0ij(2,2)+u0ij(3,3)
CellTemp = Species(1)%MassIC * u2 / (3.0*BoltzmannConst)

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  dens = totalWeight / NodeVolume
ELSE
  dens = totalWeight * Species(1)%MacroParticleFactor / NodeVolume
END IF

! Calculation of the rotational and vibrational degrees of freedom for molecules
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  Xi_vib = 0.0
  IF(BGKDoVibRelaxation) THEN
    IF(SpecDSMC(1)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(1)%SpecToPolyArray
      ALLOCATE(Xi_vib_DOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
      Xi_vib_DOF(:) = 0.
      TVib = CalcTVibPoly(Evib/totalWeight, 1)
      IF (TVib.GT.0.0) THEN
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          Xi_vib = Xi_vib + 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVib &
                              /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVib) - 1.)
        END DO
      END IF
    ELSE
      TVib=Evib / (totalWeight*BoltzmannConst*SpecDSMC(1)%CharaTVib)
      IF (TVib.GT.0.0) THEN
        Tvib= SpecDSMC(1)%CharaTVib/LOG(1. + 1./(TVib))
        Xi_vib = 2.* Evib / (totalWeight*BoltzmannConst*Tvib)
      END IF
    END IF
    Xi_Vib_old = Xi_Vib
  ELSE
    Xi_Vib_old = 0.0
    TVib = 0.0
  END IF
  Xi_rot = SpecDSMC(1)%Xi_Rot
  TRot = 2.*ERot/(Xi_rot*totalWeight*BoltzmannConst)
  InnerDOF = Xi_rot + Xi_Vib
ELSE
  InnerDOF = 0.
  Xi_rot = 0.
  Xi_vib = 0.
END IF

! Storing the relevant variables for the moving average. In the case of a limited moving average (BGKMovingAverageLength.GT.1),
! the most recent value overwrites the first entry.
IF (BGKMovingAverage) THEN
  IF (BGKMovingAverageLength.GT.1) THEN
    CorrectStep = CorrectStep + 1
    IF (CorrectStep.GT.BGKMovingAverageLength) CorrectStep = 1
    AveragingPara(1:3,CorrectStep) = vBulkAll(1:3)*nPart
    AveragingPara(4,CorrectStep) = u2Aver
    AveragingPara(5,CorrectStep) = nPart
    vBulkAver(1) = SUM(AveragingPara(1,:))
    vBulkAver(2) = SUM(AveragingPara(2,:))
    vBulkAver(3) = SUM(AveragingPara(3,:))
    u2Aver = SUM(AveragingPara(4,:))
    nPartAver = SUM(AveragingPara(5,:))
    CellTempRelax = Species(1)%MassIC / (3.*BoltzmannConst) &
                  * 1./(nPartAver-1.)*(u2Aver - (vBulkAver(1)**2. &
                  + vBulkAver(2)**2. + vBulkAver(3)**2.)/nPartAver)
  ELSE
    AveragingPara(1:3,1) = AveragingPara(1:3,1) + vBulkAll(1:3)*nPart
    AveragingPara(4,1) = AveragingPara(4,1) + u2Aver
    AveragingPara(5,1) = AveragingPara(5,1) + nPart
    CellTempRelax = Species(1)%MassIC / (3.*BoltzmannConst) &
                  * 1./(AveragingPara(5,1)-1.)*(AveragingPara(4,1) - (AveragingPara(1,1)**2. &
                  + AveragingPara(2,1)**2. + AveragingPara(3,1)**2.)/AveragingPara(5,1))
  END IF
ELSE
  CellTempRelax = CellTemp
END IF

! 2.) Calculate the reference dynamic viscosity, Prandtl number and the resulting relaxation frequency of the distribution function
dynamicvis = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*CollInf%Tref(1,1)/Pi) &
           / (4.*(4.- 2.*CollInf%omegaLaux(1,1)) * (6. - 2.*CollInf%omegaLaux(1,1))* CollInf%dref(1,1)**2.)
Prandtl =2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
CShak= Prandtl*(1.-BGKUnifiedCes)
IF (BGKCollModel.EQ.1) THEN
  relaxfreq = Prandtl*dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omegaLaux(1,1) + 0.5) &
      /dynamicvis*CellTempRelax**(-CollInf%omegaLaux(1,1) +0.5)
ELSE IF (BGKCollModel.EQ.4) THEN
  relaxfreq = dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omegaLaux(1,1) + 0.5) &
      /(dynamicvis*(1.-BGKUnifiedCes))*CellTempRelax**(-CollInf%omegaLaux(1,1) +0.5)
ELSE
  relaxfreq = dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omegaLaux(1,1) + 0.5) &
      /dynamicvis*CellTempRelax**(-CollInf%omegaLaux(1,1) +0.5)
END IF

IF(DSMC%CalcQualityFactors) THEN
  BGK_MeanRelaxFactor         = BGK_MeanRelaxFactor + relaxfreq * dtCell
  BGK_MeanRelaxFactorCounter  = BGK_MeanRelaxFactorCounter + 1
  BGK_MaxRelaxFactor          = MAX(BGK_MaxRelaxFactor,relaxfreq*dtCell)
END IF

! 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency using the collision frequency,
!     which is not the same as the relaxation frequency of distribution function, calculated above.
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  collisionfreq = SpecBGK(1)%CollFreqPreFactor(1) * Dens *CellTempRelax**(-CollInf%omegaLaux(1,1) +0.5)
  rotrelaxfreq = collisionfreq * DSMC%RotRelaxProb
  vibrelaxfreq = collisionfreq * DSMC%VibRelaxProb
  IF(SpecDSMC(1)%PolyatomicMol) THEN
    CALL CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_vib_DOF, Xi_Vib_old, RotExp, VibExp, &
                        TEqui, rotrelaxfreq, vibrelaxfreq, dtCell)
    Xi_vib = SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ELSE
    CALL CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp,  &
                    TEqui, rotrelaxfreq, vibrelaxfreq, dtCell)
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    BGK_MaxRotRelaxFactor          = MAX(BGK_MaxRotRelaxFactor,rotrelaxfreq*dtCell)
  END IF
END IF

vBulk(1:3) = 0.0; nRelax = 0; nNotRelax = 0; nRotRelax = 0; nVibRelax = 0
ALLOCATE(iPartIndx_NodeRelax(nPart), iPartIndx_NodeRelaxTemp(nPart))
iPartIndx_NodeRelaxTemp = 0

! 4.) Determine the number of particles undergoing a relaxation (including rotational and vibrational relaxation for molecules)
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  ALLOCATE(iPartIndx_NodeRelaxRot(nPart),iPartIndx_NodeRelaxVib(nPart))
  DO iLoop = 1, nPart
    partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-exp(-relaxfreq*dtCell)
    IF (ProbAddPart.GT.iRan) THEN
      nRelax = nRelax + 1
      iPartIndx_NodeRelax(nRelax) = iPartIndx_Node(iLoop)
    ELSE
      nNotRelax = nNotRelax + 1
      iPartIndx_NodeRelaxTemp(nNotRelax) = iPartIndx_Node(iLoop)
      vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
    END IF
    !Rotation
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-RotExp
    IF (ProbAddPart.GT.iRan) THEN
      nRotRelax = nRotRelax + 1
      iPartIndx_NodeRelaxRot(nRotRelax) = iPartIndx_Node(iLoop)
      OldEnRot = OldEnRot + PartStateIntEn(iPartIndx_Node(iLoop),2) * partWeight
    END IF
    ! Vibration
    IF(BGKDoVibRelaxation) THEN
      CALL RANDOM_NUMBER(iRan)
      ProbAddPart = 1.-VibExp
      IF (ProbAddPart.GT.iRan) THEN
        nVibRelax = nVibRelax + 1
        iPartIndx_NodeRelaxVib(nVibRelax) = iPartIndx_Node(iLoop)
        OldEn = OldEn + (PartStateIntEn(iPartIndx_NodeRelaxVib(nVibRelax),1) - SpecDSMC(1)%EZeroPoint) * partWeight
      END IF
    END IF
  END DO
  IF ((nRelax.EQ.0).AND.(nRotRelax.EQ.0).AND.(nVibRelax.EQ.0)) RETURN
    !! VIB RElaxation
    IF(BGKDoVibRelaxation) THEN
      IF(SpecDSMC(1)%PolyatomicMol) THEN
        ALLOCATE(VibEnergyDOF(nVibRelax,PolyatomMolDSMC(iPolyatMole)%VibDOF))
        DO iLoop = 1, nVibRelax
          partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1) = 0.0
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            CALL RANDOM_NUMBER(iRan)
            VibEnergyDOF(iLoop,iDOF) = - LOG(iRan)*Xi_vib_DOF(iDOF)/2.*TEqui*BoltzmannConst
            PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
                                                                + VibEnergyDOF(iLoop,iDOF)
          END DO
          NewEnVib = NewEnVib + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1) * partWeight
        END DO
      ELSE
        DO iLoop = 1, nVibRelax
          partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
          CALL RANDOM_NUMBER(iRan)
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = -LOG(iRan)*Xi_vib/2.*TEqui*BoltzmannConst
          NewEnVib = NewEnVib + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1) * partWeight
        END DO
      END IF
    END IF
    !! ROT Relaxation
    DO iLoop = 1, nRotRelax
      partWeight = GetParticleWeight(iPartIndx_NodeRelaxRot(iLoop))
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = -Xi_Rot / 2. * BoltzmannConst*TEqui*LOG(iRan)
      NewEnRot = NewEnRot + PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) * partWeight
    END DO
ELSE ! Atoms
  DO iLoop = 1, nPart
    partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-exp(-relaxfreq*dtCell)
    IF (ProbAddPart.GT.iRan) THEN
      nRelax = nRelax + 1
      iPartIndx_NodeRelax(nRelax) = iPartIndx_Node(iLoop)
      IF (SBGKEnergyConsMethod.EQ.2) THEN
        vBulkRelaxOld(1:3) = vBulkRelaxOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
        totalWeightRelax =  totalWeightRelax + partWeight
      END IF
    ELSE
      nNotRelax = nNotRelax + 1
      iPartIndx_NodeRelaxTemp(nNotRelax) = iPartIndx_Node(iLoop)
      vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
    END IF
  END DO
  IF (nRelax.EQ.0) RETURN
  IF (SBGKEnergyConsMethod.EQ.2) THEN
    vBulkRelaxOld = vBulkRelaxOld / totalWeightRelax
    IF (nRelax.GT.2) THEN
      DO iLoop = 1, nRelax
        partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
        OldEn = OldEn + 0.5*Species(1)%MassIC*((PartState(iPartIndx_NodeRelax(iLoop),4)-vBulkRelaxOld(1))**2.0 &
                                             + (PartState(iPartIndx_NodeRelax(iLoop),5)-vBulkRelaxOld(2))**2.0 &
                                             + (PartState(iPartIndx_NodeRelax(iLoop),6)-vBulkRelaxOld(3))**2.0) * partWeight
      END DO
    ELSE
      DO iLoop = 1, nPart
        partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
        OldEn = OldEn + 0.5*Species(1)%MassIC*((PartState(iPartIndx_Node(iLoop),4)-vBulkAll(1))**2.0 &
                                             + (PartState(iPartIndx_Node(iLoop),5)-vBulkAll(2))**2.0 &
                                             + (PartState(iPartIndx_Node(iLoop),6)-vBulkAll(3))**2.0) * partWeight
      END DO
    END IF
  END IF
END IF ! (SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)

totalWeightRelax = 0.0

! 5.) Sample new particle velocities from the target distribution function, depending on the chosen model
IF (nRelax.GT.0) THEN
  ALLOCATE(iRanPart(3,nRelax))
  SELECT CASE(BGKCollModel)
  CASE (1)  ! Ellipsoidal Statistical
    IF (ESBGKModel.EQ.1) THEN
      !! Approximated Solution
      DO fillMa1 =1, 3
        DO fillMa2 =fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandtl)/(2.*Prandtl) &
            *(Species(1)%MassIC/(BoltzmannConst*CellTemp)*nPart/(nPart-1.) &
            *(u0ij(fillMa1, fillMa2)-u0i(fillMa1)*u0i(fillMa2))-KronDelta)
!          SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandtl)/(2.*Prandtl) &
!            *(Species(1)%MassIC/(BoltzmannConst*CellTemp) &
!            *(u0ij(fillMa1, fillMa2)-u0i(fillMa1)*u0i(fillMa2))-KronDelta)
        END DO
      END DO
      SMat(2,1)=SMat(1,2)
      SMat(3,1)=SMat(1,3)
      SMat(3,2)=SMat(2,3)
      CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
    ELSE
      !! Exact Solution
      DO fillMa1 =1, 3
        DO fillMa2 =fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          A(fillMa1, fillMa2) = KronDelta - (1.-Prandtl)/Prandtl*(3.*u0ij(fillMa1, fillMa2)/trace - KronDelta)
        END DO
      END DO
      IF (ESBGKModel.EQ.2) THEN
        CALL DSYEV('V','U',3,A,3,W,Work,1000,INFO)
        SMat = 0.0
        IF (W(1).LT.0.0) THEN
          W(1) = 0.0
          IF (W(2).LT.0) W(2) = 0.0
        END IF
        IF (W(3).LT.0) THEN
          W(3) = 0.0
          DO fillMa1 =1, 3
            DO fillMa2 =fillMa1, 3
              IF (fillMa1.EQ.fillMa2) THEN
                KronDelta = 1.0
              ELSE
                KronDelta = 0.0
              END IF
              SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandtl)/(2.*Prandtl) &
                *(Species(1)%MassIC/(BoltzmannConst*CellTemp)*nPart/(nPart-1.) &
                *(u0ij(fillMa1, fillMa2)-u0i(fillMa1)*u0i(fillMa2))-KronDelta)
            END DO
          END DO
          SMat(2,1)=SMat(1,2)
          SMat(3,1)=SMat(1,3)
          SMat(3,2)=SMat(2,3)
        ELSE
          SMat(1,1) = SQRT(W(1))
          SMat(2,2) = SQRT(W(2))
          SMat(3,3) = SQRT(W(3))
          SMat = MATMUL(A, SMat)
          SMat = MATMUL(SMat, TRANSPOSE(A))
        END IF
        CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
      ELSE IF (ESBGKModel.EQ.3) THEN
        A(2,1)=A(1,2)
        A(3,1)=A(1,3)
        A(3,2)=A(2,3)
        CALL MetropolisES(nRelax, iRanPart, A)
      END IF
    END IF
  CASE (2)  ! Shakov
!    CALL MetropolisShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
    CALL ARShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
  CASE (3)  ! Standard BGK (Maxwell target distribution)
    CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
  CASE (4)  ! Unified BGK
      DO fillMa1 =1, 3
        DO fillMa2 =fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          A(fillMa1, fillMa2) = KronDelta + BGKUnifiedCes*(3.*u0ij(fillMa1, fillMa2)/trace - KronDelta)
        END DO
      END DO
      A(2,1)=A(1,2)
      A(3,1)=A(1,3)
      A(3,2)=A(2,3)
      CALL MetropolisUnified(nRelax, iRanPart,A, u2/3.,  u2i, CShak)
  END SELECT
  DO iLoop = 1, nRelax
    IF ((BGKCollModel.EQ.1).AND.(ESBGKModel.NE.3)) THEN
      tempVelo(1:3) = SQRT(BoltzmannConst*CellTemp/Species(1)%MassIC)*iRanPart(1:3,iLoop)
      DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkAll(1:3) + MATMUL(SMat,tempVelo)
    ELSE
      IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
        DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkRelaxOld(1:3) &
          + SQRT(BoltzmannConst*CellTemp/Species(1)%MassIC)*iRanPart(1:3,iLoop)
      ELSE
        DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkAll(1:3) &
          + SQRT(BoltzmannConst*CellTemp/Species(1)%MassIC)*iRanPart(1:3,iLoop)
      END IF
    END IF
    partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
    vBulkRelax(1:3) = vBulkRelax(1:3) + DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) * partWeight
    totalWeightRelax = totalWeightRelax + partWeight
  END DO
END IF ! nRelax.GT.0

! 6.) Determine the new bulk velocity and the new relative velocity of the particles that underwent relaxation
IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  vBulkRelax = vBulkRelax / totalWeightRelax
ELSE
  vBulkRelax = (vBulkRelax + vBulk) / totalWeight
  vBulk = vBulkRelax
END IF

IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  DO iLoop = 1, nRelax
    partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
    V_rel(1:3) = DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) - vBulkRelax(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC*partWeight
  END DO
ELSE
  DO iLoop = 1, nRelax
    partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
    V_rel(1:3) = DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) - vBulk(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC*partWeight
  END DO
  DO iLoop = 1, nPart-nRelax
    partWeight = GetParticleWeight(iPartIndx_NodeRelaxTemp(iLoop))
    V_rel(1:3) = PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6) - vBulk(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC*partWeight
  END DO
END IF

! 7.) Vibrational energy of the molecules: Determine the new state (either quantized or continuous) and ensure energy conservation
!     by scaling the new vibrational states with the factor alpha
IF(BGKDoVibRelaxation) THEN
  IF ((NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)) THEN
    IF (BGKUseQuantVibEn) THEN
      alpha = OldEn/NewEnVib*(Xi_Vib*nVibRelax/(3.*(nPart-1.)+Xi_Vib*nVibRelax))
      IF(SpecDSMC(1)%PolyatomicMol) THEN
        DO iLoop = 1, nVibRelax
          partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = 0.0
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            betaV = alpha*VibEnergyDOF(iLoop,iDOF)/(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst)
            CALL RANDOM_NUMBER(iRan)
            iQuant = INT(betaV+iRan)
            IF(iQuant.GT.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)) iQuant=PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)
            IF ((OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight).LT.0.0) THEN
              MaxColQua = OldEn/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*partWeight)
              IF (INT(MaxColQua).EQ.0) THEN
                iQuant = 0
              ELSE
                CALL RANDOM_NUMBER(iRan)
                iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                iQuaMax = MIN(INT(MaxColQua)+1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
                DO WHILE (iQuant.GE.iQuaMax)
                  CALL RANDOM_NUMBER(iRan)
                  iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                END DO
              END IF
            END IF
            PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1)  = PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
               + iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            VibQuantsPar(iPartIndx_NodeRelaxVib(iLoop))%Quants(iDOF) = iQuant
            OldEn = OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight
          END DO
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1)  = PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
               + SpecDSMC(1)%EZeroPoint
        END DO
      ELSE  ! Diatomic molecules
        DO iLoop = 1, nVibRelax
          partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
          betaV = alpha*PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1)/(SpecDSMC(1)%CharaTVib*BoltzmannConst)
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT(betaV+iRan)
          IF (iQuant.GT.SpecDSMC(1)%MaxVibQuant) iQuant = SpecDSMC(1)%MaxVibQuant
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1)  = (iQuant + DSMC%GammaQuant)*SpecDSMC(1)%CharaTVib*BoltzmannConst
          IF ((OldEn - (PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) - SpecDSMC(1)%EZeroPoint)*partWeight).LT.0.0) THEN
            MaxColQua = OldEn/(BoltzmannConst*SpecDSMC(1)%CharaTVib*partWeight)
            IF (INT(MaxColQua).EQ.0) THEN
              iQuant = 0
            ELSE
              CALL RANDOM_NUMBER(iRan)
              iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
              iQuaMax = MIN(INT(MaxColQua)+1, SpecDSMC(1)%MaxVibQuant)
              DO WHILE (iQuant.GE.iQuaMax)
                CALL RANDOM_NUMBER(iRan)
                iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
              END DO
            END IF
            PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1)  = (iQuant + DSMC%GammaQuant)*SpecDSMC(1)%CharaTVib*BoltzmannConst
          END IF
          OldEn = OldEn - (PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) - SpecDSMC(1)%EZeroPoint)*partWeight
        END DO
      END IF ! SpecDSMC(1)%PolyatomicMol
    ELSE ! Continuous treatment of vibrational energy
      alpha = OldEn/NewEnVib*(Xi_Vib*nVibRelax/(3.*(nPart-1.)+Xi_Vib*nVibRelax))
      DO iLoop = 1, nVibRelax
        partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
        PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = alpha*PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
          + SpecDSMC(1)%EZeroPoint
        OldEn = OldEn - (PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) - SpecDSMC(1)%EZeroPoint)*partWeight
      END DO
    END IF ! BGKUseQuantVibEn
  ELSE IF (nVibRelax.GT.0) THEN ! Relaxation towards the vibrational ground-state (new state is simply the zero-point energy)
    DO iLoop = 1, nVibRelax
      PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = SpecDSMC(1)%EZeroPoint
    END DO
  END IF ! (NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)
END IF ! BGKDoVibRelaxation

! 8.) Determine the new particle state (for molecules including rotational energy) and ensure energy conservation by scaling the new
!     velocities with the factor alpha. The actual update of particle velocity happens in the TimeDisc through the change in the
!     velocity (DSMC_RHS), to enable an easier coupling with existing routines and DSMC)
OldEn = OldEn + OldEnRot
IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  alpha = SQRT(OldEn/NewEn*(3.*(nRelax-1.))/(Xi_rot*nRotRelax+3.*(nRelax-1.)))
  DO iLoop = 1, nRelax
    DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkRelaxOld(1:3) &
                        + alpha*(DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)-vBulkRelax(1:3)) &
                        - PartState(iPartIndx_NodeRelax(iLoop),4:6)
  END DO
ELSE
  alpha = SQRT(OldEn/NewEn*(3.*(nPart-1.))/(Xi_rot*nRotRelax+3.*(nPart-1.)))
  DO iLoop = 1, nRelax
    DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkAll(1:3) + alpha*(DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)-vBulk(1:3)) &
                        - PartState(iPartIndx_NodeRelax(iLoop),4:6)
  END DO
  DO iLoop = 1, nPart-nRelax
    DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),1:3) = vBulkAll(1:3) &
                        + alpha*(PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6)-vBulk(1:3)) &
                        - PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6)
  END DO
END IF

! 9.) Rotation: Scale the new rotational state of the molecules to ensure energy conservation
IF ( (nRotRelax.GT.0)) alpha = OldEn/NewEnRot*(Xi_rot*nRotRelax/(Xi_rot*nRotRelax+3.*(nPart-1.)))
DO iLoop = 1, nRotRelax
  PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = alpha*PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2)
END DO

! CODE ANALYZE: Compare the old momentum and energy of the cell with the new, abort if relative difference is above the limits
#ifdef CODE_ANALYZE
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  Momentum_new(1:3) = Momentum_new(1:3) + (DSMC_RHS(iPartIndx_Node(iLoop),1:3) + PartState(iPartIndx_Node(iLoop),4:6))*partWeight
  Energy_new = Energy_new &
          + ((DSMC_RHS(iPartIndx_Node(iLoop),1) + PartState(iPartIndx_Node(iLoop),4))**2. &
          +  (DSMC_RHS(iPartIndx_Node(iLoop),2) + PartState(iPartIndx_Node(iLoop),5))**2. &
          +  (DSMC_RHS(iPartIndx_Node(iLoop),3) + PartState(iPartIndx_Node(iLoop),6))**2.)*0.5*Species(1)%MassIC*partWeight
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    Energy_new = Energy_new + (PartStateIntEn(iPartIndx_Node(iLoop),1) + PartStateIntEn(iPartIndx_Node(iLoop),2))*partWeight
  END IF
END DO
! Check for energy difference
IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,1.0e-12)) THEN
  WRITE(UNIT_StdOut,*) '\n'
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_old-Energy_new
  ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
    IF(energy.GT.0.0)THEN
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_old-Energy_new)/energy
    END IF
  END ASSOCIATE
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",1.0e-12
  IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha           : ", OldEn, alpha
  IPWRITE(UNIT_StdOut,*)                     " nPart, nRelax, nRotRelax, nVibRelax: ", nPart, nRelax, nRotRelax, nVibRelax
  CALL abort(&
      __STAMP__&
      ,'CODE_ANALYZE: BGK_CollisionOperator is not energy conserving!')
END IF
! Check for momentum difference
DO iMom=1,3
  IF (.NOT.ALMOSTEQUALRELATIVE(Momentum_old(iMom),Momentum_new(iMom),1.0e-9)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Direction (x,y,z)        : ",iMom
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_old             : ",Momentum_old(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_new             : ",Momentum_new(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Momentum difference : ",Momentum_old(iMom)-Momentum_new(iMom)
    ASSOCIATE( Momentum => MAX(ABS(Momentum_old(iMom)),ABS(Momentum_new(iMom))) )
      IF(Momentum.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Momentum difference : ",(Momentum_old(iMom)-Momentum_new(iMom))/Momentum
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance        : ",1.0e-9
    IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha             : ", OldEn, alpha
    IPWRITE(UNIT_StdOut,*)                     " nPart, nRelax, nRotRelax, nVibRelax: ", nPart, nRelax, nRotRelax, nVibRelax
    CALL abort(&
        __STAMP__&
        ,'CODE_ANALYZE: BGK_CollisionOperator is not momentum conserving!')
  END IF
END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE BGK_CollisionOperator


SUBROUTINE MetropolisES(nPart, iRanPart, A)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Basis ,ONLY: INV33
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: A(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: iRanPartTemp(3), V2, iRan, NewProb, OldProb, NormProb
INTEGER                        :: iLoop, iPart, iRun
LOGICAL                        :: Changed
REAL                           :: AC(3), AInvers(3,3), detA
!===================================================================================================================================
iRanPart(1,1) = rnor()
iRanPart(2,1) = rnor()
iRanPart(3,1) = rnor()
CALL INV33(A,AInvers, detA)
AC(1:3) = MATMUL(AInvers, iRanPart(1:3,1))
V2 = iRanPart(1,1)*AC(1) + iRanPart(2,1)*AC(2) + iRanPart(3,1)*AC(3)
OldProb = EXP(-0.5*V2)
!Burn in
DO iLoop = 1, 35 !50
  iRanPartTemp(1) = rnor()
  iRanPartTemp(2) = rnor()
  iRanPartTemp(3) = rnor()
  AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
  V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
  NewProb = EXP(-0.5*V2)
  NormProb = MIN(1.,NewProb/OldProb)
  CALL RANDOM_NUMBER(iRan)
  IF (NormProb.GT.iRan) THEN
    iRanPart(1:3,1) = iRanPartTemp(1:3)
    OldProb = NewProb
  END IF
END DO
! All the others
DO iPart = 2, nPart
  iRanPart(1,iPart) = iRanPart(1,iPart-1)
  iRanPart(2,iPart) = iRanPart(2,iPart-1)
  iRanPart(3,iPart) = iRanPart(3,iPart-1)
  iRun = 0
  Changed = .FALSE.
  DO WHILE ((iRun.LT.10).OR.(.NOT.Changed))
    iRun = iRun + 1
    iRanPartTemp(1) = rnor()
    iRanPartTemp(2) = rnor()
    iRanPartTemp(3) = rnor()
    AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
    V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
    NewProb = EXP(-0.5*V2)
    NormProb = MIN(1.,NewProb/OldProb)
    CALL RANDOM_NUMBER(iRan)
    IF (NormProb.GT.iRan) THEN
     Changed = .TRUE.
      iRanPart(1:3,iPart) = iRanPartTemp(1:3)
      OldProb = NewProb
    END IF
  END DO
END DO

END SUBROUTINE MetropolisES

SUBROUTINE ARShakhov(nPart, iRanPart, Vtherm, HeatVec, Prandtl)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Prandtl, Vtherm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vheat, V2, iRan, OldProb, Envelope
INTEGER                        :: iPart
!===================================================================================================================================
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope =  1.+4.*Envelope

DO iPart = 1, nPart
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARShakhov

SUBROUTINE MetropolisUnified(nPart, iRanPart,A, Vtherm, HeatVec, CShak)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Basis ,ONLY: INV33
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Vtherm
REAL, INTENT(IN)              :: A(3,3), CShak
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: iRanPartTemp(3), Vheat, V2, iRan, NewProb, OldProb, NormProb, V2ES
INTEGER                        :: iLoop, iPart, iRun
LOGICAL                        :: Changed
REAL                           :: AC(3), AInvers(3,3), detA
!===================================================================================================================================
iRanPart(1,1) = rnor()
iRanPart(2,1) = rnor()
iRanPart(3,1) = rnor()
CALL INV33(A,AInvers, detA)
AC(1:3) = MATMUL(AInvers, iRanPart(1:3,1))
V2ES = iRanPart(1,1)*AC(1) + iRanPart(2,1)*AC(2) + iRanPart(3,1)*AC(3)
V2 = iRanPart(1,1)*iRanPart(1,1) + iRanPart(2,1)*iRanPart(2,1) + iRanPart(3,1)*iRanPart(3,1)
Vheat = iRanPart(1,1)*HeatVec(1) + iRanPart(2,1)*HeatVec(2) + iRanPart(3,1)*HeatVec(3)
OldProb = EXP(-0.5*V2ES) +  EXP(-0.5*V2) * (1.-CShak)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.)
!Burn in
DO iLoop = 1, 35
  iRanPartTemp(1) = iRanPart(1,1)
  iRanPartTemp(2) = iRanPart(2,1)
  iRanPartTemp(3) = iRanPart(3,1)
  AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
  V2ES = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
  V2 = iRanPartTemp(1)*iRanPartTemp(1) + iRanPartTemp(2)*iRanPartTemp(2) + iRanPartTemp(3)*iRanPartTemp(3)
  Vheat = iRanPartTemp(1)*HeatVec(1) + iRanPartTemp(2)*HeatVec(2) + iRanPartTemp(3)*HeatVec(3)
  NewProb = EXP(-0.5*V2ES) +  EXP(-0.5*V2) * (1.-CShak)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.)
  NormProb = MIN(1.,NewProb/OldProb)
  CALL RANDOM_NUMBER(iRan)
  IF (NormProb.GT.iRan) THEN
    iRanPart(1:3,1) = iRanPartTemp(1:3)
    OldProb = NewProb
  END IF
END DO
! All the others
DO iPart = 2, nPart
  iRanPart(1,iPart) = iRanPart(1,iPart-1)
  iRanPart(2,iPart) = iRanPart(2,iPart-1)
  iRanPart(3,iPart) = iRanPart(3,iPart-1)
  iRun = 0
  Changed = .FALSE.
  DO WHILE ((iRun.LT.10).OR.(.NOT.Changed))
    iRun = iRun + 1
    iRanPartTemp(1) = iRanPart(1,iPart)
    iRanPartTemp(2) = iRanPart(2,iPart)
    iRanPartTemp(3) = iRanPart(3,iPart)
    AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
    V2ES = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
    V2 = iRanPartTemp(1)*iRanPartTemp(1) + iRanPartTemp(2)*iRanPartTemp(2) + iRanPartTemp(3)*iRanPartTemp(3)
    Vheat = iRanPartTemp(1)*HeatVec(1) + iRanPartTemp(2)*HeatVec(2) + iRanPartTemp(3)*HeatVec(3)
    NewProb = EXP(-0.5*V2ES) +  EXP(-0.5*V2) * (1.-CShak)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.)
    NormProb = MIN(1.,NewProb/OldProb)
    CALL RANDOM_NUMBER(iRan)
    IF (NormProb.GT.iRan) THEN
     Changed = .TRUE.
      iRanPart(1:3,iPart) = iRanPartTemp(1:3)
      OldProb = NewProb
    END IF
  END DO
END DO

END SUBROUTINE MetropolisUnified


SUBROUTINE BGK_BuildTransGaussNums(nPart, iRanPart)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iLoop
!===================================================================================================================================
DO iLoop = 1, nPart
  iRanPart(1,iLoop) = rnor()
  iRanPart(2,iLoop) = rnor()
  iRanPart(3,iLoop) = rnor()
END DO

END SUBROUTINE BGK_BuildTransGaussNums

SUBROUTINE CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp,  &
      TEqui, rotrelaxfreq, vibrelaxfreq, dtCell, DoVibRelaxIn)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_DSMC_Vars,              ONLY: SpecDSMC
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq, dtCell
INTEGER, INTENT(IN)             :: nPart
LOGICAL, OPTIONAL, INTENT(IN)   :: DoVibRelaxIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_vib, TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaR, betaV, RotFrac, VibFrac, TEqui_Old2
REAL                            :: eps_prec=1.0E-0
REAL                            :: correctFac, correctFacRot, maxexp    !, Xi_rel
LOGICAL                         :: DoVibRelax
!===================================================================================================================================
  IF (PRESENT(DoVibRelaxIn)) THEN
    DoVibRelax = DoVibRelaxIn
  ELSE
    DoVibRelax = BGKDoVibRelaxation
  END IF
  maxexp = LOG(HUGE(maxexp))
!  Xi_rel = 2.*(2. - CollInf%omegaLaux(1,1))
!  correctFac = 1. + (2.*SpecDSMC(1)%CharaTVib / (CellTemp*(EXP(SpecDSMC(1)%CharaTVib / CellTemp)-1.)))**2. &
!        * EXP(SpecDSMC(1)%CharaTVib /CellTemp) / (2.*Xi_rel)
!  correctFacRot = 1. + 2./Xi_rel

correctFac = 1.
correctFacRot = 1.
RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
RotFrac = nPart*(1.-RotExp)
IF(DoVibRelax) THEN
  VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
  VibFrac = nPart*(1.-VibExp)
ELSE
  VibExp = 0.0
  VibFrac = 0.0
  Xi_vib = 0.0
END IF
TEqui_Old = 0.0
TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib_old*VibFrac)
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  IF (ABS(TRot-TEqui).LT.1E-3) THEN
    RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
  ELSE
    betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dtCell/correctFacRot
    IF (-betaR.GT.0.0) THEN
      RotExp = 0.
    ELSE IF (betaR.GT.maxexp) THEN
      RotExp = 0.
    ELSE
      RotExp = exp(-betaR)
    END IF
  END IF
  RotFrac = nPart*(1.-RotExp)
  IF(DoVibRelax) THEN
    IF (ABS(TVib-TEqui).LT.1E-3) THEN
      VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
    ELSE
      betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dtCell/correctFac
      IF (-betaV.GT.0.0) THEN
        VibExp = 0.
      ELSE IF (betaV.GT.maxexp) THEN
        VibExp = 0.
      ELSE
        VibExp = exp(-betaV)
      END IF
    END IF
    IF ((SpecDSMC(1)%CharaTVib/TEqui).GT.maxexp) THEN
      Xi_Vib = 0.0
    ELSE
      Xi_vib = 2.*SpecDSMC(1)%CharaTVib/TEqui/(EXP(SpecDSMC(1)%CharaTVib/TEqui)-1.)
    END IF
    VibFrac = nPart*(1.-VibExp)
  END IF
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui
  TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib*VibFrac)
  IF(DoVibRelax) THEN
    DO WHILE( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      TEqui =(TEqui + TEqui_Old2)*0.5
      IF ((SpecDSMC(1)%CharaTVib/TEqui).GT.maxexp) THEN
        Xi_Vib = 0.0
      ELSE
        Xi_vib = 2.*SpecDSMC(1)%CharaTVib/TEqui/(EXP(SpecDSMC(1)%CharaTVib/TEqui)-1.)
      END IF
      TEqui_Old2 = TEqui
      TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib) / (3.*(nPart-1.)+2.*RotFrac+Xi_vib*VibFrac)
    END DO
  END IF
END DO

END SUBROUTINE CalcTEqui


SUBROUTINE CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_Vib_DOF, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq, &
      dtCell, DoVibRelaxIn)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_DSMC_Vars,              ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq
INTEGER, INTENT(IN)             :: nPart
REAL, INTENT(IN)                :: dtCell
LOGICAL, OPTIONAL, INTENT(IN)   :: DoVibRelaxIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_vib_DOF(:), TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaR, betaV, RotFrac, VibFrac, Xi_Rot, TEqui_Old2
REAL                            :: eps_prec=1.0
REAL                            :: correctFac, correctFacRot, maxexp
INTEGER                         :: iDOF, iPolyatMole
LOGICAL                         :: DoVibRelax
!===================================================================================================================================
IF (PRESENT(DoVibRelaxIn)) THEN
  DoVibRelax = DoVibRelaxIn
ELSE
  DoVibRelax = BGKDoVibRelaxation
END IF

maxexp = LOG(HUGE(maxexp))
Xi_Rot =   SpecDSMC(1)%Xi_Rot
iPolyatMole = SpecDSMC(1)%SpecToPolyArray
!  Xi_rel = 2.*(2. - CollInf%omegaLaux(1,1))
!  correctFac = 0.0
!  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!    correctFac = correctFac &
!        + (2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (CellTemp           &
!        *(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / CellTemp)-1.)))**2.  &
!        * EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / CellTemp) / 2.
!  END DO
!  correctFac = 1. + correctFac/Xi_rel
!  correctFacRot = 1. + Xi_Rot/Xi_rel

correctFac = 1.
correctFacRot = 1.
RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
RotFrac = nPart*(1.-RotExp)
IF(DoVibRelax) THEN
  VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
  VibFrac = nPart*(1.-VibExp)
ELSE
  VibExp = 0.0
  VibFrac = 0.0
  Xi_vib_DOF = 0.0
END IF
TEqui_Old = 0.0
TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib_old*VibFrac)
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  IF (ABS(TRot-TEqui).LT.1E-3) THEN
    RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
  ELSE
    betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dtCell/correctFacRot
    IF (-betaR.GT.0.0) THEN
      RotExp = 0.
    ELSE IF (betaR.GT.maxexp) THEN
      RotExp = 0.
    ELSE
      RotExp = exp(-betaR)
    END IF
  END IF
  RotFrac = nPart*(1.-RotExp)
  IF(DoVibRelax) THEN
    IF (ABS(TVib-TEqui).LT.1E-3) THEN
      VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
    ELSE
      betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dtCell/correctFac
      IF (-betaV.GT.0.0) THEN
        VibExp = 0.
      ELSE IF (betaV.GT.maxexp) THEN
        VibExp = 0.
      ELSE
        VibExp = exp(-betaV)
      END IF
    END IF
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      IF ((PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui).LT.maxexp) THEN
        Xi_vib_DOF(iDOF) = 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui &
                                    /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui)-1.)
      ELSE
        Xi_vib_DOF(iDOF) = 0.0
      END IF
    END DO
    VibFrac = nPart*(1.-VibExp)
  END IF
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui
  TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)  &
          / (3.*(nPart-1.)+2.*RotFrac+SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))*VibFrac)
  IF(DoVibRelax) THEN
    DO WHILE( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      TEqui =(TEqui + TEqui_Old2)*0.5
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF ((PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui).LT.maxexp) THEN
          Xi_vib_DOF(iDOF) = 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui &
                                      /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui)-1.)
        ELSE
          Xi_vib_DOF(iDOF) = 0.0
        END IF
      END DO
      TEqui_Old2 = TEqui
      TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)  &
          / (3.*(nPart-1.)+2.*RotFrac+SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))*VibFrac)
    END DO
  END IF
END DO

END SUBROUTINE CalcTEquiPoly

END MODULE MOD_BGK_CollOperator
