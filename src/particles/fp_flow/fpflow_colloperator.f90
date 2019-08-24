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

MODULE MOD_FP_CollOperator
!===================================================================================================================================
! Module containing the Fokker-Planck-based approximation of the collision operator
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE FP_CollisionOperator
  MODULE PROCEDURE FP_CollisionOperator
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: FP_CollisionOperator
!===================================================================================================================================

CONTAINS

SUBROUTINE FP_CollisionOperator(iPartIndx_Node, nPart, NodeVolume, vBulkAll)
!===================================================================================================================================
! Cell-local collision operator using different Fokker-Planck-based approximations
!===================================================================================================================================
! MODULES
#ifdef CODE_ANALYZE
USE MOD_Globals                 ,ONLY: abort,unit_stdout,myrank
#endif /* CODE_ANALYZE */
USE MOD_Globals_Vars            ,ONLY: Pi, BoltzmannConst
USE MOD_FPFlow_Vars             ,ONLY: FPCollModel, ESFPModel, SpecFP, FPUseQuantVibEn, FPDoVibRelaxation, FP_PrandtlNumber
USE MOD_FPFlow_Vars             ,ONLY: FP_MaxRelaxFactor, FP_MaxRotRelaxFactor, FP_MeanRelaxFactor, FP_MeanRelaxFactorCounter
USE MOD_Particle_Vars           ,ONLY: Species, PartState, VarTimeStep, usevMPF
USE MOD_TimeDisc_Vars           ,ONLY: dt
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, DSMC, PartStateIntEn, PolyatomMolDSMC, DSMC_RHS, VibQuantsPar, RadialWeighting
USE Ziggurat
USE MOD_FPFlow_Init             ,ONLY: FP_BuildTransGaussNums
USE MOD_DSMC_Analyze            ,ONLY: CalcTVibPoly
USE MOD_BGK_CollOperator        ,ONLY: CalcTEquiPoly, CalcTEqui
USE MOD_part_tools              ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(IN)                        :: vBulkAll(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: FPCoeffMatr(9,9), FPSolVec(9,1)
REAL                  ::  u0ij(6), u2ij(6), u0ijk(10), u2, u4i(3), OldEn, NewEn, u2i(3), u4, vBulk(3), u0ijMat(3,3),u0i(3)
REAL                  :: dens, vmag2, relaxtime, Lambda, alpha, CellTemp
REAL                  :: V_rel(3), FP_FakA, FP_FakB, FP_FakC, Prandtl, InnerDOF, SMat(3,3), tempVelo(3), KronDelta
INTEGER               :: iLoop, fillMa1, fillMa2, fillMa3, iLoop2, iPolyatMole, iDOF, IPIV(9)
REAL,ALLOCATABLE      :: Ni(:,:), iRanPart(:,:)
REAL                  :: dynamicvis, relaxfreq
REAL                  :: MaxColQua, ERot, EVib, Xi_vib, collisionfreq, vibrelaxfreq, rotrelaxfreq
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelaxRot(:),iPartIndx_NodeRelaxVib(:)
REAL, ALLOCATABLE     :: Xi_vib_DOF(:), VibEnergyDOF(:,:)
REAL                  :: betaV
REAL                  :: A(3,3), Work(1000), W(3), trace, nu, Theta, iRan, NewEnRot, NewEnVib, OldEnRot
REAL                  :: ProbAddPart,  RotExp, VibExp, TEqui, TRot, TVib, xi_rot, xi_vib_old, dtCell
INTEGER               :: INFO, iQuant, iQuaMax, nRotRelax, nVibRelax, info_dgesv
REAL                  :: partWeight, totalWeight
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

FPCoeffMatr= 0.0; dens = 0.0; u0ij=0.0; u2ij=0.0; u0ijk=0.0; u2=0.0; u4i=0.0; u2i=0.0; u4=0.0; u0ijMat=0.0; u0i=0.0
ALLOCATE(Ni(3,nPart))
Ni=0.0; ERot=0.0; EVib=0.0; OldEn = 0.0; NewEn = 0.0; OldEnRot = 0.0; NewEnRot = 0.0; NewEnVib=0.0; nRotRelax = 0; nVibRelax = 0
totalWeight = 0.0; dtCell = 0.0

DO iLoop2 = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop2))
  totalWeight = totalWeight + partWeight
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop2),4:6)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  u2= u2 + vmag2 * partWeight
  IF (FPCollModel.EQ.1) THEN
    u4i(1:3)=u4i(1:3) + vmag2 * vmag2 * V_rel(1:3) * partWeight
    iLoop = 1
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ij(iLoop)= u0ij(iLoop) + V_rel(fillMa1)*V_rel(fillMa2) * partWeight
        u2ij(iLoop)= u2ij(iLoop) + vmag2*V_rel(fillMa1)*V_rel(fillMa2) * partWeight
        iLoop = iLoop + 1
      END DO
    END DO
    iLoop = 1
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        DO fillMa3 = fillMa2 ,3
          u0ijk(iLoop)= u0ijk(iLoop) + V_rel(fillMa1)*V_rel(fillMa2)*V_rel(fillMa3) * partWeight
          iLoop = iLoop + 1
        END DO
      END DO
    END DO
  ELSE IF (FPCollModel.EQ.2) THEN
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ijMat(fillMa1, fillMa2)= u0ijMat(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2) * partWeight
      END DO
    END DO
    u0i(1:3) = u0i(1:3) + V_rel(1:3) * partWeight
  END IF
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    IF(FPDoVibRelaxation) Evib = Evib + (PartStateIntEn(iPartIndx_Node(iLoop2),1) - SpecDSMC(1)%EZeroPoint) * partWeight
    ERot = ERot + PartStateIntEn(iPartIndx_Node(iLoop2),2) * partWeight
  END IF
  OldEn = OldEn + 0.5*Species(1)%MassIC * vmag2 * partWeight
  IF(VarTimeStep%UseVariableTimeStep) THEN
    dtCell = dtCell + VarTimeStep%ParticleTimeStep(iPartIndx_Node(iLoop))
  END IF
END DO

IF(VarTimeStep%UseVariableTimeStep) THEN
  dtCell = dt * dtCell / nPart
ELSE
  dtCell = dt
END IF

u2 = u2*nPart/((nPart-1.)*totalWeight)
IF (FPCollModel.EQ.1) THEN
  u4i = u4i/totalWeight
  u0ij = u0ij/totalWeight
  u2ij = u2ij/totalWeight
  u0ijk = u0ijk/totalWeight
  u2i(1) = u0ijk(1)+u0ijk(4)+u0ijk(6)
  u2i(2) = u0ijk(2)+u0ijk(7)+u0ijk(9)
  u2i(3) = u0ijk(3)+u0ijk(8)+u0ijk(10)
  u4=u2ij(1)+u2ij(4)+u2ij(6)
ELSE IF (FPCollModel.EQ.2) THEN
  u0ijMat = u0ijMat / totalWeight
  A = u0ijMat
  CALL DSYEV('N','U',3,A,3,W,Work,1000,INFO)
  trace = u0ijMat(1,1)+u0ijMat(2,2)+u0ijMat(3,3)
  u0i = u0i / totalWeight
END IF

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  Xi_vib = 0.0
  IF(FPDoVibRelaxation) THEN
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
  Xi_rot = 0.0
  Xi_Vib = 0.0
  InnerDOF = 0.
END IF

Prandtl =2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
CellTemp = Species(1)%MassIC * u2/(3.0*BoltzmannConst)
TEqui = CellTemp
nu= 1.-3./(2.*Prandtl)
Theta = BoltzmannConst*CellTemp/Species(1)%MassIC
nu= MAX(nu,-Theta/(W(3)-Theta))

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  dens = totalWeight / NodeVolume
ELSE
  dens = totalWeight * Species(1)%MacroParticleFactor / NodeVolume
END IF

dynamicvis = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*SpecDSMC(1)%TrefVHS/Pi) &
        /(4.*(4.- 2.*SpecDSMC(1)%omegaVHS) * (6. - 2.*SpecDSMC(1)%omegaVHS)* SpecDSMC(1)%DrefVHS**2.)
relaxfreq = dens*BoltzmannConst*SpecDSMC(1)%TrefVHS**(SpecDSMC(1)%omegaVHS + 0.5) &
      /dynamicvis*CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
IF (FPCollModel.EQ.2) THEN
!  relaxtime = 2.0*(1.-nu)/relaxfreq
  relaxtime = 3.0/(Prandtl*relaxfreq)
ELSE
  relaxtime = 2.0/relaxfreq
END IF

IF(DSMC%CalcQualityFactors) THEN
  FP_MeanRelaxFactor         = FP_MeanRelaxFactor + dtCell / relaxtime
  FP_MeanRelaxFactorCounter  = FP_MeanRelaxFactorCounter + 1
  FP_MaxRelaxFactor          = MAX(FP_MaxRelaxFactor,dtCell / relaxtime)
  IF(FPCollModel.EQ.2) THEN
    FP_PrandtlNumber = FP_PrandtlNumber + (3./2.) / (1.-nu)
  END IF
END IF

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
! 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency using the collision frequency,
!     which is not the same as the relaxation frequency of distribution function, calculated above.
  collisionfreq = SpecFP(1)%CollFreqPreFactor(1) * dens *CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
  rotrelaxfreq = collisionfreq * DSMC%RotRelaxProb
  vibrelaxfreq = collisionfreq * DSMC%VibRelaxProb
  IF(SpecDSMC(1)%PolyatomicMol) THEN
    CALL CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_vib_DOF, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq, &
                        dtCell, DoVibRelaxIn=FPDoVibRelaxation)
      Xi_vib = SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ELSE
    CALL CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq, &
                    dtCell, DoVibRelaxIn=FPDoVibRelaxation)
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    FP_MaxRotRelaxFactor          = MAX(FP_MaxRotRelaxFactor,rotrelaxfreq*dtCell)
  END IF
! 4.) Determine the number of particles undergoing a relaxation (including rotational and vibrational relaxation for molecules)
  ALLOCATE(iPartIndx_NodeRelaxRot(nPart),iPartIndx_NodeRelaxVib(nPart))
  DO iLoop = 1, nPart
    !Rotation
    partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-RotExp
    IF (ProbAddPart.GT.iRan) THEN
      nRotRelax = nRotRelax + 1
      iPartIndx_NodeRelaxRot(nRotRelax) = iPartIndx_Node(iLoop)
      OldEnRot = OldEnRot + PartStateIntEn(iPartIndx_Node(iLoop),2) * partWeight
    END IF
    ! Vibration
    IF(FPDoVibRelaxation) THEN
      CALL RANDOM_NUMBER(iRan)
      ProbAddPart = 1.-VibExp
      IF (ProbAddPart.GT.iRan) THEN
        nVibRelax = nVibRelax + 1
        iPartIndx_NodeRelaxVib(nVibRelax) = iPartIndx_Node(iLoop)
        OldEn = OldEn + (PartStateIntEn(iPartIndx_NodeRelaxVib(nVibRelax),1) - SpecDSMC(1)%EZeroPoint) * partWeight
      END IF
    END IF
  END DO
  !!!! VIB RElaxation
  IF(FPDoVibRelaxation) THEN
    IF (nVibRelax.GT.0) THEN
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
  END IF
  !! ROT Relaxation
  DO iLoop = 1, nRotRelax
    partWeight = GetParticleWeight(iPartIndx_NodeRelaxRot(iLoop))
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = -Xi_Rot / 2. * BoltzmannConst*TEqui*LOG(iRan)
    NewEnRot = NewEnRot + PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) * partWeight
  END DO
END IF

ALLOCATE(iRanPart(3,nPart))
!IF ((FPCollModel.EQ.1).AND.(nPart.GE.5)) THEN
IF (FPCollModel.EQ.1) THEN
  Lambda = (u0ij(1)-1./3.*u2)**2.+2.0*u0ij(2)**2.+2.*u0ij(3)**2.+(u0ij(4)-1./3.*u2)**2.+2.*u0ij(5)**2.+(u0ij(6)-1./3.*u2)**2
  Lambda = -1.*Lambda/(relaxtime*u2**4.0)

  FPCoeffMatr(1,1) = 2.0 * u0ij(1)
  FPCoeffMatr(1,2) = 2.0 * u0ij(2)
  FPCoeffMatr(1,3) = 2.0 * u0ij(3)
  FPCoeffMatr(1,4) = 0.0
  FPCoeffMatr(1,5) = 0.0
  FPCoeffMatr(1,6) = 0.0
  FPCoeffMatr(1,7) = 2.0 * u2i(1)
  FPCoeffMatr(1,8) = 0.0
  FPCoeffMatr(1,9) = 0.0

  FPCoeffMatr(2,1) = u0ij(2)
  FPCoeffMatr(2,2) = u0ij(1) + u0ij(4)
  FPCoeffMatr(2,3) = u0ij(5)
  FPCoeffMatr(2,4) = u0ij(2)
  FPCoeffMatr(2,5) = u0ij(3)
  FPCoeffMatr(2,6) = 0.0
  FPCoeffMatr(2,7) = u2i(2)
  FPCoeffMatr(2,8) = u2i(1)
  FPCoeffMatr(2,9) = 0.0

  FPCoeffMatr(3,1) = u0ij(3)
  FPCoeffMatr(3,2) = u0ij(5)
  FPCoeffMatr(3,3) = u0ij(1) + u0ij(6)
  FPCoeffMatr(3,4) = 0.0
  FPCoeffMatr(3,5) = u0ij(2)
  FPCoeffMatr(3,6) = u0ij(3)
  FPCoeffMatr(3,7) = u2i(3)
  FPCoeffMatr(3,8) = 0.0
  FPCoeffMatr(3,9) = u2i(1)

  FPCoeffMatr(4,1) = 0.0
  FPCoeffMatr(4,2) = 2.0*u0ij(2)
  FPCoeffMatr(4,3) = 0.0
  FPCoeffMatr(4,4) = 2.0*u0ij(4)
  FPCoeffMatr(4,5) = 2.0*u0ij(5)
  FPCoeffMatr(4,6) = 0.0
  FPCoeffMatr(4,7) = 0.0
  FPCoeffMatr(4,8) = 2.0*u2i(2)
  FPCoeffMatr(4,9) = 0.0

  FPCoeffMatr(5,1) = 0.0
  FPCoeffMatr(5,2) = u0ij(3)
  FPCoeffMatr(5,3) = u0ij(2)
  FPCoeffMatr(5,4) = u0ij(5)
  FPCoeffMatr(5,5) = u0ij(4) + u0ij(6)
  FPCoeffMatr(5,6) = u0ij(5)
  FPCoeffMatr(5,7) = 0.0
  FPCoeffMatr(5,8) = u2i(3)
  FPCoeffMatr(5,9) = u2i(2)

  FPCoeffMatr(6,1) = 0.0
  FPCoeffMatr(6,2) = 0.0
  FPCoeffMatr(6,3) = 2.0 * u0ij(3)
  FPCoeffMatr(6,4) = 0.0
  FPCoeffMatr(6,5) = 2.0 * u0ij(5)
  FPCoeffMatr(6,6) = 2.0 * u0ij(6)
  FPCoeffMatr(6,7) = 0.0
  FPCoeffMatr(6,8) = 0.0
  FPCoeffMatr(6,9) = 2.0 *u2i(3)

  FPCoeffMatr(7,1) = u2i(1) + 2.0*u0ijk(1)
  FPCoeffMatr(7,2) = u2i(2) + 4.0*u0ijk(2)
  FPCoeffMatr(7,3) = u2i(3) + 4.0*u0ijk(3)
  FPCoeffMatr(7,4) = 2.0*u0ijk(4)
  FPCoeffMatr(7,5) = 4.0*u0ijk(5)
  FPCoeffMatr(7,6) = 2.0*u0ijk(6)
  FPCoeffMatr(7,7) = u4 - u2**2.0 + 2.0*u2ij(1) - 2.0*u0ij(1)*u2
  FPCoeffMatr(7,8) = 2.0*u2ij(2) - 2.0*u0ij(2)*u2
  FPCoeffMatr(7,9) = 2.0*u2ij(3) - 2.0*u0ij(3)*u2

  FPCoeffMatr(8,1) = 2.0*u0ijk(2)
  FPCoeffMatr(8,2) = u2i(1) + 4.0*u0ijk(4)
  FPCoeffMatr(8,3) = 4.0*u0ijk(5)
  FPCoeffMatr(8,4) = u2i(2) + 2.0*u0ijk(7)
  FPCoeffMatr(8,5) = u2i(3) + 4.0*u0ijk(8)
  FPCoeffMatr(8,6) = 2.0*u0ijk(9)
  FPCoeffMatr(8,7) = 2.0*u2ij(2) - 2.0*u0ij(2)*u2
  FPCoeffMatr(8,8) = u4 -u2**2.0 + 2.0*u2ij(4) - 2.0*u0ij(4)*u2
  FPCoeffMatr(8,9) = 2.0*u2ij(5) - 2.0*u0ij(5)*u2

  FPCoeffMatr(9,1) = 2.0*u0ijk(3)
  FPCoeffMatr(9,2) = 4.0*u0ijk(5)
  FPCoeffMatr(9,3) = u2i(1) + 4.0*u0ijk(6)
  FPCoeffMatr(9,4) = 2.0*u0ijk(8)
  FPCoeffMatr(9,5) = u2i(2) + 4.0*u0ijk(9)
  FPCoeffMatr(9,6) = u2i(3) + 2.0*u0ijk(10)
  FPCoeffMatr(9,7) = 2.0*u2ij(3) - 2.0*u0ij(3)*u2
  FPCoeffMatr(9,8) = 2.0*u2ij(5) - 2.0*u0ij(5)*u2
  FPCoeffMatr(9,9) = u4 -u2**2.0 + 2.0*u2ij(6) - 2.0*u0ij(6)*u2

  FPSolVec(1:6,1)=-2.0*Lambda*u2ij(1:6)
  FPSolVec(7,1)=5.0/(3.0*relaxtime)*u2i(1) + Lambda*(-3.0*u4i(1)+u2*u2i(1)+2.0*(u0ij(1)*u2i(1)+u0ij(2)*u2i(2)+u0ij(3)*u2i(3)))
  FPSolVec(8,1)=5.0/(3.0*relaxtime)*u2i(2) + Lambda*(-3.0*u4i(2)+u2*u2i(2)+2.0*(u0ij(2)*u2i(1)+u0ij(4)*u2i(2)+u0ij(5)*u2i(3)))
  FPSolVec(9,1)=5.0/(3.0*relaxtime)*u2i(3) + Lambda*(-3.0*u4i(3)+u2*u2i(3)+2.0*(u0ij(3)*u2i(1)+u0ij(5)*u2i(2)+u0ij(6)*u2i(3)))

  CALL DGESV(9, 1, FPCoeffMatr, 9, IPIV, FPSolVec, 9, info_dgesv)

  IF(info_dgesv.NE.0) THEN
    RETURN
  END IF
  DO iLoop = 1, nPart
    iRanPart(1,iLoop) = rnor()
    iRanPart(2,iLoop) = rnor()
    iRanPart(3,iLoop) = rnor()
  END DO
ELSE IF (FPCollModel.EQ.2) THEN
  IF (ESFPModel.EQ.1) THEN
    W = 0.
    A = 0.
    !! Exact Solution
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        IF (fillMa1.EQ.fillMa2) THEN
          KronDelta = 1.0
        ELSE
          KronDelta = 0.0
        END IF
        A(fillMa1, fillMa2) = (1.-nu)*Theta*KronDelta + nu*u0ijMat(fillMa1, fillMa2)
      END DO
    END DO
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
          SMat(fillMa1, fillMa2)= KronDelta - (3.-2.*Prandtl)/(4.*Prandtl)&
            *(Species(1)%MassIC/(BoltzmannConst*CellTemp)*nPart/(nPart-1.) &
            *(u0ijMat(fillMa1, fillMa2)-u0i(fillMa1)*u0i(fillMa2))-KronDelta)
        END DO
      END DO
      SMat(2,1)=SMat(1,2)
      SMat(3,1)=SMat(1,3)
      SMat(3,2)=SMat(2,3)
      SMat = SQRT(BoltzmannConst*CellTemp/Species(1)%MassIC)*SMat
    ELSE
      SMat(1,1) = SQRT(W(1))
      SMat(2,2) = SQRT(W(2))
      SMat(3,3) = SQRT(W(3))
      SMat = MATMUL(A, SMat)
      SMat = MATMUL(SMat, TRANSPOSE(A))
    END IF
  ELSE
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        IF (fillMa1.EQ.fillMa2) THEN
          KronDelta = 1.0
        ELSE
          KronDelta = 0.0
        END IF
        SMat(fillMa1, fillMa2)= KronDelta - (3.-2.*Prandtl)/(4.*Prandtl)&
          *(Species(1)%MassIC/(BoltzmannConst*CellTemp)*nPart/(nPart-1.) &
          *(u0ijMat(fillMa1, fillMa2)-u0i(fillMa1)*u0i(fillMa2))-KronDelta)
      END DO
    END DO
    SMat(2,1)=SMat(1,2)
    SMat(3,1)=SMat(1,3)
    SMat(3,2)=SMat(2,3)
  END IF

  DO iLoop = 1, nPart
    iRanPart(1,iLoop) = rnor()
    iRanPart(2,iLoop) = rnor()
    iRanPart(3,iLoop) = rnor()
  END DO
END IF

IF (FPCollModel.EQ.1) THEN
  FP_FakA = EXP(-dtCell/relaxtime)
  FP_FakB = (1.-EXP(-dtCell/relaxtime))
  FP_FakC = SQRT(BoltzmannConst*CellTemp/Species(1)%MassIC*(1.-EXP(-2.0*dtCell/relaxtime)))
ELSE
  FP_FakA = EXP(-dtCell/relaxtime)
  IF ((FPCollModel.EQ.2).AND.(ESFPModel.EQ.1)) THEN
    FP_FakC = SQRT(TEqui/CellTemp*(1.-EXP(-2.0*dtCell/relaxtime)))
  ELSE
    FP_FakC = SQRT(BoltzmannConst*TEqui/Species(1)%MassIC*(1.-EXP(-2.0*dtCell/relaxtime)))
  END IF
END IF

IF(FPDoVibRelaxation) THEN
  IF ((NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)) THEN
    IF (FPUseQuantVibEn) THEN
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
      ELSE
        DO iLoop = 1, nVibRelax
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
          OldEn = OldEn - PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) +  SpecDSMC(1)%EZeroPoint
        END DO
      END IF
    ELSE
      alpha = OldEn/NewEnVib*(Xi_Vib*nVibRelax/(3.*(nPart-1.)+Xi_Vib*nVibRelax))
      DO iLoop = 1, nVibRelax
        PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = alpha*PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
          + SpecDSMC(1)%EZeroPoint
        OldEn = OldEn - (PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) - SpecDSMC(1)%EZeroPoint)*partWeight
      END DO
    END IF
  ELSE IF (nVibRelax.GT.0) THEN
    DO iLoop = 1, nVibRelax
      PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = SpecDSMC(1)%EZeroPoint
    END DO
  END IF
END IF
OldEn = OldEn + OldEnRot

vBulk(1:3) = 0.0

DO iLoop = 1, nPart
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  DSMC_RHS(iPartIndx_Node(iLoop),1:3) = 0.0
!  IF ((FPCollModel.EQ.1).AND.(nPart.GE.5)) THEN
  IF (FPCollModel.EQ.1) THEN
    Ni(1,iLoop)  = FPSolVec(1,1)*V_rel(1)+FPSolVec(2,1)*V_rel(2) +FPSolVec(3,1)*V_rel(3) &
            + FPSolVec(7,1) * (vmag2 - u2) &
            + Lambda*(V_rel(1)*vmag2 - u2i(1))
    Ni(2,iLoop)  = FPSolVec(2,1)*V_rel(1)+FPSolVec(4,1)*V_rel(2) +FPSolVec(5,1)*V_rel(3) &
            + FPSolVec(8,1) * (vmag2 - u2) &
            + Lambda*(V_rel(2)*vmag2 - u2i(2))
    Ni(3,iLoop)  = FPSolVec(3,1)*V_rel(1)+FPSolVec(5,1)*V_rel(2) +FPSolVec(6,1)*V_rel(3) &
            + FPSolVec(9,1) * (vmag2 - u2) &
            + Lambda*(V_rel(3)*vmag2 -  u2i(3))
    DSMC_RHS(iPartIndx_Node(iLoop),1:3) = relaxtime*FP_FakB*Ni(1:3,iLoop)
  END IF

  IF (FPCollModel.EQ.2) THEN
    tempVelo(1:3) = FP_FakC*iRanPart(1:3,iLoop)
    DSMC_RHS(iPartIndx_Node(iLoop),1:3) = DSMC_RHS(iPartIndx_Node(iLoop),1:3) &
              + V_rel(1:3)*FP_FakA + MATMUL(SMat,tempVelo)
  ELSE
    DSMC_RHS(iPartIndx_Node(iLoop),1:3) = DSMC_RHS(iPartIndx_Node(iLoop),1:3) &
              + V_rel(1:3)*FP_FakA + FP_FakC*iRanPart(1:3,iLoop)
  END IF
  vBulk(1:3) = vBulk(1:3) + DSMC_RHS(iPartIndx_Node(iLoop),1:3)
END DO

vBulk(1:3) = vBulk(1:3) / totalWeight
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  V_rel(1:3) = DSMC_RHS(iPartIndx_Node(iLoop),1:3) - vBulk(1:3)
  NewEn = NewEn + (V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2 )*0.5*Species(1)%MassIC*partWeight
END DO

alpha = SQRT(OldEn/NewEn*(3.*(nPart-1.))/(Xi_rot*nRotRelax+3.*(nPart-1.)))
DO iLoop = 1, nPart
  DSMC_RHS(iPartIndx_Node(iLoop),1:3) = alpha*(DSMC_RHS(iPartIndx_Node(iLoop),1:3)-vBulk(1:3)) + vBulkAll(1:3) &
    - PartState(iPartIndx_Node(iLoop),4:6)
END DO
IF ( (nRotRelax.GT.0)) alpha = OldEn/NewEnRot*(Xi_rot*nRotRelax/(Xi_rot*nRotRelax+3.*(nPart-1.)))
DO iLoop = 1, nRotRelax
  PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = alpha*PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2)
END DO

DEALLOCATE(Ni)

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
  IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha, nPart    : ", OldEn, alpha, nPart
  CALL abort(&
      __STAMP__&
      ,'CODE_ANALYZE: FP_CollisionOperator is not energy conserving!')
END IF
! Check for momentum difference
DO iMom=1,3
  IF (.NOT.ALMOSTEQUALRELATIVE(Momentum_old(iMom),Momentum_new(iMom),1.0e-8)) THEN
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
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance        : ",1.0e-8
    IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha, nPart      : ", OldEn, alpha, nPart
    CALL abort(&
        __STAMP__&
        ,'CODE_ANALYZE: FP_CollisionOperator is not momentum conserving!')
  END IF
END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE FP_CollisionOperator

END MODULE MOD_FP_CollOperator
