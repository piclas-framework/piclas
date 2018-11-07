#include "piclas.h"

MODULE MOD_ESBGK_CollOperator
!===================================================================================================================================
! Module solving collision operator using BGK
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ESBGK_CollisionOperatorOctree
  MODULE PROCEDURE ESBGK_CollisionOperator
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ESBGK_CollisionOperatorOctree, ESBGK_Euler
!===================================================================================================================================

CONTAINS

SUBROUTINE ESBGK_CollisionOperator(iPartIndx_Node, nPart, iElem, NodeVolume, vBulkAll, AveragingPara, CorrectStep)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars ,ONLY: PartState, Species
USE MOD_DSMC_Vars     ,ONLY: DSMC_RHS, SpecDSMC, DSMC, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar
USE MOD_DSMC_Analyze  ,ONLY: CalcTVibPoly
USE MOD_TimeDisc_Vars ,ONLY: dt, TEnd, Time
USE MOD_Globals_Vars  ,ONLY: Pi, BoltzmannConst
USE MOD_ESBGK_Vars    ,ONLY: SpecESBGK, ESBGKTempCorrectFact, ESBGKModel, BGKCollModel, BGKUnifiedCes
USE MOD_ESBGK_Vars    ,ONLY: BGKDiffEn, BGKTest, BGKDiffEn2, BGKDiffEn3, BGKDiffEn4, BGKAveragingLength, BGKDoAveraging
USE MOD_ESBGK_Vars    ,ONLY: BGKDoAveragingCorrect, BGKUseQuantVibEn, BGKDoVibRelaxation, ElemSplitCells,SBGKEnergyConsMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(IN)                     :: iElem
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(IN)                        :: vBulkAll(3)
REAL, INTENT(INOUT), OPTIONAL           :: AveragingPara(5,BGKAveragingLength)
INTEGER, INTENT(INOUT), OPTIONAL        :: CorrectStep
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: KronDelta, tempVelo(3), vBulk(3), u0ij(3,3), SMat(3,3), u2, V_rel(3), vmag2, u0i(3), u2i(3)
REAL                  :: alpha, CellTemp, dens, InnerDOF, dynamicvis, iRan, NewEn, OldEn, Prandtl, relaxfreq
REAL                  :: testEnNew, testEnOld, testMomNew(3), testMomOld(3)
REAL                  :: rotrelaxfreq, vibrelaxfreq, collisionfreq, ProbAddPart, Evib, Tvib, Xi_vib, TEqui, Xi_Vib_old, Xi_rot
REAL                  :: MaxColQua, ERot, TEquiV, TEquiR
INTEGER               :: iLoop, nRelax, fillMa1, fillMa2, iQuant, iQuaMax, iDOF, iPolyatMole
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelax(:),iPartIndx_NodeRelaxTemp(:),iPartIndx_NodeRelaxRot(:),iPartIndx_NodeRelaxVib(:)
REAL, ALLOCATABLE     :: iRanPart(:,:), Xi_vib_DOF(:), VibEnergyDOF(:,:)
REAL                  :: A(3,3), Work(1000), W(3), ATemp(3,3), trace, CShak, RMax
INTEGER               :: INFO, nNotRelax, nRotRelax, nVibRelax, localBGKModel
REAL                  :: TRot, betaV, OldEnRot, RotExp, VibExp, NewEnRot, NewEnVib, OldEnVib, vBulkRelaxOld(3),vBulkRelax(3)
REAL                  :: CellTempRelax, vBulkAver(3), u2Aver, nPartAver, meanV, TimeStep, correctFac,correctFacRot
!===================================================================================================================================
!siehe paper An efficient particle Fokker–Planck algorithm for rarefied gas flows, Gorji, JCP 2014  
!testMomNew = 0.0; testMomOld = 0.0; testEnNew = 0.0; testEnOld = 0.0
NewEn = 0.; OldEn = 0.
OldEnRot = 0.
NewEnRot = 0.
NewEnVib = 0.
u2 = 0.0
u0ij = 0.0
u0i = 0.0
u2i = 0.0
Evib = 0.0
ERot = 0.0
u2Aver = 0.0
vBulkRelax = 0.0
vBulkRelaxOld = 0.0
DO iLoop = 1, nPart
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulkAll(1:3)
  IF (BGKDoAveraging) u2Aver = u2Aver + PartState(iPartIndx_Node(iLoop),4)**2. &
      +PartState(iPartIndx_Node(iLoop),5)**2. + PartState(iPartIndx_Node(iLoop),6)**2.
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  u2= u2 + vmag2
  IF ((BGKCollModel.EQ.1).OR.(BGKCollModel.EQ.4).OR.(BGKCollModel.EQ.5)) THEN
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ij(fillMa1, fillMa2)= u0ij(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2)
      END DO
    END DO
    u0i(1:3) = u0i(1:3) + V_rel(1:3)
  END IF  
  IF ((BGKCollModel.EQ.2).OR.(BGKCollModel.EQ.4).OR.(BGKCollModel.EQ.5)) u2i(1:3) = u2i(1:3) + V_rel(1:3)*vmag2
  IF (SBGKEnergyConsMethod.NE.2) OldEn = OldEn + 0.5*Species(1)%MassIC * vmag2
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN 
    IF(BGKDoVibRelaxation) Evib = Evib + PartStateIntEn(iPartIndx_Node(iLoop),1) - SpecDSMC(1)%EZeroPoint
    ERot = ERot + PartStateIntEn(iPartIndx_Node(iLoop),2)
  END IF
END DO

u2i = u2i*nPart/((nPart-1.)*(nPart-2.))
u2 = u2/nPart
u0ij = u0ij/nPart
u0i(1:3) = u0i(1:3)/nPart
trace = u0ij(1,1)+u0ij(2,2)+u0ij(3,3)

CellTemp = Species(1)%MassIC * u2/(3.0*BoltzmannConst) *nPart/(nPart-1.)

u2 = u2*nPart/(nPart-1.)
IF (BGKCollModel.EQ.5) THEN
!  IF ((MAX(ABS(u2i(1)),ABS(u2i(2)),ABS(u2i(3)))/(u2/3.)**(3./2.)).GT.4.) THEN
!  IF ((ABS(u2i(1))/(u2/3.)**(3./2.)).GT.4.) THEN
  IF (CellTemp.GT.3000.) THEN
    localBGKModel = 1
  ELSE
    localBGKModel = 2
  END IF
ELSE
  localBGKModel = BGKCollModel
END IF

dens = nPart * Species(1)%MacroParticleFactor / NodeVolume
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  Xi_vib = 0.0
  IF(BGKDoVibRelaxation) THEN
    IF(SpecDSMC(1)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(1)%SpecToPolyArray
      ALLOCATE(Xi_vib_DOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
      Xi_vib_DOF(:) = 0.
      TVib = CalcTVibPoly(Evib/nPart, 1)
      IF (TVib.GT.0.0) THEN
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          Xi_vib = Xi_vib + 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVib &
                              /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVib) - 1.)
        END DO
      END IF
    ELSE
      TVib=Evib / (nPart*BoltzmannConst*SpecDSMC(1)%CharaTVib)
      IF (TVib.GT.0.0) THEN
        Tvib= SpecDSMC(1)%CharaTVib/LOG(1. + 1./(TVib))
        Xi_vib = 2.* Evib / (nPart*BoltzmannConst*Tvib)
      END IF
    END IF
    Xi_Vib_old = Xi_Vib
  ELSE  
    Xi_Vib_old = 0.0
    TVib = 0.0
  END IF
  Xi_rot = SpecDSMC(1)%Xi_Rot
  TRot = 2.*ERot/(Xi_rot*nPart*BoltzmannConst)
!  TEquiV = CalcTEquilibVib(nPart, CellTemp,TVib, Xi_Vib)
!  TEquiR = (3.*CellTemp* (nPart-1.)/nPart+Xi_rot*TRot)/(3.* (nPart-1.)/nPart+Xi_rot)

  IF (BGKDoAveraging) THEN
    IF (BGKDoAveragingCorrect) THEN
      CorrectStep = CorrectStep + 1
      IF (CorrectStep.GT.BGKAveragingLength) CorrectStep = 1
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
  InnerDOF = Xi_rot + Xi_Vib
ELSE
  InnerDOF = 0.
END IF

dynamicvis = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*SpecDSMC(1)%TrefVHS/Pi) &
        /(4.*(4.- 2.*SpecDSMC(1)%omegaVHS) * (6. - 2.*SpecDSMC(1)%omegaVHS)* SpecDSMC(1)%DrefVHS**2.)
Prandtl =2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
CShak= Prandtl*(1.-BGKUnifiedCes)
IF (localBGKModel.EQ.1) THEN
  relaxfreq = Prandtl*dens*BoltzmannConst*SpecDSMC(1)%TrefVHS**(SpecDSMC(1)%omegaVHS + 0.5) &
      /dynamicvis*CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
ELSE IF (localBGKModel.EQ.4) THEN
  relaxfreq = dens*BoltzmannConst*SpecDSMC(1)%TrefVHS**(SpecDSMC(1)%omegaVHS + 0.5) &
      /(dynamicvis*(1.-BGKUnifiedCes))*CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
ELSE
  relaxfreq = dens*BoltzmannConst*SpecDSMC(1)%TrefVHS**(SpecDSMC(1)%omegaVHS + 0.5) &
      /dynamicvis*CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
END IF

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  collisionfreq = SpecESBGK(1)%CollFreqPreFactor(1) * Dens *CellTempRelax**(-SpecDSMC(1)%omegaVHS +0.5)
  rotrelaxfreq = collisionfreq * DSMC%RotRelaxProb
  vibrelaxfreq = collisionfreq * DSMC%VibRelaxProb
  IF(SpecDSMC(1)%PolyatomicMol) THEN
    CALL CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_vib_DOF, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq)
    Xi_vib = SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ELSE
    CALL CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp,  &
      TEqui, rotrelaxfreq, vibrelaxfreq)
  END IF
END IF

vBulk(1:3) = 0.0
nRelax = 0
nNotRelax = 0
nRotRelax = 0
nVibRelax = 0
ALLOCATE(iPartIndx_NodeRelax(nPart), iPartIndx_NodeRelaxTemp(nPart))
iPartIndx_NodeRelaxTemp = 0
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  ALLOCATE(iPartIndx_NodeRelaxRot(nPart),iPartIndx_NodeRelaxVib(nPart))
  DO iLoop = 1, nPart
!    testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6) 
!    testEnOld = testEnOld + (PartState(iPartIndx_Node(iLoop),4)**2. + PartState(iPartIndx_Node(iLoop),5)**2. &
!            + PartState(iPartIndx_Node(iLoop),6)**2.)*0.5*Species(1)%MassIC
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-exp(-relaxfreq*dt)
    IF (ProbAddPart.GT.iRan) THEN
      nRelax = nRelax + 1
      iPartIndx_NodeRelax(nRelax) = iPartIndx_Node(iLoop)
    ELSE
      nNotRelax = nNotRelax + 1
      iPartIndx_NodeRelaxTemp(nNotRelax) = iPartIndx_Node(iLoop)
      vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6)
    END IF
    !Rotation
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-RotExp
    IF (ProbAddPart.GT.iRan) THEN
      nRotRelax = nRotRelax + 1
      iPartIndx_NodeRelaxRot(nRotRelax) = iPartIndx_Node(iLoop)
      OldEnRot = OldEnRot + PartStateIntEn(iPartIndx_Node(iLoop),2)
!      testEnOld = testEnOld + PartStateIntEn(iPartIndx_Node(iLoop),2)
    END IF
    ! Vibration
    IF(BGKDoVibRelaxation) THEN
      CALL RANDOM_NUMBER(iRan)
      ProbAddPart = 1.-VibExp
      IF (ProbAddPart.GT.iRan) THEN
        nVibRelax = nVibRelax + 1
        iPartIndx_NodeRelaxVib(nVibRelax) = iPartIndx_Node(iLoop) 
        OldEn = OldEn + PartStateIntEn(iPartIndx_NodeRelaxVib(nVibRelax),1) - SpecDSMC(1)%EZeroPoint
  !      testEnOld = testEnOld + PartStateIntEn(iPartIndx_NodeRelaxVib(nVibRelax),1)
      END IF
    END IF
  END DO
  IF ((nRelax.EQ.0).AND.(nRotRelax.EQ.0).AND.(nVibRelax.EQ.0)) RETURN
  !!!! VIB RElaxation
!   DO iLoop = 1, nVibRelax
!      MaxColQua = OldEn/(BoltzmannConst*SpecDSMC(1)%CharaTVib)
!      IF (INT(MaxColQua).EQ.0) THEN
!        iQuant = 0
!      ELSE
!        CALL RANDOM_NUMBER(iRan)
!        iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
!        iQuaMax = MIN(INT(MaxColQua)+1, SpecDSMC(1)%MaxVibQuant)
!        DO WHILE (iQuant.GE.iQuaMax)
!          CALL RANDOM_NUMBER(iRan)
!          iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
!        END DO
!      END IF
!      PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
!          = (iQuant + DSMC%GammaQuant)*SpecDSMC(1)%CharaTVib*BoltzmannConst
!!      testEnNew = testEnNew + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1)
!      OldEn = OldEn - PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) +  SpecDSMC(1)%EZeroPoint
!    END DO
    IF(BGKDoVibRelaxation) THEN
      IF(SpecDSMC(1)%PolyatomicMol) THEN
        ALLOCATE(VibEnergyDOF(nVibRelax,PolyatomMolDSMC(iPolyatMole)%VibDOF))
        DO iLoop = 1, nVibRelax
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1) = 0.0
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            CALL RANDOM_NUMBER(iRan)
            VibEnergyDOF(iLoop,iDOF) = - LOG(iRan)*Xi_vib_DOF(iDOF)/2.*TEqui*BoltzmannConst
            PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
                                                                + VibEnergyDOF(iLoop,iDOF)
          END DO
          NewEnVib = NewEnVib + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1)
        END DO
      ELSE
        DO iLoop = 1, nVibRelax
          CALL RANDOM_NUMBER(iRan)
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = -LOG(iRan)*Xi_vib/2.*TEqui*BoltzmannConst
          NewEnVib = NewEnVib + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1) 
        END DO
      END IF
    END IF
    !! ROT Relaxation
    DO iLoop = 1, nRotRelax
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = -Xi_Rot / 2. * BoltzmannConst*TEqui*LOG(iRan)
      NewEnRot = NewEnRot + PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2)
    END DO
ELSE
  DO iLoop = 1, nPart
  !  testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6) 
  !  testEnOld = testEnOld + (PartState(iPartIndx_Node(iLoop),4)**2. + PartState(iPartIndx_Node(iLoop),5)**2. &
  !          + PartState(iPartIndx_Node(iLoop),6)**2.)*0.5*Species(1)%MassIC
    CALL RANDOM_NUMBER(iRan)
    ProbAddPart = 1.-exp(-relaxfreq*dt)
    IF (ProbAddPart.GT.iRan) THEN
      nRelax = nRelax + 1
      iPartIndx_NodeRelax(nRelax) = iPartIndx_Node(iLoop)
      IF (SBGKEnergyConsMethod.EQ.2) vBulkRelaxOld(1:3) = vBulkRelaxOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6)
    ELSE
      nNotRelax = nNotRelax + 1
      iPartIndx_NodeRelaxTemp(nNotRelax) = iPartIndx_Node(iLoop)
      vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6)
    END IF
  END DO
  IF (nRelax.EQ.0) RETURN

  IF (SBGKEnergyConsMethod.EQ.2) THEN
    vBulkRelaxOld = vBulkRelaxOld / nRelax
    IF (nRelax.GT.2) THEN
      DO iLoop = 1, nRelax
        OldEn = OldEn + 0.5*Species(1)%MassIC & 
              *((PartState(iPartIndx_NodeRelax(iLoop),4)-vBulkRelaxOld(1))**2.0 &
              + (PartState(iPartIndx_NodeRelax(iLoop),5)-vBulkRelaxOld(2))**2.0 &
              + (PartState(iPartIndx_NodeRelax(iLoop),6)-vBulkRelaxOld(3))**2.0)
  !    testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_NodeRelax(iLoop),4:6) 
  !    testEnOld = testEnOld + (PartState(iPartIndx_NodeRelax(iLoop),4)**2. + PartState(iPartIndx_NodeRelax(iLoop),5)**2. &
  !            + PartState(iPartIndx_NodeRelax(iLoop),6)**2.)*0.5*Species(1)%MassIC
      END DO
    ELSE
      DO iLoop = 1, nPart
        OldEn = OldEn + 0.5*Species(1)%MassIC & 
              *((PartState(iPartIndx_Node(iLoop),4)-vBulkAll(1))**2.0 &
              + (PartState(iPartIndx_Node(iLoop),5)-vBulkAll(2))**2.0 &
              + (PartState(iPartIndx_Node(iLoop),6)-vBulkAll(3))**2.0)
  !      testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6) 
  !      testEnOld = testEnOld + (PartState(iPartIndx_Node(iLoop),4)**2. + PartState(iPartIndx_Node(iLoop),5)**2. &
  !              + PartState(iPartIndx_Node(iLoop),6)**2.)*0.5*Species(1)%MassIC
      END DO
    END IF
  END IF
END IF

IF (nRelax.GT.0) THEN
  ALLOCATE(iRanPart(3,nRelax))

  SELECT CASE(localBGKModel)
  CASE (1)
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
      CALL ESBGK_BuildTransGaussNums(nRelax, iRanPart)
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
        SMat(1,1) = SQRT(W(1))
        SMat(2,2) = SQRT(W(2))
        SMat(3,3) = SQRT(W(3))
        SMat = MATMUL(A, SMat)
        SMat = MATMUL(SMat, TRANSPOSE(A))
        CALL ESBGK_BuildTransGaussNums(nRelax, iRanPart)
      ELSE IF (ESBGKModel.EQ.3) THEN
        A(2,1)=A(1,2)
        A(3,1)=A(1,3)
        A(3,2)=A(2,3)
        CALL MetropolisES(nRelax, iRanPart, A)
      END IF
    END IF
  CASE (2)
!    CALL MetropolisShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
    CALL ARShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
  CASE (3)
    CALL ESBGK_BuildTransGaussNums(nRelax, iRanPart)
  CASE (4)
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
    IF ((localBGKModel.EQ.1).AND.(ESBGKModel.NE.3)) THEN
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
    vBulkRelax(1:3) = vBulkRelax(1:3) + DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)
  END DO
END IF
IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  vBulkRelax = vBulkRelax / nRelax
ELSE
  vBulkRelax = (vBulkRelax + vBulk)/nPart
  vBulk = vBulkRelax
END IF

IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  DO iLoop = 1, nRelax 
    V_rel(1:3) = DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) - vBulkRelax(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC
  END DO
ELSE
  DO iLoop = 1, nRelax 
    V_rel(1:3) = DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) - vBulk(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC
  END DO
  DO iLoop = 1, nPart-nRelax 
    V_rel(1:3) = PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6) - vBulk(1:3)
    NewEn = NewEn + (V_rel(1)**2. + V_rel(2)**2. + V_rel(3)**2.)*0.5*Species(1)%MassIC
  END DO
END IF

IF(BGKDoVibRelaxation) THEN
  IF ((NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)) THEN
    IF (BGKUseQuantVibEn) THEN
      alpha = OldEn/NewEnVib*(Xi_Vib*nVibRelax/(3.*(nPart-1.)+Xi_Vib*nVibRelax))

      IF(SpecDSMC(1)%PolyatomicMol) THEN
        DO iLoop = 1, nVibRelax
          PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = 0.0
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            betaV = alpha*VibEnergyDOF(iLoop,iDOF)/(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst)
            CALL RANDOM_NUMBER(iRan)
            iQuant = INT(betaV+iRan)
            IF(iQuant.GT.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)) iQuant=PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)
            IF ((OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst).LT.0.0) THEN
              MaxColQua = OldEn/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
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
            OldEn = OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
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
          IF ((OldEn - PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) +  SpecDSMC(1)%EZeroPoint).LT.0.0) THEN
            MaxColQua = OldEn/(BoltzmannConst*SpecDSMC(1)%CharaTVib)
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
        !  testEnNew = testEnNew + PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop),1)
        END DO
      END IF
    ELSE
      alpha = OldEn/NewEnVib*(Xi_Vib*nVibRelax/(3.*(nPart-1.)+Xi_Vib*nVibRelax)) 
      DO iLoop = 1, nVibRelax
        PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = alpha*PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) &
          + SpecDSMC(1)%EZeroPoint
        OldEn = OldEn - PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) +  SpecDSMC(1)%EZeroPoint
      END DO
    END IF
  ELSE IF (nVibRelax.GT.0) THEN
    DO iLoop = 1, nVibRelax
      PartStateIntEn(iPartIndx_NodeRelaxVib(iLoop), 1) = SpecDSMC(1)%EZeroPoint
    END DO 
  END IF
END IF
OldEn = OldEn + OldEnRot
IF ((SBGKEnergyConsMethod.EQ.2).AND.(nRelax.GT.2)) THEN
  alpha = SQRT(OldEn/NewEn*(3.*(nRelax-1.))/(Xi_rot*nRotRelax+3.*(nRelax-1.)))
  DO iLoop = 1, nRelax
    DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkRelaxOld(1:3) &
                        + alpha*(DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)-vBulkRelax(1:3)) &
                        - PartState(iPartIndx_NodeRelax(iLoop),4:6)
!    testMomNew(1:3) = testMomNew(1:3) + DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) + PartState(iPartIndx_NodeRelax(iLoop),4:6)
!    testEnNew = testEnNew &
!            + ((DSMC_RHS(iPartIndx_NodeRelax(iLoop),1) + PartState(iPartIndx_NodeRelax(iLoop),4))**2. &
!            + (DSMC_RHS(iPartIndx_NodeRelax(iLoop),2) + PartState(iPartIndx_NodeRelax(iLoop),5))**2. &
!            + (DSMC_RHS(iPartIndx_NodeRelax(iLoop),3) + PartState(iPartIndx_NodeRelax(iLoop),6))**2.)*0.5*Species(1)%MassIC
  END DO
ELSE
  alpha = SQRT(OldEn/NewEn*(3.*(nPart-1.))/(Xi_rot*nRotRelax+3.*(nPart-1.))) 
  DO iLoop = 1, nRelax
    DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulkAll(1:3) + alpha*(DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)-vBulk(1:3)) &
                        - PartState(iPartIndx_NodeRelax(iLoop),4:6)
  !  testMomNew(1:3) = testMomNew(1:3) + DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) + PartState(iPartIndx_NodeRelax(iLoop),4:6)
  !  testEnNew = testEnNew &
  !          + ((DSMC_RHS(iPartIndx_NodeRelax(iLoop),1) + PartState(iPartIndx_NodeRelax(iLoop),4))**2. &
  !          + (DSMC_RHS(iPartIndx_NodeRelax(iLoop),2) + PartState(iPartIndx_NodeRelax(iLoop),5))**2. &
  !          + (DSMC_RHS(iPartIndx_NodeRelax(iLoop),3) + PartState(iPartIndx_NodeRelax(iLoop),6))**2.)*0.5*Species(1)%MassIC
  END DO
  DO iLoop = 1, nPart-nRelax
    DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),1:3) = vBulkAll(1:3) &
                        + alpha*(PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6)-vBulk(1:3)) &
                        - PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6)
  !  testMomNew(1:3) = testMomNew(1:3) + DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),1:3) + PartState(iPartIndx_NodeRelaxTemp(iLoop),4:6)
  !  testEnNew = testEnNew &
  !          + ((DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),1) + PartState(iPartIndx_NodeRelaxTemp(iLoop),4))**2. &
  !          + (DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),2) + PartState(iPartIndx_NodeRelaxTemp(iLoop),5))**2. &
  !          + (DSMC_RHS(iPartIndx_NodeRelaxTemp(iLoop),3) + PartState(iPartIndx_NodeRelaxTemp(iLoop),6))**2.)*0.5*Species(1)%MassIC
  END DO
END IF
IF ( (nRotRelax.GT.0)) alpha = OldEn/NewEnRot*(Xi_rot*nRotRelax/(Xi_rot*nRotRelax+3.*(nPart-1.)))
DO iLoop = 1, nRotRelax
  PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2) = alpha*PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop), 2)
!   testEnNew = testEnNew + PartStateIntEn(iPartIndx_NodeRelaxRot(iLoop),2)
END DO

IF(DSMC%CalcQualityFactors) THEN
  IF(Time.GE.(1-DSMC%TimeFracSamp)*TEnd) THEN
    DSMC%CollSepDist = DSMC%CollSepDist + 1.
    ! Calculation of the lowest required timestep in the cell
    meanV = SQRT(vBulkAll(1)*vBulkAll(1)+vBulkAll(2)*vBulkAll(2)+vBulkAll(3)*vBulkAll(3)) &
              + SQRT(8.*BoltzmannConst*CellTemp/(Pi*Species(1)%MassIC))
    IF(meanV.GT.0.0) THEN
      TimeStep = NodeVolume**(1./3.) / meanV
      DSMC%CollProbMax = MIN(DSMC%CollProbMax,TimeStep)
    END IF
  END IF
END IF

!print*, OldEn, alpha
!print*, nPart, nRelax, nRotRelax, nVibRelax
!print*, testEnOld, testEnNew, testEnOld-testEnNew
!print*, testMomOld
!print*, testMomNew
!print*, testMomOld-testMomNew
!read*

END SUBROUTINE ESBGK_CollisionOperator



SUBROUTINE ESBGK_Euler(iPartIndx_Node, nPart, iElem, NodeVolume, vBulkAll)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars ,ONLY: PartState, Species
USE MOD_DSMC_Vars     ,ONLY: DSMC_RHS, SpecDSMC, PartStateIntEn
USE MOD_TimeDisc_Vars ,ONLY: dt
USE MOD_Globals_Vars  ,ONLY: Pi
USE MOD_ESBGK_Vars    ,ONLY: SpecESBGK, ESBGKTempCorrectFact, ESBGKModel, BGKCollModel, BGKUnifiedCes
USE MOD_ESBGK_Init    ,ONLY: ESBGK_BuildTransGaussNumsEnCon
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(IN)                     :: iElem
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(IN)                        :: vBulkAll(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: RelaxationFreq, KronDelta, tempVelo(3), vBulk(3), u0ij(3,3), SMat(3,3), u2, V_rel(3), vmag2, u0i(3), u2i(3)
REAL                  :: alpha, CellTemp, cp, dens, InnerDOF, dynamicvis, iRan, NewEn, OldEn, Prandtl, pressure, relaxfreq
REAL                  :: thermalconduct, relaxfreq2, dynamicvis2, Prandtl2, testEnNew, testEnOld, testMomNew(3), testMomOld(3)
REAL                  :: rotrelaxfreq, vibrelaxfreq, collisionfreq, ProbAddPart, Evib, Tvib, Xi_vib, TEqui
REAL                  :: MaxColQua, FakXi, OldEn2, iRan2, MomCor(3), PartTemp(nPart), NewEnRot, OldEnRot, Xi_Rot
INTEGER               :: iLoop
REAL, ALLOCATABLE     :: iRanPart(:,:)
!===================================================================================================================================
!siehe paper An efficient particle Fokker–Planck algorithm for rarefied gas flows, Gorji, JCP 2014  
!testMomNew = 0.0; testMomOld = 0.0; testEnNew = 0.0; testEnOld = 0.0
NewEn = 0.; OldEn = 0.
OldEnRot = 0.0; NewEnRot = 0.0
DO iLoop = 1, nPart
!  testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_Node(iLoop),4:6) 
!  testEnOld = testEnOld + (PartState(iPartIndx_Node(iLoop),4)**2. + PartState(iPartIndx_Node(iLoop),5)**2. &
!          + PartState(iPartIndx_Node(iLoop),6)**2.)*0.5*Species(1)%MassIC
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  OldEn = OldEn +  vmag2
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) OldEnRot = OldEnRot + PartStateIntEn(iPartIndx_Node(iLoop),2)
END DO
OldEn = 0.5*Species(1)%MassIC * OldEn
!testEnOld = testEnOld + OldEnRot
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  Xi_rot = SpecDSMC(1)%Xi_Rot
!  DO iLoop = 1, nPart
!    CALL RANDOM_NUMBER(iRan)
!    PartStateIntEn(iPartIndx_Node(iLoop), 2) = -LOG(iRan)
!    NewEnRot = NewEnRot + PartStateIntEn(iPartIndx_Node(iLoop), 2)
!  END DO
END IF

ALLOCATE(iRanPart(3,nPart))
CALL ESBGK_BuildTransGaussNumsEnCon(nPart, iRanPart)
DO iLoop = 1, nPart
  DSMC_RHS(iPartIndx_Node(iLoop),1:3) = iRanPart(1:3,iLoop)
  vmag2 = DSMC_RHS(iPartIndx_Node(iLoop),1)**2. + DSMC_RHS(iPartIndx_Node(iLoop),2)**2. &
        + DSMC_RHS(iPartIndx_Node(iLoop),3)**2.
  NewEn =  NewEn + vmag2
END DO
NewEn = NewEn*0.5*Species(1)%MassIC

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  OldEn = OldEn + OldEnRot
  alpha = SQRT(OldEn/NewEn*(3.*(nPart-1.))/(Xi_rot*nPart+3.*(nPart-1.)))
ELSE
  alpha = SQRT(OldEn/NewEn)                       ! create particle index list for pairin
END IF

DO iLoop = 1, nPart
  DSMC_RHS(iPartIndx_Node(iLoop),1:3) = vBulkAll(1:3) + alpha*DSMC_RHS(iPartIndx_Node(iLoop),1:3) &
                      - PartState(iPartIndx_Node(iLoop),4:6)
!  testMomNew(1:3) = testMomNew(1:3) + DSMC_RHS(iPartIndx_Node(iLoop),1:3) + PartState(iPartIndx_Node(iLoop),4:6)
!  testEnNew = testEnNew &
!          + ((DSMC_RHS(iPartIndx_Node(iLoop),1) + PartState(iPartIndx_Node(iLoop),4))**2. &
!          + (DSMC_RHS(iPartIndx_Node(iLoop),2) + PartState(iPartIndx_Node(iLoop),5))**2. &
!          + (DSMC_RHS(iPartIndx_Node(iLoop),3) + PartState(iPartIndx_Node(iLoop),6))**2.)*0.5*Species(1)%MassIC
END DO

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  alpha = OldEn/OldEnRot*(Xi_rot*nPart/(Xi_rot*nPart+3.*(nPart-1.)))
  DO iLoop = 1, nPart
    PartStateIntEn(iPartIndx_Node(iLoop), 2) = alpha*PartStateIntEn(iPartIndx_Node(iLoop), 2)
!    testEnNew = testEnNew + PartStateIntEn(iPartIndx_Node(iLoop),2)
  END DO
END IF

!print*, testEnOld, testEnNew, testEnOld-testEnNew
!print*, testMomOld
!print*, testMomNew
!print*, testMomOld-testMomNew
!read*

END SUBROUTINE ESBGK_Euler


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
REAL                           :: iRanPartTemp(3), Vheat, V2, Vtherm, iRan, NewProb, OldProb, NormProb
INTEGER                        :: iLoop, iPart, iRun
LOGICAL                        :: Changed
REAL                           :: AC(3), CAC, AInvers(3,3), detA
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

SUBROUTINE MetropolisShakhov(nPart, iRanPart, Vtherm, HeatVec, Prandtl)
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
REAL, INTENT(IN)              :: HeatVec(3), Prandtl, Vtherm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: iRanPartTemp(3), Vheat, V2, MeanVelo(3), iRan, NewProb, OldProb, NormProb
INTEGER                        :: iLoop, iPart, iRun, iCount
LOGICAL                        :: Changed
!===================================================================================================================================
iRanPart(1,1) = rnor()
iRanPart(2,1) = rnor()
iRanPart(3,1) = rnor()
V2 = iRanPart(1,1)*iRanPart(1,1) + iRanPart(2,1)*iRanPart(2,1) + iRanPart(3,1)*iRanPart(3,1)
Vheat = iRanPart(1,1)*HeatVec(1) + iRanPart(2,1)*HeatVec(2) + iRanPart(3,1)*HeatVec(3)
OldProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
IF (OldProb.LT.0.0) THEN
  iCount = 0
  DO WHILE (OldProb.LT.0.0)
    iCount = iCount + 1
    IF (iCount.EQ.100) EXIT
    iRanPart(1,1) = rnor()
    iRanPart(2,1) = rnor()
    iRanPart(3,1) = rnor()
    V2 = iRanPart(1,1)*iRanPart(1,1) + iRanPart(2,1)*iRanPart(2,1) + iRanPart(3,1)*iRanPart(3,1)
    Vheat = iRanPart(1,1)*HeatVec(1) + iRanPart(2,1)*HeatVec(2) + iRanPart(3,1)*HeatVec(3)
    OldProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
  END DO
END IF
!OldProb = EXP(-0.5*V2) *(1. + 4./5.*(1.-Prandtl)*VHeat/(Vtherm**(3./2.))*(V2-5./2.))
!Burn in
DO iLoop = 1, 35
  iRanPartTemp(1) = rnor()
  iRanPartTemp(2) = rnor()
  iRanPartTemp(3) = rnor()
  V2 = iRanPartTemp(1)*iRanPartTemp(1) + iRanPartTemp(2)*iRanPartTemp(2) + iRanPartTemp(3)*iRanPartTemp(3)
  Vheat = iRanPartTemp(1)*HeatVec(1) + iRanPartTemp(2)*HeatVec(2) + iRanPartTemp(3)*HeatVec(3)
  NewProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
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
    V2 = iRanPartTemp(1)*iRanPartTemp(1) + iRanPartTemp(2)*iRanPartTemp(2) + iRanPartTemp(3)*iRanPartTemp(3)
    Vheat = iRanPartTemp(1)*HeatVec(1) + iRanPartTemp(2)*HeatVec(2) + iRanPartTemp(3)*HeatVec(3)
    NewProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
    NormProb = MIN(1.,NewProb/OldProb)
    CALL RANDOM_NUMBER(iRan)
    IF (NormProb.GT.iRan) THEN
     Changed = .TRUE.
      iRanPart(1:3,iPart) = iRanPartTemp(1:3)
      OldProb = NewProb
    END IF
  END DO
END DO

END SUBROUTINE MetropolisShakhov

SUBROUTINE ARShakhov(nPart, iRanPart, Vtherm, HeatVec, Prandtl)
!===================================================================================================================================
! Performs FP Momentum Evaluation
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species
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
REAL, PARAMETER                  :: PI=3.14159265358979323846  
REAL                           :: iRanPartTemp(3), Vheat, V2, MeanVelo(3), iRan, NewProb, OldProb, NormProb, Envelope, Phi
REAL                           :: VeloCrad
INTEGER                        :: iLoop, iPart, iRun, iCount
LOGICAL                        :: Changed
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
REAL                           :: iRanPartTemp(3), Vheat, V2, MeanVelo(3), iRan, NewProb, OldProb, NormProb, V2ES
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


SUBROUTINE ESBGK_BuildTransGaussNums(nPart, iRanPart)
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

END SUBROUTINE ESBGK_BuildTransGaussNums

SUBROUTINE CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp,  &
      TEqui, rotrelaxfreq, vibrelaxfreq)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC
  USE MOD_ESBGK_Vars,             ONLY : BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq
  INTEGER, INTENT(IN)             :: nPart
  REAL, INTENT(OUT)               :: Xi_vib, TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL                    :: LowerTemp, UpperTemp, TEqui_Old, betaR, betaV, RotFrac, VibFrac, TEqui_Old2
  REAL                    :: eps_prec=1.0E-0
  REAL                    :: correctFac, correctFacRot, maxexp, Xi_rel
!===================================================================================================================================
  maxexp = LOG(HUGE(maxexp))
!  Xi_rel = 2.*(2. - SpecDSMC(1)%omegaVHS)
!  correctFac = 1. + (2.*SpecDSMC(1)%CharaTVib / (CellTemp*(EXP(SpecDSMC(1)%CharaTVib / CellTemp)-1.)))**2. & 
!        * EXP(SpecDSMC(1)%CharaTVib /CellTemp) / (2.*Xi_rel)
!  correctFacRot = 1. + 2./Xi_rel
  correctFac = 1.
  correctFacRot = 1.
  RotExp = exp(-rotrelaxfreq*dt/correctFacRot) 
  RotFrac = nPart*(1.-RotExp)
  IF(BGKDoVibRelaxation) THEN
    VibExp = exp(-vibrelaxfreq*dt/correctFac) 
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
      RotExp = exp(-rotrelaxfreq*dt/correctFacRot) 
    ELSE
      betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dt/correctFacRot
      IF (-betaR.GT.0.0) THEN
        RotExp = 0.
      ELSE IF (betaR.GT.maxexp) THEN
        RotExp = 0.
      ELSE
        RotExp = exp(-betaR) 
      END IF
    END IF
    RotFrac = nPart*(1.-RotExp)
    IF(BGKDoVibRelaxation) THEN
      IF (ABS(TVib-TEqui).LT.1E-3) THEN
        VibExp = exp(-vibrelaxfreq*dt/correctFac) 
      ELSE
        betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dt/correctFac
        IF (-betaV.GT.0.0) THEN
          VibExp = 0.
    !      VibExp = exp(-betaV) 
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
    IF(BGKDoVibRelaxation) THEN
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
! print*, betaR/(rotrelaxfreq*dt*correctFacRot), betaV/(vibrelaxfreq*dt*correctFac)
END SUBROUTINE CalcTEqui

SUBROUTINE CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_Vib_DOF, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PolyatomMolDSMC
  USE MOD_ESBGK_Vars,             ONLY : BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq
  INTEGER, INTENT(IN)             :: nPart
  REAL, INTENT(OUT)               :: Xi_vib_DOF(:), TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL                    :: LowerTemp, UpperTemp, TEqui_Old, betaR, betaV, RotFrac, VibFrac, Xi_Rot
  REAL                    :: eps_prec=1.0
  REAL                    :: correctFac, correctFacRot, maxexp, TEqui_Old2
  INTEGER                 :: iDOF, iPolyatMole
!===================================================================================================================================

  maxexp = LOG(HUGE(maxexp))
  Xi_Rot =   SpecDSMC(1)%Xi_Rot
  iPolyatMole = SpecDSMC(1)%SpecToPolyArray
!  Xi_rel = 2.*(2. - SpecDSMC(1)%omegaVHS)
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
  RotExp = exp(-rotrelaxfreq*dt/correctFacRot) 
  RotFrac = nPart*(1.-RotExp)
  IF(BGKDoVibRelaxation) THEN
    VibExp = exp(-vibrelaxfreq*dt/correctFac) 
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
      RotExp = exp(-rotrelaxfreq*dt/correctFacRot) 
    ELSE
      betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dt/correctFacRot
      IF (-betaR.GT.0.0) THEN
        RotExp = 0.
      ELSE IF (betaR.GT.maxexp) THEN
        RotExp = 0.
      ELSE
        RotExp = exp(-betaR) 
      END IF
    END IF
    RotFrac = nPart*(1.-RotExp)
    IF(BGKDoVibRelaxation) THEN
      IF (ABS(TVib-TEqui).LT.1E-3) THEN
        VibExp = exp(-vibrelaxfreq*dt/correctFac) 
      ELSE
        betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dt/correctFac
        IF (-betaV.GT.0.0) THEN
          VibExp = 0.
    !      VibExp = exp(-betaV) 
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
    IF(BGKDoVibRelaxation) THEN
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


!REAL FUNCTION CalcTEquilibVib(nPart, CellTemp,TVib, Xi_Vib)
!!===================================================================================================================================
!! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!!===================================================================================================================================
!! MODULES
!  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
!  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PolyatomMolDSMC
!! IMPLICIT VARIABLE HANDLING
!  IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!  REAL, INTENT(IN)                :: CellTemp, TVib  ! Charak TVib, mean vibrational Energy of all molecules
!  INTEGER, INTENT(IN)             :: nPart
!  REAL, INTENT(IN)             :: Xi_vib  ! Charak TVib, mean vibrational Energy of all molecules
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!  REAL                            :: LowerTemp, UpperTemp, TEqui_Old
!  REAL                            :: eps_prec=1.0E-5
!  REAL                            :: Xi_Vib_new, maxexp
!  INTEGER                         :: iDOF, iPolyatMole
!!===================================================================================================================================
!  maxexp = LOG(HUGE(maxexp))
!  LowerTemp = 0.9*MIN(CellTemp, TVib)
!  UpperTemp = 1.1*MAX(CellTemp, TVib)
!  IF(SpecDSMC(1)%PolyatomicMol) THEN
!    iPolyatMole = SpecDSMC(1)%SpecToPolyArray
!    DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
!      CalcTEquilibVib = 0.5*( LowerTemp + UpperTemp)
!      Xi_vib_new  = 0.0
!      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!        IF ((PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/CalcTEquilibVib).LT.maxexp) THEN
!          Xi_vib_new = Xi_vib_new + 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/CalcTEquilibVib &
!                                      /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/CalcTEquilibVib)-1.)
!        END IF
!      END DO
!      TEqui_Old = 3.*(nPart-1.)/nPart*CellTemp+Xi_Vib*TVib
!      IF (CalcTEquilibVib*(3.*(nPart-1.)/nPart+Xi_Vib_new).GT.TEqui_Old) THEN
!        UpperTemp = CalcTEquilibVib
!      ELSE
!        LowerTemp = CalcTEquilibVib
!      END IF
!    END DO
!  ELSE
!    DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
!      CalcTEquilibVib = 0.5*( LowerTemp + UpperTemp)
!      IF ((SpecDSMC(1)%CharaTVib/CalcTEquilibVib).GT.maxexp) THEN
!        Xi_vib_new  = 0.0
!      ELSE
!        Xi_vib_new = 2.*SpecDSMC(1)%CharaTVib/CalcTEquilibVib/(EXP(SpecDSMC(1)%CharaTVib/CalcTEquilibVib)-1.)
!      END IF
!      TEqui_Old = 3.*(nPart-1.)/nPart*CellTemp+Xi_Vib*TVib
!      IF (CalcTEquilibVib*(3.*(nPart-1.)/nPart+Xi_Vib_new).GT.TEqui_Old) THEN
!        UpperTemp = CalcTEquilibVib
!      ELSE
!        LowerTemp = CalcTEquilibVib
!      END IF
!    END DO
!  END IF
!  RETURN
!END FUNCTION CalcTEquilibVib

END MODULE MOD_ESBGK_CollOperator
