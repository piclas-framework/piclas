#include "boltzplatz.h"

MODULE MOD_ESBGK_MOLEC
!===================================================================================================================================
! Module for ESBGK Flow
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ESBGK_main
  MODULE PROCEDURE ESBGK_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ESBGK_main
!===================================================================================================================================

CONTAINS


SUBROUTINE ESBGK_main()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars ,ONLY: nElems
USE MOD_DSMC_Vars ,ONLY: DSMC_RHS
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem
!===================================================================================================================================
DSMC_RHS = 0.0
DO iElem = 1, nElems
  CALL ESBGK_CollisionOperator(iElem)
END DO

END SUBROUTINE ESBGK_main

SUBROUTINE ESBGK_CollisionOperator(iElem)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_Particle_Vars      ,ONLY: PEM, PartState, Species, BoltzmannConst
USE MOD_DSMC_Vars          ,ONLY: DSMC_RHS, SpecDSMC, DSMC, PartStateIntEn
USE MOD_TimeDisc_Vars      ,ONLY: dt
USE MOD_Globals_Vars       ,ONLY: Pi
USE MOD_ESBGK_Init         ,ONLY: ESBGK_BuildTransGaussNums
USE MOD_ESBGK_Vars         ,ONLY: SpecESBGK
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: RelaxationFreq, KronDelta, tempVelo(3), vBulk(3), u0ij(3,3), SMat(3,3), u2, V_rel(3), vmag2
REAL                  :: alpha, CellTemp, cp, dens, InnerDOF, dynamicvis, iRan, NewEn, OldEn, Prandl, pressure, relaxfreq
REAL                  :: thermalconduct, relaxfreq2, dynamicvis2, Prandl2, testEnNew, testEnOld, testMomNew(3), testMomOld(3)
REAL                  :: rotrelaxfreq, vibrelaxfreq, collisionfreq, ProbAddPart, Evib, Tvib, Xi_vib, TEqui
INTEGER               :: nPart, iPart, iLoop, nRelax, fillMa1, fillMa2, nRelaxRot, nRelaxVib, iQuant, newnRelaxVib, newnRelaxRot
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:), iPartIndx_NodeRelax(:)
REAL, ALLOCATABLE     :: iRanPart(:,:)
LOGICAL, ALLOCATABLE  :: DoVibRelax(:), DoRotRelax(:)
!===================================================================================================================================

nPart = PEM%pNumber(iElem)
IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

vBulk = 0.0; Evib = 0.0
ALLOCATE(iPartIndx_Node(nPart))
iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
DO iLoop = 1, nPart
  iPartIndx_Node(iLoop) = iPart
  vBulk(1:3)  =  vBulk(1:3) + PartState(iPart,4:6)
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) Evib = Evib + PartStateIntEn(iPart,1) - SpecDSMC(1)%EZeroPoint
  iPart = PEM%pNext(iPart)
END DO
vBulk(1:3) = vBulk(1:3)/nPart
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) Evib = Evib/nPart

u0ij = 0.0
u2 = 0.0
DO iLoop = 1, nPart
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  DO fillMa1 =1, 3
    DO fillMa2 =fillMa1, 3
      u0ij(fillMa1, fillMa2)= u0ij(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2)
    END DO
  END DO
  u2= u2 + vmag2
END DO
u2 = u2/nPart
u0ij = u0ij/nPart

CellTemp = Species(1)%MassIC * u2/(3.0*BoltzmannConst)
dens = nPart * Species(1)%MacroParticleFactor / GEO%Volume(iElem)

!pressure = dens*BoltzmannConst*CellTemp
!dynamicvis = 15. * SQRT(PI * Species(1)%MassIC * BoltzmannConst * SpecDSMC(1)%TrefVHS) &
!                                  / ( 8. * PI * SpecDSMC(1)%DrefVHS**2. * (2.- SpecDSMC(1)%omegaVHS) & 
!                                  * (3. - SpecDSMC(1)%omegaVHS) ) &
!                                  * ( CellTemp / SpecDSMC(1)%TrefVHS)**(SpecDSMC(1)%omegaVHS + 0.5)
!InnerDOF = 0.
!thermalconduct = 0.25 * (15. + 2. * InnerDOF) * dynamicvis 
!cp = (5.+InnerDOF)/2.
!Prandl = cp * dynamicvis / thermalconduct
!relaxfreq = cp*pressure/thermalconduct !Prandl*pressure/dynamic viscousity
IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  IF(SpecDSMC(1)%PolyatomicMol) THEN
    ! Calculation of the vibrational degree of freedom for the particle 
      print*, 'Not implemented yet!'
      STOP
  ELSE
    InnerDOF = 2. !Rot
    TVib=Evib / (BoltzmannConst*SpecDSMC(1)%CharaTVib)
    IF (TVib.GT.0.0) THEN
      Tvib= SpecDSMC(1)%CharaTVib/LOG(1. + 1./(TVib))
      Xi_vib = 2.* Evib / (BoltzmannConst*Tvib)
      InnerDOF = InnerDOF + Xi_vib
    END IF
  END IF
ELSE
  InnerDOF = 0.
END IF
dynamicvis2 = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*SpecDSMC(1)%TrefVHS/Pi) &
        /(4.*(4.- 2.*SpecDSMC(1)%omegaVHS) * (6. - 2.*SpecDSMC(1)%omegaVHS)* SpecDSMC(1)%DrefVHS**2.)
Prandl =2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
relaxfreq = Prandl * dens*BoltzmannConst*SpecDSMC(1)%TrefVHS**(SpecDSMC(1)%omegaVHS + 0.5) &
    /dynamicvis2*CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)

CALL RANDOM_NUMBER(iRan)
nRelax = INT(nPart*(1.-exp(-relaxfreq*dt)))
ProbAddPart = nPart*(1.-exp(-relaxfreq*dt)) - REAL(INT(nPart*(1.-exp(-relaxfreq*dt))))
IF (ProbAddPart.GT.iRan) nRelax = nRelax + 1
IF (nRelax.EQ.1) THEN
  CALL RANDOM_NUMBER(iRan)
  IF (iRan.GT.0.5) THEN
    nRelax = 2
  ELSE
    nRelax = 0
  END IF
END IF
IF (nRelax.EQ.0) RETURN

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  collisionfreq = SpecESBGK(1)%CollFreqPreFactor(1) * dens *CellTemp**(-SpecDSMC(1)%omegaVHS +0.5)
  rotrelaxfreq = collisionfreq * DSMC%RotRelaxProb
  vibrelaxfreq = collisionfreq * DSMC%RotRelaxProb * DSMC%VibRelaxProb
  ALLOCATE(DoRotRelax(nRelax), DoVibRelax(nRelax))
  DoRotRelax = .FALSE.
  DoVibRelax = .FALSE.

  CALL RANDOM_NUMBER(iRan)
  nRelaxRot = INT(nPart*(1.-exp(-rotrelaxfreq*dt)))
  ProbAddPart = nPart*(1.-exp(-rotrelaxfreq*dt)) - REAL(INT(nPart*(1.-exp(-rotrelaxfreq*dt))))
  IF (ProbAddPart.GT.iRan) nRelaxRot = nRelaxRot + 1
  CALL RANDOM_NUMBER(iRan)
  nRelaxVib = INT(nPart*(1.-exp(-vibrelaxfreq*dt)))
  ProbAddPart = nPart*(1.-exp(-vibrelaxfreq*dt)) - REAL(INT(nPart*(1.-exp(-vibrelaxfreq*dt))))
  IF (ProbAddPart.GT.iRan) nRelaxVib = nRelaxVib + 1
  IF ((nRelaxRot.GT.nRelax).OR.(nRelaxVib.GT.nRelax)) THEN
    print*, 'nRelax inner > nrelax', nRelax, nRelaxRot, nRelaxVib
    STOP
  END IF
END IF

DO fillMa1 =1, 3
  DO fillMa2 =fillMa1, 3
    IF (fillMa1.EQ.fillMa2) THEN
      KronDelta = 1.0
    ELSE
      KronDelta = 0.0
    END IF
!    SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandl)/Prandl*(3.*u0ij(fillMa1, fillMa2)/(u0ij(1,1)+u0ij(2,2)+u0ij(3,3))-KronDelta) 
    SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandl)/(2.*Prandl) &
      *(Species(1)%MassIC/(BoltzmannConst*CellTemp)*(u0ij(fillMa1, fillMa2)-vBulk(fillMa1)*vBulk(fillMa2))-KronDelta)
  END DO
END DO
SMat(2,1)=SMat(1,2)
SMat(3,1)=SMat(1,3)
SMat(3,2)=SMat(2,3)

ALLOCATE(iRanPart(3,nRelax),iPartIndx_NodeRelax(nRelax))
CALL ESBGK_BuildTransGaussNums(nRelax, iRanPart)

vBulk(1:3) = 0.0
NewEn = 0.

testMomNew = 0.0; testMomOld = 0.0; testEnNew = 0.0; testEnOld = 0.0

OldEn = 0.; newnRelaxRot = 0; newnRelaxVib=0
DO iLoop = 1, nRelax
  CALL RANDOM_NUMBER(iRan)
  iPart = 1 + INT(nPart * iRan)
  iPartIndx_NodeRelax(iLoop) = iPartIndx_Node(iPart)
  iPartIndx_Node(iPart) =  iPartIndx_Node(nPart)
  nPart = nPart - 1
  vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_NodeRelax(iLoop),4:6)

  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    OldEn = OldEn + 0.5*Species(1)%MassIC*(PartState(iPartIndx_NodeRelax(iLoop),4)**2.0 &
          + PartState(iPartIndx_NodeRelax(iLoop),5)**2.0 + PartState(iPartIndx_NodeRelax(iLoop),6)**2.0)
    CALL RANDOM_NUMBER(iRan)
    IF (nRelaxRot/nRelax.GT.iRan) THEN
      newnRelaxRot = newnRelaxRot + 1
      OldEn = OldEn + PartStateIntEn(iPartIndx_NodeRelax(iLoop),2)
      DoRotRelax(iLoop) = .TRUE.
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF (nRelaxVib/nRelax.GT.iRan) THEN
      newnRelaxVib = newnRelaxVib + 1
      OldEn = OldEn + PartStateIntEn(iPartIndx_NodeRelax(iLoop),1) - SpecDSMC(1)%EZeroPoint
      DoVibRelax(iLoop) = .TRUE.
    END IF
    testEnOld = testEnOld + PartStateIntEn(iPartIndx_NodeRelax(iLoop),1) + PartStateIntEn(iPartIndx_NodeRelax(iLoop),2)
  END IF

  testMomOld(1:3) = testMomOld(1:3) + PartState(iPartIndx_NodeRelax(iLoop),4:6)
  testEnOld = testEnOld + (PartState(iPartIndx_NodeRelax(iLoop),4)**2. + PartState(iPartIndx_NodeRelax(iLoop),5)**2. &
          + PartState(iPartIndx_NodeRelax(iLoop),6)**2.)
END DO
vBulk = vBulk / nRelax

IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
  TEqui = CalcTEquilib(OldEn, nRelax, newnRelaxRot, newnRelaxVib, 1)
  OldEn = 0.
  DO iLoop = 1, nRelax
    IF (DoVibRelax(iLoop)) THEN
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
      DO WHILE (iQuant.GE.SpecDSMC(1)%MaxVibQuant)
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(1)%CharaTVib)
      END DO
      OldEn = OldEn + PartStateIntEn(iPartIndx_NodeRelax(iLoop), 1) - SpecDSMC(1)%EZeroPoint
      PartStateIntEn(iPartIndx_NodeRelax(iLoop), 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(1)%CharaTVib*BoltzmannConst
      OldEn = OldEn - PartStateIntEn(iPartIndx_NodeRelax(iLoop), 1) - SpecDSMC(1)%EZeroPoint
    END IF
    IF (DoRotRelax(iLoop)) THEN
      CALL RANDOM_NUMBER(iRan)
      OldEn = OldEn + PartStateIntEn(iPartIndx_NodeRelax(iLoop), 2)
      PartStateIntEn(iPartIndx_NodeRelax(iLoop), 2) = -BoltzmannConst*TEqui*LOG(iRan)
      OldEn = OldEn - PartStateIntEn(iPartIndx_NodeRelax(iLoop), 2)
    END IF
  END DO
ELSE
  TEqui = CellTemp
  OldEn = 0.
END IF

DO iLoop = 1, nRelax
  tempVelo(1:3) = SQRT(BoltzmannConst*TEqui/Species(1)%MassIC)*iRanPart(1:3,iLoop)
  DO fillMa1 = 1, 3
    DO fillMa2 = 1, 3
      DSMC_RHS(iPartIndx_NodeRelax(iLoop),fillMa1) = DSMC_RHS(iPartIndx_NodeRelax(iLoop),fillMa1) &
            + SMat(fillMa1, fillMa2)*tempVelo(fillMa2)
    END DO
  END DO
  NewEn = NewEn + (DSMC_RHS(iPartIndx_NodeRelax(iLoop),1)**2. + DSMC_RHS(iPartIndx_NodeRelax(iLoop),2)**2. &
        + DSMC_RHS(iPartIndx_NodeRelax(iLoop),3)**2.)
  V_rel(1:3)=PartState(iPartIndx_NodeRelax(iLoop),4:6)-vBulk(1:3)
  OldEn = OldEn + (V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2)
END DO

alpha = SQRT(OldEn/NewEn)                       ! create particle index list for pairing
alpha = 1.
DO iLoop = 1, nRelax
  DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) = vBulk(1:3) + alpha*DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3) &
                      - PartState(iPartIndx_NodeRelax(iLoop),4:6)
  testMomNew(1:3) = testMomNew(1:3) +vBulk(1:3) + alpha*DSMC_RHS(iPartIndx_NodeRelax(iLoop),1:3)
  testEnNew = testEnNew + ((vBulk(1) + alpha*DSMC_RHS(iPartIndx_NodeRelax(iLoop),1))**2. &
          + (vBulk(2) + alpha*DSMC_RHS(iPartIndx_NodeRelax(iLoop),2))**2. &
          + (vBulk(3) + alpha*DSMC_RHS(iPartIndx_NodeRelax(iLoop),3))**2.)
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20))  testEnNew = testEnNew &
    + PartStateIntEn(iPartIndx_NodeRelax(iLoop),1) + PartStateIntEn(iPartIndx_NodeRelax(iLoop),2)
END DO

print*, testEnOld, testEnNew, testEnOld-testEnNew
print*, testMomOld
print*, testMomNew
print*, testMomOld-testMomNew
read*

END SUBROUTINE ESBGK_CollisionOperator


REAL FUNCTION CalcTEquilib(ETotal, Ntrans, Nrot, Nvib, iSpec)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: ETotal  ! Charak TVib, mean vibrational Energy of all molecules
INTEGER, INTENT(IN)             :: iSpec, Ntrans, Nrot, Nvib      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iDOF,iPolyatMole
REAL(KIND=8)                    :: LowerTemp, UpperTemp, MiddleTemp ! upper and lower value of modified zero point search
REAl(KIND=8)                    :: eps_prec=1.0E-5
REAL(KIND=8)                    :: Etemp    ! both summs
REAL                            :: Xi_Vib
!===================================================================================================================================
!  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
!    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!    IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN

LowerTemp = 2.*ETotal / (BoltzmannConst*(3.*Ntrans + 2.*Nrot + 2.*Nvib))
UpperTemp = 2.*ETotal / (3.*BoltzmannConst*Ntrans)
DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
  MiddleTemp = 0.5*( LowerTemp + UpperTemp)
  Xi_Vib = (2.0*SpecDSMC(iSpec)%CharaTVib / MiddleTemp) / (EXP(SpecDSMC(iSpec)%CharaTVib / MiddleTemp) - 1.0)
  Etemp = (3.*Ntrans + 2.*Nrot + Xi_Vib*Nvib) / 2.* BoltzmannConst * MiddleTemp
  IF ( Etemp .GT. ETotal) THEN
    UpperTemp = MiddleTemp
  ELSE
    LowerTemp = MiddleTemp
  END IF
END DO
CalcTEquilib = UpperTemp ! or 0.5*( Tmax + Tmin)
RETURN
END FUNCTION CalcTEquilib

END MODULE MOD_ESBGK_MOLEC
