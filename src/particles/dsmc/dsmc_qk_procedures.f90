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

MODULE MOD_DSMC_QK_Chemistry
!===================================================================================================================================
! module including qk procedures
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! Public Part
PUBLIC :: QK_Init, QK_TestReaction, QK_CalcAnalyticRate, QK_GetAnalyticRate
!===================================================================================================================================
CONTAINS

SUBROUTINE QK_Init()
!===================================================================================================================================
!> Calculate the collision constant per collision case (for every reaction, as recombination probability through equilibrium
!> constant also utilizes RColl) and the analytical QK reaction rate
!===================================================================================================================================
! MODULES
USE MOD_Globals         ,ONLY: abort,StringBeginsWith
USE MOD_Globals_Vars    ,ONLY: Pi, BoltzmannConst
USE MOD_DSMC_Vars       ,ONLY: DSMC, ChemReac, CollInf, QKChemistry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec1, iSpec2, iCase, iReac, iInter, PartitionArraySize
REAL                          :: omega, Tref, dref, TCollExponent, Temp
!===================================================================================================================================
ALLOCATE(ChemReac%QKRColl(CollInf%NumCase),ChemReac%QKTCollCorrFac(CollInf%NumCase))
ChemReac%QKRColl = 0.
ChemReac%QKTCollCorrFac = 0.

DO iReac = 1, ChemReac%NumOfReact
  ! Skip the special case of photo ionization
  IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) CYCLE
  iSpec1 = ChemReac%Reactants(iReac,1)
  iSpec2 = ChemReac%Reactants(iReac,2)
  iCase = CollInf%Coll_Case(iSpec1, iSpec2)
  Tref = CollInf%Tref(iSpec1,iSpec2)
  omega = CollInf%omega(iSpec1,iSpec2)
  dref = CollInf%dref(iSpec1,iSpec2)
  ! Reduced mass is omitted; a weighted mass is required in case of an axisymmetric simulation and included in CalcReactionProb
  ChemReac%QKRColl(iCase) = 2. * SQRT(Pi) / (1 + CollInf%KronDelta(iCase)) * (dref)**2 * SQRT(2. * BoltzmannConst * Tref)
  TCollExponent = (0.5 - omega)
  ChemReac%QKTCollCorrFac(iCase) = (2. - omega)**TCollExponent * GAMMA(2. - omega) / GAMMA(2. - omega + TCollExponent)
END DO


IF(ChemReac%AnyQKReaction) THEN
  IF(MOD(DSMC%PartitionMaxTemp,DSMC%PartitionInterval).EQ.0.0) THEN
    PartitionArraySize = NINT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)
  ELSE
    CALL abort(&
      __STAMP__&
      ,'ERROR in Chemistry Init: Partition temperature limit must be multiple of partition interval!')
  END IF
  ALLOCATE(QKChemistry(ChemReac%NumOfReact))
  DO iReac = 1, ChemReac%NumOfReact
    IF(TRIM(ChemReac%ReactModel(iReac)).EQ.'QK')THEN
      ! Calculation of the analytical rate, to be able to calculate the backward rate with partition function
      ALLOCATE(QKChemistry(iReac)%ForwardRate(1:PartitionArraySize))
      DO iInter = 1, PartitionArraySize
        Temp = iInter * DSMC%PartitionInterval
        QKChemistry(iReac)%ForwardRate(iInter) = QK_CalcAnalyticRate(iReac,Temp)
      END DO
      ! For backward rates, the corresponding forward rate is stored to be used with the equilibrium constant
      IF(iReac.GT.ChemReac%NumOfReactWOBackward) THEN
        QKChemistry(iReac)%ForwardRate = QKChemistry(ChemReac%BackwardReacForwardIndx(iReac))%ForwardRate
      END IF
    END IF
  END DO
END IF

END SUBROUTINE QK_Init


SUBROUTINE QK_TestReaction(iPair,iReac,PerformReaction)
!===================================================================================================================================
! Decide whether a dissociation/ionization occurs: Check if the relative translation kinetic energy plus vibrational/electronic
! energy of the molecule under consideration is larger than the dissociation/ionization energy
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac, DSMC, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
LOGICAL, INTENT(INOUT)        :: PerformReaction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: MaxQua, React1Inx, React2Inx, iCase, iPath, PathIndex
REAL                          :: Weight1, Weight2, ReducedMass, IonizationEnergy, Ec
!===================================================================================================================================
! Determining, which collision partner is the dissociating particle (always the first molecule in the DefinedReact array)
IF (ChemReac%Reactants(iReac,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  React1Inx = Coll_pData(iPair)%iPart_p1
  React2Inx = Coll_pData(iPair)%iPart_p2
ELSE
  React1Inx = Coll_pData(iPair)%iPart_p2
  React2Inx = Coll_pData(iPair)%iPart_p1
END IF

iCase = CollInf%Coll_Case(PartSpecies(React1Inx), PartSpecies(React2Inx))
Weight1 = GetParticleWeight(React1Inx)
Weight2 = GetParticleWeight(React2Inx)

! Determine the number of the reaction path for this collision pair with the help of the reaction index
DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  IF(iReac.EQ.ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)) PathIndex = iPath
END DO

IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  ReducedMass = (Species(PartSpecies(React1Inx))%MassIC*Weight1 * Species(PartSpecies(React2Inx))%MassIC*Weight2) &
              / (Species(PartSpecies(React1Inx))%MassIC*Weight1 + Species(PartSpecies(React2Inx))%MassIC*Weight2)
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF

! Determine the collision energy (relative translational + vibrational energy of dissociating molecule)
Ec = 0.5 * ReducedMass*Coll_pData(iPair)%CRela2

SELECT CASE(TRIM(ChemReac%ReactType(iReac)))
CASE('I')
  IF(DSMC%ElectronicModel.GT.0) Ec = Ec + PartStateIntEn(3,React1Inx)*Weight1
  ! ionization level is last known energy level of species
  MaxQua=SpecDSMC(PartSpecies(React1Inx))%MaxElecQuant - 1
  IonizationEnergy=SpecDSMC(PartSpecies(React1Inx))%ElectronicState(2,MaxQua)*BoltzmannConst*(Weight1 + Weight2)/2.
  ! if you have electronic levels above the ionization limit, such limits should be used instead of
  ! the pure energy comparison
  IF(Ec .GT. IonizationEnergy) THEN
    PerformReaction = .TRUE.
  END IF
CASE('D')
  Ec = Ec + PartStateIntEn(1,React1Inx)*Weight1
  ! Correction for second collision partner
  IF ((SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.20)) THEN
    Ec = Ec - SpecDSMC(PartSpecies(React2Inx))%EZeroPoint*Weight2
  END IF
  ! Determination of the quantum number corresponding to the collision energy
  MaxQua = INT(Ec / ((Weight1 + Weight2)/2.*BoltzmannConst * SpecDSMC(PartSpecies(React1Inx))%CharaTVib ) - DSMC%GammaQuant)
  ! Comparing the collision quantum number with the dissociation quantum number
  IF (MaxQua.GT.SpecDSMC(PartSpecies(React1Inx))%DissQuant) THEN
    PerformReaction = .TRUE.
  END IF
CASE DEFAULT
  CALL Abort(__STAMP__,&
    'Given reaction type is not supported with the Quantum-Kinetic (QK) model. Reaction: ',iReac)
END SELECT

END SUBROUTINE QK_TestReaction

REAL FUNCTION QK_CalcAnalyticRate(iReac,Temp)
!===================================================================================================================================
! Calculation of the forward reaction rate through the analytical expression for QK
!===================================================================================================================================
! MODULES
USE MOD_Globals        ,ONLY: abort
USE MOD_DSMC_Vars      ,ONLY: SpecDSMC, ChemReac, CollInf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iReac
REAL, INTENT(IN)              :: Temp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec1, iSpec2, MaxElecQua, iQua, iCase
REAL                          :: z ! contribution of the relevant mode to the electronic or vibrational partition function
REAL                          :: Q ! incomplete gamma function
INTEGER                       :: MaxVibQuant ! highest vibrational quantum state
REAL                          :: Tref, TempRatio, omega
!===================================================================================================================================
QK_CalcAnalyticRate = 0.0
z = 0.0
iSpec1 = ChemReac%Reactants(iReac,1)
iSpec2 = ChemReac%Reactants(iReac,2)
iCase = CollInf%Coll_Case(iSpec1, iSpec2)
Tref = CollInf%Tref(iSpec1,iSpec2)
omega = CollInf%omega(iSpec1,iSpec2)

SELECT CASE (ChemReac%ReactType(iReac))
CASE('I')
  MaxElecQua=SpecDSMC(iSpec1)%MaxElecQuant - 1
  DO iQua = 0, MaxElecQua
    Q = gammainc([2.-omega,(SpecDSMC(iSpec1)%ElectronicState(2,MaxElecQua)- &
        SpecDSMC(iSpec1)%ElectronicState(2,iQua))/Temp])
    TempRatio = SpecDSMC(iSpec1)%ElectronicState(2,iQua) / Temp
    IF(CHECKEXP(TempRatio)) THEN
      QK_CalcAnalyticRate= QK_CalcAnalyticRate + Q * SpecDSMC(iSpec1)%ElectronicState(1,iQua) * EXP(-TempRatio)
      z = z + SpecDSMC(iSpec1)%ElectronicState(1,iQua) * EXP(-TempRatio)
    END IF
  END DO
  QK_CalcAnalyticRate = QK_CalcAnalyticRate*(Temp / Tref)**(0.5 - omega)*ChemReac%QKRColl(iCase) &
                        / (z * SQRT(CollInf%MassRed(iCase)))
CASE('D')
  MaxVibQuant = SpecDSMC(iSpec1)%DissQuant
  TempRatio = SpecDSMC(iSpec1)%CharaTVib / Temp
  DO iQua = 0, MaxVibQuant - 1
    Q = gammainc([2.-omega,((MaxVibQuant-iQua)*SpecDSMC(iSpec1)%CharaTVib)/Temp])
    IF(CHECKEXP(iQua*TempRatio)) THEN
      QK_CalcAnalyticRate= QK_CalcAnalyticRate + Q * EXP(-iQua*TempRatio)
    END IF
  END DO
  IF(CHECKEXP(TempRatio)) THEN
    z = 1. / (1. - EXP(-TempRatio))
  ELSE
    z = 1.
  END IF
  QK_CalcAnalyticRate = QK_CalcAnalyticRate*(Temp / Tref)**(0.5 - omega)*ChemReac%QKRColl(iCase) &
                        / (z * SQRT(CollInf%MassRed(iCase)))
END SELECT

END FUNCTION QK_CalcAnalyticRate


REAL FUNCTION QK_GetAnalyticRate(iReac,LocalTemp)
!===================================================================================================================================
! Interpolate the analytic QK reaction rate from the stored QKChemistry array
! CalcBackwardRate: for the backward rate, the value of the respective forward rate was copied during the initialization
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,              ONLY: DSMC, QKChemistry
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iReac
  REAL, INTENT(IN)              :: LocalTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: LowerLevel, UpperLevel
!===================================================================================================================================
! Determination of the lower and upper value of the temperature interval
LowerLevel = INT(LocalTemp/DSMC%PartitionInterval)
UpperLevel = LowerLevel + 1

IF((UpperLevel.GT.INT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)).OR.(LowerLevel.EQ.0)) THEN
! Instantaneous calculation of the reaction rate
  QK_GetAnalyticRate = QK_CalcAnalyticRate(iReac,LocalTemp)
ELSE
! Linear interpolation of the backward rate coefficient at the actual temperature
  QK_GetAnalyticRate = QKChemistry(iReac)%ForwardRate(LowerLevel) &
              + (QKChemistry(iReac)%ForwardRate(UpperLevel) - QKChemistry(iReac)%ForwardRate(LowerLevel))  &
              / (DSMC%PartitionInterval) * (LocalTemp - LowerLevel * DSMC%PartitionInterval)
END IF

END FUNCTION QK_GetAnalyticRate


FUNCTION gammainc( arg )
!===================================================================================================================================
! Program to test the incomplete gamma function
! the following gamma function is the one of Birds Q-K rate code
! ev. take another gamma function implementation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,PARAMETER              :: real_kind=8
  REAL(KIND=real_kind),DIMENSION(1:2), INTENT(IN) :: arg
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: n
  REAL(KIND=real_kind)           :: gamser, gln, ap, del, summ, an, ser, tmp, x,y, b,c,d,h, exparg
  REAL(KIND=real_kind)           :: gammainc
  ! parameters
  REAL(KIND=real_kind),PARAMETER,DIMENSION(6) :: &
                                      cof= [ 76.18009172947146      , &
                                            -86.50532032941677     , &
                                             24.01409824083091      , &
                                             -1.231739572450155     , &
                                              0.1208650973866179e-2  , &
                                            -0.5395239384953e-5 ]
  REAL(KIND=real_kind),PARAMETER :: stp=2.5066282746310005        , &
                                    fpmin=1.e-30
!===================================================================================================================================

  x=arg(1)
  y=x
  tmp=x+5.5
  tmp=(x+0.5)*log(tmp)-tmp
  ser=1.000000000190015
  DO n = 1, 6
    y=y+1.
    ser=ser+cof(n)/y
  END DO
  gln=tmp+log(stp*ser/x)
  IF (arg(2) < arg(1)+1.) THEN
    IF (arg(2) <= 0.) THEN
      gamser=0.
    ELSE
      ap=arg(1)
      summ=1./arg(1)
      del=summ
      DO WHILE (abs(del) > abs(summ)*1.e-8 )
        ap=ap+1.
        del=del*arg(2)/ap
        summ=summ+del
      END DO
      gamser=summ*exp(-arg(2)+arg(1)*log(arg(2))-gln)
    END IF
    gammainc=1.-gamser
  ELSE
    b =arg(2)+1.-arg(1)
    c=1./fpmin
    d=1./b
    h=d
    del=d*c
    n=0
    DO WHILE ( abs(del-1.) >= 1.e-8 )
      n=n+1
      an=-n*(n-arg(1))
      b=b+2.
      d=an*d+b
      IF ( abs(d) < fpmin ) THEN
        d=fpmin
      END IF
      c=b+an/c
      IF ( abs(c) < fpmin ) THEN
        c=fpmin
      END IF
      d=1./d
      del=d*c
      h=h*del
    END DO
    exparg = -arg(2)+arg(1)*log(arg(2))-gln
    IF(CHECKEXP(exparg))THEN
      gammainc = exp(exparg) * h
    ELSE
      gammainc = 0.
    END IF ! CHECKEXP(exparg)
  END IF
END FUNCTION gammainc

END MODULE MOD_DSMC_QK_Chemistry
