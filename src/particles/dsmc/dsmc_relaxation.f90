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

MODULE MOD_DSMC_Relaxation
!===================================================================================================================================
! Module including collisions, relaxation and reaction decision
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_VibRelaxDiatomic, CalcMeanVibQuaDiatomic, CalcXiVib, CalcXiTotalEqui, DSMC_calc_P_rot, DSMC_calc_var_P_vib
PUBLIC :: InitCalcVibRelaxProb, DSMC_calc_P_vib, SumVibRelaxProb, FinalizeCalcVibRelaxProb, DSMC_calc_P_elec
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_VibRelaxDiatomic(iPair, iPart, FakXi)
!===================================================================================================================================
! Performs the vibrational relaxation of diatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, PartStateIntEn, Coll_pData, RadialWeighting
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart, iPair
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: MaxColQua, iRan, Ec
INTEGER                       :: iQuaMax, iQua
!===================================================================================================================================
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF

MaxColQua = Ec/(BoltzmannConst*SpecDSMC(PartSpecies(iPart))%CharaTVib) - DSMC%GammaQuant
iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(iPart))%MaxVibQuant)
CALL RANDOM_NUMBER(iRan)
iQua = INT(iRan * iQuaMax)
CALL RANDOM_NUMBER(iRan)
DO WHILE (iRan.GT.(1 - REAL(iQua)/REAL(MaxColQua))**FakXi)
  !laux diss page 31
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(iRan * iQuaMax)
  CALL RANDOM_NUMBER(iRan)
END DO

PartStateIntEn(1,iPart) = (iQua + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(PartSpecies(iPart))%CharaTVib

END SUBROUTINE DSMC_VibRelaxDiatomic


SUBROUTINE CalcMeanVibQuaDiatomic()
!===================================================================================================================================
! Computes the mean vibrational quantum number of diatomic species in a cell each iteration;
! ChemReac%MeanEVibQua_PerIter is required for the determination of the vibrational degree of freedom, only used for diatomic
! molecules. For the polyatomic case, the actual vibrational degree of freedom of each molecule is utilized and not a mean value.
! The values for polyatomic molecules can have a greater spread, thus a mean value can prohibit reactions of highly excited
! molecules at lower average temperatures.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : BoltzmannConst
USE MOD_DSMC_Vars,             ONLY : DSMC, CollInf, SpecDSMC, ChemReac, BGGas
USE MOD_Particle_Vars,         ONLY : nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSpec
REAL            :: iRan, VibQuaTemp
!===================================================================================================================================

DO iSpec = 1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(.NOT.SpecDSMC(iSpec)%PolyatomicMol) THEN
      ! Skip the background gas species (value initialized in dsmc_chemical_init.f90)
      IF(BGGas%BackgroundSpecies(iSpec)) CYCLE
      ! Only treat species present in the cell
      IF(CollInf%Coll_SpecPartNum(iSpec).GT.0.) THEN
        ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
        VibQuaTemp = ChemReac%MeanEVib_PerIter(iSpec) / (BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
        CALL RANDOM_NUMBER(iRan)
        IF((VibQuaTemp-INT(VibQuaTemp)).GT.iRan) THEN
          ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(VibQuaTemp) + 1, SpecDSMC(iSpec)%MaxVibQuant-1)
        ELSE
          ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(VibQuaTemp), SpecDSMC(iSpec)%MaxVibQuant-1)
        END IF
        IF(ChemReac%MeanEVibQua_PerIter(iSpec).GT.0) THEN
          ChemReac%MeanXiVib_PerIter(iSpec) = 2. * ChemReac%MeanEVibQua_PerIter(iSpec) &
                                            * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0 )
        ELSE
          ChemReac%MeanXiVib_PerIter(iSpec) = 0.
        END IF
      ELSE
        ChemReac%MeanEVibQua_PerIter(iSpec) = 0
        ChemReac%MeanXiVib_PerIter(iSpec) = 0.
      END IF  ! CollInf%Coll_SpecPartNum(iSpec).GT.0
    END IF    ! .NOT.SpecDSMC(iSpec)%PolyatomicMol
  END IF      ! (SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)
END DO        ! iSpec = 1, nSpecies

END SUBROUTINE CalcMeanVibQuaDiatomic


SUBROUTINE CalcXiVib(TVib, iSpec, XiVibDOF, XiVibTotal)
!===================================================================================================================================
! Calculation of the vibrational degrees of freedom for each characteristic vibrational temperature, used for chemical reactions
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: TVib  !
INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT),OPTIONAL      :: XiVibDOF(:), XiVibTotal
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iDOF, iPolyatMole, VibDOF
REAL                            :: TempRatio
REAL,ALLOCATABLE                :: XiVibPart(:)
!===================================================================================================================================

IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  VibDOF = PolyatomMolDSMC(iPolyatMole)%VibDOF
  ALLOCATE(XiVibPart(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  XiVibPart = 0.0
  DO iDOF = 1 , VibDOF
    ! If the temperature is very small compared to the characteristic temperature, set the vibrational degree of freedom to zero
    ! to avoid overflows in the exponential function
    TempRatio = PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/ TVib
    IF(CHECKEXP(TempRatio)) THEN
      XiVibPart(iDOF) = (2.0*TempRatio) / (EXP(TempRatio) - 1.0)
    END IF
  END DO
ELSE
  VibDOF = 1
  ALLOCATE(XiVibPart(VibDOF))
  XiVibPart = 0.0
  TempRatio = SpecDSMC(iSpec)%CharaTVib / TVib
  IF(CHECKEXP(TempRatio)) THEN
    XiVibPart(1) = (2.0*TempRatio) / (EXP(TempRatio) - 1.0)
  END IF
END IF

IF(PRESENT(XiVibDOF)) THEN
  XiVibDOF = 0.
  XiVibDOF(1:VibDOF) = XiVibPart(1:VibDOF)
END IF

IF(PRESENT(XiVibTotal)) THEN
  XiVibTotal = SUM(XiVibPart)
END IF

RETURN

END SUBROUTINE CalcXiVib


SUBROUTINE CalcXiTotalEqui(iReac, iPair, nProd, Xi_Total, Weight, XiVibPart, XiElecPart)
!===================================================================================================================================
! Calculation of the vibrational degrees of freedom for each characteristic vibrational temperature as well as the electronic
! degrees of freedom at a common temperature, used for chemical reactions
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars              ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars                 ,ONLY: SpecDSMC, ChemReac, Coll_pData, DSMC
USE MOD_part_tools                ,ONLY: CalcXiElec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iReac, iPair, nProd
REAL, INTENT(IN)                :: Xi_Total, Weight(1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT), OPTIONAL     :: XiVibPart(:,:), XiElecPart(1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iProd, iSpec
REAL                            :: ETotal, EZeroPoint, EGuess, LowerTemp, UpperTemp, MiddleTemp, Xi_TotalTemp, XiVibTotal
REAL,PARAMETER                  :: eps_prec=1E-3
!===================================================================================================================================

ASSOCIATE( ProductReac => ChemReac%Products(iReac,1:4) )

  ! Weighted total collision energy
  ETotal = Coll_pData(iPair)%Ec

  EZeroPoint = 0.0
  DO iProd = 1, nProd
    EZeroPoint = EZeroPoint + SpecDSMC(ProductReac(iProd))%EZeroPoint * Weight(iProd)
  END DO

  LowerTemp = 1.0
  UpperTemp = 2.*(ETotal - EZeroPoint) * nProd / SUM(Weight) / (Xi_Total * BoltzmannConst)
  MiddleTemp = LowerTemp
  DO WHILE (.NOT.ALMOSTEQUALRELATIVE(0.5*(LowerTemp + UpperTemp),MiddleTemp,eps_prec))
    MiddleTemp = 0.5*( LowerTemp + UpperTemp)
    Xi_TotalTemp = Xi_Total
    DO iProd = 1, nProd
      iSpec = ProductReac(iProd)
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        CALL CalcXiVib(MiddleTemp, iSpec, XiVibDOF=XiVibPart(iProd,:), XiVibTotal=XiVibTotal)
        Xi_TotalTemp = Xi_TotalTemp + XiVibTotal
      ELSE
        IF(PRESENT(XiVibPart)) XiVibPart(iProd,:) = 0.0
      END IF
      IF((DSMC%ElectronicModel.EQ.1).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.4)) THEN
        IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          XiElecPart(iProd) = CalcXiElec(MiddleTemp, iSpec)
          Xi_TotalTemp = Xi_TotalTemp + XiElecPart(iProd)
        ELSE
          IF(PRESENT(XiElecPart)) XiElecPart(iProd) = 0.0
        END IF
      END IF
    END DO
    EGuess = EZeroPoint + Xi_TotalTemp / 2. * BoltzmannConst * MiddleTemp * SUM(Weight) / nProd
    IF (EGuess .GT. ETotal) THEN
      UpperTemp = MiddleTemp
    ELSE
      LowerTemp = MiddleTemp
    END IF
  END DO
END ASSOCIATE

RETURN

END SUBROUTINE CalcXiTotalEqui


SUBROUTINE InitCalcVibRelaxProb()
!===================================================================================================================================
! Initialize the calculation of the variable vibrational relaxation probability in the cell for each iteration
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, VarVibRelaxProb 
USE MOD_Particle_Vars      ,ONLY: nSpecies

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iSpec
!===================================================================================================================================

IF(DSMC%VibRelaxProb.EQ.2.0) THEN ! Set summs for variable vibrational relaxation to zero
  DO iSpec=1,nSpecies
    VarVibRelaxProb%ProbVibAvNew(iSpec) = 0
    VarVibRelaxProb%nCollis(iSpec) = 0
  END DO
END IF

END SUBROUTINE InitCalcVibRelaxProb


SUBROUTINE DSMC_calc_P_rot(iSpec1, iSpec2, iPair, iPart, Xi_rel, ProbRot, ProbRotMax)
!===================================================================================================================================
! Calculation of probability for rotational relaxation. Different Models implemented:
! 0 - Constant Probability
! 1 - No rotational relaxation. RotRelaxProb = 0
! 2 - Boyd
! 3 - Zhang (Nonequilibrium Direction Dependent)
!===================================================================================================================================
! MODULES
USE MOD_Globals            ,ONLY : Abort
USE MOD_Globals_Vars       ,ONLY : Pi, BoltzmannConst
USE MOD_DSMC_Vars          ,ONLY : SpecDSMC, Coll_pData, PartStateIntEn, DSMC, useRelaxProbCorrFactor, CollInf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec1, iSpec2, iPair, iPart
REAL, INTENT(IN)            :: Xi_rel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbRot, ProbRotMax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: TransEn, RotEn, RotDOF, CorrFact           ! CorrFact: To correct sample Bias
                                                                        ! (fewer DSMC particles than natural ones)
!===================================================================================================================================

TransEn    = Coll_pData(iPair)%Ec ! notice that during probability calculation,Collision energy only contains translational part
RotDOF     = SpecDSMC(iSpec1)%Xi_Rot
RotEn      = PartStateIntEn(2,iPart)
ProbRot    = 0.
ProbRotMax = 0.

! calculate correction factor according to Lumpkin et al.
! - depending on selection procedure. As only one particle undergoes relaxation
! - only one RotDOF is needed (of considered species)
IF(useRelaxProbCorrFactor) THEN
  CorrFact = 1. + RotDOF/Xi_rel
ELSE
  CorrFact = 1.
END IF

! calculate corrected probability for rotational relaxation
IF(DSMC%RotRelaxProb.GE.0.0.AND.DSMC%RotRelaxProb.LE.1.0) THEN
  ProbRot = DSMC%RotRelaxProb * CorrFact
ELSEIF(DSMC%RotRelaxProb.EQ.2.0) THEN ! P_rot according to Boyd (based on Parker's model)

  RotDOF = RotDOF*0.5 ! Only half of the rotational degree of freedom, because the other half is used in the relaxation
                      ! probability of the collision partner, see Boyd (doi:10.1063/1.858531)

  ProbRot = 1./SpecDSMC(iSpec1)%CollNumRotInf * (1. + GAMMA(RotDOF+2.-CollInf%omega(iSpec1,iSpec2)) &
          / GAMMA(RotDOF+1.5-CollInf%omega(iSpec1,iSpec2)) * (PI**(3./2.)/2.)*(BoltzmannConst*SpecDSMC(iSpec1)%TempRefRot &
          / (TransEn + RotEn) )**(1./2.) + GAMMA(RotDOF+2.-CollInf%omega(iSpec1,iSpec2))  &
          / GAMMA(RotDOF+1.-CollInf%omega(iSpec1,iSpec2)) * (BoltzmannConst*SpecDSMC(iSpec1)%TempRefRot &
          / (TransEn + RotEn) ) * (PI**2./4. + PI)) &
          * CorrFact

ELSEIF(DSMC%RotRelaxProb.EQ.3.0) THEN ! P_rot according to Zhang (NDD)
  ! if model is used for further species but N2, it should be checked if factors n = 0.5 and Cn = 1.92 are still valid
  ! (see original eq of Zhang)
  ProbRot = 1.92 * GAMMA(Xi_rel/2.) * GAMMA(RotDOF/2.) / GAMMA(Xi_rel/2.+0.5) / GAMMA(RotDOF/2.-0.5) &
          * (1 + (Xi_rel/2-0.5)*BoltzmannConst*SpecDSMC(iSpec1)%TempRefRot/TransEn) * (TransEn/RotEn)**0.5 &
          * CorrFact
  ProbRotMax = MAX(ProbRot, 0.5) ! BL energy redistribution correction factor
  ProbRot    = MIN(ProbRot, 0.5)
ELSE
  CALL Abort(__STAMP__,'Error! Model for rotational relaxation undefined:',RealInfoOpt=DSMC%RotRelaxProb)
END IF

END SUBROUTINE DSMC_calc_P_rot


SUBROUTINE DSMC_calc_P_elec(iSpec1, iSpec2, ProbElec)
!===================================================================================================================================
! Calculation of probability for electronic relaxation. 
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars          ,ONLY : SpecDSMC, useRelaxProbCorrFactor, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec1, iSpec2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbElec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: CorrFact
!===================================================================================================================================
IF(useRelaxProbCorrFactor.AND.(DSMC%ElectronicModel.EQ.1)) THEN
  CorrFact = SpecDSMC(iSpec1)%ElecRelaxCorrectFac(iSpec2)
ELSE
  CorrFact = 1.
END IF
ProbElec = SpecDSMC(iSpec1)%ElecRelaxProb*CorrFact

END SUBROUTINE DSMC_calc_P_elec



SUBROUTINE DSMC_calc_P_vib(iPair, iSpec, jSpec, Xi_rel, iElem, ProbVib)
!===================================================================================================================================
! Calculation of probability for vibrational relaxation. Different Models implemented:
! 0 - Constant Probability
! 1 - No vibrational relaxation. VibRelaxProb = 0
! 2 - Boyd with correction of Abe
!===================================================================================================================================
! MODULES
USE MOD_Globals            ,ONLY: Abort
USE MOD_DSMC_Vars          ,ONLY: SpecDSMC, DSMC, VarVibRelaxProb, useRelaxProbCorrFactor, PolyatomMolDSMC, CollInf, Coll_pData
USE MOD_MCC_Vars           ,ONLY: XSec_Relaxation, SpecXSec
USE MOD_MCC_XSec           ,ONLY: XSec_CalcVibRelaxProb
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iPair, iSpec, jSpec, iElem
REAL, INTENT(IN)          :: Xi_rel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: CorrFact       ! CorrFact: To correct sample Bias
                                            ! (fewer DSMC particles than natural ones)
INTEGER                   :: iPolyatMole, iDOF, iCase
!===================================================================================================================================

ProbVib = 0.

! calculate correction factor according to Gimelshein et al.
! - depending on selection procedure. As only one particle undergoes relaxation
! - only one VibDOF (GammaVib) is needed (of considered species)
IF(useRelaxProbCorrFactor) THEN
  CorrFact = 1. + SpecDSMC(iSpec)%GammaVib/Xi_rel
ELSE
  CorrFact = 1.
END IF

IF((DSMC%VibRelaxProb.GE.0.0).AND.(DSMC%VibRelaxProb.LE.1.0)) THEN
  IF (SpecDSMC(iSpec)%PolyatomicMol.AND.(DSMC%PolySingleMode)) THEN
    iPolyatMole  = SpecDSMC(iSpec)%SpecToPolyArray
    PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(1) = DSMC%VibRelaxProb * (1. + PolyatomMolDSMC(iPolyatMole)%GammaVib(1)/Xi_rel)
    DO iDOF = 2, PolyatomMolDSMC(iPolyatMole)%VibDOF
      PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF) = PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF - 1) + DSMC%VibRelaxProb &
                                                      * (1. + PolyatomMolDSMC(iPolyatMole)%GammaVib(1)/Xi_rel)
    END DO
  ELSE
    ProbVib = DSMC%VibRelaxProb * CorrFact
  END IF
  IF(XSec_Relaxation) THEN
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    IF(SpecXSec(iCase)%UseVibXSec) THEN
      IF(SpecXSec(iCase)%SpeciesToRelax.EQ.iSpec) THEN
        IF(SpecXSec(iCase)%UseCollXSec) THEN
          CALL XSec_CalcVibRelaxProb(iPair,iElem)
          ! Cross-section is stored in the VibProb variable
          ProbVib = SpecXSec(iCase)%VibProb / SpecXSec(iCase)%CrossSection
        ELSE
          ProbVib = SpecXSec(iCase)%VibProb / Coll_pData(iPair)%Prob
        END IF
      END IF
    END IF
  END IF
ELSE IF(DSMC%VibRelaxProb.EQ.2.0) THEN
  ! Calculation of Prob Vib in function DSMC_calc_var_P_vib.
  ! This has to average over all collisions according to Boyd (doi:10.1063/1.858495)
  ! The average value of the cell is only taken from the vector
  ProbVib = VarVibRelaxProb%ProbVibAv(iElem, iSpec) * CorrFact
ELSE
  CALL Abort(__STAMP__,'Error! Model for vibrational relaxation undefined:',RealInfoOpt=DSMC%VibRelaxProb)
END IF

IF(DSMC%CalcQualityFactors) THEN
  DSMC%CalcVibProb(iSpec,1) = DSMC%CalcVibProb(iSpec,1) + ProbVib
  DSMC%CalcVibProb(iSpec,3) = DSMC%CalcVibProb(iSpec,3) + 1
END IF

END SUBROUTINE DSMC_calc_P_vib


SUBROUTINE DSMC_calc_var_P_vib(iSpec, jSpec, iPair, ProbVib)
!===================================================================================================================================
! Calculation of probability for vibrational relaxation for variable relaxation rates. This has to average over all collisions!
! No instantanious variable probability calculateable
!===================================================================================================================================
! MODULES
USE MOD_Globals      ,ONLY: Abort
USE MOD_Globals_Vars ,ONLY: Pi, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, Coll_pData, CollInf
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iPair, iSpec, jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: TempCorr, cRela
!===================================================================================================================================
! (i) dref changed from   DrefVHS = 0.5 * (SpecDSMC(iSpec)%DrefVHS + SpecDSMC(jSpec)%DrefVHS)
!                  to   dref(iSpec,jSpec) which is identical to old definition (for averagedCollisionParameters=TRUE (DEFAULT))
! in case of averagedCollisionParameter=FALSE dref(iSpec,jSpec) contains collision specific dref see --help for details

! P_vib according to Boyd, corrected by Abe, only V-T transfer
! determine joint omega and Dref factor and rel velo
cRela=SQRT(Coll_pData(iPair)%cRela2)
! calculate non-corrected probabilities
ProbVib = 1. /SpecDSMC(iSpec)%CollNumVib(jSpec)* cRela**(3.+2.*CollInf%omega(iSpec,jSpec)) &
        * EXP(-1.*SpecDSMC(iSpec)%CharaVelo(jSpec)/cRela)
! calculate high temperature correction
TempCorr = SpecDSMC(iSpec)%VibCrossSec / (SQRT(2.)*PI*CollInf%dref(iSpec,jSpec)**2.) &
         * (  CollInf%MassRed(Coll_pData(iPair)%PairType)*cRela & !**2
         / (2.*(2.-CollInf%omega(iSpec,jSpec))*BoltzmannConst*CollInf%Tref(iSpec,jSpec)))**CollInf%omega(iSpec,jSpec)
! determine corrected probabilities
ProbVib = ProbVib * TempCorr / (ProbVib + TempCorr)        ! TauVib = TauVibStd + TauTempCorr
IF(ProbVib.NE.ProbVib) THEN !If is NAN
  ProbVib=0.
  WRITE(*,*) 'WARNING: Vibrational relaxation probability is NAN and is set to zero. cRela:', cRela
  ! CALL Abort(&
  ! __STAMP__&
  ! ,'Error! Vibrational relaxation probability is NAN (cRela);',RealInfoOpt=cRela)!, jSpec, cRela
END IF

END SUBROUTINE DSMC_calc_var_P_vib


SUBROUTINE SumVibRelaxProb(iPair)
!===================================================================================================================================
! summes up the variable vibrational realaxation probabilities
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, VarVibRelaxProb, Coll_pData, SpecDSMC
USE MOD_Particle_Vars      ,ONLY: PartSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: VibProb
INTEGER                   :: cSpec1, cSpec2
!===================================================================================================================================

  ! variable vibrational relaxation probability has to average of all collisions
IF(Coll_pData(iPair)%cRela2.EQ.0) RETURN
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

END SUBROUTINE SumVibRelaxProb


SUBROUTINE FinalizeCalcVibRelaxProb(iElem)
!===================================================================================================================================
  ! Finalize the calculation of the variable vibrational relaxation probability in the cell for each iteration
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, VarVibRelaxProb 
USE MOD_Particle_Vars      ,ONLY: nSpecies

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iSpec
!===================================================================================================================================

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

END SUBROUTINE FinalizeCalcVibRelaxProb

END MODULE MOD_DSMC_Relaxation
