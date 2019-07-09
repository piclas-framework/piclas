!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
PUBLIC :: DSMC_VibRelaxDiatomic, SetMeanVibQua, CalcXiVibPart, CalcXiTotalEqui
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_VibRelaxDiatomic(iPair, iPart, FakXi)
!===================================================================================================================================
! Performs the vibrational relaxation of diatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, PartStateIntEn, Coll_pData, RadialWeighting
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartMPF, usevMPF, PEM, VarTimeStep
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
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
REAL                          :: MaxColQua, iRan, Ec, Phi, PartStateIntEnTemp, DeltaPartStateIntEn
INTEGER                       :: iQuaMax, iQua, iElem
!===================================================================================================================================
IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF

MaxColQua = Ec/(BoltzmannConst*SpecDSMC(PartSpecies(iPart))%CharaTVib)  &
          - DSMC%GammaQuant
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
IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
  IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
    iElem = PEM%Element(iPart)
    Phi = PartMPF(Coll_pData(iPair)%iPart_p2) / PartMPF(Coll_pData(iPair)%iPart_p1)
    PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                  * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
    PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + Phi * PartStateIntEnTemp
    ! search for new vib quant
    iQua = INT(PartStateIntEnTemp/ &
            (BoltzmannConst*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib) - DSMC%GammaQuant)
    CALL RANDOM_NUMBER(iRan)
    IF(iRan .LT. PartStateIntEnTemp/(BoltzmannConst &
                * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib) - DSMC%GammaQuant &
                - DBLE(iQua)) THEN
      iQua = iQua + 1
    END IF
    PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                  * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib
    DeltaPartStateIntEn = PartMPF(Coll_pData(iPair)%iPart_p1) &
                        * (PartStateIntEnTemp - PartStateIntEn(Coll_pData(iPair)%iPart_p1,1))
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
  END IF
ELSE
  PartStateIntEn(iPart,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
              * SpecDSMC(PartSpecies(iPart))%CharaTVib
END IF

END SUBROUTINE DSMC_VibRelaxDiatomic


SUBROUTINE SetMeanVibQua()
!===================================================================================================================================
! Computes the vibrational quant of species or the mean vibrational temp/energy/dofs (polyatomic)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,             ONLY : DSMC, CollInf, SpecDSMC, ChemReac, BGGas, PolyatomMolDSMC
  USE MOD_Particle_Vars,         ONLY : nSpecies
  USE MOD_DSMC_Analyze,          ONLY : CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER         :: iSpec, iPolyatMole
  REAL            :: iRan, VibQuaTemp
!===================================================================================================================================

  DO iSpec =1, nSpecies
    ! describe evib as quantum number
    IF (iSpec.EQ.BGGas%BGGasSpecies) THEN
      ChemReac%MeanEVibQua_PerIter(iSpec) = BGGas%BGMeanEVibQua
    ELSE
      IF(CollInf%Coll_SpecPartNum(iSpec).GT.0) THEN
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
            IF(ChemReac%MeanEVib_PerIter(iSpec).GT.SpecDSMC(iSpec)%EZeroPoint) THEN
              PolyatomMolDSMC(iPolyatMole)%TVib = CalcTVibPoly(ChemReac%MeanEVib_PerIter(iSpec), iSpec)
              PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean = 2*(ChemReac%MeanEVib_PerIter(iSpec)-SpecDSMC(iSpec)%EZeroPoint) &
                                                            / (BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%TVib)
            ELSEIF(ABS(ChemReac%MeanEVib_PerIter(iSpec)-SpecDSMC(iSpec)%EZeroPoint)/SpecDSMC(iSpec)%EZeroPoint.LT.1E-9) THEN
              ! Check relative difference between vibrational energy and zero-point energy
              PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean = 0.0
              PolyatomMolDSMC(iPolyatMole)%TVib = 0.0
            ELSE
              IPWRITE(*,*) ChemReac%MeanEVib_PerIter(iSpec), CollInf%Coll_SpecPartNum(iSpec), SpecDSMC(iSpec)%EZeroPoint
              CALL abort(&
                __STAMP__&
                ,'ERROR in SetMeanVibQua, energy less than zero-point energy, Species: ',iSpec)
            END IF
          ELSE
            ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
            VibQuaTemp = ChemReac%MeanEVib_PerIter(iSpec) / (BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
            CALL RANDOM_NUMBER(iRan)
            IF((VibQuaTemp-INT(VibQuaTemp)).GT.iRan) THEN
              ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(VibQuaTemp) + 1, SpecDSMC(iSpec)%MaxVibQuant-1)
            ELSE
              ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(VibQuaTemp), SpecDSMC(iSpec)%MaxVibQuant-1)
            END IF
          END IF
        ELSE
          ChemReac%MeanEVibQua_PerIter(iSpec) = 0
        END IF
      END IF
    END IF
  END DO
END SUBROUTINE SetMeanVibQua


SUBROUTINE CalcXiVibPart(TVib, iSpec, XiVibPart)
!===================================================================================================================================
! Calculation of the vibrational degrees of freedom for each characteristic vibrational temperature, used for chemical reactions
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: TVib  !
  INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)               :: XiVibPart(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: iDOF, iPolyatMole
!===================================================================================================================================

  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    XiVibPart = 0.0
    ! The vibrational energy of the dissociating molecule and the char. vib. temps of the product are used to determine a
    ! first guess for the vibrational degree of freedom
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      XiVibPart(iDOF) = (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / TVib) &
                / (EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/ TVib) - 1.0)
    END DO
  ELSE
    XiVibPart(1) = (2.0*SpecDSMC(iSpec)%CharaTVib / TVib) &
              / (EXP(SpecDSMC(iSpec)%CharaTVib / TVib) - 1.0)
  END IF
  RETURN

END SUBROUTINE CalcXiVibPart


SUBROUTINE CalcXiTotalEqui(iReac, iPair, Xi_rel, Weight1, Weight2, WeightProd, XiVibPart, XiElecPart)
!===================================================================================================================================
! Calculation of the vibrational degrees of freedom for each characteristic vibrational temperature, used for chemical reactions
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PolyatomMolDSMC, ChemReac, Coll_pData, DSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)             :: iReac, iPair      ! Reaction Number, Grow a pair number
  REAL, INTENT(IN)                :: Xi_rel, Weight1, Weight2, WeightProd
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT), OPTIONAL     :: XiVibPart(:,:), XiElecPart(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: iDOF, iPolyatMole, nProd, iProd, iQua
  INTEGER                         :: ProductReac(1:3)
  REAL                            :: ETotal, EZeroPoint, EGuess, Xi_Total, LowerTemp, UpperTemp, MiddleTemp, Xi_TotalTemp
  REAL                            :: SumOne, SumTwo, Weight(1:3)
  REAL                            :: eps_prec=0.1
!===================================================================================================================================

  ProductReac(1:3) = ChemReac%DefinedReact(iReac,2,1:3)

  IF(ProductReac(3).EQ.0) THEN
    Xi_Total = Xi_rel + SpecDSMC(ProductReac(1))%Xi_Rot + SpecDSMC(ProductReac(2))%Xi_Rot
    nProd = 2
  ELSE
    Xi_Total = Xi_rel + SpecDSMC(ProductReac(1))%Xi_Rot + SpecDSMC(ProductReac(2))%Xi_Rot + SpecDSMC(ProductReac(3))%Xi_Rot
    nProd = 3
  END IF

  Weight(1) = Weight1; Weight(2) = Weight2; Weight(3) = WeightProd

  ! Weighted total collision energy
  ETotal = Coll_pData(iPair)%Ec

  EZeroPoint = 0.0
  DO iProd = 1, nProd
    EZeroPoint = EZeroPoint + SpecDSMC(ProductReac(iProd))%EZeroPoint * Weight(iProd)
  END DO

  LowerTemp = 1.0
  UpperTemp = 2.*(ETotal - EZeroPoint) * nProd / SUM(Weight) / (Xi_Total * BoltzmannConst)
  DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
    MiddleTemp = 0.5*( LowerTemp + UpperTemp)
    Xi_TotalTemp = Xi_Total
    DO iProd = 1, nProd
      IF((SpecDSMC(ProductReac(iProd))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(iProd))%InterID.EQ.20)) THEN
        IF(SpecDSMC(ProductReac(iProd))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(ProductReac(iProd))%SpecToPolyArray
          XiVibPart(iProd,:) = 0.0
          ! The vibrational energy of the dissociating molecule and the char. vib. temps of the product are used to determine a
          ! first guess for the vibrational degree of freedom
          DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
            XiVibPart(iProd,iDOF) = (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / MiddleTemp) &
                      / (EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/ MiddleTemp) - 1.0)
            Xi_TotalTemp = Xi_TotalTemp + XiVibPart(iProd,iDOF)
          END DO
        ELSE
          XiVibPart(iProd,1) = (2.0*SpecDSMC(ProductReac(iProd))%CharaTVib / MiddleTemp) &
                    / (EXP(SpecDSMC(ProductReac(iProd))%CharaTVib / MiddleTemp) - 1.0)
          Xi_TotalTemp = Xi_TotalTemp + XiVibPart(iProd,1)
        END IF
      ELSE
        IF(PRESENT(XiVibPart)) XiVibPart(iProd,:) = 0.0
      END IF
      IF(DSMC%ElectronicModel) THEN
        IF((SpecDSMC(ProductReac(iProd))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(iProd))%FullyIonized)) THEN
          SumOne = 0.0
          SumTwo = 0.0
          DO iQua = 0, SpecDSMC(ProductReac(iProd))%MaxElecQuant-1
            SumOne = SumOne + SpecDSMC(ProductReac(iProd))%ElectronicState(1,iQua)*BoltzmannConst &
                            * SpecDSMC(ProductReac(iProd))%ElectronicState(2,iQua)  &
                            * EXP(-SpecDSMC(ProductReac(iProd))%ElectronicState(2,iQua)/MiddleTemp)
            SumTwo = SumTwo + SpecDSMC(ProductReac(iProd))%ElectronicState(1,iQua) &
                            * EXP(-SpecDSMC(ProductReac(iProd))%ElectronicState(2,iQua)/MiddleTemp)
          END DO
          IF((SumOne.GT.0.0).AND.(SumTwo*BoltzmannConst.GT.0.0)) THEN
            XiElecPart(iProd) = 2. * SumOne / (SumTwo * BoltzmannConst * MiddleTemp)
            Xi_TotalTemp = Xi_TotalTemp + XiElecPart(iProd)
          ELSE
            XiElecPart(iProd) = 0.0
          END IF
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
  RETURN

END SUBROUTINE CalcXiTotalEqui

END MODULE MOD_DSMC_Relaxation
