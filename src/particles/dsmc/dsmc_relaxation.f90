#include "boltzplatz.h"

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
PUBLIC :: DSMC_VibRelaxDiatomic, SetMeanVibQua, CalcXiVibPart
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_VibRelaxDiatomic(iPair, iPart, FakXi)
!===================================================================================================================================
! Performs the vibrational relaxation of diatomic molecules
!===================================================================================================================================
! MODULES  
  USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC, PartStateIntEn, Coll_pData
  USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst, PartMPF, usevMPF, GEO, PEM
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iPair
  REAL(KIND=8), INTENT(IN)      :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: MaxColQua, iRan, Ec, Phi, PartStateIntEnTemp, DeltaPartStateIntEn
  INTEGER                       :: iQuaMax, iQua, iElem
!===================================================================================================================================
  Ec = Coll_pData(iPair)%Ec
  MaxColQua = Ec/(BoltzmannConst*SpecDSMC(PartSpecies(iPart))%CharaTVib)  &
            - DSMC%GammaQuant
  iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(iPart))%MaxVibQuant)
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(iRan * iQuaMax)
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
   !laux diss page 31
   CALL RANDOM_NUMBER(iRan)
   iQua = INT(iRan * iQuaMax)
   CALL RANDOM_NUMBER(iRan)
  END DO
  IF (usevMPF) THEN    
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
  USE MOD_DSMC_Vars,             ONLY : DSMC, CollInf, SpecDSMC, ChemReac, BGGas, PolyatomMolDSMC
  USE MOD_Particle_Vars,         ONLY : BoltzmannConst, nSpecies
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
            ELSEIF(ALMOSTEQUAL(ChemReac%MeanEVib_PerIter(iSpec),SpecDSMC(iSpec)%EZeroPoint)) THEN
              PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean = 0.0
              PolyatomMolDSMC(iPolyatMole)%TVib = 0.0
            ELSE
              IPWRITE(*,*) ChemReac%MeanEVib_PerIter(iSpec), CollInf%Coll_SpecPartNum(iSpec), SpecDSMC(iSpec)%EZeroPoint
              CALL abort(&
                __STAMP__&
                ,'ERROR in Calc_XiVib_Poly, energy less than zero-point energy, Species: ',iSpec)
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
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
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

END MODULE MOD_DSMC_Relaxation
