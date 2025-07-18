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

!===================================================================================================================================
!> Module with routines for relaxation
!===================================================================================================================================
MODULE MOD_DSMC_Relaxation
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

ABSTRACT INTERFACE
  SUBROUTINE RotRelaxDiaRoutine(iPair,iPart,FakXi)
    INTEGER,INTENT(IN)          :: iPair, iPart               ! index of collision pair
    REAL,INTENT(IN)             :: FakXi
  END SUBROUTINE
END INTERFACE

PROCEDURE(RotRelaxDiaRoutine),POINTER :: RotRelaxDiaRoutineFuncPTR !< pointer defining the function called for rotational relaxation
                                                                   !  depending on the RotRelaxModel (continuous or quantized)
                                                                   !  for diatomic molecules

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_VibRelaxDiatomic, CalcMeanVibQuaDiatomic, DSMC_calc_P_rot, DSMC_calc_var_P_vib
PUBLIC :: InitCalcVibRelaxProb, DSMC_calc_P_vib, SumVibRelaxProb, FinalizeCalcVibRelaxProb, DSMC_calc_P_elec
PUBLIC :: RotRelaxDiaRoutineFuncPTR, DSMC_RotRelaxDiaContinuous, DSMC_RotRelaxDiaQuant

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Rotational relaxation of diatomic molecules using quantized energy levels after Boyd
!> Physics of Fluids A: Fluid Dynamics (1989-1993) 5, 2278 (1993); doi: 10.1063/1.858531
!===================================================================================================================================
SUBROUTINE DSMC_RotRelaxDiaQuant(iPair,iPart,FakXi)
! MODULES
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC
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
INTEGER                       :: iSpec
INTEGER                       :: iQuant, J2, jIter
LOGICAL                       :: ARM
REAL                          :: fNorm, iRan, MaxValue, CurrentValue, Ec
!===================================================================================================================================
IF (usevMPF.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iSpec = PartSpecies(iPart)
! calculate maximum allowed quantum number if all of the collision energy would be transfered to rotational energy
J2 = INT(0.5 * (-1.+SQRT(1.+(4.*Ec)/(BoltzmannConst * SpecDSMC(iSpec)%CharaTRot))))
! Find max value of distribution numerically
MaxValue = 0.
DO jIter=0, J2
  CurrentValue = (2.*REAL(jIter) + 1.)*(Ec - REAL(jIter)*(REAL(jIter) + 1.)*BoltzmannConst*SpecDSMC(iSpec)%CharaTRot)**FakXi
  IF (CurrentValue .GT. MaxValue) THEN
      MaxValue = CurrentValue
  END IF
END DO
ARM = .TRUE.
CALL RANDOM_NUMBER(iRan)
iQuant = INT((1+J2)*iRan)
DO WHILE (ARM)
  fNorm = (2.*REAL(iQuant) + 1.)*(Ec - REAL(iQuant)*(REAL(iQuant) + 1.)*BoltzmannConst*SpecDSMC(iSpec)%CharaTRot)**FakXi / MaxValue
  CALL RANDOM_NUMBER(iRan)
  IF(fNorm .LT. iRan) THEN
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+J2)*iRan)
  ELSE
    ARM = .FALSE.
  END IF
END DO
PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * SpecDSMC(iSpec)%CharaTRot
END SUBROUTINE DSMC_RotRelaxDiaQuant


!===================================================================================================================================
!> Rotational relaxation of diatomic molecules using continuous energy levels after Pfeiffer/Nizenkov
!> Physics of Fluids 28, 027103 (2016); doi: 10.1063/1.4940989
!> Only separate routine for function pointer with RotRelaxModel
!===================================================================================================================================
SUBROUTINE DSMC_RotRelaxDiaContinuous(iPair,iPart,FakXi)
! MODULES
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC
USE MOD_Particle_Vars         ,ONLY: UseVarTimeStep, usevMPF
USE MOD_Particle_Vars         ,ONLY: PartSpecies
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
REAL                          :: LocalFakXi, iRan, CollEnergy
INTEGER                       :: iSpec
!===================================================================================================================================
IF (usevMPF.OR.UseVarTimeStep) THEN
  CollEnergy = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  CollEnergy = Coll_pData(iPair)%Ec
END IF

iSpec = PartSpecies(iPart)

! fix for changed FakXi for polyatomic
LocalFakXi = FakXi + 0.5*SpecDSMC(iSpec)%Xi_Rot
CALL RANDOM_NUMBER(iRan)
PartStateIntEn(2, iPart) = CollEnergy * (1.0 - iRan**(1.0/LocalFakXi))

END SUBROUTINE DSMC_RotRelaxDiaContinuous


!===================================================================================================================================
!> Performs the vibrational relaxation of diatomic molecules
!===================================================================================================================================
SUBROUTINE DSMC_VibRelaxDiatomic(iPair, iPart, FakXi)
! MODULES
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, PartStateIntEn, Coll_pData, AHO
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, c
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
REAL                          :: MaxColQua, iRan, Ec, ProbAccept
INTEGER                       :: iQuaMax, iQua, iSpec
!===================================================================================================================================

iRan = 1.
ProbAccept = 0.

IF (usevMPF.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF

iSpec = PartSpecies(iPart)

IF(DSMC%VibAHO) THEN ! AHO
  ! maximum quantum number
  iQuaMax = 2
  DO WHILE (Ec.GE.AHO%VibEnergy(iSpec,iQuaMax))
    ! collision energy is larger than vib energy for this quantum number --> increase quantum number and try again
    iQuaMax = iQuaMax + 1
    ! exit if this quantum number is larger as the table length (dissociation level is reached)
    IF (iQuaMax.GT.AHO%NumVibLevels(iSpec)) EXIT
  END DO
  ! accept iQuaMax - 1 as quantum number
  iQuaMax = iQuaMax - 1
  ! Ec without zero-point energy for calc of ProbAccept
  Ec = Ec - AHO%VibEnergy(iSpec,1)
  ! compare to random number and repeat until accepted
  DO WHILE (iRan.GT.ProbAccept)
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax + 1)
    ProbAccept = (PlanckConst * c * AHO%omegaE(iSpec) * (iQua-1) * (1. - AHO%chiE(iSpec) * iQua)) / Ec
    ! Avoid negative values in following probability calculation
    IF(ProbAccept.GT.1.) THEN
      ProbAccept = 0.
    ELSE
      ProbAccept = (1. - ProbAccept)**FakXi
      CALL RANDOM_NUMBER(iRan)
    END IF
  END DO
  ! vibrational energy is table value of the accepted quantum number
  PartStateIntEn(1,iPart) = AHO%VibEnergy(iSpec,iQua)
ELSE ! SHO
  MaxColQua = Ec/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
  iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(iSpec)%MaxVibQuant)
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(iRan * iQuaMax)
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (iRan.GT.(1 - REAL(iQua)/REAL(MaxColQua))**FakXi)
    !laux diss page 31
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
  END DO
  PartStateIntEn(1,iPart) = (iQua + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
END IF

END SUBROUTINE DSMC_VibRelaxDiatomic


!===================================================================================================================================
!> Computes the mean vibrational quantum number of diatomic species in a cell each iteration;
!> ChemReac%MeanEVibQua_PerIter is required for the determination of the vibrational degree of freedom, only used for diatomic
!> molecules. For the polyatomic case, the actual vibrational degree of freedom of each molecule is utilized and not a mean value.
!> The values for polyatomic molecules can have a greater spread, thus a mean value can prohibit reactions of highly excited
!> molecules at lower average temperatures.
!===================================================================================================================================
SUBROUTINE CalcMeanVibQuaDiatomic()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_DSMC_Vars,              ONLY : DSMC, CollInf, SpecDSMC, ChemReac, BGGas
USE MOD_Particle_Vars,          ONLY : nSpecies, Species
USE MOD_Particle_Analyze_Tools, ONLY : CalcTVibAHO, CalcXiVib
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
  IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
    IF(.NOT.SpecDSMC(iSpec)%PolyatomicMol) THEN
      ! Skip the background gas species (value initialized in dsmc_chemical_init.f90)
      IF(BGGas%BackgroundSpecies(iSpec)) CYCLE
      ! Only treat species present in the cell
      IF(CollInf%Coll_SpecPartNum(iSpec).GT.0.) THEN
        ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
        IF(DSMC%VibAHO) THEN ! AHO
          IF(ChemReac%MeanEVib_PerIter(iSpec).GT.0) THEN
            CALL CalcTVibAHO(iSpec, ChemReac%MeanEVib_PerIter(iSpec), VibQuaTemp)
            CALL CalcXiVib(VibQuaTemp, iSpec, XiVibTotal=ChemReac%MeanXiVib_PerIter(iSpec))
          ELSE
            ChemReac%MeanXiVib_PerIter(iSpec) = 0.
          END IF
        ELSE ! SHO
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
        END IF
      ELSE
        ChemReac%MeanEVibQua_PerIter(iSpec) = 0
        ChemReac%MeanXiVib_PerIter(iSpec) = 0.
      END IF  ! CollInf%Coll_SpecPartNum(iSpec).GT.0
    END IF    ! .NOT.SpecDSMC(iSpec)%PolyatomicMol
  END IF      ! (Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)
END DO        ! iSpec = 1, nSpecies

END SUBROUTINE CalcMeanVibQuaDiatomic


!===================================================================================================================================
!> Initialize the calculation of the variable vibrational relaxation probability in the cell for each iteration
!===================================================================================================================================
SUBROUTINE InitCalcVibRelaxProb()
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, VarVibRelaxProb
USE MOD_Particle_Vars      ,ONLY: nSpecies

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iSpec
!===================================================================================================================================

IF(DSMC%VibRelaxProb.EQ.2.0) THEN ! Set sums for variable vibrational relaxation to zero
  DO iSpec=1,nSpecies
    VarVibRelaxProb%ProbVibAvNew(iSpec) = 0
    VarVibRelaxProb%nCollis(iSpec) = 0
  END DO
END IF

END SUBROUTINE InitCalcVibRelaxProb


!===================================================================================================================================
!> Calculation of probability for rotational relaxation. Different Models implemented:
!> 0 - Constant Probability
!> 1 - No rotational relaxation. RotRelaxProb = 0
!> 2 - Boyd for N2
!> 3 - Zhang (Nonequilibrium Direction Dependent)
!> 4 - Boyd for H2, J . Fluid Mech. (1994), uol. 280, p p . 41-67
!===================================================================================================================================
SUBROUTINE DSMC_calc_P_rot(iSpec1, iSpec2, iPair, iPart, Xi_rel, ProbRot, ProbRotMax)
! MODULES
USE MOD_Globals            ,ONLY: Abort
USE MOD_Globals_Vars       ,ONLY: Pi, BoltzmannConst
USE MOD_Particle_Vars      ,ONLY: UseVarTimeStep, usevMPF
USE MOD_part_tools         ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars          ,ONLY: SpecDSMC, Coll_pData, PartStateIntEn, DSMC, useRelaxProbCorrFactor, CollInf
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
REAL                      :: TransEn, RotEn, RotDOF
REAL                      :: CorrFact           !> CorrFact: To correct sample bias (fewer DSMC particles than natural ones)
REAL                      :: LumpkinCorr, VHSCorr
!===================================================================================================================================
! Note that during probability calculation, collision energy only contains translational part
IF (usevMPF.OR.UseVarTimeStep) THEN
  TransEn = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  TransEn = Coll_pData(iPair)%Ec
END IF

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
ELSE IF(DSMC%RotRelaxProb.EQ.2.0) THEN ! P_rot according to Boyd (based on Parker's model)

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
ELSEIF(DSMC%RotRelaxProb.EQ.4.0) THEN ! P_rot according to Boyd for H2, J . Fluid Mech. (1994), uol. 280, p p . 41-67

  LumpkinCorr = (5.-2.*(CollInf%omega(iSpec1,iSpec2)+0.5)+2.*RotDOF)/(5.-2.*(CollInf%omega(iSpec1,iSpec2)+0.5))

  VHSCorr = 8.*GAMMA(9./2.-(CollInf%omega(iSpec1,iSpec2)+0.5))/(15.*PI*GAMMA(5./2.-(CollInf%omega(iSpec1,iSpec2)+0.5))) &
          * (208.-12.*(CollInf%omega(iSpec1,iSpec2)+0.5))/(211.-12*(CollInf%omega(iSpec1,iSpec2)+0.5)*(5./2.-(CollInf%omega(iSpec1,iSpec2)+0.5)))

  ProbRot = LumpkinCorr*VHSCorr*1./10480.*GAMMA(RotDOF+5./2.-(CollInf%omega(iSpec1,iSpec2)+0.5))/GAMMA(RotDOF+5./2.) &
          * ((TransEn+RotEn)/BoltzmannConst)**(CollInf%omega(iSpec1,iSpec2)+0.5)
ELSE
  CALL Abort(__STAMP__,'Error! Model for rotational relaxation undefined:',RealInfoOpt=DSMC%RotRelaxProb)
END IF

END SUBROUTINE DSMC_calc_P_rot


!===================================================================================================================================
!> Calculation of probability for electronic relaxation.
!===================================================================================================================================
SUBROUTINE DSMC_calc_P_elec(iSpec1, iSpec2, ProbElec)
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


!===================================================================================================================================
!> Calculation of probability for vibrational relaxation. Different Models implemented:
!> 0 - Constant Probability
!> 1 - No vibrational relaxation. VibRelaxProb = 0
!> 2 - Boyd with correction of Abe
!===================================================================================================================================
SUBROUTINE DSMC_calc_P_vib(iPair, iPart, iSpec, jSpec, Xi_rel, iElem, ProbVib)
! MODULES
USE MOD_Globals            ,ONLY: Abort
USE MOD_Globals_Vars       ,ONLY: BoltzmannConst, PlanckConst, c, eV2Joule
USE MOD_DSMC_Vars          ,ONLY: SpecDSMC, DSMC, VarVibRelaxProb, useRelaxProbCorrFactor, PolyatomMolDSMC, CollInf, Coll_pData, AHO
USE MOD_DSMC_Vars          ,ONLY: PartStateIntEn
USE MOD_MCC_Vars           ,ONLY: XSec_Relaxation, SpecXSec
USE MOD_MCC_XSec           ,ONLY: XSec_CalcVibRelaxProb
USE MOD_Part_Tools         ,ONLY: GetParticleWeight
USE MOD_Particle_Vars      ,ONLY: UseVarTimeStep, usevMPF
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iPair, iPart, iSpec, jSpec, iElem
REAL, INTENT(IN)          :: Xi_rel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: CorrFact       ! CorrFact: To correct sample Bias (fewer DSMC particles than natural ones)
REAL                      :: Ec, DissTemp, MaxColQua, Tcoll, Tref, Zref, Zvib
INTEGER                   :: iPolyatMole, iDOF, iCase, iQuantMax
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

ELSE IF(DSMC%VibRelaxProb.EQ.3.0) THEN
  ! Calculation of vibrational relaxation probability according to Bird, can be used with SHO and AHO
  ! collision energy = relative translational energy of iPair + pre-collision vibrational energy of iPart
  IF (usevMPF.OR.UseVarTimeStep) THEN
    Ec = (Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart)*GetParticleWeight(iPart)) / GetParticleWeight(iPart)
  ELSE
    Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart)*GetParticleWeight(iPart)
  END IF
  ! dissociation temperature
  DissTemp = SpecDSMC(iSpec)%Ediss_eV * eV2Joule / BoltzmannConst

  IF(DSMC%VibAHO) THEN
    ! maximum quantum number
    iQuantMax = 2
    DO WHILE (Ec.GE.AHO%VibEnergy(iSpec,iQuantMax))
      ! collision energy is larger than vib energy for this quantum number --> increase quantum number and try again
      iQuantMax = iQuantMax + 1
      ! exit if this quantum number is larger as the table length (dissociation level is reached)
      IF (iQuantMax.GT.AHO%NumVibLevels(iSpec)) EXIT
    END DO
    ! accept iQuantMax - 1 as quantum number
    iQuantMax = iQuantMax - 1
    ! collision temperature
    Tcoll = (PlanckConst * c * AHO%omegaE(iSpec) * (iQuantMax-1) * (1. - AHO%chiE(iSpec) * iQuantMax)) &
      / (BoltzmannConst * (3.-SpecDSMC(iSpec)%omega))

  ELSE
    ! maximum quantum number
    MaxColQua = Ec/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
    iQuantMax = MIN(INT(MaxColQua) + 1, SpecDSMC(iSpec)%DissQuant)
    ! collision temperature
    Tcoll = iQuantMax * SpecDSMC(iSpec)%CharaTVib / (3.-SpecDSMC(iSpec)%omega)
  END IF

  ! Programmer-defined value for the reference temperature. Rule of thumb: characteristic vibrational temperature (SHO).
  ! Can be adapted in future to
  ! Tref = SpecDSMC(iSpec)%CharaTVib
  ! but requires read-in of characteristic vibrational temperature also for AHO model
  Tref = 2500.
  ! vibrational collision number at reference temperature
  Zref = SpecDSMC(iSpec)%C1(jSpec) * EXP(SpecDSMC(iSpec)%C2(jSpec) * Tref**(-1./3.)) / (Tref**(SpecDSMC(iSpec)%omega+0.5))
  ! vibrational collision number
  Zvib = (DissTemp/Tcoll)**(SpecDSMC(iSpec)%omega+0.5) * (Zref * (DissTemp/Tref)**(-SpecDSMC(iSpec)%omega-0.5)) &
    **(((DissTemp/Tcoll)**(1./3.) - 1.) / ((DissTemp/Tref)**(1./3.) - 1.))
  IF (Zvib.LT.1.) Zvib = 1.
  ! vibrational relaxation probability
  ProbVib = 1. / Zvib

ELSE
  CALL Abort(__STAMP__,'Error! Model for vibrational relaxation undefined:',RealInfoOpt=DSMC%VibRelaxProb)
END IF

IF(DSMC%CalcQualityFactors) THEN
  DSMC%CalcVibProb(iSpec,1) = DSMC%CalcVibProb(iSpec,1) + ProbVib
  DSMC%CalcVibProb(iSpec,3) = DSMC%CalcVibProb(iSpec,3) + 1
END IF

END SUBROUTINE DSMC_calc_P_vib


!===================================================================================================================================
!> Calculation of probability for vibrational relaxation for variable relaxation rates. This has to average over all collisions!
!> No instantaneous variable probability calculateable
!===================================================================================================================================
SUBROUTINE DSMC_calc_var_P_vib(iSpec, jSpec, iPair, ProbVib)
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


!===================================================================================================================================
!> Sums up the variable vibrational relaxation probabilities
!===================================================================================================================================
SUBROUTINE SumVibRelaxProb(iPair)
! MODULES
USE MOD_DSMC_Vars          ,ONLY: DSMC, VarVibRelaxProb, Coll_pData
USE MOD_Particle_Vars      ,ONLY: PartSpecies, Species
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
  IF((Species(cSpec1)%InterID.EQ.2).OR.(Species(cSpec1)%InterID.EQ.20)) THEN
    CALL DSMC_calc_var_P_vib(cSpec1,cSpec2,iPair,VibProb)
    VarVibRelaxProb%ProbVibAvNew(cSpec1) = VarVibRelaxProb%ProbVibAvNew(cSpec1) + VibProb
    VarVibRelaxProb%nCollis(cSpec1) = VarVibRelaxProb%nCollis(cSpec1) + 1
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CalcVibProb(cSpec1,2) = MAX(DSMC%CalcVibProb(cSpec1,2),VibProb)
    END IF
  END IF
  IF((Species(cSpec2)%InterID.EQ.2).OR.(Species(cSpec2)%InterID.EQ.20)) THEN
    CALL DSMC_calc_var_P_vib(cSpec2,cSpec1,iPair,VibProb)
    VarVibRelaxProb%ProbVibAvNew(cSpec2) = VarVibRelaxProb%ProbVibAvNew(cSpec2) + VibProb
    VarVibRelaxProb%nCollis(cSpec2) = VarVibRelaxProb%nCollis(cSpec2) + 1
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%CalcVibProb(cSpec2,2) = MAX(DSMC%CalcVibProb(cSpec2,2),VibProb)
    END IF
  END IF
END IF

END SUBROUTINE SumVibRelaxProb


!===================================================================================================================================
!> Finalize the calculation of the variable vibrational relaxation probability in the cell for each iteration
!===================================================================================================================================
SUBROUTINE FinalizeCalcVibRelaxProb(iElem)
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
