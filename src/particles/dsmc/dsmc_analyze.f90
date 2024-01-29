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

MODULE MOD_DSMC_Analyze
!===================================================================================================================================
! Module for DSMC Sampling and Output
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
PUBLIC :: DSMC_data_sampling, CalcMeanFreePath,WriteDSMCToHDF5
PUBLIC :: CalcTVib, CalcGammaVib
PUBLIC :: CalcInstantTransTemp, SummarizeQualityFactors, DSMCMacroSampling
PUBLIC :: SamplingRotVibRelaxProb, CalcInstantElecTempXi
!===================================================================================================================================

CONTAINS

REAL FUNCTION CalcTVib(ChaTVib,MeanEVib,nMax)
!===================================================================================================================================
!> Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(IN)                :: ChaTVib,MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
INTEGER, INTENT(IN)             :: nMax              ! INT(CharaTDisss/CharaTVib) + 1
REAL(KIND=8)                    :: LowerVal, UpperVal, MiddleVal, MaxPosiVal  ! upper and lower value of zero point search
REAl(KIND=8)                    :: eps_prec=0.1   ! precision of zero point search
REAL(KIND=8)                    :: ZeroVal1, ZeroVal2 ! both fuction values to compare
!===================================================================================================================================

IF (MeanEVib.GT.0) THEN
  !.... Initial limits for a: lower limit = very small value
  !                           upper limit = max. value allowed by system
  !     zero point = CharaTVib / TVib
  LowerVal  = 1.0/(2.0*nMax)                                    ! Tvib is max for nMax => lower limit = 1.0/nMax
  UpperVal  = LOG(HUGE(MiddleVal*nMax))/nMax-1.0/(2.0 * nMax)   ! upper limit = for max possible EXP(nMax*MiddleVal)-value
  MaxPosiVal = LOG(HUGE(MaxPosiVal))  ! maximum value possible in system
  DO WHILE (ABS(LowerVal-UpperVal).GT.eps_prec)                      !  Let's search the zero point by bisection
    MiddleVal = 0.5*(LowerVal+UpperVal)

    IF ((LowerVal.GT.MaxPosiVal).OR.(MiddleVal.GT.MaxPosiVal)) THEN
       CALL Abort(&
__STAMP__&
,'Cannot find zero point in TVib Calculation Function! CharTVib:',RealInfoOpt=ChaTVib)
    END IF

    ! Calc of actual function values
    ZeroVal1 = DSMC%GammaQuant + 1/(EXP(LowerVal)-1) - nMax/(EXP(nMax*LowerVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
    ZeroVal2 = DSMC%GammaQuant + 1/(EXP(MiddleVal)-1) - nMax/(EXP(nMax*MiddleVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
    ! decision of direction of bisection
    IF (ZeroVal1*ZeroVal2.LT.0) THEN
      UpperVal = MiddleVal
    ELSE
      LowerVal = MiddleVal
    END IF
  END DO
  CalcTVib = ChaTVib/LowerVal ! LowerVal = CharaTVib / TVib
ELSE
  CalcTVib = 0
END IF

RETURN

END FUNCTION CalcTVib


REAL FUNCTION CalcMeanFreePath(SpecPartNum, nPart, Volume, opt_temp)
!===================================================================================================================================
!> Calculation of the mean free path for the hard sphere and variable hard sphere (if a temperature is given)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: Pi
USE MOD_Particle_Vars ,ONLY: Species, nSpecies, usevMPF
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting, CollInf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Volume,SpecPartNum(:),nPart
REAL, OPTIONAL, INTENT(IN)      :: opt_temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iSpec, jSpec
REAL                            :: Temp, MFP_Tmp, MacroParticleFactor
!===================================================================================================================================
CalcMeanFreePath = 0.0

IF (nPart.LE.1 .OR. ALL(SpecPartNum.EQ.0.) .OR.Volume.EQ.0) RETURN
! Calculation of mixture reference diameter
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  MacroParticleFactor = 1.
ELSE
  MacroParticleFactor = Species(1)%MacroParticleFactor ! assumption: weighting factor of all species are identical!!!
END IF

! Calculation of mean free path for a gas mixture (Bird 1986, p. 96, Eq. 4.77)
! (only defined for a single weighting factor, if a temperature is present calculation of the mean free path with the VHS model)
IF(PRESENT(opt_temp)) THEN
  Temp = opt_temp
  IF(Temp.LE.0.0) RETURN
    DO iSpec = 1, nSpecies
      MFP_Tmp = 0.0
      IF(SpecPartNum(iSpec).GT.0.0) THEN ! skipping species not present in the cell
        DO jSpec = 1, nSpecies
          IF(SpecPartNum(jSpec).GT.0.0) THEN ! skipping species not present in the cell
            MFP_Tmp = MFP_Tmp + (Pi*CollInf%dref(iSpec,jSpec)**2.*SpecPartNum(jSpec)*MacroParticleFactor / Volume &
                                  * (CollInf%Tref(iSpec,jSpec)/Temp)**(CollInf%omega(iSpec,jSpec)) &
                                  * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
          END IF
        END DO
        CalcMeanFreePath = CalcMeanFreePath + (SpecPartNum(iSpec) / nPart) / MFP_Tmp
      END IF
    END DO
ELSE
  DO iSpec = 1, nSpecies
    MFP_Tmp = 0.0
    IF(SpecPartNum(iSpec).GT.0.0) THEN ! skipping species not present in the cell
      DO jSpec = 1, nSpecies
        IF(SpecPartNum(jSpec).GT.0.0) THEN ! skipping species not present in the cell
          MFP_Tmp = MFP_Tmp + (Pi*CollInf%dref(iSpec,jSpec)**2.*SpecPartNum(jSpec)*MacroParticleFactor / Volume &
                                * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
        END IF
      END DO
      CalcMeanFreePath = CalcMeanFreePath + (SpecPartNum(iSpec) / nPart) / MFP_Tmp
    END IF
  END DO
END IF

RETURN

END FUNCTION CalcMeanFreePath


SUBROUTINE CalcGammaVib()
!===================================================================================================================================
!> calculate Gamma_vib factor necessary for correction of vibrational relaxation according to Gimelshein et al.
!> -> 'Vibrational Relaxation Rates in the DSMC Method', Physics of Fluids V14 No12, 2002
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iDOF, iPolyatMole
REAL                  :: CharaTVib, TempTrans, GammaVib
!===================================================================================================================================

! Calculate GammaVib Factor  = Xi_Vib² * exp(CharaTVib/T_trans) / 2
DO iSpec = 1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    ! First, reset the GammaVib array/value
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (DSMC%PolySingleMode) THEN
        PolyatomMolDSMC(iPolyatMole)%GammaVib = 0.
      ELSE
        SpecDSMC(iSpec)%GammaVib = 0.
      END IF
    ELSE
      SpecDSMC(iSpec)%GammaVib = 0.
    END IF
    TempTrans = DSMC%InstantTransTemp(iSpec)
    IF(TempTrans.GT.0.0) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        CharaTVib = MAXVAL(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
      ELSE
        CharaTVib = SpecDSMC(iSpec)%CharaTVib
      END IF
      IF(CharaTVib/TempTrans.LT.80) THEN
        ! If CharaTVib/TempTrans is too high the exp function can produce NAN
        ! CharaTVib/TempTrans=80 results in a of GammaVib=2.31020977644213E-31
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            CharaTVib = PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            GammaVib = (2.*CharaTVib / (TempTrans *(EXP(CharaTVib/TempTrans)-1.)))**2. * EXP(CharaTVib/TempTrans) / 2.
            IF (DSMC%PolySingleMode) THEN
              PolyatomMolDSMC(iPolyatMole)%GammaVib(iDOF) = GammaVib
            ELSE
              SpecDSMC(iSpec)%GammaVib = SpecDSMC(iSpec)%GammaVib + GammaVib
            END IF
          END DO
        ELSE
          CharaTVib = SpecDSMC(iSpec)%CharaTVib
          SpecDSMC(iSpec)%GammaVib = (2.*CharaTVib / (TempTrans *(EXP(CharaTVib/TempTrans)-1.)))**2. * EXP(CharaTVib/TempTrans) / 2.
        END IF
      END IF  ! CharaTVib/TempTrans.LT.80
    END IF    ! TempTrans.GT.0.0
  END IF      ! (SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)
END DO        ! iSpec = 1, nSpecies

END SUBROUTINE CalcGammaVib


SUBROUTINE CalcInstantElecTempXi(iPartIndx,PartNum)
!===================================================================================================================================
!> Calculation of the instantaneous translational temperature for the cell
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars              ,ONLY: DSMC, CollInf, PartStateIntEn
USE MOD_Particle_Vars          ,ONLY: PartSpecies, nSpecies
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Analyze_Tools ,ONLY: CalcTelec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: PartNum
INTEGER, INTENT(IN)   :: iPartIndx(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iPart, SpecPartNum_Simu(nSpecies), PartID, SpecID
REAL                  :: ElecEnergy(nSpecies), partWeight
!===================================================================================================================================
! Actual number of particles, required to avoid calculation of temperature from one particle of the species
SpecPartNum_Simu = 0
! Sum of particle number, might be weighted/multiplied with PartMPF and/or VariableTimeStep
! Setting temperature to zero
DSMC%InstantTXiElec = 0.
ElecEnergy=0.0

DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  ElecEnergy(SpecID) = ElecEnergy(SpecID) + PartStateIntEn(3,PartID) * partWeight
  SpecPartNum_Simu(SpecID) = SpecPartNum_Simu(SpecID) + 1
END DO

DO iSpec=1, nSpecies
  IF(SpecPartNum_Simu(iSpec).GT.1) THEN
    ElecEnergy(iSpec) = ElecEnergy(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
    ! Compute temperatures
    DSMC%InstantTXiElec(1,iSpec) = CalcTelec(ElecEnergy(iSpec), iSpec)
    IF (DSMC%InstantTXiElec(1,iSpec).GT.0.0) THEN
      DSMC%InstantTXiElec(2,iSpec) = 2.*ElecEnergy(iSpec) /(BoltzmannConst*DSMC%InstantTXiElec(1,iSpec))
    ELSE
      DSMC%InstantTXiElec(2,iSpec) = 0.0
    END IF
  ELSE
    DSMC%InstantTXiElec(1:2,iSpec) = 0.0
  END IF
END DO
END SUBROUTINE CalcInstantElecTempXi


SUBROUTINE CalcInstantTransTemp(iPartIndx,PartNum,vBulk)
!===================================================================================================================================
!> Calculation of the instantaneous translational temperature for the cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_DSMC_Vars     ,ONLY: DSMC
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species, nSpecies
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: PartNum
INTEGER, INTENT(IN)   :: iPartIndx(:)
REAL, OPTIONAL, INTENT(INOUT) :: vBulk(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iPart, SpecPartNum_Simu(nSpecies), PartID, SpecID
REAL                  :: PartV(nSpecies,3), PartV2(nSpecies), SumSpecPartNum, partWeight, vmag2, V_rel(3), tempweighttotal
REAL                  :: MeanPartV(nSpecies,3), totalweight2(nSpecies),totalweight(nSpecies), Ener(nSpecies), EnerTotal, tempmass
REAL                  :: vBulkTemp(3), TotalMass
!===================================================================================================================================

PartV = 0.
PartV2 = 0.
! Actual number of particles, required to avoid calculation of temperature from one particle of the species
SpecPartNum_Simu = 0
! Sum of particle number, might be weighted/multiplied with PartMPF and/or VariableTimeStep
SumSpecPartNum = 0.
! Setting temperature to zero
DSMC%InstantTransTemp = 0.
totalweight2 = 0.; totalweight = 0.
IF (PRESENT(vBulk)) vBulk(1:3) = 0.; TotalMass =0.0
DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  PartV(SpecID,1:3) = PartV(SpecID,1:3) + PartState(4:6,PartID) * partWeight
  IF (PRESENT(vBulk))  THEN
    vBulk(1:3) = vBulk(1:3) + PartState(4:6,PartID)*Species(SpecID)%MassIC*partWeight
    TotalMass = TotalMass + Species(SpecID)%MassIC*partWeight
  END IF
  totalweight(SpecID) = totalweight(SpecID) + partWeight
  totalweight2(SpecID) = totalweight2(SpecID) + partWeight*partWeight
  SpecPartNum_Simu(SpecID) = SpecPartNum_Simu(SpecID) + 1
END DO
IF (PRESENT(vBulk)) vBulk(1:3) = vBulk(1:3) / TotalMass

DO iSpec=1, nSpecies
  IF(SpecPartNum_Simu(iSpec).GT.1) THEN
    ! Compute velocity averages
    MeanPartV(iSpec,1:3)  = PartV(iSpec,1:3) / totalWeight(iSpec)
  ELSE
    MeanPartV(iSpec,1:3) = 0.
  END IF
END DO

DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  V_rel(1:3) = PartState(4:6,PartID)-MeanPartV(SpecID,1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  PartV2(SpecID) = PartV2(SpecID) + vmag2 * partWeight
END DO

Ener=0.0; EnerTotal=0.0
tempweighttotal=0.0; tempmass = 0.0; vBulkTemp=0.0
DO iSpec = 1, nSpecies
  IF ((SpecPartNum_Simu(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(PartV2(iSpec)))) THEN
    DSMC%InstantTransTemp(iSpec) = Species(iSpec)%MassIC * PartV2(iSpec) &
        /(3.0*BoltzmannConst*(totalweight(iSpec) - totalweight2(iSpec)/totalweight(iSpec)))
    Ener(iSpec) =  3./2.*BoltzmannConst*DSMC%InstantTransTemp(iSpec) * totalweight(iSpec)
    vmag2 = MeanPartV(iSpec,1)**(2.) + MeanPartV(iSpec,2)**(2.) + MeanPartV(iSpec,3)**(2.)
    Ener(iSpec) = Ener(iSpec) + totalweight(iSpec) * Species(iSpec)%MassIC / 2. * vmag2
    EnerTotal = EnerTotal + Ener(iSpec)
    tempweighttotal = tempweighttotal + totalweight(iSpec)
    tempmass = tempmass +  totalweight(iSpec) * Species(iSpec)%MassIC
    vBulkTemp(1:3) = vBulkTemp(1:3) + MeanPartV(iSpec,1:3)*totalweight(iSpec) * Species(iSpec)%MassIC
  END IF
END DO
IF ((tempmass.GT.0.0).AND.(tempweighttotal.GT.0.0)) THEN
  vBulkTemp(1:3) = vBulkTemp(1:3) / tempmass
  vmag2 = vBulkTemp(1)*vBulkTemp(1) + vBulkTemp(2)*vBulkTemp(2) + vBulkTemp(3)*vBulkTemp(3)
  EnerTotal = EnerTotal -  tempmass / 2. * vmag2
  DSMC%InstantTransTemp(nSpecies+1) = 2. * EnerTotal / (3.*tempweighttotal*BoltzmannConst)
END IF

END SUBROUTINE CalcInstantTransTemp

SUBROUTINE DSMC_data_sampling()
!===================================================================================================================================
!> Sampling of variables velocity and energy for DSMC
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars              ,ONLY: useDSMC, PartStateIntEn, DSMC, CollisMode, SpecDSMC, DSMC_Solution, AmbipolElecVelo
USE MOD_Part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, PEM, Species, DoVirtualCellMerge, VirtMergedCells
USE MOD_Mesh_Vars              ,ONLY: offSetElem
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iElem, iSpec
REAL                          :: partWeight
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
DSMC%SampNum = DSMC%SampNum + 1
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    iSpec = PartSpecies(iPart)
    iElem = PEM%LocalElemID(iPart)
    IF (DoVirtualCellMerge) THEN
      IF (VirtMergedCells(iElem)%isMerged) iElem = VirtMergedCells(iElem)%MasterCell - offSetElem
    END IF
    partWeight = GetParticleWeight(iPart)
    DSMC_Solution(1:3,iElem,iSpec) = DSMC_Solution(1:3,iElem,iSpec) + PartState(4:6,iPart)*partWeight
    DSMC_Solution(4:6,iElem,iSpec) = DSMC_Solution(4:6,iElem,iSpec) + PartState(4:6,iPart)**2*partWeight
    DSMC_Solution(7,iElem,iSpec) = DSMC_Solution(7,iElem, iSpec) + partWeight  !density number
    IF(useDSMC)THEN
      IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
        IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.20)) THEN
          DSMC_Solution(8,iElem, iSpec) = DSMC_Solution(8,iElem, iSpec) &
            + (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
          DSMC_Solution(9,iElem, iSpec) = DSMC_Solution(9,iElem, iSpec)+PartStateIntEn(2,iPart)*partWeight
        END IF
        IF (DSMC%ElectronicModel.GT.0) THEN
          IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
            DSMC_Solution(10,iElem,iSpec)=DSMC_Solution(10,iElem,iSpec)+PartStateIntEn(3,iPart)*partWeight
          END IF
        END IF
      END IF
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
          DSMC_Solution(1:3,iElem,DSMC%AmbiDiffElecSpec) = DSMC_Solution(1:3,iElem,DSMC%AmbiDiffElecSpec) &
            + AmbipolElecVelo(iPart)%ElecVelo(1:3)*partWeight
          DSMC_Solution(4:6,iElem,DSMC%AmbiDiffElecSpec) = DSMC_Solution(4:6,iElem,DSMC%AmbiDiffElecSpec) &
            + AmbipolElecVelo(iPart)%ElecVelo(1:3)**2*partWeight
          DSMC_Solution(7,iElem,DSMC%AmbiDiffElecSpec) = DSMC_Solution(7,iElem, DSMC%AmbiDiffElecSpec) + partWeight
          DSMC_Solution(11,iElem, DSMC%AmbiDiffElecSpec) = DSMC_Solution(11,iElem, DSMC%AmbiDiffElecSpec) + 1.0
        END IF
      END IF
    END IF
    DSMC_Solution(11,iElem, iSpec) = DSMC_Solution(11,iElem, iSpec) + 1.0 !simpartnum
  END IF
END DO
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE DSMC_data_sampling


SUBROUTINE DSMC_output_calc(nVar,nVar_quality,nVarloc,DSMC_MacroVal)
!===================================================================================================================================
!> Calculation of the macroscopic output from sampled particles properties including the mixture values (Total_).
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_BGK_Vars               ,ONLY: BGKInitDone, BGK_QualityFacSamp
USE MOD_DSMC_Vars              ,ONLY: DSMC_Solution, CollisMode, SpecDSMC, DSMC, useDSMC, RadialWeighting, BGGas
USE MOD_FPFlow_Vars            ,ONLY: FPInitDone, FP_QualityFacSamp
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, WriteMacroVolumeValues, usevMPF, Symmetry
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, VarTimeStep
USE MOD_Particle_Vars          ,ONLY: DoVirtualCellMerge, VirtMergedCells
USE MOD_Particle_TimeStep      ,ONLY: GetParticleTimeStep
USE MOD_Restart_Vars           ,ONLY: RestartTime
USE MOD_TimeDisc_Vars          ,ONLY: time,TEnd,iter,dt
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemMidPoint_Shared, ElemVolume_Shared
USE MOD_Mesh_Vars              ,ONLY: offSetElem
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Tools ,ONLY: CalcTelec,CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: nVar,nVar_quality,nVarloc
REAL,INTENT(INOUT)      :: DSMC_MacroVal(1:nVar+nVar_quality,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem, iSpec, nVarCount, nSpecTemp, nVarCountRelax, bgSpec
REAL                    :: TVib_TempFac, iter_loc
REAL                    :: MolecPartNum, HeavyPartNum
!===================================================================================================================================
! nullify
DSMC_MacroVal = 0.0

nVarCount=0
DO iElem = 1, nElems ! element/cell main loop
  IF (DoVirtualCellMerge) THEN
    IF (VirtMergedCells(iElem)%isMerged) CYCLE
  END IF
  MolecPartNum = 0.0
  HeavyPartNum = 0.0
  ! Avoid the output and calculation of total values for a single species, associate construct for Total_ points to the same
  ! array as the Macro_ link for a single species and total values are only calculated for nSpecies > 1. Workaround required
  ! since associate construct is around the DO iSpec=1,nSpecies loop.
  IF(nSpecies.EQ.1) THEN
    nSpecTemp = 0
  ELSE
    nSpecTemp = nSpecies
  END IF
  ASSOCIATE ( Total_Velo     => DSMC_MacroVal(nVarLoc*nSpecTemp+1:nVarLoc*nSpecTemp+3,iElem) ,&
              Total_Temp     => DSMC_MacroVal(nVarLoc*nSpecTemp+4:nVarLoc*nSpecTemp+6,iElem) ,&
              Total_TempMean => DSMC_MacroVal(nVarLoc*nSpecTemp+12,iElem)            ,&
              Total_Density  => DSMC_MacroVal(nVarLoc*nSpecTemp+7,iElem)             ,&
              Total_TempVib  => DSMC_MacroVal(nVarLoc*nSpecTemp+8,iElem)             ,&
              Total_TempRot  => DSMC_MacroVal(nVarLoc*nSpecTemp+9,iElem)             ,&
              Total_Tempelec => DSMC_MacroVal(nVarLoc*nSpecTemp+10,iElem)            ,&
              Total_PartNum  => DSMC_MacroVal(nVarLoc*nSpecTemp+11,iElem)            ,&
              SimVolume      => ElemVolume_Shared(GetCNElemID(iElem+offSetElem)) &
              )
    ! compute simulation cell volume
    DO iSpec = 1, nSpecies
      ASSOCIATE ( PartVelo   => DSMC_Solution(1:3,iElem,iSpec) ,&
                  PartVelo2  => DSMC_Solution(4:6,iElem,iSpec) ,&
                  PartNum    => DSMC_Solution(7,iElem,iSpec)   ,&
                  PartEvib   => DSMC_Solution(8,iElem,iSpec)   ,&
                  PartErot   => DSMC_Solution(9,iElem,iSpec)   ,&
                  PartEelec  => DSMC_Solution(10,iElem,iSpec)  ,&
                  SimPartNum => DSMC_Solution(11,iElem,iSpec)  ,&
                  Macro_Velo     => DSMC_MacroVal(nVarLoc*(iSpec-1)+1:nVarLoc*(iSpec-1)+3,iElem) ,&
                  Macro_Temp     => DSMC_MacroVal(nVarLoc*(iSpec-1)+4:nVarLoc*(iSpec-1)+6,iElem) ,&
                  Macro_TempMean => DSMC_MacroVal(nVarLoc*(iSpec-1)+12,iElem)                    ,&
                  Macro_Density  => DSMC_MacroVal(nVarLoc*(iSpec-1)+7,iElem)                     ,&
                  Macro_TempVib  => DSMC_MacroVal(nVarLoc*(iSpec-1)+8,iElem)                     ,&
                  Macro_TempRot  => DSMC_MacroVal(nVarLoc*(iSpec-1)+9,iElem)                     ,&
                  Macro_Tempelec => DSMC_MacroVal(nVarLoc*(iSpec-1)+10,iElem)                    ,&
                  Macro_PartNum  => DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem) &
                  )
        IF (PartNum.GT.0.0) THEN
          ! simulation particle number
          Macro_PartNum = PartNum / REAL(DSMC%SampNum)
          ! compute flow velocity
          Macro_Velo = PartVelo / PartNum
          ! compute flow Temperature
          Macro_Temp = Species(iSpec)%MassIC/BoltzmannConst * ( (PartVelo2/PartNum)- (PartVelo/PartNum)**2 )
          ! mean flow Temperature
          Macro_TempMean = (Macro_Temp(1) + Macro_Temp(2) + Macro_Temp(3)) / 3.
          ! compute number density
          IF (SimVolume.GT.0) THEN
            IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
              ! PartNum contains the weighted particle number
              Macro_Density = Macro_PartNum / SimVolume
            ELSE
              Macro_Density = Macro_PartNum*Species(iSpec)%MacroParticleFactor /SimVolume
            END IF
          ELSE
            Macro_Density = 0.
          END IF
          IF (DoVirtualCellMerge) THEN
            IF (VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
              IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
                ! PartNum contains the weighted particle number
                Macro_Density = Macro_PartNum / VirtMergedCells(iElem)%MergedVolume
              ELSE
                Macro_Density = Macro_PartNum*Species(iSpec)%MacroParticleFactor /VirtMergedCells(iElem)%MergedVolume
              END IF
            END IF
          END IF
          ! Compute total values for a gas mixture (nSpecies > 1)
          IF(nSpecies.GT.1) THEN
            Total_PartNum   = Total_PartNum + Macro_PartNum
            Total_Velo      = Total_Velo + Macro_Velo*Macro_PartNum
            Total_Temp      = Total_Temp + Macro_Temp*Macro_PartNum
            Total_TempMean  = Total_TempMean + Macro_TempMean*Macro_PartNum
            Total_Density   = Total_Density + Macro_Density
          END IF
          ! compute internal energies / has to be changed for vfd
          IF(useDSMC)THEN
            IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
              IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
                IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                  IF( (PartEvib/PartNum) .GT. 0.0 ) THEN
                    Macro_TempVib = CalcTVibPoly(PartEvib/PartNum + SpecDSMC(iSpec)%EZeroPoint, iSpec)
                  ELSE
                    Macro_TempVib = 0.0
                  END IF
                ELSE
                  TVib_TempFac = PartEvib / (PartNum * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib)
                  IF ((PartEvib /PartNum).LE.0.0) THEN
                    Macro_TempVib = 0.0
                  ELSE
                    Macro_TempVib = SpecDSMC(iSpec)%CharaTVib / LOG(1. + 1./(TVib_TempFac))
                  END IF
                END IF
                Macro_TempRot = 2. * PartERot / (PartNum*BoltzmannConst*REAL(SpecDSMC(iSpec)%Xi_Rot))
                MolecPartNum = MolecPartNum + Macro_PartNum
                IF(nSpecies.GT.1) THEN
                  Total_TempVib  = Total_TempVib  + Macro_TempVib*Macro_PartNum
                  Total_TempRot  = Total_TempRot  + Macro_TempRot*Macro_PartNum
                END IF
              END IF
              IF (DSMC%ElectronicModel.GT.0) THEN
                IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
                  Macro_TempElec = CalcTelec(PartEelec/PartNum, iSpec)
                  HeavyPartNum = HeavyPartNum + Macro_PartNum
                END IF
                IF(nSpecies.GT.1) Total_TempElec = Total_TempElec + Macro_TempElec*Macro_PartNum
              END IF
            END IF
          END IF
        END IF
        ! Write out background gas distribution but do not include it in the totals
        IF(BGGas%UseDistribution) THEN
          IF(BGGas%BackgroundSpecies(iSpec)) THEN
            bgSpec = BGGas%MapSpecToBGSpec(iSpec)
            Macro_Velo     = BGGas%Distribution(bgSpec,1:3,iElem)
            Macro_Temp     = BGGas%Distribution(bgSpec,4:6,iElem)
            Macro_TempMean = SUM(BGGas%Distribution(bgSpec,4:6,iElem)) / 3.
            Macro_Density  = BGGas%Distribution(bgSpec,7,iElem)
            Macro_TempVib  = BGGas%Distribution(bgSpec,8,iElem)
            Macro_TempRot  = BGGas%Distribution(bgSpec,9,iElem)
            Macro_Tempelec = BGGas%Distribution(bgSpec,10,iElem)
          END IF
        END IF
      END ASSOCIATE
    END DO
    ! Compute total values for a gas mixture (nSpecies > 1)
    IF(nSpecies.GT.1) THEN
      IF (Total_PartNum.GT.0.0) THEN
        Total_Velo = Total_Velo / Total_PartNum
        Total_Temp = Total_Temp / Total_PartNum
        Total_TempMean = Total_TempMean / Total_PartNum
        IF(useDSMC)THEN
          IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.(MolecpartNum.GT.0))THEN
            Total_TempVib = Total_TempVib / MolecPartNum
            Total_TempRot = Total_TempRot / MolecPartNum
          END IF
          IF ((DSMC%ElectronicModel.GT.0).AND.(HeavyPartNum.GT.0)) THEN
            Total_TempElec = Total_TempElec / HeavyPartNum
          END IF
        END IF
      END IF
    END IF
    ! Radial weighting, vMPF, variable timestep: Getting the actual number of simulation particles without weighting factors
    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
      Total_PartNum = 0.0
      DO iSpec = 1, nSpecies
        IF(DSMC%SampNum.GT.0) DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem) = DSMC_Solution(11,iElem, iSpec) / REAL(DSMC%SampNum)
        IF(nSpecies.GT.1) Total_PartNum = Total_PartNum + DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem)
      END DO
    END IF
  END ASSOCIATE
END DO
! write dsmc quality values
IF (DSMC%CalcQualityFactors) THEN
  IF(WriteMacroVolumeValues) THEN
    iter_loc = REAL(DSMC%SampNum)
  ELSE
    IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
      iter_loc = REAL(iter)
    ELSE
      iter_loc = (Time-(1-DSMC%TimeFracSamp)*TEnd) / dt
    END IF
  END IF
  DO iElem=1,nElems
    IF (DoVirtualCellMerge) THEN
      IF (VirtMergedCells(iElem)%isMerged) CYCLE
    END IF
    nVarCount = nVar
    IF(DSMC%QualityFacSamp(iElem,4).GT.0.0) THEN
      DSMC_MacroVal(nVarCount+1:nVarCount+3,iElem) = DSMC%QualityFacSamp(iElem,1:3) / DSMC%QualityFacSamp(iElem,4)
    END IF
    nVarCount = nVar + 3
    IF(UseVarTimeStep) THEN
      IF(VarTimeStep%UseLinearScaling.OR.VarTimeStep%UseDistribution) THEN
        IF(VarTimeStep%UseLinearScaling.AND.(Symmetry%Order.EQ.2)) THEN
          ! 2D/Axisymmetric uses a scaling of the time step per particle, no element values are used. For the output simply the cell
          ! midpoint is used to calculate the time step
          VarTimeStep%ElemFac(iElem) = GetParticleTimeStep(ElemMidPoint_Shared(1,GetCNElemID(iElem + offsetElem)), &
                                                      ElemMidPoint_Shared(2,GetCNElemID(iElem + offsetElem)))
        END IF
        DSMC_MacroVal(nVarCount+1,iElem) = VarTimeStep%ElemFac(iElem)
        nVarCount = nVarCount + 1
      END IF
    END IF
    IF(DoVirtualCellMerge)THEN
      DSMC_MacroVal(nVarCount+1,iElem) = VirtMergedCells(iElem)%MasterCell
      nVarCount = nVarCount + 1
    END IF
    IF(RadialWeighting%PerformCloning) THEN
      IF(DSMC%QualityFacSamp(iElem,4).GT.0.0) THEN
        DSMC_MacroVal(nVarCount+1:nVarCount+2,iElem)=DSMC%QualityFacSamp(iElem,5:6) / DSMC%QualityFacSamp(iElem,4)
      END IF
      nVarCount = nVarCount + 2
    END IF
    IF(FPInitDone) THEN
      IF(FP_QualityFacSamp(2,iElem).GT.0) THEN
        ! Mean relaxation factor (mean over all octree subcells)
        DSMC_MacroVal(nVarCount+1,iElem) = FP_QualityFacSamp(1,iElem) / FP_QualityFacSamp(2,iElem)
        ! Mean Prandtl number
        DSMC_MacroVal(nVarCount+2,iElem) = FP_QualityFacSamp(6,iElem) / FP_QualityFacSamp(2,iElem)
      END IF
      IF(FP_QualityFacSamp(4,iElem).GT.0) THEN
        ! Max relaxation factor (maximal value of all octree subcells)
        DSMC_MacroVal(nVarCount+3,iElem) = FP_QualityFacSamp(3,iElem) / FP_QualityFacSamp(4,iElem)
        ! Max rotational relaxation factor
        DSMC_MacroVal(nVarCount+4,iElem) = FP_QualityFacSamp(5,iElem) / FP_QualityFacSamp(4,iElem)
      END IF
      ! Ratio between FP and DSMC usage per cell
      DSMC_MacroVal(nVarCount+5,iElem) = FP_QualityFacSamp(4,iElem) / iter_loc
      nVarCount = nVarCount + 5
    END IF
    IF(BGKInitDone) THEN
      IF(BGK_QualityFacSamp(2,iElem).GT.0) THEN
        ! Mean relaxation factor (mean over all octree subcells)
        DSMC_MacroVal(nVarCount+1,iElem) = BGK_QualityFacSamp(1,iElem) / BGK_QualityFacSamp(2,iElem)
        ! Mean Prandtl number
        DSMC_MacroVal(nVarCount+2,iElem) = BGK_QualityFacSamp(6,iElem) / BGK_QualityFacSamp(2,iElem)
        ! Mean expected Prandtl number
        DSMC_MacroVal(nVarCount+3,iElem) = BGK_QualityFacSamp(7,iElem) / BGK_QualityFacSamp(2,iElem)
        ! Mean viscosity
        DSMC_MacroVal(nVarCount+4,iElem) = BGK_QualityFacSamp(8,iElem) / BGK_QualityFacSamp(2,iElem)
        ! Mean thermal conductivity
        DSMC_MacroVal(nVarCount+5,iElem) = BGK_QualityFacSamp(9,iElem) / BGK_QualityFacSamp(2,iElem)
      END IF
      IF(BGK_QualityFacSamp(4,iElem).GT.0) THEN
        ! Max relaxation factor (maximal value of all octree subcells)
        DSMC_MacroVal(nVarCount+6,iElem) = BGK_QualityFacSamp(3,iElem) / BGK_QualityFacSamp(4,iElem)
        ! Max rotational relaxation factor
        DSMC_MacroVal(nVarCount+7,iElem) = BGK_QualityFacSamp(5,iElem) / BGK_QualityFacSamp(4,iElem)
      END IF
      ! Ratio between BGK and DSMC usage per cell
      DSMC_MacroVal(nVarCount+8,iElem) = BGK_QualityFacSamp(4,iElem) / iter_loc
      nVarCount = nVarCount + 8
    END IF
    ! variable rotation and vibration relaxation
    IF(Collismode.GT.1) THEN
      IF((DSMC%RotRelaxProb.GE.2).OR.(DSMC%VibRelaxProb.EQ.2)) THEN
        IF(nSpecies.EQ.1) THEN
          nSpecTemp = 0
        ELSE
          nSpecTemp = nSpecies
        END IF
        DO iSpec=0,nSpecTemp
          nVarCountRelax = 13
          IF(DSMC%RotRelaxProb.GE.2) THEN
            IF(DSMC%QualityFacSampRotSamp(iElem,iSpec+1).GT.0) THEN
              DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax,iElem) = DSMC%QualityFacSampRot(iElem,iSpec+1,2) &
                                                                / REAL(DSMC%QualityFacSampRotSamp(iElem,iSpec+1))
              DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax+1,iElem) = DSMC%QualityFacSampRot(iElem,iSpec+1,1) &
                                                                / REAL(DSMC%QualityFacSampRotSamp(iElem,iSpec+1))
              nVarCountRelax = nVarCountRelax + 2
            END IF
          END IF
          IF((DSMC%VibRelaxProb.EQ.2)) THEN
            IF(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,2).GT.0) DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax,iElem) = &
                                        DSMC%QualityFacSampVib(iElem,iSpec+1,2) / REAL(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,2))
            IF(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,1).GT.0) DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax+1,iElem)=&
                                        DSMC%QualityFacSampVib(iElem,iSpec+1,1) / REAL(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,1))
          END IF
        END DO
      END IF
    END IF
  END DO
  IF (ALLOCATED(DSMC%QualityFacSampRot)) DSMC%QualityFacSampRot = 0.
  IF (ALLOCATED(DSMC%QualityFacSampRotSamp)) DSMC%QualityFacSampRotSamp = 0
  IF (ALLOCATED(DSMC%QualityFacSampVib)) DSMC%QualityFacSampVib = 0.
  IF (ALLOCATED(DSMC%QualityFacSampVibSamp)) DSMC%QualityFacSampVibSamp = 0
END IF

IF (DoVirtualCellMerge) THEN
  DO iElem = 1, nElems
    IF (VirtMergedCells(iElem)%isMerged) THEN
      DSMC_MacroVal(:,iElem) = DSMC_MacroVal(:,VirtMergedCells(iElem)%MasterCell-offSetElem)
    END IF
  END DO
END IF

END SUBROUTINE DSMC_output_calc


SUBROUTINE CalcMacroElecExcitation(MacroElecExcitation)
!===================================================================================================================================
!> Calculation of the excitation rate [1/s] based on the sum of excitation events (sum of particle weights) and the sampling time
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Vars         ,ONLY: WriteMacroVolumeValues,WriteMacroSurfaceValues,MacroValSampTime
USE MOD_Particle_Vars         ,ONLY: ExcitationSampleData,ExcitationLevelCounter
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_Timedisc_Vars         ,ONLY: time
USE MOD_DSMC_Vars             ,ONLY: DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)      :: MacroElecExcitation(1:ExcitationLevelCounter,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! INTEGER                 ::
REAL                    :: TimeSample
!===================================================================================================================================

! Determine the sampling time for the calculation of the rate
IF (WriteMacroVolumeValues) THEN
  ! Elapsed time since last sampling (variable dt's possible!)
  TimeSample = Time - MacroValSampTime
  ! Set MacroValSampTime to the current time for the next output, BUT only if it is not used in CalcSurfaceValues, otherwise CalcSurfaceValues will be skipped
  ! TODO: Have a global calculation of the sample time before the output regardless of surface and/or volume output
  IF(.NOT.WriteMacroSurfaceValues) MacroValSampTime = Time
ELSE IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
  ! Sampling at the end of the simulation: When a restart is performed and the sampling starts immediately, determine the correct sampling time
  ! (e.g. sampling is set to 20% of tend = 1s, and restart is performed at 0.9s, sample time = 0.1s)
  TimeSample = Time - RestartTime
ELSE
  ! Sampling at the end of the simulation: calculated from the user given input
  TimeSample = (Time-(1-DSMC%TimeFracSamp)*TEnd)
END IF

! Rate
MacroElecExcitation = ExcitationSampleData / TimeSample

END SUBROUTINE CalcMacroElecExcitation

SUBROUTINE WriteDSMCToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
!> Subroutine to write the sampled macroscopic solution to the HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars     ,ONLY: DSMC, RadialWeighting, CollisMode, CollInf, SpecDSMC
USE MOD_MCC_Vars      ,ONLY: SpecXSec
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: ProjectName, ElementaryCharge
USE MOD_Mesh_Vars     ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_io_HDF5
USE MOD_HDF5_Output   ,ONLY: WriteAttributeToHDF5, WriteHDF5Header, WriteArrayToHDF5
USE MOD_Particle_Vars ,ONLY: nSpecies, UseVarTimeStep, SampleElecExcitation, ExcitationLevelCounter, DoVirtualCellMerge
USE MOD_BGK_Vars      ,ONLY: BGKInitDone
USE MOD_FPFlow_Vars   ,ONLY: FPInitDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: SpecID, LevelID, SpecID2
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNamesElecExci(:)
INTEGER                        :: nVar,nVar_quality,nVarloc,nVarCount,ALLOCSTAT, iSpec, nVarRelax, nSpecOut, iCase, iLevel
INTEGER                        :: jSpec
REAL,ALLOCATABLE               :: DSMC_MacroVal(:,:), MacroElecExcitation(:,:)
REAL                           :: StartT,EndT
!===================================================================================================================================
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DSMCState',OutputTime))//'.h5'
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DSMC TO HDF5 FILE ['//TRIM(FileName)//'] ...'
GETTIME(StartT)

IF(nSpecies.EQ.1) THEN
  nSpecOut = 1
ELSE
  nSpecOut = nSpecies + 1
END IF

! Create dataset attribute "VarNames"
nVarloc=DSMC_NVARS
nVarRelax=0
IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
  IF(DSMC%RotRelaxProb.GE.2) nVarRelax = nVarRelax + 2
  IF(DSMC%VibRelaxProb.EQ.2) nVarRelax = nVarRelax + 2
END IF

nVar=(nVarloc+nVarRelax)*nSpecOut

IF (DSMC%CalcQualityFactors) THEN
  nVar_quality=3
  IF(UseVarTimeStep) nVar_quality = nVar_quality + 1
  IF(DoVirtualCellMerge) nVar_quality = nVar_quality + 1
  IF(RadialWeighting%PerformCloning) nVar_quality = nVar_quality + 2
  IF(BGKInitDone) nVar_quality = nVar_quality + 8
  IF(FPInitDone) nVar_quality = nVar_quality + 5
ELSE
  nVar_quality=0
END IF
ALLOCATE(StrVarNames(1:nVar+nVar_quality))
nVarCount=0
IF(nSpecies.GT.1) THEN
  DO iSpec=1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    StrVarNames(nVarCount+DSMC_VELOX      )='Spec'//TRIM(SpecID)//'_VeloX'
    StrVarNames(nVarCount+DSMC_VELOY      )='Spec'//TRIM(SpecID)//'_VeloY'
    StrVarNames(nVarCount+DSMC_VELOZ      )='Spec'//TRIM(SpecID)//'_VeloZ'
    StrVarNames(nVarCount+DSMC_TEMPX      )='Spec'//TRIM(SpecID)//'_TempTransX'
    StrVarNames(nVarCount+DSMC_TEMPY      )='Spec'//TRIM(SpecID)//'_TempTransY'
    StrVarNames(nVarCount+DSMC_TEMPZ      )='Spec'//TRIM(SpecID)//'_TempTransZ'
    StrVarNames(nVarCount+DSMC_NUMDENS    )='Spec'//TRIM(SpecID)//'_NumberDensity'
    StrVarNames(nVarCount+DSMC_TVIB       )='Spec'//TRIM(SpecID)//'_TempVib'
    StrVarNames(nVarCount+DSMC_TROT       )='Spec'//TRIM(SpecID)//'_TempRot'
    StrVarNames(nVarCount+DSMC_TELEC      )='Spec'//TRIM(SpecID)//'_TempElec'
    StrVarNames(nVarCount+DSMC_SIMPARTNUM )='Spec'//TRIM(SpecID)//'_SimPartNum'
    StrVarNames(nVarCount+DSMC_TEMPMEAN   )='Spec'//TRIM(SpecID)//'_TempTransMean'
    nVarCount=nVarCount+nVarloc
    IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
      IF(DSMC%RotRelaxProb.GE.2) THEN
        StrVarNames(nVarCount+1              )='Spec'//TRIM(SpecID)//'_DSMC_MaxRotRelaxProb'
        StrVarNames(nVarCount+2              )='Spec'//TRIM(SpecID)//'_DSMC_MeanRotRelaxProb'
        nvarcount=nvarcount+2
      END IF
      IF((DSMC%VibRelaxProb.EQ.2)) THEN
        StrVarNames(nVarCount+1              )='Spec'//TRIM(SpecID)//'_DSMC_MaxVibRelaxProb'
        StrVarNames(nVarCount+2              )='Spec'//TRIM(SpecID)//'_DSMC_MeanVibRelaxProb'
        nvarcount=nvarcount+2
      END IF
    END IF
  END DO ! iSpec=1,nSpecies
END IF
! fill varnames for total values
StrVarNames(nVarCount+DSMC_VELOX      )='Total_VeloX'
StrVarNames(nVarCount+DSMC_VELOY      )='Total_VeloY'
StrVarNames(nVarCount+DSMC_VELOZ      )='Total_VeloZ'
StrVarNames(nVarCount+DSMC_TEMPX      )='Total_TempTransX'
StrVarNames(nVarCount+DSMC_TEMPY      )='Total_TempTransY'
StrVarNames(nVarCount+DSMC_TEMPZ      )='Total_TempTransZ'
StrVarNames(nVarCount+DSMC_NUMDENS    )='Total_NumberDensity'
StrVarNames(nVarCount+DSMC_TVIB       )='Total_TempVib'
StrVarNames(nVarCount+DSMC_TROT       )='Total_TempRot'
StrVarNames(nVarCount+DSMC_TELEC      )='Total_TempElec'
StrVarNames(nVarCount+DSMC_SIMPARTNUM )='Total_SimPartNum'
StrVarNames(nVarCount+DSMC_TEMPMEAN   )='Total_TempTransMean'
nVarCount=nVarCount+nVarloc
IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
  IF(DSMC%RotRelaxProb.GE.2) THEN
    StrVarNames(nVarCount+1              )='Total_DSMC_MaxRotRelaxProb'
    StrVarNames(nVarCount+2              )='Total_DSMC_MeanRotRelaxProb'
    nvarcount=nvarcount+2
  END IF
  IF((DSMC%VibRelaxProb.EQ.2)) THEN
    StrVarNames(nVarCount+1              )='Total_DSMC_MaxVibRelaxProb'
    StrVarNames(nVarCount+2              )='Total_DSMC_MeanVibRelaxProb'
    nvarcount=nvarcount+2
  END IF
END IF

IF (DSMC%CalcQualityFactors) THEN
  StrVarNames(nVarCount+1) ='DSMC_MaxCollProb'
  StrVarNames(nVarCount+2) ='DSMC_MeanCollProb'
  StrVarNames(nVarCount+3) ='DSMC_MCS_over_MFP'
  nVarCount=nVarCount+3
  IF(UseVarTimeStep) THEN
    StrVarNames(nVarCount+1) ='VariableTimeStep'
    nVarCount = nVarCount + 1
  END IF
  IF(DoVirtualCellMerge) THEN
    StrVarNames(nVarCount+1) ='MergeMasterCell'
    nVarCount = nVarCount + 1
  END IF
  IF(RadialWeighting%PerformCloning) THEN
    StrVarNames(nVarCount+1) = '2D_ClonesInCell'
    StrVarNames(nVarCount+2) = '2D_IdenticalParticles'
    nVarCount=nVarCount+2
  END IF
  IF(BGKInitDone) THEN
    StrVarNames(nVarCount+1) ='BGK_MeanRelaxationFactor'
    StrVarNames(nVarCount+2) ='BGK_MeanPrandtlNumber'
    StrVarNames(nVarCount+3) ='BGK_ExpectedPrandtlNumber'
    StrVarNames(nVarCount+4) ='BGK_Viscosity'
    StrVarNames(nVarCount+5) ='BGK_ThermalConductivity'
    StrVarNames(nVarCount+6) ='BGK_MaxRelaxationFactor'
    StrVarNames(nVarCount+7) ='BGK_MaxRotationRelaxFactor'
    StrVarNames(nVarCount+8) ='BGK_DSMC_Ratio'
    nVarCount=nVarCount+8
  END IF
  IF(FPInitDone) THEN
    StrVarNames(nVarCount+1) ='FP_MeanRelaxationFactor'
    StrVarNames(nVarCount+2) ='FP_MeanPrandtlNumber'
    StrVarNames(nVarCount+3) ='FP_MaxRelaxationFactor'
    StrVarNames(nVarCount+4) ='FP_MaxRotationRelaxFactor'
    StrVarNames(nVarCount+5) ='FP_DSMC_Ratio'
    nVarCount=nVarCount+5
  END IF
END IF

! Sampling of electronic excitation: Construct the variables name based on first and second species as well as the level threshold
IF(SampleElecExcitation) THEN
  ! Number of excitation outputs (currently only electronic)
  ALLOCATE(StrVarNamesElecExci(1:ExcitationLevelCounter))
  nVarCount = 1
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      IF(.NOT.SpecXSec(iCase)%UseElecXSec) CYCLE
      ! Output of the non-electron species as first and electron species as the second index, in case multiple electron species are defined
      IF(SpecDSMC(iSpec)%InterID.EQ.4) THEN
        WRITE(SpecID,'(I3.3)') jSpec
        WRITE(SpecID2,'(I3.3)') iSpec
      ELSE
        WRITE(SpecID,'(I3.3)') iSpec
        WRITE(SpecID2,'(I3.3)') jSpec
      END IF
      DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
        WRITE(LevelID,'(F0.2)') SpecXSec(iCase)%ElecLevel(iLevel)%Threshold/ElementaryCharge
        StrVarNamesElecExci(nVarCount)='Spec'//TRIM(SpecID)//'_Spec'//TRIM(SpecID2)//'_ExcitationRate_Elec_'//TRIM(LevelID)
        nVarCount = nVarCount + 1
      END DO
    END DO
  END DO
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteHDF5Header(TRIM('DSMCState'),File_ID)
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'NSpecies',1,IntegerScalar=nSpecies)
  ! Standard variable names
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar+nVar_quality,StrArray=StrVarNames)
  ! Additional variable names: electronic excitation rate output
  IF(SampleElecExcitation) CALL WriteAttributeToHDF5(File_ID,'VarNamesExci',ExcitationLevelCounter,StrArray=StrVarNamesElecExci)
  CALL CloseDataFile()
END IF

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! Open data file for parallel output
CALL OpenDataFile(FileName,create=.false.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

ALLOCATE(DSMC_MacroVal(1:nVar+nVar_quality,nElems), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__,' Cannot allocate output array: DSMC_MacroVal!')
END IF
CALL DSMC_output_calc(nVar,nVar_quality,nVarloc+nVarRelax,DSMC_MacroVal)

! Calculation of electronic excitation rates
IF(SampleElecExcitation) THEN
  ALLOCATE(MacroElecExcitation(1:ExcitationLevelCounter,nElems), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__,' Cannot allocate output array: MacroElecExcitation!')
  END IF
  CALL CalcMacroElecExcitation(MacroElecExcitation)
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVarX        => INT(nVar+nVar_quality,IK)       ,&
  nVarExci     => INT(ExcitationLevelCounter,IK)  ,&
  PP_nElems    => INT(PP_nElems,IK)               ,&
  offsetElem   => INT(offsetElem,IK)              ,&
  nGlobalElems => INT(nGlobalElems,IK)         )
  CALL WriteArrayToHDF5(DataSetName='ElemData' , rank=2         , &
                        nValGlobal =(/nVarX    , nGlobalElems/) , &
                        nVal       =(/nVarX    , PP_nElems/)    , &
                        offset     =(/0_IK     , offsetElem/)   , &
                        collective =.false.,  RealArray=DSMC_MacroVal(:,:))
  ! Output of electronic excitation rates in a separate container
  IF(SampleElecExcitation) THEN
    CALL WriteArrayToHDF5(DataSetName='ExcitationData' , rank=2         , &
                          nValGlobal =(/nVarExci , nGlobalElems/) , &
                          nVal       =(/nVarExci , PP_nElems/)    , &
                          offset     =(/0_IK     , offsetElem/)   , &
                          collective =.false.,  RealArray=MacroElecExcitation(:,:))
  END IF
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(DSMC_MacroVal)
IF(SampleElecExcitation) DEALLOCATE(MacroElecExcitation)
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

END SUBROUTINE WriteDSMCToHDF5


SUBROUTINE SamplingRotVibRelaxProb(iElem)
!===================================================================================================================================
!> Update sampling arrays for rotational and vibrational relaxation probability for DSMC Output
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Vars          ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: MeanProb, MaxProb, PartNum
INTEGER                       :: iSpec
!===================================================================================================================================
  IF(DSMC%RotRelaxProb.EQ.2) THEN
    MeanProb = 0
    MaxProb = 0
    PartNum = 0.
    DO iSpec=1,nSpecies
      IF(DSMC%CalcRotProb(iSpec,3).GT.0) THEN
        DSMC%QualityFacSampRot(iElem,iSpec,2) = DSMC%QualityFacSampRot(iElem,iSpec,2) + DSMC%CalcRotProb(iSpec,2)
        DSMC%QualityFacSampRot(iElem,iSpec,1) = DSMC%QualityFacSampRot(iElem,iSpec,1) &
                                              + DSMC%CalcRotProb(iSpec,1) / DSMC%CalcRotProb(iSpec,3)
        DSMC%QualityFacSampRotSamp(iElem,iSpec) = DSMC%QualityFacSampRotSamp(iElem,iSpec) + 1
        MaxProb = MAX(MaxProb,DSMC%CalcRotProb(iSpec,2))
        MeanProb = MeanProb + DSMC%CalcRotProb(iSpec,1) * DSMC%CalcRotProb(iSpec,3)
        PartNum = PartNum + DSMC%CalcRotProb(iSpec,3)
      END IF
    END DO
    IF((nSpecies.GT.1).AND.(PartNum.GT.0)) THEN
      DSMC%QualityFacSampRot(iElem,nSpecies+1,2) = DSMC%QualityFacSampRot(iElem,nSpecies+1,2) + MaxProb
      DSMC%QualityFacSampRot(iElem,nSpecies+1,1) = DSMC%QualityFacSampRot(iElem,nSpecies+1,1) + MeanProb / PartNum
      DSMC%QualityFacSampRotSamp(iElem,nSpecies+1) = DSMC%QualityFacSampRotSamp(iElem,nSpecies+1) + 1
    END IF
  END IF
  ! Sample vibration relaxation probability
  IF(DSMC%VibRelaxProb.EQ.2) THEN
    MeanProb = 0
    MaxProb = 0
    PartNum = 0
    DO iSpec=1,nSpecies
      IF(DSMC%CalcVibProb(iSpec,2).GT.0) THEN
        DSMC%QualityFacSampVib(iElem,iSpec,2) = DSMC%QualityFacSampVib(iElem,iSpec,2) + DSMC%CalcVibProb(iSpec,2)
        DSMC%QualityFacSampVibSamp(iElem,iSpec,2) = DSMC%QualityFacSampVibSamp(iElem,iSpec,2) + 1
        MaxProb = MAX(MaxProb,DSMC%CalcVibProb(iSpec,2))
      END IF
      IF(DSMC%CalcVibProb(iSpec,3).GT.0) THEN
        DSMC%QualityFacSampVib(iElem,iSpec,1) = DSMC%QualityFacSampVib(iElem,iSpec,1) &
                                              + DSMC%CalcVibProb(iSpec,1) / DSMC%CalcVibProb(iSpec,3)
        DSMC%QualityFacSampVibSamp(iElem,iSpec,1) = DSMC%QualityFacSampVibSamp(iElem,iSpec,1) + 1
        MeanProb = MeanProb + DSMC%CalcVibProb(iSpec,1) * DSMC%CalcVibProb(iSpec,3)
        PartNum = PartNum + DSMC%CalcVibProb(iSpec,3)
      END IF
    END DO
    IF((nSpecies.GT.1).AND.(PartNum.GT.0)) THEN
      IF(MaxProb.GT.0.) THEN
        DSMC%QualityFacSampVib(iElem,nSpecies+1,2) = DSMC%QualityFacSampVib(iElem,nSpecies+1,2) + MaxProb
        DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,2) = DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,2) + 1
      END IF
      DSMC%QualityFacSampVib(iElem,nSpecies+1,1) = DSMC%QualityFacSampVib(iElem,nSpecies+1,1) + MeanProb / PartNum
      DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,1) = DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,1) + 1
    END IF
  END IF
END SUBROUTINE SamplingRotVibRelaxProb


SUBROUTINE SummarizeQualityFactors(iElem)
!===================================================================================================================================
!> Sample quality factors. MCS over MFP, max and mean collision probability
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: time, TEnd
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: DSMC
USE MOD_Particle_Vars         ,ONLY: WriteMacroVolumeValues
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)    :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
  ! mean collision probability of all collision pairs
  IF(DSMC%CollProbMeanCount.GT.0) THEN
    DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
    DSMC%QualityFacSamp(iElem,2) = DSMC%QualityFacSamp(iElem,2) + DSMC%CollProbMean / REAL(DSMC%CollProbMeanCount)
  END IF
  ! mean collision separation distance of actual collisions
  IF(DSMC%CollSepCount.GT.0) DSMC%QualityFacSamp(iElem,3) = DSMC%QualityFacSamp(iElem,3) + DSMC%MCSoverMFP
  ! Counting sample size
  DSMC%QualityFacSamp(iElem,4) = DSMC%QualityFacSamp(iElem,4) + 1.
  ! Sample rotation relaxation probability
  IF((DSMC%RotRelaxProb.EQ.2).OR.(DSMC%VibRelaxProb.EQ.2)) CALL SamplingRotVibRelaxProb(iElem)
END IF

END SUBROUTINE SummarizeQualityFactors


SUBROUTINE DSMCMacroSampling()
!===================================================================================================================================
!> Check if sampling should be activated and perform sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Sampling  ,ONLY: CalcSurfaceValues
USE MOD_DSMC_Vars                   ,ONLY: DSMC
USE MOD_Restart_Vars                ,ONLY: RestartTime
USE MOD_Mesh_Vars                   ,ONLY: MeshFile
USE MOD_TimeDisc_Vars               ,ONLY: time, TEnd
USE MOD_DSMC_Vars                   ,ONLY: SamplingActive
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nOutput
!-----------------------------------------------------------------------------------------------------------------------------------
! Do not perform sampling in the case of a reservoir simulation
IF (DSMC%ReservoirSimu) RETURN

! Use user given TimeFracSamp
IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).AND.(.NOT.SamplingActive))  THEN
  SamplingActive=.TRUE.
  SWRITE(*,*)'Sampling active'
ENDIF

IF(SamplingActive) THEN
  CALL DSMC_data_sampling()
  IF(DSMC%NumOutput.NE.0) THEN
    nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
    IF(Time.GE.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
      DSMC%NumOutput = DSMC%NumOutput - 1
      ! Skipping outputs immediately after the first few iterations
      IF(RestartTime.LT.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * REAL(nOutput))) THEN
        CALL WriteDSMCToHDF5(TRIM(MeshFile),time)
        IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues(during_dt_opt=.TRUE.)
      END IF
    END IF
  END IF
END IF

END SUBROUTINE DSMCMacroSampling

END MODULE MOD_DSMC_Analyze
