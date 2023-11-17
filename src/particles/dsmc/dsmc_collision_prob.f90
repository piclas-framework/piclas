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

MODULE MOD_DSMC_CollisionProb
!===================================================================================================================================
! Module calculating the collision probability
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_prob_calc
  MODULE PROCEDURE DSMC_prob_calc
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_prob_calc
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_prob_calc(iElem, iPair, NodeVolume)
!===================================================================================================================================
! Routine calculating the collision probability
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, Coll_pData, CollInf, DSMC, BGGas, ChemReac, RadialWeighting, CollisMode
USE MOD_MCC_Vars                ,ONLY: UseMCC, SpecXSec, XSec_NullCollision
USE MOD_Particle_Vars           ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF, PartTimeStep, PartState
USE MOD_TimeDisc_Vars           ,ONLY: dt
USE MOD_MCC_XSec                ,ONLY: XSec_CalcCollisionProb, XSec_CalcReactionProb, XSec_CalcVibRelaxProb, XSec_CalcElecRelaxProb
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iElem, iPair
REAL, INTENT(IN), OPTIONAL          :: NodeVolume
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iPType, NbrOfReaction, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, PairType
REAL                                :: SpecNum1, SpecNum2, Weight1, Weight2, Volume, CollProb
REAL                                :: aCEX, bCEX, aMEX, bMEX, aEL, bEL, sigma_tot, MacroParticleFactor, dtCell, CollCaseNum
!===================================================================================================================================

PairType   = Coll_pData(iPair)%PairType
iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)

iPType = SpecDSMC(iSpec_p1)%InterID + SpecDSMC(iSpec_p2)%InterID !definition of collision case

IF (PRESENT(NodeVolume)) THEN
  Volume = NodeVolume
ELSE
  Volume = ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
END IF

SpecNum1 = CollInf%Coll_SpecPartNum(iSpec_p1)
SpecNum2 = CollInf%Coll_SpecPartNum(iSpec_p2)

Weight1 = GetParticleWeight(iPart_p1)
Weight2 = GetParticleWeight(iPart_p2)

! Determing the particle weight (2D/VTS: Additional scaling of the weighting according to the position within the cell)
IF (usevMPF) THEN
  IF(RadialWeighting%DoRadialWeighting) THEN
    ! Correction factor: Collision pairs above the mean MPF within the cell will get a higher collision probability
    ! Not the actual weighting factor, since the weighting factor is included in SpecNum
    MacroParticleFactor = 0.5*(Weight1 + Weight2) * CollInf%Coll_CaseNum(PairType) / CollInf%SumPairMPF(PairType)
  ELSE
    MacroParticleFactor = 1.
  END IF
  ! Sum over the mean weighting factor of all collision pairs, is equal to the number of collision pairs
  ! (incl. weighting factor)
  CollCaseNum = CollInf%SumPairMPF(PairType)
ELSE IF (UseVarTimeStep) THEN
  ! Not the actual weighting factor, since the weighting factor is included in SpecNum
  MacroParticleFactor = 0.5*(Weight1 + Weight2) * CollInf%Coll_CaseNum(PairType) / CollInf%SumPairMPF(PairType)
  ! Sum over the mean variable time step factors (NO particle weighting factor included during SumPairMPF summation)
  CollCaseNum = CollInf%SumPairMPF(PairType) * Species(1)%MacroParticleFactor
  ! Weighting factor has to be included
  SpecNum1 = SpecNum1 * Species(1)%MacroParticleFactor
  SpecNum2 = SpecNum2 * Species(1)%MacroParticleFactor
ELSE
  MacroParticleFactor = Species(1)%MacroParticleFactor
  CollCaseNum = REAL(CollInf%Coll_CaseNum(PairType))
END IF
IF (UseVarTimeStep) THEN
  dtCell = dt * (PartTimeStep(iPart_p1) + PartTimeStep(iPart_p2))*0.5
ELSE
  dtCell = dt
END IF

IF (Volume.EQ.0.) iPType=-1
SELECT CASE(iPType)

  CASE(2,3,4,5,11,12,21,22,20,30,40,6,14,24)
  ! Atom-Atom,  Atom-Mol, Mol-Mol, Atom-Atomic (non-CEX/MEX) Ion, Molecule-Atomic Ion, Atom-Molecular Ion, Molecule-Molecular Ion
  ! 5: Atom - Electron, 6: Molecule - Electron, 14: Electron - Atomic Ion, 24: Molecular Ion - Electron
    IF(UseMCC) THEN
      ! Coll_pData(iPair)%Prob is set inside the routine
      CALL XSec_CalcCollisionProb(iPair,iElem,SpecNum1,SpecNum2,CollCaseNum,MacroParticleFactor,Volume,dtCell)
      IF(CollisMode.EQ.3) THEN
        ! Chemical reaction with cross-section based probability
        IF(ChemReac%CollCaseInfo(iCase)%HasXSecReaction) THEN
          IF(.NOT.SpecXSec(iCase)%UseCollXSec) THEN
          ! If standard collision modelling is used, the reaction probability is added to the collision probability
            CALL XSec_CalcReactionProb(iPair,iCase,iElem,SpecNum1,SpecNum2,MacroParticleFactor,Volume)
            Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob + SUM(ChemReac%CollCaseInfo(iCase)%ReactionProb(:))
          END IF
        END IF
      END IF
      IF(SpecXSec(iCase)%UseVibXSec) THEN
        IF(.NOT.SpecXSec(iCase)%UseCollXSec) THEN
          CALL XSec_CalcVibRelaxProb(iPair,iElem,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
          Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob + SpecXSec(iCase)%VibProb
        END IF
      END IF
      IF(SpecXSec(iCase)%UseElecXSec) THEN
        IF(.NOT.SpecXSec(iCase)%UseCollXSec) THEN
          CALL XSec_CalcElecRelaxProb(iPair,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
          Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob + SUM(SpecXSec(iCase)%ElecLevel(:)%Prob)
        END IF
      END IF
    ELSE
      Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(PairType))  &
              * CollInf%Cab(PairType)                                               & ! Cab species comb fac
              * MacroParticleFactor / CollCaseNum                                                     &
              * Coll_pData(iPair)%CRela2 ** (0.5-CollInf%omega(iSpec_p1,iSpec_p2)) &
              * dtCell / Volume
    END IF
  CASE(8) !Electron - Electron
    Coll_pData(iPair)%Prob = 0
  CASE(16) !Atom-Atomic CEX/MEX Ion
    NbrOfReaction = ChemReac%ReactNum(iSpec_p1,iSpec_p2,1)
    aCEX = ChemReac%CEXa(NbrOfReaction)
    bCEX = ChemReac%CEXb(NbrOfReaction)
    aMEX = ChemReac%MEXa(NbrOfReaction)
    bMEX = ChemReac%MEXb(NbrOfReaction)
    aEL = ChemReac%ELa(NbrOfReaction)
    bEL = ChemReac%ELb(NbrOfReaction)
    IF (ChemReac%DoScat(NbrOfReaction)) THEN
      IF (Coll_pData(iPair)%CRela2.EQ.0) THEN
        sigma_tot = 0
      ELSE
        !Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(PairType)*Coll_pData(iPair)%CRela2
        sigma_tot = ((aCEX+0.5*aEL)*0.5*LOG10(Coll_pData(iPair)%CRela2)+bCEX+0.5*bEL)
        !write (*,*) "CRela2=", Coll_pData(iPair)%CRela2
        !write (*,*) "CASE(16) !Atom-Atomic CEX/MEX Ion", sigma_tot , ((aCEX+aMEX)*0.5*LOG10(Coll_pData(iPair)%CRela2) + bCEX+bMEX)
      END IF
    ELSE
      IF (Coll_pData(iPair)%CRela2.EQ.0) THEN
        sigma_tot = 0
      ELSE
        ! equation based on empirical equation, probably by J.S. Miller et al (J.Appl.Phys.91,984-991 (2002))
        ! probably: sigma_tot = sigma_CEX+sigma_MEX
        ! originally: sigma_CEX = 87.3 - 13.6*LOG10(E(in eV) with LAB energy E in electronvolt.
        ! here: equation is modified in order to use only CRela2 --> has influence on the input values
        ! --> E/1.6021766208E-19 --> sigma_CEX = -168.316 - 13.6*LOG10(E(in Joule))
        ! --> E=0.5*m_Xe_ion*CRela2 --> sigma_CEX = 171.2 - 13.6*LOG10(CRela2)
        ! --> due to the formulation of the euation in the next line, the value of aCEX, 13.6, has to be doubled in the input file
        sigma_tot = ((aCEX+aMEX)*0.5*LOG10(Coll_pData(iPair)%CRela2) + bCEX+bMEX)
      END IF
    END IF
    SpecNum1 = NINT(CollInf%Coll_SpecPartNum(iSpec_p1),8) !number of particles of spec 1
    SpecNum2 = NINT(CollInf%Coll_SpecPartNum(iSpec_p2),8) !number of particles of spec 2
    IF (Coll_pData(iPair)%CRela2.eq.0.) THEN !avoid log(0)
      Coll_pData(iPair)%Prob=0.
    ELSE
      Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(PairType))  &
        !* CollInf%Cab(PairType)                                               & ! Cab species comb fac
        * MacroParticleFactor                  &
          ! weighting Fact, here only one MPF is used!!!
        / CollInf%Coll_CaseNum(PairType)                                      & ! sum of coll cases Sab
        * 1.0E-20 * SQRT(Coll_pData(iPair)%CRela2) * sigma_tot  &
          ! CEX/MEX-relation to relative velo
        * dt / Volume                   ! timestep (should be sclaed in time disc)  divided by cell volume
    END IF !avoid log(0)
  CASE(19) !Electron - Atomic CEX/MEX Ion
    Coll_pData(iPair)%Prob = 0
  CASE (-1)
    Coll_pData(iPair)%Prob = 0.
  CASE DEFAULT
    CALL Abort(__STAMP__,'ERROR in DSMC_collis: Wrong iPType case! = ',iPType)
END SELECT

IF (ISNAN(Coll_pData(iPair)%Prob)) THEN
  IPWRITE(UNIT_errOut,*)iPair,'in',iElem,'is NaN!'
  IPWRITE(UNIT_errOut,*) '-----1.PartState--of--',iPart_p1,'-----------'
  IPWRITE(UNIT_errOut,*)PartState(1:6,iPart_p1)
  IPWRITE(UNIT_errOut,*) '-----2.PartState--of--',iPart_p2,'-----------'
  IPWRITE(UNIT_errOut,*)PartState(1:6,iPart_p2)
  CALL Abort(__STAMP__,'Collision probability is NaN! CRela:',RealInfoOpt=SQRT(Coll_pData(iPair)%CRela2))
END IF
IF(DSMC%CalcQualityFactors) THEN
  CollProb = Coll_pData(iPair)%Prob
  DSMC%CollProbMax = MAX(CollProb, DSMC%CollProbMax)
  ! Remove the correction factor for the mean collision probability
  IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
    IF(BGGas%BackgroundSpecies(iSpec_p2)) THEN
      IF(XSec_NullCollision) THEN
        IF(BGGas%UseDistribution) THEN
          CollProb = CollProb * SpecXSec(iCase)%ProbNullElem(iElem)
        ELSE
          CollProb = CollProb * SpecXSec(iCase)%ProbNull
        END IF
      ELSE
        IF(BGGas%UseDistribution)THEN
          CollProb = CollProb * BGGas%SpeciesFractionElem(BGGas%MapSpecToBGSpec(iSpec_p2),iElem)
        ELSE
          CollProb = CollProb * BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(iSpec_p2))
        END IF ! BGGas%UseDistribution
      END IF
    END IF
  END IF
  DSMC%CollProbMean = DSMC%CollProbMean + CollProb
  DSMC%CollProbMeanCount = DSMC%CollProbMeanCount + 1
END IF

IF (DSMC%ReservoirSimu) THEN
  ! Sum of collision probabilities for the collision pair, required for the correct reaction rate
  IF(ChemReac%NumOfReact.GT.0) THEN
    IF (ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths.GT.0) THEN
      CollProb = Coll_pData(iPair)%Prob
      IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
        ! Calculate the collision probability for the null collision probability case
        IF(BGGas%BackgroundSpecies(iSpec_p2)) THEN
          IF(XSec_NullCollision) THEN
            IF(BGGas%UseDistribution) THEN
              CollProb = CollProb * SpecXSec(iCase)%ProbNullElem(iElem)
            ELSE
              CollProb = CollProb * SpecXSec(iCase)%ProbNull
            END IF
          ELSE
            IF(BGGas%UseDistribution)THEN
              CollProb = CollProb * BGGas%SpeciesFractionElem(BGGas%MapSpecToBGSpec(iSpec_p2),iElem)
            ELSE
              CollProb = CollProb * BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(iSpec_p2))
            END IF ! BGGas%UseDistribution
          END IF
        END IF
      END IF
      ChemReac%ReacCollMean(iCase) = ChemReac%ReacCollMean(iCase) + CollProb
    END IF
  END IF
END IF

END SUBROUTINE DSMC_prob_calc

END MODULE MOD_DSMC_CollisionProb
