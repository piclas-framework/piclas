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
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, Coll_pData, CollInf, DSMC, BGGas, ChemReac, RadialWeighting
  USE MOD_DSMC_Vars,              ONLY : ConsiderVolumePortions
  USE MOD_Particle_Vars,          ONLY : PartSpecies, Species, PartState, VarTimeStep
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_DSMC_SpecXSec,          ONLY: InterpolateCrossSection
  USE MOD_part_tools,             ONLY : GetParticleWeight
#if (PP_TimeDiscMethod==42)
  USE MOD_Particle_Vars,          ONLY : nSpecies
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                 :: iElem, iPair
  REAL(KIND=8), INTENT(IN), OPTIONAL  :: NodeVolume
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                             :: iPType, NbrOfReaction, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2
  REAL                                :: SpecNum1, SpecNum2, Weight1, Weight2, Volume
  REAL                                :: aCEX, bCEX, aMEX, bMEX, aEL, bEL, sigma_tot, MacroParticleFactor, dtCell, CollCaseNum
  REAL                                :: CollEnergy, CollProb, VeloSquare
#if (PP_TimeDiscMethod==42)
  INTEGER                             :: iReac, iSpec
#endif
!===================================================================================================================================

  iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
  iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)

  iPType = SpecDSMC(iSpec_p1)%InterID + SpecDSMC(iSpec_p2)%InterID !definition of collision case

  IF (PRESENT(NodeVolume)) THEN
    Volume = NodeVolume
  ELSE
    IF (ConsiderVolumePortions) THEN
      Volume = GEO%Volume(iElem)*(1.-GEO%MPVolumePortion(iElem))
    ELSE
      Volume = GEO%Volume(iElem)
    END IF
  END IF

  SpecNum1 = CollInf%Coll_SpecPartNum(iSpec_p1)
  SpecNum2 = CollInf%Coll_SpecPartNum(iSpec_p2)

  Weight1 = GetParticleWeight(iPart_p1)
  Weight2 = GetParticleWeight(iPart_p2)
  ! Determing the particle weight (2D/VTS: Additional scaling of the weighting according to the position within the cell)
  IF (RadialWeighting%DoRadialWeighting) THEN
    ! Correction factor: Collision pairs above the mean MPF within the cell will get a higher collision probability
    ! Not the actual weighting factor, since the weighting factor is included in SpecNum
      MacroParticleFactor = 0.5*(Weight1 + Weight2) * CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) &
                            / CollInf%MeanMPF(Coll_pData(iPair)%PairType)
    ! Sum over the mean weighting factor of all collision pairs, is equal to the number of collision pairs
    ! (incl. weighting factor)
    CollCaseNum = CollInf%MeanMPF(Coll_pData(iPair)%PairType)
  ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
    ! Not the actual weighting factor, since the weighting factor is included in SpecNum
    MacroParticleFactor = 0.5*(Weight1 + Weight2) * CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) &
                          / CollInf%MeanMPF(Coll_pData(iPair)%PairType)
    ! Sum over the mean variable time step factors (NO particle weighting factor included during MeanMPF summation)
    CollCaseNum = CollInf%MeanMPF(Coll_pData(iPair)%PairType) * Species(1)%MacroParticleFactor
    ! Weighting factor has to be included
    SpecNum1 = SpecNum1 * Species(1)%MacroParticleFactor
    SpecNum2 = SpecNum2 * Species(1)%MacroParticleFactor
  ELSE
    MacroParticleFactor = Species(1)%MacroParticleFactor
    CollCaseNum = REAL(CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType))
  END IF

  IF (VarTimeStep%UseVariableTimeStep) THEN
    dtCell = dt * (VarTimeStep%ParticleTimeStep(iPart_p1) + VarTimeStep%ParticleTimeStep(iPart_p2))*0.5
  ELSE
    dtCell = dt
  END IF

  IF (Volume.EQ.0.) iPType=-1
  SELECT CASE(iPType)

    CASE(2,3,4,5,11,12,21,22,20,30,40,6,14,24)
    ! Atom-Atom,  Atom-Mol, Mol-Mol, Atom-Atomic (non-CEX/MEX) Ion, Molecule-Atomic Ion, Atom-Molecular Ion, Molecule-Molecular Ion
    ! 5: Atom - Electron, 6: Molecule - Electron, 14: Electron - Atomic Ion, 24: Molecular Ion - Electron
      IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
        ! Using the kinetic energy of the particle (as is described in Vahedi1995 and Birdsall1991)
        VeloSquare = DOT_PRODUCT(PartState(4:6,iPart_p1),PartState(4:6,iPart_p1))
        CollEnergy = 0.5 * Species(iSpec_p1)%MassIC * VeloSquare
        ! Determining whether a real collision or a "null" collisions happens by comparing the current cross-section with the
        ! maximal collision cross section
        Coll_pData(iPair)%Prob = SQRT(VeloSquare)*InterpolateCrossSection(iSpec_p1,CollEnergy)*BGGas%NumberDensity &
                                  / SpecDSMC(iSpec_p1)%MaxCollFreq
      ELSE
        Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType))  &
                * CollInf%Cab(Coll_pData(iPair)%PairType)                                               & ! Cab species comb fac
                * MacroParticleFactor / CollCaseNum                                                     &
                * Coll_pData(iPair)%CRela2 ** (0.5-SpecDSMC(iSpec_p1)%omegaVHS) &
                        ! relative velo to the power of (1 -2omega) !! only one omega is used!!
                * dtCell / Volume
      END IF
    CASE(8) !Electron - Electron
      Coll_pData(iPair)%Prob = 0
    CASE(16) !Atom-Atomic CEX/MEX Ion
      NbrOfReaction = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p2),1)
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
          !Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2
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
      SpecNum1 = NINT(CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p1)),8) !number of particles of spec 1
      SpecNum2 = NINT(CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p2)),8) !number of particles of spec 2
      IF (Coll_pData(iPair)%CRela2.eq.0.) THEN !avoid log(0)
        Coll_pData(iPair)%Prob=0.
      ELSE
        Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType))  &
          !* CollInf%Cab(Coll_pData(iPair)%PairType)                                               & ! Cab species comb fac
          * Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MacroParticleFactor                  &
            ! weighting Fact, here only one MPF is used!!!
          / CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType)                                      & ! sum of coll cases Sab
          * 1.0E-20 * SQRT(Coll_pData(iPair)%CRela2) * sigma_tot  &
            ! CEX/MEX-relation to relative velo
          * dt / Volume                   ! timestep (should be sclaed in time disc)  divided by cell volume
      END IF !avoid log(0)
    CASE(19) !Electron - Atomic CEX/MEX Ion
      Coll_pData(iPair)%Prob = 0
    CASE (-1)
      Coll_pData(iPair)%Prob = 0.
    CASE DEFAULT
      CALL Abort(&
__STAMP__&
,'ERROR in DSMC_collis: Wrong iPType case! = ',iPType)
  END SELECT
  IF ( SQRT(DOT_PRODUCT(PartState(4:6,Coll_pData(iPair)%iPart_p1),PartState(4:6,Coll_pData(iPair)%iPart_p1)))&
    .LT.DSMC%veloMinColl(PartSpecies(Coll_pData(iPair)%iPart_p1)) .OR. &
       SQRT(DOT_PRODUCT(PartState(4:6,Coll_pData(iPair)%iPart_p2),PartState(4:6,Coll_pData(iPair)%iPart_p2)))&
    .LT.DSMC%veloMinColl(PartSpecies(Coll_pData(iPair)%iPart_p2))) Coll_pData(iPair)%Prob = 0.
  IF (ISNAN(Coll_pData(iPair)%Prob)) THEN
    IPWRITE(UNIT_errOut,*)iPair,'in',iElem,'is NaN!'
    CALL Abort(&
__STAMP__&
,'Collision probability is NaN! CRela:',RealInfoOpt=SQRT(Coll_pData(iPair)%CRela2))
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
      ! Calculate the collision probability for cross section case
      CollProb = 1. - EXP(-SQRT(VeloSquare)*InterpolateCrossSection(iSpec_p1,CollEnergy)*BGGas%NumberDensity*dt)
    ELSE
      CollProb = Coll_pData(iPair)%Prob
    END IF
    DSMC%CollProbMax = MAX(CollProb, DSMC%CollProbMax)
    DSMC%CollProbMean = DSMC%CollProbMean + CollProb
    DSMC%CollProbMeanCount = DSMC%CollProbMeanCount + 1
  END IF
#if (PP_TimeDiscMethod==42)
  ! Sum of collision probabilities for the collision pair and the corresponding reaction, required for the correct reaction rate
  IF(ChemReac%NumOfReact.GT.0) THEN
    DO iSpec=1, nSpecies
      iReac=ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p2),iSpec)
      IF (iReac.NE.0) THEN
        IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
          ! Calculate the collision probability for cross section case
          CollProb = 1. - EXP(-SQRT(VeloSquare)*InterpolateCrossSection(iSpec_p1,CollEnergy)*BGGas%NumberDensity*dt)
        ELSE
          CollProb = Coll_pData(iPair)%Prob
        END IF
        ChemReac%ReacCollMean(iReac) = ChemReac%ReacCollMean(iReac) + CollProb
        ChemReac%ReacCollMeanCount(iReac) = ChemReac%ReacCollMeanCount(iReac) + 1
      END IF
    END DO
  END IF
#endif
END SUBROUTINE DSMC_prob_calc

END MODULE MOD_DSMC_CollisionProb