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

MODULE MOD_vmpf_collision
!===================================================================================================================================
! module including collisions with different particle weighting factors
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DSMC_vmpf_prob
  MODULE PROCEDURE DSMC_vmpf_prob
END INTERFACE
INTERFACE vMPF_PostVelo
  MODULE PROCEDURE vMPF_PostVelo
END INTERFACE
INTERFACE vMPF_AfterSplitting
  MODULE PROCEDURE vMPF_AfterSplitting
END INTERFACE
INTERFACE AtomRecomb_vMPF
  MODULE PROCEDURE AtomRecomb_vMPF
END INTERFACE

PUBLIC :: DSMC_vmpf_prob, vMPF_PostVelo, vMPF_AfterSplitting, AtomRecomb_vMPF
!===================================================================================================================================

CONTAINS


SUBROUTINE DSMC_vmpf_prob(iElem, iPair, NodeVolume)
!===================================================================================================================================
! collision probability calculation with different particle weighting factors
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  USE MOD_Globals
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, Coll_pData, CollInf, DSMC, BGGas
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartMPF, Species
  USE MOD_Particle_Mesh_Vars,     ONLY : Geo
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_DSMC_SpecXSec

!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                 !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                             :: iPType, SpecToExec
REAL                                :: SpecNum1, SpecNum2
REAL                                :: Volume
REAL                                :: MaxMPF, BGGasDensity_new
!--------------------------------------------------------------------------------------------------!
! input variable declaration                                                                       !
INTEGER, INTENT(IN)                :: iElem, iPair
REAL,INTENT(IN), OPTIONAL         :: NodeVolume
!--------------------------------------------------------------------------------------------------!

print*,'1'
read*

  iPType = SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID &
         + SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%InterID !definition of col case
  IF (PRESENT(NodeVolume)) THEN
    Volume = NodeVolume
  ELSE
    Volume = GEO%Volume(iElem)
  END IF
  IF (BGGas%BGGasSpecies.NE.0) THEN
    MaxMPF = PartMPF(Coll_pData(iPair)%iPart_p1)
  ELSE
    MaxMPF = MAX(PartMPF(Coll_pData(iPair)%iPart_p1), PartMPF(Coll_pData(iPair)%iPart_p2))
  END IF
  SELECT CASE(iPType)
    CASE(2,3,4,11) !Atom-Atom,  Atom-Mol, Mol-Mol, Atom-Atomic Ion
      SpecNum1 = CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p1)) !number of particles of spec 1
      SpecNum2 = CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p2)) !number of particles of spec 2
      IF (BGGas%BGGasSpecies.NE.0) THEN
        IF (Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensityFileID.GT.0) THEN
          BGGasDensity_new=Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensity(iElem)
        ELSE
          BGGasDensity_new=BGGas%BGGasDensity
        END IF
        Coll_pData(iPair)%Prob = BGGasDensity_new/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType)) &
                * CollInf%Cab(Coll_pData(iPair)%PairType)                                               & ! Cab species comb fac
                * Coll_pData(iPair)%CRela2 ** (0.5-SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS) &
                * dt
      ELSE
        Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType))  &
                * CollInf%Cab(Coll_pData(iPair)%PairType)                                               & ! Cab species comb fac
                * MaxMPF                                                                                &
                / CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType)                                      & ! sum of coll cases Sab
                * Coll_pData(iPair)%CRela2 ** (0.5-SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS) &
                * dt / Volume                     ! timestep (should be sclaed in time disc)  divided by cell volume
      END IF
    CASE(5) !Atom - Electron
      ALLOCATE(Coll_pData(iPair)%Sigma(0:3))  ! Cross Section of Collision of this pair
      Coll_pData(iPair)%Sigma = 0
      Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2

      ! Define what spezies is the atom and is execuded
      IF (SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.eq.1) THEN
        SpecToExec = PartSpecies(Coll_pData(iPair)%iPart_p1)
      ELSE
        SpecToExec = PartSpecies(Coll_pData(iPair)%iPart_p2)
      END IF

      SELECT CASE(SpecDSMC(SpecToExec)%NumOfPro)        ! Number of protons, which element
        CASE (18)                                      ! Argon
          CALL XSec_Argon_DravinLotz(SpecToExec, iPair)
        CASE DEFAULT
          CALL abort(&
          __STAMP__&
          ,'ERROR: spec proton not defined!')
      END SELECT

      SpecNum1 = CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p1)) !number of particles of spec 1
      SpecNum2 = CollInf%Coll_SpecPartNum(PartSpecies(Coll_pData(iPair)%iPart_p2)) !number of particles of spec 2
      ! generally this is only a HS calculation of the prob
      IF (BGGas%BGGasSpecies.NE.0) THEN
        IF (Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensityFileID.GT.0) THEN
          BGGasDensity_new=Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensity(iElem)
        ELSE
          BGGasDensity_new=BGGas%BGGasDensity
        END IF
        Coll_pData(iPair)%Prob = BGGasDensity_new  &
                * SQRT(Coll_pData(iPair)%CRela2)*Coll_pData(iPair)%Sigma(0) &
                * dt                     ! timestep (should be sclaed in time disc)
      ELSE
        Coll_pData(iPair)%Prob = SpecNum1*SpecNum2                                                      &
                * MaxMPF                                                                                &
                / CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType)                                      & ! sum of coll cases Sab
                * SQRT(Coll_pData(iPair)%CRela2)*Coll_pData(iPair)%Sigma(0)                             &
                        ! relative velo to the power of (1 -2omega) !! only one omega is used!!
                * dt / Volume                     ! timestep (should be sclaed in time disc)  divided by cell volume
      END IF
    CASE(8) !Electron - Electron
      Coll_pData(iPair)%Prob = 0
    CASE(14) !Electron - Atomic Ion
      Coll_pData(iPair)%Prob = 0
    CASE(20) !Atomic Ion - Atomic Ion
      Coll_pData(iPair)%Prob = 0
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'ERROR in DSMC_collis: Wrong iPType case! =', iPType)
  END SELECT
  IF(DSMC%CalcQualityFactors) THEN
    DSMC%CollProbMax = MAX(Coll_pData(iPair)%Prob, DSMC%CollProbMax)
    DSMC%CollProbMean = DSMC%CollProbMean + Coll_pData(iPair)%Prob
    DSMC%CollProbMeanCount = DSMC%CollProbMeanCount + 1
  END IF
END SUBROUTINE DSMC_vmpf_prob

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE vMPF_PostVelo(iPair, iElem)
!===================================================================================================================================
! determination of velocities after collision
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : DSMC_RHS, Coll_pData
  USE MOD_Particle_Vars,          ONLY : PartMPF, PartSpecies, Species
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                 !
!--------------------------------------------------------------------------------------------------!
! Local variable declaration
  REAL                                ::  Phi
!--------------------------------------------------------------------------------------------------!
! input variable declaration                                                                       !
  INTEGER, INTENT(IN)                 :: iPair, iElem                                            !
!--------------------------------------------------------------------------------------------------!

print*,'2'
read*

  IF(PartMPF(Coll_pData(iPair)%iPart_p1).gt.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
    Phi = PartMPF(Coll_pData(iPair)%iPart_p2) / PartMPF(Coll_pData(iPair)%iPart_p1)
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + 0.5 * PartMPF(Coll_pData(iPair)%iPart_p1) &
                          * Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC &
                          * Phi * (1 - Phi) &
                          * ( DSMC_RHS(Coll_pData(iPair)%iPart_p1,1)**2 &
                            + DSMC_RHS(Coll_pData(iPair)%iPart_p1,2)**2 &
                            + DSMC_RHS(Coll_pData(iPair)%iPart_p1,3)**2 )
    DSMC_RHS(Coll_pData(iPair)%iPart_p1,1) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p1,1)
    DSMC_RHS(Coll_pData(iPair)%iPart_p1,2) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p1,2)
    DSMC_RHS(Coll_pData(iPair)%iPart_p1,3) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p1,3)
  ELSE ! IF(PartMPF(Coll_pData(iPair)%iPart_p1).lt.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
    Phi = PartMPF(Coll_pData(iPair)%iPart_p1) / PartMPF(Coll_pData(iPair)%iPart_p2)
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + 0.5 * PartMPF(Coll_pData(iPair)%iPart_p2) &
                          * Species(PartSpecies(Coll_pData(iPair)%iPart_p2))%MassIC &
                          * Phi * (1 - Phi) &
                          * ( DSMC_RHS(Coll_pData(iPair)%iPart_p2,1)**2 &
                            + DSMC_RHS(Coll_pData(iPair)%iPart_p2,2)**2 &
                            + DSMC_RHS(Coll_pData(iPair)%iPart_p2,3)**2 )
    DSMC_RHS(Coll_pData(iPair)%iPart_p2,1) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p2,1)
    DSMC_RHS(Coll_pData(iPair)%iPart_p2,2) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p2,2)
    DSMC_RHS(Coll_pData(iPair)%iPart_p2,3) = Phi * DSMC_RHS(Coll_pData(iPair)%iPart_p2,3)
  END IF

END SUBROUTINE vMPF_PostVelo

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE vMPF_AfterSplitting(OrgPartIndex, W_Part, W_Spec)
!===================================================================================================================================
! Optional splitting of particle after collision arcording its new MPF
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,              ONLY : DSMC_RHS, DSMCSumOfFormedParticles, PartStateIntEn, DSMC
  USE MOD_Particle_Vars,          ONLY : PartMPF, PartSpecies, PartState, PDM, PEM
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                  !
!--------------------------------------------------------------------------------------------------!
! Local variable declaration
  INTEGER                             :: NumOfPart, iPart, PositionNbr
!--------------------------------------------------------------------------------------------------!
! input variable declaration                                                                       !
  REAL, INTENT(IN)                    :: W_Part, W_Spec                                          !
  INTEGER, INTENT(IN)                 :: OrgPartIndex
!--------------------------------------------------------------------------------------------------!


print*,'3'
read*

  NumOfPart = INT(W_Part / W_Spec)
  PartMPF(OrgPartIndex) = W_Part + (1 - NumOfPart) * W_Spec
  NumOfPart = NumOfPart - 1 ! first Part already exist (OrgPartIndex)
  DO iPart = 1, NumOfPart
  !.... Get free particle index for new part
    DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
    PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
    IF (PositionNbr.EQ.0) THEN
      CALL abort(&
      __STAMP__&
      ,'New Particle Number greater max Part Num')
    END IF
  ! Copy molecule data for non-reacting particle part
    PDM%ParticleInside(PositionNbr) = .true.
    PartSpecies(PositionNbr)        = PartSpecies(OrgPartIndex)
    PartState(PositionNbr,1:6)      = PartState(OrgPartIndex,1:6)
    PartStateIntEn(PositionNbr, 1)  = PartStateIntEn(OrgPartIndex, 1)
    PartStateIntEn(PositionNbr, 2)  = PartStateIntEn(OrgPartIndex, 2)
    PEM%Element(PositionNbr)        = PEM%Element(OrgPartIndex)
    PartMPF(PositionNbr)            = W_Spec
    IF (DSMC%ElectronicModel) THEN
      PartStateIntEn(PositionNbr, 3)  = PartStateIntEn(OrgPartIndex, 3)
    END IF

    DSMC_RHS(PositionNbr,1) = DSMC_RHS(OrgPartIndex,1)
    DSMC_RHS(PositionNbr,2) = DSMC_RHS(OrgPartIndex,2)
    DSMC_RHS(PositionNbr,3) = DSMC_RHS(OrgPartIndex,3)
  END DO

END SUBROUTINE vMPF_AfterSplitting

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE AtomRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
!===================================================================================================================================
! atom recombination routine           A + B + X -> AB + X
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles
  USE MOD_DSMC_Vars,             ONLY : ChemReac, PartStateIntEn
  USE MOD_Particle_Vars,         ONLY : PartSpecies, PartState, PDM, PEM
  USE MOD_Particle_Vars,         ONLY : usevMPF, PartMPF, Species
  USE MOD_Particle_Mesh_Vars,    ONLY : GEO
  USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange
  USE MOD_part_tools,             ONLY : DiceUnitVector
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE                                                                                      !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                      !
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  REAL                          :: FakXi, Xi, iRan
  INTEGER                       :: iQuaMax, iQua, React1Inx, React2Inx, NonReacPart, NonReacPart2
  REAL                          :: MaxColQua, Phi
  REAL                          :: ReacMPF, PartStateIntEnTemp, DeltaPartStateIntEn
  REAL                          :: MPartStateVibEnOrg, MPartStateElecEnOrg ! inner energy of M molec before Reaction
  REAL                          :: DSMC_RHS_M_Temp(3)
  REAL                          :: RanVec(3)
!--------------------------------------------------------------------------------------------------!
! input variable declaration
  INTEGER, INTENT(IN)           :: iPair, iReac, iPart_p3
  INTEGER, INTENT(IN)           :: iElem
!--------------------------------------------------------------------------------------------------!

print*,'4'
read*


  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF

  DSMC_RHS(React1Inx,1) = 0.0
  DSMC_RHS(React1Inx,2) = 0.0
  DSMC_RHS(React1Inx,3) = 0.0
  DSMC_RHS(iPart_p3,1) = 0.0
  DSMC_RHS(iPart_p3,2) = 0.0
  DSMC_RHS(iPart_p3,3) = 0.0

  ReacMPF = MIN(PartMPF(React1Inx), PartMPF(React2Inx), PartMPF(iPart_p3))
  MPartStateVibEnOrg = PartStateIntEn(iPart_p3,1)
  MPartStateElecEnOrg = PartStateIntEn(iPart_p3,3)

  IF (PartMPF(React1Inx).GT.ReacMPF) THEN ! just a part of the atom A recomb
  !.... Get free particle index for the non-reacting particle part
    DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
    NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
    IF (NonReacPart.EQ.0) THEN
      CALL abort(&
      __STAMP__&
      ,'New Particle Number greater max Part Num')
    END IF
  ! Copy molecule data for non-reacting particle part
    PDM%ParticleInside(NonReacPart) = .true.
    PartSpecies(NonReacPart)        = PartSpecies(React1Inx)
    PartState(NonReacPart,1:6)      = PartState(React1Inx,1:6)
    PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React1Inx, 1)
    PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React1Inx, 2)
    IF (DSMC%ElectronicModel) THEN
      PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React1Inx, 3)
    END IF
    PEM%Element(NonReacPart)        = PEM%Element(React1Inx)
    PartMPF(NonReacPart)            = PartMPF(React1Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
    PartMPF(React1Inx)              = ReacMPF ! reacting part MPF = ReacMPF
    IF (PartMPF(React2Inx).GT.ReacMPF) THEN ! just a part of the atom B recomb
    !.... Get free particle index for the non-reacting particle part
      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
      NonReacPart2 = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
      IF (NonReacPart2.EQ.0) THEN
        CALL abort(&
        __STAMP__&
        ,'New Particle Number greater max Part Num')
      END IF
    ! Copy molecule data for non-reacting particle part
      PDM%ParticleInside(NonReacPart2) = .true.
      PartSpecies(NonReacPart2)        = PartSpecies(React2Inx)
      PartState(NonReacPart2,1:6)      = PartState(React2Inx,1:6)
      PartStateIntEn(NonReacPart2, 1)  = PartStateIntEn(React2Inx, 1)
      PartStateIntEn(NonReacPart2, 2)  = PartStateIntEn(React2Inx, 2)
      IF (DSMC%ElectronicModel) THEN
        PartStateIntEn(NonReacPart2, 3)  = PartStateIntEn(React2Inx, 3)
      END IF
      PEM%Element(NonReacPart2)        = PEM%Element(React2Inx)
      PartMPF(NonReacPart2)            = PartMPF(React2Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
      PartMPF(React2Inx)               = ReacMPF ! reacting part MPF = ReacMPF
      CALL DSMC_RelaxForNonReacPart(iPair, NonReacPart, NonReacPart2, iElem)
    ELSE IF (PartMPF(iPart_p3).GT.ReacMPF) THEN ! just a part of the M take place in recomb
      CALL DSMC_RelaxForNonReacPart(iPair, NonReacPart, iPart_p3, iElem)
    END IF
  ELSE IF (PartMPF(React2Inx).GT.ReacMPF) THEN ! just a part of the atom B recomb
  !.... Get free particle index for the non-reacting particle part
    DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
    NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
    IF (NonReacPart.EQ.0) THEN
      CALL abort(&
      __STAMP__&
      ,'New Particle Number greater max Part Num')
    END IF
  ! Copy molecule data for non-reacting particle part
    PDM%ParticleInside(NonReacPart) = .true.
    PartSpecies(NonReacPart)        = PartSpecies(React2Inx)
    PartState(NonReacPart,1:6)      = PartState(React2Inx,1:6)
    PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React2Inx, 1)
    PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React2Inx, 2)
    IF (DSMC%ElectronicModel) THEN
      PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React2Inx, 3)
    END IF
    PEM%Element(NonReacPart)        = PEM%Element(React2Inx)
    PartMPF(NonReacPart)            = PartMPF(React2Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
    PartMPF(React2Inx)              = ReacMPF ! reacting part MPF = ReacMPF
    IF (PartMPF(iPart_p3).GT.ReacMPF) THEN ! just a part of the M take place in recomb
      CALL DSMC_RelaxForNonReacPart(iPair, NonReacPart, iPart_p3, iElem)
    END IF
  END IF

! The input particle 1 is replaced by the product molecule, the
!     second input particle is deleted
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PDM%ParticleInside(React2Inx) = .FALSE.

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation of AB and X (if X is a molecule)
!--------------------------------------------------------------------------------------------------!
  Xi = 2.0 * (2.0 - SpecDSMC(PartSpecies(iPart_p3))%omegaVHS) + SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot &
     + SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

  ! check if atomic electron shell is modelled
  IF ( DSMC%ElectronicModel ) THEN
  ! Add heat of formation to collision energy
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac) - MPartStateVibEnOrg + &
                           PartStateIntEn(Coll_pData(iPair)%iPart_p1,3) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)
!--------------------------------------------------------------------------------------------------!
! electronic relaxation  of AB and X (if X is not an electron)
!--------------------------------------------------------------------------------------------------!
    CALL ElectronicEnergyExchange(iPair,React1Inx,FakXi,iPart_p3 )
    Coll_pData(iPair)%Ec =  Coll_pData(iPair)%Ec + MPartStateElecEnOrg
    CALL ElectronicEnergyExchange(iPair,iPart_p3,FakXi,React1Inx )
  ELSE
! Add heat of formation to collision energy
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac) - MPartStateVibEnOrg
  END IF

  MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(React1Inx))%CharaTVib)  &
            - DSMC%GammaQuant
  iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(React1Inx))%MaxVibQuant)
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(iRan * iQuaMax)
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
    !laux diss page 31
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
  END DO
  PartStateIntEn(React1Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                * SpecDSMC(PartSpecies(React1Inx))%CharaTVib
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1)+ MPartStateVibEnOrg
! X particle
  IF((SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.20)) THEN
    MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(iPart_p3))%CharaTVib)  &
              - DSMC%GammaQuant
    iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(iPart_p3))%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
     !laux diss page 31
     CALL RANDOM_NUMBER(iRan)
     iQua = INT(iRan * iQuaMax)
     CALL RANDOM_NUMBER(iRan)
    END DO
    IF (PartMPF(iPart_p3).GT.ReacMPF) THEN
!      DeltaPartStateIntEn = 0.0
      Phi = ReacMPF / PartMPF(iPart_p3)
      PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpecies(iPart_p3))%CharaTVib
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
      PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn(iPart_p3,1) + Phi * PartStateIntEnTemp
      ! searche for new vib quant
      iQua = INT(PartStateIntEnTemp/(BoltzmannConst*SpecDSMC(PartSpecies(iPart_p3))%CharaTVib) - DSMC%GammaQuant)
      CALL RANDOM_NUMBER(iRan)
      IF(iRan .LT. PartStateIntEnTemp/(BoltzmannConst &
                 * SpecDSMC(PartSpecies(iPart_p3))%CharaTVib) - DSMC%GammaQuant &
                 - DBLE(iQua)) THEN
        iQua = iQua + 1
      END IF
      PartStateIntEn(iPart_p3,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpecies(iPart_p3))%CharaTVib
      DeltaPartStateIntEn = PartMPF(iPart_p3) &
                          * (PartStateIntEnTemp - PartStateIntEn(iPart_p3,1))
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2))
      GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
    ELSE
      PartStateIntEn(iPart_p3,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpecies(iPart_p3))%CharaTVib
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,1)
    END IF
  END IF
!--------------------------------------------------------------------------------------------------!
! rotational Relaxation of AB and X (if X is a molecule)
!--------------------------------------------------------------------------------------------------!
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
  IF((SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    IF (PartMPF(iPart_p3).GT.ReacMPF) THEN
      Phi = ReacMPF / PartMPF(iPart_p3)
      PartStateIntEnTemp = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
      PartStateIntEn(iPart_p3,2) = (1-Phi) * PartStateIntEn(iPart_p3,2) &
                                                   + Phi * PartStateIntEnTemp
    ELSE
      PartStateIntEn(iPart_p3,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,2)
    END IF
  END IF
!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), &
                  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(iPart_p3)))
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(iPart_p3), &
                  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(iPart_p3)))

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(iPart_p3, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(iPart_p3, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(iPart_p3, 6)

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * Coll_pData(iPair)%Ec/ &
            CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(iPart_p3)))
  RanVec(1:3) = DiceUnitVector()
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RanVec(1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RanVec(2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RanVec(3)

  ! deltaV particle 1
  DSMC_RHS(React1Inx,1) = DSMC_RHS(React1Inx,1) + VeloMx + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = DSMC_RHS(React1Inx,2) + VeloMy + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = DSMC_RHS(React1Inx,3) + VeloMz + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

 ! deltaV particle 2
  DSMC_RHS_M_Temp(1) = VeloMx - FracMassCent1*RanVelox - PartState(iPart_p3, 4)
  DSMC_RHS_M_Temp(2) = VeloMy - FracMassCent1*RanVeloy - PartState(iPart_p3, 5)
  DSMC_RHS_M_Temp(3) = VeloMz - FracMassCent1*RanVeloz - PartState(iPart_p3, 6)

  IF(PartMPF(iPart_p3).GT.ReacMPF) THEN
    Phi = ReacMPF / PartMPF(iPart_p3)
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + 0.5 * PartMPF(iPart_p3) &
                          * Species(PartSpecies(iPart_p3))%MassIC &
                          * Phi * (1 - Phi) &
                          * ( DSMC_RHS_M_Temp(1)**2 &
                            + DSMC_RHS_M_Temp(2)**2 &
                            + DSMC_RHS_M_Temp(3)**2 )
    DSMC_RHS_M_Temp(1) = Phi * DSMC_RHS_M_Temp(1)
    DSMC_RHS_M_Temp(2) = Phi * DSMC_RHS_M_Temp(2)
    DSMC_RHS_M_Temp(3) = Phi * DSMC_RHS_M_Temp(3)
  END IF
  DSMC_RHS(iPart_p3,1) = DSMC_RHS(iPart_p3,1) + DSMC_RHS_M_Temp(1)
  DSMC_RHS(iPart_p3,2) = DSMC_RHS(iPart_p3,2) + DSMC_RHS_M_Temp(2)
  DSMC_RHS(iPart_p3,3) = DSMC_RHS(iPart_p3,3) + DSMC_RHS_M_Temp(3)

  IF((usevMPF).AND.(ReacMPF.GT.(Species(PartSpecies(React1Inx))%MacroParticleFactor))) THEN
    CALL vMPF_AfterSplitting(React1Inx, ReacMPF, Species(PartSpecies(React1Inx))%MacroParticleFactor)
  END IF

END SUBROUTINE AtomRecomb_vMPF

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_RelaxForNonReacPart(iPair, Part_1, Part_2, iElem)
!===================================================================================================================================
! perform collision for non reactive particles
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : CollInf, DSMC_RHS, DSMC, SpecDSMC, PartStateIntEn, DSMC
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, PartState, PartMPF
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_DSMC_ElectronicModel,   ONLY : ElectronicEnergyExchange
  USE MOD_part_tools,             ONLY : DiceUnitVector
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                        !
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: CRela2
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  REAL                          :: CollisionEnergy
  INTEGER                       :: PartSpec1, PartSpec2, iCase
  REAL                          :: iRan
  LOGICAL                       :: DoRot2, DoVib2, DoElec1, DoElec2  ! Check whether rot or vib relax is performed
  REAL                          :: Xi_rel, Xi, FakXi                ! Factors of DOF
  INTEGER                       :: iQuaMax, iQua                    ! Quantum Numbers
  REAL                          :: MaxColQua                        ! Max. Quantum Number
  REAL                          :: PartStateIntEnTemp, Phi, DeltaPartStateIntEn ! temp. var for inertial energy (needed for vMPF)
  REAL                          :: RanVec(3)
!--------------------------------------------------------------------------------------------------!
! input variable declaration                                                                       !
  INTEGER, INTENT(IN)           :: iPair, Part_1, Part_2
!--------------------------------------------------------------------------------------------------!
! output variable declaration
  INTEGER, INTENT(IN)           :: iElem
!--------------------------------------------------------------------------------------------------!

print*,'6'
read*


  DoRot2 = .FALSE.
  DoVib2 = .FALSE.
  DoElec1 = .FALSE.
  DoElec2 = .FALSE.

  PartSpec1 = PartSpecies(Part_1)
  PartSpec2 = PartSpecies(Part_2)
  iCase = CollInf%Coll_Case(PartSpec1, PartSpec2)
  CRela2 = (PartState(Part_1,4) - PartState(Part_2,4))**2 &
         + (PartState(Part_1,5) - PartState(Part_2,5))**2 &
         + (PartState(Part_1,6) - PartState(Part_2,6))**2
  Xi_rel = 2*(2 - SpecDSMC(PartSpec1)%omegaVHS)
    ! DOF of relative motion in VHS model, only for one omega!!
  CollisionEnergy = 0.5 * CollInf%MassRed(iCase)*CRela2
  Xi = Xi_rel ! Xi are all DOF in the collision

!--------------------------------------------------------------------------------------------------!
! Decision if Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!

  IF ( DSMC%ElectronicModel ) THEN
    ! step as TRUE
    IF ((SpecDSMC(PartSpec1)%InterID.NE.4).AND.(.NOT.SpecDSMC(PartSpec1)%FullyIonized)) THEN
      CALL RANDOM_NUMBER(iRan)
      IF ( SpecDSMC(PartSpec1)%ElecRelaxProb .GT. iRAN ) THEN
        DoElec1 = .TRUE.
      END IF
    END IF
  END IF

  IF((SpecDSMC(PartSpec2)%InterID.EQ.2).OR.(SpecDSMC(PartSpec2)%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    IF(SpecDSMC(PartSpec2)%RotRelaxProb.GT.iRan) THEN
      DoRot2 = .TRUE.
      CollisionEnergy = CollisionEnergy + PartStateIntEn(Part_2,2) ! adding rot energy to coll energy
      Xi = Xi + SpecDSMC(PartSpec2)%Xi_Rot
      IF(SpecDSMC(PartSpec1)%VibRelaxProb.GT.iRan) DoVib2 = .TRUE.
    END IF
  END IF
  IF ( DSMC%ElectronicModel ) THEN
    ! step as TRUE
    IF ((SpecDSMC(PartSpec2)%InterID.NE.4).AND.(.NOT.SpecDSMC(PartSpec2)%FullyIonized)) THEN
      IF ( SpecDSMC(PartSpec2)%ElecRelaxProb .GT. iRAN ) THEN
        DoElec2 = .TRUE.
      END IF
    END IF
  END IF
!--------------------------------------------------------------------------------------------------!
! Electronic Relaxation / Transition
!--------------------------------------------------------------------------------------------------!
                       ! moved from "Vibrational Relaxation" up here as it is used in function calls
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40
  IF ( DSMC%ElectronicModel ) THEN
    ! Relaxation of first particle
    IF ( DoElec1 ) THEN
      ! calculate energy for electronic relaxation of particle 1
      CollisionEnergy = CollisionEnergy + PartStateIntEn(Part_1,3)
      CALL ElectronicEnergyExchange(iPair,Part_1,FakXi,Part_2, iElem)
    END IF
    ! Electronic relaxation of second particle
    IF ( DoElec2 ) THEN
      ! calculate energy for electronic relaxation of particle 1
      CollisionEnergy = CollisionEnergy + PartStateIntEn(Part_2,3)
      CALL ElectronicEnergyExchange(iPair,Part_2,FakXi,Part_1, iElem)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40
  IF(DoVib2) THEN
    CollisionEnergy = CollisionEnergy + PartStateIntEn(Part_2,1) ! adding vib energy to coll energy
    MaxColQua = CollisionEnergy/(BoltzmannConst*SpecDSMC(PartSpec2)%CharaTVib)  &
              - DSMC%GammaQuant
    iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpec2)%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
     !laux diss page 31
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
    END DO
    IF (PartMPF(Part_2).GT.PartMPF(Part_1)) THEN
!      DeltaPartStateIntEn = 0.0
      Phi = PartMPF(Part_1) / PartMPF(Part_2)
      PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpec2)%CharaTVib
      CollisionEnergy = CollisionEnergy - PartStateIntEnTemp
      PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn(Part_2,1) + Phi * PartStateIntEnTemp
      ! searche for new vib quant
      iQua = INT(PartStateIntEnTemp/(BoltzmannConst*SpecDSMC(PartSpec2)%CharaTVib) - DSMC%GammaQuant)
      CALL RANDOM_NUMBER(iRan)
      IF(iRan .LT. PartStateIntEnTemp/(BoltzmannConst &
                 * SpecDSMC(PartSpec2)%CharaTVib) - DSMC%GammaQuant &
                 - DBLE(iQua)) THEN
        iQua = iQua + 1
      END IF
      PartStateIntEn(Part_2,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpec2)%CharaTVib
      DeltaPartStateIntEn = PartMPF(Part_2) &
                          * (PartStateIntEnTemp - PartStateIntEn(Part_2,1))
!      CollisionEnergy = CollisionEnergy + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2))
      GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
    ELSE
      PartStateIntEn(Part_2,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpec2)%CharaTVib
      CollisionEnergy = CollisionEnergy - PartStateIntEn(Part_2,1)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoRot2) THEN
    CALL RANDOM_NUMBER(iRan)
    IF (PartMPF(Part_2).GT.PartMPF(Part_1)) THEN
      Phi = PartMPF(Part_1) / PartMPF(Part_2)
      PartStateIntEnTemp = CollisionEnergy * (1.0 - iRan**(1.0/FakXi))
      CollisionEnergy = CollisionEnergy - PartStateIntEnTemp
      PartStateIntEn(Part_2,2) = (1-Phi) * PartStateIntEn(Part_2,2) &
                                                   + Phi * PartStateIntEnTemp
    ELSE
      PartStateIntEn(Part_2,2) = CollisionEnergy * (1.0 - iRan**(1.0/FakXi))
      CollisionEnergy = CollisionEnergy - PartStateIntEn(Part_2,2)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  FracMassCent1 = CollInf%FracMassCent(PartSpec1, iCase)
  FracMassCent2 = CollInf%FracMassCent(PartSpec2, iCase)

  ! Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(Part_1, 4) &
         + FracMassCent2 * PartState(Part_2, 4)
  VeloMy = FracMassCent1 * PartState(Part_1, 5) &
         + FracMassCent2 * PartState(Part_2, 5)
  VeloMz = FracMassCent1 * PartState(Part_1, 6) &
         + FracMassCent2 * PartState(Part_2, 6)

  ! Calculate random vec and new squared velocities
  CRela2 = 2 * CollisionEnergy/CollInf%MassRed(iCase)
  RanVec(1:3) = DiceUnitVector()
  RanVelox = SQRT(CRela2) * RanVec(1)
  RanVeloy = SQRT(CRela2) * RanVec(2)
  RanVeloz = SQRT(CRela2) * RanVec(3)

  ! deltaV particle 1
  DSMC_RHS(Part_1,1) = VeloMx + FracMassCent2*RanVelox &
          - PartState(Part_1, 4)
  DSMC_RHS(Part_1,2) = VeloMy + FracMassCent2*RanVeloy &
          - PartState(Part_1, 5)
  DSMC_RHS(Part_1,3) = VeloMz + FracMassCent2*RanVeloz &
          - PartState(Part_1, 6)
 ! deltaV particle 2
  DSMC_RHS(Part_2,1) = VeloMx - FracMassCent1*RanVelox &
          - PartState(Part_2, 4)
  DSMC_RHS(Part_2,2) = VeloMy - FracMassCent1*RanVeloy &
          - PartState(Part_2, 5)
  DSMC_RHS(Part_2,3) = VeloMz - FracMassCent1*RanVeloz &
          - PartState(Part_2, 6)

  IF(PartMPF(Part_1).gt.PartMPF(Part_2)) THEN
    Phi = PartMPF(Part_2) / PartMPF(Part_1)
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + 0.5 * PartMPF(Part_1) &
                          * Species(PartSpecies(Part_1))%MassIC &
                          * Phi * (1 - Phi) &
                          * ( DSMC_RHS(Part_1,1)**2 &
                            + DSMC_RHS(Part_1,2)**2 &
                            + DSMC_RHS(Part_1,3)**2 )
    DSMC_RHS(Part_1,1) = Phi * DSMC_RHS(Part_1,1)
    DSMC_RHS(Part_1,2) = Phi * DSMC_RHS(Part_1,2)
    DSMC_RHS(Part_1,3) = Phi * DSMC_RHS(Part_1,3)
  ELSE ! IF(PartMPF(Part_1).lt.PartMPF(Part_2)) THEN
    Phi = PartMPF(Part_1) / PartMPF(Part_2)
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + 0.5 * PartMPF(Part_2) &
                          * Species(PartSpecies(Part_2))%MassIC &
                          * Phi * (1 - Phi) &
                          * ( DSMC_RHS(Part_2,1)**2 &
                            + DSMC_RHS(Part_2,2)**2 &
                            + DSMC_RHS(Part_2,3)**2 )
    DSMC_RHS(Part_2,1) = Phi * DSMC_RHS(Part_2,1)
    DSMC_RHS(Part_2,2) = Phi * DSMC_RHS(Part_2,2)
    DSMC_RHS(Part_2,3) = Phi * DSMC_RHS(Part_2,3)
  END IF

END SUBROUTINE DSMC_RelaxForNonReacPart

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_vmpf_collision
