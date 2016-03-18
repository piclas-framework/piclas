#include "boltzplatz.h"

MODULE MOD_DSMC_ChemReact
!===================================================================================================================================
! Module for chemical reactions including calculation of probabilities and collisions
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ElecImpactIoni
  MODULE PROCEDURE ElecImpactIoni
END INTERFACE

INTERFACE ElecImpactIoniQK
  MODULE PROCEDURE ElecImpactIoniQk
END INTERFACE

INTERFACE SetMeanVibQua
  MODULE PROCEDURE SetMeanVibQua
END INTERFACE

INTERFACE MolecDissoc
  MODULE PROCEDURE MolecDissoc
END INTERFACE

INTERFACE MolecExch
  MODULE PROCEDURE MolecExch
END INTERFACE

INTERFACE AtomRecomb
  MODULE PROCEDURE AtomRecomb
END INTERFACE

INTERFACE simpleCEX
  MODULE PROCEDURE simpleCEX
END INTERFACE

INTERFACE AssoIonization
  MODULE PROCEDURE AssoIonization
END INTERFACE

INTERFACE CalcReactionProb
  MODULE PROCEDURE CalcReactionProb
END INTERFACE

INTERFACE CalcBackwardRate
  MODULE PROCEDURE CalcBackwardRate
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ElecImpactIoni, SetMeanVibQua, MolecDissoc, MolecExch, AtomRecomb, simpleCEX, AssoIonization, CalcReactionProb
PUBLIC :: CalcBackwardRate, IonRecomb, ElecImpactIoniQK
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcReactionProb(iPair,iReac,ReactionProb,iPart_p3,nPartNode,Volume)
!===================================================================================================================================
! Calculates the reaction probability for dissociation, exchange, recombination and associative ionization reactions
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_PolyAtomicModel,   ONLY : Calc_Beta_Poly
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, DSMC, SpecDSMC, PartStateIntEn, ChemReac, PolyatomMolDSMC
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, BoltzmannConst, nSpecies
  USE MOD_DSMC_Analyze,           ONLY : CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair, iReac
  INTEGER, INTENT(IN), OPTIONAL :: iPart_p3, nPartNode
  REAL, INTENT(IN), OPTIONAL    :: Volume
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)             :: ReactionProb
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: PartToExec, PartReac2, iPolyatMole
  REAL                          :: EZeroPoint_Educt, EZeroPoint_Prod, EReact 
  REAL                          :: Xi_vib1, Xi_vib2, Xi_vib3, Xi_Total
  REAL(KIND=8)                 :: BetaReaction, BackwardRate
!===================================================================================================================================

  IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
    PartToExec = Coll_pData(iPair)%iPart_p1
    PartReac2 = Coll_pData(iPair)%iPart_p2
  ELSE
    PartToExec = Coll_pData(iPair)%iPart_p2
    PartReac2 = Coll_pData(iPair)%iPart_p1
  END IF

  IF((TRIM(ChemReac%ReactType(iReac)).EQ.'R').AND.(.NOT.PRESENT(iPart_p3))) THEN
    CALL abort(&
     __STAMP__&
      ,'Optional argument (iPart_p3) is missing for the recombination reaction. Reaction: ',iReac)
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calculation of the zero-point-energies
  !---------------------------------------------------------------------------------------------------------------------------------
  EZeroPoint_Educt = 0.0
  EZeroPoint_Prod = 0.0
  ! Testing if the first reacting particle is an atom or molecule, if molecule: is it polyatomic?
  IF((SpecDSMC(PartSpecies(PartToExec))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartToExec))%InterID.EQ.20)) THEN 
    IF(SpecDSMC(PartSpecies(PartToExec))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(PartSpecies(PartToExec))%SpecToPolyArray
      EZeroPoint_Educt = EZeroPoint_Educt + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
      ! Using the mean vibrational degree of freedom within the cell
!      Xi_vib1 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
      ! Calculation of the vibrational degree of freedom for the particle 
      IF (PartStateIntEn(PartToExec,1).GT.PolyatomMolDSMC(iPolyatMole)%EZeroPoint) THEN
        Xi_vib1 = 2.*(PartStateIntEn(PartToExec,1)-PolyatomMolDSMC(iPolyatMole)%EZeroPoint)                                  &
                / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(PartToExec,1), PartSpecies(PartToExec)))
      ELSE
        Xi_vib1 = 0.0
      END IF
    ELSE
      EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartToExec))%CharaTVib
      IF(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)).GT.0.0) THEN
        Xi_vib1 = 2.*ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
              * LOG(1.0/ ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) + 1.0 )
      ELSE
        Xi_vib1 = 0.0
      END IF
    END IF
  ELSE
    Xi_vib1 = 0.0
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Testing if the second particle (particle not dissociating) is an atom or molecule, if molecule: is it polyatomic?
  IF((SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.20)) THEN  
    IF(SpecDSMC(PartSpecies(PartReac2))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(PartSpecies(PartReac2))%SpecToPolyArray
      EZeroPoint_Educt = EZeroPoint_Educt + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
        ! Adding the zero-point energy of the non-reacting molecule in the dissociation case
        EZeroPoint_Prod = EZeroPoint_Prod + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
      END IF
      ! Using the mean vibrational degree of freedom within the cell
!      Xi_vib2 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
      ! Calculation of the vibrational degree of freedom for the particle 
      IF (PartStateIntEn(PartReac2,1).GT.PolyatomMolDSMC(iPolyatMole)%EZeroPoint) THEN
        Xi_vib2 = 2.*(PartStateIntEn(PartReac2,1)-PolyatomMolDSMC(iPolyatMole)%EZeroPoint)                                  &
                / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(PartReac2,1), PartSpecies(PartReac2)))
      ELSE
        Xi_vib2 = 0.0
      END IF 
    ELSE
      EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
        ! Adding the zero-point energy of the non-reacting molecule in the dissociation case
        EZeroPoint_Prod = EZeroPoint_Prod + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
      END IF
      IF(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)).GT.0.0) THEN
        Xi_vib2 = 2.*ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)) &
        * LOG(1.0/ ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)) + 1.0 )
      ELSE
        Xi_vib2 = 0.0
      END IF
    END IF
  ELSE
    Xi_vib2 = 0.0
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    ! Testing if the third particle is an atom or molecule, if molecule: is it polyatomic?
    IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2) THEN
      IF(SpecDSMC(PartSpecies(iPart_p3))%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(PartSpecies(iPart_p3))%SpecToPolyArray
        EZeroPoint_Prod = EZeroPoint_Prod + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
        EZeroPoint_Educt = EZeroPoint_Educt + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
        ! Using the mean vibrational degree of freedom within the cell
!        Xi_vib3 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
        ! Calculation of the vibrational degree of freedom for the particle 
        IF (PartStateIntEn(iPart_p3,1).GT.PolyatomMolDSMC(iPolyatMole)%EZeroPoint) THEN
          Xi_vib3 = 2.*(PartStateIntEn(iPart_p3,1)-PolyatomMolDSMC(iPolyatMole)%EZeroPoint)                                  &
                  / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(iPart_p3,1), PartSpecies(iPart_p3)))
        ELSE
          Xi_vib3 = 0.0
        END IF
      ELSE
        EZeroPoint_Prod = EZeroPoint_Prod + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(iPart_p3))%CharaTVib
        EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(iPart_p3))%CharaTVib
        IF(ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart_p3)).GT.0.0) THEN
          Xi_vib3 = 2.0*ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart_p3)) &
                  * LOG(1.0/ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart_p3)) + 1.0)
        ELSE
          Xi_vib3 = 0.0
        END IF
      END IF
    ELSE 
       Xi_vib3 = 0.0
    END IF
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Testing if the first produced particle is an atom or molecule, if molecule: is it polyatomic?
  IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%SpecToPolyArray
      EZeroPoint_Prod = EZeroPoint_Prod + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
    ELSE
      EZeroPoint_Prod = EZeroPoint_Prod + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
    END IF
  END IF
  ! Testing if the second produced particle is an atom or molecule, if molecule: is it polyatomic?
  IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%SpecToPolyArray
      EZeroPoint_Prod = EZeroPoint_Prod + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
    ELSE
      EZeroPoint_Prod = EZeroPoint_Prod + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib
    END IF
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calculation of the reaction probability, if collision energy minus the zero-point energy of the EDUCTS is greater than the
  ! activation energy AND collision energy minus the zero-point energy of the PRODUCTS is greater than the heat of formation
  !---------------------------------------------------------------------------------------------------------------------------------
  IF(((Coll_pData(iPair)%Ec-EZeroPoint_Educt).GE.ChemReac%EActiv(iReac)) .AND. &
    ((Coll_pData(iPair)%Ec-EZeroPoint_Prod).GE.(-1*ChemReac%EForm(iReac)))) THEN
    ! Determination of the total degree of freedom
    Xi_Total = Xi_vib1 + Xi_vib2  &
            + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
            + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
            + 2.*(2.-SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
    IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
      Xi_Total = Xi_Total + 3. + Xi_vib3 + SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot
    END IF
    ! Zero-point energy of educts is removed from the collision energy utilized for the calculation of the reaction probability
    EReact = Coll_pData(iPair)%Ec - EZeroPoint_Educt
    ! Determination of the Beta coefficient (array for diatomic molecules, calculation for polyatomic)
    IF (SpecDSMC(PartSpecies(PartReac2))%PolyatomicMol                &
        .OR.SpecDSMC(PartSpecies(PartToExec))%PolyatomicMol           &
        .OR.SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol  &
        .OR.SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
      BetaReaction = Calc_Beta_Poly(iReac,Xi_Total)
    ELSE
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
        BetaReaction = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(                                                         &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)),                                          &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))
      ELSE IF(TRIM(ChemReac%ReactType(iReac)).EQ.'E') THEN
        BetaReaction = ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(                                                         &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)),                                          &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))
      ELSE IF(TRIM(ChemReac%ReactType(iReac)).EQ.'i') THEN
        BetaReaction = ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(                                                          &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)),                                          &
                              ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))
      ELSE IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        IF(SpecDSMC(PartSpecies(iPart_p3))%PolyatomicMol) THEN
          BetaReaction = Calc_Beta_Poly(iReac,Xi_Total)
        ELSE
          BetaReaction = &
            ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(PartSpecies(iPart_p3),ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart_p3)))
        END IF
      ELSE
        CALL abort(&
       __STAMP__&
        ,'Reaction Type is not properly specified. Reaction: ',iReac)
      END IF
    END IF
    ! Calculation of the backward reaction rate coefficient and applying to Beta coefficient after Boyd "Modeling backward chemical
    ! rate processes in the direct simulation Monte Carlo method", Phys. Fluids 19, 1261103 (2007)
    IF(DSMC%BackwardReacRate.AND.((iReac.GT.ChemReac%NumOfReact/2))) THEN
      CALL CalcBackwardRate(iReac,DSMC%InstantTransTemp(nSpecies+1),BackwardRate)
      BetaReaction = BetaReaction * BackwardRate &
                / (ChemReac%Arrhenius_Prefactor(iReac) * DSMC%InstantTransTemp(nSpecies+1)**ChemReac%Arrhenius_Powerfactor(iReac))
    END IF
    ! Actual calculation of the reaction probability, different equation for recombination reaction
    IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
      ReactionProb = BetaReaction * (nPartNode*Species(PartSpecies(iPart_p3))%MacroParticleFactor/Volume)    &
               * EReact**(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 + SpecDSMC(PartSpecies(iPart_p3))%omegaVHS)
    ELSE
      IF(SpecDSMC(PartSpecies(PartReac2))%PolyatomicMol.OR.SpecDSMC(PartSpecies(PartToExec))%PolyatomicMol) THEN
        ! Energy is multiplied by a factor to increase the resulting exponent and avoid floating overflows for high vibrational
        ! degree of freedom, later the reaction probability is scaled again with the same factor and the respective exponents
        ReactionProb = BetaReaction * ((EReact - ChemReac%EActiv(iReac))*1E6)                                                   &
              ** (ChemReac%Arrhenius_Powerfactor(iReac)-1.5+SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS+Xi_Total/2.)    &
               * (EReact * 1E6)**(1.0 - Xi_Total/2.)
        ReactionProb = ReactionProb / ((1E6)**(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5                                      &
              + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS))
      ELSE
        ReactionProb = BetaReaction * ((EReact - ChemReac%EActiv(iReac)))                                                       &
              ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS             &
              + Xi_Total/2.) * (EReact) ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor - Xi_Total/2.)
      END IF
    END IF
  ELSE
    ReactionProb = 0.0
  END IF
END SUBROUTINE CalcReactionProb


SUBROUTINE ElecImpactIoniQK(iReac, iPair)
!===================================================================================================================================
! Perfoms the electron impact ionization with Q-K
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, PartStateIntEn !, Debug_Energy
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec
USE MOD_Particle_Vars,         ONLY : usevMPF, Species,PartPosRef
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange, TVEEnergyExchange
USE MOD_Eval_xyz,              ONLY : eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars,ONLY : DoRefmapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec
  REAL                          :: JToEv, iRan, FacEtraDistri
  REAL                          :: ERel_React1_React2, ERel_React2_Elec
  INTEGER                       :: PositionNbr, React1Inx, ElecInx
  REAL                          :: VxPseuAtom, VyPseuAtom, VzPseuAtom
  INTEGER                       :: MaxElecQua
  REAL                          :: IonizationEnergy
  REAL                          :: FakXi, Xi_rel,Xi
  INTEGER                       :: iQuaMax, iQua, MaxColQua
  !REAL                          :: ElecTransfer
  INTEGER                       :: newSpecies
  LOGICAL                       :: DoVib, DoRot
!===================================================================================================================================


!..Get the index of react1 and the electron
IF (SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.eq.4) THEN
  ElecInx = Coll_pData(iPair)%iPart_p1
  React1Inx = Coll_pData(iPair)%iPart_p2 
ELSE
  React1Inx = Coll_pData(iPair)%iPart_p1
  ElecInx = Coll_pData(iPair)%iPart_p2
END IF

!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) +&
!      0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
!    + 0.5* Species(PartSpecies(ElecInx))%MassIC*(PartState(ElecInx,4)**2+PartState(ElecInx,5)**2+PartState(ElecInx,6)**2)&
!    + PartStateIntEn(React1Inx,3)

! ionization level is last known energy level of species
MaxElecQua=SpecDSMC(PartSpecies(React1Inx))%MaxElecQuant - 1
IonizationEnergy=SpecDSMC(PartSpecies(React1Inx))%ElectronicState(2,MaxElecQua)*BoltzmannConst
! nullify
!ElecTransfer = 0.

IF(usevMPF) CALL abort(&
       __STAMP__&
        ,' Reaction not implemented with vMPF ',iReac)

IF (SpecDSMC(PartSpecies(ElecInx))%InterID.NE.4) CALL abort(&
       __STAMP__&
        ,' Only electron impact ionization. Further collision partner not implemented!  ',iReac)


! remove ionization energy from collision
Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - IonizationEnergy

Xi_rel = 4.*(2. - SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS) 
Xi     = Xi_rel !Xi are all DOF in the collision

! if first collison partner is a molecule add its internal energy to the collision energy
PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
IF((SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.20)) THEN
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec +  PartStateIntEn(React1Inx,1) + PartStateIntEn(React1Inx,2)
  Xi = Xi + SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
  FakXi = 0.5*Xi  - 1 
  CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi)
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2) 
ELSE
  ! store the electronic energy of the ionized particle
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40  
  CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi)
  !ElecTransfer = PartStateIntEn(React1Inx,3)
END IF

! distribute Etra of pseudo neutral particle (i+e) and old electron
CALL RANDOM_NUMBER(iRan)
FacEtraDistri = iRan
CALL RANDOM_NUMBER(iRan)
! laux diss page 40, omegaVHS only of one species
DO WHILE ((4 *FacEtraDistri*(1-FacEtraDistri))**(1-SpecDSMC(PartSpecies(React1Inx))%omegaVHS).LT.iRan)
  CALL RANDOM_NUMBER(iRan)
  FacEtraDistri = iRan
  CALL RANDOM_NUMBER(iRan)
END DO
ERel_React1_React2 = Coll_pData(iPair)%Ec * FacEtraDistri
ERel_React2_Elec = Coll_pData(iPair)%Ec - ERel_React1_React2

!.... Get free particle index for the 3rd particle produced
DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
IF (PositionNbr.EQ.0) THEN
   CALL Abort(&
    __STAMP__,&
    ' New Particle Number greater than max particle number!')
END IF

!Set new Species of electron
PDM%ParticleInside(PositionNbr) = .true.
PartSpecies(PositionNbr) = ChemReac%DefinedReact(iReac,2,2)
PartState(PositionNbr,1:3) = PartState(React1Inx,1:3)
IF(DoRefMapping)THEN ! here Nearst-GP is missing
  PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,React1Inx)
END IF

PartStateIntEn(PositionNbr, 1) = 0.
PartStateIntEn(PositionNbr, 2) = 0.
PartStateIntEn(PositionNbr, 3) = 0.
PEM%Element(PositionNbr) = PEM%Element(React1Inx)

!Scattering of pseudo atom (e-i) and collision partner e (scattering of e)
FracMassCent1 = CollInf%FracMassCent(PartSpecies(ChemReac%DefinedReact(iReac,1,1)), &
    CollInf%Coll_Case(PartSpecies(ChemReac%DefinedReact(iReac,1,1)),PartSpecies(ElecInx)))
FracMassCent2 = CollInf%FracMassCent(PartSpecies(ElecInx), &
    CollInf%Coll_Case(PartSpecies(ChemReac%DefinedReact(iReac,1,1)),PartSpecies(ElecInx)))

!Calculation of velo from center of mass
VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
       + FracMassCent2 * PartState(ElecInx, 4)
VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
       + FracMassCent2 * PartState(ElecInx, 5)
VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
       + FracMassCent2 * PartState(ElecInx, 6)

!calculate random vec and new squared velocities
Coll_pData(iPair)%CRela2 = 2 * ERel_React1_React2 /CollInf%MassRed(Coll_pData(iPair)%PairType)
CALL RANDOM_NUMBER(iRan)
iVec = INT(NumRanVec * iRan + 1)
RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

! deltaV particle 2
DSMC_RHS(ElecInx,1) = VeloMx - FracMassCent1*RanVelox - PartState(ElecInx, 4)
DSMC_RHS(ElecInx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(ElecInx, 5)
DSMC_RHS(ElecInx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(ElecInx, 6)

!Set velocity of pseudo atom (i+e)
VxPseuAtom = (VeloMx + FracMassCent2*RanVelox)* Species(ChemReac%DefinedReact(iReac,1,1))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)  
VyPseuAtom = (VeloMy + FracMassCent2*RanVeloy)* Species(ChemReac%DefinedReact(iReac,1,1))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)   
VzPseuAtom = (VeloMz + FracMassCent2*RanVeloz)* Species(ChemReac%DefinedReact(iReac,1,1))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)   

!Scattering of i + e
FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), &
              &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
FracMassCent2 = CollInf%FracMassCent(PartSpecies(PositionNbr), & 
              &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))

!calculate random vec and new squared velocities
Coll_pData(iPair)%CRela2 = 2 *  ERel_React2_Elec / & 
        CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
CALL RANDOM_NUMBER(iRan)
iVec = INT(NumRanVec * iRan + 1)
RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

!deltaV particle 1
DSMC_RHS(React1Inx,1) = VxPseuAtom + FracMassCent2*RanVelox - PartState(React1Inx, 4)
DSMC_RHS(React1Inx,2) = VyPseuAtom + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
DSMC_RHS(React1Inx,3) = VzPseuAtom + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

!deltaV new formed particle
PartState(PositionNbr,4:6) = 0
DSMC_RHS(PositionNbr,1) = VxPseuAtom - FracMassCent1*RanVelox 
DSMC_RHS(PositionNbr,2) = VyPseuAtom - FracMassCent1*RanVeloy 
DSMC_RHS(PositionNbr,3) = VzPseuAtom - FracMassCent1*RanVeloz 

!Debug_Energy(2) = Debug_Energy(2)&
!      + 0.5* Species(PartSpecies(ElecInx))%MassIC * (&
!           (VeloMx - FracMassCent1*RanVelox)**2 &
!          +(VeloMy - FracMassCent1*RanVeloy)**2 &
!          +(VeloMz - FracMassCent1*RanVeloz)**2) &
!      + 0.5* Species(PartSpecies(React1Inx))%MassIC  * (&
!           (VxPseuAtom + FracMassCent2*RanVelox)**2 &
!          +(VyPseuAtom + FracMassCent2*RanVeloy)**2 &
!          +(VzPseuAtom + FracMassCent2*RanVeloz)**2) &
!      + 0.5* Species(PartSpecies(PositionNbr))%MassIC  * (&
!           (VxPseuAtom - FracMassCent1*RanVelox)**2 &
!          +(VyPseuAtom - FracMassCent1*RanVeloy)**2 &
!          +(VzPseuAtom - FracMassCent1*RanVeloz)**2) &
!        + PartStateIntEn(React1Inx,3)

!print*, Debug_Energy(1),Debug_Energy(2), Debug_Energy(2)-Debug_Energy(1)
!read*
END SUBROUTINE ElecImpactIoniQK


SUBROUTINE IonRecomb(iReac, iPair, iPart_p3)
!===================================================================================================================================
! performe three-body ion-recombination to neutral or less charged ion
! ion recombination routine           A+ + e + X -> A + X
!===================================================================================================================================
USE MOD_Globals
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles!, Debug_Energy
USE MOD_DSMC_Vars,             ONLY : ChemReac, CollisMode, PartStateIntEn
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec,Species
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange, TVEEnergyExchange
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac, iPart_p3
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz, VeloMxOld, VeloMyOld, VeloMzOld ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec
  REAL                          :: FakXi, Xi, iRan, Xi_rel
  INTEGER                       :: iQuaMax, iQua, React1Inx, React2Inx
  REAL                          :: MaxColQua
! additional for Q-K theory
  REAL                          :: ksum, Tcoll
  INTEGER                       :: ii
  INTEGER                       :: newSpeciesID,oldSpeciesID
  INTEGER                       :: MaxElecQua
  REAL                          :: IonizationEnergy, evor, enach
  LOGICAL                       :: DoVib, DoRot
!===================================================================================================================================
IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
  React1Inx = Coll_pData(iPair)%iPart_p1
  React2Inx = Coll_pData(iPair)%iPart_p2
ELSE
  React2Inx = Coll_pData(iPair)%iPart_p1
  React1Inx = Coll_pData(iPair)%iPart_p2
END IF

!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) +&
!      0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
!    + 0.5* Species(PartSpecies(React2Inx))%MassIC*(PartState(React2Inx,4)**2+PartState(React2Inx,5)**2+PartState(React2Inx,6)**2)&
!    + 0.5* Species(PartSpecies(iPart_p3))%MassIC* (PartState(iPart_p3, 4)**2+PartState(iPart_p3, 5)**2+PartState(iPart_p3, 6)**2)&
!    + PartStateIntEn(React1Inx,3)

  ! Calculation of the centre of mass of the product molecule
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(React2Inx), Coll_pData(iPair)%PairType)

  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(React2Inx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(React2Inx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(React2Inx, 6)  

! The input particle 1 is replaced by the product molecule, the
!     second input particle is deleted
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PDM%ParticleInside(React2Inx) = .FALSE.
  PartState(React1Inx, 4) = VeloMx
  PartState(React1Inx, 5) = VeloMy
  PartState(React1Inx, 6) = VeloMz

  ! has to be calculated earlier because of setting of electronic energy
  Xi = 2.0 * (2.0 - SpecDSMC(PartSpecies(iPart_p3))%omegaVHS) 
 

! ionization level is last known energy level of species
  MaxElecQua=SpecDSMC(PartSpecies(React1Inx))%MaxElecQuant - 1
  IonizationEnergy=SpecDSMC(PartSpecies(React1Inx))%ElectronicState(2,MaxElecQua)*BoltzmannConst

  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                       &
                 + 0.5 * CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx), PartSpecies(iPart_p3)))                  &
                 * ((VeloMx-PartState(iPart_p3,4))**2+(VeloMy-PartState(iPart_p3,5))**2+(VeloMz-PartState(iPart_p3,6))**2)  &
                 + IonizationEnergy + PartStateIntEn(React1Inx,3)
        
  IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ. 2) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart_p3,2)
    Xi = Xi + SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot
  END IF
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

!--------------------------------------------------------------------------------------------------!
! electronic relaxation  of AB and X (if X is not an electron)  and vibrational relaxation
!--------------------------------------------------------------------------------------------------!
  CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi)

  IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ. 2) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart_p3,3) + PartStateIntEn(iPart_p3,1)
    CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,iPart_p3,FakXi)
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(iPart_p3,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,2) 
  ELSE IF ((SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.1).OR.(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.10)) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart_p3,3)
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,iPart_p3,FakXi)
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
CALL RANDOM_NUMBER(iRan)
iVec = INT(NumRanVec * iRan + 1)
RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

! deltaV particle 1
DSMC_RHS(React1Inx,1) = VeloMx + FracMassCent2*RanVelox - PartState(React1Inx, 4)
DSMC_RHS(React1Inx,2) = VeloMy + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
DSMC_RHS(React1Inx,3) = VeloMz + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

! deltaV particle 2
DSMC_RHS(iPart_p3,1) = VeloMx - FracMassCent1*RanVelox - PartState(iPart_p3, 4)
DSMC_RHS(iPart_p3,2) = VeloMy - FracMassCent1*RanVeloy - PartState(iPart_p3, 5)
DSMC_RHS(iPart_p3,3) = VeloMz - FracMassCent1*RanVeloz - PartState(iPart_p3, 6)

!Debug_Energy(2) = Debug_Energy(2)&
!      + 0.5* Species(PartSpecies(React1Inx))%MassIC * (&
!           (VeloMx + FracMassCent2*RanVelox)**2 &
!          +(VeloMy + FracMassCent2*RanVeloy)**2 &
!          +(VeloMz + FracMassCent2*RanVeloz)**2) &
!      + 0.5* Species(PartSpecies(iPart_p3))%MassIC  * (&
!           (VeloMx - FracMassCent1*RanVelox)**2 &
!          +(VeloMy - FracMassCent1*RanVeloy)**2 &
!          +(VeloMz - FracMassCent1*RanVeloz)**2) &
!        + PartStateIntEn(React1Inx,3)

!print*, Debug_Energy(1),Debug_Energy(2), Debug_Energy(2)-Debug_Energy(1)
END SUBROUTINE IonRecomb

SUBROUTINE ElecImpactIoni(iReac, iPair)
!===================================================================================================================================
! Perfoms the electron impact ionization
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, PartStateIntEn
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec, PartPosRef
USE MOD_Particle_Tracking_Vars,ONLY : DoRefmapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec
  REAL                          :: JToEv, iRan, FacEtraDistri
  REAL                          :: ERel_React1_React2, ERel_React2_Elec
  INTEGER                       :: PositionNbr, React1Inx, ElecInx
  REAL                          :: VxPseuAtom, VyPseuAtom, VzPseuAtom
!  REAL                           :: FakXi, Xi_rel
!  INTEGER                        :: iQuaMax, iQua, MaxColQua
!===================================================================================================================================

JToEv = 1.602176565E-19

!..Get the index of react1 and the electron
  IF (SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.EQ.4) THEN
    ElecInx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2 
  ELSE
    React1Inx = Coll_pData(iPair)%iPart_p1
    ElecInx = Coll_pData(iPair)%iPart_p2
  END IF

Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - SpecDSMC(PartSpecies(React1Inx))%Eion_eV*JToEv 

!spectoexec muss evtl noch geändert werden, ist ja kein spec sondern ein partikel
!page 31 diss laux DOF

! ! first, the internal energies relax if SpecToExec is a molecule
! IF(molekül) THEN
! ! Vibrational Relaxation if molec
!   Xi_rel = 2*(2 - SpecDSMC(PartSpecies(SpecToExec))%omegaVHS) ! DOF of relative motion in VHS model, only for one omega!!
!           ! this is a result of the mean value of the relative energy in the vhs model, laux diss page 31
!   FakXi = 0.5*(Xi_rel + SpecDSMC(PartSpecies(SpecToExec))%Xi_Rot) - 1
!           ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40
!   MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(SpecToExec))%CharaTVib)  - DSMC%GammaQuant
!   iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(SpecToExec))%MaxVibQuant)
!   CALL RANDOM_NUMBER(iRan)
!   iQua = INT(iRan * iQuaMax)
!   CALL RANDOM_NUMBER(iRan)
!   DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi) 
!     !GammaQuant was added, laux diss page 31, this was not in eq in LasVegas
!     CALL RANDOM_NUMBER(iRan)
!     iQua = INT(iRan * iQuaMax)    
!     CALL RANDOM_NUMBER(iRan)
!   END DO
! !spectoexec_evib muss noch genau definiert werden!!!!!!!!!!!!!!
!    spectoexec_evib = (iQua + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(PartSpecies(SpecToExec))%CharaTVib 
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - spectoexec_evib 

! ! Rotational Relaxation if molec
! !spectoexec_erot muss noch genau definiert werden!!!!!!!!!!!!!!
!    CALL RANDOM_NUMBER(iRan)
!    spectoexec_erot = Coll_pData(iPair)%Ec * (1 - iRan**(1/FakXi))
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - spectoexec_erot 
! END IF
  
  ! distribute Etra of pseudo neutral particle (i+e) and old electron
  CALL RANDOM_NUMBER(iRan)
  FacEtraDistri = iRan
  CALL RANDOM_NUMBER(iRan)
  ! laux diss page 40, omegaVHS only of one species
  DO WHILE ((4 *FacEtraDistri*(1-FacEtraDistri))**(1-SpecDSMC(PartSpecies(React1Inx))%omegaVHS).LT.iRan)
    CALL RANDOM_NUMBER(iRan)
    FacEtraDistri = iRan
    CALL RANDOM_NUMBER(iRan)
  END DO
  ERel_React1_React2 = Coll_pData(iPair)%Ec * FacEtraDistri
  ERel_React2_Elec = Coll_pData(iPair)%Ec - ERel_React1_React2

  !.... Get free particle index for the 3rd particle produced
  DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
  PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
  IF (PositionNbr.EQ.0) THEN
    PRINT*, 'New Particle Number greater max Part Num'
    STOP
  END IF

  !Set new Species of electron
  PDM%ParticleInside(PositionNbr) = .true.
  PartSpecies(PositionNbr) = ChemReac%DefinedReact(iReac,2,2)
  PartState(PositionNbr,1:3) = PartState(React1Inx,1:3)
  IF(DoRefMapping)THEN ! here Nearst-GP is missing
    PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,React1Inx)
  END IF
  PartStateIntEn(PositionNbr, 1) = 0
  PartStateIntEn(PositionNbr, 2) = 0
  PEM%Element(PositionNbr) = PEM%Element(React1Inx)

  !Scattering of pseudo atom (e-i) and collision partner e (scattering of e)
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(ElecInx), Coll_pData(iPair)%PairType)

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(ElecInx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(ElecInx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(ElecInx, 6)

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * ERel_React1_React2 /CollInf%MassRed(Coll_pData(iPair)%PairType)
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)
  
 ! deltaV particle 2
  DSMC_RHS(ElecInx,1) = VeloMx - FracMassCent1*RanVelox - PartState(ElecInx, 4)
  DSMC_RHS(ElecInx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(ElecInx, 5)
  DSMC_RHS(ElecInx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(ElecInx, 6)
  
  !Set velocity of pseudo atom (i+e)
  VxPseuAtom = VeloMx + FracMassCent2*RanVelox  
  VyPseuAtom = VeloMy + FracMassCent2*RanVeloy 
  VzPseuAtom = VeloMz + FracMassCent2*RanVeloz 

  !Set new Species of formed ion
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)

  !Scattering of i + e
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), &
                &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(PositionNbr), & 
                &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 *  ERel_React2_Elec / & 
          CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  !deltaV particle 1
  DSMC_RHS(React1Inx,1) = VxPseuAtom + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = VyPseuAtom + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = VzPseuAtom + FracMassCent2*RanVeloz - PartState(React1Inx, 6)
  
  !deltaV new formed particle
  PartState(PositionNbr,4:6) = 0
  DSMC_RHS(PositionNbr,1) = VxPseuAtom - FracMassCent1*RanVelox 
  DSMC_RHS(PositionNbr,2) = VyPseuAtom - FracMassCent1*RanVeloy 
  DSMC_RHS(PositionNbr,3) = VzPseuAtom - FracMassCent1*RanVeloz 

END SUBROUTINE ElecImpactIoni


SUBROUTINE MolecDissoc(iReac, iPair)
!===================================================================================================================================
! Perfom the molecular dissociation
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY : abort
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, CollisMode, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar, Debug_Energy
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec
USE MOD_Particle_Vars,         ONLY : usevMPF, PartMPF, RandomVec, Species, PartPosRef
USE MOD_Particle_Mesh_Vars,    ONLY : GEO
USE MOD_vmpf_collision,        ONLY : vMPF_AfterSplitting
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange
USE MOD_DSMC_PolyAtomicModel,  ONLY : DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, DSMC_InsertPolyProduct
USE MOD_DSMC_Analyze,          ONLY : CalcTVibPoly
USE MOD_Particle_Tracking_Vars,ONLY : DoRefMapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec, iPolyatMole, iDOF
  REAL                          :: JToEv, FakXi, Xi_rel, iRan, FacEtraDistri, Xi_vib1, Xi_vib2
  REAL                          :: ERel_React1_React2, ERel_React1_React3
  INTEGER                       :: iQuaMax, iQua, PositionNbr, React1Inx, React2Inx, NonReacPart
  REAL                          :: MaxColQua
  REAL                          :: VxPseuMolec, VyPseuMolec, VzPseuMolec
  REAL                          :: DeltaPartStateIntEn, PartStateIntEnTemp, Phi, ReacMPF
  REAL                          :: ElecTransfer
  REAL                          :: EZeroTempToExec1, EZeroTempToExec2, TVibTemp
!===================================================================================================================================

JToEv = 1.602176565E-19
ElecTransfer = 0.
Xi_vib1 = 0.0
Xi_vib2 = 0.0
EZeroTempToExec1 = 0.0
EZeroTempToExec2 = 0.0

!..Get the index of react1 and the react2
  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2 
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF


!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) +&
!      0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
!    + 0.5* Species(PartSpecies(React2Inx))%MassIC*(PartState(React2Inx,4)**2+PartState(React2Inx,5)**2+PartState(React2Inx,6)**2)&
!    + PartStateIntEn(React1Inx,1) + PartStateIntEn(React2Inx,1) &
!    + PartStateIntEn(React1Inx,2) + PartStateIntEn(React2Inx,2)


  IF (usevMPF) THEN ! reaction MPF definition
    ReacMPF = MIN(PartMPF(React1Inx), PartMPF(React2Inx))
    IF (PartMPF(React1Inx).GT.ReacMPF) THEN ! just a part of the molecule diss
    !.... Get free particle index for the non-reacting particle part
      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
      NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
      IF (NonReacPart.EQ.0) THEN
        CALL abort(&
       __STAMP__&
        ,'New Particle Number greater max Part Num in MolecDiss. Reaction: ',iReac)
      END IF
    ! Copy molecule data for non-reacting particle part
      PDM%ParticleInside(NonReacPart) = .true.
      PartSpecies(NonReacPart)        = PartSpecies(React1Inx)
      PartState(NonReacPart,1:6)      = PartState(React1Inx,1:6)
      IF(DoRefMapping)THEN ! here Nearst-GP is missing
        PartPosRef(1:3,NonReacPart)=PartPosRef(1:3,React1Inx)
      END IF
      PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React1Inx, 1)
      PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React1Inx, 2)
      IF (DSMC%ElectronicState) THEN
        PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React1Inx, 3)
      END IF
      PEM%Element(NonReacPart)        = PEM%Element(React1Inx)
      PartMPF(NonReacPart)            = PartMPF(React1Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
      PartMPF(React1Inx)              = ReacMPF ! reacting part MPF = ReacMPF
    END IF
  END IF
  
  ! Add heat of formation to collision energy
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac)

  Xi_rel = 4.*(2. - SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS)

  IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN
    ! If the dissociating molecule is polyatomic, additional rotational DOFs of the products of the dissociation have to be
    ! included in FakXi
    FakXi = 0.5*(Xi_rel + SpecDSMC(PartSpecies(React2Inx))%Xi_Rot + SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot &
    + SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot) - 1.0
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Vibrational degree of freedom of first product
    !-------------------------------------------------------------------------------------------------------------------------------
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2) THEN
      IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%SpecToPolyArray
        IF(PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean.GT.0) THEN
          Xi_vib1 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
        ELSE
          Xi_vib1 = 0.0
          ! The vibrational energy of the dissociating molecule and the char. vib. temps of the product are used to determine a
          ! first guess for the vibrational degree of freedom
          IF(PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib.GT.0) THEN
            DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
              Xi_vib1 = Xi_vib1 + (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) & 
                        / PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib) &
                        / (exp(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
                        / PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib) - 1.0)
            END DO
          END IF 
        END IF
      ELSE
        IF(ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,1)).GT.0) THEN
          Xi_vib1 = 2.0*ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,1)) &
                * LOG(1.0/ ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,1)) + 1.0)
        ELSE
          iPolyatMole = SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray
          IF(PolyatomMolDSMC(iPolyatMole)%TVib.GT.0) THEN
            Xi_vib1 = (2.0*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib / PolyatomMolDSMC(iPolyatMole)%TVib) &
                      / (exp(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib / PolyatomMolDSMC(iPolyatMole)%TVib) - 1.0)
          ELSE
            Xi_vib1 = 0.0
          END IF
        END IF
      END IF
      FakXi = FakXi + 0.5*Xi_vib1
    END IF

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Vibrational degree of freedom of second product
    !-------------------------------------------------------------------------------------------------------------------------------
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2) THEN
      IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%SpecToPolyArray
        IF(CollInf%Coll_SpecPartNum(ChemReac%DefinedReact(iReac,2,2)).NE.0) THEN
          Xi_vib2 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
        ELSE
          Xi_vib2 = 0.0
          IF(PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib.GT.0) THEN
            DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
              Xi_vib2 = Xi_vib2 + (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) & 
                        / PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib) &
                        / (exp(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
                        / PolyatomMolDSMC(SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray)%TVib) - 1.0)
            END DO
          END IF
        END IF
      ELSE
        IF(ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,2)).GT.0) THEN
          Xi_vib2 = 2.0*ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,2)) &
                * LOG(1.0/ ChemReac%MeanEVibQua_PerIter(ChemReac%DefinedReact(iReac,2,2)) + 1.0)
        ELSE
          iPolyatMole = SpecDSMC(PartSpecies(React1Inx))%SpecToPolyArray
          IF(PolyatomMolDSMC(iPolyatMole)%TVib.GT.0) THEN
            Xi_vib2 = (2.0*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib / PolyatomMolDSMC(iPolyatMole)%TVib) &
                      / (exp(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib / PolyatomMolDSMC(iPolyatMole)%TVib) - 1.0)
          ELSE
            Xi_vib2 = 0.0
          END IF
        END IF
      END IF
      FakXi = FakXi + 0.5*Xi_vib2
    END IF
  ELSE
    FakXi = 0.5*(Xi_rel + SpecDSMC(PartSpecies(React2Inx))%Xi_Rot) - 1.0 
  END IF

   ! check if electronic model is used
  IF ( DSMC%ElectronicState ) THEN
    ! add electronic energy to collision energy
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(React1Inx,3) + &
                                                  PartStateIntEn(React2Inx,3)
    IF (SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - &
            DSMC%GammaQuant *BoltzmannConst*SpecDSMC(PartSpecies(React2Inx))%CharaTVib 
    END IF
    IF (SpecDSMC(PartSpecies(React2Inx))%InterID.NE.4) THEN
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React2Inx,FakXi,React1Inx,PEM%Element(React1Inx))
      ! store the electronic energy of the dissociating molecule
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi,React2Inx,PEM%Element(React2Inx))
      ElecTransfer = PartStateIntEn(React1Inx,3)
    END IF
    IF (SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + &
            DSMC%GammaQuant *BoltzmannConst*SpecDSMC(PartSpecies(React2Inx))%CharaTVib 
    END IF
  END IF
  ! Zero-point energy of the products of dissociation is substracted from total collision energy before the relaxation
  ! of non-reacting partner and added again after the relaxation
  IF (SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
      EZeroTempToExec1 = PolyatomMolDSMC(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%SpecToPolyArray)%EZeroPoint
    ELSE
      EZeroTempToExec1 = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
    END IF
  END IF
  IF (SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
      EZeroTempToExec2 = PolyatomMolDSMC(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%SpecToPolyArray)%EZeroPoint
    ELSE
      EZeroTempToExec2 = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib
    END IF
  END IF
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - EZeroTempToExec1 - EZeroTempToExec2
  ! Relaxation of non-reacting collision partner
  IF (SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2) THEN
    IF(SpecDSMC(PartSpecies(React2Inx))%PolyatomicMol) THEN
!      CALL DSMC_VibRelaxPoly(Coll_pData(iPair)%Ec,PartSpecies(React2Inx),React2Inx,FakXi)
      CALL DSMC_InsertPolyProduct(PartSpecies(React2Inx),PolyatomMolDSMC(SpecDSMC(PartSpecies(React2Inx))%SpecToPolyArray)%TVib, &
                                                                                                                  React2Inx,iPair)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,1)
    ELSE
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(React2Inx))%CharaTVib)  &
                - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(React2Inx))%MaxVibQuant)
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
        IF (PartMPF(React2Inx).GT.ReacMPF) THEN
        !Vibrational Relaxation of React2Inx
          DeltaPartStateIntEn = 0.0
          Phi = ReacMPF / PartMPF(React2Inx)
          PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(React2Inx))%CharaTVib
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
          PartStateIntEnTemp = (1-Phi) * PartStateIntEn(React2Inx,1) + Phi * PartStateIntEnTemp
          ! searche for new vib quant
          MaxColQua = PartStateIntEnTemp/(BoltzmannConst*SpecDSMC(PartSpecies(React2Inx))%CharaTVib)  &
                    - DSMC%GammaQuant
          iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(React2Inx))%MaxVibQuant)
          CALL RANDOM_NUMBER(iRan)
          iQua = INT(iRan * iQuaMax)
          CALL RANDOM_NUMBER(iRan)
          DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
           !laux diss page 31
           CALL RANDOM_NUMBER(iRan)
           iQua = INT(iRan * iQuaMax)
           CALL RANDOM_NUMBER(iRan)
          END DO
          PartStateIntEn(React2Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(React2Inx))%CharaTVib
          DeltaPartStateIntEn = PartMPF(React2Inx) &
                              * (PartStateIntEnTemp - PartStateIntEn(React2Inx,1))
          ! adding in-energy lost due to vMPF
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + DeltaPartStateIntEn / ReacMPF
        END IF
      ELSE ! no vMPF or MPF of React2Inx .eq. ReacMPF
      !Vibrational Relaxation of React2Inx
        PartStateIntEn(React2Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(PartSpecies(React2Inx))%CharaTVib
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,1)
      END IF ! (usevMPF).AND.(PartMPF(React2Inx).GT.ReacMPF)
    END IF
  END IF
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + EZeroTempToExec1
  ! Relaxation of first dissociation product
  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2) THEN
    FakXi = FakXi - 0.5*Xi_vib1
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib1)*Xi_vib1)/(Xi_vib1*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,1),TVibTemp,React1Inx,iPair)
    ELSE
      IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) DEALLOCATE(VibQuantsPar(React1Inx)%Quants)
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib)  &
      - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(React1Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
      * SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
      DO WHILE ((iRan.GT.(1 - iQua/MaxColQua)**FakXi).OR.((Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1)).LT.0.0)) 
        !laux diss page 31
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEn(React1Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
        * SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
      END DO
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1)
  END IF

  !.... Get free particle index for the 3rd particle produced
  DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
  PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
  IF (PositionNbr.EQ.0) THEN
    CALL abort(&
        __STAMP__,&
    'New Particle Number greater max Part Num in MolecDissoc. Reaction: ',iReac)
  END IF  

  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + EZeroTempToExec2
  !Set new Species of new particle
  PDM%ParticleInside(PositionNbr) = .true.
  PartSpecies(PositionNbr) = ChemReac%DefinedReact(iReac,2,2)
  PartState(PositionNbr,1:3) = PartState(React1Inx,1:3)
  IF(DoRefMapping)THEN ! here Nearst-GP is missing
    PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,React1Inx)
  END IF
  PartStateIntEn(PositionNbr, 1) = 0
  PartStateIntEn(PositionNbr, 2) = 0
  PEM%Element(PositionNbr) = PEM%Element(React1Inx)
  IF(usevMPF) PartMPF(PositionNbr) = ReacMPF
  ! Relaxation of new particle (second dissociation product)
  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2) THEN
    FakXi = FakXi - 0.5*Xi_vib2
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib2)*Xi_vib2)/(Xi_vib2*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,2),TVibTemp,PositionNbr,iPair)
    ELSE
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib)  &
      - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(PositionNbr,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
      * SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib
      DO WHILE ((iRan.GT.(1 - iQua/MaxColQua)**FakXi).OR.((Coll_pData(iPair)%Ec - PartStateIntEn(PositionNbr,1)).LT.0.0)) 
        !laux diss page 31
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)    
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEn(PositionNbr,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
        * SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib
      END DO
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(PositionNbr,1)
  END IF

  IF (SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2) THEN
      IF(SpecDSMC(PartSpecies(React2Inx))%Xi_Rot.EQ.3) THEN
        FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(React2Inx))%Xi_Rot
        CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React2Inx, FakXi)
      ELSE
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEn(React2Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
        FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(React2Inx))%Xi_Rot
      END IF    
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,2)
  END IF

  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React1Inx, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
    END IF
  END IF

  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, PositionNbr, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(PositionNbr,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(PositionNbr,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(PositionNbr,2)
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot
    END IF
  END IF

  ! distribute Etra of pseudo neutral particle (i+e) and old electron
  CALL RANDOM_NUMBER(iRan)
  FacEtraDistri = iRan
  CALL RANDOM_NUMBER(iRan)
  ! laux diss page 40, omegaVHS only of one species
  DO WHILE ((4 *FacEtraDistri*(1-FacEtraDistri))**(1-SpecDSMC(PartSpecies(React1Inx))%omegaVHS).LT.iRan)
    CALL RANDOM_NUMBER(iRan)
    FacEtraDistri = iRan
    CALL RANDOM_NUMBER(iRan)
  END DO
  ERel_React1_React2 = Coll_pData(iPair)%Ec * FacEtraDistri
  ERel_React1_React3 = Coll_pData(iPair)%Ec - ERel_React1_React2

  !Scattering of pseudo molecule (a-b) and collision partner e (scattering of e)
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(React2Inx), Coll_pData(iPair)%PairType)

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(React2Inx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(React2Inx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(React2Inx, 6)

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * ERel_React1_React2 /CollInf%MassRed(Coll_pData(iPair)%PairType)
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)
  
 ! deltaV particle 2
  DSMC_RHS(React2Inx,1) = VeloMx - FracMassCent1*RanVelox - PartState(React2Inx, 4)
  DSMC_RHS(React2Inx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(React2Inx, 5)
  DSMC_RHS(React2Inx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(React2Inx, 6)  

!Debug_Energy(2)  = Debug_Energy(2) &
!      +0.5* Species(PartSpecies(React2Inx))%MassIC&
!                                            *((VeloMx - FracMassCent1*RanVelox)**2   &
!                                             +(VeloMy - FracMassCent1*RanVeloy)**2   &
!                                             +(VeloMz - FracMassCent1*RanVeloz)**2 )


  IF (usevMPF) THEN
    IF (PartMPF(React2Inx).GT.ReacMPF) THEN
      Phi = ReacMPF / PartMPF(React2Inx)
      GEO%DeltaEvMPF(PEM%Element(React2Inx)) = GEO%DeltaEvMPF(PEM%Element(React2Inx)) + 0.5 * PartMPF(React2Inx) &
                                             * Species(PartSpecies(React2Inx))%MassIC &
                                             * Phi * (1 - Phi) &
                                             * ( DSMC_RHS(React2Inx,1)**2 &
                                               + DSMC_RHS(React2Inx,2)**2 &
                                               + DSMC_RHS(React2Inx,3)**2 ) 
      DSMC_RHS(React2Inx,1) = Phi * DSMC_RHS(React2Inx,1)
      DSMC_RHS(React2Inx,2) = Phi * DSMC_RHS(React2Inx,2)
      DSMC_RHS(React2Inx,3) = Phi * DSMC_RHS(React2Inx,3) 
    END IF
  END IF

  !Set velocity of pseudo molec (a+b) and calculate the centre of mass frame velocity: m_pseu / (m_3 + m_4) * v_pseu
  !(Velocity of pseudo molecule is NOT equal to the COM frame velocity)
  VxPseuMolec = (VeloMx + FracMassCent2*RanVelox) * Species(PartSpecies(React1Inx))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)
  VyPseuMolec = (VeloMy + FracMassCent2*RanVeloy) * Species(PartSpecies(React1Inx))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)
  VzPseuMolec = (VeloMz + FracMassCent2*RanVeloz) * Species(PartSpecies(React1Inx))%MassIC &
                  / (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(ChemReac%DefinedReact(iReac,2,2))%MassIC)

  !Set new Species of dissoc atom
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)

  ! here set the electronic level of the 2 atoms
  IF ( DSMC%ElectronicState ) THEN
    Coll_pData(iPair)%Ec = ElecTransfer + ERel_React1_React3
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec, React1Inx   , FakXi, PositionNbr,PEM%Element(React1Inx) )
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec, PositionNbr , FakXi, React1Inx,PEM%Element(PositionNbr) )
    ERel_React1_React3 = Coll_pData(iPair)%Ec
  END IF

  !Scattering of a + b
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), &
                  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(PositionNbr), & 
                  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 *  ERel_React1_React3 / & 
          CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  !deltaV particle 1
  DSMC_RHS(React1Inx,1) = VxPseuMolec + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = VyPseuMolec + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = VzPseuMolec + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

  IF(SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.1) THEN
    PartStateIntEn(React1Inx, 1) = 0
    PartStateIntEn(React1Inx, 2) = 0
  END IF
  !deltaV new formed particle
  PartState(PositionNbr,4:6) = 0
  DSMC_RHS(PositionNbr,1) = VxPseuMolec - FracMassCent1*RanVelox 
  DSMC_RHS(PositionNbr,2) = VyPseuMolec - FracMassCent1*RanVeloy 
  DSMC_RHS(PositionNbr,3) = VzPseuMolec - FracMassCent1*RanVeloz 

  IF(usevMPF) THEN
    IF (ReacMPF.GT.(Species(PartSpecies(React1Inx))%MacroParticleFactor)) THEN
      CALL vMPF_AfterSplitting(React1Inx, ReacMPF, Species(PartSpecies(React1Inx))%MacroParticleFactor)
    END IF
    IF (ReacMPF.GT.(Species(PartSpecies(PositionNbr))%MacroParticleFactor)) THEN
      CALL vMPF_AfterSplitting(PositionNbr, ReacMPF, Species(PartSpecies(PositionNbr))%MacroParticleFactor)
    END IF
  END IF

!Debug_Energy(2)  = Debug_Energy(2) &
!      +0.5* Species(PartSpecies(React1Inx))%MassIC&
!                                            *((VxPseuMolec + FracMassCent2*RanVelox )**2   &
!                                             +(VyPseuMolec + FracMassCent2*RanVeloy )**2   &
!                                             +(VzPseuMolec + FracMassCent2*RanVeloz )**2 ) &
!      +0.5* Species(PartSpecies(PositionNbr))%MassIC&
!                                            *((VxPseuMolec - FracMassCent1*RanVelox )**2   &
!                                             +(VyPseuMolec - FracMassCent1*RanVeloy )**2   &
!                                             +(VzPseuMolec - FracMassCent1*RanVeloz )**2 ) &
!      + PartStateIntEn(React2Inx,1) &
!      + PartStateIntEn(React2Inx,2)

!print*,Debug_Energy(1), Debug_Energy(2), "    Diff=",Debug_Energy(2)-Debug_Energy(1)
!read*
END SUBROUTINE MolecDissoc


SUBROUTINE MolecExch(iReac, iPair)
!===================================================================================================================================
! Perform molecular exchange reaction
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY : abort
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, CollisMode, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar, Debug_Energy
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec
USE MOD_vmpf_collision,        ONLY : vMPF_AfterSplitting
USE MOD_Particle_Vars,         ONLY : usevMPF, PartMPF, RandomVec, Species, PartPosRef
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange
USE MOD_DSMC_PolyAtomicModel,  ONLY : DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, FakXiPoly, DSMC_InsertPolyProduct
USE MOD_DSMC_Analyze,          ONLY : CalcTVib, CalcTVibPoly
USE MOD_Particle_Tracking_Vars,ONLY : DoRefmapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec, iPolyatMole
  REAL                          :: JToEv, FakXi, Xi_rel, iRan
  INTEGER                       :: iQuaMax, iQua, React1Inx, React2Inx, NonReacPart
  REAL                          :: MaxColQua
  REAL                          :: ReacMPF
  REAL                          :: Xi_vib1, Xi_vib2
  REAL                          :: EZeroPoint, TVibTemp
!===================================================================================================================================

  JToEv = 1.602176565E-19
  Xi_vib1 = 0.
  Xi_vib2 = 0.
  EZeroPoint = 0.0

!..Get the index of react1 and the react2
  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2 
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF

!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) &
    !+ 0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
    !+ 0.5* Species(PartSpecies(React2Inx))%MassIC*(PartState(React2Inx,4)**2+PartState(React2Inx,5)**2+PartState(React2Inx,6)**2)&
      !+ PartStateIntEn(React1Inx,1)+ PartStateIntEn(React2Inx,1) &
      !+ PartStateIntEn(React1Inx,2)+ PartStateIntEn(React2Inx,2)

  IF (usevMPF) THEN ! reaction MPF definition
    ReacMPF = MIN(PartMPF(React1Inx), PartMPF(React2Inx))
    IF (PartMPF(React1Inx).GT.ReacMPF) THEN ! just a part of the molecule 1 react
    !.... Get free particle index for the non-reacting particle part
      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
      NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
      IF (NonReacPart.EQ.0) THEN
        CALL abort(&
            __STAMP__,&
        'New Particle Number greater max Part Num in MolecExchange. Reaction: ',iReac)
      END IF
    ! Copy molecule data for non-reacting particle part
      PDM%ParticleInside(NonReacPart) = .true.
      PartSpecies(NonReacPart)        = PartSpecies(React1Inx)
      PartState(NonReacPart,1:6)      = PartState(React1Inx,1:6)
      IF(DoRefMapping)THEN ! here Nearst-GP is missing
        PartPosRef(1:3,NonReacPart)=PartPosRef(1:3,React1Inx)
      END IF
      PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React1Inx, 1)
      PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React1Inx, 2)
      IF (DSMC%ElectronicState) THEN
        PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React1Inx, 3)
      END IF
      PEM%Element(NonReacPart)        = PEM%Element(React1Inx)
      PartMPF(NonReacPart)            = PartMPF(React1Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
      PartMPF(React1Inx)              = ReacMPF ! reacting part MPF = ReacMPF
    ELSE IF (PartMPF(React2Inx).GT.ReacMPF) THEN ! just a part of the molecule 2 react
    !.... Get free particle index for the non-reacting particle part
      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
      NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
      IF (NonReacPart.EQ.0) THEN
        CALL abort(&
            __STAMP__,&
        'New Particle Number greater max Part Num in MolecExchange. Reaction: ',iReac)
      END IF
    ! Copy molecule data for non-reacting particle part
      PDM%ParticleInside(NonReacPart) = .true.
      PartSpecies(NonReacPart)        = PartSpecies(React2Inx)
      PartState(NonReacPart,1:6)      = PartState(React2Inx,1:6)
      IF(DoRefMapping)THEN ! here Nearst-GP is missing
        PartPosRef(1:3,NonReacPart)=PartPosRef(1:3,React2Inx)
      END IF
      PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React2Inx, 1)
      PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React2Inx, 2)
      IF (DSMC%ElectronicState) THEN
        PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React2Inx, 3)
      END IF
      PEM%Element(NonReacPart)        = PEM%Element(React2Inx)
      PartMPF(NonReacPart)            = PartMPF(React2Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
      PartMPF(React2Inx)              = ReacMPF ! reacting part MPF = ReacMPF
    END IF
  END IF

  ! Add heat of formation to collision energy
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac)

  ! DOF of relative motion in VHS model, only for one omega!!
  ! this is a result of the mean value of the relative energy in the vhs model, laux diss page 31

  Xi_rel = 2.*(2. - SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%omegaVHS)

  FakXi = 0.5*(Xi_rel + SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot + SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot) - 1
  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol.OR.SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) &
                                                                  CALL FakXiPoly(iReac, FakXi, Xi_vib1, Xi_vib2)
  ! Set new Species of molec and atom
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PartSpecies(React2Inx) = ChemReac%DefinedReact(iReac,2,2)

  ! check if electronic model is used
  IF ( DSMC%ElectronicState ) THEN
    ! add electronic energy to collision energy
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,3) + &
                                                  PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi,React2Inx,PEM%Element(React1Inx) )
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React2Inx,FakXi,React1Inx,PEM%Element(React2Inx) )
  END IF

  ! Vibrational Relaxation of React1Inx
  IF((SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.20)) THEN
    FakXi = FakXi - 0.5*Xi_vib1
    IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.20)) THEN
      IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%SpecToPolyArray
        EZeroPoint = PolyatomMolDSMC(iPolyatMole)%EZeroPoint
      ELSE
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%CharaTVib
      END IF
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - EZeroPoint
    END IF
    IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib1)*Xi_vib1)/(Xi_vib1*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,1),TVibTemp,React1Inx,iPair)
    ELSE
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
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1) + EZeroPoint
  END IF
  ! Vibrational Relaxation of React2Inx
  IF((SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(React2Inx))%InterID.EQ.20)) THEN
    FakXi = FakXi - 0.5*Xi_Vib2
    IF(SpecDSMC(PartSpecies(React2Inx))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib2)*Xi_vib2)/(Xi_vib2*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,2),TVibTemp,React2Inx,iPair)
    ELSE
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(React2Inx))%CharaTVib)  &
                - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(React2Inx))%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
      DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi) 
        !laux diss page 31
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)    
        CALL RANDOM_NUMBER(iRan)
      END DO
      PartStateIntEn(React2Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                    * SpecDSMC(PartSpecies(React2Inx))%CharaTVib 
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,1)
  END IF
  ! Rotational Relaxation 1
  IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React1Inx, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
    END IF
  ELSE
    ! If new particle is an atom, internal energies are zero
    PartStateIntEn(React1Inx,1) = 0.0
    PartStateIntEn(React1Inx,2) = 0.0
  END IF
  ! Rotational Relaxation 2
  IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React2Inx, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(React2Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React2Inx,2)
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,2))%Xi_Rot
    END IF
  ELSE
    ! If new particle is an atom, internal energies are zero
    PartStateIntEn(React2Inx,1) = 0.0
    PartStateIntEn(React2Inx,2) = 0.0
  END IF
  
!--------------------------------------------------------------------------------------------------! 
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------! 
  FracMassCent1 = CollInf%FracMassCent(ChemReac%DefinedReact(iReac,1,1), &
                CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))
  FracMassCent2 = CollInf%FracMassCent(ChemReac%DefinedReact(iReac,1,2), & 
                CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))

  !Calculation of velo from center of mass from old particle pair
  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(React2Inx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(React2Inx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(React2Inx, 6)

  ! FracMassCent with the masses of products for calc of CRela2 and velo distribution
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(React2Inx), CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * Coll_pData(iPair)%Ec/ &
            CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  !deltaV particle 1  
  DSMC_RHS(React1Inx,1) = VeloMx + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = VeloMy + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = VeloMz + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

  ! deltaV particle 2
  DSMC_RHS(React2Inx,1) = VeloMx - FracMassCent1*RanVelox - PartState(React2Inx, 4)
  DSMC_RHS(React2Inx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(React2Inx, 5)
  DSMC_RHS(React2Inx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(React2Inx, 6)

  IF(usevMPF) THEN
    IF (ReacMPF.GT.(Species(PartSpecies(React1Inx))%MacroParticleFactor)) THEN
      CALL vMPF_AfterSplitting(React1Inx, ReacMPF, Species(PartSpecies(React1Inx))%MacroParticleFactor)
    END IF
  END IF
  IF(usevMPF) THEN
    IF (ReacMPF.GT.(Species(PartSpecies(React2Inx))%MacroParticleFactor)) THEN
      CALL vMPF_AfterSplitting(React2Inx, ReacMPF, Species(PartSpecies(React2Inx))%MacroParticleFactor)
    END IF
  END IF


!Debug_Energy(2)  = Debug_Energy(2) &
      !+0.5* Species(PartSpecies(React1Inx))%MassIC&
                                            !*((VeloMx + FracMassCent2*RanVelox )**2   &
                                             !+(VeloMx + FracMassCent2*RanVeloy )**2   &
                                             !+(VeloMx + FracMassCent2*RanVeloz )**2 ) &
      !+0.5* Species(PartSpecies(React2Inx))%MassIC&
                                            !*((VeloMx - FracMassCent1*RanVelox )**2   &
                                             !+(VeloMx - FracMassCent1*RanVeloy )**2   &
                                             !+(VeloMx - FracMassCent1*RanVeloz )**2 ) &
      !+ PartStateIntEn(React1Inx,1)+ PartStateIntEn(React1Inx,2) &
      !+ PartStateIntEn(React2Inx,1)+ PartStateIntEn(React2Inx,2) ! is zero if atom
!print*,Debug_Energy(1), Debug_Energy(2), "    Diff=",Debug_Energy(2)-Debug_Energy(1)
!print*,"press enter"
!read*
END SUBROUTINE MolecExch


SUBROUTINE AtomRecomb(iReac, iPair, iPart_p3)
!===================================================================================================================================
! Performs recombination between two atoms to one molecule
! atom recombination routine           A + B + X -> AB + X
!===================================================================================================================================
! MODULES
  USE MOD_Globals,               ONLY : abort
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles
  USE MOD_DSMC_Vars,             ONLY : ChemReac, CollisMode, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar, Debug_Energy
  USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec, Species
  USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange
  USE MOD_DSMC_PolyAtomicModel,  ONLY : DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, DSMC_InsertPolyProduct
  USE MOD_DSMC_Analyze,          ONLY : CalcTVib, CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac, iPart_p3
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec
  REAL                          :: FakXi, Xi, iRan
  INTEGER                       :: iQuaMax, iQua, React1Inx, React2Inx
  REAL                          :: MaxColQua
! additional for Q-K theory
  REAL                          :: ksum, Tcoll
  INTEGER                       :: ii
  REAL                          :: TVibTemp, Xi_vib1
  INTEGER                       :: iDOF,iPolyatMole
!===================================================================================================================================

  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF

!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) +&
!      0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
!    + 0.5* Species(PartSpecies(React2Inx))%MassIC*(PartState(React2Inx,4)**2+PartState(React2Inx,5)**2+PartState(React2Inx,6)**2)&
!    + 0.5* Species(PartSpecies(iPart_p3))%MassIC* (PartState(iPart_p3, 4)**2+PartState(iPart_p3, 5)**2+PartState(iPart_p3, 6)**2)

  ! Calculation of the centre of mass of the product molecule
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(React2Inx), Coll_pData(iPair)%PairType)

  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(React2Inx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(React2Inx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(React2Inx, 6)  

  ! The input particle 1 is replaced by the product molecule, the
  !     second input particle is deleted
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PDM%ParticleInside(React2Inx) = .FALSE.

  PartState(React1Inx, 4) = VeloMx
  PartState(React1Inx, 5) = VeloMy
  PartState(React1Inx, 6) = VeloMz

  ! has to be calculated earlier because of setting of electronic energy
  Xi = 2.0 * (2.0 - SpecDSMC(PartSpecies(iPart_p3))%omegaVHS) + SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot &
      + SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%SpecToPolyArray
    IF(PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean.GT.0.0) THEN  ! does the species already exist in the cell? 
      Xi_vib1 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
    ELSE                                                      ! if not, use the calculated vib temp
      Xi_vib1 = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Xi_vib1 = Xi_vib1 + (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)    &
                            / CalcTVibPoly((ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,1)) &
                            + ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,2))),ChemReac%DefinedReact(iReac,2,1))) &
                            / (exp(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)  &
                            / CalcTVibPoly((ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,1)) &
                            + ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,2))),ChemReac%DefinedReact(iReac,2,1))) - 1)
      END DO
    END IF
  END IF

  ! check if atomic electron shell is modelled
  IF ( DSMC%ElectronicState ) THEN
  ! Add heat of formation to collision energy
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac) - PartStateIntEn(iPart_p3,1) + &
                           PartStateIntEn(Coll_pData(iPair)%iPart_p1,3) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)
!--------------------------------------------------------------------------------------------------!
! electronic relaxation  of AB and X (if X is not an electron) 
!--------------------------------------------------------------------------------------------------!
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi )
    Coll_pData(iPair)%Ec =  Coll_pData(iPair)%Ec + PartStateIntEn(iPart_p3,3)
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,iPart_p3,FakXi )
  ELSE
    ! Collision energy is the sum of the relative translational energy of the recombining educts (first term) 
    ! and the product and third non-reacting particle (second term), heat of formation and internal energy of both educts and
    ! rotational energy of the third non-reacting particle
    Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                       &
                 + 0.5 * CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx), PartSpecies(iPart_p3)))                  &
                 * ((VeloMx-PartState(iPart_p3,4))**2+(VeloMy-PartState(iPart_p3,5))**2+(VeloMz-PartState(iPart_p3,6))**2)  &
                 + ChemReac%EForm(iReac) + PartStateIntEn(React1Inx,1) + PartStateIntEn(React2Inx,1)                        &
                 + PartStateIntEn(React1Inx,2) + PartStateIntEn(React2Inx,2) + PartStateIntEn(iPart_p3,2)
  END IF

!--------------------------------------------------------------------------------------------------! 
! Vibrational Relaxation of AB and X (if X is a molecule) 
!--------------------------------------------------------------------------------------------------!
  ! chose between Q-K and Arrhenius for detailed balancing
  IF (ChemReac%QKProcedure(iReac)) THEN ! Q-K
    MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(React1Inx))%CharaTVib)  &
              - DSMC%GammaQuant
    iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(React1Inx))%MaxVibQuant)
    ksum = 0.
    Tcoll = CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2  &
          / ( 2 * BoltzmannConst * ( 2 - SpecDSMC(PartSpecies(React1Inx))%omegaVHS ) )
    DO ii = 0, iQuaMax -1
      ksum = ksum + gammainc( [2. - SpecDSMC(PartSpecies(React1Inx))%omegaVHS,                          &
                               (iQuaMax - ii)*SpecDSMC(PartSpecies(React1Inx))%CharaTVib / Tcoll ] ) * &
             exp( - ii * SpecDSMC(PartSpecies(React1Inx))%CharaTVib / Tcoll )
    END DO
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT. ( gammainc([2. - SpecDSMC(PartSpecies(React1Inx))%omegaVHS,                          &
                               (iQuaMax - iQua)*SpecDSMC(PartSpecies(React1Inx))%CharaTVib / Tcoll ] ) *   &
                          exp( - iQua * SpecDSMC(PartSpecies(React1Inx))%CharaTVib/ Tcoll ) / ksum ) )
      ! diss page 31
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)    
      CALL RANDOM_NUMBER(iRan)
    END DO
    PartStateIntEn(React1Inx,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                  * SpecDSMC(PartSpecies(React1Inx))%CharaTVib
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1)+ PartStateIntEn(iPart_p3,1)
  ELSE
    IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib1)*Xi_vib1)/(Xi_vib1*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,1),TVibTemp,React1Inx,iPair)
    ELSE
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
    END IF

    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1) + PartStateIntEn(iPart_p3,1)

    ! X particle
    IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ. 2) THEN
      IF(SpecDSMC(PartSpecies(iPart_p3))%PolyatomicMol) THEN
        CALL DSMC_VibRelaxPoly(Coll_pData(iPair)%Ec,PartSpecies(iPart_p3),iPart_p3,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,1)
      ELSE
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
        PartStateIntEn(iPart_p3,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(PartSpecies(iPart_p3))%CharaTVib 
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,1) 
      END IF
    END IF
  END IF
!--------------------------------------------------------------------------------------------------! 
! rotational Relaxation of AB and X (if X is a molecule) 
!--------------------------------------------------------------------------------------------------!
  IF(SpecDSMC(PartSpecies(React1Inx))%Xi_Rot.EQ.3) THEN
    FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
    CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React1Inx, FakXi)      
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
  ELSE  
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
  END IF

  IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2) THEN
    IF(SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(iPart_p3))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, iPart_p3, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart_p3,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
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
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  ! deltaV particle 1
  DSMC_RHS(React1Inx,1) = VeloMx + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = VeloMy + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = VeloMz + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

 ! deltaV particle 2
  DSMC_RHS(iPart_p3,1) = VeloMx - FracMassCent1*RanVelox - PartState(iPart_p3, 4)
  DSMC_RHS(iPart_p3,2) = VeloMy - FracMassCent1*RanVeloy - PartState(iPart_p3, 5)
  DSMC_RHS(iPart_p3,3) = VeloMz - FracMassCent1*RanVeloz - PartState(iPart_p3, 6)

!Debug_Energy(2) = Debug_Energy(2)&
!      + 0.5* Species(PartSpecies(React1Inx))%MassIC * (&
!           (VeloMx + FracMassCent2*RanVelox)**2 &
!          +(VeloMy + FracMassCent2*RanVeloy)**2 &
!          +(VeloMz + FracMassCent2*RanVeloz)**2) &
!      + 0.5* Species(PartSpecies(iPart_p3))%MassIC  * (&
!           (VeloMx - FracMassCent1*RanVelox)**2 &
!          +(VeloMy - FracMassCent1*RanVeloy)**2 &
!          +(VeloMz - FracMassCent1*RanVeloz)**2) &
!          +PartStateIntEn(React1Inx,1)&
!          +PartStateIntEn(React1Inx,2)&
!          -ChemReac%EForm(iReac)
!
!print*, Debug_Energy(1),Debug_Energy(2), Debug_Energy(2)-Debug_Energy(1)
!read*
END SUBROUTINE AtomRecomb


SUBROUTINE simpleCEX(iReac, iPair)
!===================================================================================================================================
! simple charge exchange interaction     
! ION(v1) + ATOM(v2) -> ATOM(v1) + ION(v2)
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS
  USE MOD_DSMC_Vars,             ONLY : ChemReac
  USE MOD_Particle_Vars,         ONLY : PartSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: React1Inx, React2Inx
!===================================================================================================================================

  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF
  ! change species
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PartSpecies(React2Inx) = ChemReac%DefinedReact(iReac,2,2)
  ! deltaV particle 1
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,1) = 0.
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,2) = 0.
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,3) = 0.
  ! deltaV particle 2
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,1) = 0.
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,2) = 0.
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,3) = 0.

END SUBROUTINE simpleCEX

SUBROUTINE AssoIonization(iReac, iPair)
!===================================================================================================================================
! Perform molecular exchange reaction
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY : abort
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, DSMC, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, CollisMode, PartStateIntEn, PolyatomMolDSMC, VibQuantsPar
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec
USE MOD_vmpf_collision,        ONLY : vMPF_AfterSplitting
USE MOD_Particle_Vars,         ONLY : usevMPF, PartMPF, RandomVec, Species
USE MOD_Particle_Mesh_Vars,    ONLY : GEO
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange
USE MOD_DSMC_PolyAtomicModel,  ONLY : DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, FakXiPoly, DSMC_InsertPolyProduct
USE MOD_DSMC_Analyze,          ONLY : CalcTVib, CalcTVibPoly
USE MOD_Particle_Tracking_Vars,ONLY : DoRefmapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec, iPolyatMole, iDOF
  REAL                          :: JToEv, FakXi, Xi_rel, iRan
  INTEGER                       :: iQuaMax, iQua, React1Inx, React2Inx
  REAL                          :: MaxColQua
  REAL                          :: Xi_vib1
  REAL                          :: EZeroPoint, TVibTemp
!--- Variables for vMPF
!  REAL                          :: ReacMPF
!  INTEGER                       :: NonReacPart
!===================================================================================================================================

JToEv = 1.602176565E-19
Xi_vib1 = 0.
EZeroPoint = 0.0

!..Get the index of react1 and the react2
  IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%DefinedReact(iReac,1,1)) THEN
    React1Inx = Coll_pData(iPair)%iPart_p1
    React2Inx = Coll_pData(iPair)%iPart_p2 
  ELSE
    React2Inx = Coll_pData(iPair)%iPart_p1
    React1Inx = Coll_pData(iPair)%iPart_p2
  END IF

!  IF (usevMPF) THEN ! reaction MPF definition
!    ReacMPF = MIN(PartMPF(React1Inx), PartMPF(React2Inx))
!    IF (PartMPF(React1Inx).GT.ReacMPF) THEN ! just a part of the molecule 1 react
!    !.... Get free particle index for the non-reacting particle part
!      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
!      NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
!      IF (NonReacPart.EQ.0) THEN
!        CALL abort(__STAMP__,&
!        'New Particle Number greater max Part Num in MolecExchange. Reaction: ',iReac)
!      END IF
!    ! Copy molecule data for non-reacting particle part
!      PDM%ParticleInside(NonReacPart) = .true.
!      PartSpecies(NonReacPart)        = PartSpecies(React1Inx)
!      PartState(NonReacPart,1:6)      = PartState(React1Inx,1:6)
!      IF(DoRefMapping)THEN ! here Nearst-GP is missing
!        PartPosRef(1:3,NonReacPart)=PartPosRef(1:3,React1Inx)
!      END IF
!      PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React1Inx, 1)
!      PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React1Inx, 2)
!      IF (DSMC%ElectronicState) THEN
!        PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React1Inx, 3)
!      END IF
!      PEM%Element(NonReacPart)        = PEM%Element(React1Inx)
!      PartMPF(NonReacPart)            = PartMPF(React1Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
!      PartMPF(React1Inx)              = ReacMPF ! reacting part MPF = ReacMPF
!    ELSE IF (PartMPF(React2Inx).GT.ReacMPF) THEN ! just a part of the molecule 2 react
!    !.... Get free particle index for the non-reacting particle part
!      DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
!      NonReacPart = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
!      IF (NonReacPart.EQ.0) THEN
!        CALL abort(__STAMP__,&
!        'New Particle Number greater max Part Num in MolecExchange. Reaction: ',iReac)
!      END IF
!    ! Copy molecule data for non-reacting particle part
!      PDM%ParticleInside(NonReacPart) = .true.
!      PartSpecies(NonReacPart)        = PartSpecies(React2Inx)
!      PartState(NonReacPart,1:6)      = PartState(React2Inx,1:6)
!      IF(DoRefMapping)THEN ! here Nearst-GP is missing
!        PartPosRef(1:3,NonReacPart)=PartPosRef(1:3,React2Inx)
!      END IF
!      PartStateIntEn(NonReacPart, 1)  = PartStateIntEn(React2Inx, 1)
!      PartStateIntEn(NonReacPart, 2)  = PartStateIntEn(React2Inx, 2)
!      IF (DSMC%ElectronicState) THEN
!        PartStateIntEn(NonReacPart, 3)  = PartStateIntEn(React2Inx, 3)
!      END IF
!      PEM%Element(NonReacPart)        = PEM%Element(React2Inx)
!      PartMPF(NonReacPart)            = PartMPF(React2Inx) - ReacMPF ! MPF of non-reacting particle part = MPF Diff
!      PartMPF(React2Inx)              = ReacMPF ! reacting part MPF = ReacMPF
!    END IF
!  END IF

  ! Add heat of formation to collision energy
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + ChemReac%EForm(iReac)

  ! DOF of relative motion in VHS model, only for one omega!!
  ! this is a result of the mean value of the relative energy in the vhs model, laux diss page 31

  Xi_rel = 2.*(2. - SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%omegaVHS)

  FakXi = 0.5*(Xi_rel + SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot) - 1
  
  IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%SpecToPolyArray
    IF(PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean.GT.0.0) THEN  ! does the species already exist in the cell? 
      Xi_vib1 = PolyatomMolDSMC(iPolyatMole)%Xi_Vib_Mean
    ELSE                                                      ! if not, use the calculated vib temp
      Xi_vib1 = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Xi_vib1 = Xi_vib1 + (2.0*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)    &
                            / CalcTVibPoly((ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,1)) &
                            + ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,2))),ChemReac%DefinedReact(iReac,2,1))) &
                            / (exp(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)  &
                            / CalcTVibPoly((ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,1)) &
                            + ChemReac%MeanEVib_PerIter(ChemReac%DefinedReact(iReac,1,2))),ChemReac%DefinedReact(iReac,2,1))) - 1)
      END DO
    END IF
  END IF

  ! Set new Species of molec and atom
  PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
  PartSpecies(React2Inx) = ChemReac%DefinedReact(iReac,2,2)

!  ! check if electronic model is used
!  IF ( DSMC%ElectronicState ) THEN
!    ! add electronic energy to collision energy
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,3) + &
!                                                  PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)
!    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi,React2Inx,PEM%Element(React1Inx) )
!    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React2Inx,FakXi,React1Inx,PEM%Element(React2Inx) )
!  END IF

  ! Vibrational Relaxation of React1Inx
  IF(SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.20) THEN
    FakXi = FakXi - 0.5*Xi_vib1
    IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN
      TVibTemp = 2.0*(Coll_pData(iPair)%Ec/((FakXi+1.0)*2.0+Xi_vib1)*Xi_vib1)/(Xi_vib1*BoltzmannConst)
      CALL DSMC_InsertPolyProduct(ChemReac%DefinedReact(iReac,2,1),TVibTemp,React1Inx,iPair)
    ELSE
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
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1) + EZeroPoint
  END IF

  ! Rotational Relaxation 1
  IF((SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.2).OR.(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot.EQ.3) THEN
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, React1Inx, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2)
      FakXi = FakXi - 0.5*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
    END IF
  ELSE
    ! If new particle is an atom, internal energies are zero
    PartStateIntEn(React1Inx,1) = 0.0
    PartStateIntEn(React1Inx,2) = 0.0
  END IF

!--------------------------------------------------------------------------------------------------! 
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------! 
  FracMassCent1 = CollInf%FracMassCent(ChemReac%DefinedReact(iReac,1,1), &
                CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))
  FracMassCent2 = CollInf%FracMassCent(ChemReac%DefinedReact(iReac,1,2), & 
                CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))

  !Calculation of velo from center of mass from old particle pair
  VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
         + FracMassCent2 * PartState(React2Inx, 4)
  VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
         + FracMassCent2 * PartState(React2Inx, 5)
  VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
         + FracMassCent2 * PartState(React2Inx, 6)

  ! FracMassCent with the masses of products for calc of CRela2 and velo distribution
  FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(React2Inx), CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * Coll_pData(iPair)%Ec/ &
            CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(React2Inx)))
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  !deltaV particle 1  
  DSMC_RHS(React1Inx,1) = VeloMx + FracMassCent2*RanVelox - PartState(React1Inx, 4)
  DSMC_RHS(React1Inx,2) = VeloMy + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
  DSMC_RHS(React1Inx,3) = VeloMz + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

  ! deltaV particle 2
  DSMC_RHS(React2Inx,1) = VeloMx - FracMassCent1*RanVelox - PartState(React2Inx, 4)
  DSMC_RHS(React2Inx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(React2Inx, 5)
  DSMC_RHS(React2Inx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(React2Inx, 6)

!  IF(usevMPF) THEN
!    IF (ReacMPF.GT.(Species(PartSpecies(React1Inx))%MacroParticleFactor)) THEN
!      CALL vMPF_AfterSplitting(React1Inx, ReacMPF, Species(PartSpecies(React1Inx))%MacroParticleFactor)
!    END IF
!  END IF
!  IF(usevMPF) THEN
!    IF (ReacMPF.GT.(Species(PartSpecies(React2Inx))%MacroParticleFactor)) THEN
!      CALL vMPF_AfterSplitting(React2Inx, ReacMPF, Species(PartSpecies(React2Inx))%MacroParticleFactor)
!    END IF
!  END IF

END SUBROUTINE AssoIonization

SUBROUTINE SetMeanVibQua()
!===================================================================================================================================
! Computes the vibrational quant of species or the mean vibrational energy (polyatomic)
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,             ONLY : DSMC, CollInf, SpecDSMC, ChemReac, BGGas
  USE MOD_Particle_Vars,         ONLY : BoltzmannConst, nSpecies
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

  DO iSpec =1, nSpecies
    ! describe evib as quantum number
    IF (iSpec.EQ.BGGas%BGGasSpecies) THEN
      ChemReac%MeanEVibQua_PerIter(iSpec) = BGGas%BGMeanEVibQua
    ELSE
      IF(CollInf%Coll_SpecPartNum(iSpec).GT.0) THEN
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN 
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
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


SUBROUTINE CalcBackwardRate(iReacTmp,LocalTemp,BackwardRate)
!===================================================================================================================================
! Calculation of the backward reaction rate with partition sums, interpolation within the given temperature interval
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,             ONLY : DSMC, SpecDSMC, ChemReac
  USE MOD_Particle_Vars,         ONLY : nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iReacTmp
  REAL, INTENT(IN)              :: LocalTemp
  REAL, INTENT(OUT)             :: BackwardRate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: iReac, iSpec, LowerLevel, UpperLevel, iChemDir
  REAL                            :: PartFuncProduct(2), k_b_lower, k_b_upper
!===================================================================================================================================

  ! Determination of the lower and upper value of the temperature interval
  LowerLevel = INT(LocalTemp/DSMC%PartitionInterval)
  UpperLevel = LowerLevel + 1

  ! Reading the stoichiometric coefficients from the reactants
  iReac = iReacTmp - ChemReac%NumOfReact/2

  IF(UpperLevel.GT.INT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)) THEN
    CALL abort(&
     __STAMP__&
      ,'Temperature limit for the backward reaction rate calculation exceeds the given value! Temp: ',RealInfoOpt=LocalTemp)
  END IF

  ! Calculation of the backward reaction rate at the lower temperature value (using the equilibrium constant)
  PartFuncProduct(1:2) = 1.
  DO iSpec = 1, nSpecies
    DO iChemDir = 1,2
      IF(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir).NE.0) THEN
        PartFuncProduct(iChemDir) = PartFuncProduct(iChemDir)   &
          * SpecDSMC(iSpec)%PartitionFunction(LowerLevel)**(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir))
      END IF
    END DO
  END DO
  IF((PartFuncProduct(1).NE.0.).AND.(PartFuncProduct(2).NE.0.)) THEN
    k_b_lower = ChemReac%Arrhenius_Prefactor(iReac) * (LowerLevel * DSMC%PartitionInterval)**ChemReac%Arrhenius_Powerfactor(iReac) &
              * (PartFuncProduct(1)/PartFuncProduct(2))
  ELSE
    k_b_lower = 0.0
  END IF
! Calculation of the backward reaction rate at the upper temperature value (using the equilibrium constant)
  PartFuncProduct(1:2) = 1.
  DO iSpec = 1, nSpecies
    DO iChemDir = 1,2
      IF(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir).NE.0) THEN
        PartFuncProduct(iChemDir) = PartFuncProduct(iChemDir)   &
          * SpecDSMC(iSpec)%PartitionFunction(UpperLevel)**(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir))
      END IF
    END DO
  END DO
  IF((PartFuncProduct(1).NE.0.).AND.(PartFuncProduct(2).NE.0.)) THEN
    k_b_upper = ChemReac%Arrhenius_Prefactor(iReac) * (UpperLevel * DSMC%PartitionInterval)**ChemReac%Arrhenius_Powerfactor(iReac) &
            * (PartFuncProduct(1)/PartFuncProduct(2))
  ELSE
    k_b_upper = 0.0
  END IF
! Linear interpolation of the backward rate coefficient at the actual temperature
  BackwardRate = k_b_lower &
            + (k_b_upper - k_b_lower)  &
            / (DSMC%PartitionInterval) * (LocalTemp - LowerLevel * DSMC%PartitionInterval)


END SUBROUTINE CalcBackwardRate


RECURSIVE FUNCTION lacz_gamma(a) RESULT(g)
!===================================================================================================================================
! gamma function taken from
! http://rosettacode.org/wiki/Gamma_function#Fortran
! variefied against build in and compiled with double precision
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  REAL(KIND=8), INTENT(IN) :: a
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(KIND=8) :: g
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(KIND=8), PARAMETER :: pi = 3.14159265358979324
  INTEGER, PARAMETER :: cg = 7
  ! these precomputed values are taken by the sample code in Wikipedia,
  ! and the sample itself takes them from the GNU Scientific Library
  REAL(KIND=8), DIMENSION(0:8), PARAMETER :: p = &
       (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
       771.32342877765313, -176.61502916214059, 12.507343278686905, &
       -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
  REAL(KIND=8) :: t, w, x
  INTEGER :: i
!===================================================================================================================================

  x = a

  IF ( x < 0.5 ) THEN
     g = pi / ( sin(pi*x) * lacz_gamma(1.0-x) )
  ELSE
     x = x - 1.0
     t = p(0)
     DO i=1, cg+2
        t = t + p(i)/(x+real(i))
     END DO
     w = x + real(cg) + 0.5
     g = sqrt(2.0*pi) * w**(x+0.5) * exp(-w) * t
  END IF
END FUNCTION lacz_gamma


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
  REAL(KIND=real_kind)           :: gamser, gln, ap, del, summ, an, ser, tmp, x,y, b,c,d,h
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
    gammainc=exp(-arg(2)+arg(1)*log(arg(2))-gln) * h
  END IF
END FUNCTION gammainc

END MODULE MOD_DSMC_ChemReact
