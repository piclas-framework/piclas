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

MODULE MOD_DSMC_ChemReact
!===================================================================================================================================
! Module for chemical reactions including calculation of probabilities and collisions
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
PUBLIC :: DSMC_Chemistry, CalcReactionProb, CalcBackwardRate
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Calculates the reaction probability for dissociation, exchange, recombination and associative ionization reactions
!> using the Arrhenius-type reaction rate
!===================================================================================================================================
SUBROUTINE CalcReactionProb(iPair,iReac,ReactionProb,nPair,NumDens)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, maxEXP
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, DSMC, SpecDSMC, PartIntEn, ChemReac, CollInf, ReactionProbGTUnityCounter
USE MOD_Particle_Vars           ,ONLY: PartState, Species, PartSpecies, nSpecies, UseVarTimeStep, usevMPF
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcTVibPoly, CalcTelec
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_DSMC_QK_Chemistry       ,ONLY: QK_GetAnalyticRate
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, nPair
REAL, INTENT(IN)              :: NumDens
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: ReactionProb
!-----------------------------------------------------------------------------------------------------------------------------------
! IN-OUTPUT VARIABLES
INTEGER, INTENT(INOUT)        :: iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ReactInx(1:3), ProductReac(1:4), EductReac(1:3), iReacForward, iCase, iPart
REAL                          :: EZeroPoint_Educt, EZeroPoint_Prod, EReact, ReducedMass, ReducedMassUnweighted, omega, Tref
REAL                          :: Xi_vib(1:3), Xi_elec(1:3), Xi_Total
REAL                          :: BetaReaction, BackwardRate
REAL                          :: Rcoll, Tcoll, Telec, NumWeightEduct, NumWeightProd
INTEGER                       :: iPath, PathIndex
REAL                          :: Weight(1:4), SumWeightEduct, SumWeightProd, betaEXP
!===================================================================================================================================

Weight = 0.; ReactInx = 0
NumWeightEduct = 2.; NumWeightProd = 2.
SumWeightEduct = 0.; SumWeightProd = 0.

EductReac(1:3) = ChemReac%Reactants(iReac,1:3)
iCase = CollInf%Coll_Case(EductReac(1), EductReac(2))

IF(EductReac(3).NE.0) THEN
  IF((iPair.GE.ChemReac%LastPairForRec).AND.(ChemReac%LastPairForRec.GT.0)) THEN
    ReactInx(3) = 0
  ELSE IF(ChemReac%RecombParticle.EQ.0) THEN
    IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
      ChemReac%LastPairForRec = nPair - ChemReac%nPairForRec
      ReactInx(3) = Coll_pData(ChemReac%LastPairForRec)%iPart_p1
      ChemReac%RecombParticle = ReactInx(3)
    ELSE
      ReactInx(3) = 0
    END IF
  ELSE
    ReactInx(3) = ChemReac%RecombParticle
  END IF
  IF(ReactInx(3).GT.0) THEN
    NumWeightEduct = 3.
    EductReac(3) = PartSpecies(ReactInx(3))
    ! Determine the number of the reaction path for this collision pair with the help of the reaction index
    DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
      IF(iReac.EQ.ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)) PathIndex = iPath
    END DO
    iReac = ChemReac%ReactNumRecomb(EductReac(1), EductReac(2), EductReac(3))
    ! Save the new reaction index (depending on the third partner) for the case to be used in DSMC_Chemistry
    ChemReac%CollCaseInfo(iCase)%ReactionIndex(PathIndex) = iReac
  ELSE
    ! If no third collision partner can be found e.g. last (available) pair, set reaction probability to zero and leave routine
    ReactionProb = 0.
    RETURN
  END IF
END IF

! Saving the product species, might be different for the recombination case as reaction index could have changed depending on the
! third reaction partner
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

IF (ChemReac%Reactants(iReac,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  ReactInx(1) = Coll_pData(iPair)%iPart_p1
  ReactInx(2) = Coll_pData(iPair)%iPart_p2
ELSE
  ReactInx(1) = Coll_pData(iPair)%iPart_p2
  ReactInx(2) = Coll_pData(iPair)%iPart_p1
END IF

DO iPart = 1, NINT(NumWeightEduct)
  Weight(iPart) = GetParticleWeight(ReactInx(iPart))
  SumWeightEduct = SumWeightEduct + Weight(iPart)
END DO

SumWeightProd = Weight(1) + Weight(2)
IF(ProductReac(3).NE.0) THEN
  NumWeightProd = 3.
  IF(EductReac(3).EQ.0) Weight(3) = Weight(1)
  SumWeightProd = SumWeightProd + Weight(3)
END IF
IF(ProductReac(4).NE.0) THEN
  NumWeightProd = 4.
  Weight(4) = Weight(1)
  SumWeightProd = SumWeightProd + Weight(4)
END IF

IF (UseVarTimeStep.OR.usevMPF) THEN
  ReducedMass = (Species(EductReac(1))%MassIC *Weight(1)  * Species(EductReac(2))%MassIC * Weight(2)) &
    / (Species(EductReac(1))%MassIC * Weight(1)+ Species(EductReac(2))%MassIC * Weight(2))
  ReducedMassUnweighted = ReducedMass * 2./(Weight(1)+Weight(2))
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF

!---------------------------------------------------------------------------------------------------------------------------------
! Calculation of the collision energy
!---------------------------------------------------------------------------------------------------------------------------------
Coll_pData(iPair)%Ec = 0.5 * ReducedMass*Coll_pData(iPair)%cRela2
DO iPart = 1, NINT(NumWeightEduct)
  IF((Species(EductReac(iPart))%InterID.EQ.2).OR.(Species(EductReac(iPart))%InterID.EQ.20)) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(iPart))%EVib(1)*Weight(iPart) &
                                              + PartIntEn(ReactInx(iPart))%ERot(1)*Weight(iPart)
  END IF
  IF (DSMC%ElectronicModel.GT.0) THEN
    IF((Species(EductReac(iPart))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(iPart))%FullyIonized)) &
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(iPart))%EElec(1)*Weight(iPart)
  END IF
END DO
! Add the translational energy of the third particle in case of a recombination reaction
IF(EductReac(3).NE.0) THEN
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + 0.5 * Species(EductReac(3))%MassIC * DOTPRODUCT(PartState(4:6,ReactInx(3)))
END IF
!---------------------------------------------------------------------------------------------------------------------------------
! Calculation of the zero-point-energies
!---------------------------------------------------------------------------------------------------------------------------------
EZeroPoint_Educt = 0.0; EZeroPoint_Prod = 0.0
DO iPart = 1, NINT(NumWeightEduct)
  IF((Species(EductReac(iPart))%InterID.EQ.2).OR.(Species(EductReac(iPart))%InterID.EQ.20)) THEN
    EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(iPart))%EZeroPoint * Weight(iPart)
  END IF
END DO
DO iPart = 1, NINT(NumWeightProd)
  IF((Species(ProductReac(iPart))%InterID.EQ.2).OR.(Species(ProductReac(iPart))%InterID.EQ.20)) THEN
    EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(iPart))%EZeroPoint * Weight(iPart)
  END IF
END DO
!---------------------------------------------------------------------------------------------------------------------------------
! Calculation of the reaction probability, if collision energy minus the zero-point energy of the EDUCTS is greater than the
! activation energy AND collision energy minus the zero-point energy of the PRODUCTS is greater than the heat of formation
!---------------------------------------------------------------------------------------------------------------------------------
IF(((Coll_pData(iPair)%Ec-EZeroPoint_Educt).GE.(SumWeightEduct/NumWeightEduct*ChemReac%EActiv(iReac))) .AND. &
    ((Coll_pData(iPair)%Ec-EZeroPoint_Prod).GE.(-1.*SumWeightProd/NumWeightProd*ChemReac%EForm(iReac)))) THEN
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calculation of the vibrational degrees of freedom and electronic shell
  !---------------------------------------------------------------------------------------------------------------------------------
  Xi_vib = 0.0; Xi_elec = 0.0
  DO iPart = 1, NINT(NumWeightEduct)
    IF((Species(EductReac(iPart))%InterID.EQ.2).OR.(Species(EductReac(iPart))%InterID.EQ.20)) THEN
      IF(SpecDSMC(EductReac(iPart))%PolyatomicMol) THEN
        IF (PartIntEn(ReactInx(iPart))%EVib(1).GT.SpecDSMC(EductReac(iPart))%EZeroPoint) THEN
          Xi_vib(iPart) = 2.*(PartIntEn(ReactInx(iPart))%EVib(1)-SpecDSMC(EductReac(iPart))%EZeroPoint)  &
                  / (BoltzmannConst*CalcTVibPoly(PartIntEn(ReactInx(iPart))%EVib(1), EductReac(iPart)))
        END IF
      ELSE
        Xi_vib(iPart) = ChemReac%MeanXiVib_PerIter(EductReac(iPart))
      END IF
    END IF
    IF (DSMC%ElectronicModel.GT.0) THEN
      IF((Species(EductReac(iPart))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(iPart))%FullyIonized)) THEN
        IF(PartIntEn(ReactInx(iPart))%EElec(1).GT.0.0)THEN
          Telec=CalcTelec( PartIntEn(ReactInx(iPart))%EElec(1), EductReac(iPart))
          Xi_elec(iPart)=2.*PartIntEn(ReactInx(iPart))%EElec(1)/(BoltzmannConst*Telec)
        END IF
      END IF
    END IF
  END DO
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Determination of the total degree of freedom
  !---------------------------------------------------------------------------------------------------------------------------------
  omega = CollInf%omega(EductReac(1),EductReac(2))
  Tref = CollInf%Tref(EductReac(1),EductReac(2))
  Xi_Total = 2.*(2.-omega)
  DO iPart = 1, NINT(NumWeightEduct)
    Xi_Total = Xi_Total + Xi_elec(iPart) + Xi_vib(iPart) + SpecDSMC(EductReac(iPart))%Xi_Rot
  END DO
  IF(EductReac(3).NE.0) Xi_Total = Xi_Total + 3.
  ! Calculation of the beta factor (requires the collision pair-specific vibrational and electronic degrees of freedom)
  ! Not calculated when treating recombination reactions with the automatic backward reactions.
  IF((TRIM(ChemReac%ReactType(iReac)).NE.'R').OR.(iReac.LE.ChemReac%NumOfReactWOBackward)) BetaReaction = Calc_Beta_TCE(iReac,Xi_Total)
  ! Zero-point energy of educts is removed from the collision energy utilized for the calculation of the reaction probability
  EReact = NumWeightEduct*(Coll_pData(iPair)%Ec - EZeroPoint_Educt)/ SumWeightEduct
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calculation of the backward reaction rate coefficient and applying to Beta coefficient after Boyd "Modeling backward chemical
  ! rate processes in the direct simulation Monte Carlo method", Phys. Fluids 19, 1261103 (2007)
  !---------------------------------------------------------------------------------------------------------------------------------
  IF(iReac.GT.ChemReac%NumOfReactWOBackward) THEN
    iReacForward = ChemReac%BackwardReacForwardIndx(iReac)
    IF(DSMC%InstantTransTemp(nSpecies+1).GT.0.0) THEN
      CALL CalcBackwardRate(iReac,DSMC%InstantTransTemp(nSpecies+1),BackwardRate)
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'E') THEN
        betaEXP = (ChemReac%EActiv(iReacForward)-ChemReac%EActiv(iReac))/(BoltzmannConst*DSMC%InstantTransTemp(nSpecies+1))
        betaEXP = MIN(maxexp,betaEXP)
        BetaReaction = BetaReaction * BackwardRate / EXP(-betaEXP)  &
          / (ChemReac%Arrhenius_Prefactor(iReac) * DSMC%InstantTransTemp(nSpecies+1)**ChemReac%Arrhenius_Powerfactor(iReac))
      END IF
    ELSE
      BackwardRate = 0.0
      BetaReaction = 0.0
    END IF
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Actual calculation of the reaction probability, different equation for recombination reaction
  !---------------------------------------------------------------------------------------------------------------------------------
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    IF(iReac.GT.ChemReac%NumOfReactWOBackward) THEN
      Tcoll =ReducedMassUnweighted*Coll_pData(iPair)%CRela2 / (BoltzmannConst * 2.*(2.-omega))
      Rcoll = (Tcoll / Tref)**(0.5 - omega) &
              * ChemReac%QKRColl(iCase) / SQRT(ReducedMassUnweighted) * ChemReac%QKTCollCorrFac(iCase)
      ReactionProb = BackwardRate / Rcoll * NumDens
    ELSE
      ! Reaction probability after regular TCE-model
      ReactionProb = BetaReaction * NumDens &
                * EReact**(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 + CollInf%omega(EductReac(3),EductReac(3)))
    END IF
  ELSE
    IF(SpecDSMC(EductReac(2))%PolyatomicMol.OR.SpecDSMC(EductReac(1))%PolyatomicMol) THEN
      ! Energy is multiplied by a factor to increase the resulting exponent and avoid floating overflows for high vibrational
      ! degree of freedom, later the reaction probability is scaled again with the same factor and the respective exponents
      ReactionProb = BetaReaction * ((EReact - ChemReac%EActiv(iReac))*1E6) &
            ** (ChemReac%Arrhenius_Powerfactor(iReac)-1.5+omega+Xi_Total/2.) * (EReact * 1E6)**(1.0 - Xi_Total/2.)
      ReactionProb = ReactionProb / ((1E6)**(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 + omega))
    ELSE
      ReactionProb = BetaReaction * ((EReact - ChemReac%EActiv(iReac))) &
            ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + omega + Xi_Total/2.) * (EReact) ** (1.0 - Xi_Total/2.)
    END IF
  END IF
ELSE
  ReactionProb = 0.0
END IF

IF(DSMC%ReservoirSimu) THEN
  IF(DSMC%ReservoirRateStatistic) THEN
    IF((ReactionProb.GT.1).AND.(ReactionProbGTUnityCounter.LT.100)) THEN
      ReactionProbGTUnityCounter=ReactionProbGTUnityCounter+1
      IPWRITE(*,*) 'Warning: ReactionProb greater than unity! ReacNbr:', iReac,'    ReactionProb:',ReactionProb
      IF(ReactionProbGTUnityCounter.EQ.100)THEN
        IPWRITE(*,*) ' Counted 100 ReactionProb greater than unity. Turning this warning off.'
      END IF
    END IF
  END IF
END IF
! ReactionProb should not be gt 1 to avoid meaningless high weighting of a single reaction
IF (ReactionProb.GT.1) ReactionProb = 1.0
IF (DSMC%ReservoirSimu.AND..NOT.DSMC%ReservoirRateStatistic) THEN
  ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb * SumWeightEduct / NumWeightEduct
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
END IF

END SUBROUTINE CalcReactionProb


!===================================================================================================================================
!> Calculates the Beta coefficient for polyatomic reactions
!===================================================================================================================================
PPURE REAL FUNCTION Calc_Beta_TCE(iReac,Xi_Total)
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars       ,ONLY: ChemReac, CollInf
USE MOD_Globals_Vars    ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: iReac
REAL, INTENT(IN)        :: Xi_Total
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: omega
!===================================================================================================================================

omega = CollInf%omega(ChemReac%Reactants(iReac,1),ChemReac%Reactants(iReac,2))

IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
  Calc_Beta_TCE = ChemReac%Arrhenius_Prefactor(iReac) * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) - omega) &
    * GAMMA(Xi_total/2. - omega + 2.) / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) + Xi_total/2. + 1.5))
ELSE
  IF((ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 + omega + Xi_Total/2.).GT.0.0) THEN
    Calc_Beta_TCE = ChemReac%Arrhenius_Prefactor(iReac) * (BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) - omega)) &
    * GAMMA(Xi_Total/2.) / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) + Xi_Total/2. - 0.5 + omega))
  ELSE
    Calc_Beta_TCE = 0.0
  END IF
END IF

END FUNCTION Calc_Beta_TCE


!===================================================================================================================================
!> Routine performs a reaction of the type A + B + C -> D + E + F + G, where A, B, C, D, E, F can be anything
!===================================================================================================================================
SUBROUTINE DSMC_Chemistry(iPair, iReac)
! MODULES
USE MOD_Globals                ,ONLY: abort,DOTPRODUCT,StringBeginsWith,UNIT_StdOut
USE MOD_Globals_Vars
USE MOD_DSMC_Vars              ,ONLY: Coll_pData, DSMC, CollInf, SpecDSMC
USE MOD_DSMC_Vars              ,ONLY: ChemReac, PartIntEn, PolyatomMolDSMC, BGGas, ElecRelaxPart
USE MOD_DSMC_Vars              ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi, newElecRelaxParts, iPartIndx_NodeNewElecRelax
USE MOD_DSMC_Vars              ,ONLY: iPartIndx_NodeElecRelaxChem, nElecRelaxChemParts
USE MOD_Particle_Vars          ,ONLY: PartSpecies, PartState, PDM, PEM, PartPosRef, Species, PartMPF, usevMPF
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_DSMC_ElectronicModel   ,ONLY: ElectronicEnergyExchange
USE MOD_DSMC_PolyAtomicModel   ,ONLY: RotRelaxPolyRoutineFuncPTR, DSMC_RelaxVibPolyProduct
USE MOD_DSMC_Relaxation        ,ONLY: DSMC_VibRelaxDiatomic, RotRelaxDiaRoutineFuncPTR
USE MOD_DSMC_CollisVec         ,ONLY: PostCollVec
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Analyze_Tools ,ONLY: CalcXiTotalEqui
USE MOD_part_tools             ,ONLY: GetParticleWeight, GetNextFreePosition
USE MOD_part_operations        ,ONLY: RemoveParticle
#ifdef CODE_ANALYZE
USE MOD_Globals                ,ONLY: unit_stdout
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
#endif /* CODE_ANALYZE */
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance,nPartIn,nPartOut,PartEkinIn,PartEkinOut
USE MOD_Particle_Analyze_Pure   ,ONLY: CalcEkinPart
#if USE_MPI
USE MOD_Globals                ,ONLY: myrank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ReactInx(1:4), EductReac(1:3), ProductReac(1:4), nDOFMAX, iProd, iPolyatMole, iSpec
INTEGER                       :: SpecToDelete, iPart, NumEduct, NumProd, LocalElemID
REAL                          :: FracMassCent1, FracMassCent2, MassRed      ! mx/(mx+my)
REAL                          :: VeloCOM(1:3)                               !> Centre of mass velocity
REAL                          :: omega, FakXi, Xi_total, iRan, FacEtraDistri
REAL                          :: ERel_React1_React2, ERel_React1_React3, ERel_React2_React4
REAL                          :: Xi_elec(1:4), EZeroTempToExec(1:4)
REAL, ALLOCATABLE             :: XiVibPart(:,:)
REAL                          :: Weight(1:4), SumWeightProd, WeightProd(1:4)
REAL                          :: cRelaNew(3), TempVelo(3)
#ifdef CODE_ANALYZE
REAL,PARAMETER                :: RelMomTol=5e-9  ! Relative tolerance applied to conservation of momentum before/after reaction
REAL,PARAMETER                :: RelEneTol=5e-12 ! Relative tolerance applied to conservation of energy before/after reaction
!                                                ! 2024/1: Increased from 2e-12 to 5e-12 due to reggies failing when using CORE_SPLIT
REAL                          :: Energy_old,Energy_new,Momentum_old(3),Momentum_new(3)
INTEGER                       :: iMom, iMomDim
#endif /* CODE_ANALYZE */
!===================================================================================================================================

EductReac(1:3) = ChemReac%Reactants(iReac,1:3)
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

IF(EductReac(3).NE.0) THEN
  IF(ChemReac%RecombParticle.NE.0) THEN
    ReactInx(3) = ChemReac%RecombParticle
    ! If a pair ahead in the list was broken, check which of the particles is used
    IF(ChemReac%LastPairForRec.GT.0) THEN
      Coll_pData(ChemReac%LastPairForRec)%NeedForRec = .TRUE.
      IF(ChemReac%RecombParticle.EQ.Coll_pData(ChemReac%LastPairForRec)%iPart_p2) THEN
        ! Both particles of the pair were already used, set the RecombParticle index to zero and advance the nPairForRec counter
        ! -> New pair will be broken up (in CalcReactionProb routine)
        ChemReac%RecombParticle = 0
        ChemReac%nPairForRec = ChemReac%nPairForRec + 1
      ELSE
        ! Utilize the second particle of the pair for the recombination
        ChemReac%RecombParticle = Coll_pData(ChemReac%LastPairForRec)%iPart_p2
      END IF
    ELSE
      ChemReac%RecombParticle = 0
    END IF
  ELSE
    CALL Abort(__STAMP__,'New Particle Number greater max Part Num in DSMC_Chemistry. Reaction: ',IntInfoOpt=iReac)
  END IF
END IF

Xi_elec = 0.
EZeroTempToExec = 0.

Weight = 0.
WeightProd = 0.
NumEduct = 2; NumProd = 2; SumWeightProd = 0.

! Get the index of reactants from the collision pair
IF((    PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%Reactants(iReac,1)) &
  .AND.(PartSpecies(Coll_pData(iPair)%iPart_p2).EQ.ChemReac%Reactants(iReac,2))) THEN
  ReactInx(1) = Coll_pData(iPair)%iPart_p1
  ReactInx(2) = Coll_pData(iPair)%iPart_p2
ELSEIF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%Reactants(iReac,2)&
  .AND.(PartSpecies(Coll_pData(iPair)%iPart_p2).EQ.ChemReac%Reactants(iReac,1))) THEN
  ReactInx(2) = Coll_pData(iPair)%iPart_p1
  ReactInx(1) = Coll_pData(iPair)%iPart_p2
ELSE
  CALL abort(__STAMP__,'ERROR in DSMC_Chemistry: Chosen pair does not correspond to the reactants!')
END IF

! Species-specific time step: store the velocity of species at the lower time step
IF(VarTimeStep%UseSpeciesSpecific) THEN
  IF(VarTimeStep%DisableForMCC) THEN
    IF(Species(PartSpecies(ReactInx(1)))%TimeStepFactor.LT.1.) THEN
      TempVelo(1:3) = PartState(4:6,ReactInx(1))
    ELSE IF(Species(PartSpecies(ReactInx(2)))%TimeStepFactor.LT.1.) THEN
      TempVelo(1:3) = PartState(4:6,ReactInx(2))
    END IF
  END IF
END IF

! Third reactant is available: set the species and number of reactants and products in case of recombination and exchange (3 -> 3)
IF(EductReac(3).NE.0) THEN
  NumEduct = 3
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    EductReac(3) = PartSpecies(ReactInx(3))
    ProductReac(2) = PartSpecies(ReactInx(3))
  END IF
  IF(ProductReac(3).NE.0) THEN
    PartSpecies(ReactInx(3)) = ProductReac(3)
    NumProd = 3
  END IF
END IF

! Set the particle weights of the reactants
DO iPart = 1, NumEduct
  Weight(iPart) = GetParticleWeight(ReactInx(iPart))
END DO

! Do not perform the reaction in case the reaction is to be calculated at a constant gas composition (DSMC%ReservoirSimuRate = T)
IF (DSMC%ReservoirSimu.AND.DSMC%ReservoirSimuRate) THEN
  ! Count the number of reactions to determine the actual reaction rate
  IF (DSMC%ReservoirRateStatistic) ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + SUM(Weight(:)) / NumEduct
  ! Leave the routine again
  RETURN
END IF

! Set the particle weights of the products (they are equal to the weights of the reactants except in the case of a recombination)
IF ((usevMPF).AND.(TRIM(ChemReac%ReactType(iReac)).EQ.'R')) THEN
  ! Set the weight of the non-reacting collision partner
  WeightProd(2) = Weight(3)
  ! Calculate the new mass-weighted weight
  PartMPF(ReactInx(1)) = (Species(EductReac(1))%MassIC*PartMPF(ReactInx(1)) + Species(EductReac(2))%MassIC*PartMPF(ReactInx(2))) &
    /(Species(EductReac(1))%MassIC+Species(EductReac(2))%MassIC)
  ! Get the weight in case of variable time step in combination with vMPF
  WeightProd(1) = GetParticleWeight(ReactInx(1))
ELSE
  WeightProd = Weight
END IF

! Add reactants to the counter of outgoing particles
IF(CalcPartBalance) THEN
  IF(.NOT.StringBeginsWith(ChemReac%ReactModel(iReac),'phIon'))THEN
    DO iPart = 1, NumEduct
      IF(BGGas%BackgroundSpecies(EductReac(iPart))) CYCLE
      nPartOut(EductReac(iPart))=nPartOut(EductReac(iPart)) + 1
      PartEkinOut(EductReac(iPart))=PartEkinOut(EductReac(iPart))+CalcEkinPart(ReactInx(iPart))
    END DO
  END IF
END IF

! Calculate the reduced mass (Special case: photo-ionization, where the other reactant is a photon)
IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon'))THEN
  MassRed = 0.
ELSE
  IF (UseVarTimeStep.OR.usevMPF) THEN
    MassRed = Species(EductReac(1))%MassIC*Weight(1) * Species(EductReac(2))%MassIC*Weight(2) &
                / (Species(EductReac(1))%MassIC*Weight(1) + Species(EductReac(2))%MassIC*Weight(2))
  ELSE
    MassRed = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF
END IF

! Calculate the sum of the weights of the products (at least 2 products in any reaction)
SumWeightProd = WeightProd(1) + WeightProd(2)

IF(EductReac(3).EQ.0) THEN
  IF(ProductReac(3).NE.0) THEN
    ! === Get free particle index for the 3rd product
    ReactInx(3) = GetNextFreePosition()
    PDM%ParticleInside(ReactInx(3)) = .true.
    PDM%IsNewPart(ReactInx(3)) = .true.
    PDM%dtFracPush(ReactInx(3)) = .FALSE.
    ! Set species index of new particle
    PartSpecies(ReactInx(3)) = ProductReac(3)
    PartState(1:3,ReactInx(3)) = PartState(1:3,ReactInx(1))
    IF(TrackingMethod.EQ.REFMAPPING) THEN
      PartPosRef(1:3,ReactInx(3))=PartPosRef(1:3,ReactInx(1))
    END IF
    IF((Species(ProductReac(3))%InterID.EQ.2).OR.(Species(ProductReac(3))%InterID.EQ.20)) THEN
      ALLOCATE(PartIntEn(ReactInx(3))%ERot(1), PartIntEn(ReactInx(3))%EVib(1))
      PartIntEn(ReactInx(3))%EVib(1)= 0.
      PartIntEn(ReactInx(3))%ERot(1) = 0.
    END IF
    IF(DSMC%ElectronicModel.GT.0) THEN
      IF((Species(ProductReac(3))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(3))%FullyIonized)) THEN
        ALLOCATE(PartIntEn(ReactInx(3))%EElec(1))
        PartIntEn(ReactInx(3))%EElec(1) = 0.
      END IF
    END IF
    PEM%GlobalElemID(ReactInx(3)) = PEM%GlobalElemID(ReactInx(1))
    PEM%LastGlobalElemID(ReactInx(3)) = PEM%GlobalElemID(ReactInx(3))
    IF(usevMPF) PartMPF(ReactInx(3)) = PartMPF(ReactInx(1))
    IF(UseVarTimeStep) PartTimeStep(ReactInx(3)) = PartTimeStep(ReactInx(1))
    IF (DSMC%DoAmbipolarDiff) THEN
      newAmbiParts = newAmbiParts + 1
      iPartIndx_NodeNewAmbi(newAmbiParts) = ReactInx(3)
    END IF
    IF (DSMC%ElectronicModel.EQ.4) THEN
      newElecRelaxParts = newElecRelaxParts + 1
      iPartIndx_NodeNewElecRelax(newElecRelaxParts) = ReactInx(3)
    END IF
    WeightProd(3) = Weight(1)
    NumProd = 3
    SumWeightProd = SumWeightProd + WeightProd(3)
    ! Particle index linked list (required for merge)
    LocalElemID = PEM%LocalElemID(ReactInx(3))
    PEM%pNext(PEM%pEnd(LocalElemID)) = ReactInx(3)
    PEM%pEnd(LocalElemID) = ReactInx(3)
    PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
  END IF
END IF

IF(ProductReac(4).NE.0) THEN
  ! === Get free particle index for the 4th product
  ReactInx(4) = GetNextFreePosition()
  PDM%ParticleInside(ReactInx(4)) = .true.
  PDM%IsNewPart(ReactInx(4)) = .true.
  PDM%dtFracPush(ReactInx(4)) = .FALSE.
  ! Set species index of new particle
  PartSpecies(ReactInx(4)) = ProductReac(4)
  PartState(1:3,ReactInx(4)) = PartState(1:3,ReactInx(1))
  IF(TrackingMethod.EQ.REFMAPPING) THEN ! here Nearst-GP is missing
    PartPosRef(1:3,ReactInx(4))=PartPosRef(1:3,ReactInx(1))
  END IF
  IF((Species(ProductReac(4))%InterID.EQ.2).OR.(Species(ProductReac(4))%InterID.EQ.20)) THEN
    ALLOCATE(PartIntEn(ReactInx(4))%ERot(1), PartIntEn(ReactInx(4))%EVib(1))
    PartIntEn(ReactInx(4))%EVib(1)= 0.
    PartIntEn(ReactInx(4))%ERot(1) = 0.
  END IF
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(ProductReac(4))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(4))%FullyIonized)) THEN
      ALLOCATE(PartIntEn(ReactInx(4))%EElec(1))
      PartIntEn(ReactInx(4))%EElec(1) = 0.
    END IF
  END IF
  PEM%GlobalElemID(ReactInx(4)) = PEM%GlobalElemID(ReactInx(1))
  PEM%LastGlobalElemID(ReactInx(4)) = PEM%GlobalElemID(ReactInx(4))
  IF(usevMPF) PartMPF(ReactInx(4)) = PartMPF(ReactInx(1))
  IF(UseVarTimeStep) PartTimeStep(ReactInx(4)) = PartTimeStep(ReactInx(1))
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = ReactInx(4)
  END IF
  IF (DSMC%ElectronicModel.EQ.4) THEN
    newElecRelaxParts = newElecRelaxParts + 1
    iPartIndx_NodeNewElecRelax(newElecRelaxParts) = ReactInx(4)
  END IF
  WeightProd(4) = Weight(1)
  NumProd = 4
  SumWeightProd = SumWeightProd + WeightProd(4)
  ! Particle index linked list (required for merge)
  LocalElemID = PEM%LocalElemID(ReactInx(4))
  PEM%pNext(PEM%pEnd(LocalElemID)) = ReactInx(4)
  PEM%pEnd(LocalElemID) = ReactInx(4)
  PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
END IF

#ifdef CODE_ANALYZE
! Energy conservation
Energy_old=0.5*Species(PartSpecies(ReactInx(1)))%MassIC*DOTPRODUCT(PartState(4:6,ReactInx(1))) * Weight(1) &
          +0.5*Species(PartSpecies(ReactInx(2)))%MassIC*DOTPRODUCT(PartState(4:6,ReactInx(2))) * Weight(2) &
          + ChemReac%EForm(iReac)*SumWeightProd/REAL(NumProd)
IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) &
  Energy_old = Energy_old + (PartIntEn(ReactInx(1))%EVib(1) + PartIntEn(ReactInx(1))%ERot(1)) * Weight(1)
IF((Species(EductReac(2))%InterID.EQ.2).OR.(Species(EductReac(2))%InterID.EQ.20)) &
  Energy_old = Energy_old + (PartIntEn(ReactInx(2))%EVib(1) + PartIntEn(ReactInx(2))%ERot(1)) * Weight(2)
IF(DSMC%ElectronicModel.GT.0) THEN
  IF((Species(EductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(1))%FullyIonized)) &
    Energy_old=Energy_old + PartIntEn(ReactInx(1))%EElec(1)*Weight(1)
  IF((Species(EductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(2))%FullyIonized)) &
    Energy_old=Energy_old + PartIntEn(ReactInx(2))%EElec(1)*Weight(2)
END IF
IF (EductReac(3).NE.0) THEN
  Energy_old=Energy_old+0.5*Species(PartSpecies(ReactInx(3)))%MassIC*DOTPRODUCT(PartState(4:6,ReactInx(3)))* Weight(3)
  IF((Species(EductReac(3))%InterID.EQ.2).OR.(Species(EductReac(3))%InterID.EQ.20)) &
     Energy_old=Energy_old + (PartIntEn(ReactInx(3))%EVib(1) + PartIntEn(ReactInx(3))%ERot(1)) * Weight(3)
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(EductReac(3))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(3))%FullyIonized)) &
      Energy_old=Energy_old + PartIntEn(ReactInx(3))%EElec(1) * Weight(3)
  END IF
END IF
! Momentum conservation
Momentum_old(1:3) = Species(PartSpecies(ReactInx(1)))%MassIC * PartState(4:6,ReactInx(1)) * Weight(1) &
                  + Species(PartSpecies(ReactInx(2)))%MassIC * PartState(4:6,ReactInx(2)) * Weight(2)
IF (EductReac(3).NE.0) Momentum_old(1:3) = Momentum_old(1:3) + Species(PartSpecies(ReactInx(3)))%MassIC &
                                                                * PartState(4:6,ReactInx(3)) * Weight(3)
#endif /* CODE_ANALYZE */

! This might fall below zero
Coll_pData(iPair)%Ec = 0.5 * MassRed * Coll_pData(iPair)%CRela2 + ChemReac%EForm(iReac)*SumWeightProd/REAL(NumProd)

!-------------------------------------------------------------------------------------------------------------------------------
! Relative translational degrees of freedom
!-------------------------------------------------------------------------------------------------------------------------------
IF(ProductReac(4).NE.0) THEN
  ! 4 Products
  ! Determine the omega of the pseudo pair (1-3) and (2-4)
  omega = (CollInf%omega(ProductReac(1),ProductReac(3)) + CollInf%omega(ProductReac(2),ProductReac(4))) / 2.
  ! Consider the relative degrees of freedom between the 1. and 3. product, 2. and 4. product and the between the pseudo pairs
  Xi_total = 2.*(2. - CollInf%omega(ProductReac(1),ProductReac(3))) + 2.*(2. - CollInf%omega(ProductReac(2),ProductReac(4))) &
              + 2.*(2. - omega)
ELSE IF(ProductReac(3).NE.0) THEN
  ! 3 Products
  ! Determine the omega of the pseudo pair (1-3) and the 2. product
  omega = (CollInf%omega(ProductReac(1),ProductReac(3)) + CollInf%omega(ProductReac(2),ProductReac(2))) / 2.
  ! Consider the relative degrees of freedom between the 1. and 3. product and the 2. product and the pseudo pair (1-3)
  Xi_total = 2.*(2. - CollInf%omega(ProductReac(1),ProductReac(3))) + 2.*(2. - omega)
ELSE
  ! 2 Products
  Xi_total = 2.*(2. - CollInf%omega(ProductReac(1),ProductReac(2)))
END IF
!-------------------------------------------------------------------------------------------------------------------------------
! Rotational degrees of freedom
!-------------------------------------------------------------------------------------------------------------------------------
DO iProd = 1, NumProd
  Xi_total = Xi_total + SpecDSMC(ProductReac(iProd))%Xi_Rot
END DO
!-------------------------------------------------------------------------------------------------------------------------------
! Adding the vibrational and rotational energy to the collision energy
!-------------------------------------------------------------------------------------------------------------------------------
IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) &
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (PartIntEn(ReactInx(1))%EVib(1) + PartIntEn(ReactInx(1))%ERot(1))*Weight(1)
IF((Species(EductReac(2))%InterID.EQ.2).OR.(Species(EductReac(2))%InterID.EQ.20)) &
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (PartIntEn(ReactInx(2))%EVib(1) + PartIntEn(ReactInx(2))%ERot(1))*Weight(2)
!-------------------------------------------------------------------------------------------------------------------------------
! Addition of the electronic energy to the collision energy)
!-------------------------------------------------------------------------------------------------------------------------------
IF (DSMC%ElectronicModel.GT.0) THEN
  IF((Species(EductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(1))%FullyIonized)) &
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(1))%EElec(1)*Weight(1)
  IF((Species(EductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(2))%FullyIonized)) &
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(2))%EElec(1)*Weight(2)
END IF
! Sanity check
IF(Coll_pData(iPair)%Ec.LT.0)THEN
  IPWRITE(UNIT_StdOut,*) "DSMC Chemistry: Added vibrational and electronic energy"
  IPWRITE(UNIT_StdOut,*) "iReac       =", iReac
  IPWRITE(UNIT_StdOut,*) "ReactInx(1) =", ReactInx(1), " ReactInx(2) =", ReactInx(2)
  IPWRITE(UNIT_StdOut,*) "Weight(1)   =", Weight(1)  , " Weight(2)   =", Weight(2)
  IPWRITE(UNIT_StdOut,*) "Coll_pData(iPair)%Ec     = ", Coll_pData(iPair)%Ec
  IPWRITE(UNIT_StdOut,*) "MassRed                  = ", MassRed
  IPWRITE(UNIT_StdOut,*) "ChemReac%EForm(iReac)    =", ChemReac%EForm(iReac)
  IPWRITE(UNIT_StdOut,*) "Coll_pData(iPair)%CRela2 = ", Coll_pData(iPair)%CRela2
  IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) &
    IPWRITE(UNIT_StdOut,*) "PartIntEn(ReactInx(1))%EVib,ERot =", PartIntEn(ReactInx(1))%EVib, PartIntEn(ReactInx(1))%ERot
  IF((Species(EductReac(2))%InterID.EQ.2).OR.(Species(EductReac(2))%InterID.EQ.20)) &
    IPWRITE(UNIT_StdOut,*) "PartIntEn(ReactInx(2))%EVib,ERot =", PartIntEn(ReactInx(2))%EVib, PartIntEn(ReactInx(2))%ERot
  IF(DSMC%ElectronicModel.GT.0)THEN
    IF((Species(EductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(1))%FullyIonized)) &
      IPWRITE(UNIT_StdOut,*) "PartIntEn(ReactInx(1))%EElec =", PartIntEn(ReactInx(1))%EElec
    IF((Species(EductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(2))%FullyIonized)) &
      IPWRITE(UNIT_StdOut,*) "PartIntEn(ReactInx(2))%EElec =", PartIntEn(ReactInx(2))%EElec
  END IF ! DSMC%ElectronicModel.GT.0
  CALL abort(__STAMP__,'Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
END IF ! Coll_pData(iPair)%Ec.LT.0

IF(EductReac(3).NE.0) THEN
  ! If a third collision partner exists (recombination/exchange reactions with defined third educt, A + B+ C), calculation of
  ! the centre of mass of a pseudo-molecule consisting of the first two educts -> (AB) + C
  IF (usevMPF) THEN
    FracMassCent1 = Species(EductReac(1))%MassIC *Weight(1) &
        /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
    FracMassCent2 = Species(EductReac(2))%MassIC *Weight(2) &
        /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(EductReac(1), Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(EductReac(2), Coll_pData(iPair)%PairType)
  END IF

  ! Centre of mass velocity between the first two products
  VeloCOM(1:3) = FracMassCent1 * PartState(4:6,ReactInx(1)) + FracMassCent2 * PartState(4:6,ReactInx(2))

  ! Overwriting the PartState of the first particle with the new PartState of the pseudo-molecule (AB)
  PartState(4:6,ReactInx(1)) = VeloCOM(1:3)

  ! Calculation of the reduced mass of the pseudo-molecule and third collision partner
  CALL CalcPseudoScatterVars(EductReac(1),EductReac(2),EductReac(3),FracMassCent1,FracMassCent2,MassRed,Weight(1:3))
  ! Addition of the relative translation energy between (AB) and C, rotational and vibrational energy of the third
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + 0.5 * MassRed * ((VeloCOM(1)-PartState(4,ReactInx(3)))**2 &
                                                               + (VeloCOM(2)-PartState(5,ReactInx(3)))**2 &
                                                               + (VeloCOM(3)-PartState(6,ReactInx(3)))**2)
  IF((Species(EductReac(3))%InterID.EQ.2).OR.(Species(EductReac(3))%InterID.EQ.20)) &
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (PartIntEn(ReactInx(3))%EVib(1) + PartIntEn(ReactInx(3))%ERot(1)) * Weight(3)
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(EductReac(3))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(3))%FullyIonized)) &
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(3))%EElec(1) * Weight(3)
  END IF
END IF

!-------------------------------------------------------------------------------------------------------------------------------
! Redistribution of collisional energy according to the equipartion theorem
!-------------------------------------------------------------------------------------------------------------------------------
! Determining the maximal number of vibrational SHOs for allocation of the XiVibPart array
nDOFMAX = 0
DO iProd = 1, NumProd
  IF((Species(ProductReac(iProd))%InterID.EQ.2).OR.(Species(ProductReac(iProd))%InterID.EQ.20)) THEN
    IF(SpecDSMC(ProductReac(iProd))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ProductReac(iProd))%SpecToPolyArray
      nDOFMAX = MAX(nDOFMAX,PolyatomMolDSMC(iPolyatMole)%VibDOF)
    ELSE
      nDOFMAX = MAX(nDOFMAX,1)
    END IF
  END IF
END DO
!-------------------------------------------------------------------------------------------------------------------------------
! Root-finding algorithm to determine the vibrational and electronic degrees of freedom
IF((nDOFMAX.GT.0).AND.((DSMC%ElectronicModel.EQ.1).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.4))) THEN
  ! Electronic and vibrational energy is considered
  ALLOCATE(XiVibPart(NumProd,nDOFMAX))
  XiVibPart = 0.
  CALL CalcXiTotalEqui(iReac,iPair,NumProd,Xi_total,WeightProd,XiVibPart=XiVibPart,XiElecPart=Xi_elec)
ELSEIF((DSMC%ElectronicModel.EQ.1).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.4)) THEN
  ! Only electronic energy is considered
  CALL CalcXiTotalEqui(iReac,iPair,NumProd,Xi_total,WeightProd,XiElecPart=Xi_elec)
ELSEIF(nDOFMAX.GT.0) THEN
  ! Only vibrational energy is considered (in case of ElectronicModel = 3, products are always in the ground-state)
  ALLOCATE(XiVibPart(NumProd,nDOFMAX))
  XiVibPart = 0.
  CALL CalcXiTotalEqui(iReac,iPair,NumProd,Xi_total,WeightProd,XiVibPart=XiVibPart)
END IF
!-------------------------------------------------------------------------------------------------------------------------------
! Addition of the vibrational and electronic degrees of freedom to Xi_total
!-------------------------------------------------------------------------------------------------------------------------------
IF(nDOFMAX.GT.0) THEN
  DO iProd = 1, NumProd
    IF((Species(ProductReac(iProd))%InterID.EQ.2).OR.(Species(ProductReac(iProd))%InterID.EQ.20)) THEN
      Xi_total = Xi_total + SUM(XiVibPart(iProd,:))
      EZeroTempToExec(iProd) = SpecDSMC(ProductReac(iProd))%EZeroPoint*WeightProd(iProd)
    END IF
  END DO
END IF
IF (DSMC%ElectronicModel.GT.0) THEN
  DO iProd = 1, NumProd
    Xi_total = Xi_total + Xi_elec(iProd)
  END DO
END IF
FakXi = 0.5*Xi_total - 1.0
!-------------------------------------------------------------------------------------------------------------------------------
! Subtracting the zero-point energy of the products (is added back later)
Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - SUM(EZeroTempToExec(:))
! Sanity check
IF(Coll_pData(iPair)%Ec.LT.0.)THEN
  IPWRITE(UNIT_StdOut,*) "DSMC Chemistry: Zero-point energy"
  IPWRITE(UNIT_StdOut,*) "iReac =", iReac
  CALL abort(__STAMP__,'Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
END IF ! Coll_pData(iPair)%Ec.LT.0.
!-------------------------------------------------------------------------------------------------------------------------------
! Set the new species of the products
DO iProd = 1, NumProd
  PartSpecies(ReactInx(iProd)) = ProductReac(iProd)
END DO
!--------------------------------------------------------------------------------------------------
! Electronic energy exchange
!--------------------------------------------------------------------------------------------------
IF (DSMC%ElectronicModel.GT.0) THEN
  DO iProd = 1, NumProd
    IF((Species(ProductReac(iProd))%InterID.EQ.4).OR.SpecDSMC(ProductReac(iProd))%FullyIonized) THEN
      IF (DSMC%ElectronicModel.EQ.2) THEN
        IF(ALLOCATED(PartIntEn(ReactInx(iProd))%DistriFunc)) DEALLOCATE(PartIntEn(ReactInx(iProd))%DistriFunc)
      END IF
      SDEALLOCATE(PartIntEn(ReactInx(iProd))%EElec)
      IF (DSMC%ElectronicModel.EQ.4) THEN
        nElecRelaxChemParts = nElecRelaxChemParts + 1
        iPartIndx_NodeElecRelaxChem(nElecRelaxChemParts) = ReactInx(iProd)
        ElecRelaxPart(ReactInx(iProd)) = .FALSE.
      END IF
    ELSE
      IF (.NOT.ALLOCATED(PartIntEn(ReactInx(iProd))%EElec)) ALLOCATE(PartIntEn(ReactInx(iProd))%EElec(1))
      IF (DSMC%ElectronicModel.EQ.4) THEN
        PartIntEn(ReactInx(iProd))%EElec = 0.0
        nElecRelaxChemParts = nElecRelaxChemParts + 1
        iPartIndx_NodeElecRelaxChem(nElecRelaxChemParts) = ReactInx(iProd)
        ElecRelaxPart(ReactInx(iProd)) = .FALSE.
      ELSE
        IF (DSMC%ElectronicModel.EQ.2) THEN
          IF(ALLOCATED(PartIntEn(ReactInx(iProd))%DistriFunc)) DEALLOCATE(PartIntEn(ReactInx(iProd))%DistriFunc)
          ALLOCATE(PartIntEn(ReactInx(iProd))%DistriFunc(1:SpecDSMC(ProductReac(iProd))%MaxElecQuant))
          PartIntEn(ReactInx(iProd))%DistriFunc = 0.0
          PartIntEn(ReactInx(iProd))%EElec = 0.0
        END IF
        FakXi = FakXi - 0.5*Xi_elec(iProd)
        CALL ElectronicEnergyExchange(iPair,ReactInx(iProd),FakXi, NewPart = .TRUE., XSec_Level = 0)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartIntEn(ReactInx(iProd))%EElec(1)*WeightProd(iProd)
        IF(Coll_pData(iPair)%Ec.LT.0.) THEN
          IPWRITE(UNIT_StdOut,*) "DSMC Chemistry: Electronic energy exchange"
          IPWRITE(UNIT_StdOut,*) "iProd =", iProd
          IPWRITE(UNIT_StdOut,*) "iReac =", iReac
          CALL abort(__STAMP__,'Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
        END IF
      END IF
    END IF
  END DO
END IF
!--------------------------------------------------------------------------------------------------
! Vibrational energy exchange
!--------------------------------------------------------------------------------------------------
DO iProd = 1, NumProd
  IF((Species(ProductReac(iProd))%InterID.EQ.2).OR.(Species(ProductReac(iProd))%InterID.EQ.20)) THEN
    IF (.NOT.ALLOCATED(PartIntEn(ReactInx(iProd))%EVib)) ALLOCATE(PartIntEn(ReactInx(iProd))%EVib(1))
    FakXi = FakXi - 0.5*XiVibPart(iProd,1)
    IF(SpecDSMC(ProductReac(iProd))%PolyatomicMol) THEN
      ! Zero-point energy is added (for every vibrational dof separately) and new vibrational state is substracted
      ! from the collision energy within the routine
      CALL DSMC_RelaxVibPolyProduct(iPair, ReactInx(iProd), FakXi, XiVibPart(iProd,:), WeightProd(iProd))
    ELSE
      IF(EductReac(iProd).NE.0) THEN
        IF(SpecDSMC(EductReac(iProd))%PolyatomicMol) DEALLOCATE(PartIntEn(ReactInx(iProd))%QVib)
      END IF
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + EZeroTempToExec(iProd)
      CALL DSMC_VibRelaxDiatomic(iPair,ReactInx(iProd),FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartIntEn(ReactInx(iProd))%EVib(1)*WeightProd(iProd)
      ! Sanity check
      IF(Coll_pData(iPair)%Ec.LT.0.) THEN
        IPWRITE(UNIT_StdOut,*) "DSMC Chemistry: Vibrational energy exchange"
        IPWRITE(UNIT_StdOut,*) "iProd =", iProd
        IPWRITE(UNIT_StdOut,*) "iReac =", iReac
        CALL abort(__STAMP__,'Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
      END IF
    END IF
  ELSE
    SDEALLOCATE(PartIntEn(ReactInx(iProd))%EVib)
  END IF
END DO
!--------------------------------------------------------------------------------------------------
! Rotational energy exchange (additional check: If new particle is an atom, internal energies are zero)
!--------------------------------------------------------------------------------------------------
DO iProd = 1, NumProd
  IF ((Species(ProductReac(iProd))%InterID.EQ.2).OR.(Species(ProductReac(iProd))%InterID.EQ.20)) THEN
    IF (.NOT.ALLOCATED(PartIntEn(ReactInx(iProd))%ERot)) ALLOCATE(PartIntEn(ReactInx(iProd))%ERot(1))
    FakXi = FakXi - 0.5*SpecDSMC(ProductReac(iProd))%Xi_Rot
    IF(SpecDSMC(ProductReac(iProd))%PolyatomicMol) THEN
      CALL RotRelaxPolyRoutineFuncPTR(iPair, ReactInx(iProd), FakXi)
    ELSE
      CALL RotRelaxDiaRoutineFuncPTR(iPair, ReactInx(iProd), FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartIntEn(ReactInx(iProd))%ERot(1)*WeightProd(iProd)
    ! Sanity check
    IF(Coll_pData(iPair)%Ec.LT.0.) THEN
      IPWRITE(UNIT_StdOut,*) "DSMC Chemistry: Rotational energy exchange"
      IPWRITE(UNIT_StdOut,*) "iProd =", iProd
      IPWRITE(UNIT_StdOut,*) "iReac =", iReac
      CALL abort(__STAMP__,'Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
    END IF
  ELSE
    SDEALLOCATE(PartIntEn(ReactInx(iProd))%ERot)
  END IF
END DO
!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!
IF(ProductReac(3).NE.0) THEN
  ! === 3/4 Products ============================================================================= !
  ! Determine the centre of mass velocity of the reactants
  IF(EductReac(3).NE.0) THEN
    ! Scattering 3 -> 3/4: Utilizing the FracMassCent's from above, calculated for the pseudo-molecule and the third educt,
    ! PartState(ReactInx(1)) is the centre of mass of the pseudo-molecule
    VeloCOM(1:3) = FracMassCent1 * PartState(4:6,ReactInx(1)) + FracMassCent2 * PartState(4:6,ReactInx(3))
  ELSE
    ! Scattering 2 -> 3/4
    IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) THEN
    ! Do not consider the momentum of the photon
      FracMassCent1 = 1.
      FracMassCent2 = 0.
    ELSE
      IF (usevMPF) THEN
        FracMassCent1 = Species(EductReac(1))%MassIC *Weight(1) &
            /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
        FracMassCent2 = Species(EductReac(2))%MassIC *Weight(2) &
            /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
      ELSE
        FracMassCent1 = CollInf%FracMassCent(EductReac(1),Coll_pData(iPair)%PairType)
        FracMassCent2 = CollInf%FracMassCent(EductReac(2),Coll_pData(iPair)%PairType)
      END IF
    END IF
    ! Calculation of velocity from center of mass
    VeloCOM(1:3) = FracMassCent1 * PartState(4:6,ReactInx(1)) + FracMassCent2 * PartState(4:6,ReactInx(2))
  END IF

  IF(ProductReac(4).NE.0) THEN
    ! === 4 Products ============================================================================= !
    ! FracMassCent's and reduced mass are calculated for the pseudo-molecule (1-3) and the pseudo-molecule (2-4)
    CALL CalcPseudoScatterVars_4Prod(ProductReac(1),ProductReac(3),ProductReac(2),ProductReac(4),FracMassCent1,FracMassCent2, &
          MassRed, WeightProd)

    ! Distribute the remaining collision energy
    Coll_pData(iPair)%CRela2 = 2. * Coll_pData(iPair)%Ec / MassRed

    cRelaNew(1:3) = PostCollVec(iPair)
    ! Calculate the energy value (only the relative part) to be distributed onto the products 1 and 3
    ERel_React1_React3 = 0.5 * (Species(ProductReac(1))%MassIC* WeightProd(1) + Species(ProductReac(3))%MassIC * WeightProd(3)) * DOTPRODUCT(FracMassCent2*cRelaNew(1:3))
    ! Calculate the energy value (only the relative part) to be distributed onto the products 2 and 4
    ERel_React2_React4 = 0.5 * (Species(ProductReac(2))%MassIC* WeightProd(2) + Species(ProductReac(4))%MassIC * WeightProd(4)) * DOTPRODUCT(FracMassCent1*cRelaNew(1:3))
    ! Distribute the energy onto the components of the pseudo molecule (2-4)
    IF (usevMPF) THEN
      FracMassCent1 = Species(ProductReac(2))%MassIC * WeightProd(2) &
          /(Species(ProductReac(2))%MassIC* WeightProd(2) + Species(ProductReac(4))%MassIC * WeightProd(4))
      FracMassCent2 = Species(ProductReac(4))%MassIC * WeightProd(4) &
          /(Species(ProductReac(2))%MassIC* WeightProd(2) + Species(ProductReac(4))%MassIC * WeightProd(4))
      MassRed = Species(ProductReac(2))%MassIC *WeightProd(2)* Species(ProductReac(4))%MassIC *WeightProd(4) &
          / (Species(ProductReac(2))%MassIC*WeightProd(2) + Species(ProductReac(4))%MassIC *WeightProd(4))
    ELSE
      FracMassCent1 = CollInf%FracMassCent(ProductReac(2),CollInf%Coll_Case(ProductReac(2),ProductReac(4)))
      FracMassCent2 = CollInf%FracMassCent(ProductReac(4),CollInf%Coll_Case(ProductReac(2),ProductReac(4)))
      MassRed = CollInf%MassRed(CollInf%Coll_Case(ProductReac(2),ProductReac(4)))
    END IF

    Coll_pData(iPair)%CRela2 = 2. * ERel_React2_React4 / MassRed

    ! Set the energy for the 2-4 pair
    cRelaNew(1:3) = PostCollVec(iPair)
    ! deltaV particle 2
    PartState(4:6,ReactInx(2)) = VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)
    ! deltaV particle 4
    PartState(4:6,ReactInx(4)) = VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)
#ifdef CODE_ANALYZE
    Energy_new=0.5*Species(ProductReac(2))%MassIC*DOTPRODUCT(VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)) * WeightProd(2)
    Momentum_new(1:3) = Species(ProductReac(2))%MassIC* (VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)) * WeightProd(2)
    Energy_new = Energy_new &
                  + 0.5*Species(ProductReac(4))%MassIC*DOTPRODUCT(VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)) * WeightProd(4)
    Momentum_new(1:3) = Momentum_new(1:3) &
                        + Species(ProductReac(4))%MassIC*(VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)) * WeightProd(4)
    IF((Species(ProductReac(2))%InterID.EQ.2).OR.(Species(ProductReac(2))%InterID.EQ.20)) &
      Energy_new = Energy_new + (PartIntEn(ReactInx(2))%EVib(1) + PartIntEn(ReactInx(2))%ERot(1)) * WeightProd(2)
    IF((Species(ProductReac(4))%InterID.EQ.2).OR.(Species(ProductReac(4))%InterID.EQ.20)) &
      Energy_new = Energy_new + (PartIntEn(ReactInx(4))%EVib(1) + PartIntEn(ReactInx(4))%ERot(1)) * WeightProd(4)
    IF(DSMC%ElectronicModel.GT.0) THEN
      IF((Species(ProductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(2))%FullyIonized)) &
        Energy_new=Energy_new + PartIntEn(ReactInx(2))%EElec(1)*WeightProd(2)
      IF((Species(ProductReac(4))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(4))%FullyIonized)) &
        Energy_new=Energy_new + PartIntEn(ReactInx(4))%EElec(1)*WeightProd(4)
    END IF
#endif /* CODE_ANALYZE */
  ELSE
    ! === 3 Products ============================================================================= !
    ! The remaining collision energy has to distributed onto three particles (pseudo molecule requires more energy due to additional
    ! translational degrees of freedom)
    CALL RANDOM_NUMBER(iRan)
    FacEtraDistri = iRan
    CALL RANDOM_NUMBER(iRan)
    DO WHILE ((4 *FacEtraDistri*(1-FacEtraDistri))**(1-CollInf%omega(ProductReac(1),ProductReac(2))).LT.iRan)
      CALL RANDOM_NUMBER(iRan)
      FacEtraDistri = iRan
      CALL RANDOM_NUMBER(iRan)
    END DO
    ERel_React1_React2 = Coll_pData(iPair)%Ec * FacEtraDistri
    ERel_React1_React3 = Coll_pData(iPair)%Ec - ERel_React1_React2
    ! FracMassCent's and reduced mass are calculated for the pseudo-molecule 1-3 and the second product, in the case of dissociation
    ! this is the non-reacting collision partner
    CALL CalcPseudoScatterVars(ProductReac(1),ProductReac(3),ProductReac(2),FracMassCent1,FracMassCent2,MassRed &
          , (/WeightProd(1),WeightProd(3),WeightProd(2)/))

    Coll_pData(iPair)%cRela2 = 2 * ERel_React1_React2 / MassRed

    cRelaNew(1:3) = PostCollVec(iPair)
    PartState(4:6,ReactInx(2)) = VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)
#ifdef CODE_ANALYZE
    Energy_new=0.5*Species(ProductReac(2))%MassIC*DOTPRODUCT(VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3))* WeightProd(2)
    Momentum_new(1:3) = Species(ProductReac(2))%MassIC* (VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)) * WeightProd(2)

    IF((Species(ProductReac(2))%InterID.EQ.2).OR.(Species(ProductReac(2))%InterID.EQ.20)) &
      Energy_new = Energy_new + (PartIntEn(ReactInx(2))%EVib(1) + PartIntEn(ReactInx(2))%ERot(1)) * WeightProd(2)
    IF(DSMC%ElectronicModel.GT.0) THEN
      IF((Species(ProductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(2))%FullyIonized)) &
        Energy_new=Energy_new + PartIntEn(ReactInx(2))%EElec(1)*WeightProd(2)
    END IF
#endif /* CODE_ANALYZE */
    ! Set velocity of pseudo molec (AB) and calculate the centre of mass frame velocity: m_pseu / (m_3 + m_4) * v_pseu
    ! (Velocity of pseudo molecule is NOT equal to the COM frame velocity)
    VeloCOM(1:3) = (VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3))
  END IF
  ! === 3/4 Products ============================================================================= !
  ! Treatment of the components of the pseudo-molecule (1-3)
  IF (usevMPF) THEN
    FracMassCent1 = Species(ProductReac(1))%MassIC *WeightProd(1) &
        /(Species(ProductReac(1))%MassIC* WeightProd(1) + Species(ProductReac(3))%MassIC * WeightProd(3))
    FracMassCent2 = Species(ProductReac(3))%MassIC *WeightProd(3) &
        /(Species(ProductReac(1))%MassIC* WeightProd(1) + Species(ProductReac(3))%MassIC * WeightProd(3))
    MassRed = Species(ProductReac(1))%MassIC *WeightProd(1)* Species(ProductReac(3))%MassIC *WeightProd(3) &
        / (Species(ProductReac(1))%MassIC*WeightProd(1) + Species(ProductReac(3))%MassIC *WeightProd(3))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(ProductReac(1),CollInf%Coll_Case(ProductReac(1),ProductReac(3)))
    FracMassCent2 = CollInf%FracMassCent(ProductReac(3),CollInf%Coll_Case(ProductReac(1),ProductReac(3)))
    MassRed = CollInf%MassRed(CollInf%Coll_Case(ProductReac(1),ProductReac(3)))
  END IF

  Coll_pData(iPair)%cRela2 = 2 * ERel_React1_React3 / MassRed

  cRelaNew(1:3) = PostCollVec(iPair)

  !deltaV particle 1
  PartState(4:6,ReactInx(1)) = VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)

  !deltaV particle 3
  PartState(4:6,ReactInx(3)) = VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)

#ifdef CODE_ANALYZE
  ! New total energy
  Energy_new=Energy_new + 0.5*Species(ProductReac(1))%MassIC*DOTPRODUCT(VeloCOM(1:3)+FracMassCent2*cRelaNew(1:3))*WeightProd(1) &
                        + 0.5*Species(ProductReac(3))%MassIC*DOTPRODUCT(VeloCOM(1:3)-FracMassCent1*cRelaNew(1:3))*WeightProd(3)
  IF((Species(ProductReac(1))%InterID.EQ.2).OR.(Species(ProductReac(1))%InterID.EQ.20)) &
    Energy_new = Energy_new + (PartIntEn(ReactInx(1))%EVib(1) + PartIntEn(ReactInx(1))%ERot(1)) * WeightProd(1)
  IF((Species(ProductReac(3))%InterID.EQ.2).OR.(Species(ProductReac(3))%InterID.EQ.20)) &
    Energy_new = Energy_new + (PartIntEn(ReactInx(3))%EVib(1) + PartIntEn(ReactInx(3))%ERot(1)) * WeightProd(3)
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(ProductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(1))%FullyIonized)) &
      Energy_new=Energy_new + PartIntEn(ReactInx(1))%EElec(1)*WeightProd(1)
    IF((Species(ProductReac(3))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(3))%FullyIonized)) &
      Energy_new=Energy_new + PartIntEn(ReactInx(3))%EElec(1)*WeightProd(3)
  END IF
  ! New total momentum
  Momentum_new(1:3) = Momentum_new(1:3) &
                    + Species(ProductReac(1))%MassIC * (VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)) * WeightProd(1) &
                    + Species(ProductReac(3))%MassIC * (VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)) * WeightProd(3)
#endif /* CODE_ANALYZE */

ELSEIF(ProductReac(3).EQ.0) THEN
  ! === 2 Products =============================================================================== !
  IF(EductReac(3).NE.0) THEN
    ! Scattering 3 -> 2
    VeloCOM(1:3) = FracMassCent1 * PartState(4:6,ReactInx(1)) + FracMassCent2 * PartState(4:6,ReactInx(3))
    ! When RHS is set, ReactInx(2) is utilized, not an error as the old state cancels out after the particle push in the time disc,
    ! therefore, there is no need to set change the index as the proper species, ProductReac(2), was utilized for the relaxation
    CALL RemoveParticle(ReactInx(3))
  ELSE
    ! Scattering 2 -> 2
    IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) THEN
    ! Do not consider the momentum of the photon
      FracMassCent1 = 1.
      FracMassCent2 = 0.
    ELSE
      IF (usevMPF) THEN
        FracMassCent1 = Species(EductReac(1))%MassIC *Weight(1) &
            /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
        FracMassCent2 = Species(EductReac(2))%MassIC *Weight(2) &
            /(Species(EductReac(1))%MassIC* Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
      ELSE
        FracMassCent1 = CollInf%FracMassCent(EductReac(1),CollInf%Coll_Case(EductReac(1),EductReac(2)))
        FracMassCent2 = CollInf%FracMassCent(EductReac(2),CollInf%Coll_Case(EductReac(1),EductReac(2)))
      END IF
    END IF
    ! Centre of mass velocity
    VeloCOM(1:3) = FracMassCent1 * PartState(4:6,ReactInx(1)) + FracMassCent2 * PartState(4:6,ReactInx(2))
  END IF
  ERel_React1_React3 = Coll_pData(iPair)%Ec

  ! Sanity check
  IF(ERel_React1_React3.LT.0.) THEN
    IPWRITE(UNIT_StdOut,*) "iReac =", iReac
    CALL abort(__STAMP__,'DSMC Chemistry: Coll_pData(iPair)%Ec < 0.)',RealInfoOpt = Coll_pData(iPair)%Ec)
  END IF

  IF (usevMPF) THEN
    FracMassCent1 = Species(ProductReac(1))%MassIC *WeightProd(1) &
        /(Species(ProductReac(1))%MassIC* WeightProd(1) + Species(ProductReac(2))%MassIC * WeightProd(2))
    FracMassCent2 = Species(ProductReac(2))%MassIC *WeightProd(2) &
        /(Species(ProductReac(1))%MassIC* WeightProd(1) + Species(ProductReac(2))%MassIC * WeightProd(2))
    MassRed = Species(ProductReac(1))%MassIC *WeightProd(1)* Species(ProductReac(2))%MassIC *WeightProd(2) &
        / (Species(ProductReac(1))%MassIC*WeightProd(1) + Species(ProductReac(2))%MassIC *WeightProd(2))
  ELSE
    ! Scattering of (AB)
    FracMassCent1 = CollInf%FracMassCent(ProductReac(1),CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
    FracMassCent2 = CollInf%FracMassCent(ProductReac(2),CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
    MassRed = CollInf%MassRed(CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
  END IF

  Coll_pData(iPair)%cRela2 = 2 * ERel_React1_React3 / MassRed

  cRelaNew(1:3) = PostCollVec(iPair)

  !deltaV particle 1
  PartState(4:6,ReactInx(1)) = VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)
  !deltaV particle 2
  PartState(4:6,ReactInx(2)) = VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)

#ifdef CODE_ANALYZE
  ! New total energy of remaining products (here, recombination: 2 products)
  Energy_new = 0.5*Species(ProductReac(1))%MassIC*((VeloCOM(1) + FracMassCent2*cRelaNew(1))**2  &
                                                 + (VeloCOM(2) + FracMassCent2*cRelaNew(2))**2  &
                                                 + (VeloCOM(3) + FracMassCent2*cRelaNew(3))**2) * WeightProd(1) &
              +0.5*Species(ProductReac(2))%MassIC*((VeloCOM(1) - FracMassCent1*cRelaNew(1))**2  &
                                                 + (VeloCOM(2) - FracMassCent1*cRelaNew(2))**2  &
                                                 + (VeloCOM(3) - FracMassCent1*cRelaNew(3))**2) * WeightProd(2)
  IF((Species(ProductReac(1))%InterID.EQ.2).OR.(Species(ProductReac(1))%InterID.EQ.20)) &
    Energy_new = Energy_new + (PartIntEn(ReactInx(1))%EVib(1) + PartIntEn(ReactInx(1))%ERot(1)) * WeightProd(1)
  IF((Species(ProductReac(2))%InterID.EQ.2).OR.(Species(ProductReac(2))%InterID.EQ.20)) &
    Energy_new = Energy_new + (PartIntEn(ReactInx(2))%EVib(1) + PartIntEn(ReactInx(2))%ERot(1)) * WeightProd(2)
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(ProductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(1))%FullyIonized)) &
      Energy_new=Energy_new + PartIntEn(ReactInx(1))%EElec(1)*WeightProd(1)
    IF((Species(ProductReac(2))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(2))%FullyIonized)) &
      Energy_new=Energy_new + PartIntEn(ReactInx(2))%EElec(1)*WeightProd(2)
  END IF
  ! New total momentum
    Momentum_new(1:3) = Species(ProductReac(1))%MassIC * (/VeloCOM(1) + FracMassCent2*cRelaNew(1),  &
                                                           VeloCOM(2) + FracMassCent2*cRelaNew(2),  &
                                                           VeloCOM(3) + FracMassCent2*cRelaNew(3)/) * WeightProd(1) &
                      + Species(ProductReac(2))%MassIC * (/VeloCOM(1) - FracMassCent1*cRelaNew(1),  &
                                                           VeloCOM(2) - FracMassCent1*cRelaNew(2),  &
                                                           VeloCOM(3) - FracMassCent1*cRelaNew(3)/) * WeightProd(2)
#endif /* CODE_ANALYZE */
END IF

IF(ChemReac%NumDeleteProducts.GT.0) THEN
  DO iSpec = 1, ChemReac%NumDeleteProducts
    SpecToDelete = ChemReac%DeleteProductsList(iSpec)
    DO iProd = 1, NumProd
      IF(ProductReac(iProd).EQ.SpecToDelete) THEN
        ! Remove the respective particle
        CALL RemoveParticle(ReactInx(iProd))
      END IF
    END DO
  END DO
END IF

IF(CalcPartBalance) THEN
  DO iPart = 1, NumProd
    nPartIn(ProductReac(iPart))=nPartIn(ProductReac(iPart)) + 1
    PartEkinIn(ProductReac(iPart))=PartEkinIn(ProductReac(iPart))+CalcEkinPart(ReactInx(iPart))
  END DO
END IF

! Species-specific time step: check whether the product species with the lower time step is retained or deleted and the velocity
! vector of the reactant restored
IF(VarTimeStep%UseSpeciesSpecific) THEN
  IF(VarTimeStep%DisableForMCC) THEN
    DO iProd = 3, NumProd
      iSpec = ProductReac(iProd)
      IF(Species(iSpec)%TimeStepFactor.LT.1.) THEN
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.GT.Species(iSpec)%TimeStepFactor) THEN
          CALL RemoveParticle(ReactInx(iProd))
          IF(Species(PartSpecies(ReactInx(1)))%TimeStepFactor.LT.1.) THEN
            PartState(4:6,ReactInx(1)) = TempVelo(1:3)
          ELSE IF(Species(PartSpecies(ReactInx(2)))%TimeStepFactor.LT.1.) THEN
            PartState(4:6,ReactInx(2)) = TempVelo(1:3)
          END IF
        END IF
      END IF
    END DO
  END IF
END IF

#ifdef CODE_ANALYZE
! Check for energy difference
IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,RelEneTol)) THEN
  WRITE(UNIT_StdOut,*) '\n'
  IPWRITE(UNIT_StdOut,*)    " Species, Reactants     : ",EductReac
  IPWRITE(UNIT_StdOut,*)    " Species, Products      : ",ProductReac
  IPWRITE(UNIT_StdOut,*)    " Weights                : ",Weight
  IPWRITE(UNIT_StdOut,*)    " Weights, Products      : ",WeightProd
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_old-Energy_new
  ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
    IF(energy.GT.0.0)THEN
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_old-Energy_new)/energy
    END IF
  END ASSOCIATE
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelEneTol
  CALL abort(&
      __STAMP__&
      ,'CODE_ANALYZE: DSMC_Chemistry is not energy conserving for chemical reaction:', IntInfoOpt=iReac)
END IF
! Check for momentum difference
IF(Symmetry%Order.EQ.3) THEN
  ! Do not check the momentum in z as it can be very small (close to machine precision), leading to greater relative errors
  iMomDim = 3
ELSE IF(Symmetry%Order.EQ.2) THEN
  iMomDim = 2
ELSE
  iMomDim = 1
END IF
DO iMom=1,iMomDim
  IF (.NOT.ALMOSTEQUALRELATIVE(Momentum_old(iMom),Momentum_new(iMom),RelMomTol)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Direction (x,y,z)        : ",iMom
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_old             : ",Momentum_old(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_new             : ",Momentum_new(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Momentum difference : ",Momentum_old(iMom)-Momentum_new(iMom)
    ASSOCIATE( Momentum => MAX(ABS(Momentum_old(iMom)),ABS(Momentum_new(iMom))) )
      IF(Momentum.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Momentum difference : ",(Momentum_old(iMom)-Momentum_new(iMom))/Momentum
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelMomTol
    CALL abort(&
        __STAMP__&
        ,'CODE_ANALYZE: DSMC_Chemistry is not momentum conserving for chemical reaction:', IntInfoOpt=iReac)
  END IF
END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE DSMC_Chemistry


!===================================================================================================================================
!> Calculation of the backward reaction rate with partition sums, interpolation within the given temperature interval
!===================================================================================================================================
SUBROUTINE CalcBackwardRate(iReacTmp,LocalTemp,BackwardRate)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: maxEXP
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, ChemReac, QKChemistry
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_DSMC_QK_Chemistry     ,ONLY: QK_CalcAnalyticRate
USE MOD_Particle_Analyze_Tools,ONLY: CalcPartitionFunction
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
INTEGER                         :: iReac, iSpec, LowerLevel, UpperLevel, iChemDir, MaxElecQua
REAL                            :: PartFuncProduct(2), k_b_lower, k_b_upper, ActivationEnergy_K, PartitionFunction
REAL                            :: Qtra, Qrot, Qvib, Qelec, expVal
!===================================================================================================================================
  ! Determination of the lower and upper value of the temperature interval
  LowerLevel = INT(LocalTemp/DSMC%PartitionInterval)
  UpperLevel = LowerLevel + 1

  ! Reading the stoichiometric coefficients from the reactants
  iReac = ChemReac%BackwardReacForwardIndx(iReacTmp)
  IF (TRIM(ChemReac%ReactModel(iReac)).EQ.'QK') THEN
    IF (TRIM(ChemReac%ReactType(iReac)).EQ.'I') THEN
      MaxElecQua=SpecDSMC(ChemReac%Reactants(iReac,1))%MaxElecQuant - 1
      ActivationEnergy_K = SpecDSMC(ChemReac%Reactants(iReac,1))%ElectronicState(2,MaxElecQua)
    ELSEIF(TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
      ActivationEnergy_K = SpecDSMC(ChemReac%Reactants(iReac,1))%Ediss_eV * 11604.52500617 ! eV -> K
    END IF
  END IF

  ! Calculation of the backward reaction rate using the equilibrium constant)
  IF((UpperLevel.GT.INT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)).OR.(LowerLevel.EQ.0)) THEN
  ! Direct calculation at given temperature
    PartFuncProduct(1:2) = 1.
    DO iSpec = 1, nSpecies
      DO iChemDir = 1,2
        IF(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir).NE.0) THEN
          CALL CalcPartitionFunction(iSpec, LocalTemp, Qtra, Qrot, Qvib, Qelec)
          PartitionFunction = Qtra * Qrot * Qvib * Qelec
          PartFuncProduct(iChemDir) = PartFuncProduct(iChemDir)   &
            * PartitionFunction**(ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir))
        END IF
      END DO
    END DO
    IF (TRIM(ChemReac%ReactModel(iReac)).EQ.'QK') THEN
      expVal = MIN(maxexp,ActivationEnergy_K/LocalTemp)
      BackwardRate = QK_CalcAnalyticRate(iReac,LocalTemp)*(PartFuncProduct(1)/PartFuncProduct(2))*EXP(expVal)
    ELSE
      BackwardRate = ChemReac%Arrhenius_Prefactor(iReac)  &
              * (LocalTemp)**ChemReac%Arrhenius_Powerfactor(iReac) &
              * (PartFuncProduct(1)/PartFuncProduct(2))
    END IF
  ELSE
  ! Interpolation between tabulated lower and upper values
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
      IF (TRIM(ChemReac%ReactModel(iReac)).EQ.'QK') THEN
        expVal = MIN(maxexp,ActivationEnergy_K/(LowerLevel * DSMC%PartitionInterval))
        k_b_lower = QKChemistry(iReac)%ForwardRate(LowerLevel)* (PartFuncProduct(1)/PartFuncProduct(2))* EXP(expVal)
      ELSE
        k_b_lower = ChemReac%Arrhenius_Prefactor(iReac)  &
                * (LowerLevel * DSMC%PartitionInterval)**ChemReac%Arrhenius_Powerfactor(iReac) &
                * (PartFuncProduct(1)/PartFuncProduct(2))
      END IF
    ELSE
      k_b_lower = 0.0
    END IF

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
      IF (TRIM(ChemReac%ReactModel(iReac)).EQ.'QK') THEN
        expVal = MIN(maxexp,ActivationEnergy_K/(UpperLevel * DSMC%PartitionInterval))
        k_b_upper = QKChemistry(iReac)%ForwardRate(UpperLevel)* (PartFuncProduct(1)/PartFuncProduct(2)) * EXP(expVal)
      ELSE
        k_b_upper = ChemReac%Arrhenius_Prefactor(iReac) &
              * (UpperLevel * DSMC%PartitionInterval)**ChemReac%Arrhenius_Powerfactor(iReac) &
              * (PartFuncProduct(1)/PartFuncProduct(2))
      END IF
    ELSE
      k_b_upper = 0.0
    END IF
  ! Linear interpolation of the backward rate coefficient at the actual temperature
    BackwardRate = k_b_lower &
              + (k_b_upper - k_b_lower)  &
              / (DSMC%PartitionInterval) * (LocalTemp - LowerLevel * DSMC%PartitionInterval)
  END IF

END SUBROUTINE CalcBackwardRate


!===================================================================================================================================
!> Routine determines the reduced mass and the mass fraction between a pseudo-molecule and third species
!===================================================================================================================================
SUBROUTINE CalcPseudoScatterVars(PseuSpec1, PseuSpec2, ScatterSpec3, FracMassCent1, FracMassCent2, MassRed, Weight)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: PseuSpec1, PseuSpec2, ScatterSpec3
REAL, INTENT(IN)              :: Weight(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: FracMassCent1, FracMassCent2, MassRed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Mass
!===================================================================================================================================
Mass = Species(PseuSpec1)%MassIC*Weight(1) + Species(PseuSpec2)%MassIC*Weight(2)
FracMassCent1 = Mass / (Mass + Species(ScatterSpec3)%MassIC*Weight(3))
FracMassCent2 = Species(ScatterSpec3)%MassIC*Weight(3) / (Mass + Species(ScatterSpec3)%MassIC*Weight(3))
MassRed = (Mass*Species(ScatterSpec3)%MassIC*Weight(3)) &
                        / (Mass+Species(ScatterSpec3)%MassIC*Weight(3))
END SUBROUTINE CalcPseudoScatterVars


!===================================================================================================================================
!> Routine determines the reduced mass and the mass fraction between a pseudo-molecule and third species
!===================================================================================================================================
SUBROUTINE CalcPseudoScatterVars_4Prod(Spec1, Spec2, Spec3, Spec4, FracMassCent1, FracMassCent2, MassRed, Weight)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: Spec1, Spec2, Spec3, Spec4
REAL, INTENT(IN)              :: Weight(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: FracMassCent1, FracMassCent2, MassRed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: MassPseu1, MassPseu2
!===================================================================================================================================
MassPseu1 = Species(Spec1)%MassIC*Weight(1) + Species(Spec2)%MassIC*Weight(2)
MassPseu2 = Species(Spec3)%MassIC*Weight(3) + Species(Spec4)%MassIC*Weight(4)
FracMassCent1 = MassPseu1 / (MassPseu1 + MassPseu2)
FracMassCent2 = MassPseu2 / (MassPseu1 + MassPseu2)
MassRed = (MassPseu1*MassPseu2) / (MassPseu1+MassPseu2)
END SUBROUTINE CalcPseudoScatterVars_4Prod

END MODULE MOD_DSMC_ChemReact