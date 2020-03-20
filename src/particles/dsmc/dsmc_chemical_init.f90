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

MODULE MOD_DSMC_ChemInit
!===================================================================================================================================
! Initialization of chemical module
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_chemical_init
  MODULE PROCEDURE DSMC_chemical_init
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_chemical_init
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_chemical_init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY: ChemReac,CollisMode, DSMC, QKAnalytic, SpecDSMC, BGGas
  USE MOD_ReadInTools
  USE MOD_Globals
  USE MOD_Globals_Vars,           ONLY: BoltzmannConst
  USE MOD_PARTICLE_Vars,          ONLY: nSpecies
  USE MOD_Particle_Analyze_Vars,  ONLY: ChemEnergySum
  USE MOD_DSMC_ChemReact,         ONLY: CalcPartitionFunction, CalcQKAnalyticRate
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(LEN=3)         :: hilf
  INTEGER               :: iReac, iReac2, iReac3, iReac4, iSpec, iChemDir, iReacForward, iReacDiss, PartitionArraySize, iInter
  INTEGER, ALLOCATABLE  :: PairCombID(:,:), DummyRecomb(:,:)
  LOGICAL, ALLOCATABLE  :: YetDefined_Help(:)
  LOGICAL               :: DoScat
  INTEGER               :: Reactant1, Reactant2, Reactant3, MaxSpecies, MaxElecQua, ReadInNumOfReact
  REAL                  :: Temp, Qtra, Qrot, Qvib, Qelec, BGGasEVib
!===================================================================================================================================

! reading reaction values
  ChemReac%NumOfReact = GETINT('DSMC-NumOfReactions','0')
  ReadInNumOfReact = ChemReac%NumOfReact
  IF(CollisMode.EQ.3)THEN
    IF(ChemReac%NumOfReact.EQ.0)THEN
      CALL Abort(&
__STAMP__&
,' Collisions with chemical reactions require chemical reaction database! ',CollisMode,REAL(ChemReac%NumOfReact))
    END IF
  END IF
  ALLOCATE(ChemReac%ArbDiss(ChemReac%NumOfReact))
  ! Allowing unspecified non-reactive collision partner (CH4 + M -> CH3 + H + M, e.g. (/1,0,0/) -> (/2,0,3/)
  iReacDiss = ChemReac%NumOfReact
  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    ChemReac%ArbDiss(iReac)%NumOfNonReactives = GETINT('DSMC-Reaction'//TRIM(hilf)//'-NumberOfNonReactives','0')
    IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.GT.0) THEN
      ALLOCATE(ChemReac%ArbDiss(iReac)%NonReactiveSpecies(ChemReac%ArbDiss(iReac)%NumOfNonReactives))
      ChemReac%ArbDiss(iReac)%NonReactiveSpecies = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-NonReactiveSpecies', &
                                                ChemReac%ArbDiss(iReac)%NumOfNonReactives)
      ! First reaction is saved within the dummy input reaction, thus "- 1"
      ChemReac%NumOfReact = ChemReac%NumOfReact + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
    END IF
  END DO
  ! Delete products if they belong to a certain species
  ChemReac%NumDeleteProducts = GETINT('Particles-Chemistry-NumDeleteProducts')
  IF(ChemReac%NumDeleteProducts.GT.0) THEN
    ALLOCATE(ChemReac%DeleteProductsList(ChemReac%NumDeleteProducts))
    ChemReac%DeleteProductsList = GETINTARRAY('Particles-Chemistry-DeleteProductsList', ChemReac%NumDeleteProducts)
  END IF
  ! Calculation of the backward reaction rates
  IF(DSMC%BackwardReacRate) THEN
   ChemReac%NumOfReact = 2 * ChemReac%NumOfReact
  END IF
  ChemEnergySum = 0.
  SWRITE(*,*) '| Number of considered reaction paths (including dissociation and recombination): ', ChemReac%NumOfReact
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (ChemReac%NumOfReact.GT.0) THEN
    ALLOCATE(ChemReac%NumReac(ChemReac%NumOfReact))
    ChemReac%NumReac = 0
#if (PP_TimeDiscMethod==42)
    ALLOCATE(ChemReac%ReacCount(ChemReac%NumOfReact))
    ChemReac%ReacCount = 0
    ALLOCATE(ChemReac%ReacCollMean(ChemReac%NumOfReact))
    ChemReac%ReacCollMean = 0.0
    ALLOCATE(ChemReac%ReacCollMeanCount(ChemReac%NumOfReact))
    ChemReac%ReacCollMeanCount = 0
#endif
    ALLOCATE(ChemReac%QKProcedure(ChemReac%NumOfReact))
    ChemReac%QKProcedure = .FALSE.
    ALLOCATE(ChemReac%QKMethod(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%QKCoeff(ChemReac%NumOfReact,2))
    ALLOCATE(YetDefined_Help(ChemReac%NumOfReact))
    YetDefined_Help = .FALSE.
    ALLOCATE(ChemReac%ReactType(ChemReac%NumOfReact))
    ChemReac%ReactType = '0'
    ALLOCATE(ChemReac%DefinedReact(ChemReac%NumOfReact,2,3))
    ChemReac%DefinedReact = 0
    ALLOCATE(ChemReac%ReactCase(nSpecies,nSpecies))
    ChemReac%ReactCase = 0
    ALLOCATE(ChemReac%ReactNum(nSpecies, nSpecies, nSpecies))
    ChemReac%ReactNum = 0
    ALLOCATE(ChemReac%ReactNumRecomb(nSpecies, nSpecies, nSpecies))
    ChemReac%ReactNumRecomb = 0
    ALLOCATE(ChemReac%Arrhenius_Prefactor(ChemReac%NumOfReact),&
             ChemReac%Arrhenius_Powerfactor(ChemReac%NumOfReact),&
             ChemReac%EActiv(ChemReac%NumOfReact),&
             ChemReac%EForm(ChemReac%NumOfReact),&
             ChemReac%Hab(ChemReac%NumOfReact))
    ChemReac%Hab=0.0
    ALLOCATE(ChemReac%MeanEVibQua_PerIter(nSpecies))
    ChemReac%MeanEVibQua_PerIter = 0
    ALLOCATE(ChemReac%MeanEVib_PerIter(nSpecies))
    ChemReac%MeanEVib_PerIter = 0.0
    ALLOCATE(ChemReac%MeanXiVib_PerIter(nSpecies))
    ChemReac%MeanXiVib_PerIter = 0.0
    ALLOCATE(DummyRecomb(nSpecies,nSpecies))
    DummyRecomb = 0
    ALLOCATE(ChemReac%CEXa(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%CEXb(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MEXa(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MEXb(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%ELa(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%ELb(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%DoScat(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%ReactInfo(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%TLU_FileName(ChemReac%NumOfReact))

    IF (BGGas%NumberOfSpecies.GT.0) THEN
      DO iSpec = 1, nSpecies
        IF(BGGas%BackgroundSpecies(iSpec)) THEN
          ! Background gas: Calculation of the mean vibrational quantum number of diatomic molecules
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            IF(.NOT.SpecDSMC(iSpec)%PolyatomicMol) THEN
              BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib &
                + BoltzmannConst * SpecDSMC(iSpec)%CharaTVib / (EXP(SpecDSMC(iSpec)%CharaTVib / SpecDSMC(iSpec)%Init(0)%TVib) - 1)
              BGGasEVib = BGGasEVib/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
              ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(BGGasEVib) + 1, SpecDSMC(iSpec)%MaxVibQuant)
              ChemReac%MeanXiVib_PerIter(iSpec) = 2. * ChemReac%MeanEVibQua_PerIter(iSpec) &
                                                * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0 )
            END IF
          ELSE
            ChemReac%MeanEVibQua_PerIter(iSpec) = 0
            ChemReac%MeanXiVib_PerIter(iSpec) = 0.
          END IF
        END IF
      END DO
    END IF

    DoScat = .false.
    DO iReac = 1, ReadInNumOfReact
      WRITE(UNIT=hilf,FMT='(I0)') iReac
      ChemReac%ReactType(iReac)             = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-ReactionType','0'))
      ChemReac%QKProcedure(iReac)           = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-QKProcedure','.FALSE.')
      CHemReac%QKMethod(iReac)       = GETINT('DSMC-Reaction'//TRIM(hilf)//'-QK-Method','0')
      ChemReac%QKCoeff(iReac,1)      = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-QK-Coeff1','0')
      ChemReac%QKCoeff(iReac,2)      = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-QK-Coeff2','0')
      ChemReac%DefinedReact(iReac,1,:)      = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Reactants',3,'0,0,0')
      ChemReac%DefinedReact(iReac,2,:)      = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Products',3,'0,0,0')
      ChemReac%Arrhenius_Prefactor(iReac)   = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Prefactor','0')
      ChemReac%Arrhenius_Powerfactor(iReac) = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Powerfactor','0')
      ChemReac%EActiv(iReac)                = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Activation-Energy_K','0')*BoltzmannConst
      ChemReac%CEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXa','-27.2')
      ChemReac%CEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXb','175.269')
      ChemReac%DoScat(iReac)               = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-DoScat','.FALSE.')
      IF (ChemReac%DoScat(iReac)) THEN
        DoScat = .true.
        ChemReac%ELa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-ELa','-26.8')        ! Passen diese Zeilen so?
        ChemReac%ELb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-ELb','148.975')
        ChemReac%TLU_FileName(iReac)        = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-TLU_FileName','0'))
        IF(TRIM(ChemReac%TLU_FileName(iReac)).EQ.'0') THEN
          CALL abort(__STAMP__,&
            'Input Error: No file containing table lookup data for scattering simulation has been found.')
        END IF
      ELSE
        ChemReac%MEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXa','-27.2')
        ChemReac%MEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXb','175.269')
      END IF
      ! ChemReac%MEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXa','0')
      ! ChemReac%MEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXb','0')

      ! Filling up ChemReac-Array for the given non-reactive dissociation/electron-impact ionization partners
      IF((TRIM(ChemReac%ReactType(iReac)).EQ.'D').OR.(TRIM(ChemReac%ReactType(iReac)).EQ.'iQK')) THEN
        IF((ChemReac%DefinedReact(iReac,1,2).EQ.0).AND.(ChemReac%DefinedReact(iReac,2,2).EQ.0)) THEN
          IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.EQ.0) THEN
            CALL abort(__STAMP__,&
            'Error in Definition: Non-reacting partner(s) has to be defined!',iReac)
          END IF
          DO iReac2 = 1, ChemReac%ArbDiss(iReac)%NumOfNonReactives
            IF(iReac2.EQ.1) THEN
              ! The first non-reacting partner is written into the reaction, which was read-in.
              ChemReac%DefinedReact(iReac,1,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
              ChemReac%DefinedReact(iReac,2,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
            ELSE
              ! The following reaction are added after the number of originally read-in reactions (counter: iReacDiss)
              ChemReac%ReactType(iReacDiss+iReac2-1)             = ChemReac%ReactType(iReac)
              ChemReac%QKProcedure(iReacDiss+iReac2-1)           = ChemReac%QKProcedure(iReac)
              ChemReac%DefinedReact(iReacDiss+iReac2-1,1,:)      = ChemReac%DefinedReact(iReac,1,:)
              ChemReac%DefinedReact(iReacDiss+iReac2-1,1,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
              ChemReac%DefinedReact(iReacDiss+iReac2-1,2,:)      = ChemReac%DefinedReact(iReac,2,:)
              ChemReac%DefinedReact(iReacDiss+iReac2-1,2,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
              ChemReac%Arrhenius_Prefactor(iReacDiss+iReac2-1)   = ChemReac%Arrhenius_Prefactor(iReac)
              ChemReac%Arrhenius_Powerfactor(iReacDiss+iReac2-1) = ChemReac%Arrhenius_Powerfactor(iReac)
              ChemReac%EActiv(iReacDiss+iReac2-1)                = ChemReac%EActiv(iReac)
            END IF
          END DO
          iReacDiss = iReacDiss + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
        ELSE IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.NE.0) THEN
          CALL abort(__STAMP__,&
            'Dissociation - Error in Definition: Non-reacting partner(s) has to be zero!',iReac)
        END IF
      END IF
    END DO

    IF (DoScat) CALL Init_TLU_Data()

    ! Calculation of stoichiometric coefficients and calculation of the heat of formation
    DO iReac = 1, ChemReac%NumOfReact
      ALLOCATE(ChemReac%ReactInfo(iReac)%StoichCoeff(nSpecies,2))
      ChemReac%ReactInfo(iReac)%StoichCoeff(1:nSpecies,1:2) = 0
      ChemReac%EForm(iReac) = 0.0
      DO iSpec=1, nSpecies
        DO iChemDir=1,2
          IF(ChemReac%DefinedReact(iReac,iChemDir,1).EQ.iSpec) THEN
            ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
          END IF
          IF(ChemReac%DefinedReact(iReac,iChemDir,2).EQ.iSpec) THEN
            ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
          END IF
          IF(ChemReac%DefinedReact(iReac,iChemDir,3).EQ.iSpec) THEN
            ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
          END IF
        END DO
        ! Calculation of the enthalpy of reaction by the summation of the enthalpies of formation of the respective species
        ! (ionization energy of ionized species was already added in dsmc_init.f90)
        ChemReac%EForm(iReac) = ChemReac%EForm(iReac) &
                              - ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,2)*SpecDSMC(iSpec)%HeatOfFormation  &
                              + ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,1)*SpecDSMC(iSpec)%HeatOfFormation
        ! For the impact-ionization, the heat of reaction is equal to the ionization energy
        IF(TRIM(ChemReac%ReactType(iReac)).EQ.'iQK') THEN
          IF(.NOT.ALLOCATED(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%ElectronicState)) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR: QK reactions require the definition of at least the ionization energy as electronic level!',iReac)
          END IF
          MaxElecQua=SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxElecQuant - 1
          ChemReac%EForm(iReac) = - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%ElectronicState(2,MaxElecQua)*BoltzmannConst
        END IF
      END DO
    END DO

    ! Initialize partition functions required for automatic backward rates
    IF(DSMC%BackwardReacRate) THEN
      IF(MOD(DSMC%PartitionMaxTemp,DSMC%PartitionInterval).EQ.0.0) THEN
        PartitionArraySize = NINT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)
      ELSE
        CALL abort(&
__STAMP__&
,'ERROR in Chemistry Init: Partition temperature limit must be multiple of partition interval!')
      END IF
      DO iSpec = 1, nSpecies
        ALLOCATE(SpecDSMC(iSpec)%PartitionFunction(1:PartitionArraySize))
        DO iInter = 1, PartitionArraySize
          Temp = iInter * DSMC%PartitionInterval
          CALL CalcPartitionFunction(iSpec, Temp, Qtra, Qrot, Qvib, Qelec)
          SpecDSMC(iSpec)%PartitionFunction(iInter) = Qtra * Qrot * Qvib * Qelec
        END DO
      END DO
    END IF
    ! Initialize analytic QK reaction rate (required for calculation of backward rate with QK and if multiple QK reactions can occur
    ! during a single collision, e.g. N2+e -> ionization or dissociation)
    IF(ANY(ChemReac%QKProcedure)) THEN
      IF(MOD(DSMC%PartitionMaxTemp,DSMC%PartitionInterval).EQ.0.0) THEN
        PartitionArraySize = NINT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)
      ELSE
        CALL abort(&
__STAMP__&
,'ERROR in Chemistry Init: Partition temperature limit must be multiple of partition interval!')
      END IF
      ALLOCATE(QKAnalytic(ChemReac%NumOfReact))
      DO iReac = 1, ChemReac%NumOfReact
        IF(ChemReac%QKProcedure(iReac)) THEN
          ! Calculation of the analytical rate, to be able to calculate the backward rate with partition function
          ALLOCATE(QKAnalytic(iReac)%ForwardRate(1:PartitionArraySize))
          DO iInter = 1, PartitionArraySize
            Temp = iInter * DSMC%PartitionInterval
            QKAnalytic(iReac)%ForwardRate(iInter) = CalcQKAnalyticRate(iReac,Temp)
          END DO
        END IF
      END DO
    END IF
    ! Filling up ChemReac-Array with forward rate coeff., switching reactants with products and setting new energies
    IF(DSMC%BackwardReacRate) THEN
      DO iReac = ChemReac%NumOfReact/2+1, ChemReac%NumOfReact
        iReacForward = iReac - ChemReac%NumOfReact/2
        IF(ChemReac%QKProcedure(iReacForward)) THEN
          IF((TRIM(ChemReac%ReactType(iReacForward)).EQ.'iQK').OR.&
             (TRIM(ChemReac%ReactType(iReacForward)).EQ.'D'  )) THEN
            ChemReac%ReactType(iReac) = 'r'
            IF(TRIM(ChemReac%ReactType(iReacForward)).EQ.'D') ChemReac%ReactType(iReac) = 'R'
            ChemReac%DefinedReact(iReac,1,1)      = ChemReac%DefinedReact(iReacForward,2,1)
            ! Products of the dissociation (which are the educts of the recombination) have to be swapped in order to comply with
            ! definition of the recombination reaction (e.g. CH3 + H + M -> CH4 + M but CH4 + M -> CH3 + M + H)
            ChemReac%DefinedReact(iReac,1,2)      = ChemReac%DefinedReact(iReacForward,2,3)
            ChemReac%DefinedReact(iReac,1,3)      = ChemReac%DefinedReact(iReacForward,2,2)
            ChemReac%DefinedReact(iReac,2,:)      = ChemReac%DefinedReact(iReacForward,1,:)
            ChemReac%EForm(iReac)                 = -ChemReac%EForm(iReacForward)
            ChemReac%EActiv(iReac) = 0.0
            QKAnalytic(iReac)%ForwardRate         = QKAnalytic(iReacForward)%ForwardRate
          ELSE
            CALL abort(__STAMP__,&
            'Other reaction types than iQK are not implemented with the automatic backward rate determination, Reaction:', iReac)
          END IF
        ELSE
          IF(TRIM(ChemReac%ReactType(iReacForward)).EQ.'D') THEN
            ! Analogous to the iQK case
            ChemReac%ReactType(iReac) = 'R'
            ChemReac%DefinedReact(iReac,1,1)      = ChemReac%DefinedReact(iReacForward,2,1)
            ChemReac%DefinedReact(iReac,1,2)      = ChemReac%DefinedReact(iReacForward,2,3)
            ChemReac%DefinedReact(iReac,1,3)      = ChemReac%DefinedReact(iReacForward,2,2)
            ChemReac%DefinedReact(iReac,2,:)      = ChemReac%DefinedReact(iReacForward,1,:)
            ChemReac%EActiv(iReac) = 0.0
          ELSEIF(TRIM(ChemReac%ReactType(iReacForward)).EQ.'E') THEN
            ChemReac%ReactType(iReac) = 'E'
            ChemReac%DefinedReact(iReac,1,:)      = ChemReac%DefinedReact(iReacForward,2,:)
            ChemReac%DefinedReact(iReac,2,:)      = ChemReac%DefinedReact(iReacForward,1,:)
            ChemReac%EActiv(iReac)                = ChemReac%EForm(iReacForward) + ChemReac%EActiv(iReacForward)
            IF(ChemReac%EActiv(iReac).LT.0.0) THEN
              ! The absolute value of the heat of formation cannot be larger than the activation energy but Arrhenius fits require
              ! sometimes a different value to better reproduce the experimental results. Doesnt matter for backward rate.
              ChemReac%EActiv(iReac) = 0.0
            END IF
          ELSE
            CALL abort(&
            __STAMP__&
            ,'Automatic calculation of backward reaction rate not supported with the chosen react type:',iReac)
          END IF
          ChemReac%Arrhenius_Prefactor(iReac)     = ChemReac%Arrhenius_Prefactor(iReacForward)
          ChemReac%Arrhenius_Powerfactor(iReac)   = ChemReac%Arrhenius_Powerfactor(iReacForward)
          ChemReac%EForm(iReac)                   = -ChemReac%EForm(iReacForward)
        END IF
      END DO
    END IF

    DO iReac = 1, ChemReac%NumOfReact
      ! Proof of reactant definition
      IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)*ChemReac%DefinedReact(iReac,1,3)).EQ.0) THEN
          CALL abort(__STAMP__,&
          'Recombination - Error in Definition: Not all reactant species are defined! ReacNbr: ',iReac)
        END IF
      ELSE
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)).EQ.0) THEN
          CALL abort(__STAMP__,&
          'Chemistry - Error in Definition: Reactant species not properly defined. ReacNbr:',iReac)
        END IF
      END IF
      ! Proof of product definition
      IF (TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)*ChemReac%DefinedReact(iReac,2,3)).EQ.0) THEN
          CALL abort(__STAMP__,&
          'Dissociation - Error in Definition: Not all product species are defined!  ReacNbr: ',iReac)
        END IF
        IF(ChemReac%DefinedReact(iReac,1,2).NE.ChemReac%DefinedReact(iReac,2,2)) THEN
          CALL abort(__STAMP__,&
          'Dissociation - Error in Definition: Non-reacting partner has to remain second product (/1,2,3/)  ReacNbr: ',iReac)
        END IF
      ELSE
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)).EQ.0) THEN
          CALL abort(__STAMP__,&
          'Chemistry - Error in Definition: Product species not properly defined. ReacNbr:',iReac)
        END IF
      END IF
      MaxSpecies = MAXVAL(ChemReac%DefinedReact(iReac,1:2,1:3))
      IF(MaxSpecies.GT.nSpecies) THEN
        CALL abort(__STAMP__,&
          'Chemistry - Error in Definition: Defined species does not exist, check number of species. ReacNbr:',iReac)
      END IF
    END DO

    ALLOCATE(PairCombID(nSpecies, nSpecies))
    PairCombID = 0
    CALL DSMC_BuildChem_IDArray(PairCombID)


    ! Case 6: One ionization and one ion recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'iQK').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF (((TRIM(ChemReac%ReactType(iReac2)).EQ.'r').AND.(.NOT.YetDefined_Help(iReac2))).OR.&
              ((TRIM(ChemReac%ReactType(iReac2)).EQ.'rQK').AND.(.NOT.YetDefined_Help(iReac2)))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 6
                  ChemReac%ReactCase(Reactant2, Reactant1) = 6
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 19: only ion recombination
    DO iReac = 1, ChemReac%NumOfReact
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'r') THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        Reactant3 = ChemReac%DefinedReact(iReac,1,3)
        IF(.NOT.YetDefined_Help(iReac)) THEN
          ChemReac%ReactCase(Reactant1, Reactant2) = 19
          ChemReac%ReactCase(Reactant2, Reactant1) = 19
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
        END IF
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
        DummyRecomb(Reactant1,Reactant2) = iReac
      END IF
    END DO
    ! Case 19: only ion recombination
    DO iReac = 1, ChemReac%NumOfReact
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'rQK') THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        Reactant3 = ChemReac%DefinedReact(iReac,1,3)
        IF(.NOT.YetDefined_Help(iReac)) THEN
          ChemReac%ReactCase(Reactant1, Reactant2) = 19
          ChemReac%ReactCase(Reactant2, Reactant1) = 19
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
        END IF
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
        DummyRecomb(Reactant1,Reactant2) = iReac
      END IF
    END DO
    ! Case 18: only electron impact ionization
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'iQK').AND.(.NOT.YetDefined_Help(iReac))) THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        ChemReac%ReactCase(Reactant1, Reactant2) = 18
        ChemReac%ReactCase(Reactant2, Reactant1) = 18
        ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
        ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
              Reactant1 = ChemReac%DefinedReact(iReac,1,1)
              Reactant2 = ChemReac%DefinedReact(iReac,1,2)
              ChemReac%ReactCase(Reactant1, Reactant2) = 20
              ChemReac%ReactCase(Reactant2, Reactant1) = 20
              ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
              ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
              ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
              ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
              YetDefined_Help(iReac2) = .TRUE.
            END IF
          END IF
        END DO
        YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 16: only simple CEX/MEX possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'x').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 16
          ChemReac%ReactCase(Reactant2, Reactant1) = 16
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 5/7/8/9/10/11/12: At least two dissociations possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF (iReac.NE.iReac2) THEN
            IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac2))) THEN
              IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                  PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                    ! Case 5: Two dissociations only
                    Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                    Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                    ChemReac%ReactCase(Reactant1, Reactant2) = 5
                    ChemReac%ReactCase(Reactant2, Reactant1) = 5
                    ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                    ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                    ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                    ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                    DO iReac3 = 1, ChemReac%NumOfReact
                      IF ((iReac3.NE.iReac2).AND.(iReac3.NE.iReac)) THEN
                        IF(DSMC%NumPolyatomMolecs.GT.0) THEN ! 3 diss only possible with polyatomic molecules
                          IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                            IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                  ! Case 7: Three dissociations
                                  Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                  Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                  ChemReac%ReactCase(Reactant1, Reactant2) = 7
                                  ChemReac%ReactCase(Reactant2, Reactant1) = 7
                                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                  ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                  ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                  ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                  ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                  YetDefined_Help(iReac3) = .TRUE.
                              DO iReac4 = 1, ChemReac%NumOfReact
                                IF ((iReac4.NE.iReac3).AND.(iReac4.NE.iReac2).AND.(iReac4.NE.iReac)) THEN
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 8: Four dissociations
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 8
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 8
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 9: Three dissociations and one exchange
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 9
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 9
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
                                END IF
                              END DO
                            END IF
                          END IF
                        END IF
                        IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                          IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                              PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                ! Case 10: Two dissociations and one exchange
                                Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                ChemReac%ReactCase(Reactant1, Reactant2) = 10
                                ChemReac%ReactCase(Reactant2, Reactant1) = 10
                                ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                YetDefined_Help(iReac3) = .TRUE.
                                DO iReac4 = 1, ChemReac%NumOfReact
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 11: Two dissociations, one exchange and recombination
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 11
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 11
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
                                END DO
                          END IF
                        END IF
                        IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                          IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                              PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                ! Case 12: Two dissociations and one recombination
                                Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                ChemReac%ReactCase(Reactant1, Reactant2) = 12
                                ChemReac%ReactCase(Reactant2, Reactant1) = 12
                                ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                YetDefined_Help(iReac3) = .TRUE.
                          END IF
                        END IF
                      END IF
                    END DO
                    YetDefined_Help(iReac) = .TRUE.
                    YetDefined_Help(iReac2) = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Cases 4/13/14: Dissociation, recombination and exchange reactions possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  ! Case 4: One dissociation and one exchange
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 4
                  ChemReac%ReactCase(Reactant2, Reactant1) = 4
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                  ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
                  DO iReac3 = 1, ChemReac%NumOfReact
                    IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                      IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                          PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                            ! Case 13: One dissociation, one recombination and one exchange
                            Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                            Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                            ChemReac%ReactCase(Reactant1, Reactant2) = 13
                            ChemReac%ReactCase(Reactant2, Reactant1) = 13
                            ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                            ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                            ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                            ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                            ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                            ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                            YetDefined_Help(iReac3) = .TRUE.
                      END IF
                    END IF
                  END DO
                  DO iReac3 = 1, ChemReac%NumOfReact
                    IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                      IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                          PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                            ! Case 17: One dissociation and two exchange reactions
                            Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                            Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                            ChemReac%ReactCase(Reactant1, Reactant2) = 17
                            ChemReac%ReactCase(Reactant2, Reactant1) = 17
                            ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                            ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                            ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                            ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                            ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                            ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                            YetDefined_Help(iReac3) = .TRUE.
                      END IF
                    END IF
                  END DO
            END IF
          END IF
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  ! Case 14: One dissociation and one recombination
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 14
                  ChemReac%ReactCase(Reactant2, Reactant1) = 14
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 15: One exchange and one recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 15
                  ChemReac%ReactCase(Reactant2, Reactant1) = 15
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 3: only exchange reaction possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 3
          ChemReac%ReactCase(Reactant2, Reactant1) = 3
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 2: only dissociation possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 2
          ChemReac%ReactCase(Reactant2, Reactant1) = 2
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 1: only molecular recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        Reactant3 = ChemReac%DefinedReact(iReac,1,3)
        IF(.NOT.YetDefined_Help(iReac)) THEN
          ChemReac%ReactCase(Reactant1, Reactant2) = 1
          ChemReac%ReactCase(Reactant2, Reactant1) = 1
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
        END IF
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
        DummyRecomb(Reactant1,Reactant2) = iReac
      END IF
    END DO
    ! Filling empty values of ReactNumRecomb with the last recombination reaction for that collision pair
    DO iReac = 1, ChemReac%NumOfReact
      IF((TRIM(ChemReac%ReactType(iReac)).EQ.'R').OR.(TRIM(ChemReac%ReactType(iReac)).EQ.'rQK') &
        .OR.(TRIM(ChemReac%ReactType(iReac)).EQ.'r')) THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        DO iSpec = 1, nSpecies
          IF(ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec).EQ.0) THEN
            ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec) = DummyRecomb(Reactant1,Reactant2)
            ChemReac%ReactNumRecomb(Reactant2, Reactant1, iSpec) = DummyRecomb(Reactant1,Reactant2)
          END IF
        END DO
      END IF
    END DO

    DEALLOCATE(PairCombID)
    DEALLOCATE(DummyRecomb)
    DEALLOCATE(ChemReac%ArbDiss)
    CALL Calc_Arrhenius_Factors()
  ELSE
    SWRITE(*,'(A)') 'NO REACTIONS DEFINED!'
  END IF

END SUBROUTINE DSMC_chemical_init


SUBROUTINE DSMC_BuildChem_IDArray(PairCombID)
!===================================================================================================================================
! Builds array for the identification of different species combinations
!===================================================================================================================================
! MODULES
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(INOUT)  :: PairCombID(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: iSpec, jSpec, iComb
!===================================================================================================================================
  iComb = 1
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec,  nSpecies
        PairCombID(iSpec, jSpec) = iComb
        PairCombID(jSpec, iSpec) = iComb
        iComb = iComb + 1
    END DO
  END DO

END SUBROUTINE DSMC_BuildChem_IDArray


SUBROUTINE Calc_Arrhenius_Factors()
!===================================================================================================================================
! Calculates variables for Arrhenius calculation (H_ab, Beta) for diatomic chemical reactions
!===================================================================================================================================
! MODULES
  USE MOD_Globals,            ONLY : Abort
  USE MOD_Globals_Vars,       ONLY: BoltzmannConst
  USE MOD_DSMC_Vars,          ONLY : ChemReac, SpecDSMC, CollInf
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies
  USE MOD_Globals_Vars,       ONLY : Pi
  USE MOD_DSMC_Analyze,       ONLY : CalcTVib
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: iQuaMax1, iQuaMax2, iQua1, iQua2, iReac, iQuaMax1_temp, iQuaMax2_temp
  INTEGER                 :: iQuaMax3, iQua3, iSpec
  REAL                    :: Xi_vib1, Xi_vib2, Xi_vib3
!  REAL                    :: EVib1, TVib1, EVib2, TVib2              ! needed for TSHO
!===================================================================================================================================

  DO iReac = 1, ChemReac%NumOfReact
    ! Calculate the Arrhenius arrays only if the reaction is not a QK reaction
    IF (.NOT.ChemReac%QKProcedure(iReac)) THEN
      ! Compute VHS Factor H_ab necessary for reaction probs, only defined for one omega for all species (laux diss page 24)
      ChemReac%Hab(iReac) = GAMMA(2.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) * 2.0 &
              * CollInf%Cab(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2))) &
              / ((1 + CollInf%KronDelta(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))) &
              * SQRT(Pi)) * (2.0 * BoltzmannConst &
              / CollInf%MassRed(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1), ChemReac%DefinedReact(iReac,1,2)))) &
              ** (0.5 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
      IF((.NOT.SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%PolyatomicMol) &
          .AND.(.NOT.SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%PolyatomicMol)) THEN
! ----------------------------------------------------------------------------------------------------------------------------------
! Recombination
! ----------------------------------------------------------------------------------------------------------------------------------
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
          IF(SpecDSMC(ChemReac%DefinedReact(iReac,1,3))%PolyatomicMol &
            .OR.SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) CYCLE
          iQuaMax3 = MAXVAL(SpecDSMC(:)%MaxVibQuant)
          ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(nSpecies,0:iQuaMax3-1))
          ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(:,:) = 0.0
          DO iSpec=1,nSpecies
            IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
              DO iQua3 = 0 , iQuaMax3-1
                IF (iQua3.NE.0) THEN
                  Xi_vib3 = 2 * iQua3 * LOG(1.0/ iQua3 + 1.0 )
                ELSE
                  Xi_vib3 = 0
                END IF
                ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iSpec,iQua3) = ChemReac%Arrhenius_Prefactor(iReac) &
                                         * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                         * GAMMA(( 3.0 + Xi_vib3 + SpecDSMC(iSpec)%Xi_Rot)/2 &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 ) &
                                         / (ChemReac%Hab(iReac) * GAMMA(( 3.0 + Xi_vib3 &
                                         + SpecDSMC(iSpec)%Xi_Rot)/2 &
                                         + 1.5 + ChemReac%Arrhenius_Powerfactor(iReac) ))
              END DO
            ELSE
              ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iSpec,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                                         * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                         * GAMMA(1.5 &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 ) &
                                         / ( ChemReac%Hab(iReac) * GAMMA(3.0 + ChemReac%Arrhenius_Powerfactor(iReac)))
            END IF
          END DO
        END IF
! ----------------------------------------------------------------------------------------------------------------------------------
! Dissociation
! ----------------------------------------------------------------------------------------------------------------------------------
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
          iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
          iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
          iQuaMax1_temp = 100
          iQuaMax2_temp = 100
          !!!!!!would be needed in TSHO Case!!!!!
          ! iQuaMax1_temp = iQuaMax1 - 1
          ! iQuaMax2_temp = iQuaMax2 - 1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(iQuaMax2.EQ.0) THEN
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(0:iQuaMax1-1,0:1))
          ELSE
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:iQuaMax2-1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
          END IF
          DO iQua1 = 0 , iQuaMax1_temp
            !Calculation of the vibrational DOF in TSHO Model
            IF (iQua1.NE.0) THEN
              Xi_vib1 = 2 * iQua1 * LOG(1.0/ iQua1 + 1.0 )
              !!!!!!would be needed in TSHO Case!!!!!
              !  EVib1 = (iQua1 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib
              !  TVib1 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib &
              !          , EVib1, iQuaMax1)
              !  Xi_vib1 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib / TVib1) - 1)) &
              !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1 / TVib1) - 1))
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ELSE
              Xi_vib1 = 0
            END IF

            IF(iQuaMax2.EQ.0) THEN
              Xi_vib2 = 0
              ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
              IF (Xi_vib1.GT.0) THEN
                IF (SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor.EQ.0) THEN
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2))
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(Xi_vib1/2) * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                      + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                      + Xi_vib1/2))
                END IF
              ELSE
                ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = 0
              END IF
            ELSE
              DO iQua2 = 0, iQuaMax2_temp
                !Calculation of the vibrational DOF in TSHO Model
                IF (iQua2.NE.0) THEN
                  Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                  !!!!!!would be needed in TSHO Case!!!!!
                  !  EVib2 = (iQua2 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib
                  ! TVib2 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib &
                  !          , EVib2, iQuaMax2)
                  !  Xi_vib2 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib / TVib2) - 1)) &
                  !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2 / TVib2) - 1))
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ELSE
                  Xi_vib2 = 0
                END IF
                !also only defined for one omega VHS for each species
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2) = Xi_vib1 + Xi_vib2  &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
                IF (Xi_vib1.GT.0) THEN
                  IF (SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor.EQ.0) THEN
                    ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) &
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2))
                  ELSE
                    ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) &
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(Xi_vib1/2) * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                        + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                        + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                        + Xi_vib1/2))
                  END IF
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = 0
                END IF
              END DO
            END IF
          END DO
          !! would not be needed in TSHO Case
          IF (iQuaMax2.EQ.0) THEN
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
          ELSE
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(100,100)
          END IF
        END IF
! ----------------------------------------------------------------------------------------------------------------------------------
! Exchange
! ----------------------------------------------------------------------------------------------------------------------------------
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'E') THEN
          iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
          iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
          iQuaMax1_temp = 100
          iQuaMax2_temp = 100
          !!!!!!would be needed in TSHO Case!!!!!
          ! iQuaMax1_temp = iQuaMax1 - 1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(iQuaMax1.EQ.0) THEN
            IF(iQuaMax2.EQ.0) THEN
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:1, 0:1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:1,0:1))
            ELSE
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:1, 0:iQuaMax2-1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:1,0:iQuaMax2-1))
            END IF
          ELSE
            IF(iQuaMax2.EQ.0) THEN
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:iQuaMax1-1,0:1))
            ELSE
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:iQuaMax2-1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
            END IF
          END IF
          IF(iQuaMax1.EQ.0) THEN
            IF(iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(0,0) = 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
              ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(0,0)/2) &
                    / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(0,0)/2))
            ELSE
              DO iQua2 = 0, iQuaMax2_temp
                ! Calculation of the vibrational DOF in TSHO Model
                IF (iQua2.NE.0) THEN
                  Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                ELSE
                  Xi_vib2 = 0
                END IF
                ! also only defined for one omega VHS for each species
                ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2) = Xi_vib2  &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
                ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) &
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2)/2) &
                    / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2)/2))
              END DO
            END IF
          ELSE
            DO iQua1 = 0 , iQuaMax1_temp
              !Calculation of the vibrational DOF in TSHO Model
              IF (iQua1.NE.0) THEN
                Xi_vib1 = 2 * iQua1 * LOG(1.0/ iQua1 + 1.0 )
                !!!!!!would be needed in TSHO Case!!!!!
                !  EVib1 = (iQua1 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib
                !  TVib1 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib &
                !          , EVib1, iQuaMax1)
                !  Xi_vib1 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib  &
                !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib / TVib1) - 1)) &
                !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1  &
                !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1 / TVib1) - 1))
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ELSE
                Xi_vib1 = 0
              END IF
              IF(iQuaMax2.EQ.0) THEN
                Xi_vib2 = 0
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                          + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
                IF (Xi_vib1.GT.0) THEN
                  ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2))
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,0) = 0
                END IF
              ELSE
                DO iQua2 = 0, iQuaMax2_temp
                  !Calculation of the vibrational DOF in TSHO Model
                  IF (iQua2.NE.0) THEN
                    Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                    !!!!!!would be needed in TSHO Case!!!!!
                    !  EVib2 = (iQua2 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib
                    ! TVib2 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib &
                    !          , EVib2, iQuaMax2)
                    !  Xi_vib2 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib  &
                    !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib / TVib2) - 1)) &
                    !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2  &
                    !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2 / TVib2) - 1))
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ELSE
                    Xi_vib2 = 0
                  END IF
                  !also only defined for one omega VHS for each species
                  ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2) = Xi_vib1 + Xi_vib2  &
                          + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                          + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                          + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
                  IF (Xi_vib1.GT.0) THEN
                    ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) &
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2))
                  ELSE
                    ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,iQua2) = 0
                  END IF
                END DO
              END IF
            END DO
          END IF
          !! would not be needed in TSHO Case
          IF (iQuaMax1.EQ.0) THEN
            IF (iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(0,0)
            ELSE
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(0,100)
            END IF
          ELSE
            IF (iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
            ELSE
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(100,100)
            END IF
          END IF
        END IF ! Exchange reaction
      END IF ! First and second educt not polyatomic
    END IF ! Not Q-K
  END DO

END SUBROUTINE Calc_Arrhenius_Factors


SUBROUTINE Init_TLU_Data
!===================================================================================================================================
! Reads Scattering Angle Lookup Table from Test_Lookup_komplett.txt
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,          ONLY : TLU_Data, ChemReac
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

!-- parameters
INTEGER,PARAMETER             :: unit1=20
!DOUBLE PRECISION, PARAMETER   :: mass_ion=2.180d-25 !Xenon
!DOUBLE PRECISION, PARAMETER   :: mass_neutral=2.180d-25 !Xenon

!-- local variables
INTEGER                       :: io_error1,read_error,iLine
CHARACTER(LEN=1000000)        :: string1
INTEGER                       :: N_b, N_E
DOUBLE PRECISION              :: real1, real2
!===================================================================================================================================

OPEN(UNIT=unit1,file=TRIM(ChemReac%TLU_FileName(ChemReac%NumOfReact)),STATUS='old',ACTION='READ',IOSTAT=io_error1)
IF ( io_error1 .EQ. 0) THEN
  !----------------------------------Schleife ueber alle Zeilen in file1----------------------------------!
  iLine=0
  DO
    !----------------------------------Einlesen----------------------------------!
    READ(unit1,'(A)',IOSTAT=read_error) string1
    IF ( read_error .GT. 0) THEN
      CALL abort(__STAMP__,&
        'Chemistry - Error in Init_TLU_Data, Error:',read_error)
    ELSE IF (read_error .LT. 0 ) THEN
      EXIT ! Dateiende erreicht
    ELSE
      iLine=iLine+1
      IF (iLine.EQ.1) THEN
        READ(string1,*,IOSTAT=read_error) real1, real2
        N_b=NINT(real1)
        N_E=NINT(real2)
        ALLOCATE( TLU_Data%Chitable(1:N_E,1:N_b))
        ALLOCATE( TLU_Data%deltabj(1:N_E))
      ELSE IF (iLine.EQ.2) THEN
        READ(string1,*,IOSTAT=read_error) TLU_Data%Emin, TLU_Data%Emax, TLU_Data%deltaE
      ELSE IF (iLine.EQ.3) THEN
        READ(string1,*,IOSTAT=read_error) TLU_Data%deltabj(1:N_E)
      ELSE IF (iLine.EQ.4 .OR. iLine.EQ.5) THEN
        !dummy lines...
      ELSE IF (iLine.GE.6 .AND. iLine.LE.N_b+5) THEN
        READ(string1,*,IOSTAT=read_error) TLU_Data%Chitable(:,iLine-5)
      ELSE
        CALL abort(__STAMP__,&
          'Chemistry - Error in Init_TLU_Data, File too long, Error:',read_error)
      END IF
      IF ( read_error .NE. 0) THEN
        !STOP "Datenfehler im gelesenen String!"
        CALL abort(__STAMP__,&
          'Chemistry - Error in Init_TLU_Data, Data error in loaded string,Error:',read_error)
      END IF
    END IF
  END DO
ELSE
  !STOP "Datenfehler im gelesenen String!"
        CALL abort(__STAMP__,&
          'Chemistry - Error in Init_TLU_Data, Error while opening of Test_Lookup_komplett.txt, Error:',io_error1)
  !WRITE(*,'(A,I0,A)') 'Beim Oeffenen der Datei Test_Lookup_komplett.txt ist Fehler Nr. ', io_error1,' aufgetreten!'
END IF

! Force Chi at N_b to be 0
TLU_Data%Chitable(:,N_b) = 0

CLOSE(unit=unit1)
END SUBROUTINE Init_TLU_Data

END MODULE MOD_DSMC_ChemInit
