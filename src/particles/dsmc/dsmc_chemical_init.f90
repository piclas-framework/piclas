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

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DefineParametersChemistry
PUBLIC :: DSMC_chemical_init
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for DSMC (Direct Simulation Monte Carlo)
!==================================================================================================================================
SUBROUTINE DefineParametersChemistry()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection("DSMC Chemistry")
CALL prms%CreateIntOption(      'DSMC-NumOfReactions'  &
                                           ,'Number of reactions.')
CALL prms%CreateIntOption(      'DSMC-Reaction[$]-NumberOfNonReactives'  &
                                           ,'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-NonReactiveSpecies'  &
                                           ,'Array with the non-reactive collision partners for dissociation'&
                                           ,numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'DSMC-Reaction[$]-ReactionType'  &
                                           ,'Used reaction type\n'//&
                                            'I: electron impact ionization\n'//&
                                            'R: molecular recombination\n'//&
                                            'D: molecular dissociation\n'//&
                                            'E: molecular exchange reaction\n'//&
                                            'X: simple charge exchange reaction)', 'none', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DSMC-Reaction[$]-QKProcedure'  &
                                           ,'Flag to use quantum-kinetic model', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DSMC-Reaction[$]-QK-Method'  &
                                           ,'Recombination Method for Q-K model\n'//&
                                            '1: by Bird\n'//&
                                            '2: by Gallis)\n'//&
                                            'If using bird, define the variables:\n'//&
                                            'DSMC-Reaction[$]-QK-Coeff1\n'//&
                                            'DSMC-Reaction[$]-QK-Coeff2 ', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-QK-Coeff1'  &
                                           ,'First Q-K coefficient for Birds method.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-QK-Coeff2'  &
                                           ,'Second Q-K coefficient for Birds method.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2,\n'//&
                                            'SpecNumOfReactant3)', '0 , 0 , 0' , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Products'  &
                                           ,'Products of Reaction[j] (Product1, Product2, Product3)', '0 , 0 , 0' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Prefactor'  &
                                           , 'TODO-DEFINE-PARAMETER ', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Powerfactor'  &
                                           , 'TODO-DEFINE-PARAMETER', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Activation-Energy_K'  &
                                           , 'Activation energy (relativ to k_Boltzmann) for Reaction[$].', '0.' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CEXa'  &
                                , 'CEX log-factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '-27.2' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CEXb'  &
                                , 'CEX const. factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '175.269' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DSMC-Reaction[$]-DoScat'  &
                                , 'Perform scattering-based charge-exchange instead of isotropic '//&
                                '(model of Samuel Araki by lookup table)', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-ELa'  &
                                , 'with DoScat=T: EL log-factor '//&
                                '(g&cut-off-angle-dep. cs in Angstrom, def.: value for Xe+)', '-26.8' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-ELb'  &
                                , 'with DoScat=T: EL const. factor '//&
                                '(g&cut-off-angle-dep. cs in Angstrom, def.: value for Xe+)', '148.975' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-MEXa'  &
                                , 'with DoScat=F: MEX log-factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '-27.2' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-MEXb'  &
                                , 'with DoScat=F: MEX const. factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '175.269' , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(     'DSMC-Reaction[$]-TLU_FileName'  &
                                , 'with DoScat=F: No TLU-File needed '//&
                                '(def.: )', '0' , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Particles-Chemistry-NumDeleteProducts','Number of species, which should be deleted if they are '//&
                                'a product of chemical reactions', '0')
CALL prms%CreateIntArrayOption( 'Particles-Chemistry-DeleteProductsList','List of the species indices to be deleted if they are '//&
                                'a product of chemical reactions')

CALL prms%CreateLogicalOption(  'Part-Species[$]-UseCollXSec'  &
                                           ,'Utilize collision cross sections for the determination of collision probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseVibXSec'  &
                                           ,'Utilize vibrational cross sections for the determination of relaxation probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Particles-CollXSec-Database', 'File name for the collision cross section database. Container '//&
                                                               'should be named with species pair (e.g. "Ar-electron"). The '//&
                                                               'first column shall contain the energy in eV and the second '//&
                                                               'column the cross-section in m^2', 'none')
CALL prms%CreateLogicalOption(  'Particles-CollXSec-NullCollision'  &
                                  ,'Utilize the null collision method for the determination of the number of pairs '//&
                                  'based on the maximum collision frequency and time step (only with a background gas)' &
                                  ,'.TRUE.')
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CrossSection'  &
                                , 'Photon-ionization cross-section', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersChemistry

SUBROUTINE DSMC_chemical_init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars               ,ONLY: ChemReac, DSMC, SpecDSMC, BGGas
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies, Species
USE MOD_Particle_Analyze_Vars   ,ONLY: ChemEnergySum
USE MOD_DSMC_ChemReact          ,ONLY: CalcPartitionFunction
USE MOD_part_emission_tools     ,ONLY: CalcPhotonEnergy
USE MOD_DSMC_QK_PROCEDURES      ,ONLY: QK_Init
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iReac, iReac2, iReac3, iReac4, iSpec, iChemDir, iReacForward, iReacDiss, PartitionArraySize, iInter
INTEGER, ALLOCATABLE  :: PairCombID(:,:), DummyRecomb(:,:)
LOGICAL, ALLOCATABLE  :: YetDefined_Help(:)
LOGICAL               :: DoScat
INTEGER               :: Reactant1, Reactant2, Reactant3, MaxSpecies, MaxElecQua, ReadInNumOfReact
REAL                  :: Temp, Qtra, Qrot, Qvib, Qelec, BGGasEVib, PhotonEnergy
INTEGER               :: iInit
!===================================================================================================================================

ChemReac%NumOfReact = GETINT('DSMC-NumOfReactions')
ReadInNumOfReact = ChemReac%NumOfReact
IF(ChemReac%NumOfReact.LE.0) THEN
  CALL Abort(&
    __STAMP__&
    ,' CollisMode = 3 requires a chemical reaction database. DSMC-NumOfReactions cannot be zero!')
END IF
ALLOCATE(ChemReac%ArbDiss(ChemReac%NumOfReact))
! Allowing unspecified non-reactive collision partner (CH4 + M -> CH3 + H + M, e.g. (/1,0,0/) -> (/2,0,3/)
iReacDiss = ChemReac%NumOfReact
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  ChemReac%ArbDiss(iReac)%NumOfNonReactives = GETINT('DSMC-Reaction'//TRIM(hilf)//'-NumberOfNonReactives')
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
!----------------------------------------------------------------------------------------------------------------------------------
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
ALLOCATE(ChemReac%CrossSection(ChemReac%NumOfReact))
ChemReac%CrossSection = 0.

IF (BGGas%NumberOfSpecies.GT.0) THEN
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      ! Background gas: Calculation of the mean vibrational quantum number of diatomic molecules
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        IF(.NOT.SpecDSMC(iSpec)%PolyatomicMol) THEN
          BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib &
            + BoltzmannConst * SpecDSMC(iSpec)%CharaTVib / (EXP(SpecDSMC(iSpec)%CharaTVib / SpecDSMC(iSpec)%Init(1)%TVib) - 1)
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
  ChemReac%ReactType(iReac)             = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-ReactionType'))
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
    ChemReac%ELa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-ELa','-26.8')
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
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'phIon') THEN
    ChemReac%CrossSection(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CrossSection')
    ! Check if species 3 is an electron and abort (this is not implemented yet)
    IF(ChemReac%DefinedReact(iReac,2,3).GT.0)THEN
      IF(SpecDSMC(ChemReac%DefinedReact(iReac,2,3))%InterID.EQ.4) CALL abort(&
        __STAMP__&
        ,'Chemical reaction with electron as 3rd product species. This is not implemented yet for photoionization! iReac=',&
        IntInfoOpt=iReac)
    END IF ! ChemReac%DefinedReact(iReac,2,3).GT.0
  END IF
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
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'phIon') THEN
    PhotonEnergy = 0.
    DO iSpec = 1, nSpecies
      DO iInit = 1, Species(iSpec)%NumberOfInits
        IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'photon_cylinder') THEN
          PhotonEnergy = CalcPhotonEnergy(Species(iSpec)%Init(iInit)%WaveLength)
          EXIT
        END IF
      END DO
    END DO
    ChemReac%EForm(iReac) = ChemReac%EForm(iReac) + PhotonEnergy
    IF(ChemReac%EForm(iReac).LE.0.0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Photon energy is not sufficient for the given ionization reaction: ',iReac)
    END IF
  END IF
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

! Initialize analytic QK reaction rate (required for calculation of backward rate with QK and if multiple QK reactions can occur
! during a single collision, e.g. N2+e -> ionization or dissociation)
IF(ANY(ChemReac%QKProcedure).OR.DSMC%BackwardReacRate) THEN
  CALL QK_Init()
END IF

DO iReac = 1, ChemReac%NumOfReact
  ! Proof of reactant definition
  IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)*ChemReac%DefinedReact(iReac,1,3)).EQ.0) THEN
      CALL abort(__STAMP__,&
      'Recombination - Error in Definition: Not all reactant species are defined! ReacNbr: ',iReac)
    END IF
  ELSE IF (TRIM(ChemReac%ReactType(iReac)).NE.'phIon') THEN
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

! Case: photo-ionization (NOT to be included in regular chemistry)
DO iReac = 1, ChemReac%NumOfReact
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'phIon') THEN
    IF(.NOT.YetDefined_Help(iReac)) THEN
      YetDefined_Help(iReac) = .TRUE.
    END IF
  END IF
END DO

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
USE MOD_Globals         ,ONLY: Abort
USE MOD_Globals_Vars    ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars       ,ONLY: ChemReac, SpecDSMC, CollInf
USE MOD_PARTICLE_Vars   ,ONLY: nSpecies
USE MOD_Globals_Vars    ,ONLY: Pi
USE MOD_DSMC_Analyze    ,ONLY: CalcTVib
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iQuaMax1, iQuaMax2, iQuaMax3, iQua1, iQua2, iQua3, iReac, iSpec, iSpec1, iSpec2
REAL                    :: omega, Xi_rot1, Xi_rot2, Xi_vib1, Xi_vib2, Xi_vib3, Xi_rel_rot, Xi_total
!===================================================================================================================================

DO iReac = 1, ChemReac%NumOfReact
  ! Skip photo-ionization reactions
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'phIon') CYCLE
  ! Skip reactions modelled with QK
  IF(ChemReac%QKProcedure(iReac)) CYCLE
  iSpec1 = ChemReac%DefinedReact(iReac,1,1)
  iSpec2 = ChemReac%DefinedReact(iReac,1,2)
  omega = CollInf%omega(iSpec1,iSpec2)
  ! Compute VHS Factor H_ab necessary for reaction probs, only defined for one omega for all species (laux diss page 24)
  ChemReac%Hab(iReac) = GAMMA(2.0 - omega) * 2.0 * CollInf%Cab(CollInf%Coll_Case(iSpec1,iSpec2)) &
                        / ((1 + CollInf%KronDelta(CollInf%Coll_Case(iSpec1,iSpec2))) * SQRT(Pi)) &
                        * (2.0 * BoltzmannConst / CollInf%MassRed(CollInf%Coll_Case(iSpec1, iSpec2))) ** (0.5 - omega)
  ! Skip reactions with polyatomic educts, currently beta is calculated on the fly
  IF((SpecDSMC(iSpec1)%PolyatomicMol).OR.(SpecDSMC(iSpec2)%PolyatomicMol)) CYCLE
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
        iQuaMax3 = SpecDSMC(iSpec)%MaxVibQuant
      ELSE
        iQuaMax3 = 1
      END IF
      DO iQua3 = 0 , iQuaMax3-1
        IF (iQua3.NE.0) THEN
          Xi_vib3 = 2. * iQua3 * LOG(1.0/iQua3 + 1.0)
        ELSE
          Xi_vib3 = 0
        END IF
        Xi_total = 3. + SpecDSMC(iSpec)%Xi_Rot + Xi_vib3
        ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iSpec,iQua3) = ChemReac%Arrhenius_Prefactor(iReac) &
                                  * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) - omega) &
                                  * GAMMA(Xi_total/2. - omega + 2.) / (ChemReac%Hab(iReac) * GAMMA(Xi_total/2. &
                                  + 1.5 + ChemReac%Arrhenius_Powerfactor(iReac) ))
      END DO
    END DO
  END IF
! ----------------------------------------------------------------------------------------------------------------------------------
! Dissociation & Exchange
! ----------------------------------------------------------------------------------------------------------------------------------
  IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').OR.(TRIM(ChemReac%ReactType(iReac)).EQ.'E')) THEN
    IF ((SpecDSMC(iSpec1)%InterID.EQ.2).OR.(SpecDSMC(iSpec1)%InterID.EQ.20)) THEN
      iQuaMax1 = SpecDSMC(iSpec1)%MaxVibQuant
      Xi_Rot1 = SpecDSMC(iSpec1)%Xi_Rot
    ELSE
      iQuaMax1 = 1
      Xi_Rot1 = 0.
    END IF
    IF ((SpecDSMC(iSpec2)%InterID.EQ.2).OR.(SpecDSMC(iSpec2)%InterID.EQ.20)) THEN
      iQuaMax2 = SpecDSMC(iSpec2)%MaxVibQuant
      Xi_Rot2 = SpecDSMC(iSpec2)%Xi_Rot
    ELSE
      iQuaMax2 = 1
      Xi_Rot2 = 0.
    END IF
    Xi_rel_rot = 2.*(2. - omega) + Xi_Rot1 + Xi_Rot2
    ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
    DO iQua1 = 0 , iQuaMax1 - 1
      IF (iQua1.NE.0) THEN
        Xi_vib1 = 2. * iQua1 * LOG(1.0/iQua1 + 1.0)
      ELSE
        Xi_vib1 = 0.
      END IF
      DO iQua2 = 0, iQuaMax2 - 1
        IF (iQua2.NE.0) THEN
          Xi_vib2 = 2. * iQua2 * LOG(1.0/iQua2 + 1.0)
        ELSE
          Xi_vib2 = 0
        END IF
        Xi_total = Xi_rel_rot + Xi_vib1 + Xi_vib2
        ChemReac%ReactInfo(iReac)%Beta_Arrhenius(iQua1,iQua2) = Calc_Beta_Arrhenius_Diatomic(iReac,Xi_total)
      END DO
    END DO
  END IF
END DO

END SUBROUTINE Calc_Arrhenius_Factors


PURE REAL FUNCTION Calc_Beta_Arrhenius_Diatomic(iReac,Xi_total)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals           ,ONLY: Abort
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: ChemReac, CollInf
USE MOD_DSMC_Analyze      ,ONLY: CalcTVib
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iReac
REAL, INTENT(IN)          :: Xi_total
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iSpec1, iSpec2
REAL                      :: omega, A_Arrhenius, b_Arrhenius
!===================================================================================================================================

iSpec1 = ChemReac%DefinedReact(iReac,1,1)
iSpec2 = ChemReac%DefinedReact(iReac,1,2)

omega = CollInf%omega(iSpec1,iSpec2)
A_Arrhenius = ChemReac%Arrhenius_Prefactor(iReac)
b_Arrhenius = ChemReac%Arrhenius_Powerfactor(iReac)

Calc_Beta_Arrhenius_Diatomic = A_Arrhenius * (BoltzmannConst**(0.5 - b_Arrhenius - omega)) * GAMMA(Xi_total/2.) &
                                / (ChemReac%Hab(iReac) * GAMMA(b_Arrhenius - 0.5 + omega + Xi_total/2.))

END FUNCTION Calc_Beta_Arrhenius_Diatomic


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
