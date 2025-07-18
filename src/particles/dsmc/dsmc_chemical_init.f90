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
PUBLIC :: InitReactionPaths
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
CALL prms%SetSection("Chemistry")
CALL prms%CreateStringOption(   'DSMC-ChemistryModel', 'Name of available chemistry reaction model, e.g.:\n'//&
                                  'Air_5Spec_8Reac_Park1993\n'//&
                                  'Mars_11Spec_27Reac_Johnston2014\n'//&
                                  'Titan_14Spec_24Reac_Savajano2011\n'//&
                                  'A complete list can be found in the documentation.', 'none')
CALL prms%CreateIntOption(      'DSMC-NumOfReactions','Number of chemical reactions')
CALL prms%CreateIntOption(      'DSMC-Reaction[$]-NumberOfNonReactives', &
                                          'Number of non-reactive collision partners (Length of the read-in vector)', &
                                          '0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-NonReactiveSpecies'  &
                                           ,'Array with the non-reactive collision partners for dissociation' &
                                           ,numberedmulti=.TRUE., no=0)
CALL prms%CreateStringOption(   'DSMC-Reaction[$]-ReactionModel'  &
                                           ,'Used reaction model\n'//&
                                            'TCE: Total Collision Energy\n'//&
                                            'phIon: photon-ionization\n'//&
                                            'phIonXSec: photon-ionization with cross-section-based data for the reaction\n'//&
                                            'QK: quantum kinetic\n'//&
                                            'XSec: cross-section-based data for the reaction', 'TCE', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2,\n'//&
                                            'SpecNumOfReactant3)', numberedmulti=.TRUE.,no=3)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Products'  &
                                        ,'Products of Reaction[j] (Product1, Product2, Product3, Product 4)',numberedmulti=.TRUE., no=4)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Prefactor', &
                                    'Prefactor A of the extended Arrhenius equation, k = A * T^b * EXP(-E_a/T), '//&
                                    'Units: 1/s, m3/s, m6/s (depending on the type of the reaction)', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Powerfactor', &
                                    'Temperature exponent b of the extended Arrhenius equation, k = A * T^b * EXP(-E_a/T), '//&
                                    'Units: -', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Activation-Energy_K', &
                                    'Activation energy E_a of the extended Arrhenius equation, k = A * T^b * EXP(-E_a/T), '//&
                                    'Units: Kelvin', '0.', numberedmulti=.TRUE.)
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
                                , 'with DoScat=F: No TLU-File needed ', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Particles-Chemistry-NumDeleteProducts','Number of species, which should be deleted if they are '//&
                                'a product of chemical reactions', '0')
CALL prms%CreateIntArrayOption( 'Particles-Chemistry-DeleteProductsList','List of the species indices to be deleted if they are '//&
                                'a product of chemical reactions', no=0)

CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CrossSection', &
                                'Photon-ionization cross-section for the reaction type phIon', numberedmulti=.TRUE.)

CALL prms%SetSection("Chemistry - Backward Reaction Rates")
CALL prms%CreateLogicalOption(  'Particles-DSMC-BackwardReacRate', &
                                          'Set [TRUE] to enable the automatic calculation of the backward reaction rate '//&
                                          'coefficient using the equilibrium constant calculated by partition functions\n'//&
                                          '[FALSE] if they are defined as separate reactions. This option can be overwritten '//&
                                          'by the reaction-specific flag, e.g. DSMC-Reaction1-BackwardReac = T/F ' , '.FALSE.')
CALL prms%CreateLogicalOption(  'DSMC-Reaction[$]-BackwardReac', &
                                          'Consider automatic backward reaction for this specific reaction. The default value '//&
                                          'is equal the value of the global flag Particles-DSMC-BackwardReacRate.', &
                                          numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Particles-DSMC-PartitionMaxTemp', &
                                          'Define temperature limit for pre-stored partition function that are used for '//&
                                          'calculation of backwards rates', '20000.0')
CALL prms%CreateRealOption(     'Particles-DSMC-PartitionInterval', &
                                          'Define temperature interval for pre-stored partition functions that are used for '//&
                                          'calculation of backwards rates', '10.0')
CALL prms%CreateIntOption(      'Part-Species[$]-SymmetryFactor', &
                                          'Rotational symmetry factor, depending on the molecule configuration', &
                                          numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersChemistry

SUBROUTINE DSMC_chemical_init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi
USE MOD_DSMC_Vars               ,ONLY: ChemReac, DSMC, SpecDSMC, BGGas, CollInf
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies, Species, SpeciesDatabase
USE MOD_Particle_Analyze_Vars   ,ONLY: ChemEnergySum
USE MOD_DSMC_QK_Chemistry       ,ONLY: QK_Init
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcXiVib
USE MOD_MCC_Vars                ,ONLY: NbrOfPhotonXsecReactions
USE MOD_io_hdf5
USE MOD_HDF5_input              ,ONLY: ReadAttribute, DatasetExists, AttributeExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iReac, iReac2, iSpec, iPart, iReacDiss, iSpec2
INTEGER, ALLOCATABLE  :: DummyRecomb(:,:)
LOGICAL               :: DoScat
REAL                  :: BGGasEVib, omega, ChargeProducts, ChargeReactants
INTEGER               :: Reactant1, Reactant2, Reactant3, MaxSpecies, ReadInNumOfReact, CustomNumReact
!===================================================================================================================================

NbrOfPhotonXsecReactions = 0
ChemReac%AnyQKReaction = .FALSE.
ChemReac%AnyXSecReaction = .FALSE.

! Read-in of the name of the chemistry model (number of reactions is then given by the available reactions)
ChemReac%ChemistryModel  = TRIM(GETSTR('DSMC-ChemistryModel'))

IF(ChemReac%ChemistryModel.EQ.'none') THEN
  ChemReac%NumOfReact = GETINT('DSMC-NumOfReactions')
  ReadInNumOfReact = ChemReac%NumOfReact
  ChemReac%NumOfReactWOBackward = ChemReac%NumOfReact
  ALLOCATE(ChemReac%ArbDiss(ChemReac%NumOfReact))
  ! Allowing unspecified non-reactive collision partner (CH4 + M -> CH3 + H + M, e.g. (/1,0,0/) -> (/2,0,3/)
  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    ChemReac%ArbDiss(iReac)%NumOfNonReactives = GETINT('DSMC-Reaction'//TRIM(hilf)//'-NumberOfNonReactives')
    IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.GT.0) THEN
      ALLOCATE(ChemReac%ArbDiss(iReac)%NonReactiveSpecies(ChemReac%ArbDiss(iReac)%NumOfNonReactives))
      ChemReac%ArbDiss(iReac)%NonReactiveSpecies = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-NonReactiveSpecies', &
                                                ChemReac%ArbDiss(iReac)%NumOfNonReactives)
      ! First reaction is saved within the dummy input reaction, thus "- 1"
      ChemReac%NumOfReact = ChemReac%NumOfReact + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
      ChemReac%NumOfReactWOBackward = ChemReac%NumOfReactWOBackward + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
    END IF
  END DO
ELSE
  ! Number of reactions is set based on the chemistry model
  CALL GetNumReacFromDatabase(ReadInNumOfReact)
  CustomNumReact = GETINT('DSMC-NumOfReactions','0')
  IF(CustomNumReact.GT.0) THEN
    CALL Abort(__STAMP__,' Custom reactions in addition to a chemistry model are not supported yet!')
  END IF
  IF(ChemReac%NumOfReact.LE.0) THEN
    CALL Abort(__STAMP__,' No reactions found for the selected chemistry model: '//TRIM(ChemReac%ChemistryModel))
  END IF
END IF

IF(ChemReac%NumOfReact.LE.0) THEN
  CALL Abort(__STAMP__,' CollisMode = 3 requires a chemical reaction database. DSMC-NumOfReactions cannot be zero!')
END IF

! Set counter for filling-up the reaction array with the non-reactives at the end of the read-in reaction list.
iReacDiss = ReadInNumOfReact

!-----------------------------------------------------------------------------------
! Flag for the automatic calculation of the backward reaction rate with the partition functions and equilibrium constant.
DSMC%BackwardReacRate  = GETLOGICAL('Particles-DSMC-BackwardReacRate')
! Partition functions are calculated for each species during initialization and stored for values starting with the
! DSMC%PartitionInterval up to DSMC%PartitionMaxTemp, interpolation between the stored values (also used for analytic QK reactions)
DSMC%PartitionMaxTemp  = GETREAL('Particles-DSMC-PartitionMaxTemp','20000')
DSMC%PartitionInterval = GETREAL('Particles-DSMC-PartitionInterval','10')
ALLOCATE(ChemReac%BackwardReac(ChemReac%NumOfReactWOBackward))
ChemReac%BackwardReac = .FALSE.
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  ! Calculation of the backward reaction rates
  IF(DSMC%BackwardReacRate) THEN
    ChemReac%BackwardReac(iReac) = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-BackwardReac','.TRUE.')
  ELSE
    ChemReac%BackwardReac(iReac) = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-BackwardReac','.FALSE.')
  END IF
  IF (ChemReac%BackwardReac(iReac)) THEN
    ChemReac%NumOfReact = ChemReac%NumOfReact + 1
    IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.GT.0) &
      ChemReac%NumOfReact = ChemReac%NumOfReact + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
  END IF
END DO
IF (ANY(ChemReac%BackwardReac)) DSMC%BackwardReacRate = .TRUE.
!-----------------------------------------------------------------------------------
! Delete products if they belong to a certain species
ChemReac%NumDeleteProducts = GETINT('Particles-Chemistry-NumDeleteProducts')
IF(ChemReac%NumDeleteProducts.GT.0) THEN
  ALLOCATE(ChemReac%DeleteProductsList(ChemReac%NumDeleteProducts))
  ChemReac%DeleteProductsList = GETINTARRAY('Particles-Chemistry-DeleteProductsList', ChemReac%NumDeleteProducts)
END IF

ChemEnergySum = 0.
CALL PrintOption('Number of considered reaction paths (including dissociation and recombination)','INFO',IntOpt=ChemReac%NumOfReact)
!----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(ChemReac%NumReac(ChemReac%NumOfReact))
ChemReac%NumReac = 0
IF (DSMC%ReservoirSimu) THEN
  ALLOCATE(ChemReac%ReacCount(ChemReac%NumOfReact))
  ChemReac%ReacCount = 0
  ALLOCATE(ChemReac%ReacCollMean(CollInf%NumCase))
  ChemReac%ReacCollMean = 0.0
END IF
ALLOCATE(ChemReac%ReactType(ChemReac%NumOfReact))
ChemReac%ReactType = '0'
ALLOCATE(ChemReac%ReactModel(ChemReac%NumOfReact))
ChemReac%ReactModel = '0'
ALLOCATE(ChemReac%ReactionName(ChemReac%NumOfReact))
ALLOCATE(ChemReac%Reactants(ChemReac%NumOfReact,3))
ChemReac%Reactants = 0
ALLOCATE(ChemReac%Products(ChemReac%NumOfReact,4))
ChemReac%Products = 0
ALLOCATE(ChemReac%ReactCase(ChemReac%NumOfReact))
ChemReac%ReactCase = 0
ALLOCATE(ChemReac%ReactNum(nSpecies, nSpecies, nSpecies))
ChemReac%ReactNum = 0
ALLOCATE(ChemReac%ReactNumRecomb(nSpecies, nSpecies, nSpecies))
ChemReac%ReactNumRecomb = 0
ALLOCATE(ChemReac%Arrhenius_Prefactor(ChemReac%NumOfReact), ChemReac%Arrhenius_Powerfactor(ChemReac%NumOfReact),&
          ChemReac%EActiv(ChemReac%NumOfReact), ChemReac%EForm(ChemReac%NumOfReact), ChemReac%Hab(ChemReac%NumOfReact))
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
      IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
        IF(.NOT.SpecDSMC(iSpec)%PolyatomicMol) THEN
          IF(DSMC%VibAHO) THEN ! AHO
            CALL CalcXiVib(SpecDSMC(iSpec)%Init(1)%TVib, iSpec, XiVibTotal=ChemReac%MeanXiVib_PerIter(iSpec))
          ELSE ! SHO
            BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib &
              + BoltzmannConst * SpecDSMC(iSpec)%CharaTVib / (EXP(SpecDSMC(iSpec)%CharaTVib / SpecDSMC(iSpec)%Init(1)%TVib) - 1)
            BGGasEVib = BGGasEVib/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) - DSMC%GammaQuant
            ChemReac%MeanEVibQua_PerIter(iSpec) = MIN(INT(BGGasEVib) + 1, SpecDSMC(iSpec)%MaxVibQuant)
            ChemReac%MeanXiVib_PerIter(iSpec) = 2. * ChemReac%MeanEVibQua_PerIter(iSpec) &
                                              * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0 )
          END IF
        END IF
      ELSE
        ChemReac%MeanEVibQua_PerIter(iSpec) = 0
        ChemReac%MeanXiVib_PerIter(iSpec) = 0.
      END IF
    END IF
  END DO
END IF

DoScat = .false.
IF(SpeciesDatabase.EQ.'none') THEN
  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    ChemReac%ReactModel(iReac)  = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-ReactionModel'))
    ChemReac%Reactants(iReac,:) = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Reactants',3)
    ChemReac%Products(iReac,:)  = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Products',4)
    SELECT CASE (TRIM(ChemReac%ReactModel(iReac)))
      CASE('TCE')
        ! Total Collision Energy: Arrhenius-based chemistry model
        ChemReac%Arrhenius_Prefactor(iReac)   = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Prefactor')
        ChemReac%Arrhenius_Powerfactor(iReac) = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Powerfactor')
        ChemReac%EActiv(iReac)                = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Activation-Energy_K')*BoltzmannConst
      CASE('QK')
        ! Quantum Kinetic: Threshold energy based chemistry model
        ChemReac%AnyQKReaction = .TRUE.
      CASE('XSec')
        ! Chemistry model based on cross-section data
        ChemReac%AnyXSecReaction = .TRUE.
      CASE('CEX')
        ! Simple charge exchange reactions
        ChemReac%CEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXa','-27.2')
        ChemReac%CEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXb','175.269')
        ChemReac%DoScat(iReac)               = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-DoScat','.FALSE.')
        IF (ChemReac%DoScat(iReac)) THEN
          DoScat = .true.
          ChemReac%ELa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-ELa','-26.8')
          ChemReac%ELb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-ELb','148.975')
          ChemReac%TLU_FileName(iReac)        = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-TLU_FileName'))
        ELSE
          ChemReac%MEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXa','-27.2')
          ChemReac%MEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXb','175.269')
        END IF
      CASE('phIon')
        ! Photo-ionization reactions
        ChemReac%CrossSection(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CrossSection')
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE('phIonXSec')
        ! Photo-ionization reactions (data read-in from database)
        NbrOfPhotonXsecReactions = NbrOfPhotonXsecReactions + 1
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE DEFAULT
        CALL abort(__STAMP__,'Selected reaction model is not supported in reaction number: ', IntInfoOpt=iReac)
    END SELECT
  END DO
ELSE
  CALL ReadReacFromDatabase(ReadInNumOfReact)
END IF

! Filling up ChemReac-Array for the given non-reactive dissociation/electron-impact ionization partners
DO iReac = 1, ReadInNumOfReact
  IF((ChemReac%Reactants(iReac,2).EQ.0).AND.(ChemReac%Products(iReac,2).EQ.0)) THEN
    IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.EQ.0) THEN
      CALL abort(__STAMP__,'Error in Definition: Non-reacting partner(s) has to be defined!',IntInfoOpt=iReac)
    END IF
    DO iReac2 = 1, ChemReac%ArbDiss(iReac)%NumOfNonReactives
      IF(iReac2.EQ.1) THEN
        ! The first non-reacting partner is written into the reaction, which was read-in.
        ChemReac%Reactants(iReac,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
        ChemReac%Products(iReac,2)      = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
      ELSE
        ! The following reaction are added after the number of originally read-in reactions (counter: iReacDiss)
        ChemReac%BackwardReac(iReacDiss+iReac2-1)           = ChemReac%BackwardReac(iReac)
        ChemReac%ReactModel(iReacDiss+iReac2-1)             = ChemReac%ReactModel(iReac)
        ChemReac%Reactants(iReacDiss+iReac2-1,:)            = ChemReac%Reactants(iReac,:)
        ChemReac%Reactants(iReacDiss+iReac2-1,2)            = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
        ChemReac%Products(iReacDiss+iReac2-1,:)             = ChemReac%Products(iReac,:)
        ChemReac%Products(iReacDiss+iReac2-1,2)             = ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iReac2)
        ChemReac%Arrhenius_Prefactor(iReacDiss+iReac2-1)    = ChemReac%Arrhenius_Prefactor(iReac)
        ChemReac%Arrhenius_Powerfactor(iReacDiss+iReac2-1)  = ChemReac%Arrhenius_Powerfactor(iReac)
        ChemReac%EActiv(iReacDiss+iReac2-1)                 = ChemReac%EActiv(iReac)
      END IF
    END DO
    iReacDiss = iReacDiss + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
  ELSE IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.NE.0) THEN
    CALL abort(__STAMP__,'Dissociation/Ionization - Error in Definition: Non-reacting partner(s) has to be zero!',iReac)
  END IF
END DO

! Automatic determination of the reaction type based on the number of reactants, products and charge balance (for ionization)
DO iReac = 1, ChemReac%NumOfReact
  ! Skip photo-ionization reactions
  IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) CYCLE
  IF (ChemReac%Reactants(iReac,3).NE.0) THEN
    ChemReac%ReactType(iReac)           = 'R'
  ELSE IF (ChemReac%Products(iReac,3).NE.0) THEN
    ChargeReactants = ABS(Species(ChemReac%Reactants(iReac,1))%ChargeIC) + ABS(Species(ChemReac%Reactants(iReac,2))%ChargeIC)
    ChargeProducts = ABS(Species(ChemReac%Products(iReac,1))%ChargeIC) + ABS(Species(ChemReac%Products(iReac,2))%ChargeIC) &
        + ABS(Species(ChemReac%Products(iReac,3))%ChargeIC)
    IF (ChemReac%Products(iReac,4).GT.0) ChargeProducts = ChargeProducts + ABS(Species(ChemReac%Products(iReac,4))%ChargeIC)
    IF (ChargeReactants.NE.ChargeProducts) THEN
      ChemReac%ReactType(iReac)         = 'I'
    ELSE
      ChemReac%ReactType(iReac)         = 'D'
    END IF
  ELSE
    ChemReac%ReactType(iReac)           = 'E'
  END IF
END DO

IF (DoScat) CALL abort(__STAMP__,'Deactivated Init_TLU_Data() routine')
!IF (DoScat) CALL Init_TLU_Data()

! Calculation of stoichiometric coefficients and calculation of the heat of formation
DO iReac = 1, ChemReac%NumOfReact
  ALLOCATE(ChemReac%ReactInfo(iReac)%StoichCoeff(nSpecies,2))
  ChemReac%ReactInfo(iReac)%StoichCoeff(1:nSpecies,1:2) = 0
  ChemReac%EForm(iReac) = 0.0
  DO iSpec=1, nSpecies
    ! Reactants
    DO iPart = 1,3
      IF(ChemReac%Reactants(iReac,iPart).EQ.iSpec) THEN
          ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,1) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,1) + 1
      END IF
    END DO
    ! Products
    DO iPart = 1,4
      IF(ChemReac%Products(iReac,iPart).EQ.iSpec) THEN
          ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,2) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,2) + 1
      END IF
    END DO
    ! Calculation of the enthalpy of reaction by the summation of the enthalpies of formation of the respective species
    ! (ionization energy of ionized species was already added in dsmc_init.f90)
    ChemReac%EForm(iReac) = ChemReac%EForm(iReac) &
                          - ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,2)*SpecDSMC(iSpec)%HeatOfFormation  &
                          + ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,1)*SpecDSMC(iSpec)%HeatOfFormation
    ! For the impact-ionization, the heat of reaction is equal to the ionization energy
    IF(TRIM(ChemReac%ReactType(iReac)).EQ.'I') THEN
      IF(.NOT.ALLOCATED(SpecDSMC(ChemReac%Reactants(iReac,1))%ElectronicState)) CALL abort(&
        __STAMP__,'ERROR: Ionization reactions require the definition of at least the ionization energy as electronic level!',iReac)
    END IF
  END DO
END DO ! iReac = 1, ChemReac%NumOfReact

! Photoionization
CALL InitPhotonReactions()

! Populate the background reaction arrays and initialize the required partition functions
IF(DSMC%BackwardReacRate) THEN
  CALL DSMC_BackwardRate_init()
END IF

! Check for possible input errors in the chemistry definition
DO iReac = 1, ChemReac%NumOfReact
  ! Proof of recombination definition
  IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    IF ((ChemReac%Reactants(iReac,1)*ChemReac%Reactants(iReac,2)*ChemReac%Reactants(iReac,3)).EQ.0) THEN
      CALL abort(__STAMP__,'Recombination - Error in Definition: Not all reactant species are defined! ReacNbr: ',iReac)
    END IF
    IF (ChemReac%Reactants(iReac,3).NE.ChemReac%Products(iReac,2)) THEN
      CALL abort(__STAMP__,&
      'Recombination - Error in Definition: Third-collision partner does not correspond to the second product! ReacNbr: ',iReac)
    END IF
  ELSE IF (.NOT.StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) THEN
    IF ((ChemReac%Reactants(iReac,1)*ChemReac%Reactants(iReac,2)).EQ.0) THEN
      CALL abort(__STAMP__,'Chemistry - Error in Definition: Reactant species not properly defined. ReacNbr:',iReac)
    END IF
  END IF
  ! Proof of dissociation definition
  IF (TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
    ! Three product species are given
    IF ((ChemReac%Products(iReac,1)*ChemReac%Products(iReac,2)*ChemReac%Products(iReac,3)).EQ.0) THEN
      CALL abort(__STAMP__,'Dissociation - Error in Definition: Not all product species are defined!  ReacNbr: ',iReac)
    END IF
    IF(TRIM(ChemReac%ReactModel(iReac)).NE.'XSec') THEN
      ! Cross-section based chemistry does not require this definition as no backward reaction rates are implemented
      IF(ChemReac%Reactants(iReac,2).NE.ChemReac%Products(iReac,2)) THEN
        CALL abort(__STAMP__,&
        'Dissociation - Error in Definition: Non-reacting partner has to remain second product (/1,2,3/)  ReacNbr: ',iReac)
      END IF
    END IF
    ! At least a molecule is given as a reactant
    IF((Species(ChemReac%Reactants(iReac,1))%InterID.NE.2).AND.(Species(ChemReac%Reactants(iReac,1))%InterID.NE.20) &
      .AND.(Species(ChemReac%Reactants(iReac,2))%InterID.NE.2).AND.(Species(ChemReac%Reactants(iReac,2))%InterID.NE.20)) THEN
      CALL abort(__STAMP__,&
        'Dissociation - Error in Definition: None of the reactants is a molecule, check species indices and charge definition. ReacNbr: ',iReac)
    END IF
  ELSE
    IF ((ChemReac%Products(iReac,1)*ChemReac%Products(iReac,2)).EQ.0) THEN
      CALL abort(__STAMP__,'Chemistry - Error in Definition: Product species not properly defined. ReacNbr:',iReac)
    END IF
  END IF
  ! Check if the maximum species index is not greater than the number of species
  MaxSpecies = MAXVAL(ChemReac%Reactants(iReac,1:3))
  IF(MaxSpecies.GT.nSpecies) CALL abort(__STAMP__,&
      'Chemistry - Error in Definition: Defined species does not exist, check number of species. ReacNbr:',iReac)
END DO

! Initialize analytic QK reaction rate (required for calculation of backward rate with QK and if multiple QK reactions can occur
! during a single collision, e.g. N2+e -> ionization or dissociation)
IF(ChemReac%AnyQKReaction.OR.DSMC%BackwardReacRate) THEN
  CALL QK_Init()
END IF

! Count the number of possible reactions paths per collision case
ALLOCATE(ChemReac%CollCaseInfo(CollInf%NumCase))
ChemReac%CollCaseInfo(:)%NumOfReactionPaths = 0

! Initialize reaction paths
CALL InitReactionPaths()

! Recombination: saving the reaction index based on the third species
DO iReac = 1, ChemReac%NumOfReact
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    Reactant1 = ChemReac%Reactants(iReac,1)
    Reactant2 = ChemReac%Reactants(iReac,2)
    Reactant3 = ChemReac%Reactants(iReac,3)
    ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
    ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
    DummyRecomb(Reactant1,Reactant2) = iReac
  END IF
END DO
! Filling empty values of ReactNumRecomb with the last recombination reaction for that collision pair
DO iReac = 1, ChemReac%NumOfReact
  IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
    Reactant1 = ChemReac%Reactants(iReac,1)
    Reactant2 = ChemReac%Reactants(iReac,2)
    DO iSpec = 1, nSpecies
      IF(ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec).EQ.0) THEN
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec) = DummyRecomb(Reactant1,Reactant2)
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, iSpec) = DummyRecomb(Reactant1,Reactant2)
      END IF
    END DO
  END IF
END DO

DEALLOCATE(DummyRecomb)
DEALLOCATE(ChemReac%ArbDiss)

DO iReac = 1, ChemReac%NumOfReact
  ! Skip reactions that are not TCE
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'TCE') CYCLE
  iSpec = ChemReac%Reactants(iReac,1)
  iSpec2 = ChemReac%Reactants(iReac,2)
  omega = CollInf%omega(iSpec,iSpec2)
  ! Compute VHS Factor H_ab necessary for reaction probs, only defined for one omega for all species (laux diss page 24)
  ChemReac%Hab(iReac) = GAMMA(2.0 - omega) * 2.0 * CollInf%Cab(CollInf%Coll_Case(iSpec,iSpec2)) &
                        / ((1 + CollInf%KronDelta(CollInf%Coll_Case(iSpec,iSpec2))) * SQRT(Pi)) &
                        * (2.0 * BoltzmannConst / CollInf%MassRed(CollInf%Coll_Case(iSpec, iSpec2))) ** (0.5 - omega)
END DO

END SUBROUTINE DSMC_chemical_init


!===================================================================================================================================
!> Initialize all possible reaction paths
!===================================================================================================================================
SUBROUTINE InitReactionPaths()
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars ,ONLY: ChemReac,CollInf
#if USE_HDG
USE MOD_Particle_Vars ,ONLY: Species
USE MOD_HDG_Vars  ,ONLY: UseBRElectronFluid ! Used for skipping reactions involving electrons as products
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iReac,iCase,iCase2,ReacIndexCounter
LOGICAL     :: RecombAdded
#if USE_HDG
INTEGER     :: iProd
#endif /*USE_HDG*/
!===================================================================================================================================
DO iCase = 1, CollInf%NumCase
  RecombAdded = .FALSE.
  REACLOOP: DO iReac = 1, ChemReac%NumOfReact
#if USE_HDG
    ! Skip reactions involving electrons as products
    IF(UseBRElectronFluid) THEN
      DO iProd = 1,4
        IF(ChemReac%Products(iReac,iProd).NE.0)THEN
          IF(Species(ChemReac%Products(iReac,iProd))%InterID.EQ.4) CYCLE REACLOOP
        END IF ! ChemReac%Products(iReac,iProd).NE.0
      END DO
    END IF
#endif /*USE_HDG*/
    ! Skip the special case of photo ionization
    IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) CYCLE
    iCase2 = CollInf%Coll_Case(ChemReac%Reactants(iReac,1),ChemReac%Reactants(iReac,2))
    IF(iCase.EQ.iCase2) THEN
      ! Only add recombination reactions once
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        IF(RecombAdded) THEN
          CYCLE
        ELSE
          RecombAdded = .TRUE.
        END IF
      END IF
      ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths = ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths + 1
    END IF
  END DO REACLOOP
END DO

DO iCase = 1, CollInf%NumCase
  ! Allocate the case specific type with the number of the possible reaction paths
  ALLOCATE(ChemReac%CollCaseInfo(iCase)%ReactionIndex(ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths))
  ChemReac%CollCaseInfo(iCase)%ReactionIndex = 0
  ALLOCATE(ChemReac%CollCaseInfo(iCase)%ReactionProb(ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths))
  ChemReac%CollCaseInfo(iCase)%ReactionProb = 0.
  ChemReac%CollCaseInfo(iCase)%HasXSecReaction    = .FALSE.
  ReacIndexCounter = 0
  RecombAdded = .FALSE.
  REACLOOP2: DO iReac = 1, ChemReac%NumOfReact
#if USE_HDG
    ! Skip reactions involving electrons as products
    IF(UseBRElectronFluid) THEN
      DO iProd = 1,4
        IF(ChemReac%Products(iReac,iProd).NE.0)THEN
          IF(Species(ChemReac%Products(iReac,iProd))%InterID.EQ.4) CYCLE REACLOOP2
        END IF ! ChemReac%Products(iReac,iProd).NE.0
      END DO
    END IF
#endif /*USE_HDG*/
    ! Skip the special case of photo ionization
    IF(StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) CYCLE
    iCase2 = CollInf%Coll_Case(ChemReac%Reactants(iReac,1),ChemReac%Reactants(iReac,2))
    ! Save the reaction index for the specific collision case
    IF(iCase.EQ.iCase2) THEN
      ! Save the case number for the reaction
      ChemReac%ReactCase(iReac) = iCase
      ! But only add one recombination reaction to the number of reaction paths (the index of the others is stored in ReactNumRecomb
      ! and is chosen based on the third collision partner selected during the simulation)
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        IF(RecombAdded) THEN
          CYCLE
        ELSE
          RecombAdded = .TRUE.
        END IF
      END IF
      ReacIndexCounter = ReacIndexCounter + 1
      ChemReac%CollCaseInfo(iCase)%ReactionIndex(ReacIndexCounter) = iReac
      IF(TRIM(ChemReac%ReactModel(iReac)).EQ.'XSec')  THEN
        ChemReac%CollCaseInfo(iCase)%HasXSecReaction = .TRUE.
        IF(ChemReac%Reactants(iReac,3).NE.0) THEN
          CALL abort(__STAMP__,&
            'Chemistry - Error: Cross-section based chemistry for reactions with three reactants is not supported yet!')
        END IF
      END IF
    END IF
  END DO REACLOOP2
END DO

END SUBROUTINE InitReactionPaths


SUBROUTINE DSMC_BackwardRate_init()
!===================================================================================================================================
!> Initialize and read-in of variables required for the automatic backward reaction rate calculation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_DSMC_Vars               ,ONLY: ChemReac, DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies, Species, SpeciesDatabase
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcPartitionFunction
USE MOD_io_hdf5
USE MOD_HDF5_input              ,ONLY: ReadAttribute
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iSpec, iReac, iReacForward, iPolyatMole, PartitionArraySize, iInter,err
REAL                  :: Temp, Qtra, Qrot, Qvib, Qelec
CHARACTER(LEN=64)     :: dsetname
INTEGER(HID_T)        :: file_id_specdb                       ! File identifier
!===================================================================================================================================

IF(ChemReac%AnyXSecReaction) CALL abort(__STAMP__,&
  'Automatic calculation of backward reaction rates in combination with cross-section based reactions is NOT supported!')

! 1.) Read-in of species parameters for the partition function calculation
DO iSpec = 1, nSpecies
  IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    IF(SpeciesDatabase.NE.'none') THEN
      ! Initialize FORTRAN interface.
      CALL H5OPEN_F(err)
      CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
      dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
      CALL ReadAttribute(file_id_specdb,'SymmetryFactor',1,DatasetName = dsetname,IntScalar=SpecDSMC(iSpec)%SymmetryFactor,ReadFromGroup=.TRUE.)
      CALL PrintOption('SymmetryFactor '//TRIM(Species(iSpec)%Name),'DB',IntOpt=SpecDSMC(iSpec)%SymmetryFactor)
      ! Close the file.
      CALL H5FCLOSE_F(file_id_specdb, err)
      ! Close FORTRAN interface.
      CALL H5CLOSE_F(err)
    END IF

    IF(Species(iSpec)%DoOverwriteParameters) THEN
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      SpecDSMC(iSpec)%SymmetryFactor = GETINT('Part-Species'//TRIM(hilf)//'-SymmetryFactor')
    END IF

    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
        IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
          CALL abort(__STAMP__,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
        END IF
      ELSE
        IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)  &
            * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
          CALL abort(__STAMP__,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
        END IF
      END IF
    ELSE
      IF(SpecDSMC(iSpec)%CharaTRot*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
        CALL abort(__STAMP__,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
      END IF
    END IF
  END IF
  IF((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized).AND.(Species(iSpec)%InterID.NE.100)) THEN
    IF(.NOT.ALLOCATED(SpecDSMC(iSpec)%ElectronicState)) THEN
      CALL abort(__STAMP__,'ERROR: Electronic energy levels required for the calculation of backward reaction rate!',iSpec)
    END IF
  END IF
END DO

! 2.) Initialize the tabulated partition functions

IF(MOD(DSMC%PartitionMaxTemp,DSMC%PartitionInterval).EQ.0.0) THEN
  PartitionArraySize = NINT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)
ELSE
  CALL abort(__STAMP__,'ERROR in Chemistry Init: Partition temperature limit must be multiple of partition interval!')
END IF
DO iSpec = 1, nSpecies
  ALLOCATE(SpecDSMC(iSpec)%PartitionFunction(1:PartitionArraySize))
  DO iInter = 1, PartitionArraySize
    Temp = iInter * DSMC%PartitionInterval
    CALL CalcPartitionFunction(iSpec, Temp, Qtra, Qrot, Qvib, Qelec)
    SpecDSMC(iSpec)%PartitionFunction(iInter) = Qtra * Qrot * Qvib * Qelec
  END DO
END DO

! 3.) Filling up ChemReac-Array with forward rate coefficients, switching reactants with products and setting new energies
ALLOCATE(ChemReac%BackwardReacForwardIndx(ChemReac%NumOfReactWOBackward+1:ChemReac%NumOfReact))
iReac = ChemReac%NumOfReactWOBackward
DO iReacForward = 1, ChemReac%NumOfReactWOBackward
  IF (.NOT.ChemReac%BackwardReac(iReacForward)) CYCLE
  iReac = iReac + 1
  ChemReac%BackwardReacForwardIndx(iReac) = iReacForward
  IF(TRIM(ChemReac%ReactModel(iReacForward)).EQ.'QK') THEN
    IF((TRIM(ChemReac%ReactType(iReacForward)).EQ.'I').OR.&
        (TRIM(ChemReac%ReactType(iReacForward)).EQ.'D'  )) THEN
      ChemReac%ReactType(iReac) = 'R'
      ChemReac%ReactModel(iReac) = 'TCE'
      ChemReac%Reactants(iReac,1)      = ChemReac%Products(iReacForward,1)
      ! Products of the dissociation (which are the reactants of the recombination) have to be swapped in order to comply with
      ! definition of the recombination reaction (e.g. CH3 + H + M -> CH4 + M but CH4 + M -> CH3 + M + H)
      ChemReac%Reactants(iReac,2)      = ChemReac%Products(iReacForward,3)
      ChemReac%Reactants(iReac,3)      = ChemReac%Products(iReacForward,2)
      ChemReac%Products(iReac,1:3)     = ChemReac%Reactants(iReacForward,1:3)
      ChemReac%EForm(iReac)            = -ChemReac%EForm(iReacForward)
      ChemReac%EActiv(iReac) = 0.0
    ELSE
      CALL abort(__STAMP__,&
      'Other reaction types than I and D are not implemented with the automatic backward rate determination, Reaction:', iReac)
    END IF
  ELSE
    IF((TRIM(ChemReac%ReactType(iReacForward)).EQ.'I').OR.&
        (TRIM(ChemReac%ReactType(iReacForward)).EQ.'D'  )) THEN
      ! Analogous to the I case
      ChemReac%ReactType(iReac) = 'R'
      ChemReac%ReactModel(iReac) = 'TCE'
      ChemReac%Reactants(iReac,1)      = ChemReac%Products(iReacForward,1)
      ChemReac%Reactants(iReac,2)      = ChemReac%Products(iReacForward,3)
      ChemReac%Reactants(iReac,3)      = ChemReac%Products(iReacForward,2)
      ChemReac%Products(iReac,1:3)     = ChemReac%Reactants(iReacForward,1:3)
      ChemReac%EActiv(iReac) = 0.0
    ELSEIF(TRIM(ChemReac%ReactType(iReacForward)).EQ.'E') THEN
      ChemReac%ReactType(iReac) = 'E'
      ChemReac%ReactModel(iReac) = 'TCE'
      ChemReac%Reactants(iReac,1:3)      = ChemReac%Products(iReacForward,1:3)
      ChemReac%Products(iReac,1:3)       = ChemReac%Reactants(iReacForward,1:3)
      ChemReac%EActiv(iReac)             = ChemReac%EForm(iReacForward) + ChemReac%EActiv(iReacForward)
      IF(ChemReac%EActiv(iReac).LT.0.0) THEN
        ! The absolute value of the heat of formation cannot be larger than the activation energy but Arrhenius fits require
        ! sometimes a different value to better reproduce the experimental results. Doesnt matter for backward rate.
        ChemReac%EActiv(iReac) = 0.0
      END IF
    ELSE
      CALL abort(__STAMP__,'Automatic calculation of backward reaction rate not supported with the chosen react type:',iReac)
    END IF
    ChemReac%Arrhenius_Prefactor(iReac)     = ChemReac%Arrhenius_Prefactor(iReacForward)
    ChemReac%Arrhenius_Powerfactor(iReac)   = ChemReac%Arrhenius_Powerfactor(iReacForward)
    ChemReac%EForm(iReac)                   = -ChemReac%EForm(iReacForward)
  END IF
END DO

END SUBROUTINE DSMC_BackwardRate_init


SUBROUTINE GetNumReacFromDatabase(ReadInNumOfReact)
!===================================================================================================================================
!> Determine the number of reactions depending on the selected chemistry model
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Particle_Vars     ,ONLY: SpeciesDatabase, Species, nSpecies
USE MOD_DSMC_Vars         ,ONLY: ChemReac
USE MOD_HDF5_Input        ,ONLY: DatasetExists,nDims,HSize,ReadAttribute,GetDataSize,AttributeExists
USE MOD_StringTools       ,ONLY: STRICMP, SPLIT_STRING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)               :: ReadInNumOfReact
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=256)                :: dsetname, groupname, dsetname2, NonReactiveSpeciesString
CHARACTER(LEN=255),ALLOCATABLE    :: ModelNames(:), NonReactiveSpeciesNames(:)
INTEGER                           :: storage, max_corder, err
INTEGER                           :: totalNumReac, ModelNamesSize, iModel, numReac, iSpec, jSpec, iNonReac, NumChemModels, NumOfCols, LenNonReacArray
INTEGER(HID_T)                    :: file_id_specdb, group_id, dset_id             ! File/group/dataset identifiers
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iReac, ReacID
LOGICAL                           :: GroupFound, AttrExists, SpeciesFound
!===================================================================================================================================

! 1) Open the species database
! Initialize FORTRAN interface
CALL H5OPEN_F(err)
! Open the file
CALL H5FOPEN_F(TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
! Check if the REACTIONS group exists
groupname = 'Reactions'
CALL H5LEXISTS_F(file_id_specdb, TRIM(groupname), GroupFound, err)
! Abort if the group does not exist
IF(.NOT.GroupFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: No reactions group found!')
! Open the group and get number of available reactions in the database
CALL H5GOPEN_F(file_id_specdb,TRIM(groupname), group_id, err)
CALL H5Gget_info_f(group_id, storage, totalNumReac, max_corder, err)
IF(totalNumReac.EQ.0) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: No reactions found!')

! Loop over all the reactions and filter out the selected chemistry model
numReac = 0
ALLOCATE(ChemReac%totalReacToModel(totalNumReac))
DO iReac = 0, totalNumReac-1
  ! Get name and size of name
  CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, iReac, dsetname, err, size)
  dsetname2 = TRIM(groupname)//'/'//TRIM(dsetname)
  ! Open the reaction dataset
  CALL H5DOPEN_F(file_id_specdb, dsetname2, dset_id, err)
  ! Get the size of the attribute array
  CALL GetDataSize(dset_id,'ChemistryModel',nDims,HSize,attrib=.TRUE.)
  IF (nDims.EQ.1) THEN  ! Catch for old database where ChemistryModel was not an 2D array
    NumChemModels = 1
  ELSE                  ! NumChemModels of current reaction equals number of rows of attribute
    NumChemModels=INT(HSize(2),4)
  END IF
  ! NumOfCols = 1 if Non Reactives dont exist, = 2 if Non Reactives exist
  NumOfCols = INT(HSize(1),4)
  ModelNamesSize = NumChemModels*NumOfCols
  IF(ALLOCATED(ModelNames)) THEN
    DEALLOCATE(ModelNames)
    ALLOCATE(ModelNames(ModelNamesSize))
  ELSE
    ALLOCATE(ModelNames(ModelNamesSize))
  END IF
  ModelNames = ''
  CALL ReadAttribute(file_id_specdb,'ChemistryModel',ModelNamesSize,TRIM(dsetname2),StrArray=ModelNames)
  DO iModel = 1,NumChemModels
    ! ModelNames contains strings of alterning ChemistryModel and (if existent) corrensponding NonReactives, thus index of ModelNames to check if ChemistryModel matches
    IF(STRICMP(ChemReac%ChemistryModel,ModelNames((iModel - 1) * NumOfCols + 1))) THEN
      numReac = numReac + 1
      ChemReac%totalReacToModel(numReac) = INT(iReac,4)
    END IF
  END DO
  CALL H5DCLOSE_F(dset_id, err)
END DO

! Set the number of reactions as found in the database
ChemReac%NumOfReact = numReac
ChemReac%NumOfReactWOBackward = numReac
ReadInNumOfReact = numReac
ALLOCATE(ChemReac%ArbDiss(numReac))

! Look for the definition of the non-reactive species
DO iReac = 1, numReac
  ReacID = ChemReac%totalReacToModel(iReac)
  ! Get name and size of name
  CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, ReacID, dsetname, err, size)
  dsetname2 = TRIM(groupname)//'/'//TRIM(dsetname)
  ! Open the reaction dataset
  CALL H5DOPEN_F(file_id_specdb, dsetname2, dset_id, err)
  ! Get the size of the attribute array
  CALL GetDataSize(dset_id,'ChemistryModel',nDims,HSize,attrib=.TRUE.)
  IF (nDims.EQ.1) THEN  ! Catch for old database where ChemistryModel was not an 2D array
    NumChemModels = 1
    ! dummy -1 to handle old case
    NumOfCols = -1
  ELSE                  ! NumChemModels of current reaction equals number of rows of attribute
    NumChemModels=INT(HSize(2),4)
    NumOfCols = INT(HSize(1),4)
  END IF
  ! Check if NonReactives exist
  IF(NumOfCols.EQ.2)THEN    ! new format of database
    ModelNamesSize = NumChemModels*NumOfCols
    IF(ALLOCATED(ModelNames)) THEN
      DEALLOCATE(ModelNames)
      ALLOCATE(ModelNames(ModelNamesSize))
    ELSE
      ALLOCATE(ModelNames(ModelNamesSize))
    END IF
    ModelNames = ''
    CALL ReadAttribute(file_id_specdb,'ChemistryModel',ModelNamesSize,TRIM(dsetname2),StrArray=ModelNames)
    NonReactiveSpeciesString = ''
    DO iModel = 1,NumChemModels
      ! ModelNames contains string of alterning ChemistryModel and (if existent) corrensponding NonReactives, thus index of ModelNames
      IF(STRICMP(ChemReac%ChemistryModel,ModelNames((iModel - 1) * NumOfCols + 1))) THEN
        ! saving NonReactives of selected ChemistryModel in comma separated string
        NonReactiveSpeciesString = ADJUSTL(Trim(ModelNames((iModel - 1) * NumOfCols + 2)))
      END IF
    END DO
    LenNonReacArray = 0
    DO iNonReac = 1, LEN_TRIM(NonReactiveSpeciesString)
      ! count elements (commas) for array allocation
      IF(NonReactiveSpeciesString(iNonReac:iNonReac).EQ.',')THEN
        LenNonReacArray = LenNonReacArray + 1
      END IF
    END DO
    ! Allocate memory NonReactiveSpecies Array with count of commas +1 since string does not end with a comma
    ALLOCATE(NonReactiveSpeciesNames(LenNonReacArray + 1))
    ALLOCATE(ChemReac%ArbDiss(iReac)%NonReactiveSpecies(LenNonReacArray + 1))
    CALL SPLIT_STRING(NonReactiveSpeciesString,',',NonReactiveSpeciesNames,ChemReac%ArbDiss(iReac)%NumOfNonReactives)
    ! Find the index of the respective species
    DO iSpec = 1, ChemReac%ArbDiss(iReac)%NumOfNonReactives
      SpeciesFound = .FALSE.
      DO jSpec = 1, nSpecies
        IF(STRICMP(NonReactiveSpeciesNames(iSpec),Species(jSpec)%Name)) THEN
          ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iSpec) = jSpec
          SpeciesFound = .TRUE.
          EXIT
        END IF
      END DO
      IF(.NOT.SpeciesFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: Species defined as non-reactive has not beend found!')
    END DO
    DEALLOCATE(NonReactiveSpeciesNames)
    DEALLOCATE(ModelNames)
    ! First reaction is saved within the dummy input reaction, thus "- 1"
    ChemReac%NumOfReact = ChemReac%NumOfReact + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
    ChemReac%NumOfReactWOBackward = ChemReac%NumOfReactWOBackward + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
  ELSE IF(NumOfCols.EQ.-1)THEN      ! old format of database
    ! Check if the non-reactive species list exists
    CALL AttributeExists(file_id_specdb,'NonReactiveSpecies',TRIM(dsetname2), AttrExists=AttrExists)
    IF(AttrExists) THEN
      ! Get the size of the attribute array
      CALL GetDataSize(dset_id,'NonReactiveSpecies',nDims,HSize,attrib=.TRUE.)
      ChemReac%ArbDiss(iReac)%NumOfNonReactives = INT(HSize(1),4)
      IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.GT.0) THEN
        ALLOCATE(ChemReac%ArbDiss(iReac)%NonReactiveSpecies(ChemReac%ArbDiss(iReac)%NumOfNonReactives))
        ALLOCATE(NonReactiveSpeciesNames(ChemReac%ArbDiss(iReac)%NumOfNonReactives))
        ! Read-in the non-reactive species names array
        CALL ReadAttribute(file_id_specdb,'NonReactiveSpecies',ChemReac%ArbDiss(iReac)%NumOfNonReactives,TRIM(dsetname2),StrArray=NonReactiveSpeciesNames)
        ! Find the index of the respective species
        DO iSpec = 1, ChemReac%ArbDiss(iReac)%NumOfNonReactives
          SpeciesFound = .FALSE.
          DO jSpec = 1, nSpecies
            IF(STRICMP(NonReactiveSpeciesNames(iSpec),Species(jSpec)%Name)) THEN
              ChemReac%ArbDiss(iReac)%NonReactiveSpecies(iSpec) = jSpec
              SpeciesFound = .TRUE.
              EXIT
            END IF
          END DO
          IF(.NOT.SpeciesFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: Species defined as non-reactive has not beend found!')
        END DO
        DEALLOCATE(NonReactiveSpeciesNames)
        ! First reaction is saved within the dummy input reaction, thus "- 1"
        ChemReac%NumOfReact = ChemReac%NumOfReact + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
        ChemReac%NumOfReactWOBackward = ChemReac%NumOfReactWOBackward + ChemReac%ArbDiss(iReac)%NumOfNonReactives - 1
      END IF
    ELSE
      ChemReac%ArbDiss(iReac)%NumOfNonReactives = 0
    END IF
  ELSE
    ChemReac%ArbDiss(iReac)%NumOfNonReactives = 0
  END IF
END DO

! Close the group and the file
CALL H5GCLOSE_F(group_id, err)
CALL H5FCLOSE_F(file_id_specdb, err)
! Close FORTRAN interface
CALL H5CLOSE_F(err)

END SUBROUTINE GetNumReacFromDatabase


SUBROUTINE ReadReacFromDatabase(ReadInNumOfReact)
!===================================================================================================================================
!> Read-in of the reactions from the species database, depending on whether the reactions are defined in the ini or by the chemistry
!> model
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_Particle_Vars     ,ONLY: SpeciesDatabase, Species, nSpecies
USE MOD_DSMC_Vars         ,ONLY: ChemReac
USE MOD_HDF5_Input        ,ONLY: DatasetExists,nDims,HSize,ReadAttribute,GetDataSize,AttributeExists
USE MOD_StringTools       ,ONLY: STRICMP
USE MOD_MCC_Vars          ,ONLY: NbrOfPhotonXsecReactions
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: ReadInNumOfReact
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=256)                :: dsetname, groupname, dsetname2
INTEGER                           :: err
INTEGER(HID_T)                    :: file_id_specdb, group_id, dset_id             ! File/group/dataset identifiers
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: ReacID
INTEGER                           :: nSets, start
LOGICAL                           :: AttrExists, GroupFound, SpeciesFound
INTEGER                           :: iReac, iVar, indexSubset, iSpec, iProd
CHARACTER(LEN=32)                 :: hilf
CHARACTER(LEN=255)                :: reacNameReactants, reacNameProducts
CHARACTER(LEN=255)                :: ReactantNames(3), ProductNames(4)
REAL                              :: SumProdMass, SumReactMass
REAL,PARAMETER                    :: RelMassTol=1e-14 ! Relative tolerance applied to conservation of mass before/after reaction
!===================================================================================================================================

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
! Check if the REACTIONS group exists
groupname = 'Reactions'
CALL H5LEXISTS_F(file_id_specdb, TRIM(groupname), GroupFound, err)
IF(.NOT.GroupFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: No reactions group found!')
! Open the group and get number of available reactions in the database
CALL H5GOPEN_F(file_id_specdb,TRIM(groupname), group_id, err)

IF(ChemReac%ChemistryModel.EQ.'none') THEN
  DO iReac = 1, ReadInNumOfReact
    ! Read-in the reaction name from the parameter file by combining the reactant and product definition
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    ChemReac%ReactModel(iReac)  = TRIM(GETSTR('DSMC-Reaction'//TRIM(hilf)//'-ReactionModel'))
    ChemReac%Reactants(iReac,:) = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Reactants',3)
    ChemReac%Products(iReac,:)  = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Products',4)
    ! Get the first reactant and write its species name to the string
    IF(ChemReac%Reactants(iReac,1).NE.0) THEN
      reacNameReactants = TRIM(Species(ChemReac%Reactants(iReac,1))%Name)
    ELSE
      CALL abort(__STAMP__,'ERROR in parameter.ini: First reactant cannot be zero!')
    END IF
    ! Add additional reactants to the string
    DO iVar = 2, 3
      IF(ChemReac%Reactants(iReac,iVar).GT.0) THEN
        reacNameReactants = TRIM(reacNameReactants)//'+'//TRIM(Species(ChemReac%Reactants(iReac,iVar))%Name)
      ELSE IF(ChemReac%Reactants(iReac,iVar).LT.0) THEN
        ! Treatment of non-reactive species (e.g. CH4 + M/A -> CH3 + H + M/A)
        IF(ABS(ChemReac%Reactants(iReac,iVar)).EQ.2) THEN
          reacNameReactants = TRIM(reacNameReactants)//'+M'
        ELSE IF (ABS(ChemReac%Reactants(iReac,iVar)).EQ.1) THEN
          reacNameReactants = TRIM(reacNameReactants)//'+A'
        ELSE
          CALL abort(__STAMP__,'ERROR in parameter.ini: Only -1 (=A) and -2 (=M) is supported as generic non-reactive species!')
        END IF
      END IF
    END DO
    ! Get the first product and write its species name to the string
    IF(ChemReac%Products(iReac,1).NE.0) THEN
      reacNameProducts = TRIM(Species(ChemReac%Products(iReac,1))%Name)
    ELSE
      CALL abort(__STAMP__,'ERROR in parameter.ini: First product cannot be zero!')
    END IF
    ! Add additional products to the string
    DO iVar = 2, 4
      IF(ChemReac%Products(iReac,iVar).GT.0) THEN
        reacNameProducts = TRIM(reacNameProducts)//'+'//TRIM(Species(ChemReac%Products(iReac,iVar))%Name)
      ELSE IF(ChemReac%Products(iReac,iVar).LT.0) THEN
        ! Treatment of non-reactive species (e.g. CH4 + M/A -> CH3 + H + M/A)
        IF(ABS(ChemReac%Products(iReac,iVar)).EQ.2) THEN
          reacNameProducts = TRIM(reacNameProducts)//'+M'
        ELSE IF (ABS(ChemReac%Products(iReac,iVar)).EQ.1) THEN
          reacNameProducts = TRIM(reacNameProducts)//'+A'
        ELSE
          CALL abort(__STAMP__,'ERROR in parameter.ini: Only -1 (=A) and -2 (=M) is supported as generic non-reactive species!')
        END IF
      END IF
    END DO
    ! Combine both strings to the reaction name to read-in from the database
    ChemReac%ReactionName(iReac) = TRIM(reacNameReactants)//'_'//TRIM(reacNameProducts)
    ! TODO: Reactions might have multiple variants #1, #2, etc.
    dsetname2 = TRIM(groupname)//'/'//TRIM(ChemReac%ReactionName(iReac))
    SELECT CASE (TRIM(ChemReac%ReactModel(iReac)))
      CASE('TCE')
        ! Total Collision Energy: Arrhenius-based chemistry model
        CALL ReadArrheniusFromDatabase(file_id_specdb,dsetname2,iReac)
      CASE('QK')
        ! Quantum Kinetic: Threshold energy based chemistry model
        ChemReac%AnyQKReaction = .TRUE.
      CASE('XSec')
        ! Chemistry model based on cross-section data
        ChemReac%AnyXSecReaction = .TRUE.
        ! CURRENTLY READ-IN DURING MCC_Chemistry_Init()
      CASE('phIon')
        ! Photo-ionization reactions
        ChemReac%CrossSection(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CrossSection')
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE('phIonXSec')
        ! Photo-ionization reactions (data read-in from database)
        NbrOfPhotonXsecReactions = NbrOfPhotonXsecReactions + 1
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE DEFAULT
        CALL abort(__STAMP__,'Selected reaction model is not supported in reaction number: ', IntInfoOpt=iReac)
    END SELECT
  END DO
ELSE
  ! Read-in of reactions part of the chemistry model
  DO iReac = 1, ReadInNumOfReact
    ReacID = ChemReac%totalReacToModel(iReac)
    ! Get name and size of name
    CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, ReacID, dsetname, err, size)
    indexSubset = INDEX(TRIM(dsetname),'#')
    dsetname2 = TRIM(groupname)//'/'//TRIM(dsetname)
    ! Remove a potential #1,#2, etc. if reaction has multiple variants
    IF(indexSubset.GT.0) THEN
      ChemReac%ReactionName(iReac) = TRIM(dsetname(1:indexSubset-1))
    ELSE
      ChemReac%ReactionName(iReac) = TRIM(dsetname)
    END IF
    CALL PrintOption('Read-in from database for reaction','INFO',StrOpt=TRIM(ChemReac%ReactionName(iReac)))
    ! Open the reaction dataset
    CALL H5DOPEN_F(file_id_specdb, dsetname2, dset_id, err)
    ! Read-in the reaction model
    CALL ReadAttribute(file_id_specdb,'ReactionModel',1,TRIM(dsetname2), StrScalar=ChemReac%ReactModel(iReac))
    ! Read-in of reactants
    CALL AttributeExists(file_id_specdb,'Reactants',TRIM(dsetname2), AttrExists=AttrExists)
    IF(AttrExists) THEN
      ! Open the reaction dataset
      CALL H5DOPEN_F(file_id_specdb, dsetname2, dset_id, err)
      ! Get the size of the attribute array
      CALL GetDataSize(dset_id,'Reactants',nDims,HSize,attrib=.TRUE.)
      nSets = INT(HSize(1),4)
      CALL ReadAttribute(dset_id,'Reactants',nSets,StrArray=ReactantNames(1:nSets))
      DO iProd = 1,nSets
        DO iSpec = 1, nSpecies
          IF(STRICMP(ReactantNames(iProd),Species(iSpec)%Name)) THEN
            ChemReac%Reactants(iReac,iProd) = iSpec
            SpeciesFound = .TRUE.
          END IF
        END DO
        IF(.NOT.SpeciesFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: Species defined in the reaction has not been found!')
      END DO
    ELSE
      CALL abort(__STAMP__,'ERROR in reaction definition: No reactants found in the selected reaction in the database')
    END IF
    ! Read-in of products
    CALL AttributeExists(file_id_specdb,'Products',TRIM(dsetname2), AttrExists=AttrExists)
    IF(AttrExists) THEN
      ! Open the reaction dataset
      CALL H5DOPEN_F(file_id_specdb, dsetname2, dset_id, err)
      ! Get the size of the attribute array
      CALL GetDataSize(dset_id,'Products',nDims,HSize,attrib=.TRUE.)
      nSets = INT(HSize(1),4)
      CALL ReadAttribute(dset_id,'Products',nSets,StrArray=ProductNames(1:nSets))
      start = 0
      DO iProd = 1,nSets
        ! If non-reactives are present, keep the second entry empty
        IF(ChemReac%ArbDiss(iReac)%NumOfNonReactives.GT.0.AND.iProd.GT.1) start = 1
        DO iSpec = 1, nSpecies
          IF(STRICMP(ProductNames(iProd),Species(iSpec)%Name)) THEN
            ChemReac%Products(iReac,start+iProd) = iSpec
            SpeciesFound = .TRUE.
          END IF
        END DO
        IF(.NOT.SpeciesFound) CALL abort(__STAMP__,'ERROR in SpeciesDatabase: Species defined in the reaction has not been found!')
      END DO
    ELSE
      CALL abort(__STAMP__,'ERROR in reaction definition: No products found in the selected reaction in the database')
    END IF
    ! Get mass of products and reactants
    SumProdMass = 0
    DO iSpec = 1, 4
      IF(ChemReac%Products(iReac,iSpec).NE.0) SumProdMass = SumProdMass + Species(ChemReac%Products(iReac,iSpec))%MassIC
    END DO
    SumReactMass = 0
    DO iSpec = 1, 3
      IF(ChemReac%Reactants(iReac,iSpec).NE.0) SumReactMass = SumReactMass + Species(ChemReac%Reactants(iReac,iSpec))%MassIC
    END DO
    ! Sanity mass check for reactions of selected chemistry model, real compare with RelMassTol
    IF(.NOT.ALMOSTEQUALRELATIVE(SumProdMass,SumReactMass,RelMassTol)) THEN
      CALL PrintOption('DSMC_Chemistry might not be mass conserving for chemical reaction:','WARNING',StrOpt=TRIM(ChemReac%ReactionName(iReac)))
      CALL abort(__STAMP__,'DSMC_Chemistry might not be mass conserving for current chemical reaction!')
    END IF
      ! Read-in of the reaction parameters, depending on the model
    SELECT CASE (TRIM(ChemReac%ReactModel(iReac)))
      CASE('TCE')
        ! Total Collision Energy: Arrhenius-based chemistry model
        CALL ReadArrheniusFromDatabase(file_id_specdb,dsetname2,iReac)
      CASE('QK')
        ! Quantum Kinetic: Threshold energy based chemistry model
        ! Only requires the species-specific dissociation energy
        ChemReac%AnyQKReaction = .TRUE.
      CASE('XSec')
        ! CURRENTLY READ-IN DURING MCC_Chemistry_Init()
        ! Chemistry model based on cross-section data
        ChemReac%AnyXSecReaction = .TRUE.
      CASE('phIon')
        ! Photo-ionization reactions
        ChemReac%CrossSection(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CrossSection')
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE('phIonXSec')
        ! Photo-ionization reactions (data read-in from database)
        NbrOfPhotonXsecReactions = NbrOfPhotonXsecReactions + 1
        ChemReac%AnyPhIonReaction = .TRUE.
      CASE DEFAULT
        CALL abort(__STAMP__,'Selected reaction model as part of a chemistry model is not supported in reaction number: ', IntInfoOpt=iReac)
    END SELECT
  END DO
END IF

! Close the group and the file
CALL H5GCLOSE_F(group_id, err)
CALL H5FCLOSE_F(file_id_specdb, err)
! Close FORTRAN interface
CALL H5CLOSE_F(err)

END SUBROUTINE ReadReacFromDatabase


SUBROUTINE ReadArrheniusFromDatabase(file_id_specdb,dsetname2,iReac)
!===================================================================================================================================
!> Read-in of the Arrhenius parameters from the database
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: ChemReac
USE MOD_HDF5_Input        ,ONLY: DatasetExists,ReadAttribute
USE MOD_ReadInTools       ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T),INTENT(IN)         :: file_id_specdb
CHARACTER(LEN=256),INTENT(IN)     :: dsetname2
INTEGER,INTENT(IN)                :: iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: ReactionFound
!===================================================================================================================================

CALL DatasetExists(file_id_specdb,TRIM(dsetname2),ReactionFound)
IF(.NOT.ReactionFound) CALL abort(__STAMP__,'ERROR in parameter.ini: Defined reaction has not been found in the database!')

CALL ReadAttribute(file_id_specdb,'Arrhenius-Prefactor',1,DatasetName = dsetname2,RealScalar=ChemReac%Arrhenius_Prefactor(iReac))
CALL PrintOption('Arrhenius-Prefactor','DB',RealOpt=ChemReac%Arrhenius_Prefactor(iReac))
CALL ReadAttribute(file_id_specdb,'Arrhenius-Powerfactor',1,DatasetName = dsetname2,RealScalar=ChemReac%Arrhenius_Powerfactor(iReac))
CALL PrintOption('Arrhenius-Powerfactor','DB',RealOpt=ChemReac%Arrhenius_Powerfactor(iReac))
CALL ReadAttribute(file_id_specdb,'Activation-Energy_K',1,DatasetName = dsetname2,RealScalar=ChemReac%EActiv(iReac))
ChemReac%EActiv(iReac) = ChemReac%EActiv(iReac)*BoltzmannConst
CALL PrintOption('Activation-Energy_K','DB',RealOpt=ChemReac%EActiv(iReac))

END SUBROUTINE ReadArrheniusFromDatabase


!===================================================================================================================================
!> Calculation and sanity check of ChemReac%EForm(iReac) for photoionization reactions, i.e.,  ChemReac%ReactModel(iReac).EQ.'phIon'
!> The sanity check will determine if the photon energy is sufficient to trigger any reactions.
!===================================================================================================================================
SUBROUTINE InitPhotonReactions()
! MODULES
USE MOD_Globals             ,ONLY: abort
USE MOD_DSMC_Vars           ,ONLY: ChemReac,CollisMode,UseDSMC
USE MOD_part_emission_tools ,ONLY: CalcPhotonEnergy
USE MOD_PARTICLE_Vars       ,ONLY: nSpecies
USE MOD_RayTracing_Vars     ,ONLY: RayPartBound,Ray
USE MOD_Particle_Vars       ,ONLY: Species
USE MOD_MCC_Vars            ,ONLY: NbrOfPhotonXsecReactions
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: PhotonEnergy
INTEGER :: iInit, iSpec, iReac
!===================================================================================================================================
IF(.NOT.UseDSMC) RETURN
IF(CollisMode.NE.3) RETURN

DO iReac = 1, ChemReac%NumOfReact
  ! Check whether the photon energy is sufficient to trigger the chemical reaction
  IF(TRIM(ChemReac%ReactModel(iReac)).EQ.'phIon') THEN
    PhotonEnergy = 0.

    ! Const. photon energy from wavelength given in ray tracing model
    IF(RayPartBound.GT.0)THEN
      PhotonEnergy = CalcPhotonEnergy(Ray%WaveLength)
    END IF ! RayPartBound.GT.0

    ! Const. photon energy from wavelength given in particle emission model
    DO iSpec = 1, nSpecies
      DO iInit = 1, Species(iSpec)%NumberOfInits
        SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
        CASE('photon_cylinder','photon_honeycomb','photon_rectangle')
          PhotonEnergy = CalcPhotonEnergy(Species(iSpec)%Init(iInit)%WaveLength)
          EXIT
        END SELECT
      END DO
    END DO

    ChemReac%EForm(iReac) = ChemReac%EForm(iReac) + PhotonEnergy
    IF(ChemReac%EForm(iReac).LE.0.0) THEN
      CALL abort(__STAMP__,'ERROR: Photon energy is not sufficient for the given ionization reaction: ',iReac)
    END IF
    ! Abort if photon-ionization reactions using cross-sections have been defined
    IF(NbrOfPhotonXsecReactions.GT.0) CALL abort(__STAMP__,&
      'Photoionization reactions with constant cross-sections cannot be combined with XSec data cross-sections for photoionization')
  END IF ! TRIM(ChemReac%ReactModel(iReac)).EQ.'phIon'
END DO ! iReac = 1, ChemReac%NumOfReact

END SUBROUTINE InitPhotonReactions


END MODULE MOD_DSMC_ChemInit
