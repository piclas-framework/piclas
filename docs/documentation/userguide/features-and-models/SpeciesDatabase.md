(sec:SpeciesDatabase)=
# Unified Species Database

A unified database of species data, electronic states, cross-sections, and chemistry models can be used as a more convenient alternative input for the simulations. The use of the database allows to reduce the length of the input files and ensures a meaningful storage of the parameters, together with the reference from which they are taken.

An predefined database (`SpeciesDatabase.h5`) is provided at the top level of the PICLas directory. This database contains over 40 of the most common species in PICLas simulations and the most common complex chemical reaction models with 300 single reactions. In addition, several cross-section models are included.

To read in data from a given species database, the database name must be specified in the `parameter.ini` file of the simulation

    Particles-Species-Database = SpeciesDatabase.h5

When reading the database, all available parameters are taken directly from the database by default. Missing parameters can be added as usual in the `parameter.ini` or the `DSMC.ini` files. If no database is specified or the specified database is not available, the parameters are read from the `parameter.ini` and the `DSMC.ini` files as described in the sections above.

All data available in the database and how to choose between different input forms is described in the following sections for the different data types available. An example where the given database is used can be found in `regressioncheck/NIG_Reservoir/CHEM_EQUI_TCE_Air_5Spec_Database`.

For instructions on extending the existing database or creating a new one, please refer to Chapter {ref}`sec:tools-usd` for instructions.

(ssec:Species)=
## Species-specific data

To include a species in a simulation, it has to be selected by its name with the following command:

    Part-Species1-SpeciesName = CH4

The database follows the standard nomenclature for atoms and molecules. Cations are given by 

    Part-Species4-SpeciesName = CH3Ion1

where the number after 'Ion' refers to the degree of ionization. The database contains general species data such as the mass, charge and number of atoms, as well as with VHS and VSS parameters. For molecules, rotational and vibrational frequencies can be given as well. A complete list of available parameters and the corresponding variable names in the `DSMC.ini` input are given below. 

    Part-Species1-SpeciesName
    Part-Species1-MassIC
    Part-Species1-ChargeIC
    Part-Species1-PolyatomicMol
    Part-Species1-InteractionID
    Part-Species1-NumOfAtoms
    Part-Species1-LinearMolec
    Part-Species1-SymmetryFactor 
    Part-Species1-Tref
    Part-Species1-dref
    Part-Species1-omega
    Part-Species1-alphaVSS
    Part-Species1-CharaTempVibX 
    Part-Species1-CharaTempRotX
    Part-Species1-Ediss_eV
    Part-Species1-HeatOfFormation_K 
    Part-Species1-PreviousState

Per default, all available parameters are read from the database. Species for which the parameters should be taken from the input files can be specified by an overwrite command. 

    Part-Species1-DoOverwriteParameters = true

If this flag is set, all parameters for this species need to be set manually.

(ssec:Reaction)=
## Reaction data

The database contains different chemistry models including various reactions.

| Name                                       | Description | Species | React. | Reference |
| :----------------------------------------: |  :------------:  | :---------: |  :------: | :----: |
| Air_5Spec_8Reac_Laux1995                   | 2.3.0 (Nov 2021) |  5 |  8 | 1.12.1 |
| Air_5Spec_8Reac_Park1993                   | 2.2.0 (Nov 2021) |  8 |  8 | 1.10.5 |
| Air_11Spec_27Reac_Park1993                 | 2.1.0 (Nov 2021) | 11 | 27 | 1.10.6 |
| Air_11Spec_43Reac_Farbar2008               | 2.0.0 (Nov 2021) | 11 | 43 | 1.10.5 |
| Air_11Spec_51Reac_Park1993                 | 2.0.0 (Nov 2021) | 11 | 51 | 1.10.5 |
| Fe-in-Air_3Spec_7Reac_Voronov1997Plane2015 | 2.0.0 (Nov 2021) |  3 |  7 | 1.10.5 |
| CH4_7Spec_7Reac                            | 2.0.0 (Nov 2021) |  7 |  7 | 1.10.5 |
| CH4-Ar_8Spec_7Reac                         | 2.3.0 (Nov 2021) |  8 |  7 | 1.12.1 |
| CO2_6Spec_10Reac_Johnston2014              | 2.2.0 (Nov 2021) |  6 | 10 | 1.10.5 |
| Mars_11Spec_27Reac_Johnston2014            | 2.3.0 (Nov 2021) | 11 | 27 | 1.12.1 |
| Mars_16Spec_31Reac_Park1994                | 2.2.0 (Nov 2021) | 16 | 31 | 1.10.5 |
| Mars_17Spec_42Reac_Johnston2014            | 2.3.0 (Nov 2021) | 17 | 42 | 1.12.1 |
| Titan_14Spec_24Reac_Gokcen2007             | 2.2.0 (Nov 2021) | 14 | 24 | 1.10.5 |
| Titan_18Spec_30Reac_Gokcen2007             | 2.3.0 (Nov 2021) | 18 | 30 | 1.12.1 |

The database contains data for the TCE and QK model. Reactions to be included in the simualtion are specified by their name:

    ! Reaction1: CH4 + M -> CH3 + H + M 
    DSMC-Reaction1-ReactionName = DISS_CH4_TO_CH3
    
To ensure consistency and findability, a naming convention is used in the database. The available flags for the reaction are `DISS` (dissociation), `EXC` (exchange) and `RECOMB` (recombination or addition). The main reactants and products are specified as well. To further state atomic, electronic or molecular collision parameters, as well as ionization processes, the flags `_AT`, `_MOL`, `_EL` and `_ION` are used respectively. The name of reactants, products and non-reactive collision parameters are deposited in the database as well, but they are not read out to avoid restrictions in the simulation setup. These values are given manually in the `parameter.ini`, as described in Section {ref} `sec:DSMC-chemistry`.

If the reaction parameters should be given manually, the following command can be set:

    DSMC-OverwriteReacDatabase = true

(sssec:Arrh-rates)=
### Arrhenius rates

Arrhenius-type reaction rates needed for the Total-Collision-Energy model ({ref}`ssec:TCE`) can be taken directly from the database. Here, the Arrhenius prefactor and powerfactor, together with the activation energy is stored. 

(sssec:QK-model)=
### QK model

The dissociation energy used in the Quantum-Kinetic model ({ref}`ssec:QK`) is deposited in the species database.

(ssec:Catalysis)=
## Heterogenous reactions

All data necessary for the modelling of gas-surface reactions with PICLas can be found in the species database. To include parameters for a certain reaction, the reaction name needs to be specified by
     
    ! Reaction 1: CO(g) -> CO(ads) 
    Surface-Reaction1-SurfName = Adsorption_CO
    
According to the naming convention used in the database, the surface mechanism (adsorption, desorption, ER and LH reaction) is stated together with the main reactant and product (only in the case of reactions). The database contains the sticking coefficients for the Kisluik and Langmuir adsorption model, together with the Polanyi-Wigner parameters and the Arrhenius-type reaction rates. All parameters that can be found in the database are givel below. 

    Surface-Reaction1-SurfName 
    Surface-Reaction1-Type 
    Surface-Reaction1-StickingCoefficient 
    Surface-Reaction1-EqConstant
    Surface-Reaction1-DissOrder 
    Surface-Reaction1-ReactHeat 
    Surface-Reaction1-HeatScaling 
    Surface-Reaction1-Energy 
    Surface-Reaction1-Prefactor
    Surface-Reaction1-LateralInteraction 
    Surface-Reaction1-Ca 
    Surface-Reaction1-Cb

The reactants and products are deposited in the database as well, but again not read out to ensure the correct input. Reactants and products need to be defined by 

    Surface-Reaction1-Reactants = (/2,2/) 
    Surface-Reaction1-Products = (/3,0/)
    
So far, the database does not contain any surface-specific or diffusion parameters. These need to be set in the `parameter.ini` by

    Surface-Reaction1-NumOfBoundaries = 1
    Surface-Reaction1-Boundaries = (/1,0/)
    Part-Boundary1-WallTemp = 323
    Part-Boundary1-LatticeVector = 0.389E-9
    Part-Boundary1-NbrOfMol-UnitCell = 2
    Part-Boundary1-Species1-Coverage = 0.1
    Part-Boundary1-Species1-MaxCoverage = 0.333 
    
If a parameter is not defined or set to a default value, the code will abort.

To define parameters manually instead of using the database, the following option can be selected:

    OverwriteCatParameters = true

(ssec:El-states)=
## Electronic states

The modelling of electronic relaxation follows the principle as described in {ref} `sec:DSMC-electronic-relaxation`.  As for the species data, the electronic levels can be given separately by specifiying 

    Part-Species1-DoOverwriteParameters = true
    
In this case, the electronic relaxation data is taken from an additional Electronic-Database:

    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5
 

(ssec:Xsec-data)=
## Cross-section data

The use of the unififed species database for the cross-section data, follows the description given in Section {ref} `ssec:xsec-chemistry`. All reaction paths are again stored by their reaction names and can be called in the `parameter.ini`.

(ssec:Rad-data)=
## Radiation data

TO-DO

(ssec:Overview)=
## Overview

An overview of the functionalities of the unified species database is given in the scheme below:

TO-DO


