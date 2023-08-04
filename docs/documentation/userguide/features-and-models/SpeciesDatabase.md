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
| Air_5Spec_8Reac_Laux1995                   | Air without ions |  5 |  8 | M. Laux. Direkte Simulation verdünnter, reagierender Strömungen. PhD Thesis, University of Stuttgart, 1995. |
| Air_5Spec_8Reac_Park1993                   | Air without ions |  8 |  8 | C. Park. Review of chemical-kinetic problems of future NASA missions. I - Earth entries. Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993. doi:10.2514/3.431. |
| Air_11Spec_27Reac_Park1993                 | Air with ions | 11 | 27 | C. Park. Review of chemical-kinetic problems of future NASA missions. I - Earth entries. Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993. doi:10.2514/3.431. |
| Air_11Spec_43Reac_Farbar2008               | Air with ions and backward reaction rates| 11 | 43 | E. Farbar and I. Boyd. Simulation of Fire II Reentry Flow Using the Direct Simulation Monte Carlo Method. 2008. doi:10.2514/6.2008-4103. |
| Air_11Spec_51Reac_Park1993                 | Air with ions | 11 | 51 | C. Park. Review of chemical-kinetic problems of future NASA missions. I - Earth entries. Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993. doi:10.2514/3.431. |
| Fe-in-Air_3Spec_7Reac_Voronov1997Plane2015 | Outgassing of iron into air |  3 |  7 | Voronov. G. A practical fit formula for ionization rate coefficients of atoms and ions by electron impact: Z= 1–28. Atomic Data and Nuclear Data Tables, 65(1):1–35, 1997. + J. M. Plane, W. Feng, and E. C. Dawkins. The mesosphere and metals: Chemistry and changes. Chemical reviews, 115(10):4497–4541, 2015. |
| CH4_7Spec_7Reac                            |  |  7 |  7 | - |
| CH4-Ar_8Spec_7Reac                         |  |  8 |  7 | - |
| CO2_6Spec_10Reac_Johnston2014              |  |  6 | 10 | C. Johnston and A. Brandis. Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases. Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014. doi:10.1016/j.jqsrt.2014.08.025. |
| Mars_11Spec_27Reac_Johnston2014            | Mars without ions | 11 | 27 | C. Johnston and A. Brandis. Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases. Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014. doi:10.1016/j.jqsrt.2014.08.025. |
| Mars_16Spec_31Reac_Park1994                | Mars with ions | 16 | 31 | C. Park, J. T. Howe, R. L. Jaffe, and G. V. Candler. Review of chemical-kinetic problems of future NASA missions. II - Mars entries. Journal of Thermophysics and Heat Transfer, 8(1):9–23, 1994. doi:10.2514/3.496. |
| Mars_17Spec_42Reac_Johnston2014            | Mars with ions and O2+ | 17 | 42 | C. Johnston and A. Brandis. Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases. Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014. doi:10.1016/j.jqsrt.2014.08.025. |
| Titan_14Spec_24Reac_Savajano2011           | Titan without ions but with Argon | 14 | 24 | R. Savajano, R. Sobbia, M. Gaffuri, and P. Leyland. Reduced Chemical Kinetic Model for Titan Entries, International Journal of Chemical Engineering, vol. 2011, Article ID 970247, 2011. doi:10.1155/2011/970247. |
| Titan_18Spec_30Reac_Savajano2011           | Titan with ions but without Argon | 18 | 30 | R. Savajano, R. Sobbia, M. Gaffuri, and P. Leyland. Reduced Chemical Kinetic Model for Titan Entries, International Journal of Chemical Engineering, vol. 2011, Article ID 970247, 2011. doi:10.1155/2011/970247. |

he database contains data for the TCE and QK model. Reactions to be included in the simulation are specified by their reaction equation or their chemical model:

    ! Reaction1: CH4 + M -> CH3 + H + M 
    DSMC-Reaction1-Reactants = (/1,0,0/)
    DSMC-Reaction1-Products = (/2,0,3,0/)
    DSMC-Reaction1-ReactionName = CH4+M_CH3_M+H
    
The reaction name is generated automatically and follows a set convention, that is enforced inside the database. To ensure a correct read-in, all species in the chosen reactions must have a defined species name in the parameter.ini Reactants are ordered according to a predefined list, with nonreacting partners listed always at the end. The same general order is used for the products, however the nonreacting partners are given always at the second position. If one reaction appears in multiple models or with multiple parameter sets in the database an additional enumerator in the form of f.e. '#5' is given at the end of the reaction name. In these cases, the correct number needs to be supplied in the parameter.ini as well.

If a chemistry model is defined, all reactions with this model are read-in from te database and no additional reaction names need to be supplied. The use of a set of reaction equations from the database can be initialized with the following command:

    ! Reaction set 1
    DSMC-ChemistryModel = Titan_14Spec_24Reac_Savajano2011
    
If the reaction parameters should be given manually in the parameter.ini, the following command can be set: 

    DSMC-OverwriteReacDatabase = true

(sssec:Arrh-rates)=
### Arrhenius rates

Arrhenius-type reaction rates needed for the Total-Collision-Energy model ({ref}`ssec:TCE`) can be taken directly from the database. Here, the Arrhenius prefactor and powerfactor, together with the activation energy is stored. 

(sssec:QK-model)=
### QK model

The dissociation energy used in the Quantum-Kinetic model ({ref}`ssec:QK`) is deposited in the species database.

(ssec:El-states)=
## Electronic states

The modelling of electronic relaxation follows the principle as described in {ref} `sec:DSMC-electronic-relaxation`.  As for the species data, the electronic levels can be given separately by specifiying 

    Part-Species1-DoOverwriteParameters = true
    
In this case, the electronic relaxation data is taken from an additional Electronic-Database:

    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5
 

(ssec:Xsec-data)=
## Cross-section data

The use of the unififed species database for the cross-section data, follows the description given in Section {ref} `ssec:xsec-chemistry`. All reaction paths are again stored by their reaction names and can be called in the `parameter.ini`.

(ssec:Overview)=
## Overview

An overview of the functionalities of the unified species database is given in the scheme below:

TO-DO


