(sec:unified-species-database)=
# Unified Species Database

A unified database of species data, electronic states, cross-sections, and chemistry models can be used as a more convenient alternative input for the simulations. The use of the database allows to reduce the length of the input files and ensures a consistent storage of the parameters, together with the respective reference. The database (`SpeciesDatabase.h5`) is provided at the top level of the PICLas directory. It contains over 40 of the most common species in PICLas simulations and the most common complex chemical reaction models with 300 single reactions. In addition, several cross-section models are included. This feature is still under development and the provided species & reaction data should be treated carefully as we are in process of verifying the data.

The data in the Unified Species Database is grouped, as shown in the following example:

    Cross-Sections (group)
        H2-H2Ion1 (dataset)
    Reaction (group)
        CH3_CH2+H (dataset)
            Chemistry model (attribute)
            Reaction model (attribute)
            Arrhenius parameters (attribute)
            Products    (attribute)
            Reactants   (attribute)
        O2+M_O+O+M (dataset)
            Chemistry model (attribute)
            Reaction model (attribute)
            Arrhenius parameters (attribute)
            Products    (attribute)
            Reactants   (attribute)
        Fe_FeIon1+electron (dataset)
            Chemistry model (attribute)
            Reaction model (attribute)
            Products    (attribute)
            Reactants   (attribute)
    Species (group)
        H2 (group)
            Electronic levels (dataset)
            Species parameters (attribute)
        H2Ion1 (group)
            Electronic levels (dataset)
            Species parameters (attribute)
        electron (group)
            Electronic levels (dataset)
            Species parameters (attribute)

To read in data from a given species database, the database name must be specified in the `parameter.ini` file of the simulation

    Particles-Species-Database = SpeciesDatabase.h5

When reading the database, all available parameters are taken directly from the database by default. Missing parameters can be added as usual in the `parameter.ini` or the `DSMC.ini` files. If no database is specified or the specified database is not available, the parameters are read from the regular parameter files as described in the sections above. All data available in the database and how to choose between different input forms is described in the following sections for the different data types available. An example where the given database is used can be found in `regressioncheck/NIG_Reservoir/CHEM_EQUI_TCE_Air_5Spec_Database`.

For instructions on extending the existing database or creating a new one, please refer to Chapter {ref}`sec:tools-usd` for instructions.

(ssec:usd-species)=
## Species data

To include a species in a simulation, it has to be selected by its name with the following command:

    Part-Species1-SpeciesName = CH4

The database follows the PICLas nomenclature for atoms and molecules. Cations are given by 

    Part-Species4-SpeciesName = CH3Ion1

where the number after 'Ion' refers to the degree of ionization. The database contains general species data such as the mass, charge and number of atoms, as well as the VHS and VSS parameters. For molecules, rotational and vibrational frequencies can be given as well. Additionally, each species contains a dataset with the electronic energy levels. A complete list of the parameter and models is available in Section {ref}`sec:DSMC-species`. Per default, all available parameters are read from the database. While we are in the process of verifying the database, the table gives an overview over the already verified parameters and the respective data sources. In general, the source/reference will also be provided with the parameter in the database itself.

|                                    Parameter |    Unit    | Source (unless otherwise noted)                                                                  |
| -------------------------------------------: | :--------: | :----------------------------------------------------------------------------------------------- |
|               Heat of formation (at 298.15K) |   Kelvin   | [Active Thermochemical Tables (ATcT)](https://atct.anl.gov/)                                     |

<!-- | Electronic energy levels for atoms (Degeneracy,Energy) | - , Kelvin | [Atoms: NIST Atomic Spectra Database](https://physics.nist.gov/PhysRefData/ASD/levels_form.html) | -->


It possible to revert to the regular parameter read-in from parameter files per species

    Part-Species1-DoOverwriteParameters = true

If this flag is set, all parameters for this species need to be set manually and a separate electronic state database has to be provided:

    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5


(ssec:usd-reaction)=
## Reaction data

The database contains different chemistry models including various reactions, utilizing the {ref}`ssec:TCE` model or {ref}`ssec:QK` model. Each reaction is assigned to a single or more chemistry models (N). This information is stored in the chemistry model attribute along with the non reactive species of this reaction. This attribute contains an array which either has the dimension (N,1) or (N,2), where the second column contains the non reactive species, if existent.
For the TCE model, the Arrhenius prefactor, powerfactor and the activation energy are stored. For the QK model only dissociation energy of the molecule is required, which has been provided as a species parameter. 


To utilize the read-in of reaction from the database, a chemistry model has to be defined. All reactions with this model are then read-in from the database and no additional parameter input is required:

    DSMC-ChemistryModel = Titan_14Spec_24Reac_Gokcen2007

|                    Name                    |                Description                | Species | React. | Reference                                                                                                                                                                                                                                          |
| :----------------------------------------: | :---------------------------------------: | :-----: | :----: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|          Air_5Spec_8Reac_Laux1995          |             Air without ions              |    5    |   8    | M. Laux. *Direkte Simulation verdünnter, reagierender Strömungen.* PhD Thesis, University of Stuttgart, 1995.                                                                                                                                      |
|          Air_5Spec_8Reac_Park1993          |             Air without ions              |    8    |   8    | C. Park. *Review of chemical-kinetic problems of future NASA missions. I - Earth entries.* [Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993.](https://doi.org/10.2514/3.431)                                                        |
|         Air_11Spec_27Reac_Park1993         |               Air with ions               |   11    |   27   | C. Park. *Review of chemical-kinetic problems of future NASA missions. I - Earth entries.* [Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993.](https://doi.org/10.2514/3.431)                                                        |
|        Air_11Spec_43Reac_Farbar2008        | Air with ions and backward reaction rates |   11    |   43   | E. Farbar and I. Boyd. *Simulation of Fire II Reentry Flow Using the Direct Simulation Monte Carlo Method.* [40th Thermophysics Conference, 2008.](https://doi.org/10.2514/6.2008-4103)                                                            |
|         Air_11Spec_51Reac_Park1993         |               Air with ions               |   11    |   51   | C. Park. *Review of chemical-kinetic problems of future NASA missions. I - Earth entries.* [Journal of Thermophysics and Heat Transfer, 7(3):385–398, 1993.](https://doi.org/10.2514/3.431)                                                        |
| Fe-in-Air_3Spec_7Reac_Voronov1997Plane2015 |        Outgassing of iron into air        |    3    |   7    | Voronov. G. *A practical fit formula for ionization rate coefficients of atoms and ions by electron impact: Z= 1–28.* [Atomic Data and Nuclear Data Tables, 65(1):1–35, 1997.](https://doi.org/10.1006/adnd.1997.0732)                             |
|                                            |                                           |         |        | J. M. Plane, W. Feng, and E. C. Dawkins. *The Mesosphere and Metals: Chemistry and Changes.* [Chemical Reviews, 115(10):4497–4541, 2015.](https://doi.org/10.1021/cr500501m)                                                                       |
|              CH4_7Spec_7Reac               |                                           |    7    |   7    | -                                                                                                                                                                                                                                                  |
|             CH4-Ar_8Spec_7Reac             |                                           |    8    |   7    | -                                                                                                                                                                                                                                                  |
|       CO2_6Spec_10Reac_Johnston2014        |                                           |    6    |   10   | C. Johnston and A. Brandis. *Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases.* [Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014.](https://doi.org/10.1016/j.jqsrt.2014.08.025) |
|      Mars_11Spec_27Reac_Johnston2014       |             Mars without ions             |   11    |   27   | C. Johnston and A. Brandis. *Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases.* [Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014.](https://doi.org/10.1016/j.jqsrt.2014.08.025) |
|        Mars_16Spec_31Reac_Park1994         |              Mars with ions               |   16    |   31   | C. Park, J. T. Howe, R. L. Jaffe, and G. V. Candler. *Review of chemical-kinetic problems of future NASA missions. II - Mars entries.* [Journal of Thermophysics and Heat Transfer, 8(1):9–23, 1994.](https://doi.org/10.2514/3.496)               |
|      Mars_17Spec_42Reac_Johnston2014       |          Mars with ions and O2+           |   17    |   42   | C. Johnston and A. Brandis. *Modeling of nonequilibrium CO Fourth-Positive and CN Violet emission in CO2–N2 gases.* [Journal of Quantitative Spectroscopy and Radiative Transfer, 149:303–317, 2014.](https://doi.org/10.1016/j.jqsrt.2014.08.025) |
|       Titan_14Spec_24Reac_Gokcen2007       |     Titan without ions but with Argon     |   14    |   24   | R. Savajano, R. Sobbia, M. Gaffuri, and P. Leyland. *Reduced Chemical Kinetic Model for Titan Entries.* [International Journal of Chemical Engineering, vol. 2011, Article ID 970247, 2011.](https://doi.org/10.1155/2011/970247)                  |
|                                            |                                           |         |        | T. Gokcen. *N2-CH4-Ar Chemical Kinetic Model for Simulations of Atmospheric Entry to Titan.* [Journal of Thermophysics and Heat Transfer, 21(1):9–18, 2007.](https://doi.org/10.2514/1.22095)                                                      |
|       Titan_18Spec_30Reac_Gokcen2007       |     Titan with ions but without Argon     |   18    |   30   | T. Gokcen. *N2-CH4-Ar Chemical Kinetic Model for Simulations of Atmospheric Entry to Titan.* [Journal of Thermophysics and Heat Transfer, 21(1):9–18, 2007.](https://doi.org/10.2514/1.22095)                                                      |

Currently, these reactions cannot be utilized separately, but only as part of a chemistry model.

<!-- Reactions to be included in the simulation are specified by their reaction equation or their chemical model:

    ! Reaction1: CH4 + M -> CH3 + H + M 
    DSMC-Reaction1-Reactants = (/1,0,0/)
    DSMC-Reaction1-Products = (/2,0,3,0/)

A reaction name is generated automatically based on the given species names and read-in from the database. To ensure a correct read-in, all species in the chosen reactions must have a defined species name in the parameter.ini Reactants are ordered according to a predefined list, with nonreacting partners listed always at the end. The same general order is used for the products, however the nonreacting partners are given always at the second position. If one reaction appears in multiple models or with multiple parameter sets in the database an additional enumerator in the form of e.g. '#5' is given at the end of the reaction name in the database. -->

(ssec:usd-xsec-data)=
## Cross-section data

The use of the unified species database for the cross-section data follows the description given in Section {ref}`sec:background-gas-collision-xsec` for collision and {ref}`sec:xsec-chemistry` for chemistry modelling, respectively. An example for is located in `regressioncheck/NIG_Reservoir/CHEM_RATES_XSec_Chem_H2_Plasma_Database`.
