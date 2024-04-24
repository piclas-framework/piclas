# Appendix

## Tested compiler combinations

The following list summarizes all **tested combinations** of the required libraries (HDF5, OpenMPI, CMake etc.)

| Dev |  PICLas Version  |    System   |  Compiler |  HDF5  |      MPI      |   CMake  |
| --- |  :------------:  | :---------: |  :------: | :----: | :-----------: | :------: |
|  SC | 2.3.0 (Nov 2021) |      PC     | gcc11.2.0 | 1.12.1 | openmpi-4.1.1 |  3.21.3  |
|  SC | 2.2.0 (Nov 2021) |      PC     | gcc10.1.0 | 1.10.5 | openmpi-4.0.2 |  3.17.0  |
|  AM | 2.1.0 (Nov 2021) |      PC     | gcc9.3.0  | 1.10.6 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) |  boltzhawk  |  gcc9.3.0 | 1.10.5 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) | boltzreggie |  gcc9.2.0 | 1.10.5 | openmpi-4.0.2 |  3.15.0  |
|  SC | 2.0.0 (Nov 2021) |  boltzplatz |  gcc9.2.0 | 1.10.5 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) |     hawk    |  gcc9.2.0 | 1.10.5 |    mpt2.23    |  3.16.4  |
|  SC | 2.0.0 (Nov 2021) |     fh1     | intel18.1 |  1.10  |    impi2018   |   3.17   |
|  SC | 2.0.0 (Nov 2021) |     fh2     | intel19.1 |  1.10  |    impi2019   |   3.17   |
|  PN |  1.4.0 (Nov 19)  |  boltzplatz |  gnu7.4.0 | 1.10.5 | openmpi-3.1.3 | 3.15.3-d |
|  SC |  1.4.0 (Nov 19)  | boltzreggie |  gnu9.2.0 | 1.10.5 | openmpi-4.0.2 | 3.15.3-d |

Combinations that can cause problems are listed in the following table

| Dev | PICLas Version |    System   | Compiler |  HDF5  |      MPI      |   CMake  |                                            Notes                                            |
| --- | :------------: | :---------: | :------: | :----: | :-----------: | :------: |          :-----------------------------------------------------------------------:          |
|  SC | 1.4.0 (Nov 19) | boltzreggie | gnu9.2.0 | 1.10.5 | openmpi-4.0.1 | 3.15.3-d | Does not work for more than 3 processors probably due to a problem with the OpenMPI version |
|     |                |             |          |        |               |          |                                                                                             |

## Unified Species Database (USD)

This section documents the initial creation process for the database, but since it is updated all the time it is recommended that modifications and maintenance are performed using the scripts in the *tools* folder: `piclas/tools/species_database/`. These scripts are meant to provide a controlled and uniform approach to database management.
All files to create the original database were moved to `piclas/tools/archive/`.


A tool to create a database containing cross-section, electronic states, Arrhenius rates, and species data can be found in the *tools* folder: `piclas/tools/species_database/`.
The Python script (python3.7) `create_species_database.py` creates a new database or expands an existing one combining all necessary parameters (formerly read-in through ini files). The script uses the `numpy`, `h5py`, `argparse`,`datetime`, `cmath`, and `matplotlib.rcsetup` packages. To create the species database run the command:

    python3.7 create_species_database.py
    
If electronic states or cross-section data should be added to the species database, an electronic states `Electronic-State-Database.h5` and a cross-section database `XSec-Database.h5` need to be built before.

If nothing additionally is specified, the following filenames are called: `DSMC.ini` for the parameter, gas-phase reaction input and `Species_Database.h5` for the final output. For custom file names and for the inclusion of electronic and cross-section data, the following options can be added: 

    python3 create_species_database.py --parameter parameter-filename --electronic electronic_statefile --crosssection crosssection_statefile --output output_filename --reference reference-name

or

    python3 create_species_database.py -p parameter-filename -e electronic_statefile -c crosssection_statefile -o output_filename -r reference-name
    
The data is grouped in the output file, as shown in the following example:
    
    Cross-Sections (group)
        H2-H2Ion1 (dataset)
    Reaction (group)
        CH3_CH2+H (dataset)
        Chemistry model (attributes)
        Arrhenius parameters (attributes)
        O2+M_O+O+M (dataset)
        Chemistry model (attributes)
        Arrhenius parameters (attributes)
        Fe_FeIon1+electron (dataset)
        Chemistry model (attributes)
        Arrhenius parameters (attributes)
    Species (group)
        H2 (group)
            Electronic levels (dataset)
            Species parameters (attributes)
        H2Ion1 (group)
            Electronic levels (dataset)
            Species parameters (attributes)
        electron (group)
            Electronic levels (dataset)
            Species parameters (attributes)
        
For cross-sections, reactions and species data, the former `DSMC.ini` files are used to create the database. However, every species must be defined to create the database and to run simulations. 

    Part-Species1-SpeciesName=CO2
     
The name of the reaction is optional, if none is given, the reaction name is created automatically from the reactants and products, following the naming convention defined below:

    Reac1+Reac2_Prod1+Prod2+Prod3 
    Reac1+Reac2_Prod1+Prod2+Prod3#2 
    
Non-reacting partners can be given as A (atoms or similar) or M (molecules or similar) in the reaction name. If a name is defined in the input file, the programm automatically renames it according to the set convention. If multiple sets of parameters or multiple models exist for the same reaction, a counter variable f.e. '#5' is added at the end of the reaction name.

For reactions a chemistry model is defined in all cases. The name of the given parameter-file is automatically taken as the model name. To have a clear distinction, the following naming convention should be used for the parameter-filename and thus the chemistry model for the reactions:

    PlanetAtmosphere_XSpec_XReac_Source (Titan_18Spec_30Reac_Gokcen2007)
    TestCase_XSpec_XReac_Source (CO2_6Spec_10Reac_Johnston2014)
    
In addition to creating a new database, the same script can be used to extend an existing version. For this, only the name of the existing database needs to be defined in the call to create_species_database.py. The function automatically tests if the provided data is already included and adds them if not.

After creating the database the script `extend_reactions.py` in `piclas/tools/archive/` is used to allow one reaction to be part of more than one chemistry model. In this script there are some hard coded changes to ensure a consistent documentation, so no changes by hand had to be made, e.g. changing the reaction model to 'QK' for these reactions 'O2+electron_O2Ion1+electron+electron#1','O+electron_OIon1+electron+electron#1' and 'CO+electron_COIon1+electron+electron#1'.

### Electronic database

This function was also moved to the *tools* folder: `piclas/tools/species_database/`. If a new atom species is added the electronic excitation states are added to the database with the data from the NIST database (https://physics.nist.gov/cgi-bin/ASD/energy1.pl).
The folder `piclas/tools/electronic_database/` was removed to enable a uniform handling of electronic excitation states, but the section is kept for clarity.

A tool to create a database containing electronic excitation states can be found in the *tools* folder: `piclas/tools/electronic_database/`.
The Python script (python3.7) `create_electronic_database_atoms.py` can be used to populate a PICLas-compatible cross-section database, using
the `pandas`, `h5py`, `io`, `re`, `datetime` and `requests` packages. It can be excuted with

    python3.7 `create_electronic_database_atoms.py` 
    
The script gets the data from the NIST database (https://physics.nist.gov/cgi-bin/ASD/energy1.pl) and stores it in an h5-database. Additional species can be added by adapting the `species-list` parameter.

### Reactions

This function was also moved to the *tools* folder: `piclas/tools/species_database/`. If a new reaction is added the necessary data is read in from user inputs or from the name of the reaction.

Reactions can be defined in a separate group as shown above, named by the product species, and added manually using h5py or
[HDF View](https://www.hdfgroup.org/downloads/hdfview/) (Make sure to re-open the file as `Read/Write` to be able to modify and create the dataset.). The ionization process from the LXCat database is not added as a reaction,
however, the cross-section can be used. An additional source for cross-sectional data for reactions is the [NIFS database](https://dbshino.nifs.ac.jp/nifsdb/).