# Tools

This section gives an overview over the tools and scripts contained in the **PICLas** repository.
It also provides references to the tutorials where their usage is explained.
An overview of the tools is given in [TOOLS.md](https://github.com/piclas-framework/piclas/blob/master/tools/TOOLS.md).

(sec:tools-usd)=
## Unified Species Database (USD)

The unified species database was created for a more convenient alternative input for the simulations of species data, electronic states, cross-sections, and chemistry models.
For the general structure and using the unified species database, please refer to Chapter {ref}`sec:unified-species-database` for instructions.


(ssec:tools-maintain-database)=
### Maintain and edit database

A tool to modify or maintain the database it is recommended to use the `maintain_database.py` which can be found in the *tools* directory: `piclas/tools/species_database/`. Simply run

    python3 maintain_database.py

, note that the SpeciesDatabase cannot be used by another program at the same time.

The script shows 5 options to choose from:

    1 to maintain/edit species  or
    2 to maintain/edit chemical reactions  or
    3 to maintain/edit cross section data  or
    4 to maintain/edit surface chemistry  or
    5 to exit program

(sssec:tools-maintain-database-species)=
#### Species parameters


The first option loops over all atom species in the database and compares the electronic excitation states (datasets in the .h5 file) to the electronic excitation states from the [NIST database](https://physics.nist.gov/PhysRefData/ASD/levels_form.html). This option is currently only available for atom species due to the lack of data in other databases for molecules. Afterwards, the attribute values including the mass and heat of formation are compared with data obtained from the [Active Thermochemical Tables](https://atct.anl.gov/) for all species.
Currently the verification of electronic excitation states for molecules is not implemented yet due to the lack of data in other databases.


After selecting to edit/maintain species the script offers the following options

    1 check existing species  or
    2 add new data to existing species  or
    3 add new species  or
    4 to exit program

The first option loops over all existing species and compares the data to different references. Generally the excitation states (electronic excitation states for all species and additionally rotational/vibrational excitation states for molecules) and a few species attributes, specifically the `HeatOfFormation_K`, the atomic mass `MassIC`, the charge `ChargeIC` and the internal species identifier `InteractionID`, are compared.

The electronic excitation states are compared to the data from the [NIST database](https://physics.nist.gov/PhysRefData/ASD/levels_form.html), which only offers this data for atomic species. Due to the lack of data in other databases for the electronic excitation states of molecules, an automated verification of the electronic excitation states is therefore currently not possible.

Alternatively there is the option to add electronic excitation states by providing the data as comma-separated values (csv) in a separate file. For activating this feature the electronic excitation states need to be stored in a .csv file in the `piclas/tools/species_database/` directory. The name and format should follow the template `custom_electronic_levels_SPECIES.csv`, where SPECIES is replaced with the actual name of the species, e.g. 'H2'. The format should look like this

|  J  | Levelcm-1 |
|-----|-----------|
| XXX |    XX     |
| XXX |    XX     |
| 1   |    XX     |

where each row represents one electronic state and the last row contains the ionization energy with a degeneracy of 1. The first column contains the total angular momentum $J$, which is later converted to the degeneracy $g$ with $g = 2J+1$. The second column contains the energy value in 1/cm, which will be converted to Kelvin.

The electronic excitation states are kept, if the data is within a relative tolerance of 1e-05 (which is set in `piclas/tools/species_database/edit_species.py`). If the data is not within this range or the number of electronic levels do not match the differences are displayed and it is possible to choose which data should be kept

    1 to keep data from unified species database  or
    2 to save only electronic level data from https://physics.nist.gov/cgi-bin/ASD/energy1.pl  or
    3 to exit program

Similar to the electronic energy levels, it is possible to store the vibrational energy levels for molecules in the database, which is needed for the Anharmonic Oscillator model ({ref}`sec:DSMC-species`). The vibrational excitation states are provided with as comma-separated values (csv) in a separate file aswell. The template `custom_vibrational_levels_SPECIES.csv` can be used in this case. The format differs here sightly, because each vibrational state is exactly degenerated once, so this column is unnecessary, which leads to

| VibLevelJ |
|-----------|
|    XX     |
|    XX     |
|    XX     |

where the energy is given in Joule. Note that the first level here corresponds to the zero-point energy. The relative tolerance is the same as for the electronic excitation states.

The different attributes are compared for each species, where the `HeatOfFormation_K` and `MassIC` are obtained from the [Active Thermochemical Tables](https://atct.anl.gov/Thermochemical%20Data/version%201.130), while the `ChargeIC` and the `InteractionID` are derived from the species name. If no differences are determined, the output should look like this

    HeatOfFormation_K for SPECIES is equal so will be kept
    MassIC for SPECIES is equal so will be kept
    ChargeIC for SPECIES is equal so will be kept
    InteractionID for SPECIES is equal so will be kept

where SPECIES corresponds to the current species.
Otherwise the differences are displayed and it is possible to select which data should be saved

    1 to keep attributes from unified species database or
    2 to save attributes from https://atct.anl.gov/Thermochemical%20Data/version%201.130 or
    3 to exit program here


The second mode can be used to add new attributes to a species in the database. This can be useful when additional parameters are needed for a new implemented model. The species can be provided by a comma-separated string, e.g. Fe,H2. The excitation states are read in aswell and new data can be added here, if no vibrational or rotational states were found. Afterwards all existing attributes are displayed and the name of the new attribute(s) can be entered. If the attribute already exists, it can be overwritten with the new value. To keep consistancy with newly added parameters, the datatype (e.g. integer, float, str, ...) of the given parameter needs to be added to the specified aswell and is then saved to the config.py file.

The third mode is used to add new species to the database. Which attributes are needed is decided by the name of the species, which needs to be given as input. All attributes which cannot be used from other databases or references, but are needed for the new species (e.g. characteristic temperature(s) for polyatomic species) have to be provided as user input.

The script can be executed by providing the species as arguments

    python3 maintain_database.py --species H Ar

to limit the check to these species.

The other options follow the same logic. If no species is given with the --species argument, the species to add should be input comma seperated string, e.g. `Fe,Ar,H,CIon1,CIon2,C`.


(sssec:tools-maintain-database-reaction)=
#### Chemical reactions

After selecting to edit/maintain reactions the script offers the following options

    1 add new reactions  or
    2 delete reactions  or
    3 to exit program

The first option will expect a list of reactions as a comma separated string, e.g. `C+N_CNIon1+electron,C2+M_C+M+C` and loop over all given reactions. This list can also be provided in the `piclas/tools/species_database/maintain_database.py` file as a list called `code_reactions`. If the reaction already exists in the database all entries will be listed, e.g.

    Reaction: C+N_CNIon1+electron#1
    * Created :  August 02, 2023
    Activation-Energy_K :  164400.0
    Arrhenius-Powerfactor :  1.5
    Arrhenius-Prefactor :  1.66053927673551e-15

    Chemistry Model 1:  Titan_18Spec_30Reac_Gokcen2007

    Products: CNIon1,electron
    Reactants: N,C
    ReactionModel :  TCE

Currently there are the following options to choose from

    1 add a new reaction with different attributes  or
    2 add a new chemistry model to existing reaction  or
    3 to skip this reaction  or
    4 to exit program

When adding a new reaction some parameters such as the reaction model and chemistry model need to be set per user input. Other parameters such as the reactants and products are constructed from the reaction name or if the non reactive species need to be set as well.

The 'delete reactions' option will expect a list of reactions as comma separated string, e.g. `C+N_CNIon1+electron,C2+M_C+M+C`. If there is only one reaction in the database the reaction will be deleted and if there is more than one reaction stored in the database all reactions will be displayed like shown above. To choose which of the listed reactions should be deleted the number(s) of the reaction(s) to delete have to be entered as comma separated string, e.g. `1,2,3`.

(sssec:tools-maintain-database-xsec-collision)=
### Collision cross-sections

The option to create new collision cross-section data is not implemented in the `maintain_database.py` script. To add new collision cross-section data it is recommended to use the old workflow and revert to the regular parameter read-in for these species as described in {ref}`ssec:usd-species` or insert the cross-section data by hand to the unified species database.

A tool to create a database containing cross-section data can be found in the *tools* folder: `piclas/tools/crosssection_database/`.
The Python script (python3.7) `create_xsec_db_lxcat.py` can be used to populate a PICLas-compatible cross-section database, using the `numpy`, `h5py` and `lxcat_data_parser` packages.

    python3.7 create_xsec_db_lxcat.py

A database (containing multiple species and cross-section types) downloaded directly from the Plasma Data Exchange Project and the
[LXCat database](https://fr.lxcat.net/home/) and the name of output database can be supplied to the script with

    database_input = "Database.txt"
    database_output = "Database.h5"

Currently, PICLas only utilizes the elastic, effective and vibrational cross-sections, however, all excitation cross-section types
are grouped and stored in the output file. An example is given below

    CO2-electron (group)
        EFFECTIVE (dataset)
        ROTATION (group)
            0.02 (dataset)
        VIBRATION (group)
            0.29
            0.59
        REACTION (group)
            CO2Ion1-electron-electron

Datasets, which cannot be identified as rotational, vibrational or electronic excitation will grouped within an `UNDEFINED` group.
By defining a species list, only certain species can be included in the output database

    species_list = ["Ar","CO"]

Finally, the utilized cross-section data should be properly referenced by adding the information to the HDF5 database as an attribute

    reference = 'XXX database, www.lxcat.net, retrieved on MMMM DD, YYYY.'

Users of cross-section data are encouraged to download the data directly from the [LXCat project website](https://fr.lxcat.net/home/)
and to consider the guidelines regarding referencing and publication.

(sssec:tools-maintain-database-surfchem)=
### Surface chemistry

The option to create new surface chemistry data is not implemented in the `maintain_database.py` script. To add surface chemistry data it is recommended to add data by hand to the unified species database.

## Userblock

The `userblock` contains the complete information about a **PICLas** run (git branch of the
repository, differences to that branch, `cmake` configuration and parameter file) and is prepended
to every `.h5` state file. The parameter file is prepended in ASCII format, the rest is binary and
is generated automatically during the build process with the `generateuserblock.sh` script.

### `extract_userblock.py`

It can be extracted and printed using the `extract_userblock.py` script. Its basic usage is

    python2 extract_userblock.py -XXX [statefile.h5]

where `-XXX` can be replaced by

* `-s` to show all available parts of the userblock (such as `CMAKE` or `GIT BRANCH`)
* `-a` to print the complete userblock
* `-p [part]` to print one of the parts listed with the `-s` command.

### `rebuild.py`

The second python tool in this folder is `rebuild.py`. It extracts the userblock from a state file
and builds a **PICLas** repository and binary identical to the one that the state file was created
with. In order to do so, it clones a **PICLas** git repository, checks out the given branch, applies
the stored changes to the git `HEAD` and builds **PICLas** with the stored `cmake` options.
If run with the parameter file given in the `INIFILE` part of the userblock, this binary should
reproduce the same results/behaviour (possible remaining sources of different output are for example
differences in restart files, compilers, linked libraries or machines). The basic usage is

    python2 rebuild.py [dir] [statefile.h5]

where `dir` is an empty directory that the repository is cloned into and where the `piclas`
executable is built. `statefile.h5` is the state file whose userblock is used to rebuild the `piclas`
executable. Help can be shown via `-h` for both userblock scripts.


