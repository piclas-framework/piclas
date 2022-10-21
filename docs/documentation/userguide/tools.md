# Tools

This section gives an overview over the tools and scripts contained in the **PICLas** repository.
It also provides references to the tutorials where their usage is explained.
An overview of the tools is given in [TOOLS.md](https://github.com/piclas-framework/piclas/blob/master/tools/TOOLS.md).

(sec:tools-usd)=
## Unified Species Database (USD)

A tool to create a database containing cross-section, electronic states, surface reactions, Arrhenius rates, and species data can be found in the *tools* folder: `piclas/tools/species_database/`.
The Python script (python3.7) `create_species_database.py` creates a new database or expands an existing one combining all necessary parameters (formerly read-in through ini files). The script uses the `numpy`, `h5py`, `argparse`,`datetime`, `cmath`, and `matplotlib.rcsetup` packages. To create the unified database with

    python3.7 create_species_database.py
    
an electronic states and a cross-section database need to be built before ({ref}`ssec:tools-xsec-collision`, {ref}`ssec:tools-electronic-database`).
If nothing additionally is specified, the following filenames are called: `DSMC.ini` for the parameter and (gas and surface) reaction input, `Electronic-State-Database.h5` for the electronic state database, `XSec_Database.h5` for the cross-section data,`Rad.dat` for the radiation data, and `Species_Database.h5` for the final output. For custom file names, the following options can be added: 

    python3 create_species_database.py --parameter parameter-filename --electronic electronic_statefile --crosssection crosssection_statefile --radiation radiation-filename --output output_filename --reference reference-name

or

    python3 create_species_database.py -p parameter-filename -e electronic_statefile -c crosssection_statefile -s radiation-filename -o output_filename -r reference-name
    
The data is grouped in the output file, as shown in the following example:
    
    Cross-Sections (group)
        H2-H2Ion1 (dataset)
    Radiation (group)
        N (group)
            Levels (dataset)
            Lines (dataset)
    Reaction (group)
        DISS_CH3_TO_CH2 (dataset)
        Arrhenius parameters (attributes)
        EXC_CN+_N (dataset)
        Arrhenius parameters (attributes)
        RECOMB_C_N_ION (dataset)
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
    Surface-Chemistry
        Adsorption_CO (dataset)
        Reaction parameters (attributes)
        LH_Oxidation_CO (dataset)
        Reaction parameters (attributes)
        
For cross-sections, reactions, radiation and species data, the former `DSMC.ini` files are used to create the database. However, for every species and for every reaction, names must be defined to create the database and to run simulations.

    Part-Species1-SpeciesName=CO2
    DSMC-Reaction1-ReactionName = DISS_CH3_TO_CH2
    Surface-Reaction1-SurfName = Adsorption_CO
    
        

(ssec:tools-electronic-database)=
### Electronic database

A tool to create a database containing electronic excitation states can be found in the *tools* folder: `piclas/tools/electronic_database/`.
The Python script (python3.7) `create_electronic_database_atoms.py` can be used to populate a PICLas-compatible cross-section database, using
the `pandas`, `h5py`, `io`, `re`, `datetime` and `requests` packages. It can be excuted with

    python3.7 `create_electronic_database_atoms.py` 
    
The script gets the data from the NIST database (https://physics.nist.gov/cgi-bin/ASD/energy1.pl) and stores it in an h5-database. Additional species can be added by adapting the `species-list` parameter.

(ssec:tools-xsec-collision)=
### Collision cross-section database

A tool to create a database containing cross-section data can be found in the *tools* folder: `piclas/tools/crosssection_database/`.
The Python script (python3.7) `create_xsec_db_lxcat.py` can be used to populate a PICLas-compatible cross-section database, using
the `numpy`, `h5py` and `lxcat_data_parser` packages.

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

Chemical reaction can be added to the database manually using [HDF View](https://www.hdfgroup.org/downloads/hdfview/).
Make sure to re-open the file as `Read/Write` to be able to modify and create the dataset.

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


