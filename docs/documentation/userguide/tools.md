# Tools

This section gives an overview over the tools and scripts contained in the **PICLas** repository.
It also provides references to the tutorials where their usage is explained.
An overview of the tools is given in [TOOLS.md](https://github.com/piclas-framework/piclas/blob/master/tools/TOOLS.md).

(sec:tools-xsec-collision)=
## Collision cross-section database

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


