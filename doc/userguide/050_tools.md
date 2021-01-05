\hypertarget{tools}{}

# Tools Overview \label{chap:tools}

This section gives an overview over the tools and scripts contained in the **PICLas** repository. It also provides references to the tutorials where their usage is explained.

## Collision cross-section database \label{sec:tools_mcc}

A tool to create a database containing cross-section data can be found in the *tools* folder: `piclas/tools/crosssection_database/`. The Python script (python3.7) `create_xsec_db_lxcat.py` can be used to populate a PICLas-compatible cross-section database, using the `numpy`, `h5py` and `lxcat_data_parser` packages.

    python3.7 create_xsec_db_lxcat.py

A database (containing multiple species and cross-section types) downloaded directly from the Plasma Data Exchange Project and the [LXCat database](https://fr.lxcat.net/home/) and the name of output database can be supplied to the script with

    database_input = "Database.txt"
    database_output = "Database.h5"

Currently, PICLas only utilizes the elastic, effective and vibrational cross-sections, however, all excitation cross-section types are grouped and stored in the output file. An example is given below 

    CO2-electron (group)
        EFFECTIVE (dataset)
        ROTATION (group)
            0.02 (dataset)
        VIBRATION (group)
            0.29
            0.59
        REACTION (group)
            CO2Ion1-electron-electron

Datasets, which cannot be identified as rotational, vibrational or electronic excitation will grouped within an `UNDEFINED` group. By defining a species list, only certain species can be included in the output database

    species_list = ["Ar","CO"]

Finally, the utilized cross-section data should be properly referenced by adding the information to the HDF5 database as an attribute

    reference = 'XXX database, www.lxcat.net, retrieved on MMMM DD, YYYY.'

Users of cross-section data are encouraged to download the data directly from the [LXCat project website](https://fr.lxcat.net/home/) and to consider the guidelines regarding referencing and publication.

Chemical reaction can be added to the database manually using [HDF View](https://www.hdfgroup.org/downloads/hdfview/). Make sure to re-open the file as `Read/Write` to be able to modify and create the dataset.

## Visualization (NEEDS AN UPDATE)

### Installation of ParaView and Plugin

A ParaView reader based on `posti_visu` to load **PICLas** state files in ParaView. Provides the interface to adjust `posti_visu` parameters in the ParaView GUI. For this purpose the `libVisuReader.so` has to be loaded as a Plugin in ParaView.

#### Currently Tested Combinations

* GNU **v7.3.0** with Ubuntu **v18.04**
* Paraview **v5.6.0** (via cmake and GNU)
* OPENMPI **v3.1.2** (via configure and GNU)
* HDF5 **v1.10.2** parallel compiled (via cmake and GNU)
* **flexi** and **fleximultiphase** with GNU

and

* GNU **v7.3.0** with Mint **v19**
* Paraview **v5.6.0** (via cmake and GNU)
* OPENMPI **v3.1.3** (via configure and GNU)
* HDF5 **v1.10.0-patch1** parallel compiled (via cmake and GNU)
* **flexi**  with GNU

and

* INTEL **v19.0** with Ubuntu **v18.04**
* Paraview **v5.6.0** (via cmake and Intel)
* OPENMPI **v3.1.2** (via configure and Intel)
* HDF5 **v1.10.2** parallel compiled (via cmake and Intel)
* **flexi** and **fleximultiphase** with Intel

What needs to be considered when comparing with older versions, eg. **v5.3.0**?
* Paraview should be compiled with **an own HDF5** (pre-installed in /opt/hdf5/...)
* Paraview uses **by default qt5** (see Paraview with Qt5)

What else should be considered?

Use the following "cmake .." command:
```
ccmake .. -DHDF5_PARALLEL=ON
```
or
```
cmake .. -DHDF5_PARALLEL=ON ... "Flags siehe unten"
```
Set the following variables:

**new**
* VTK_USE_SYSTEM_HDF5 = **ON**
* HDF5_IS_PARALLEL = **ON**
* VTK_MODULE_vtkhdf5_IS_SHARED = **OFF**

**old**
- CMAKE_BUILD_TYPE = **Release**
- PARAVIEW_ENABLE_PYTHON = **ON**
- PARAVIEW_USE_MPI = **ON**
- PARAVIEW_INSTALL_DEVELOPMENT_FILES = **ON**
- CMAKE_INSTALL_PREFIX = **/opt/paraview/...**

In order to set all packages correctly, e.g., HDF5 (it does not matter if compiled with cmake or configure), it is assumed that all
environment variables are set correctly (e.g. .bashrc):

* **PATH**
* **LD_LIBRARY_PATH**
* **CMAKE_PREFIX_PATH**
* **CMAKE_LIRBRARY_PATH**
* **CMAKE_INCLUDE_PATH**


Afterwards, everything should work as expected.

What else must be considered when using the **INTEL compiler**?

*  PARAVIEW_ENABLE_MOTIONFX = **OFF**
*  [bugfix for plugin](https://gitlab.com./fleximultiphase/Codes/fleximultiphase/merge_requests/141)

#### Overview
For Ubuntu 16.04:
```
sudo apt-get install gfortran libpython-dev libboost-dev make libphonon-dev libphonon4 qt4-dev-tools libqt4-dev qt4-qmake libxt-dev g++ gcc cmake-curses-gui libqt4-opengl-dev mesa-common-dev git vim
```

Additional steps:
* OpenMPI (Version: 2.0.0 tested. Others should also work.) compiled with GNU (set exports in the .bashrc and reload!)
  **ATTENTION: All other Software (hdf5, ParaView, Flexi, ...) must be compiled with the same MPI version.**
* HDF5 (Version: 1.10.0-patch1 tested. Others should also work.) with cmake as shown in the following,
  but remove the last hdf5 in the export (place in the .bashrc), as follows:
   * export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake/
* Paraview (5.4.1) downloaded and extracted
* Vote for this fix and apply the patch: [FileParser-Patch fuer Timestamps](https://paraview.uservoice.com/forums/11350-general/suggestions/12591231-expand-the-name-patterns-for-temporal-file-series)
```
patch -p1 < FileParser.patch
```
im ParaView root directory.
   * PATCH for 5.3.0: [FileParser.patch](/uploads/ec2723307910f061389d83bba0c1897c/FileParser.patch), [ErrorMessage.patch](/uploads/793df4cdaee862d10a48fdc321b8df93/ErrorMessage.patch)
   * PATCH for 5.4.1: [FileParser.patch](/uploads/4f413a1cac1a04645749def3b016c28d/FileParser.patch), [ErrorMessage.patch](/uploads/412d199f264c51aafebb871e16fd730b/ErrorMessage.patch)
* Compile ParaView with the following cmake command and don't forget the exports
* Compile the plugin

#### OpenMPI installation
OpenMPI compilation (Versions 1.8.4, 1.8.8, 2.0.0, 3.1.3 tested):
```
 GNU: ./configure --prefix=/opt/openmpi/1.X.X/gnu
 INTEL: ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=/opt/openmpi/1.X.X/intel
 make && make install
```
Set the environment variables (in /etc/profile or .bashrc or .zshrc ...):
```
 export PATH=/opt/openmpi/1.X.X/YYYYY/bin:$PATH
 export LD_LIBRARY_PATH=/opt/openmpi/1.X.X/YYYYY/lib/:$LD_LIBRARY_PATH
 export CMAKE_PREFIX_PATH=/opt/openmpi/1.X.X/YYYYY/:$CMAKE_PREFIX_PATH
```
**ATTENTION: reload .bashrc or. .zshrc!!!**

For older Ubuntu versions (14.04), the OpenMPI version from the package manager can be used (version 1.6.x) but this is **NOT
ENCOURAGED**. For Ubuntu 16.04, OpenMPI 1.10.x is supplied, which should be compiled, otherwise the plugin might not work later on.
#### HDF5 installation
HDF5 compilation (Versions 1.8.13, 1.8.14, 1.8.16, 1.10.0-patch1 tested):

###### either via '''cmake''' (**ENCOURAGED**)
```
 cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_ENABLE_F2003=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/hdf5/1.X.X ..
 make && make install
```
Set the environment variables (in /etc/profile or .bashrc or .zshrc ...):
```
 export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake/hdf5
 export CMAKE_PREFIX_PATH=/opt/hdf5/1.X.X/:$CMAKE_PREFIX_PATH
```
Attention: At least for HDF5 1.8.16 the correct path does not include /hdf5:
```
 export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake
```
**ATTENTION: reload .bashrc or .zshrc!!!**

###### or via configure (**NOT ENCOURAGED**)
```
 ./configure --prefix=/opt/hdf5/1.X.X --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel
 make && make install
```
Set the environment variables (in /etc/profile or .bashrc or .zshrc ...):
```
 export HDF5_DIR=/opt/hdf5/1.X.X/
 export CMAKE_PREFIX_PATH=/opt/hdf5/1.X.X/:$CMAKE_PREFIX_PATH
```
**ATTENTION: reload .bashrc or .zshrc!!!**


#### ParaView installation
Download the ParaView source files, e.g., from github
```
git clone https://github.com/Kitware/ParaView.git ./paraview
git submodule update --init --recursive
```

ParaView compilation (Versions 4.3.1, 5.0.0, 5.3.0, 5.4.1 tested):

See the "Prerequisites" listed on the ParaView documentation [ParaView-Doku](http://www.paraview.org/Wiki/ParaView:Build_And_Install#Download_ParaView_Source_Code).
For Ubuntu 14.04 and 16.04 the following packages are required:
```
sudo apt-get install libqt4-dev libboost-dev libpython-dev qt4-dev-tools
```
If the following executables are not found during the installation in the qt4/bin folder, a symbolic link to /usr/bin/ can be set
for:
* usr/bin/xmlpatterns
* (usr/bin/qhelpgenerator)
* if somebody required VisitBridge (e.g. for CGNS files), set -DPARAVIEW_USE_VISITBRIDGE to ON
* Apply the patch here (please vote): [FileParser-Patch fuer Timestamps](https://paraview.uservoice.com/forums/11350-general/suggestions/12591231-expand-the-name-patterns-for-temporal-file-series)
   * PATCH for 5.3.0: [FileParser.patch](/uploads/ec2723307910f061389d83bba0c1897c/FileParser.patch), [ErrorMessage.patch](/uploads/793df4cdaee862d10a48fdc321b8df93/ErrorMessage.patch)
   * PATCH for 5.4.1: [FileParser.patch](/uploads/4f413a1cac1a04645749def3b016c28d/FileParser.patch), [ErrorMessage.patch](/uploads/412d199f264c51aafebb871e16fd730b/ErrorMessage.patch)

ParaView compilation:
```
 cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_VISITBRIDGE=OFF -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_QT_VERSION=4 -DCMAKE_INSTALL_PREFIX=/opt/paraview/X.X.X ..
 make && make install
```
Set the environment variables (in /etc/profile or .bashrc or .zshrc ...):
```
 export ParaView_DIR=/opt/paraview/X.X.X
```
**ATTENTION: reload .bashrc or .zshrc!!!**

##### ParaView with Qt5

The new ParaView versions are built with Qt5 by default. The following packages are required (tested with Ubuntu 18.04):
```
sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libxt-dev libgl1-mesa-dev
```

#### Compile the ParaView-Plugin
* Set the Flexi cmake to compile the plugin
```
FLEXI_BUILDPOSTI = ON
POSTI_BUILD_VISU = ON
POSTI_USE_PARAVIEW = ON
```
* set an appropriate alias for ParaView in .bash_aliases/.zshrc/...
```
 alias paraview='paraview --mpi'
```
* Multiple plugins can be used, e.g., from different `build` folders, by setting the variable `PV_PLUGIN_PATH`. The following command can be utilized
```
PV_PLUGIN_PATH=pfad_zum_flexi_build_ordner/lib paraview
```

### ParaView on Hazel Hen (ParaView 5.3)

HLRS is providing ParaView as a module. Some special steps have to be taken to connect a ParaView client on a local machine to a server on the Hazel Hen.

#### Compile client

First, you need to compile the sources provided by the HLRS, since they are not the same as the official 5.3.0 sources (because reasons). The source files can be found on Hazel Hen in the folder `/opt/hlrs/tools/paraview/5.3`. Follow the steps above to compile these source files on your local machine, but make sure that the following options are set:
* PARAVIEW_QT_VERSION=4
* VTK_RENDERING_BACKEND=OpenGl

since otherwise the handshake between client and server will not work.

#### Compile Plugin on Hazel Hen

Compiling the Plugin on Hazel Hen must be done on a mom node, not a login node since during the compile process some helper programs are called which can not be executed on a login node. So first, get an interactive job.
The Plugin requires dynamic linking, activate that by running
```
export CRAYPE_LINK_TYPE=dynamic
```
then load the required modules
```
module switch gcc/8.2.0 gcc/7.3.0
module load tools/paraview/5.3.0-git-master-parallel-Mesa
module unload cray-hdf5-parallel/1.10.2.0
module load cray-hdf5-parallel/1.10.1.1
```
ATTENTION: If you don't do it in this order, the wrong HDF5 module will be loaded and/or the environment variables not set correctly. To be sure everything worked, check if `HDF5_DIR` is pointing to `/opt/cray/pe/hdf5-parallel/1.10.1.1/GNU/5.1/`.

The older compiler version is necessary since HDF5 1.10.1.1 (which is the one used by ParaView) was build with that version.

Switch to/create your build directory and generate a Makefile - remember to turn `USE_PARAVIEW` on. Compilation must be done using the `aprun` command:
```
aprun -n 1 make visuReader
```
The Plugin should now be available in the `lib/` directory.

#### Start server

Next the server needs to be started on Hazel Hen. When the server should be started in an interactive job, run a command like
```
qsub -I -l nodes=1:ppn=24,walltime=01:00:00
```
to get an interactive job. When the job starts, change again to the older compiler and load the ParaView module
```
module switch gcc/8.2.0 gcc/7.3.0
module load tools/paraview/5.3.0-git-master-parallel-Mesa
```
and then start the ParaView server using aprun
```
aprun -n 24 pvserver
```

#### Connect to server

STEP 1:

The mom node has to be reachable from your local machine. To achieve this, there are several options:
* Option 1: Add the mom nodes to your ssh config using the login nodes as a [proxy ](https://wiki.gentoo.org/wiki/SSH_jump_host). Add something like this to your `~/.ssh/config` file:
```
Host hen
  User iagmustermann
  HostName hazelhen.hww.de

Host mom01
  User iagmustermann
  HostName mom01
  ProxyJump hen
```
Repeat the last four lines for all mom nodes. The `User` lines are only needed if your username on HazelHen differs from that on your local machine. For ssh versions less than 7.3, replace the line `ProxyJump  hen` by `ProxyCommand ssh -W %h:%p hen`.
* Option 2: Use the special queue `io` which will put the job on a mom node that is reachable from everywhere inside of the university network. (ATTENTION: The `io`-queue is currently only open for test users! The nodes are named something like `hazelhen-network10` and can be reached using `hazelhen-network10.hww.hlrs.de`.)
* Option 3: Change the script to manually tunnel twice - once from your local machine to a login node and then to the mom node.

STEP 2:

To connect to the server, a helper script `pvconnect` can be found on Hazel Hen in the folder `/sw/hazelhen/hlrs/tools/paraview/5.3.0-git-master/Pre/bin`. A slightly modified version can be found [here](uploads/15eb2c6e5c10f62b166fdd38ffab39ed/pvconnect) (adapt `/PATH/TO/PARAVIEW` in the script before running it).  This script will create a ssh tunnel to the mom node that runs the pvserver job and then starts a ParaView-client on your local machine that automatically connects to the server. The syntax is
```
pvconnect -pvs nidXXXXXX:YYYYYYY -via momZZ
```
where nidXXXXXX is the name of the pvserver and YYYYYY the port of the server, both of which are shown when the server starts. The connection is created through momZZ, where you need to insert the number of the mom node the pvserver job runs on.

If you work with the `PV_PLUGIN_PATH` environment variable: Reset it with `export PV_PLUGIN_PATH=''` before starting the script to avoid errors.

<!-- ## Swap meshes

---------------------------------------------------------------------------------------------
**posti_swapmesh**
--------------------------------- ------------------------------------------------------------
Brief description                  Interpolates state file data from one mesh to another. Uses high-order interpolation and a Newton coordinate search algorithm. Meshes do not have to be conforming. A reference state can be given for areas in the target mesh not covered by the source mesh.

Basic usage                        `posti_swapmesh [parameter.ini] [statefile.h5]`

Further info / usage example       No tutorials so far
--------------------------------------------------------------------------------------------- -->

<!-- ## Record points

----------------------------------------------------------------------------------------------
**posti_preparerecordpoints**
--------------------------------- ------------------------------------------------------------
Brief description                  Enables **PICLas** to record values at a set of physical points over time with a higher temporal sampling rate than the state file output interval. The record point coordinates and the **PICLas** mesh are defined in the parameter file. Creates an additional `.h5` file, whose path is passed to **PICLas** as a parameter.

Basic usage                        `posti_preparerecordpoints [parameter_prepareRP.ini]`

Further info / usage example       \ref{sec:postiRecordpoints}
----------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------
**posti_visualizerecordpoints**
--------------------------------- ------------------------------------------------------------
Brief description                  Performs post-processing of the `*_RP_*` files written by **PICLas**: merges several time steps and writes output such as value over time or spectra.

Basic usage                        `posti_visualizerecordpoints [parameter_visuRP.ini] [projectname_RP_*.h5]`

Further info / usage example       \ref{sec:postiRecordpoints}
----------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------
**posti_evaluaterecordpoints**
--------------------------------- ------------------------------------------------------------
Brief description                  Evaluate the values at recorpoints a posteri from existing statefiles. Can be used if the recordpoints have not been set during the simulation, but will only give coarse temporal resolution.

Basic usage                        `posti_evaluaterecordpoints [parameter.ini] [statefile.h5]`

Further info / usage example       No tutorials so far
----------------------------------------------------------------------------------------------- -->

## Initialization of An Initial Particle Distribution Function (NEEDS UPDATING)

Currently, the initial particle distributions are simple functions which do not consider great changes in any
coordinate. Here is a fast HOW-TO of creating any initial particle distribution.

Following information are given analytically: n=n(x),v=v(x). The initial distribution which are non-constant
are evaluated in each element using the DSMC cell_local sampling. The required values are evaluated at the
element-origin. This requires a fine sampling mesh with a sufficient high resolution to sample the PDF. These information
are then stored in the InitPart_DSMCHOState.h5 and used for a restart with the MacroRestartValues, resulting in the
InitPDF_State.h5. From this file, the PartData is copied to the initial file, generated with the coarse computational grid


    h5copy -i InitPDF_State.h5 -s PartData -o Coarse_State.h5 -d PartData


Finally, the last entry of PartInt of Coarse_State.h5 has to be set to the number of all particles in InitPDF_State.h5. Note,
that PICLas is executed on a single core.

## Tools Folder

PICLas comes with a collection of loose tools that perform various tasks, e.g., post-processing. An
overview of the tools is given in [TOOLS.md](https://github.com/piclas-framework/piclas/blob/master/tools/TOOLS.md).

### Userblock

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


