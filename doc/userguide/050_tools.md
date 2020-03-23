\hypertarget{tools}{}

# Tools Overview \label{chap:tools}

This section gives an overview over the tools and scripts contained in the **PICLas** repository. It also provides references to the tutorials where their usage is explained.

## Collision cross-section database \label{sec:tools_mcc}

A tool to create a database containing cross-section data can be found in the *tools* folder: `piclas/tools/crosssection_database/`. The Python script (python3.7) `create_xsec_db_lxcat.py` can be used to populate a PICLas-compatible cross-section database, using the `numpy`, `h5py` and `lxcat_data_parser` packages.

    python3.7 create_xsec_db_lxcat.py

A database (containing multiple species and cross-section types) downloaded directly from the Plasma Data Exchange Project and the [LXCat database](https://fr.lxcat.net/home/) and the name of output database can be supplied to the script with

    database_input = "Database.txt"
    database_output = "Database.h5"

Currently, PICLas only utilizes effective cross-sections between neutral species and electrons and as such only these cross-section types are stored in the output file. By defining a species list, only certain species can be included in the output database

    species_list = ["Ar","CO"]

Finally, the utilized cross-section data should be properly referenced by adding the information to the HDF5 database as an attribute

    reference = 'XXX database, www.lxcat.net, retrieved on MMMM DD, YYYY.'

Users of cross-section data are encouraged to download the data directly from the [LXCat project website](https://fr.lxcat.net/home/) and to consider the guidelines regarding referencing and publication.

## Visualization (NEEDS AN UPDATE)

### Installation of ParaView and Plugin

A ParaView reader based on `posti_visu` to load **PICLas** state files in ParaView. Provides the interface to adjust `posti_visu` parameters in the ParaView GUI. For this purpose the `libVisuReader.so` has to be loaded as a Plugin in ParaView.

#### Aktuell

Getestet mit:

* GNU **v7.3.0** mit Ubuntu **v18.04**
* Paraview **v5.6.0** (via cmake und GNU)
* OPENMPI **v3.1.2** (via configure und GNU)
* HDF5 **v1.10.2** parallel gebaut (via cmake und GNU)
* **flexi** und **fleximultiphase** mit GNU

und

* GNU **v7.3.0** mit Mint **v19**
* Paraview **v5.6.0** (via cmake und GNU)
* OPENMPI **v3.1.3** (via configure und GNU)
* HDF5 **v1.10.0-patch1** parallel gebaut (via cmake und GNU)
* **flexi**  mit GNU

und

* INTEL **v19.0** mit Ubuntu **v18.04**
* Paraview **v5.6.0** (via cmake und Intel)
* OPENMPI **v3.1.2** (via configure und Intel)
* HDF5 **v1.10.2** parallel gebaut (via cmake und Intel)
* **flexi** und **fleximultiphase** mit Intel

Was ist zu beachten im Vergleich zu den älteren Versionen, z.B. Paraview **v5.3.0**?
* Paraview sollte mit **eigenem HDF5** (vorinstalliert in /opt/hdf5/...) gebaut werden
* Paraview verwendet **standardmäßig qt5** (siehe Paraview mit Qt5)

Was ist sonst noch zu beachten?

Führe anstatt "ccmake .." folgenden Befehl aus:
```
ccmake .. -DHDF5_PARALLEL=ON
```
oder
```
cmake .. -DHDF5_PARALLEL=ON ... "Flags siehe unten"
```

Setze während des Konfigurierens:

**neu**
* VTK_USE_SYSTEM_HDF5 = **ON**
* HDF5_IS_PARALLEL = **ON**
* VTK_MODULE_vtkhdf5_IS_SHARED = **OFF**

**alt**
- CMAKE_BUILD_TYPE = **Release**
- PARAVIEW_ENABLE_PYTHON = **ON**
- PARAVIEW_USE_MPI = **ON**
- PARAVIEW_INSTALL_DEVELOPMENT_FILES = **ON**
- CMAKE_INSTALL_PREFIX = **/opt/paraview/...**

Damit alle Pakete richtig gefunden werden, z.B. HDF5 (egal ob mit cmake oder configure), wird natürlich vorausgesetzt, dass alle Umgebungsvariablen

* **PATH**
* **LD_LIBRARY_PATH**
* **CMAKE_PREFIX_PATH**
* **CMAKE_LIRBRARY_PATH**
* **CMAKE_INCLUDE_PATH**

mit den relevanten Pfaden in der ".bashrc" richtig gesetzt sind.

Anschließend sollte alles wie gehabt funktionieren.

Was ist zusätzlich mit dem **INTEL compiler** zu beachten?

*  PARAVIEW_ENABLE_MOTIONFX = **OFF**
*  [bugfix für plugin](https://gitlabext.iag.uni-stuttgart.de/fleximultiphase/Codes/fleximultiphase/merge_requests/141)

#### Übersicht
Fuer ein jungfräuliches Ubuntu 16.04 folgende Pakete installieren:
```
sudo apt-get install gfortran libpython-dev libboost-dev make libphonon-dev libphonon4 qt4-dev-tools libqt4-dev qt4-qmake libxt-dev g++ gcc cmake-curses-gui libqt4-opengl-dev mesa-common-dev git vim
```

Weitere Schritte wie folgt:
* OpenMPI (Version: 2.0.0 getestet. Andere sollten aber auch funktionieren.) mit GNU bauen (exports in der .bashrc nicht vergessen und diese auch neu laden!!!) **ATTENTION: Alle andere Software (hdf5, ParaView, Flexi, ...) sollte mit diesem MPI gebaut werden.**
* HDF5 (Version: 1.10.0-patch1 getestet. Andere sollten aber auch funktionieren.) mit cmake wie unten beschrieben bauen, aber beim export das letzte hdf5 streichen (am besten wieder in die .bashrc), also so:
   * export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake/
* Paraview (5.4.1) runterladen und entpacken
* Hier: [FileParser-Patch fuer Timestamps](https://paraview.uservoice.com/forums/11350-general/suggestions/12591231-expand-the-name-patterns-for-temporal-file-series) voten und den patch anwenden mit
```
patch -p1 < FileParser.patch
```
im paraview root directory.
   * PATCH fuer 5.3.0: [FileParser.patch](/uploads/ec2723307910f061389d83bba0c1897c/FileParser.patch), [ErrorMessage.patch](/uploads/793df4cdaee862d10a48fdc321b8df93/ErrorMessage.patch)
   * PATCH fuer 5.4.1: [FileParser.patch](/uploads/4f413a1cac1a04645749def3b016c28d/FileParser.patch), [ErrorMessage.patch](/uploads/412d199f264c51aafebb871e16fd730b/ErrorMessage.patch)
* Paraview mit dem untigem cmake-Commando bauen, wieder die exports nicht vergessen
* Plugin bauen

#### OpenMPI kompilieren
OpenMPI selbst kompilieren (Versionen 1.8.4, 1.8.8, 2.0.0, 3.1.3 getestet):
```
 GNU: ./configure --prefix=/opt/openmpi/1.X.X/gnu
 INTEL: ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=/opt/openmpi/1.X.X/intel
 make && make install
```
Setzen der Umgebungsvariablen (in /etc/profile oder .bashrc oder .zshrc ...):
```
 export PATH=/opt/openmpi/1.X.X/YYYYY/bin:$PATH
 export LD_LIBRARY_PATH=/opt/openmpi/1.X.X/YYYYY/lib/:$LD_LIBRARY_PATH
 export CMAKE_PREFIX_PATH=/opt/openmpi/1.X.X/YYYYY/:$CMAKE_PREFIX_PATH
```
**ATTENTION: .bashrc bzw. .zshrc neu laden!!!**

Für ältere Ubuntu (14.04) kann hier auch direkt das OpenMPI aus den Paketquellen (Version 1.6.x) verwendet werden (**NICHT EMPFOHLEN**). Bei Ubuntu 16.04 wird OpenMPI 1.10.x mitgeliefert, hier **selber kompilieren**! Ansonsten geht die Ausführung des Plugins später schief.

#### HDF5 kompilieren
HDF5 bauen (Versionen 1.8.13, 1.8.14, 1.8.16, 1.10.0-patch1 getestet):

###### entweder mit '''cmake''' (**EMPFOHLEN**)
```
 cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_ENABLE_F2003=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/hdf5/1.X.X ..
 make && make install
```
Setzen der Umgebungsvariablen (in /etc/profile oder .bashrc oder .zshrc ...):
```
 export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake/hdf5
 export CMAKE_PREFIX_PATH=/opt/hdf5/1.X.X/:$CMAKE_PREFIX_PATH
```
Achtung: Zumindest bei HDF5 1.8.16 lautet der richtige Pfad (ohne /hdf5):
```
 export HDF5_DIR=/opt/hdf5/1.X.X/share/cmake
```
**ATTENTION: .bashrc bzw. .zshrc neu laden!!!**

###### oder mit configure (**NICHT EMPFOHLEN**)
```
 ./configure --prefix=/opt/hdf5/1.X.X --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel
 make && make install
```
Setzen der Umgebungsvariablen (in /etc/profile oder .bashrc oder .zshrc ...):
```
 export HDF5_DIR=/opt/hdf5/1.X.X/
 export CMAKE_PREFIX_PATH=/opt/hdf5/1.X.X/:$CMAKE_PREFIX_PATH
```
**ATTENTION: .bashrc bzw. .zshrc neu laden!!!**


#### ParaView kompilieren
Paraview source files gibt es auf github unter
```
git clone https://github.com/Kitware/ParaView.git ./paraview
git submodule update --init --recursive
```

ParaView bauen (Versionen 4.3.1, 5.0.0, 5.3.0, 5.4.1 getestet):

Siehe auch "Prerequisits" auf der Paraview Dokumentationsseite [ParaView-Doku](http://www.paraview.org/Wiki/ParaView:Build_And_Install#Download_ParaView_Source_Code).
Für Ubuntu 14.04 und 16.04 werden u.a. folgende Pakete benötigt:
```
sudo apt-get install libqt4-dev libboost-dev libpython-dev qt4-dev-tools
```
Falls die folgenden Executables beim Kompilieren im qt4/bin Ordner nicht gefunden werden, kann dort ein symbolischer Link zu /usr/bin/ gesetzt werden:
* usr/bin/xmlpatterns
* (usr/bin/qhelpgenerator)
* wer die VisitBridge benoetigt (z.B. fuer CGNS Files), der stellt -DPARAVIEW_USE_VISITBRIDGE auf ON
* Hier: [FileParser-Patch fuer Timestamps](https://paraview.uservoice.com/forums/11350-general/suggestions/12591231-expand-the-name-patterns-for-temporal-file-series) voten und den patch anwenden
   * PATCH fuer 5.3.0: [FileParser.patch](/uploads/ec2723307910f061389d83bba0c1897c/FileParser.patch), [ErrorMessage.patch](/uploads/793df4cdaee862d10a48fdc321b8df93/ErrorMessage.patch)
   * PATCH fuer 5.4.1: [FileParser.patch](/uploads/4f413a1cac1a04645749def3b016c28d/FileParser.patch), [ErrorMessage.patch](/uploads/412d199f264c51aafebb871e16fd730b/ErrorMessage.patch)

ParaView bauen:
```
 cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_VISITBRIDGE=OFF -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_QT_VERSION=4 -DCMAKE_INSTALL_PREFIX=/opt/paraview/X.X.X ..
 make && make install
```

Setzen der Umgebungsvariablen (in /etc/profile oder .bashrc oder .zshrc ...):
```
 export ParaView_DIR=/opt/paraview/X.X.X
```
**ATTENTION: .bashrc bzw. .zshrc neu laden!!!**

##### ParaView mit Qt5

Die neueren ParaView-Versionen werden standardmäßig mit Qt5 gebaut, wenn dies nicht wie oben explizit ausgeschlossen wird. Wer schon auf Qt5 ist, der benötigt die folgenden Pakete, um zu kompilieren (getestet mit Ubuntu 18.04):
```
sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libxt-dev libgl1-mesa-dev
```

#### ParaView-Plugin kompilieren
* Im Flexi cmake das bauen des Plugins aktivieren:
```
FLEXI_BUILDPOSTI = ON
POSTI_BUILD_VISU = ON
POSTI_USE_PARAVIEW = ON
```
* paraview in der .bash_aliases/.zshrc/... umlenken auf:
```
 alias paraview='paraview --mpi'
```
* Mehrere Plugins kann man ganz einfach benutzen (z.b. bei mehreren `build` Ordnern), indem man auf der Commandline die `PV_PLUGIN_PATH`-Variable setzt. Zum Beispiel paraview mit folgendem Befehl ausführen um das Plugin im entsprechenden lib-Ordner zu verwenden:
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
