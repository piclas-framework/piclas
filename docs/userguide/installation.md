# Installation

The following chapter describes the installation procedure on a Linux machine requiring root access. This includes the installation
of required prerequisites, setting up MPI and HDF5. Please note that high-performance clusters usually have a module environment,
where you have to load the appropriate modules instead of compiling them yourself. The module configuration for some of the clusters
used by the research group are given in Chapter {ref}`cluster_guide:Cluster Guidelines`.
In that case, you can jump directly to the description of the download and installation procedure of PICLas in Section
{ref}`sec:optaining-the-source`.

## Prerequisites
**PICLas** has been tested on various Linux distributions. This includes Ubuntu 14.04 LTS, 16.04 LTS and 18.04 LTS, 20.04 LTS
20.10 and 21.04, OpenSUSE 42.1 and CentOS 7.
For **tested combinations** of prerequisites (HDF5, OpenMPI, CMake etc.) and known problems that may occur, visit Chapter
{ref}`appendix:Appendix`.

The suggested packages in this section can be replaced by self compiled versions. The required packages for the Ubuntu Linux
distributions are listed in {numref}`tab:installation_prereqs_ubuntu` and for Centos Linux in
{numref}`tab:installation_prereqs_centos`.
Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git

```{table} Debian/Ubuntu packages. x: required, o: optional, -: not available
---
name: tab:installation_prereqs_ubuntu
---
|      Package     | Ubuntu 14.04 | Ubuntu 16.04 | Ubuntu 18.04 | Ubuntu 20.04 |
| :--------------: | :----------: | :----------: | :----------: | :----------: |
|        git       |       x      |       x      |       x      |       x      |
|       cmake      |       x      |       x      |       x      |       x      |
| cmake-curses-gui |       o      |       o      |       o      |       o      |
|    liblapack3    |       x      |       x      |       x      |       x      |
|   liblapack-dev  |       x      |       x      |       x      |       x      |
|     gfortran     |       x      |       x      |       x      |       x      |
|        g++       |       x      |       x      |       x      |       x      |
|  mpi-default-dev |       x      |       x      |       x      |       x      |
|    zlib1g-dev    |       -      |       x      |       x      |       x      |
|  exuberant-ctags |       o      |       o      |       o      |       o      |
```

and under CentOS via

    sudo yum install git

For extra packages install EPEL and SCL

    sudo yum install epel-release centos-release-scl

```{table} Centos packages. x: required, o: optional, -: not available
---
name: tab:installation_prereqs_centos
---
|      Package     |   CentOS 7   |
| :--------------: | :----------: |
|        git       |       x      |
|       cmake      |       x      |
|      cmake3      |       x      |
|      libtool     |       x      |
|   ncurses-devel  |       x      |
|   lapack-devel   |       x      |
|  openblas-devel  |       x      |
|   devtoolset-9   |       x      |
|        gcc       |       x      |
|      gcc-c++     |       x      |
|      zlib1g      |       x      |
|   zlib1g-devel   |       o      |
|  exuberant-ctags |       o      |
|   numactl-devel  |       x      |
|  rdma-core-devel |       o      |
```

On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active
subroutines) in order to execute **PICLas** correctly. This is done using the command

    ulimit -s unlimited

from the command line. For convenience, you can add this line to your `.bashrc`.

## Required Libraries

For a list of tested library version combinations, see Chapter {ref}`appendix:Appendix`.

### Installing/setting up GCC

Additional packages are required starting at specific versions of the GCC compiler suite.

|    GCC Version   | Ubuntu 20.04 (and older) |
| :--------------: |       :----------:       |
|       9.3.0      |        libmpfr-dev       |
|                  |        libmpc-dev        |


(sec:installing-mpi)=
### Installing/setting up OpenMPI

PICLas requires the installation of MPI with properly setup environment variables for parallel computations. The recommended MPI
version to use with PICLas is *OpenMPI 2.1.5*. An example of properly set environment variable are given below

    export MPI_DIR=/opt/openmpi/2.1.5
    export PATH="/opt/openmpi/2.1.5/bin:$PATH"
    export LD_LIBRARY_PATH="/opt/openmpi/2.1.5/lib:$LD_LIBRARY_PATH"

For convenience, you can add these lines to your `.bashrc`.

(sec:hdf5-installation)=
### Installing/setting up HDF5

An available installation of HDF5 can be utilized with PICLas. This requires properly setup environment variables and the
compilation of HDF5 during the PICLas compilation has to be turned off (`LIBS_BUILD_HDF5 = OFF`). If this option is enabled,
HDF5 will be downloaded and compiled. However, this means that every time a clean compilation of PICLas is performed, HDF5 will
be recompiled. It is preferred to either install HDF5 on your system locally or utilize the packages provided on your cluster.

The recommended HDF5 version to use with PICLas is **1.10.0-patch1**. In the following a manual installation of HDF5 is described,
if HDF5 is already available on your system you can skip to the next section {ref}`sec:setting-env-vars`.

#### Manual HDF5 installation

First, download HDF5 from [HDFGroup (external website)](https://portal.hdfgroup.org/display/support/Downloads) and extract it

    tar xvf hdf5-version.tar.gz

Then create a build folder

    cd hdf-version && mkdir -p build

and configure HDF5 to install into "/opt/hdf5/1.X.X" (your choice, should be writable)

    cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_ENABLE_F2003=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/hdf5/1.X.X ..

Make and install (if you chosen a folder required root access)

    make && make install

(sec:setting-env-vars)=
#### Setting environment variables

Depending whether HDF5 was installed using *configure* or *CMake*, different settings for the HDF5_DIR variable are required

* Configure

        export HDF5_DIR = /opt/hdf5/1.X.X/

* CMake

        export HDF5_DIR = /opt/hdf5/1.X.X/shared/cmake/XXX

If your CMake version is above 3.9.X, CMake uses a new findPackage routine, requiring that **HDF5_ROOT** is set

    export HDF5_ROOT=/opt/hdf5/1.10.0-patch1/mpi/

For convenience, you can add these lines to your `.bashrc`.

(sec:optaining-the-source)=
## Obtaining the source

The **PICLas** repository is available at GitHub. To obtain the most recent version you have two possibilities:

* Clone the **PICLas** repository from Github

        git clone https://github.com/piclas-framework/piclas.git

* Download **PICLas** from Github:

        wget https://github.com/piclas-framework/piclas/archive/master.tar.gz
        tar xzf master.tar.gz

Note that cloning **PICLas** from GitHub may not be possible on some machines, as e.g. the HLRS at the University of Stuttgart
restricts internet access. Please refer to section {ref}`sec:cloning-at-hlrs` of this user guide.

## Compiling the code

* Open a terminal
* Change into the **PICLas** directory
* Create a new subdirectory and use CMake to configure and compile the code

        mkdir build; cd build
        ccmake ..
        make

For a list of all compiler options visit Section {ref}`sec:compiler-options`. Finally, the executables **PICLas** and **piclas2vtk**
are contained in your **PICLas** directory in `build/bin/`.

(sec:directory-paths)=
### Directory paths

In the following, we write `$PICLASROOT` as a substitute for the path to the **PICLas** repository. Please replace `$PICLASROOT`
in all following commands with the path to your **PICLas** repository *or* add an environment variable `$PICLASROOT`.

Furthermore, the path to executables is omitted in the following, so for example, we write `piclas` instead of
`$PICLASROOT/build/bin/piclas`.

In order to execute a file, you have to enter the full path to it in the terminal. There are two different ways to enable typing
`piclas` instead of the whole path (do not use both at the same time!)

1. You can add an alias for the path to your executable. Add a command of the form

~~~~~~~
alias piclas='$PICLASROOT/build/bin/piclas'
~~~~~~~

to the bottom of the file `~/.bashrc`. Source your `~/.bashrc` afterwards with

~~~~~~~
. ~/.bashrc
~~~~~~~

2. You can add the **PICLas** binary directory to your `$PATH` environment variable by adding

~~~~~~~
export PATH=$PATH:$PICLASROOT/build/bin
~~~~~~~

to the bottom of the file `~/.bashrc` and sourcing your `~/.bashrc` afterwards.

Now you are ready for the utilization of **PICLas**.
