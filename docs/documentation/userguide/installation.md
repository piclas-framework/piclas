# Installation

The following chapter describes the installation procedure on a Linux machine requiring root access. This includes the installation
of required prerequisites, setting up MPI and HDF5. Please note that high-performance clusters usually have a module environment,
where you have to load the appropriate modules instead of compiling them yourself. The module configuration for some of the clusters
used by the research group are given in Chapter {ref}`userguide/cluster_guide:Cluster Guidelines`.
In that case, you can jump directly to the description of the download and installation procedure of PICLas in Section
{ref}`sec:obtaining-the-source`.

## AppImage executable download

PICLas and its tools can be installed on a Linux machine without the need of compiling the source code.
Currently, PICLas executables are only available as *[AppImage](https://appimage.org/)* for Linux.
The only requirements are that [GNU C Library (glibc)](https://www.gnu.org/software/libc/) and [OpenMPI](https://www.open-mpi.org/)
are pre-installed on the target system and available when running the AppImage executables.
The minimum version requirements are listed in the following table and it is not guaranteed that more recent versions of the
libraries listed in the table work automatically

| PICLas Version |                              MPI                              | glibc |
| :------------: | :-----------------------------------------------------------: | :---: |
| 3.6.0 - 4.0.0  | [openmpi-4.1.6](https://www.open-mpi.org/software/ompi/v4.1/) | 2.18  |
| 3.3.0 - 3.5.0  | [openmpi-4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) | 2.18  |
|     <3.3.0     | [openmpi-4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) | 2.17  |

Static libraries for [OpenMPI](https://www.open-mpi.org/) are not distributed within the AppImage package because of the system-dependent optimizations
(e.g. specific InfiniBand settings).
Additional external libraries and versions that are used for compiling are the following but do not have to be installed separately
on the system where the AppImage is going to be executed

| PICLas Version |                                 GNU GCC                                  |                                          HDF5                                          |                         PETSc                          |
| :------------: | :----------------------------------------------------------------------: | :------------------------------------------------------------------------------------: | :----------------------------------------------------: |
| 3.6.0 - 4.0.0  | [gcc (GCC) 8.3.1 20190311 (Red Hat 8.3.1-3)](https://gcc.gnu.org/gcc-8/) | [HDF5 1.12.2](https://www.hdfgroup.org/2022/04/release-of-hdf5-1-12-2-newsletter-183/) | [PETSc 3.21.6](https://petsc.org/release/changes/321/) |
|    <=3.5.0     | [gcc (GCC) 8.3.1 20190311 (Red Hat 8.3.1-3)](https://gcc.gnu.org/gcc-8/) | [HDF5 1.12.2](https://www.hdfgroup.org/2022/04/release-of-hdf5-1-12-2-newsletter-183/) | [PETSc 3.18.4](https://petsc.org/release/changes/318/) |

Other operating systems, such as Windows or MacOS might be supported in the future.

Download the pre-compiled executables from the [PICLas release tag assets](https://github.com/piclas-framework/piclas/releases).
Note that versions prior to v3.0.0 are not supported for AppImage download.
Unzip the files, switch into the directory an then and check their MD5 hashes by running

    md5sum -c md5sum.txt

which should produce output looking like

    piclasDSMC: OK
    piclasLeapfrogHDG: OK
    piclas2vtk: OK
    superB: OK

indicating that everything is fine.
After downloading the binary files, it has to be checked that all files are executable and if not simply run

    chmod +x piclas*

for all files beginning with piclas (add the other files to the list if required) before they can be used.
If problems occur when executing the AppImage, check the [troubleshooting]( https://docs.appimage.org/user-guide/troubleshooting/index.html)
section for possible fixes.

## Prerequisites
**PICLas** has been used on various Linux distributions in the past. This includes different (K)Ubuntu version up to 24.04, OpenSUSE 42.1 and CentOS 7. For a list of tested library version combinations, see Section {ref}`sec:required-libs`. The suggested packages in this section can be replaced by self compiled versions.

### Ubuntu or similar

 The required packages for the Ubuntu or similar Linux distributions can be obtained using the apt environment:

    sudo apt-get install git cmake cmake-curses-gui liblapack3 liblapack-dev gfortran g++ mpi-default-dev zlib1g-dev

The following packages are optional:

    sudo apt-get install exuberant-ctags ninja pkg-config

Additionally, the most recent versions of the GCC compiler (>13) require the following packages:

    sudo apt-get install libmpfr-dev libmpc-dev libgmp-dev

### CentOS

For CentOS, the packages can be installed using the following command

    sudo yum install git cmake cmake3 libtool ncurses-devel lapack-devel openblas-devel devtoolset-9 gcc gcc-c++ zlib1g zlib1g-devel exuberant-ctags numactl-devel rdma-core-devel binutils tar epel-release centos-release-scl

On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active
subroutines) in order to execute **PICLas** correctly. This is done using the command

    ulimit -s unlimited

from the command line. For convenience, you can add this line to your `.bashrc`.

(sec:required-libs)=
## Required Libraries
The following list contains the **recommended library combinations** for the Intel and GNU compiler in combination with HDF5, OpenMPI, CMake etc.

| PICLas Version | CMake  | Compiler  |            MPI             |  HDF5  | PETSc  |
| :------------: | :----: | :-------: | :------------------------: | :----: | :----: |
|     4.1.0      | 4.2.1  | gcc14.2.0 | openmpi-5.0.8, mpich-4.3.1 | 1.14.6 | 3.22.5 |
| 3.6.0 - 4.0.0  | 3.31.1 | gcc14.2.0 | openmpi-5.0.6, mpich-4.1.2 | 1.14.0 | 3.21.6 |
|     3.4.0      | 3.31.1 | gcc13.2.0 | openmpi-4.1.5, mpich-4.1.2 | 1.14.0 | 3.19.3 |
|     3.3.0      | 3.26.4 | gcc13.2.0 |       openmpi-4.1.5        | 1.14.0 | 3.19.3 |
|     2.8.0      | 3.24.2 | gcc12.2.0 |       openmpi-4.1.4        | 1.12.2 |   -    |
| 2.3.0 - 2.7.0  | 3.21.3 | gcc11.2.0 |       openmpi-4.1.1        | 1.12.1 |   -    |
|     2.0.0      |  3.17  | intel19.1 |          impi2019          |  1.10  |   -    |
| 2.0.0 - 2.2.2  |  3.17  | intel19.1 |          impi2019          |  1.10  |   -    |

and the **minimum requirements**

| PICLas Version | Compiler  |  HDF5  |      MPI      | CMake |
| :------------: | :-------: | :----: | :-----------: | :---: |
| 2.3.0 - 2.8.0  | gnu9.2.0  | 1.10.6 | openmpi-3.1.6 | 3.17  |
| 2.0.0 - 2.2.2  | intel18.1 |  1.10  |   impi2018    | 3.17  |

A full list of all previously tested combinations is found in Chapter {ref}`userguide/appendix:Appendix`. Alternative combinations might work as well, however, have not been tested.

If you are setting-up a fresh system for the simulation with PICLas, it is recommended to use a Module environment, which can be set
up with the provided shell scripts in `piclas/tools/Setup_ModuleEnv`.
A description is available here: `piclas/tools/Setup_ModuleEnv/README.md`.
This allows installing and switching between different compiler, MPI, HDF5 and PETSc versions.

### Installing GCC

You can install a specific version of the GCC compiler from scratch, if for example the compiler provided by the repositories of your operating systems is too old. If you already have a compiler installed, make sure the environment variables are set correctly as shown at the end of this section. For this purpose, download the GCC source code directly using a mirror website (detailed information and different mirrors can be found here: [GCC mirror sites](https://gcc.gnu.org/mirrors.html))

    wget "ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/gcc-11.2.0/gcc-11.2.0.tar.gz"

After a successful download, make sure to unpack the archive

    tar -xzf gcc-11.2.0.tar.gz

Now you can enter the folder, create a build folder within and enter it

    cd gcc-11.2.0 && mkdir build && cd build

Here, you can configure the installation and create the Makefile. Make sure to adapt your installation path (`--prefix=/home/user/gcc/11.2.0`) before continuing with the following command

    ../configure -v \
    --prefix=/home/user/gcc/11.2.0 \
    --enable-languages=c,c++,objc,obj-c++,fortran \
    --enable-shared \
    --disable-multilib \
    --disable-bootstrap \
    --enable-checking=release \
    --with-sysroot=/ \
    --with-system-zlib

After the Makefile has been generated, you can compile the compiler in parallel, e.g. with 4 cores

    make -j4 2>&1 | tee make.out

and finally install it by

    make install 2>&1 | tee install.out

If you encounter any difficulties, you can submit an issue on GitHub and attach the `make.out` and/or `install.out` files, where the console output is stored. To make sure that the installed compiler is also utilized by CMake, you have to set the environment variables, again making sure to use your installation folder and the correct version

    export CC = /home/user/gcc/11.2.0/bin/gcc
    export CXX = /home/user/gcc/11.2.0/bin/g++
    export FC = /home/user/gcc/11.2.0/bin/gfortran
    export PATH="/home/user/gcc/11.2.0/include/c++/11.2.0:$PATH"
    export PATH="/home/user/gcc/11.2.0/bin:$PATH"
    export LD_LIBRARY_PATH="/home/user/gcc/11.2.0/lib64:$LD_LIBRARY_PATH"

For convenience, you can add these lines to your `.bashrc`.

(sec:installing-mpi)=
### Installing OpenMPI

You can install a specific version of OpenMPI from scratch, using the same compiler that will be used for the code installation. If you already have OpenMPI installed, make sure the environment variables are set correctly as shown at the end of this section. First, download the OpenMPI source code either from the website [https://www.open-mpi.org/](https://www.open-mpi.org/) or directly using the following command:

    wget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz"

Unpack it, enter the folder, create a build folder and enter it (procedure analogous to the GCC installation). For the configuration, utilize the following command, which specifies the compiler to use and the installation directory with `--prefix=/home/user/openmpi/4.1.1`

    ../configure --prefix=/home/user/openmpi/4.1.1 CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)

After the Makefile has been generated, you can compile OpenMPI in parallel, e.g. with 4 cores

    make -j4 2>&1 | tee make.out

and finally install it by

    make install 2>&1 | tee install.out

If you encounter any difficulties, you can submit an issue on GitHub and attach the `make.out` and/or `install.out` files, where the console output is stored. To make sure that the installed OpenMPI is also utilized by CMake, you have to set the environment variables, again making sure to use your installation folder and the correct version

    export MPI_DIR=/home/user/openmpi/4.1.1
    export PATH="/home/user/openmpi/4.1.1/bin:$PATH"
    export LD_LIBRARY_PATH="/home/user/openmpi/4.1.1/lib:$LD_LIBRARY_PATH"

For convenience, you can add these lines to your `.bashrc`.

(sec:hdf5-installation)=
### Installing HDF5

An available installation of HDF5 can be utilized with PICLas. This requires properly setup environment variables and the compilation of HDF5 during the PICLas compilation has to be turned off (`LIBS_BUILD_HDF5 = OFF`, as per default). If this option is enabled, HDF5 will be downloaded and compiled. However, this means that every time a clean compilation of PICLas is performed, HDF5 will be recompiled. It is recommended to install HDF5 within the Module environment or utilize the packages provided on your cluster.

You can install a specific version of HDF5 from scratch, using the same compiler and MPI configuration that will be used for the code installation. If you already have HDF5 installed, make sure the environment variables are set correctly as shown at the end of this section. First, download HDF5 from [HDFGroup (external website)](https://portal.hdfgroup.org/display/support/Downloads) or download it directly

    wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.gz"

and extract it, enter the folder, create a build folder and enter it

    tar -xvf hdf5-1.12.1.tar.gz && cd hdf5-1.12.1 && mkdir -p build && cd build

Afterwards, configure HDF5 and specify the installation directory with `--prefix=/home/user/hdf5/1.12.1`

    ../configure --prefix=/home/user/hdf5/1.12.1 --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel CC=$(which mpicc) CXX=$(which mpicxx) FC=$(which mpifort)

After the Makefile has been generated, you can compile HDF5 in parallel, e.g. with 4 cores

    make -j4 2>&1 | tee make.out

and finally install it by

    make install 2>&1 | tee install.out

If you encounter any difficulties, you can submit an issue on GitHub and attach the `make.out` and/or `install.out` files, where the console output is stored. To make sure that the installed HDF5 is also utilized by CMake, you have to set the environment variables, again making sure to use your installation folder and the correct version

    export HDF5_DIR = /home/user/hdf5/1.12.1
    export HDF5_ROOT = /home/user/hdf5/1.12.1
    export PATH="/home/user/hdf5/1.12.1/include:$PATH"
    export PATH="/home/user/hdf5/1.12.1/bin:$PATH"
    export LD_LIBRARY_PATH="/home/user/hdf5/1.12.1/lib:$LD_LIBRARY_PATH"

For convenience, you can add these lines to your `.bashrc`.

(sec:petsc-installation)=
### Installing PETSc

The following list contains the **recommended/working library versions** for PETSc and PICLas

| PICLas Version | PETSc Version |
| :------------: | :-----------: |
|     3.6.0      |    3.21.6     |
|     3.3.0      |    3.19.3     |
|     3.0.0      |  3.17, 3.18   |
|     2.9.0      |     3.17      |

To utilize PETSc, it is recommended to install it with the Module environment script (`tools/Setup_ModuleEnv/InstallPETSc.sh`) or the built-in installation during PICLas configuration (see Section {ref}`sec:solver-settings`). The recommended external libraries for PETSc are MUMPS, Hypre, ScaLAPACK, METIS, and ParMETIS, which are automatically downloaded when PETSc is installed with the mentioned approaches.

On a cluster with restricted internet access, the procedure described above will not work. However, you can build it directly on the cluster by downloading locally and manually uploading PETSc to your cluster. More information is provided specifically for the servers of the HLRS in Section {ref}`sec:petsc-at-hlrs`, which can be applied to other systems.

(sec:obtaining-the-source)=
## Obtaining the source

The **PICLas** repository is available at GitHub. To obtain the most recent version you have two possibilities:

* Clone the **PICLas** repository from Github

    ````
    git clone https://github.com/piclas-framework/piclas.git
    ````

* Download **PICLas** from Github:

    ````
    wget https://github.com/piclas-framework/piclas/archive/master.tar.gz
    tar xzf master.tar.gz
    ````

Note that cloning **PICLas** from GitHub may not be possible on some machines, as e.g. the HLRS at the University of Stuttgart
restricts internet access. Please refer to section {ref}`sec:cloning-at-hlrs` of this user guide.

## Compiling the code

To compile the code, open a terminal, change into the **PICLas** directory and create a new subdirectory. Use the CMake curses interface (requires the `cmake-curses-gui` package on Ubuntu) to configure and generate the Makefile

    mkdir build && cd build
    ccmake ..

Here you will be presented with the graphical user interface, where you can navigate with the arrow keys.
Before you can generate a Makefile, you will have to configure with `c`, until the generate option `g` appears at the bottom of the terminal. Alternatively, you can configure the code directly by supplying the compiler flags

    cmake -DPICLAS_TIMEDISCMETHOD=DSMC ..

For a list of all compiler options visit Section {ref}`sec:compiler-options`.
PICLas supports Unix Makefiles (default) and [Ninja](https://ninja-build.org/) as generators.
To select ninja either export the environment variable `export CMAKE_GENERATOR=Ninja` or add the cmake option `-GNinja`, e.g.

    cmake -GNinja ..
    ccmake -GNinja ..

After a successful generation of the Makefile (or ninja build files), compile the code in parallel (e.g. with 4 cores) using

    make -j4 2>&1 | tee piclas_make.out

or the corresponding Ninja command

    ninja -j4 2>&1 | tee piclas_make.out

To use all available cores for parallel compilation use either `make -j` or `ninja -j0`.
Finally, the executables **piclas** and **piclas2vtk** are contained in your **PICLas** directory in `build/bin/`.
If you encounter any difficulties, you can submit an issue on GitHub and attach the `piclas_make.out` file, where the console output is stored.

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