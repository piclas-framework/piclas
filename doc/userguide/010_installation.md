\hypertarget{installation}{}

# Installation

## Prerequisites
**PICLas** has been tested for various Linux distributions. This includes Ubuntu 14.04 LTS, 16.04 LTS and 18.04 LTS, OpenSUSE 42.1 and CentOS 7.
The suggested packages in this section can of course be replaced by self compiled versions.

The required packages for the Ubuntu Linux distributions are listed in table \ref{tab:installation_prereqs_ubuntu}. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git

| Package          | Ubuntu 14.04    | Ubuntu 16.04    | Ubuntu 18.04    |
|:----------------:|:---------------:|:---------------:|:---------------:|
| git              | x               |      x          |      x          |
| cmake            | x               |      x          |      x          |
| cmake-curses-gui | o               |      o          |      o          |
| liblapack3       | x               |      x          |      x          |
| liblapack-dev    | x               |      x          |      x          |
| gfortran         | x               |      x          |      x          |
| g++              | x               |      x          |      x          |
| mpi-default-dev  | x               |      x          |      x          |
| zlib1g-dev       | -               |      x          |      x          |
| exuberant-ctags  | o               |      o          |      o          |

Table: Debian/Ubuntu packages.\label{tab:installation_prereqs_ubuntu}
x: required, o: optional, -: not available

The required packages for OpenSUSE and CentOS are listed in table \ref{tab:installation_prereqs_redhat}.

Under OpenSUSE, packages are installed by the following command.

    sudo zypper install git

The `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/bin

Under CentOS, packages are installed by the following command.

    sudo yum install git

Additionally, the `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/openmpi/bin



| Package          | OpenSUSE 42.1 | CentOS 7 |
|:----------------:|:-------------:|:--------:|
| git              |      x        |    x     |
| cmake            |      x        |    x     |
| lapack-devel     |      x        |    x     |
| openmpi          |      x        |    x     |
| openmpi-devel    |      x        |    x     |
| zlib-devel       |      x        |    x     |
| gcc-fortran      |      x        |    x     |
| gcc              |      x        |    -     |
| gcc-c++          |      x        |    x     |
| ctags-etags      |      -        |    o     |

Table: OpenSUSE/CentOS packages.\label{tab:installation_prereqs_redhat}
x: required, o: optional, -: not available


On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active subroutines) in order to execute **PICLas** correctly. This is done using the command

~~~~~~~
ulimit -s unlimited
~~~~~~~

from the command line. For convenience, you can add this line to your `.bashrc`.

## Obtaining the source

The **PICLas** repository is available at GitHub. To obtain the most recent version you have two possibilities:

* Clone the **PICLas** repository from Github

        git clone https://github.com/piclas-framework/piclas.git

* Download **PICLas** from Github:

        wget https://github.com/piclas-framework/piclas/archive/master.tar.gz
        tar xzf master.tar.gz

Note that cloning **PICLas** from GitHub may not be possible on some machines, as e.g. the HLRS at the University of Stuttgart restricts internet access. Please refer to section \ref{sec:cloninghlrs} of this user guide.

## Compiling the code \label{sec:compilingthecode}

* Open a terminal
* Change into the **PICLas** directory
* Create a new subdirectory and use CMake to configure and compile the code

        mkdir build; cd build
        cmake ../
        make

For a list of all compiler options see Section \ref{sec:compileroptions}. The executables **PICLas** and **h5piclas2vtk** are contained in your **PICLas** directory in `build/bin/`.

### Directory paths \label{sec:installation_directory}

In the following, we write `$PICLASROOT` as a substitute for the path to the **PICLas** repository. Please replace `$PICLASROOT` in all following commands with the path to your **PICLas** repository *or* add an environment variable `$PICLASROOT`. 

Furthermore, the path to executables is omitted in the following, so for example, we write `piclas` instead of `$PICLASROOT/build/bin/piclas`. 

Here is some explanation for Linux beginners:

In order to execute a file, you have to enter the full path to it in the terminal. There are two different ways to enable typing `piclas` instead of the whole path (do not use both at the same time!)

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