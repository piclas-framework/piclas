## PICLas installation procedure

## Prerequisites

For **tested combinations** of prerequisities (HDF5, openMPI, cmake etc.) and known problems that may occur, see [BuildConfigurations.md](BuildConfigurations.md).

PICLas has been tested for various Linux distributions. This includes Ubuntu 14.04 LTS and 16.04 LTS, Chakra OS. In addition, Paraview or Tecplot can be used for visualization.

The required packages for the Ubuntu Linux distributions are listed in table \ref{tab:installation_prereqs_ubuntu}. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git
    

| Package          | Ubuntu 14.04    | Ubuntu 16.04    |
|:----------------:|:---------------:|:---------------:|
| git              | x               |      x          |
| cmake            | x               |      x          |
| cmake-curses-gui | x               |      x          |
| liblapack3       | x               |      x          |
| liplapack-dev    | x               |      x          |
| gfortran         | x               |      x          |
| g++              | x               |      x          |
|  mpi-default-dev | x               |      x          |
| zlib1g-dev       | -               |      x          |

Table: Required debian packages under Ubuntu.


## Compiling the code

* Open a terminal
* clone the GIT repository: 

		git clone git@gitlab.com:piclas/piclas.git

* Change into the PICLas directory
* Create a new subdirectory and use CMake to configure and compile the code
* Create a CMake build directory and choose the build options for PICLas

		mkdir -p build && cd build 
    ccmake ../

* Configure the build using CMake; ENABLE_ triggers the build of the corresponding tool.
* Build PICLas using CMake:

		make

The executable **piclas** is contained in your PICLas directory in `build/bin/`.

Custom configuration of compiler options may be done using

    ccmake ../

## Running the code

* Open a terminal
* Navigate to a directory and copy a case folder 

        cd temp

* Run piclas for PIC simulations

        $piclas parameter_piclas.ini

* Run piclas for DSMC simulations

        $piclas parameter_piclas.ini DSMCSpecies.ini

* Restart a simulation

        $piclas parameter_piclas.ini (DSMCSpecies.ini)  StateFile.h5

