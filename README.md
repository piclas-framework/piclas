# Welcome to Boltzplatz

## Instructions to build Boltzplatz git version

The following steps are needed to build Boltzplatz. Be sure that you have installed the following prerequisites:

* MPI installation, e.g. [Open-MPI][openmpi]
* [CMake][cmake]

### Builing steps

For a smooth building step ensure that the following environement variables are set (if the libraries are installed globally on the system):

* HDF5: HDF5_DIR (for CMake-build HDF5) or HDF5_ROOT (for HDF5 build with the configure script) or globally installed HDF5
* Lapack/Blas library
* GIT

The following libraries are installed locally within Boltzplatz and build automatically (in case the Boltzplatz needs  in case they are not found on the system):

* [HDF5][hdf5]: Library for efficient I/O on large scale systems (always needed)

  Manual HDF5-Installation

*  Download HDF5

         tar xvf hdf5-version.tar.gz
*  Create Build folder

         cd hdf-version && mkdir -p build
* Configure and make HDF5

         cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_ENABLE_F2003=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/hdf5/1.X.X ..
* Make 
         make && make install

* Setting of the variables
         
         export HDF5_ROOT = /opt/hdf5/1.X.X/

* [Tecplot IO Library][tecio]: Library to write the proprietary Tecplot binary data format

  Manual Tecio-installation
  
* Create new folder and clone  TECIO

         git clone git@129.69.43.151:libs/TECPLOT.git  TECPLOT
         
* Unzip archive with

         cd TECPLOT && tar xvf tecio-2013.tar.bz2
        
* Build tecio

         cd tecio-2013 && ./Runmake linuxg27x64.24 -tecio 

* Move tecio to be found by cmake

         mkdir -p /opt/tecio-2013
         mv * /opt/tecio-2013
   

The Tools are known to work with the following software versions:

* OpenMPI > 1.8
* HDF5 > 1.8.13


To build the Boltzplatz the following steps are needed:

* clone the GIT repository: 

		git clone git@gitlab.iag.uni-stuttgart.de:piclas/boltzplatz.git

* Create a CMake build directory and choose the build options for Boltzplatz

		mkdir -p build && cd build && ccmake ../

* Configure the build using CMake; ENABLE_ triggers the build of the corresponding tool.
* Build Boltzplatz using CMake:

		make

* Now Boltzplatz is build in the chosen CMake build directory and can be used. You may choose to install them to a different location by setting the CMAKE_INSTALL_DIRECTORY and perform:

		make install

* After this step the code can be used.

[openmpi]: https://www.open-mpi.org/
[paraview]: https://www.paraview.org
[cmake]: https://www.cmake.org
[hdf5]: https://www.hdfgroup.org/
[tecio]: http://www.tecplot.com/downloads/tecio-library/