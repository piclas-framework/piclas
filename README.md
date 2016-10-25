# Welcome to Boltzplatz

[![build status](https://gitlabext.iag.uni-stuttgart.de/piclas/boltzplatz/badges/master/build.svg)](https://gitlabext.iag.uni-stuttgart.de/piclas/boltzplatz/builds/)


## Instructions to build Boltzplatz git version

The following steps are needed to build Boltzplatz. Be sure that you have installed the following prerequisites:

* MPI installation, e.g. [Open-MPI][openmpi]
* [CMake][cmake]

### Building steps

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

* Setting variables
         
         installed with config script: export HDF5_DIR = /opt/hdf5/1.X.X/
         installed with cmake: export HDF5_DIR = /opt/hdf5/1.X.X/shared/cmake/XXX

The tools are known to work with the following software versions:

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
 
### Building with ccmake on forhlr1

For building with *CMAKE* on the forhlr1-cluster, the following steps are needed:

* include in .bashrc or .profile:
  
        module load devel/cmake/3.5
        module swap compiler/intel compiler/intel/16.0
        module swap mpi/openmpi mpi/impi/5.1
        module load lib/hdf5/1.8_intel_16.0
        export HDF5_DIR=/opt/hdf5/1.8.17_intel_16.0_impi_5.1

* ccmake: as described above, but instead of building hdf5 include the following changes:

        set HDF5F90 to TRUE
        *IGNORE* possible error messages 

### Building with ccmake on hazelhen

For building on the hazelhen-cluster, the following steps are needed:

* include in .bashrc or .profile:
        
		module unload PrgEnv-cray
        module load PrgEnv-intel
        module load cray-hdf5-parallel
        module load tools/cmake

* Profiling with Craypad

    * Compile Boltzplatz with 
    
        module load perftools-base && module load perftools-lite && export CRAYPAT_LITE=event_profile
    
    * Run boltzplatz with normal submit script
    * Program has to finish normally! Enough time during execution. Note, that the profiled version is slower, hence, 
      the testqueue is maybe to short. 
    * Visualize the *.app2 files 


[openmpi]: https://www.open-mpi.org/
[paraview]: https://www.paraview.org
[cmake]: https://www.cmake.org
[hdf5]: https://www.hdfgroup.org/