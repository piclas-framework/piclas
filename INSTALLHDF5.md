### Building HDF5

For a smooth building step ensure that the following environment variables are set (if the libraries are installed globally on the system):

* HDF5: HDF5_DIR (for CMake-build HDF5) or HDF5_ROOT (for HDF5 build with the configure script) or globally installed HDF5

The following libraries are installed locally within PICLas and build automatically (in case the PICLas needs  in case they are not found on the system):

* [HDF5][hdf5]: Library for efficient I/O on large scale systems (always needed)

  Manual HDF5-Installation

*  Download HDF5 from [www.support.hdfgroup.org](https://support.hdfgroup.org/downloads/index.html)

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

* CMAKE>3.9.+

    CMake uses a new findPackage, hence **HDF5_ROOT** is used, e.g.
    
         export HDF5_ROOT=/opt/hdf5/1.10.0-patch1/mpi/

