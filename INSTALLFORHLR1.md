### Building with ccmake on forhlr1

For building with *CMAKE* on the forhlr1-cluster, the following steps are needed:

* include in .bashrc or .profile:
  
        module load devel/cmake/3.5
        module swap compiler/intel compiler/intel/17.0
        module swap mpi/openmpi mpi/impi/2018
        module load lib/hdf5/1.10_intel_17.0_impi_2017
        export HDF5_DIR=/software/all/lib/hdf5/1.10.2_intel_17.0_impi_2017

* ccmake: as described above, but instead of building hdf5 include the following changes:

        * set HDF5F90 to TRUE is deprecated and no longer needed
        *IGNORE* possible error messages 

