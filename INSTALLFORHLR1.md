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

