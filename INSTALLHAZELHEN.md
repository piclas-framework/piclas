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
