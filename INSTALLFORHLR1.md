### Building with ccmake on forhlr1

For building with *CMAKE* on the forhlr1-cluster, the following modules (Intel compiler) should be loaded in the .bashrc or .profile:
  
    module load devel/cmake
    module load compiler/intel/18.0
    module load mpi/impi/2018
    module load lib/hdf5/1.10

Example submit script:

    #!/bin/bash
    #SBATCH --nodes=5
    #SBATCH --ntasks-per-node=20
    #SBATCH --time=04:00:00
    #SBATCH --job-name=PLACEHOLDER
    #SBATCH --partition multinode
    #SBATCH --mail-user=nizenkov@irs.uni-stuttgart.de
    #SBATCH --mail-type=ALL
    
    module load mpi/impi/2018
    module load lib/hdf5/1.10
    
    mpiexec.hydra -bootstrap slurm ./piclas parameter.ini DSMCSpecies.ini 1>log 2>log.err