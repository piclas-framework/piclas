\hypertarget{cluster_guide}{}

# Cluster guidelines \label{chap:cluster_guide}

## Simulating at HLRS \label{sec:cloninghlrs}

Unfortunately, the GitHub and GitLab servers are not available on machines at the HLRS, such as the Hawk, due to 
restricted internet access. The workaround is to use ssh tunneling and proxy jump or remote forwarding to access 
the repositories.

### Cloning with the SSH protocol

Two methods for checking out the code are described in the following.

#### Method 1 (Proxy Jump)
To clone a repository from, e.g., gitlab.com on a HLRS system, the ssh proxy jump must first be set up. Simply connect to
the system via ssh and add the following lines to the ssh configuration file under `~/.ssh/config` on the cluster

    Host gitlab.com
       HostName   gitlab.com
       ProxyJump  RemoteHost

where the *RemoteHost* has internet access and can be accessed via ssh from the HLRS system. It is also defined

    Host RemoteHost
       User         username
       HostName     hostname.something.com
       ForwardX11   yes

Then simply clone the repository via

    git clone git@gitlab.com:mygroup/myproject.git

and all the ssh traffic via `gitlab.com` is automatically re-routed over the *RemoteHost*.

#### Method 2 (Remote Forwarding)

You can use the SSH protocol to clone the repository. You have to connect to the cluster with the `RemoteForward` option

    ssh -R 7777:github.com:22 username@hazelhen.hww.de

To avoid using the above command every time, you can add the following to your `.ssh/config` file:

    host hlrs
       hostname        hazelhen.hww.de
       user            username
       RemoteForward   7777 gitlab.com:22

and login with `ssh hlrs`. Now you can clone the repository when logged onto the cluster by

    git clone ssh://git@localhost:7777/piclas/piclas.git

Note that if you experience problems when connecting, e.g., when the warning

    Warning: remote port forwarding failed for listen port 7777

is issued,
you have to choose a different port, e.g., 1827 (or any other 4-digit number) and re-connect.
If the code has already been cloned using the original port, the number of the port must be changed
in `./git/config` to the new number for git fetch/pull operations, which would then look like

    [remote "origin"]
      url = ssh://git@localhost:1827/piclas/piclas.git
      fetch = +refs/heads/*:refs/remotes/origin/*

### Cloning with the HTTPS protocol

The HLRS provides a tutorial for this case in their [https://wickie.hlrs.de](https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git). However, this method has not been verified.

### Compiling and executing PICLas

For building on the hazelhen cluster, certain modules have to be loaded and included in the .bashrc or .profile:

    module unload PrgEnv-cray
    module load PrgEnv-intel
    module load cray-hdf5-parallel

An example submit script for the test queue is then

    #!/bin/bash
    #PBS -N Testcase
    #PBS -M email@university.de
    #PBS -m abe
    #PBS -l nodes=4:ppn=24
    #PBS -l walltime=00:24:59
    #PBS -q test

    # number of cores per node
    nodecores=24

    # switch to submit dir
    cd $PBS_O_WORKDIR

    # get number of total cores
    ncores=$(cat $PBS_NODEFILE | wc -l)

    module unload craype-hugepages16M

    # restart
    aprun -n $ncores -N $nodecores -j 1 ./piclas parameter.ini DSMCSpecies.ini restart.h5 1>log 2>log.err 

More information on using the queue system can be found in the [HLRS wiki](https://wickie.hlrs.de/platforms/index.php/CRAY_XC40_Using_the_Batch_System).

Section last updated: 27.03.2019

### Profiling with Craypat

* Compile PICLas with 
       module load perftools-base && module load perftools-lite && export CRAYPAT_LITE=event_profile
* Run PICLas with normal submit script
* Program has to finish normally! Enough time during execution. Note, that the profiled version is slower, hence, the testqueue is maybe too short. 
* Visualize the *.app2 files 

## Simulating at forHLR \label{sec:forhlr}

For building with *CMake* on the forhlr1 cluster, the following modules (Intel compiler) should be loaded and included in the .bashrc or .profile:
  
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
    #SBATCH --mail-user=your@mail.de
    #SBATCH --mail-type=ALL
    
    module load mpi/impi/2018
    module load lib/hdf5/1.10
    
    mpiexec.hydra -bootstrap slurm ./piclas parameter.ini DSMCSpecies.ini 1>log 2>log.err

More information about the cluster and the batch system can be found at the [ForHLR wiki](https://wiki.scc.kit.edu/hpc/index.php/Category:ForHLR).

Section last updated: 27.03.2019

## Simulating at bwUniCluster \label{sec:bwuni}

For building with *CMake* on the bwUniCluster cluster, the following modules (Intel compiler) should be loaded and included in the .bashrc or .profile:
  
    module load devel/cmake
    module load compiler/intel/17.0
    module load mpi/openmpi/3.1-intel-17.0

During the configuration of PICLas with CMake, the option to build the HDF5 module has to be activated

    PICLAS_BUILD_HDF5 = ON

More information about the cluster and the batch system can be found at the [bwUniCluster wiki](https://www.scc.kit.edu/dienste/bwUniCluster.php).

Section last updated: 05.06.2019
