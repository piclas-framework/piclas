# Cluster Guidelines

PICLas has been tested on numerous different Linux- and CPU-based clusters and is continuously been developed to run efficiently on
high-performance systems.

## Simulating at HLRS

Unfortunately, the GitHub and GitLab servers are not available on machines at the HLRS, such as the Hawk, due to 
restricted internet access. The workaround is to use ssh tunneling and proxy jump or remote forwarding to access 
the repositories.

(sec:cloning-at-hlrs)=
### Cloning with the SSH protocol

Two methods for checking out the code are described in the following.

#### Method 1 (Proxy Jump)
To clone a repository from, e.g., github.com on a HLRS system, the ssh proxy jump must first be set up. Simply connect to
the system via ssh and add the following lines to the ssh configuration file under `~/.ssh/config` on the cluster

    Host github.com
       HostName   github.com
       ProxyJump  RemoteHost

where the *RemoteHost* has internet access and can be accessed via ssh from the HLRS system. It is also defined

    Host RemoteHost
       User         username
       HostName     hostname.something.com
       ForwardX11   yes

Then simply clone the repository via

    git clone git@github.com:mygroup/myproject.git

and all the ssh traffic via `github.com` is automatically re-routed over the *RemoteHost*.

#### Method 2 (Remote Forwarding)

You can use the SSH protocol to clone the repository. You have to connect to the cluster with the `RemoteForward` option

    ssh -R 7777:github.com:22 username@hazelhen.hww.de

To avoid using the above command every time, you can add the following to your `.ssh/config` file:

    host hlrs
       hostname        hazelhen.hww.de
       user            username
       RemoteForward   7777 github.com:22

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
    aprun -n $ncores -N $nodecores -j 1 ./piclas parameter.ini DSMC.ini restart.h5 1>log 2>log.err

More information on using the queue system can be found in the [HLRS wiki](https://wickie.hlrs.de/platforms/index.php/CRAY_XC40_Using_the_Batch_System).

Section last updated: 27.03.2019

