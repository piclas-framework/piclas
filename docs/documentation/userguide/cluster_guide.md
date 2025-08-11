# Cluster Guidelines

PICLas has been tested on numerous different Linux- and CPU-based clusters and is continuously being developed to run efficiently on
high-performance systems.

## Simulating at HLRS

Unfortunately, the GitHub and GitLab servers are not available on machines at the HLRS due to
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


(sec:modules-at-hlrs)=
### Setting up the module environment

#### Setup for Vulcan

Following modules have been loaded before compiling and within the submit script before the execution of PICLas:

    module load cmake/3.22.3 hdf5/1.12.1-openmpi-4.1.2-gcc-11.2.0 mkl/2024.1
    module load gcc/13.1.0

#### Setup for Hunter (CPU only)

    module swap HLRS/APU HLRS/CPU
    module swap PrgEnv-cray PrgEnv-gnu/8.5.0
    module load cray-hdf5-parallel

(sec:petsc-at-hlrs)=
### Installing PETSc with restricted internet access

To utilize PETSc at the HLRS, you can download it manually to your local machine and upload to the cluster. In case the provided link does not work, check the [PETSc website](https://petsc.org/release/install/download/#doc-download) for a current download link:

    wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.21.6.tar.gz
    scp petsc-3.21.6.tar.gz vulcan:~/.
    ssh vulcan
    tar xf petsc-3.21.6.tar.gz
    cd petsc-3.21.6

You will have to manually download additional libraries for the PETSc installation. For that purpose, create a folder

    mkdir externalpackages

and execute the PETSc configuration to get shown the download links in the console

    ./configure --with-mpi --with-debugging=0 --with-shared-libraries=1 --with-mpi-f90module-visibility=0 --with-bison=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --with-packages-download-dir=externalpackages/ --download-hypre --download-mumps --download-scalapack --download-metis --download-parmetis

You should see a list with git and https links to the requested packages. After downloading the packages to your local machine and uploading them to the `externalpackages` folder, re-run the command to configure PETSc and follow the commands provided in the console, where the configuration is followed by a `make` and `make check` command. Afterwards, make sure to set the `PETSC_DIR` as provided at the end in your user profile and submit script:

    export PETSC_DIR=/path/to/petsc

If you are having problems with MPI, `--with-mpi-dir=` can be used to point directly to the MPI installation you want to use. Additionally, if you are compiling for a specific architecture, which might be different from your login node, make sure to specify `-march=` and `-mtune=` or other necessary commands.

#### Configure for Vulcan

The specific PETSc configure command for the Vulcan cluster and the AMD Epyc compute nodes (`genoa` queue):

    ./configure --with-mpi-dir=/opt/hlrs/non-spack/2022-02/mpi/openmpi/4.1.2-gcc-11.2.0/ --with-debugging=0 --with-shared-libraries=1 --with-mpi-f90module-visibility=0 --with-bison=0 COPTFLAGS='-O3 -march=znver3 -mtune=znver3' CXXOPTFLAGS='-O3 -march=znver3 -mtune=znver3' FOPTFLAGS='-O3 -march=znver3 -mtune=znver3' --with-packages-download-dir=../externalpackages/ --download-hypre --download-mumps --download-scalapack --download-metis --download-parmetis

#### Configure for Hunter

    ./configure --with-mpi --with-debugging=0 --with-shared-libraries=1 --with-mpi-f90module-visibility=0 --with-bison=0 COPTFLAGS='-O3 -march=znver4 -mtune=znver4' CXXOPTFLAGS='-O3 -march=znver4 -mtune=znver4' FOPTFLAGS='-O3 -march=znver4 -mtune=znver4' --with-packages-download-dir=externalpackages/ --download-hypre --download-mumps --download-scalapack --download-metis --download-parmetis

### Compiling and executing PICLas

After loading the required modules, the installation of PICLas can proceed as usual. To perform a calculation a job must be submitted, ideally using a submit script:

    qsub submit.sh

More information on using the queue system can be found in the [HLRS wiki](https://kb.hlrs.de/platforms/index.php/Platforms).

Section last updated: 16.06.2025
