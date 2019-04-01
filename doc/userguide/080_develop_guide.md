\hypertarget{develop_guide}{}

# Development guidelines \label{chap:develop_guide}

This chapter contains information about the development process and other issues concerning Git (GitLab/GitHub).

## Development process

Naming convention for branches, workflow for development, milestones etc.

After the successful completion of all regression checks (check-in, nightly, weekly), the master.dev branch can be merged into the master.

### Style Guide

* Why do we need a style guide?
    * It creates a unified appearance and coding structure
    * It makes the code more understandable and therefore important information is understood more
        easily
    * It forces the developers to think more actively about their work
* General rules
    * Coding language: English
    * A maximum of 132 characters are allowed per line (incl. Comments)
    * Indentation: 2 spaces (no tabs!)
    * Line breaks in comments -> the following line must be indented appropriately
    * Comments of modules and input-/output variables: Doxygen style
    * Comments of preprocessor directives in C-Style

#### Header of Functions and Subroutines

Always use `USE` statements with `ONLY`

    USE MODULE, ONLY: ...
    
this accounts for variables and function/subroutines. An exception are the initilization and finalization routines.

    !==============================================================
    !> \brief Fills the solution array U with a initial solution.
    !>
    !> Fills the solution array U with a initial solution provided by the ExactFunc subroutine through interpolation. Function is
    !> specified with the IniExactFunc paramter.
    !==============================================================
    SUBROUTINE FillIni(NLoc,xGP,U)
    !--------------------------------------------------------------
    ! MODULES
    USE MOD_PreProc
    USE MOD_Equation_Vars ,ONLY: IniExactFunc
    USE MOD_Exactfunc     ,ONLY: ExactFunc
    USE MOD_Mesh_Vars     ,ONLY: nElems
    IMPLICIT NONE
    !--------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    INTEGER,INTENT(IN)              :: NLoc                                    !< Polynomial degree of solution 
    REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)  !< Coordinates of Gauss-points
    REAL,INTENT(OUT)                :: U(PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
    !--------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER                         :: i,j,k,iElem
    !==============================================================
    
    ! Evaluate the initial solution at the nodes and fill the solution vector U. 
    DO iElem=1,nElems
      DO k=0,NLoc; DO j=0,NLoc; DO i=0,NLoc
        CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
    END SUBROUTINE FillIni

The separators `!====` and `!----` are exactly 132 characters long have been shortened for visualization purposes.

#### Variables

* Proprocessor variables: ```PP_$var```
     ```
     PP_nVar
     ```
     
* Counters: the counting variable (lower case) + description (the first character is capital case)
    ```
    DO iVar=1,PP_nVar
    ```
    
* Variables generally begin with a capital letter (composite words also)
    ```
    ALLOCATE(ActualElem)
    ```
    
* When using single characters: small at the beginning when using composite words otherwise in
  capital letters. Both is possible when purely single characters are used. Exceptions are allowed in
  special cases, but they are not recommened.
    ```
    hTilde, TildeH, (Elem%U)
    ```

#### Functions and Control Structures
* FORTRAN intrinsics generally in capital letters
    ```
    ALLOCATE(), DO, MAX(), SQRT(), INT(), etc.
    ```
* END-X is to be separated by a space
    ```
    END DO, END IF, END SUBROUTINE
    ```
* For loops and `IF` statements etc. comments are to be inserted at the end (and inbetween, e.g. when
`ELSE IF` is used)

    ```
    DO iVar=1,PP_nVar
      IF (a.EQ.b) THEN
    ...
      ELSE ! a.NE.b
    ...
      END IF ! a.EQ.b
    ...
    END DO ! PP_nVar
    ```

## Release and deploy

### Collaborative Numerics Group

The master branch of development group can be merged after the successful regression check with the master of the collaborative group. For this purpose, the collaborative repository can be added as a remote

    git remote add remote_name git@gitlab.com:collaborative-numerics-group/piclas/piclas.git

Now you can checkout the most recent version of the master branch of the collaborative-numerics-group and create a local branch with that version (a simple checkout will create a detached HEAD state)

    git fetch
    git checkout -b branch_name remote_name/master

The master branch of the development repository can now be merged into the newly created branch. Make sure to have the most recent version of the master branch (of the development repository) as well.

    git merge origin/master

Finally, the changes can be pushed from the *branch_name* to the master of collaborative-numerics-group

    git push remote_name master

If a tag has also been created, it should be pushed separately.

    git push remote_name tag_name

### GitHub

Upon completion of a milestone leading to tagged version, the tag should be deployed to GitHub.

## Simulating at the HLRS \label{sec:cloninghlrs}

Unfortunately, the GitHub and GitLab servers are not available on machines at the HLRS, such as the Hazelhen, due to restricted internet access. The workaround is to use ssh tunneling and remote forwarding to access the repositories.

### Cloning with the SSH protocol

You can use the SSH protocol to clone the repository. You have to connect to the cluster with the `RemoteForward` option

    ssh -R 7777:github.com:22 username@hazelhen.hww.de

To avoid using the above command every time, you can add the following to your `.ssh/config` file:

    host hlrs
       hostname hazelhen.hww.de
       user username
       RemoteForward 7777 gitlab.com:22

and login with `ssh hlrs`. Now you can clone the repository when logged onto the cluster by

    git clone ssh://git@localhost:7777/piclas/piclas.git

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

### Profiling with Craypad

* Compile PICLas with 
       module load perftools-base && module load perftools-lite && export CRAYPAT_LITE=event_profile
* Run PICLas with normal submit script
* Program has to finish normally! Enough time during execution. Note, that the profiled version is slower, hence, the testqueue is maybe too short. 
* Visualize the *.app2 files 

## Simulating at the forHLR \label{sec:forhlr}

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