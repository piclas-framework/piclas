# Regression Testing

The purpose of regression testing is summarized by the following [Wikipedia quote](https://en.wikipedia.org/wiki/Regression_testing) {cite}`Basu2015`:

> *Regression testing (rarely non-regression testing) is re-running functional and non-functional tests to ensure that previously developed and tested software still performs after a change.*

## reggie2.0 Tool

PICLas is continuously tested by utilizing a Python based regression testing environment, which is
run by a [Gitlab Runner](https://docs.gitlab.com/runner/). Therefore, the tool *reggie2.0* is used, which is found under
[https://github.com/piclas-framework/reggie2.0](https://github.com/piclas-framework/reggie2.0).
Additionally, the different [Analyze routines](https://github.com/piclas-framework/reggie2.0#analyze-routines-for-analyzeini)
defined in the *analysis.ini* files that can be applied and the general structure of a regression test is described there.
Different tests are executed on check-in, during nightly or weekly testing. These tests are defined
in the file *.gitlab-ci.yml* that is located in the top level repository directory of PICLas.
In this file, various tests are defined, which are found under *regressioncheck* and a summary
of the different tests for PICLas are given [here](https://github.com/piclas-framework/piclas/blob/master/REGGIE.md).
The automatic execution by a [Gitlab Runner](https://docs.gitlab.com/runner/) can be performed on any machine that is connected to the
internet and in the following sections, the setup of such a machine is described.

## Local execution of reggie2.0
To quickly test regression checks locally, either to reproduce an error that has occurred during a GitLab pipeline or to test newly
developed code, the [reggie2.0](https://github.com/piclas-framework/reggie2.0) tool can be executed without the need to install a
[Gitlab Runner](https://docs.gitlab.com/runner/).

1. Install reggie from GitHub as described [here](https://github.com/piclas-framework/reggie2.0?tab=readme-ov-file#installation) and
   run the tool with `--help` to get an overview of the available options

       reggie --help

1. Build the required executable, e.g., *piclas*, either automatically using the [reggie2.0](https://github.com/piclas-framework/reggie2.0)
   tool or configure cmake and compile the executable by hand.

   Building the executable automatically requires a directory under *regressioncheck* that contains a *builds.ini* file from which
   all the compilation flags are read automatically and is described in the next step because automatic compilation is always
   accompanied by also running the code and analysing the results.

   There is also a *bash* script to extract a single *cmake* command line containing all the compiler settings from the *builds.ini* to
   compile a single executable by hand.
   Navigate to a *build* directory and run the script

       cd ~/piclas/build
       ./../tools/cmake-builds-ini.sh ../regressioncheck/NIG_convtest_poisson/builds.ini

   The output of the script will look like this

        CMAKE_BUILD_TYPE ........................ Debug,Release
        LIBS_BUILD_HDF5 ......................... OFF
        PICLAS_POLYNOMIAL_DEGREE ................ N
        PICLAS_EQNSYSNAME ....................... poisson
        PICLAS_TIMEDISCMETHOD ................... RK3
        LIBS_USE_MPI ............................ ON,OFF
        PICLAS_NODETYPE ......................... GAUSS
        PICLAS_PARTICLES ........................ ON,OFF
        PICLAS_CODE_ANALYZE ..................... ON,OFF
        LIBS_USE_PETSC .......................... OFF,ON

       Select the first set of options via

       cmake ..  -DCMAKE_BUILD_TYPE=Debug -DLIBS_BUILD_HDF5=OFF -DPICLAS_POLYNOMIAL_DEGREE=N -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=RK3 -DLIBS_USE_MPI=ON -DPICLAS_NODETYPE=GAUSS -DPICLAS_PARTICLES=ON -DPICLAS_CODE_ANALYZE=ON -DLIBS_USE_PETSC=OFF

   The last line can directly be copied into the terminal within a build directory to generate the make files for compiling PICLas with the
   first set of parameter options given in the *builds.ini* file (ignoring the *nocrosscombination* statements).
   Simply adjust the last line to have the correct flags set, execute *cmake* in the *build* directory and run *make* to compile:

       cmake ..  -DCMAKE_BUILD_TYPE=Debug -DLIBS_BUILD_HDF5=OFF -DPICLAS_POLYNOMIAL_DEGREE=N -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=RK3 -DLIBS_USE_MPI=ON -DPICLAS_NODETYPE=GAUSS -DPICLAS_PARTICLES=ON -DPICLAS_CODE_ANALYZE=ON -DLIBS_USE_PETSC=ON
       make -j

   Note that `-DLIBS_USE_PETSC=ON` has been adjusted in the above command.
   This will compile the *piclas* executable and it will be placed in the current directory under *bin*.
   Important note: Some regression tests build the hopr meshes "on-the-fly", hence, the *hopr* executable is required additionally.
   There are different possibilities to supply the *hopr* executable that is required:

   1. Set the compile flag `LIBS_DOWNLOAD_HOPR=ON`, which automatically downloads the *hopr* executable and places a link under
      *./bin* in the current build directory. This can be used if the *piclas* executable is compiled by hand.
   1. Set the environment variable via `export HOPR_PATH=/path/to/hopr` to point to the *hopr* executable on the current system

1. Run the [reggie2.0](https://github.com/piclas-framework/reggie2.0) tool either a) automatic mode or b) pre-compiled mode:

   This file is used to compile one or more *piclas* executables and the directories that accompany the *builds.ini* file will be used
   for testing. Note that not all executables might be used for those directories, as they might be excluded via the
   *excludebuild.ini* file.

   To run the automatic compilation and testing procedure, simply navigate the terminal to a location outside of the
   *regressioncheck* directory (which is found in the *piclas* repository).
   Switch to the *home* directory and run the [reggie2.0](https://github.com/piclas-framework/reggie2.0) tool there via

       cd ~
       reggie /path/to/piclas/regressioncheck/NIG_DSMC

   to start compiling and executing the resulting code.
   All output is placed under a new directory *output_dir* within the current directory.
   Never run the [reggie2.0](https://github.com/piclas-framework/reggie2.0) tool from within the */path/to/piclas/regressioncheck/*
   //as the directory tree structure is copied from there and the source path and target path cannot be the same!
   This procedure will run all the example directories under *NIG_DSMC*.

       ls /path/to/piclas/regressioncheck/NIG_DSMC

       2D_VTS_Distribution     builds.ini                  RotPeriodicBC                 SURF_PROB_DifferentProbs  VSS_VHS_SelfDiffusion
       Ambipolar_Diffusion     Macroscopic_Restart         RotPeriodicBCMulti            SURF_PROB_MultiReac
       Ambipolar_Diffusion_SF  MCC_BGG_Elec_XSec_Sampling  RotPeriodicBCMultiInterPlane  VirtualCellMerge

   To run the pre-compiled executable, navigate to the corresponding *build* directory and run the
   [reggie2.0](https://github.com/piclas-framework/reggie2.0) tool there

       cd ~/piclas/build
       reggie -e ./bin/piclas ../regressioncheck/NIG_DSMC/Ambipolar_Diffusion

   to run a specific test case, e.g., *Ambipolar_Diffusion*.

1. The *analysis.ini* file within the *Ambipolar_Diffusion* directory lists the analysis that is performed after the successful
   execution of *piclas*.
   An overview of the available analysis functions can be found [here](https://github.com/piclas-framework/reggie2.0?tab=readme-ov-file#analyze-routines-for-analyzeini).
   Also, look into the existing *regressioncheck* examples to get an idea how to construct a new *regressioncheck* setup or modify
   an existing one.


1. Reference files: They can be created automatically with the command line arguments

       -z, --rc              Create/Replace reference files that are required for analysis. After running the program, the output files are stored in the check-/example-directory.
       -i, --noMPI           Run program without "mpirun" (single thread execution).

   for analysis files, which use a reference file with which the output of a run is compared.
   Here, the flag `-i` is used to create the reference file with a single-core run, which is suggested as a best practice as the
   actual run might be performed with multiple cores and the output should ideally be the same.
   An example is given under [regressioncheck/WEK_DSMC/ChannelFlow_SurfChem_AdsorpDesorp_CO_O2](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/WEK_DSMC/ChannelFlow_SurfChem_AdsorpDesorp_CO_O2/analyze.ini)
   where *.h5* files are compared.

## Running GitLab *.gitlab-ci.yml* Tests

The GitLab CI/CD tests can either be run *locally* or *remotely* and both methods are explained in the following.
The tests are defined in the file *.gitlab-ci.yml* in the top-level directory of the piclas repository.

## Remote Testing on Gitlab

Open a browser and go to the [piclas gitlab pipelines website](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines), where the
latest pipeline jobs are displayed. To start a new pipeline, click the button *Run pipeline* and select the required branch name or
tag, which should be tested. Then, define the necessary *Variables*, which are summarized in {numref}`tab:pipeline_vars`.

```{table} Gitlab pipeline variables
---
name: tab:pipeline_vars
---
| Property                      | Description                                                                                                                              |  Value  |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | :-----: |
| DO_CHECKIN                    | short tests that are also run, when new commits are pushed                                                                               |    T    |
| DO_NIGHTLY                    | longer tests, executed every day                                                                                                         |    T    |
| DO_WEEKLY                     | very long tests, executed once a week                                                                                                    |    T    |
| DO_CODE_COVERAGE              | Generate coverage data of piclas for all jobs in current pipeline                                                                        |    T    |
| DO_REGGIE_COVERAGE            | Generate coverage data of reggie2.0 itself for all jobs in current pipeline                                                              |    T    |
| DO_NODE_SPLIT                 | MPI: virtual CPU splitting for multi-node testing, where a specific number of cores/threads are grouped in separate nodes (default is 2) |    T    |
| DO_CORE_SPLIT                 | MPI: virtual CPU splitting for multi-node testing, where each core/thread resembles a separate node                                      |    T    |
| DO_MPICH                      | MPI: Force compilation using MPICH instead of OpenMPI                                                                                    |    T    |
```

Per default, `DO_CHECKIN`, `DO_NIGHTLY`, `DO_WEEKLY`, `DO_NODE_SPLIT` and `DO_CORE_SPLIT` are tested automatically for the branch
*master.dev*. For details, see the [Pipeline schedules](https://piclas.boltzplatz.eu/piclas/piclas/-/pipeline_schedules) section in Gitlab.

## Code Coverage

Code coverage information can be inspected in multiple ways: on GitLab, by downloading GitLab artifacts and inspecting the data locally, or by running the tests locally in the first place. By default, all lines containing `Call Abort()` or `Call CollectiveStop()` are not considered when creating coverage reports.

### Coverage on GitLab

For regression testing on GitLab with code coverage, a separate flag `DO_CODE_COVERAGE` is required. This creates a separate stage, which is executed after all other stages. This way the code coverage stage is able to collect the coverage data from all (previous) jobs/runs of the current pipeline. This is done in the following way.

For each job (e.g., CHE_DSMC), the coverage data is collected and stored in a `.json` report file per build. These reports are named after the build and stored in the `Coverage` directory. If the build is reused for another job (e.g., same DSMC build for both CHE_DSMC and NIG_DSMC), the report file is updated. After all jobs in the other stages are finished, the coverage stage starts. The report files from all builds are combined into a single report containing the coverage data from all previously triggered runs/builds of the current pipeline. For example, if only `DO_CODE_COVERAGE` and `DO_CHECKIN` are set, the report will show coverage data from all builds with corresponding runs in the `DO_CHECKIN` case. The output of the coverage stage will also display which `.json` report files are found and therefore used to check whether the correct files and builds are considered.

#### Inspecting data on GitLab

The coverage stage will create a [GitLab report artifact](https://docs.gitlab.com/ci/yaml/artifacts_reports/#artifactsreportscoverage_report), which is used for the visualization. GitLab uses the last created report file of the current branch (if not expired yet) to display the coverage. Keep in mind that new pipelines might change the coverage data if a different set of tests is run. The coverage report is shown in the merge request difference view/changes. A line that was tested is indicated by a green bar to its left, otherwise a red bar appears. This allows inspection of regression tests for new features directly on GitLab. For smaller features/tests, it is recommended to check the coverage locally first to avoid triggering unnecessary tests.

To generate a full coverage report for all available regression tests of PICLas, either:
* Execute all tests in the same pipeline with `DO_CODE_COVERAGE`, or
* Combine separate reports manually

On GitLab the coverage is shown as a single number either in the output of the coverage job, on the right side when inspecting the job, or even in the merge request view. The displayed number is the line coverage (different coverage types in "Inspecting data locally"), which is set in `.gitlab-ci.yml`.

Currently the setup does not allow to test node/core split together, since the report can only be created per pipeline. Therefore, the reports have to be combined manually or inspected separately.

#### Inspecting data locally

In addition to the report artifact, the data is stored as a GitLab job artifact, which can be downloaded and inspected locally. This can be done on GitLab in the Pipeline overview on the right under Actions by choosing the Coverage artifact. To get a well-structured view of the data, it is recommended to use the `Coverage/html/coverage.html` file. Total coverage is displayed at the top right and split into three different categories: Lines, Functions and Branches. **Lines** checks how many of the compiled lines were actually executed. Lines that were covered by tests are highlighted in green, while untested lines appear in red. Yellow lines indicate some kind of branching, e.g. an if clause, where one condition was not tested. **Functions** works similarly, but only checks for function calls. So if a file contains only one function and it is called somewhere else, this will result in 100% coverage even if only 10% of lines were executed in this function. Lastly, **Branches** check how many different paths the program could take and then counts how many it did. So this mainly counts if/else statements or case select blocks.

### Local coverage tests

Since reggie supports coverage outputs, it can also be generated locally. To enable code coverage information, each executable must be compiled with additional flags via the `PICLAS_CODE_COVERAGE` option. This generates additional `.gcno` and `.gcda` files per object file, which track all compiled lines and the number of calls per line.

To enable code coverage when using reggie locally, use the `-v` option. More information can be found in the [reggie documentation](https://github.com/reggie-framework/reggie2.0).


#### Combining single reports

In some cases, it might be helpful to combine single reports of different reggie runs. This can be done using [gcovr](https://github.com/gcovr/gcovr), which is the same tool that reggie2.0 uses itself. For this case, `.json` files are used. Separate reports can be combined with
```
gcovr --root <root_dir> --add-tracefile <json_file1> --add-tracefile <json_file2> --html-nested report_name.html
```
where `--add-tracefile` takes wildcard arguments as well. The `--html-nested` option generates the same output format as reggie2.0 on GitLab, which is nicely structured analogous to the src directory. `--root` specifies the root directory containing the source files and `report_name.html` the output file.

Note that `<root_dir>` must be the same directory as the one used for the gcovr call that originally created the single report files. E.g. all single reports are created with `--root ~/some_dir`, then the `<root_dir>` for combining the reports must also be `--root ~/some_dir`. Reggie2.0 usually tries to append `src` to the found root directory to exclude UnitTests in the coverage information for piclas. Therefore combining reports generated by reggie2.0 for piclas would be done with
```
gcovr --root /path/to/piclas/src --add-tracefile Coverage/*.json --html-nested report_name.html
```

For more information on gcovr and coverage report formats, see the [gcovr documentation](https://gcovr.com/).

### Reggie coverage

Besides generating code coverage reports of PICLas, it is also possible to generate a report of the reggie tool itself. To do this, set `DO_REGGIE_COVERAGE`, which wraps each reggie call with the [Python coverage tool](https://coverage.readthedocs.io/). This generates a coverage report of all used lines in the reggie module, which is stored as a GitLab artifact. The report can be inspected using the `Coverage/reggie/index.html` file.

## Local Testing using *gitlab-ci-local*

To locally test the GitLab CI (including a YAML verification), [gitlab-ci-local](https://github.com/firecow/gitlab-ci-local) can be used.
An installation guide can be found [here](https://github.com/firecow/gitlab-ci-local#linux-based-on-debian).
After a successful installation, you can view the available parameters through
```
gitlab-ci-local --help
```
To view the stages for the default check-in pipeline, execute in the main folder of piclas:
```
gitlab-ci-local --list
```
To view all stages and tests:
```
gitlab-ci-local --list-all
```
To execute the check-in pipeline locally (that is the jobs that were shown with the `--list` command), use
```
gitlab-ci-local --shell-isolation
```
to avoid errors due to parallel writing of the ctags.txt file. An alternative is to limit the concurrent execution to one job, which
is analogous to the current configuration on the [Gitlab Runner](https://docs.gitlab.com/runner/) (requires gitlab-ci-local in version 4.42.0)
```
gitlab-ci-local --concurrency=1
```
It should be noted that currently the cache creation & utilization does not seem to represent the remote execution, meaning that some
errors might only be recognized after a push to the remote. A specific job can be executed simply by reference its name, and to also
consider the dependencies (i.e. the `needs:`), the following command can be utilized to execute, for example the DSMC check-in job:
```
gitlab-ci-local --needs CHE_DSMC
```
Another useful option to check the resulting configuration file is
```
gitlab-ci-local --preview preview.yml
```
which gives the expanded version of utilized `extends:` and `<<:` templates.

When running `gitlab-ci-local` on a system with a module environment, it is neccessary to pass the local modules that are used for compiling
```
DO_RUN_LOCAL="cmake/3.30.3   gcc/14.2.0   mpich/4.1.2/gcc/14.2.0    hdf5/1.14.0/gcc/14.2.0/mpich/4.1.2    hopr/master/gcc/14.2.0/mpich/4.1.2/hdf5/1.14.0    petsc/3.21.6/gcc/14.2.0/mpich/4.1.2"
gitlab-ci-local --variable DO_RUN_LOCAL=$DO_RUN_LOCAL
```
If multiple variables are required add them to the command
```
gitlab-ci-local --variable DO_RUN_LOCAL=$DO_RUN_LOCAL --variable CHECK_WARNINGS=True
```
to envoke additional options of the pipeline.

### Example
To run a specific reggie job, in this case a *weekly* reggie that depends on another job, the following parameters are passed
```
DO_RUN_LOCAL="cmake/3.30.3   gcc/14.2.0   mpich/4.1.2/gcc/14.2.0    hdf5/1.14.0/gcc/14.2.0/mpich/4.1.2    hopr/master/gcc/14.2.0/mpich/4.1.2/hdf5/1.14.0    petsc/3.21.6/gcc/14.2.0/mpich/4.1.2"
gitlab-ci-local --shell-isolation --needs WEK_Radiation --variable DO_RUN_LOCAL=$DO_RUN_LOCAL --variable DO_WEEKLY=T
```
where the arguments are listed and explained in the following table
```{table} gitlab-ci-local variables example
---
name: tab:gitlab_ci_local_vars
---
  | Parameter                             | Description                                                                         |
  | ------------------------------------- | ----------------------------------------------------------------------------------- |
  | --shell-isolation                     | Avoid errors due to parallel writing of the ctags.txt file                          |
  | --needs WEK_Radiation                 | Run WEK_Radiation, which also requires WEK_DSMC_Radiation                           |
  | --variable DO_RUN_LOCAL=$DO_RUN_LOCAL | Pass the locally installed and used modules                                         |
  | --variable DO_WEEKLY=T                | WEK_Radiation is a weekly reggie that requires the DO_WEEKLY to be passed           |
```

## Regression Test *Gitlab Runner* Setup for self-hosted Servers
This section describes the necessary steps to install a [Gitlab Runner](https://docs.gitlab.com/runner/) on a Ubuntu system to run *Gitlab Build Pipelines*.
In a first step, the required software packages for PICLas and [reggie2.0](https://github.com/piclas-framework/reggie2.0) are installed on a new system.
In a second step, the *gitlab-runner* program is installed and the setup of runner is described.

### Prerequisites: Installation of Software on Clean Ubuntu Setup (18.04)
Latest tests on

  * Ubuntu (18.04), 3 Jul 2019
  * Ubuntu server (18.04.3 LTS), 19 Nov 2019

The following packages can be installed automatically by using the script located at `./tools/Setup_ModuleEnv/InstallPackagesReggie.sh`.
The system inquiries can be skipped by forcing `yes` as input via

    yes | ./InstallPackagesRe

The script contains the following packages

```
# Check for updates
sudo apt-get update

# compiler
sudo apt-get install make cmake cmake-curses-gui gfortran g++ gcc

# python
sudo apt-get install python python-numpy python3 python3-numpy python3-matplotlib python-matplotlib

# editors
sudo apt-get install gedit vim vim-runtime vim-common vim-tiny gdb vim-gui-common

# git && versions
sudo apt-get install git gitg qgit subversion

# paraview
sudo apt-get install libpython-dev libboost-dev libphonon-dev libphonon4 libxt-dev mesa-common-dev
#apt-get install qt4-default qt4-dev-tools libqt4-dev qt4-qmake libqt4-opengl-dev
sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libgl1-mesa-dev

# Tecplot
sudo apt-get install libstdc++5

# tools
sudo apt-get install gzip gimp htop meld gnuplot gnuplot-x11 vlc okular ddd gmsh unzip
sudo apt-get install openvpn openssl openssh-client

# for FLEXI/PICLas
sudo apt-get install liblapack3 liblapack-dev zlib1g-dev exuberant-ctags

# for documentation
sudo apt-get install texlive-base
sudo apt-get install texlive-latex-extra

# hdf5-file viewer
sudo apt-get install hdfview

# Install libs for reggie
sudo apt-get install python-h5py

```

When no module environment is to be used on the server, the following packages are also required

```
# Further libs
sudo apt-get install hdf5-tools libhdf5-dev # this is maybe not required (do not install them if it works without these packages)


# openMPI
wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.3.tar.gz
tar -xvf openmpi-3.1.3.tar.gz
cd openmpi-3.1.3
mkdir -p build
cd build
../configure --prefix=/opt/openmpi/3.1.3
sudo make all install
export PATH="/opt/openmpi/3.1.3/bin:$PATH"
export LD_LIBRARY_PATH="/opt/openmpi/3.1.3/lib:$LD_LIBRARY_PATH"
cd ../..
rm openmpi-3.1.3.tar.gz


#HDF5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.bz2
tar -xjf hdf5-1.10.5.tar.bz2
cd hdf5-1.10.5
mkdir -p build
cd build
cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5
sudo make && sudo make install
export HDF5_DIR=/opt/hdf5/1.10.5/share/cmake
cd ../..
rm hdf5-1.10.5.tar.bz2
```

otherwise a module environment can be installed at this point, see
`~/piclas/tools/Setup_ModuleEnv/README.md`, which is explained in detail in Chapter {ref}`developerguide/tools:Developer Tools` under Section {ref}`developerguide/tools:Module Environment`.

When no module environment is to be used on the server, the following commands must be places in the
*.gitlab-ci.yml* file:

```
# Export the paths on new reggie2@reggie2 (no module env any more)
before_script:
  - ulimit -s unlimited
  - export PATH=/opt/openmpi/3.1.3/bin:$PATH
  - export LD_LIBRARY_PATH=/opt/openmpi/3.1.3/lib/:$LD_LIBRARY_PATH
  - export CMAKE_PREFIX_PATH=/opt/openmpi/3.1.3/share/cmake:$CMAKE_PREFIX_PATH
  - export CMAKE_LIBRARY_PATH=/opt/openmpi/3.1.3/lib:$CMAKE_LIBRARY_PATH
  - export HDF5_DIR=/opt/hdf5/1.10.5/share/cmake/
  - export PATH=/opt/hdf5/1.10.5/bin:$PATH
  - export LD_LIBRARY_PATH=/opt/hdf5/1.10.5/lib/:$LD_LIBRARY_PATH
  - export CMAKE_PREFIX_PATH=/opt/hdf5/1.10.5/:$CMAKE_PREFIX_PATH
  - export CMAKE_LIBRARY_PATH=/opt/hdf5/1.10.5/lib:$CMAKE_LIBRARY_PAT
```


otherwise, the correct environment must be loaded by adding the following in `/etc/profile`

    ```
    # Default modules
    module load gcc/9.2.0  cmake/3.15.3-d  openmpi/4.0.1/gcc/9.2.0  hdf5/1.10.5/gcc/9.2.0/openmpi/4.0.1
    ```

or by loading the modules directly in the gilab script file, e.g.,

    ```
    module load XX/XX
    ```

NOTE: The stack size limit has been removed here by `ulimit -s unlimited`, which might be required
by memory consuming programs

### Installation of Gitlab Runners
Latest tests on

  * Ubuntu (18.04) with gitlab-runner 10.5.0 (10.5.0), 3 Jul 2019
  * Ubuntu server (18.04.3 LTS) with gitlab-runner 10.5.0 (10.5.0), 19 Nov 2019

1. Install gitlab-runner from ubuntu packages (choose old version to avoid problems https://gitlab.com/gitlab-org/gitlab-runner/issues/1379)
   This creates the user gitlab-runner and a home directory (for 10.5 in /var/lib/gitlab-runner/)

    ```
    sudo apt-get install gitlab-runner
    ```

2. Start the runner program

    ```
    sudo gitlab-runner start
    ```
3. Register a runner using a shell executor (follow the information in the official guideline as it
   can vary slightly from version to version), on gitlab see `Settings` $\rightarrow$ `CI / CD Settings`
   $\rightarrow$ `Runners` (registration token during setup)
    ```
    sudo gitlab-runner register
    ```
4. Restart the runner
    ```
    sudo gitlab-runner restart
    ```
5. create ssh keys for normal user and set up password free access to gitlab (https://piclas.boltzplatz.eu)
    ```
    ssh-keygen -t ecdsa -b 521
    ```
    Add key to `Enabled deploy keys`. If multiple codes are on gitlab, add the key to one
    repository and select the key on the other repositories via `Privately accessible deploy keys`.

    Clone a code from each platform to create known hosts then
    ```
    sudo cp .ssh/*  /var/lib/gitlab-runner/.ssh
    ```
    Check rights in this folder. If necessary make gitlab-runner own the files with
    ```
    sudo chown -R gitlab-runner:gitlab-runner /var/lib/gitlab-runner/.ssh/
    ```
    If the runner is used to push to remote repositories, add the public key under *deploy keys*
    and execute, e.g.,
    ```
    sudo -u gitlab-runner git clone git@github.com:piclas-framework/piclas.git piclas_github
    ```
    to establish the first connection with the new repository and add the repo IP to the list of
    known hosts.
6. Start pipeline in gitlab or github for testing of reggie

NOTE: Interesting information is found in `/etc/systemd/system/gitlab-runner.service`.

### Configuration Files

The runner services can be adjusted by changing the settings in the file

    /etc/systemd/system/gitlab-runner.service

in which the runner configuration file is specified:

```
[Unit]
Description=GitLab Runner
After=syslog.target network.target
ConditionFileIsExecutable=/usr/bin/gitlab-runner

[Service]
StartLimitInterval=5
StartLimitBurst=10
ExecStart=/usr/bin/gitlab-runner "run" "--working-directory" "/var/lib/gitlab-runner/" "--config" "/etc/gitlab-runner/config.toml" "--service" "gitlab-runner" "--syslog" "--user" "gitlab-runner"


Restart=always
RestartSec=120

[Install]
WantedBy=multi-user.target
```

The runner configuration settings can be edited by changing the settings in the file

    /etc/gitlab-runner/config.toml

where the number of runners, the concurrency level and runner limits are specified:

```
concurrent = 2
check_interval = 0

[[runners]]
  name = "myrunner1"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]

[[runners]]
  name = "myrunner2"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]

[[runners]]
  name = "myrunner3"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]
```

### Automatic Deployment to Other Platforms (GitHub)

1. Add the required ssh key to the deploy keys on the respective platform (e.g. github)
1. Clone a code from the platform to update the list of known hosts. Do not forget to copy the
    information to the correct location for the runner to have access to the platform
    ```
    sudo cp~/.ssh/.ssh/known_hosts /var/lib/gitlab-runner/.ssh/known_hosts
    ```
    This might have to be performed via the gitlab-runner user, which can be accomplished by
    executing the following command
    ```
    sudo -u gitlab-runner git clone git@github.com:piclas-framework/piclas.git piclas_github
    ```
1. PICLas deployment is performed by the gitlab runner in the *deployment stage*
    ```
    github:
      stage: deploy
      tags:
        - withmodules-concurrent
      script:
        - if [ -z "${DO_DEPLOY}" ]; then exit ; fi
        - rm -rf piclas_github || true ;
        - git clone -b master --single-branch git@piclas.boltzplatz.eu:piclas/piclas.git piclas_github ;
        - cd piclas_github ;
        - git remote add piclas-framework git@github.com:piclas-framework/piclas.git ;
        - git push --force --follow-tags piclas-framework master ;
    ```

    This script clones the master branch of PICLas and deploys it on github.