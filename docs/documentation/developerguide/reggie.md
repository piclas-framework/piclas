# Regression Testing

The purpose of regression testing is summarized by the following [Wikipedia quote](https://en.wikipedia.org/wiki/Regression_testing) {cite}`Basu2015`:

> *Regression testing (rarely non-regression testing) is re-running functional and non-functional tests to ensure that previously developed and tested software still performs after a change.*

## reggie Tool

PICLas is continuously tested by utilizing a Python based regression testing environment, which is
run by a [Gitlab Runner](https://docs.gitlab.com/runner/). Therefore, the tool *reggie* is used, which is found under
[https://github.com/reggie-framework/reggie](https://github.com/reggie-framework/reggie).
Additionally, the different [Analyze routines](https://github.com/reggie-framework/reggie#analyze-routines-for-analyzeini)
defined in the *analysis.ini* files that can be applied and the general structure of a regression test is described there.
Different tests are executed on check-in, during nightly or weekly testing. These tests are defined
in the file *.gitlab-ci.yml* that is located in the top level repository directory of PICLas.
In this file, various tests are defined, which are found under *regressioncheck* and a summary
of the different tests for PICLas are given [here](https://github.com/piclas-framework/piclas/blob/master/REGGIE.md).
The automatic execution is performed by a [Gitlab Runner](https://docs.gitlab.com/runner/) on a machine that is connected to the internet.

## Local execution of reggie
To quickly test regression checks locally, either to reproduce an error that has occurred during a GitLab pipeline or to test newly
developed code, the [reggie](https://github.com/reggie-framework/reggie) tool can be executed without the need to install a
[Gitlab Runner](https://docs.gitlab.com/runner/).

1. Install reggie from GitHub as described [here](https://github.com/reggie-framework/reggie?tab=readme-ov-file#installation) and
   run the tool with `--help` to get an overview of the available options

       reggie --help

1. Build the required executable, e.g., *piclas*, either automatically using the *reggie*
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
   Important note: Some regression tests build the **PyHOPE** meshes "on-the-fly", hence, the *pyhope* tool is required additionally.
   **pyhope** needs to be pre-installed on the system and an installation guide can be found on
   [GitHub - PyHOPE](https://github.com/hopr-framework/PyHOPE).

1. Run the *reggie* tool either a) automatic mode or b) pre-compiled mode:

   This file is used to compile one or more *piclas* executables and the directories that accompany the *builds.ini* file will be used
   for testing. Note that not all executables might be used for those directories, as they might be excluded via the
   *excludebuild.ini* file.

   To run the automatic compilation and testing procedure, simply navigate the terminal to a location outside of the
   *regressioncheck* directory (which is found in the *piclas* repository).
   Switch to the *home* directory and run the *reggie* tool there via

       cd ~
       reggie /path/to/piclas/regressioncheck/NIG_DSMC

   to start compiling and executing the resulting code.
   All output is placed under a new directory *output_dir* within the current directory.
   Never run the *reggie* tool from within the */path/to/piclas/regressioncheck/*
   //as the directory tree structure is copied from there and the source path and target path cannot be the same!
   This procedure will run all the example directories under *NIG_DSMC*.

       ls /path/to/piclas/regressioncheck/NIG_DSMC

       2D_VTS_Distribution     builds.ini                  RotPeriodicBC                 SURF_PROB_DifferentProbs  VSS_VHS_SelfDiffusion
       Ambipolar_Diffusion     Macroscopic_Restart         RotPeriodicBCMulti            SURF_PROB_MultiReac
       Ambipolar_Diffusion_SF  MCC_BGG_Elec_XSec_Sampling  RotPeriodicBCMultiInterPlane  VirtualCellMerge

   To run the pre-compiled executable, navigate to the corresponding *build* directory and run the
   *reggie* tool there

       cd ~/piclas/build
       reggie -e ./bin/piclas ../regressioncheck/NIG_DSMC/Ambipolar_Diffusion

   to run a specific test case, e.g., *Ambipolar_Diffusion*.

1. The *analysis.ini* file within the *Ambipolar_Diffusion* directory lists the analysis that is performed after the successful
   execution of *piclas*.
   An overview of the available analysis functions can be found [here](https://github.com/reggie-framework/reggie?tab=readme-ov-file#analyze-routines-for-analyzeini).
   Also, look into the existing *regressioncheck* examples to get an idea how to construct a new *regressioncheck* setup or modify
   an existing one.


1. Reference files: They can be created automatically with the command line arguments

       -z, --rc              Create/Replace reference files that are required for analysis. After running the program, the output files are stored in the check-/example-directory.
       -i, --noMPI           Run program without "mpirun" (single thread execution).

   for analysis files, which use a reference file with which the output of a run is compared.
   Here, the flag `-i` is used to create the reference file with a single-core run, which is suggested as a best practice as the
   actual run might be performed with multiple cores and the output should ideally be the same.
   An example is given under [regressioncheck/WEK_DSMC/SurfChem_AdsorpDesorp_CO](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/WEK_DSMC/SurfChem_AdsorpDesorp_CO/analyze.ini)
   where *.h5* files are compared.

## Compression of HDF5 reference files
When creating .h5 reference files that are used for comparison in regression tests, it is very beneficial to compress these to a minimum
before committing and pushing them to the git repository.
Consider the following example with an `analysis.ini` file containing

    ! hdf5 diff
    h5diff_file            = TestRotatingWall_DSMCSurfState_000.00100000000000000.h5
    h5diff_reference_file  = TestRotatingWall_DSMCSurfState_000.00100000000000000_reference.h5
    h5diff_data_set        = SurfaceData
    h5diff_tolerance_value = 40E-2
    h5diff_tolerance_type  = relative
    h5diff_var_attribute   = VarNamesSurface
    h5diff_var_name        = Total_TorqueZ
    h5diff_max_differences = 15

where the only field in the .h5 file that is used for comparison in the regression test is `Total_TorqueZ`.
Everything else in the file in not required for the reggie.
As described in the previous section, most reggie analysis functions can be run with the argument `-z, --rc`, which copies a file
resulting from a piclas simulation as reference file directly into the regression test directory to be used in future runs without
the need for copying the file and renaming it by hand

    reggie -iz ../regressioncheck/WEK_DSMC/Torque_Output

which creates the file `TestRotatingWall_DSMCSurfState_000.00100000000000000_reference.h5` in the directory
`/WEK_DSMC/Torque_Output`.
The argument `-i` is only used to create the reference file from a single-core execution as the test later uses more processes, but
the same result is obviously required.
Additionally, unnecessary parts of the .h5 file can be removed if these are not required for the reggie analysis.
This can be achieved by modifying the .h5 file directly with, e.g., [hdfview](https://www.hdfgroup.org/download-hdfview/).
Open the reference file with [hdfview](https://www.hdfgroup.org/download-hdfview/), remove all datasets that are not needed and
resize the dataset `SurfaceData` so that it only contains one single property, name `Total_TorqueZ`. If most of the options,
when right-clicking on the dataset are not available, make sure to `Reload File As` -> `Read/Write` first.
After saving the file to the disk, it is necessary to remove the occupied disk space of the deleted datasets via

    h5repack -i original.h5 -o compressed.h5

and

    h5repack -v -f SHUF -f GZIP=9 original.h5 compressed-shuffled-GZIP9.h5
    h5repack -v -f GZIP=9 original.h5 compressed-GZIP9.h5

by trying both commands and using the smaller resulting .h5 file.
[Note that shuffle by itself will do nothing to compress. It simply re-orders bytes in memory in hopes of making the resulting
byte-stream easier to compress for something like GZIP than it would be able to do otherwise. This is because GZIP is a byte-level
compressor. It doesn’t know about things like shorts or ints or doubles.](https://forum.hdfgroup.org/t/h5repack-gzip-1-slow/4283/5)

## Running GitLab *.gitlab-ci.yml* Tests

The GitLab CI/CD tests can either be run *locally* or *remotely* and both methods are explained in the following.
The tests are defined in the file *.gitlab-ci.yml* in the top-level directory of the piclas repository.

## Remote Testing on Gitlab

Open a browser and go to the [piclas gitlab pipelines website](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines), where the
latest pipeline jobs are displayed. To start a new pipeline, click the button *Run pipeline* and select the required branch name or
tag, which should be tested. Then, activate the necessary *inputs* by selecting `true`, which are summarized in {numref}`tab:pipeline_inputs`.

```{table} Gitlab pipeline inputs
---
name: tab:pipeline_inputs
---
| Property                      | Description                                                                                                                              |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| DO_CHECKIN                    | short tests that are also run, when new commits are pushed                                                                               |
| DO_NIGHTLY                    | longer tests, executed every day                                                                                                         |
| DO_WEEKLY                     | very long tests, executed once a week                                                                                                    |
| DO_CODE_COVERAGE              | Generate coverage data of piclas for all jobs in current pipeline                                                                        |
| DO_REGGIE_COVERAGE            | Generate coverage data of reggie itself for all jobs in current pipeline                                                              |
| DO_NODE_SPLIT                 | MPI: virtual CPU splitting for multi-node testing, where a specific number of cores/threads are grouped in separate nodes (default is 2) |
| DO_CORE_SPLIT                 | MPI: virtual CPU splitting for multi-node testing, where each core/thread resembles a separate node                                      |
| DO_MPICH                      | MPI: Force compilation using MPICH instead of OpenMPI                                                                                    |
```

Per default, `DO_CHECKIN`, `DO_NIGHTLY`, `DO_WEEKLY`, `DO_NODE_SPLIT` and `DO_CORE_SPLIT` are tested automatically for the branch
*master.dev*. For details, see the [Pipeline schedules](https://piclas.boltzplatz.eu/piclas/piclas/-/pipeline_schedules) section in Gitlab.

## Code Coverage

Code coverage information can be inspected in multiple ways: on GitLab, by downloading GitLab artifacts and inspecting the data locally, or by running the tests locally in the first place. By default, all lines containing `Call Abort()` or `Call CollectiveStop()` are not considered when creating coverage reports.

### Coverage on GitLab

To enable code coverage when running regression tests on GitLab, set the pipeline input `DO_CODE_COVERAGE` to `true`. This creates a separate stage, which is executed after all other stages. This way the code coverage stage is able to collect the coverage data from all (previous) jobs/runs of the current pipeline. This is done in the following way.

For each job (e.g., CHE_DSMC), the coverage data is collected and stored in a `.json` report file per build. These reports are named after the build and stored in the `Coverage` directory. If the build is reused for another job (e.g., same DSMC build for both CHE_DSMC and NIG_DSMC), the report file is updated. After all jobs in the other stages are finished, the coverage stage starts. The report files from all builds are combined into a single report containing the coverage data from all previously triggered runs/builds of the current pipeline. For example, if only the inputs `DO_CODE_COVERAGE` and `DO_CHECKIN` are set, the report will show coverage data from all builds with corresponding runs in the `DO_CHECKIN` case. The output of the coverage stage will also display which `.json` report files are found and therefore used to check whether the correct files and builds are considered.

#### Inspecting data on GitLab

The coverage stage will create a [GitLab report artifact](https://docs.gitlab.com/ci/yaml/artifacts_reports/#artifactsreportscoverage_report), which is used for the visualization. GitLab currently uses the coverage report of the latest pipeline of the current branch (if not expired yet) to display the coverage. Note that if the latest pipeline did not create coverage data, nothing will be displayed in the merge request difference view. This is done to prevent a false representation for new code evaluated with older coverage information. Moreover, displaying older coverage data would not be possible if the source code has changed. Keep in mind that new pipelines might change the coverage data if a different set of tests is run. The coverage report is shown in the merge request difference view/changes. A line that was tested is indicated by a green bar to its left, otherwise a red bar appears. This allows inspection of regression tests for new features directly on GitLab. For smaller features/tests, it is recommended to check the coverage locally first to avoid triggering unnecessary tests.

To generate a full coverage report for all available regression tests of PICLas, either:
* Execute all tests in the same pipeline with `DO_CODE_COVERAGE` set to `true`, or
* Combine separate reports manually

On GitLab the coverage is shown as a single number either in the output of the coverage job, on the right side when inspecting the job, or even in the merge request view. The displayed number is the line coverage (different coverage types in "Inspecting data locally"), which is set in `.gitlab-ci.yml`.

Currently the setup does not allow to test node/core split together, since the report can only be created per pipeline. Therefore, the reports have to be combined manually or inspected separately.

#### Inspecting data locally

In addition to the report artifact, the data is stored as a GitLab job artifact, which can be downloaded and inspected locally. This can be done on GitLab in the Pipeline overview on the right under Actions by choosing the Coverage artifact. To get a well-structured view of the data, it is recommended to use the `Coverage/html/coverage.html` file. Total coverage is displayed at the top right and split into three different categories: Lines, Functions and Branches. **Lines** checks how many of the compiled lines were actually executed. Lines that were covered by tests are highlighted in green, while untested lines appear in red. Yellow lines indicate some kind of branching, e.g. an if clause, where one condition was not tested. **Functions** works similarly, but only checks for function calls. So if a file contains only one function and it is called somewhere else, this will result in 100% coverage even if only 10% of lines were executed in this function. Lastly, **Branches** check how many different paths the program could take and then counts how many it did. So this mainly counts if/else statements or case select blocks.

### Local coverage tests

Since reggie supports coverage outputs, it can also be generated locally. To enable code coverage information, each executable must be compiled with additional flags via the `PICLAS_CODE_COVERAGE` option. This generates additional `.gcno` and `.gcda` files per object file, which track all compiled lines and the number of calls per line.
For example, create a build directory and compile the Poisson solver with the Leapfrom time integration scheme

    mkdir build_poisson_code_coverage && cd build_poisson_code_coverage
    cmake .. -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=Leapfrog -DLIBS_USE_PETSC=ON -DPICLAS_CODE_COVERAGE=ON
    make -j

To enable code coverage when using reggie locally, use the `-o` option. More information can be found in the [reggie documentation](https://github.com/reggie-framework/reggie).
Note that due to the module names, e.g., "__mod_dsmc_MOD_dsmc_main", gcovr must be run with `--include-internal-functions`, otherwise all functions will be excluded.
This can be done with the additional reggie flag `--gcovr_extra`, e.g. `--gcovr_extra \'--include-internal-functions\'`.

To run the regression test and create the code coverage data for piclas

    reggie -e bin/piclas ../regressioncheck/CHE_poisson_p_adaption/Laplace_h_N1_pAdaptionType0/ -o 1

and check if the report is created in the output_dir directory

    find . -name *report.html

View the report with any browser

    firefox ./output_dir/Coverage/combined_report/html/combined_report.html

#### Combining single reports

In some cases, it might be helpful to combine single reports of different reggie runs. This can be done using [gcovr](https://github.com/gcovr/gcovr), which is the same tool that reggie uses itself. For this case, `.json` files are used. Separate reports can be combined with
```
gcovr --root <root_dir> --add-tracefile <json_file1> --add-tracefile <json_file2> --html-nested report_name.html
```
where `--add-tracefile` takes wildcard arguments as well. The `--html-nested` option generates the same output format as reggie on GitLab, which is nicely structured analogous to the src directory. `--root` specifies the root directory containing the source files and `report_name.html` the output file.

Note that `<root_dir>` must be the same directory as the one used for the gcovr call that originally created the single report files. E.g. all single reports are created with `--root ~/some_dir`, then the `<root_dir>` for combining the reports must also be `--root ~/some_dir`. Reggie usually tries to append `src` to the found root directory to exclude UnitTests in the coverage information for piclas. Therefore combining reports generated by reggie for piclas would be done with
```
gcovr --root /path/to/piclas/src --add-tracefile Coverage/*.json --html-nested report_name.html
```

For more information on gcovr and coverage report formats, see the [gcovr documentation](https://gcovr.com/).

### Reggie coverage

Besides generating code coverage reports of PICLas, it is also possible to generate a report of the reggie tool itself. To do this, set the pipeline input `DO_REGGIE_COVERAGE` to `true`, which wraps each reggie call with the [Python coverage tool](https://coverage.readthedocs.io/). This generates a coverage report of all used lines in the reggie module, which is stored as a GitLab artifact. The report can be inspected using the `Coverage/reggie/index.html` file.

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
DO_RUN_LOCAL="cmake/3.30.3   gcc/14.2.0   mpich/4.1.2/gcc/14.2.0    hdf5/1.14.0/gcc/14.2.0/mpich/4.1.2    petsc/3.21.6/gcc/14.2.0/mpich/4.1.2"
gitlab-ci-local --input DO_RUN_LOCAL=$DO_RUN_LOCAL
```
If multiple inputs are required add them to the command
```
gitlab-ci-local --input DO_RUN_LOCAL=$DO_RUN_LOCAL --input CHECK_WARNINGS=true
```
to envoke additional options of the pipeline.

### Example
To run a specific reggie job, in this case a *weekly* reggie that depends on another job, the following parameters are passed
```
DO_RUN_LOCAL="cmake/3.30.3   gcc/14.2.0   mpich/4.1.2/gcc/14.2.0    hdf5/1.14.0/gcc/14.2.0/mpich/4.1.2    petsc/3.21.6/gcc/14.2.0/mpich/4.1.2"
gitlab-ci-local --shell-isolation --needs WEK_Radiation --input DO_RUN_LOCAL=$DO_RUN_LOCAL --input DO_WEEKLY=true
```
where the arguments are listed and explained in {numref}`tab:gitlab_ci_local_inputs`
```{table} gitlab-ci-local inputs example
---
name: tab:gitlab_ci_local_inputs
---
  | Parameter                             | Description                                                                         |
  | ------------------------------------- | ----------------------------------------------------------------------------------- |
  | --shell-isolation                     | Avoid errors due to parallel writing of the ctags.txt file                          |
  | --needs WEK_Radiation                 | Run WEK_Radiation, which also requires WEK_DSMC_Radiation                           |
  | --input DO_RUN_LOCAL=$DO_RUN_LOCAL    | Pass the locally installed and used modules                                         |
  | --input DO_WEEKLY=true                | WEK_Radiation is a weekly reggie that requires the DO_WEEKLY to be passed           |
```
