## To-Do's

* [ ] ToDo

## Merge Request Checklist

### Code

* [ ] Make sure the [Style Guide](https://piclas.readthedocs.io/en/latest/developerguide/styleguide.html) is respected and the [MPI rules](https://piclas.readthedocs.io/en/latest/developerguide/mpi.html#general-rules) are followed
  * [ ] Check if newly introduced `CALL abort(...)` statements can be replaced with `CALL CollectiveStop(...)`, which can mostly be achieved during initialisation.
    For details on using this function, see the [Developer Guide: CollectiveStop](https://piclas.readthedocs.io/en/latest/developerguide/code_extension.html#collectivestop) section.
  * [ ] Are there new or changed shared memory windows (SHM)? Check if the [rules in the Developer Guide are being followed](https://piclas.readthedocs.io/en/latest/developerguide/mpi.html#shared-memory-windows).
* Maximum number of 10 compiler warnings
  * [ ] Check with specific compiler settings for the feature branch via `./tools/test_max_warnings.sh`. Number of found warnings:
  * [ ] Run [pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch and set the inputs `DO_CHECKIN` and `CHECK_WARNINGS` to `true` for automatic compiler warning tests for other compiler flag combinations
* [ ] Check file size via *./tools/test_max_file_size.sh*. Write the name and file size of the largest here: _________

### Regression testing

The new feature must be tested with at least one new or old [regression test(s)](https://github.com/piclas-framework/piclas/tree/master/regressioncheck):

* [ ] Add small test setup if the new feature is not covered by any old regression tests
* [ ] Add entry in [REGGIE.md table](https://github.com/piclas-framework/piclas/blob/master/REGGIE.md) by running the reggie table script within the corresponding reggie folder where the builds.ini file is via `cd ~/piclas/regressioncheck/MY_REGGIE_EXAMPLE && ./../../tools/reggietable.sh` and adjusting the output
* Check correct memory allocation and deallocation for the reggie test case
  * [ ] Implement automatic restart functionality check within the reggie example via load balancing. Add a suitable [load balancing
  block](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/NIG_PIC_poisson_Leapfrog/MCC_EBeam_SpeciesSpecificTimestep/parameter.ini) in the `parameter.ini` file.
  * [ ] Implement restart functionality check within the reggie example via
    [parameter-pre.ini](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/CHE_drift_diffusion_explicit-FV/periodic_box_free-stream/parameter-pre.ini) file, which in executed in
    [externals.ini](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/CHE_drift_diffusion_explicit-FV/periodic_box_free-stream/externals.ini)
    and creates a restart state file that is used in [command_line.ini](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/CHE_drift_diffusion_explicit-FV/periodic_box_free-stream/command_line.ini),
    see [example](https://github.com/piclas-framework/piclas/blob/master/regressioncheck/CHE_drift_diffusion_explicit-FV/periodic_box_free-stream/).
  * [ ] Compile PICLas with Sanitizer and `MPI=OFF` as well as `MPI=ON` and run with one process to find possible memory leaks.
        When using MPICH, the test should also be performed with multiple processes. Leaks can be identified using
        [this approach](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#possible-memory-leak-detection-when-using-mpich).
  * [ ] New memory allocation: How much memory is now allocated? Are new arrays, types, structures allocated when the new model is active or even when the model is deactivated? Can the memory footprint be improved by only allocating arrays when the model is active?
  * [ ] Are there arrays being allocated in the declaration section? See [the problems that can occur](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#seemingly-meaningless-change-in-code-triggers-segmentation-fault-or-slow-down-of-the-code) and an example where this has been fixed in [this commit](https://github.com/piclas-framework/piclas/commit/0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42).
* Test the three [shared memory modes](https://piclas.readthedocs.io/en/latest/userguide/workflow.html#compiler-options) for the reggie by hand
  * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
  * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, for example, each process yields a logical node
  * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] Check that all reggies under `regressioncheck` are compatible with `pyhope` by running the script `convertHoprToPyHopeIni.sh`
   within the `regressioncheck` directory `cd regressioncheck && ../tools/convertHoprToPyHopeIni.sh` and adding the files changed by
   the script to the MR
* [ ] When all the above points regarding the reggie have been completed, [run a Gitlab pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch with the
  inputs `DO_NIGHTLY` and `DO_CORE_SPLIT` set to `true`, which are described in the [Developer Guide: Remote Testing on Gitlab](https://piclas.readthedocs.io/en/latest/developerguide/reggie.html#remote-testing-on-gitlab) section.
* [ ] Run the regression checks, which should test the new feature (either new tests or existing tests using the added lines) with code coverage (`DO_CODE_COVERAGE` set to `true` or locally to avoid unnecessary runs) and check that
  * [ ] all new features are tested (visible as green/red bars next to each code line in merge request diff view)

### Documentation

* Descriptions for new/changed routines
  * [ ] Short header title: Do not just spell out the name of the subroutine! Add units for important variables if applicable.
  * [ ] Workflow
    * [ ] Short [header summary](https://github.com/piclas-framework/piclas/blob/790daf835fd76e24a2ca8eb2e3021e149c5f5c09/src/globals/globals.f90#L417)
    * [ ] [Inside the routine](https://github.com/piclas-framework/piclas/blob/790daf835fd76e24a2ca8eb2e3021e149c5f5c09/src/dg/fillflux.f90#L85) at the appropriate positions
* [ ] New feature description in appropriate documentation in the [User and/or Developer Guide](https://piclas.readthedocs.io/en/latest/index.html)
* [ ] Make sure to label the merge request accordingly (Improvement / Feature) and that the merge request title is appropriate and concise, since it will be automatically utilized for the release notes