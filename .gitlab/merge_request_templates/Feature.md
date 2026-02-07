## To-Do's

* [ ] ToDo

## Merge Request Checklist

* [ ] Make sure the [Style Guide](https://piclas.readthedocs.io/en/latest/developerguide/styleguide.html) is respected
* [ ] Make sure the [Best Practices](https://piclas.readthedocs.io/en/latest/developerguide/styleguide.html) guide is followed
* Maximum number of 10 compiler warnings
  * [ ] Check with specific compiler settings for the feature branch via `./tools/test_max_warnings.sh`. Number of found warnings:
  * [ ] Run [pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch and supply the variables `DO_CHECKIN=T` and `CHECK_WARNINGS=T` for automatic compiler warning tests for other compiler flag combinations
* [ ] Check file size via *./tools/test_max_file_size.sh*. Write the name and file size of the largest here: _________
* [ ] Check if newly introduced `CALL abort(...)` statements can be replaced with `CALL CollectiveStop(...)`, which can mostly be achieved during initialisation.
  For details on using this function, see the [Developer Guide: CollectiveStop](https://piclas.readthedocs.io/en/latest/developerguide/bestpractices.html#collectivestop) section.
* Descriptions for new/changed routines
  * [ ] Short header title: Do not just spell out the name of the subroutine! Add units for important variables if applicable.
  * [ ] Workflow
    * [ ] Short [header summary](https://github.com/piclas-framework/piclas/blob/790daf835fd76e24a2ca8eb2e3021e149c5f5c09/src/globals/globals.f90#L417)
    * [ ] [Inside the routine](https://github.com/piclas-framework/piclas/blob/790daf835fd76e24a2ca8eb2e3021e149c5f5c09/src/dg/fillflux.f90#L85) at the appropriate positions
* Reggie: The new feature must be tested with at least one new or old [regression test(s)](https://github.com/piclas-framework/piclas/tree/master/regressioncheck)
  * [ ] Add small test setup if the new feature is not covered by any old regression tests
  * [ ] Add entry in [REGGIE.md table](https://github.com/piclas-framework/piclas/blob/master/REGGIE.md) by running the reggie table script within the corresponding reggie folder where the builds.ini file is via `./../../tools/reggietable.sh` and adjusting the output
  * Check correct memory allocation and deallocation for the reggie test case
    * [ ] Check automatic restart functionality of reggie example via load balancing
    * [ ] Compile PICLas with Sanitizer and `MPI=OFF` as well as `MPI=ON` and run with one process to find possible memory leaks.
          When using MPICH, the test should also be performed with multiple processes. Leaks can be identified using
          [this approach](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#possible-memory-leak-detection-when-using-mpich).
    * [ ] New memory allocation: How much memory is now allocated? Are new arrays, types, structures allocated when the new model is active or even when the model is deactivated? Can the memory footprint be improved by only allocating arrays when the model is active?
    * [ ] Are there arrays being allocated in the declaration section? See [the problems that can occur](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#seemingly-meaningless-change-in-code-triggers-segmentation-fault-or-slow-down-of-the-code) and an example where this has been fixed in [this commit](https://github.com/piclas-framework/piclas/commit/0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42).
  * Test the three [shared memory modes](https://piclas.readthedocs.io/en/latest/userguide/workflow.html#compiler-options) for the reggie by hand
    * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
    * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, for example, each process yields a logical node
    * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
  * [ ] When all the above points regarding the reggie have been completed, [run a Gitlab pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch with the
    variables `DO_NIGHTLY=T` and `DO_CORE_SPLIT=T`, which are described in the [Developer Guide: Remote Testing on Gitlab](https://piclas.readthedocs.io/en/latest/developerguide/reggie.html#remote-testing-on-gitlab) section.
* [ ] New feature description in appropriate documentation in the [User and/or Developer Guide](https://piclas.readthedocs.io/en/latest/index.html)
* [ ] Are there new or changed shared memory windows (SHM)? Check if the [rules in the Developer Guide are being followed](https://piclas.readthedocs.io/en/latest/developerguide/bestpractices.html#shared-memory-windows).
* Add line to **Release Notes** for the next **Release X.X.X** under [Merge requests](https://piclas.boltzplatz.eu/piclas/piclas/-/merge_requests) with a short description of this MR. Note that the commit hash is created upon pressing the merge button and is displayed in the **Activity** section afterwards.
* [ ] Run the regression checks, which should test the new feature (either new tests or existing tests using the added lines) with code coverage (`DO_CODE_COVERAGE=T` or locally to avoid unnecessary runs) and check that
  * [ ] all new features are tested (visible as green/red bars next to each code line in merge request diff view)