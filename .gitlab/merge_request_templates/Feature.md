## To-Do's

* [ ] ToDo

## Merge Request Checklist

* [ ] Make sure the [Style Guide](https://piclas.readthedocs.io/en/latest/developerguide/styleguide.html) is respected
* Maximum number of 10 compiler warnings
  * [ ] Check with specific compiler settings for the feature branch via *./tools/test_max_warnings.sh*. Number of found warnings:
  * [ ] Run [pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch and supply the variables `DO_CHECKIN=T` and `CHECK_WARNINGS=T` for automatic compiler warning tests for other compiler flag combinations
* [ ] Check file size via *./tools/test_max_file_size.sh*. Largest file size:
* Descriptions for new/changed routines
  * [ ] Short header description: Do not just spell out the name of the subroutine and add units for important variables if applicable
  * [ ] Workflow
    * [ ] Short summary in the header
    * [ ] Inside the routine at the appropriate positions
* Reggie: The new feature must be tested with at least one new or old regression test(s)
  * [ ] Add small test setup if the new feature is not covered by any old regression tests
  * [ ] Add entry in REGGIE.md table by running the reggie table script within the reggie folder where the builds.ini file is via `./../../tools/reggietable.sh` and adjusting the output
  * [ ] Check correct memory allocation and deallocation for the reggie test case
    * [ ] Check automatic restart functionality of reggie example via Load Balance
    * [ ] Compile PICLas with Sanitizer and MPI=OFF as well as MPI=ON and run with one process to find possible memory
          leaks. When using MPICH, the test should also be performed with multiple processes. Leaks can be identified using
          [this approach](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#possible-memory-leak-detection-when-using-mpich).
    * [ ] New memory allocation: How much memory is now allocated? Are new arrays, types, structures allocated when the new model is active or always? Can the memory footprint be improved?
    * [ ] Are there arrays being allocated in the declaration section? See [the problems that can occur](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#seemingly-meaningless-change-in-code-triggers-segmentation-fault-or-slow-down-of-the-code) and an example where this has been fixed in [0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42](https://github.com/piclas-framework/piclas/commit/0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42).
  * Test the three shared memory modes for the reggie by hand
    * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
    * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, .i.e, each process yields a logical node
    * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] New feature description in appropriate documentation (user/developer guide)
* [ ] Are there new/changed shared memory windows? Check if the [rules in the dev guide are being followed](https://piclas.readthedocs.io/en/latest/developerguide/bestpractices.html#shared-memory-windows).