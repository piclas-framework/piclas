## To-Do's

* [ ] ToDo

## Merge Request Checklist

* [ ] Style Guide
* [ ] Maximum of 10 compile warnings via *./tools/test_max_warnings.sh*. How many warning were found?
* [ ] No large files via *./tools/test_max_file_size.sh*. What is the largest file?
* [ ] Descriptions for new/changed routines
  * [ ] Short header description (do not just spell out the name of the subroutine, units for important variables if applicable)
  * [ ] Workflow
    * [ ] Short summary in the header
    * [ ] Inside the routine at the appropriate positions
* [ ] Reggie: The new feature must be tested with at least one new or old regression test(s)
  * [ ] Add small test setup if the new feature is not covered by any old regression tests
  * [ ] Add entry in REGGIE.md table by running the reggie table script within the reggie folder where the builds.ini file is via `./../../tools/reggietable.sh` and adjusting the output
  * [ ] Check correct allocation and deallocation for the test case
    * [ ] Either check automatic restart functionality of reggie example via Load Balance
    * [ ] And/or compile PICLas with Sanitizer and MPI=OFF as well as MPI=ON and run with one process to find possible memory
          leaks. When using MPICH, the test should also be performed with multiple processes. Leaks can be identified using
          [this approach](https://piclas.readthedocs.io/en/latest/developerguide/troubleshooting.html#possible-memory-leak-detection-when-using-mpich).
  * [ ] Test the three shared memory modes
    * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
    * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, .i.e, each process yields a logical node
    * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] New feature description in appropriate documentation (user/developer guide)
* [ ] Replace `MPI_COMM_WORLD` with `MPI_COMM_PICLAS`
