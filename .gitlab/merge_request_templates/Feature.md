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
* [ ] Reggie
  * [ ] Add small test setup
  * [ ] Add entry in REGGIE.md table
  * [ ] Check automatic restart functionality of reggie example via Load Balance (checks correct allocation and deallocation for the test case)
  * [ ] Test the three shared memory modes
    * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` (default) for splitting shared memory domains on the physical node
    * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, .i.e, each process yields a logical node
    * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] New feature description in appropriate documentation (user/developer guide)
* [ ] Replace `MPI_COMM_WORLD` with `MPI_COMM_PICLAS`
