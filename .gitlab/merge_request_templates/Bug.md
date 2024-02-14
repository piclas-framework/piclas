## Related Issue

Closes #number

## Merge Request Checklist

* [ ] Style Guide
* [ ] Maximum of 10 compile warnings via *./tools/test_max_warnings.sh*
* [ ] No large files via *./tools/test_max_file_size.sh*. What is the largest file?
* [ ] Test the three shared memory modes
  * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
  * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, .i.e, each process yields a logical node
  * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] Replace `MPI_COMM_WORLD` with `MPI_COMM_PICLAS`
