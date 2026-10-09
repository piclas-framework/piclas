## Related Issue

Closes #number

## Merge Request Checklist

* [ ] Make sure the [Style Guide](https://piclas.readthedocs.io/en/latest/developerguide/styleguide.html) is respected and the [MPI rules](https://piclas.readthedocs.io/en/latest/developerguide/mpi.html#general-rules) are followed
  * [ ] Check if newly introduced `CALL abort(...)` statements can be replaced with `CALL CollectiveStop(...)`, which can mostly be achieved during initialisation.
    For details on using this function, see the [Developer Guide: CollectiveStop](https://piclas.readthedocs.io/en/latest/developerguide/code_extension.html#collectivestop) section.
  * [ ] Are there new or changed shared memory windows (SHM)? Check if the [rules in the Developer Guide are being followed](https://piclas.readthedocs.io/en/latest/developerguide/mpi.html#shared-memory-windows).
* Maximum number of 10 compiler warnings
  * [ ] Check with specific compiler settings for the feature branch via `./tools/test_max_warnings.sh`. Number of found warnings:
  * [ ] Run [pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) for the feature branch and set the inputs `DO_CHECKIN` and `CHECK_WARNINGS` to `true` for automatic compiler warning tests for other compiler flag combinations
* [ ] Check file size via *./tools/test_max_file_size.sh*. Write the name and file size of the largest here: _________
* [ ] Test the three shared memory modes
  * [ ] `PICLAS_SHARED_MEMORY = MPI_COMM_TYPE_SHARED` (default) for splitting shared memory domains on the physical node
  * [ ] `PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE` for splitting at process level, .i.e, each process yields a logical node
  * [ ] `PICLAS_SHARED_MEMORY = PICLAS_COMM_TYPE_NODE` for splitting at 2 processes per logical node
* [ ] Replace `MPI_COMM_WORLD` with `MPI_COMM_PICLAS`
* [ ] Make sure to label the merge request accordingly (Bug / Improvement) and that the merge request title is appropriate and concise, since it will be automatically utilized for the release notes
