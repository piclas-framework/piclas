# Troubleshooting
The following examples are the result of lengthy debugging sessions and can be checked if problems arise during the development
process. Especially when using MPI, the debugging process can be quite cumbersome.


## WriteArrayToHDF5() and the collective flag
**Error Description**
- One or more array elements in a .h5 file are corrupt, e.g., the element contains the value `1.16828195e-314` instead of the expected
value of `7960`

**Necessary Conditions**
- `LIBS_USE_MPI=ON`: the error only occurs when using more than 1 process (a multi-node run with a large number of processes
  might yield a high chance to trigger this problem)

**Sufficient Conditions**
- `CALL WriteArrayToHDF5()` with `collective=.TRUE.` but not all processes enter this function

**Finding the Bug**
- This error can be found with a regression test that runs with `LIBS_USE_MPI=OFF` or one process and again with multiple processes
at best using the multi-node feature `PICLAS_SHARED_MEMORY=OMPI_COMM_TYPE_CORE`

**Resolution**
- `CALL WriteArrayToHDF5()` with collective=.FALSE. when it is not 100% certain that all processes enter this routine

**Explanation**
- Setting collective=.TRUE. triggers the usage of `H5FD_MPIO_COLLECTIVE_F` and `collective=.FALSE.` uses `H5FD_MPIO_INDEPENDENT_F` in
  `H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)`, which configures the "transfer mode" in the hdf5 output

**git hashes**
- One of these bugs was specifically fixed in
  [0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42](https://github.com/piclas-framework/piclas/commit/0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42)


## Seemingly meaningless change in code triggers segmentation fault or slow down of the code
**Error Description**
- A user-defined read-in parameter is commented in/out which triggers the bug due to different values that are read into the
  parameter, e.g., setting `1e-14` or `10e-15` (which is basically the same value!)
  The suspected problem is that something is read from the memory, which should not be read, e.g., an uninitialized variable
  (that therefore points to a possibly random location in the memory)
  This causes an undefined state (non-deterministic outcome, different compilers/different machines yield different effects).
  The code might not crash but hang for a certain amount of time (this can be used to find the bug).

**Necessary Conditions**
- An integer variable that is used for (indirect) allocation of an array is uninitialized

**Sufficient Conditions**
- A function or subroutine is called that declares an array depending on an uninitialized integer variable

**Finding the Bug**
- When the code does not crash but instead hangs for a certain amount of time, the following print statements can be used (when not
  utilizing a debugger)

      IPWRITE(UNIT_StdOut,'(I0,A,I0)') ": v "//TRIM(__FILE__)//" +",__LINE__

  to see the exact position in the code where the code hangs (due to a possible function call and subsequent allocation process)
  for a short period of time

**Resolution**
- Nullify all integer variables that are used for allocation per default

**Explanation**
- The following variable was declared in a subroutine

      REAL :: DistriOld(SpecDSMC(PartSpecies(iPart1))%MaxElecQuant)

  but the integer `SpecDSMC(PartSpecies(iPart1))%MaxElecQuant` was not initialized because the corresponding model is not used in
  this specific case and therefore `DistriOld` is never used. It is however allocated with an undefined state, with undefined outcome!
  In this case the bug was fixed by using the "assignment to an allocatable array"

      REAL,ALLOCATABLE              :: DistriOld(:)
      ...
      DistriOld = ElectronicDistriPart(iPart1)%DistriFunc

  For gfortran, Fortran 2003's assignment to an allocatable array was introduced in 4.6, 2011-01-28.

**git hashes**
- One of these bugs was specifically fixed in
  [8d1129be95abc91bf56e94bf7f12987e4c214666](https://github.com/piclas-framework/piclas/commit/8d1129be95abc91bf56e94bf7f12987e4c214666)
