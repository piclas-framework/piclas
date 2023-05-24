# Troubleshooting
The following examples are the result of lengthy debugging sessions and can be checked if problems arise during the development
process. Especially when using MPI, the debugging process can be quite cumbersome.

## WriteArrayToHDF5() and the collective flag
**Error**
- One or more array elements in a .h5 file are corrupt, e.g., the element contains the value `1.16828195e-314` instead of the expected
value of `7960`

**Necessary Conditions**
- `LIBS_USE_MPI=ON`: the error only occurs when using more than 1 process (a multi-node run with a large number of processes
  might yield a high chance to trigger this problem)

**Sufficient Conditions**
- `CALL WriteArrayToHDF5()` with `collective=.TRUE.` but not all processes enter this function

**Resolution**
- `CALL WriteArrayToHDF5()` with collective=.FALSE. when it is not 100% certain that all processes enter this routine

**git hashes**
- one of these bugs was specifically fixed in
  [0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42](https://github.com/piclas-framework/piclas/commit/0b2f7b12ecdf84d095caeb8c4b35e08a8484ce42)

**Explanation**
- collective=.TRUE. triggers the usage of `H5FD_MPIO_COLLECTIVE_F` and `collective=.FALSE.` uses `H5FD_MPIO_INDEPENDENT_F` in
  `H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)`
