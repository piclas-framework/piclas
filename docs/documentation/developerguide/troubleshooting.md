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
- Setting `collective=.TRUE.` triggers the usage of `H5FD_MPIO_COLLECTIVE_F` (`collective=.FALSE.` uses `H5FD_MPIO_INDEPENDENT_F`) in
  `H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)`, which configures the "transfer mode" in the hdf5 output.
  Collective MPI output requires that all processes take part in the operation!

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


## Possible memory leak detection when using MPICH
**Error Description**
- The error output looks like this

      ====================================================================================================================================
      PICLAS FINISHED! [ 2.41 sec ] [ 0:00:00:02 ]
      ====================================================================================================================================
      Abort(810645775): Fatal error in internal_Finalize: Other MPI error, error stack:
      internal_Finalize(50)...........: MPI_Finalize failed
      MPII_Finalize(400)..............:
      MPID_Finalize(652)..............:
      MPIDI_OFI_mpi_finalize_hook(856):
      destroy_vni_context(1094).......: OFI domain close failed (ofi_init.c:1094:destroy_vni_context:Device or resource busy)
      [WARNING] yaksa: 4 leaked handle pool objects

  and shows that piclas finishes successfully, but an MPI error is invoked afterwards.
  Note that last line containing "[WARNING] yaksa: 4 leaked handle pool objects" might not be there and sometimes reflects the
  number of processes minus one.

**Necessary Conditions**
- MPICH must be used instead of OpenMPI. The problem even occurs when only one single process is used.


**Finding the Bug**
- Activate `PICLAS_DEBUG_MEMORY=ON` and check all the pairs of, e.g.,

      myrank=      0               Allocated ElemBaryNGeo_Shared_Win with WIN_SIZE =                  240

  with

      myrank=      0                          Unlocking ElemBaryNGeo_Shared_Win with MPI_WIN_UNLOCK_ALL()
      myrank=      0                     Freeing window ElemBaryNGeo_Shared_Win with       MPI_WIN_FREE()

  to find the missing `CALL UNLOCK_AND_FREE(ElemBaryNGeo_Shared_Win)` by running piclas and storing the output in, e.g., `std.out`
  and then running the following command

      STDOUT='std.out'; dashes='----------------------------------------'; for line in $(grep -o -P '(?<=Allocated).*(?=with)' ${STDOUT} | sort -u | xargs); do printf 'Checking [%s] %s' "$line" "${dashes:${#line}}"; NbrOfA=$(grep -iw "${line}" ${STDOUT} | grep -ic "Allocated"); printf ' allocated %sx' "$NbrOfA"; NbrOfDA=$(grep -iw "${line}" ${STDOUT} | grep -ic "Unlocking"); printf ' deallocated %sx' "$NbrOfDA"; if [[ $NbrOfDA -lt $NbrOfA ]]; then echo " ---> Could not find MPI_WIN_UNLOCK_ALL() and MPI_WIN_FREE() for this variable"; else echo "... okay"; fi; done

  Replace `std.out` in the command if a different file name is used.
  If a variable is not correctly freed, the output of the script should look like this

      Checking [ElemSideNodeID_Shared_Win] --------------- allocated 2x deallocated 2x... okay
      Checking [ElemMidPoint_Shared_Win] ----------------- allocated 2x deallocated 2x... okay
      Checking [ElemNodeID_Shared_Win] ------------------- allocated 2x deallocated 1x ---> Could not find all required MPI_WIN_UNLOCK_ALL() and MPI_WIN_FREE() for this variable
      Checking [ElemBaryNGeo_Shared_Win] ----------------- allocated 2x deallocated 2x... okay
      Checking [ElemRadius2NGeo_Shared_Win] -------------- allocated 2x deallocated 2x... okay
      Checking [ElemCurved_Shared_Win] ------------------- allocated 2x deallocated 2x... okay

**Resolution**
- Add the missing `CALL UNLOCK_AND_FREE(MYNAME_Win)`, where `MYNAME` is the name of the shared memory window.

**Explanation**
- `MPI_WIN_UNLOCK_ALL()` and `MPI_WIN_FREE()` must be applied to shared memory windows before `CALL MPI_FINALIZE(iError)` is called.

**git hashes**
- One of these bugs was specifically fixed in
  [1a151c24bab7ea22809d3d7554ff5ddf18379cf1](https://github.com/piclas-framework/piclas/commit/1a151c24bab7ea22809d3d7554ff5ddf18379cf1)

### PETSc

To check whether any PETSc objects have not been freed after the `PetscFinalize` call, prepend the following parameters:

    PETSC_OPTIONS=-malloc_debug ./piclas parameter.ini

## Compiler Options

* Release: optimized with -O3 for execution runs
* Debug: debugger options
* Sanitize: GNU sanitizer for further debugging
  
|            Compiler-Flag           |     Options,List     |                                                                                                                                                                                                  What does it do?                                                                                                                                                                                                  |
| ---------------------------------- | -------------------- |                                                                                                                                                 -------------------------------------------------------------------------------------------------------------------                                                                                                                                                |
|          --ffpe-trap=list          |       *invalid*      |                                                                                                                                                                                invalid floating point operation, such as SQRT(-1.0)                                                                                                                                                                                |
|                                    |        *zero*        |                                                                                                                                                                                                  division by zero                                                                                                                                                                                                  |
|                                    |      *overflow*      |                                                                                                                                                                                       overflow in a floating point operation                                                                                                                                                                                       |
|                                    |      *underflow*     |                                                                                                                                                   underflow in a floating point. **DO NOT USE**. Because a small value can occure, such as exp(-766.2). operation                                                                                                                                                  |
|                                    |      *precision*     |                                                                                                                                                                                         loss of precision during operation                                                                                                                                                                                         |
|                                    |                      |                                                                                            Some of the routines in the Fortran runtime library, like **CPU_TIME**, are likely to trigger floating point exceptions when ffpe-trap=precision is used. For this reason, the use of ffpe-trap=precision is not recommended.                                                                                           |
|                                    |      *denormal*      |                                                                                                                                                                                         operation produced a denormal value                                                                                                                                                                                        |
|             -fbacktrace            |                      |                                                                                                                                                                                runtime error should lead to a backtrace of the error                                                                                                                                                                               |
|           -fcheck=keyword          |         *all*        |                                                                                                                                                                                              enable all run-time check                                                                                                                                                                                             |
|                                    |     *array-temps*    |                                                                                                   Warns at run time when for passing an actual argument a temporary array had to be generated. The information generated by this warning is sometimes useful in optimization, in order to avoid such temporaries.                                                                                                  |
|                                    |       *bounds*       |                                         Enable generation of run-time checks for array subscripts and against the declared minimum and maximum values. It also checks array indices for assumed and deferred shape arrays against the actual allocated bounds and ensures that all string lengths are equal for character array constructors without an explicit typespec.                                         |
|                                    |         *do*         |                                                                                                                                                              Enable generation of run-time checks for invalid modification of loop iteration variables                                                                                                                                                             |
|                                    |         *mem*        |                                                                                                                    Enable generation of run-time checks for memory allocation. Note: This option does not affect explicit allocations using theALLOCATE statement, which will be always checked.                                                                                                                   |
|                                    |       *pointer*      |                                                                                                                                                                         Enable generation of run-time checks for pointers and allocatables.                                                                                                                                                                        |
|                                    |      *recursion*     |                                                                          Enable generation of run-time checks for recursively called subroutines and functions which are not marked as recursive. See also -frecursive. Note: This check does not work for OpenMP programs and is disabled if used together with -frecursiveand -fopenmp.                                                                          |
|             -fdump-core            |                      |                                                                                                       Request that a core-dump file is written to disk when a runtime error is encountered on systems that support core dumps. This option is only effective for the compilation of the Fortran main program                                                                                                       |
|           -fstack-arrays           |                      |                                    Adding this option will make the Fortran compiler put all local arrays, even those of unknown size onto stack memory. If your program uses very large local arrays it is possible that you will have to extend your runtime limits for stack memory on some operating systems. This flag is enabled by default at optimization level -Ofast.                                    |
|           -frepack-arrays          |                      | In some circumstances GNU Fortran may pass assumed shape array sections via a descriptor describing a noncontiguous area of memory. This option adds code to the function prologue to repack the data into a contiguous block at runtime.This should result in faster accesses to the array. However it can introduce significant overhead to the function call, especially when the passed data is noncontiguous. |
|       -finline-matmul-limit=n      |                      |                                                                                                                                                                                                                                                                                                                                                                                                                    |
|          -finit-local-zero         |                      |                                                                                                         The -finit-local-zero option instructs the compiler to initialize local INTEGER, REAL, and COMPLEX variables to zero, LOGICALvariables to false, and CHARACTER variables to a string of null bytes                                                                                                         |

## Maximum Number of Warnings

The maximum number of allowed warnings is limited to 10 and can be checked by using the script *test_max_warnings.sh*.
Navigate to the top-level directory of the repository and execute

```
./tools/test_max_warnings.sh
```

This creates a build directory *build_test_max_warnings* and launches *cmake*, which must be
configured by the user (hit *c* for configure). The user must then supply the desired compilation flags and must 
complete configuration by hitting *c* again and then continue the script by generating the make
files (hit *g* for generate). The
number of warnings is then supplied.
