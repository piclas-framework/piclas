# Best Practices

The following collection of best practice guidelines are intended to prevent bugs and improve the computational performance.

## MPI

The general rules can be summarized as follows:

1. **The first rule of MPI is**: You do not send subsets of arrays, only complete continuous data ranges.
2. **The second rule of MPI is**: You do not send subsets of arrays, only complete continuous data ranges.
3. **Third rule of MPI**: Someone sends non-continuous data, the simulation is over.
4. **Fourth rule**: Only two processors to a single send-receive message. One sender and one receiver, no more and no less.
5. **Fifth rule**: Only one processor access (read or write) to a shared memory region.
6. **Sixth rule**: After nullification of a shared array follows WIN_SYNC and BARRIER because until the sync is complete the status of the memory is undefined, i.e., old or new value or utter nonsense.

Please also read the general implementation information and, e.g., mappings used for elements, sides and nodes in the chapter
{ref}`developerguide/mpi:MPI Implementation`.

If you want to break the first/second rule, remember that in FORTRAN the last dimension can be used for slicing.

## Shared Memory Windows

The following principles should always be considered when using shared memory windows

- Only the node root process initializes the shared memory array

      ! Allocate the shared memory window
      CALL Allocate_Shared((/nUniqueGlobalNodes/), NodeVolume_Shared_Win, NodeVolume_Shared)

      ! Lock the window
      CALL MPI_WIN_LOCK_ALL(0, NodeVolume_Shared_Win, IERROR)

      ! Set pointer
      NodeVolume => NodeVolume_Shared

      ! Only CN root nullifies
      IF (myComputeNodeRank.EQ.0) NodeVolume = 0.0

      ! This sync/barrier is required as it cannot be guaranteed that the zeros have been
      ! written to memory by the time the MPI_REDUCE is executed (see MPI specification).
      ! Until the Sync is complete, the status is undefined, i.e., old or new value or utter
      ! nonsense.
      CALL BARRIER_AND_SYNC(NodeVolume_Shared_Win, MPI_COMM_SHARED)

- When all processes on a node write exclusively to their separate region in the shared memory array, using designated elements IDs which are
  assigned to a single process only

      ! Get offset
      ! J_N is only built for local DG elements. Therefore, array is only filled for elements on the same compute node
      offsetElemCNProc = offsetElem - offsetComputeNodeElem

      ! Allocate shared array
      CALL Allocate_Shared((/nComputeNodeElems/),ElemVolume_Shared_Win,ElemVolume_Shared)
      ...

      ! Calculate element volumes
      DO iElem = 1,nElems
        CNElemID=iElem+offsetElemCNProc
        !--- Calculate and save volume of element iElem
        J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
          ElemVolume_Shared(CNElemID) = ElemVolume_Shared(CNElemID) + wGP(i)*wGP(j)*wGP(k)*J_N(1,i,j,k)
        END DO; END DO; END DO
      END DO

- When all processes on a node write to all regions in the shared memory array, an additional local array is required, which has to be reduced to the shared array at the end

      CALL Allocate_Shared((/nSpecies,4,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallImpactEnergy_Shared_Win,SampWallImpactEnergy_Shared)
      CALL MPI_WIN_LOCK_ALL(0,SampWallImpactEnergy_Shared_Win,IERROR)
      IF (myComputeNodeRank.EQ.0) SampWallImpactEnergy_Shared = 0.
      CALL BARRIER_AND_SYNC(SampWallImpactEnergy_Shared_Win,MPI_COMM_SHARED)
      ALLOCATE(SampWallImpactEnergy(1:nSpecies,1:4,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
      SampWallImpactEnergy = 0.
      SampWallImpactEnergy(SpecID,1,SubP,SubQ,SurfSideID) = SampWallImpactEnergy(SpecID,1,SubP,SubQ,SurfSideID) + ETrans * MPF
      CALL MPI_REDUCE(SampWallImpactEnergy,SampWallImpactEnergy_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
      CALL BARRIER_AND_SYNC(SampWallImpactEnergy_Shared_Win,MPI_COMM_SHARED)

- When possible, never read from the shared memory array in a round robin manner, as shown in this [commit [eaff78c]](https://github.com/piclas-framework/piclas/commit/eaff78c158884e0bab05c555bf72b4ff6198e42f).
  Split the work and then use `MPI_REDUCE` or `MPI_ALLREDUCE`.
  Instead of

      CNVolume = SUM(ElemVolume_Shared(:))

  where all processes traverse over the same memory addresses, which slows down the computation, use

      offsetElemCNProc = offsetElem - offsetComputeNodeElem
      CNVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,CNVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SHARED,iError)

  to split the operations and use MPI to distribute the information among the processes.

- Atomic MPI operations on shared memory
  - Example 1: [Store the min/max extent when building the CN FIBGM [6350cc2]](https://github.com/piclas-framework/piclas/commit/6350cc2575d15c7ceb804bc8d839ca5ef2b33dbb?diff=split#diff-aa2cf11ef2c11ce88cdefcf02fe06b643771c968021311ea356c428bbb20d041L1214)
  - Example 2: [Use atomic MPI operations to read/write from contested shared memory [772c371]](https://github.com/piclas-framework/piclas/commit/772c3711bbb0c935659b2d08fccd18c80e6b72dc)
  - The main idea is to access and change parts of a shared array with multiple processes to, e.g., sum up numbers from different
    processes and guarantee that in the end the sum is correct without having a predefined order in which the numbers are added to the
    entry in the shared array.

    In the example in [772c371], get the memory window while bypassing local caches

        CALL MPI_FETCH_AND_OP(ElemDone,ElemDone,MPI_INTEGER,0,INT(posElem*SIZE_INT,MPI_ADDRESS_KIND),MPI_NO_OP,ElemInfo_Shared_Win,iError)

    Flush only performs the pending operations (getting the value)

        CALL MPI_WIN_FLUSH(0,ElemInfo_Shared_Win,iError)

    Using `MPI_REPLACE` makes sure that the correct value is written in the end by one of the processes in an undefined order.

        MPI_FETCH_AND_OP(haloChange,dummyInt,MPI_INTEGER,0,INT(posElem*SIZE_INT,MPI_ADDRESS_KIND),MPI_REPLACE,ElemInfo_Shared_Win,iError)
        CALL MPI_WIN_FLUSH(0,ElemInfo_Shared_Win,iError)

## Hawk

Before running a simulation, check out the HLRS Wiki pages [Batch System PBSPro (Hawk)](https://kb.hlrs.de/platforms/index.php/Batch_System_PBSPro_(Hawk)).

### Striping
Always use user-defined striping in the simulation case folders that are on the work spaces as the default stiping setting (dynamic
striping) has caused massive problems in the past. Add the following code to your submit script

    # Set fixed striping to avoid problems with the progressive Lustre file layout
    # - Region 1 [0, 1GiB): Stripe-Size=1 MiB, Stripe-Count=1
    #lfs setstripe -c 1 -S 1M $PBS_O_WORKDIR
    # - Region 2 [1GiB, 4GiB): Stripe-Size=1 MiB, Stripe-Count=4
    #lfs setstripe -c 4 -S 1M $PBS_O_WORKDIR
    # - Region 3 [4 GiB, EOF): Stripe-Size=4 MiB, Stripe-Count=8
    lfs setstripe -c 8 -S 4M $PBS_O_WORKDIR

Note that the correct line should be commented in and the other lines should be commented out, all depending on the size of your
output files.
Also consider the stripe settings for large mesh files just to be sure.

### Species-zero bug
It has repeatedly occurred that particles with species index zero have been produced on hawk.
This might be due to the output to .h5, which could reflect the previous section regarding the striping settings, but could also lie
deeper the Lustre file system itself.
If this problem occurs, the corrupted particles must be removed from the .h5 file by hand if a restart from such a corrupted file is
performed in order to prevent piclas from crashing.

## CollectiveStop
When using the `CALL abort(__STAMP__,'ERROR ...')` subroutine, each MPI process that encounters this
call emits an output to std.out, which can result in a huge amount of output when running on hundreds
or thousands of processes, especially if an error in the .ini files produces the error in the initialization step.
To prevent excessive output to std.out, the `CollectiveStop` subroutine has been created.
```
!==================================================================================================================================
!> \brief Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!>
!> Safely terminate program using a soft MPI_FINALIZE in the MPI case and writes the error message only on the root.
!> Terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> This routine can only be used if ALL processes are guaranteed to generate the same error at the same time!
!> Prime use is to exit FLEXI without MPI errors and with a single error message if some parameters are not set in the init
!> routines or a file is not found.
!>
!> Criteria where CollectiveStop may be used:
!> 0. In case of doubt stick with Abort, which is always safe!
!> 1. A routine is BY DESIGN (!) called by all processes, i.e. does not permit to be called by single processes or subgroups.
!> 2. The criteria for the CollectiveStop must be identical among all processors.
!> 3. The routine is only used during the init phase.
!> 4. The error must not originate from MPI errors (e.g. during MPI init)
!> 5. The error must not originate from checking roundof errors (e.g. accuracy of interpolation matrices)
!>
!==================================================================================================================================
```
An example where `CollectiveStop` should be used instead of `abort` is in `gradients.f90`

```
GradLimiterType=GETINT('Grad-LimiterType')
GradLimVktK=GETREAL('Grad-VktK')
SELECT CASE(GradLimiterType)
CASE(0)
  LBWRITE(UNIT_stdOut,*)'Limiter = 0 -> first order FV'
CASE(1) !minmax
  LBWRITE(UNIT_stdOut,*)'Using Barth-Jespersen Limiter'
CASE(4) !venkatakrishnan
  LBWRITE(UNIT_stdOut,*)'Using Venkatakrishnan limiter with K =', GradLimVktK
CASE(9) ! no limiter (central)
  LBWRITE(UNIT_stdOut,*)'Not using any limiter'
CASE DEFAULT
  CALL abort(__STAMP__,'Limiter type not implemented.')
END SELECT
```
which should read
```
GradLimiterType=GETINT('Grad-LimiterType')
GradLimVktK=GETREAL('Grad-VktK')
SELECT CASE(GradLimiterType)
CASE(0)
  LBWRITE(UNIT_stdOut,*)'Limiter = 0 -> first order FV'
CASE(1) !minmax
  LBWRITE(UNIT_stdOut,*)'Using Barth-Jespersen Limiter'
CASE(4) !venkatakrishnan
  LBWRITE(UNIT_stdOut,*)'Using Venkatakrishnan limiter with K =', GradLimVktK
CASE(9) ! no limiter (central)
  LBWRITE(UNIT_stdOut,*)'Not using any limiter'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Limiter type not implemented.')
END SELECT
```
