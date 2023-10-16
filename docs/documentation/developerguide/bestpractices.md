# Best Practices

The following collection of best practice guidelines are intended to prevent bugs and improve the computational performance.

## MPI

The general rules can be summarized as follows:

> The first rule of MPI is: You do not send subsets of arrays, only complete continuous data ranges.
> The second rule of MPI is: You do not send subsets of arrays, only complete continuous data ranges.
> Third rule of MPI: Someone sends non-continuous data, the simulation is over.
> Fourth rule: Only two procs to a single send-receive message.
> Fifth rule: Only one proc access (read or write) toa shared memory region.

## Shared Memory Windows

The following principals should always be considered when using shared memory windows

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

- When all processes on a node write to their separate region in the shared memory array, e.g., designated elements IDs which are
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

- When all processes on a node write to all regions in the shared memory array, an additional local array is required

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

- [Use atomic MPI operations to read/write from contested shared memory [772c371]](https://github.com/piclas-framework/piclas/commit/772c3711bbb0c935659b2d08fccd18c80e6b72dc)

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
