# MPI Implementation

This chapter describes how PICLas subroutines and functions are parallelized.

## General Rules

The general rules can be summarized as follows:

1. **The first rule of MPI is**: You do not send subsets of arrays, only complete continuous data ranges.
2. **The second rule of MPI is**: You do not send subsets of arrays, only complete continuous data ranges.
3. **Third rule of MPI**: Someone sends non-continuous data, the simulation is over.
4. **Fourth rule**: Only two processors to a single send-receive message. One sender and one receiver, no more and no less.
5. **Fifth rule**: Only one processor access (read or write) to a shared memory region.

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

## Debugging MPI

Debug MPI

    mpirun --mca mpi_abort_print_stack 1 -np 3 ~/piclas/build/bin/piclas parameter.ini DSMC.ini

## Construction of Halo Region (MPI 3.0 Shared Memory)
The general idea is to have the geometry information in a shared memory array on each node. This
reduces the memory overhead and removes the need for halo regions between processors on the same
node. An exemplary unstructured mesh is given in {numref}`fig:dev_mpi_shared_mesh`.
The workflow is as follows:

1. Load the complete mesh into a shared memory window that is solely allocated by the compute-node root processor of size
   `1:nTotalElems` (only node coords, not derived properties, such as metric terms).
1. Load the elements that are assigned to processors on a single compute-node into an array of size `1:nSharedElems`
   (red elements in {numref}`fig:dev_mpi_shared_mesh`).
1. Find the compute-node halo elements and store in an array of size `1:nSharedHaloElems` (blue elements in
   {numref}`fig:dev_mpi_shared_mesh`).
   1. Background mesh (BGM) reduction
   2. Shared sides reduction
   3. Halo communicator

```{figure} figures/mpi_shared_mesh/dev_mpi_shared_mesh.png
---
name: fig:dev_mpi_shared_mesh
width: 200px
align: center
---

Node local (red), halo (blue) and complete mesh (white). The red and blue geometry information is
allocated once on each node and available to all processors on that node
```

### Mesh Geometry

The complete mesh geometry is read-in by each compute-node and available through the ElemInfo_Shared, SideInfo_Shared, NodeInfo_Shared, NodeCoords_Shared arrays. Other derived variables and mappings (e.g. in `MOD_Particle_Mesh_Vars`) are usually only available on the compute-node.

#### ElemInfo

Element information is read-in and built in `ReadMeshElems`

    ElemInfo_Shared(1:ELEMINFOSIZE,1:nGlobalElems)
    ! Possible preprocessor variables for first entry
    ! ElemInfo from mesh file
    ELEM_TYPE         1    ! PyHOPE classification of element (currently not used)
    ELEM_ZONE         2    ! Zone as defined in PyHOPE (currently not used)
    ELEM_FIRSTSIDEIND 3    ! Index of the first side belonging to the element
    ELEM_LASTSIDEIND  4    ! Index of the last side belonging to the element
    ELEM_FIRSTNODEIND 5    ! Index of the first node belonging to the element
    ELEM_LASTNODEIND  6    ! Index of the last node belonging to the element
    ! ElemInfo for shared array (only with MPI)
    ELEM_RANK         7    ! Processor rank to which the element is assigned
    ELEM_HALOFLAG     8    ! Element type:  1 = Regular element on compute node
                           !                2 = Halo element (non-periodic)
                           !                3 = Halo element (periodic)

#### SideInfo

Side information is read-in and built in `ReadMeshSides`

    SideInfo_Shared(1:SIDEINFOSIZE+1,1:nNonUniqueGlobalSides)
    ! Possible preprocessor variables for first entry
    SIDE_TYPE         1    ! PyHOPE classification of side
    SIDE_ID           2    ! tbd
    SIDE_NBELEMID     3    ! Neighbouring element ID
    SIDE_FLIP         4    ! tbd
    SIDE_BCID         5    ! Index of the boundary as given parameter.ini
    SIDE_ELEMID       6    ! Element ID
    SIDE_LOCALID      7    ! Local side ID
    SIDE_NBSIDEID     8    ! Neighbouring side ID
    SIDE_NBELEMTYPE   9    ! Neighbouring element type, 0 = No connection to another element
                           !                            1 = Neighbour element is on compute-node
                           !                            2 = Neighbour element is outside of compute node

#### NodeInfo

Node information is read-in and built in `ReadMeshNodes`

    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
    ! Non-unique node coordinates per element
    Coords(1:3) = NodeCoords_Shared(3,1:8*nGlobalElems)

### Element/Side Mappings

To get the required element/side ID, utilize one of these functions from the module `MOD_Mesh_Tools`

    ! Global element
    GlobalElemID = GetGlobalElemID(CNElemID)
    ! Compute-node element
    CNElemID = GetCNElemID(GlobalElemID)
    ! Global side
    GlobalSideID = GetGlobalSideID(CNSideID)
    ! Compute-node side
    CNSideID = GetCNSideID(GlobalSideID)

If a loop is over the local number of elements, the compute-node and global element can be determined with the respective offsets.

    DO iElem = 1, nElems
        CNElemID = iElem + offsetComputeNodeElem
        GlobalElemID = iElem + offsetElem
    END DO

#### SurfSide

Additionally to conventional sides, mappings for the sides that belong to a boundary condition are available.

    DO iSurfSide = 1,nComputeNodeSurfSides
      GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    END DO
### Particle Element Mappings

    PEM%GlobalElemID(iPart)         ! Global element ID
    PEM%CNElemID(iPart)             ! Compute-node local element ID (GlobalElem2CNTotalElem(PEM%GlobalElemID(iPart)))
    PEM%LocalElemID(iPart)          ! Core local element ID (PEM%GlobalElemID(iPart) - offsetElem)

## Custom communicators

To limit the number of communicating processors, feature specific communicators can be built. In the following, an example is given
for a communicator, which only contains processors with a local surface side (part of the `InitParticleBoundarySurfSides` routine). First, a global variable `SurfCOMM`, which is based on the `tMPIGROUP` type, is required:
```
TYPE tMPIGROUP
  INTEGER         :: UNICATOR=MPI_COMM_NULL     !< MPI communicator for surface sides
  INTEGER         :: nProcs                     !< number of MPI processes
  INTEGER         :: MyRank                     !< MyRank, ordered from 0 to nProcs - 1
END TYPE
TYPE (tMPIGROUP)    :: SurfCOMM
```
To create a subset of processors, a condition is required, which is defined by the `color` variable:
```
color = MERGE(1337, MPI_UNDEFINED, nSurfSidesProc.GT.0)
```
Here, every processor with the same `color` will be part of the same communicator. The condition `nSurfSidesProc.GT.0` in this case includes every processor with a surface side. Every other processor will be set to `MPI_UNDEFINED` and consequently be part of `MPI_COMM_NULL`. Now, the communicator itself can be created:
```
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, SurfCOMM%UNICATOR, iError)
```
`MPI_COMM_PICLAS` denotes the global PICLas communicator containing every processor (but can also be a previously created subset) and the `MPI_INFO_NULL` entry denotes the rank assignment within the new communicator (default: numbering from 0 to nProcs - 1). Additional information can be stored within the created variable:
```
IF(SurfCOMM%UNICATOR.NE.MPI_COMM_NULL) THEN
  ! Stores the rank within the given communicator as MyRank
  CALL MPI_COMM_RANK(SurfCOMM%UNICATOR, SurfCOMM%MyRank, iError)
  ! Stores the total number of processors of the given communicator as nProcs
  CALL MPI_COMM_SIZE(SurfCOMM%UNICATOR, SurfCOMM%nProcs, iError)
END IF
```
Through the IF clause, only processors that are part of the communicator can be addressed. And finally, it is important to free the communicator during the finalization routine:
```
IF(SurfCOMM%UNICATOR.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(SurfCOMM%UNICATOR,iERROR)
```
This works for communicators, which have been initialized with MPI_COMM_NULL, either initially during the variable definition or during the split call.
If not initialized initially, you have to make sure that the freeing call is only performed, if the respective split routine has been called to guarantee
that either a communicator exists and/or every (other) processor has been set to MPI_COMM_NULL.

### Available communicators

| Handle                  | Description                                                                                   | Derived from            |
| ----------------------- | --------------------------------------------------------------------------------------------- | ----------------------- |
| MPI_COMM_WORLD          | Default global communicator                                                                   | -                       |
| MPI_COMM_PICLAS         | Duplicate of MPI_COMM_WORLD                                                                   | MPI_COMM_PICLAS         |
| MPI_COMM_SHARED         | Processors on a node                                                                          | MPI_COMM_PICLAS         |
| MPI_COMM_LEADERS_SHARED | Group of node leaders (myComputeNodeRank = 0, and the root is a node leader as well)          | MPI_COMM_PICLAS         |
| MPI_COMM_LEADERS_SURF   | Node leaders with surface sides (including halo sides) on a node (might not include the root) | MPI_COMM_LEADERS_SHARED |

#### Feature-specific

| Handle                              | Description                                                            | Derived from    |
| ----------------------------------- | ---------------------------------------------------------------------- | --------------- |
| PartMPIInitGroup(nInitRegions)%COMM | Emission groups                                                        | MPI_COMM_PICLAS |
| SurfCOMM%UNICATOR                   | Processors with a surface side (e.g. reflective), including halo sides | MPI_COMM_PICLAS |
| CPPCOMM%UNICATOR                    | Coupled power potential                                                | MPI_COMM_PICLAS |
| EDC%COMM(iEDCBC)%UNICATOR           | Electric displacement current (per BC)                                 | MPI_COMM_PICLAS |
| FPC%COMM(iUniqueFPCBC)%UNICATOR     | Floating potential (per BC)                                            | MPI_COMM_PICLAS |
| EPC%COMM(iUniqueEPCBC)%UNICATOR     | Electric potential (per BC)                                            | MPI_COMM_PICLAS |
| BiasVoltage%COMM%UNICATOR           | Bias voltage                                                           | MPI_COMM_PICLAS |

#### Flexi-specific

Following communicators are from flexi and not utilized except in `GatheredWriteArray` (`hdf5_output.f90`):

| Handle           | Description                                                                              | Derived from    |
| ---------------- | ---------------------------------------------------------------------------------------- | --------------- |
| MPI_COMM_NODE    | Processors on a node (identical to MPI_COMM_SHARED in PICLas, but not in flexi)          | MPI_COMM_PICLAS |
| MPI_COMM_LEADERS | Group of node leaders (identical to MPI_COMM_LEADERS_SHARED in PICLas, but not in flexi) | MPI_COMM_PICLAS |
| MPI_COMM_WORKERS | All remaining processors, who are not leaders                                            | MPI_COMM_PICLAS |
