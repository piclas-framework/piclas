# MPI Implementation

This chapter describes how PICLas subroutines and functions are parallelized.

## General Remarks: Things to consider

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
    ELEM_TYPE         1    ! HOPR classification of element (currently not used)
    ELEM_ZONE         2    ! Zone as defined in HOPR (currently not used)
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
    SIDE_TYPE         1    ! HOPR classification of side
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
  INTEGER         :: UNICATOR            !< MPI communicator for surface sides
  INTEGER         :: nProcs              !< number of MPI processes
  INTEGER         :: MyRank              !< MyRank, ordered from 0 to nProcs - 1
END TYPE
TYPE (tMPIGROUP)    :: SurfCOMM
```
To create a subset of processors, a condition is required, which is defined by the `color` variable:
```
color = MERGE(1337, MPI_UNDEFINED, nSurfSidesProc.GT.0)
```
Here, every processor with the same `color` will be part of the same communicator. The condition `nSurfSidesProc.GT.0` in this case includes every processor with a surface side. Every other processor will be set to `MPI_UNDEFINED` and consequently be part of `MPI_COMM_NULL`. Now, the communicator itself can be created:
```
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, MPI_INFO_NULL, SurfCOMM%UNICATOR, iError)
```
`MPI_COMM_WORLD` denotes the original communicator containing every processor (but can also be a previously created subset) and the `MPI_INFO_NULL` entry denotes the rank assignment within the new communicator (default: numbering from 0 to nProcs - 1). Additional information can be stored within the created variable:
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
Otherwise load balancing might produce undefined errors that are very difficult to find.