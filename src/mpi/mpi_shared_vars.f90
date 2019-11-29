!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!===================================================================================================================================
!> Contains variables to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_MPI_Shared_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
#if USE_MPI
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: MPISharedInitIsDone=.FALSE.

! Communication
INTEGER            :: myComputeNodeRank               !> Rank of current proc on current compute-node
INTEGER            :: myLeaderGroupRank               !> Rank of compute-node root in compute-node-root comm
INTEGER,ALLOCATABLE:: MPIRankGlobal(:)                !> Array of size nProcessors holding the global rank of each proc
INTEGER,ALLOCATABLE:: MPIRankShared(:)                !> Array of size nProcessors holding the shared rank of each proc
INTEGER            :: nComputeNodeProcessors          !> Number of procs on current compute-node
INTEGER            :: nLeaderGroupProcs               !> Number of nodes
INTEGER            :: nProcessors_Global              !> Number of total procs
INTEGER            :: MPI_COMM_SHARED                 !> Communicator on current compute-node
INTEGER            :: MPI_COMM_LEADERS_SHARED         !> Communicator compute-node roots (my_rank_shared=0)

! Mesh
!> Counters
INTEGER            :: nNonUniqueGlobalSides           !> total nb. of non-unique sides of mesh (hexahedral: 6*nElems)
INTEGER            :: nNonUniqueGlobalNodes           !> total nb. of non-unique nodes of mesh (hexahedral: 8**NGeo * nElems)
INTEGER            :: nNonUniqueGlobalTrees           !> total nb. of trees
INTEGER            :: nComputeNodeElems               !> Number of elems on current compute-node
INTEGER            :: nComputeNodeSides               !> Number of sides on current compute-node
INTEGER            :: nComputeNodeNodes               !> Number of nodes on current compute-node
INTEGER            :: nComputeNodeTotalElems          !> Number of elems on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalSides          !> Number of sides on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalNodes          !> Number of nodes on current compute-node (including halo region)
INTEGER            :: offsetComputeNodeElem           !> elem offset of compute-node root
INTEGER            :: offsetComputeNodeSide           !> side offset of compute-node root
INTEGER            :: offsetComputeNodeNode           !> node offset of compute-node root

INTEGER, ALLOCATABLE :: CNTotalElem2GlobalElem(:) !> Compute Nodes mapping 1:nTotal -> 1:nGlobal
INTEGER, ALLOCATABLE :: GlobalElem2CNTotalElem(:) !> Reverse Mapping


INTEGER,POINTER :: ElemInfo_Shared(:,:)
INTEGER         :: ElemInfo_Shared_Win

INTEGER,POINTER :: SideInfo_Shared(:,:)
INTEGER         :: SideInfo_Shared_Win

INTEGER,POINTER :: NodeInfo_Shared(:)
INTEGER         :: NodeInfo_Shared_Win
REAL,POINTER    :: NodeCoords_Shared(:,:)
INTEGER         :: NodeCoords_Shared_Win

REAL,POINTER    :: xiMinMax_Shared(:,:,:)
INTEGER         :: xiMinMax_Shared_Win
REAL,POINTER    :: TreeCoords_Shared(:,:,:,:,:)
INTEGER         :: TreeCoords_Shared_Win
INTEGER,POINTER :: ElemToTree_Shared(:)
INTEGER         :: ElemToTree_Shared_Win

INTEGER,POINTER :: FIBGM_nElems_Shared(:,:,:)           !> FastInitBackgroundMesh of compute node
INTEGER         :: FIBGM_nElems_Shared_Win
INTEGER,POINTER :: FIBGM_Element_Shared(:)             !> FastInitBackgroundMesh of compute node
INTEGER         :: FIBGM_Element_Shared_Win

REAL,POINTER    :: BoundsOfElem_Shared(:,:,:)          !> Cartesian bounding box around element
INTEGER         :: BoundsOfElem_Shared_Win
INTEGER,POINTER :: ElemToBGM_Shared(:,:)               !> BGM Bounding box around element (respective BGM indices) of compute node
INTEGER         :: ElemToBGM_Shared_Win
INTEGER,POINTER :: FIBGM_offsetElem_Shared(:,:,:)
INTEGER         :: FIBGM_offsetElem_Shared_Win

REAL,POINTER    :: XCL_NGeo_Shared(:,:,:,:,:)
INTEGER         :: XCL_NGeo_Shared_Win
REAL,POINTER    :: dXCL_NGeo_Shared(:,:,:,:,:,:)
INTEGER         :: dXCL_NGeo_Shared_Win
REAL,POINTER    :: BezierControlPoints3D_Shared(:,:,:,:)
INTEGER         :: BezierControlPoints3D_Shared_Win

REAL,POINTER    :: ElemBaryNGeo_Shared(:,:)
INTEGER         :: ElemBaryNGeo_Shared_Win
REAL,POINTER    :: ElemRadiusNGeo_Shared(:)
INTEGER         :: ElemRadiusNGeo_Shared_Win
REAL,POINTER    :: ElemRadius2NGeo_Shared(:)
INTEGER         :: ElemRadius2NGeo_Shared_Win
REAL,POINTER    :: XiEtaZetaBasis_Shared(:,:,:)
INTEGER         :: XiEtaZetaBasis_Shared_Win
REAL,POINTER    :: slenXiEtaZetaBasis_Shared(:,:)
INTEGER         :: slenXiEtaZetaBasis_Shared_Win

REAL,POINTER    :: SideSlabNormals_Shared(:,:,:)
INTEGER         :: SideSlabNormals_Shared_Win
REAL,POINTER    :: SideSlabIntervals_Shared(:,:)
INTEGER         :: SideSlabIntervals_Shared_Win
LOGICAL,POINTER :: BoundingBoxIsEmpty_Shared(:)
INTEGER         :: BoundingBoxIsEmpty_Shared_Win

INTEGER,POINTER :: SideType_Shared(:)
INTEGER         :: SideType_Shared_Win
REAL,POINTER    :: SideDistance_Shared(:)
INTEGER         :: SideDistance_Shared_Win
REAL,POINTER    :: SideNormVec_Shared(:,:)
INTEGER         :: SideNormVec_Shared_Win

#endif /* USE_MPI */
END MODULE
