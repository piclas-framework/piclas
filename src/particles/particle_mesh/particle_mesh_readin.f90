!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

!===================================================================================================================================
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Particle_Mesh_Readin
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
INTEGER,ALLOCATABLE       :: SideInfo_Shared_tmp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! Interfaces
INTERFACE ReadMeshBasics
  MODULE PROCEDURE ReadMeshBasics
END INTERFACE

INTERFACE ReadMeshElems
  MODULE PROCEDURE ReadMeshElems
END INTERFACE

INTERFACE ReadMeshSides
  MODULE PROCEDURE ReadMeshSides
END INTERFACE

INTERFACE ReadMeshSideNeighbors
  MODULE PROCEDURE ReadMeshSideNeighbors
END INTERFACE

INTERFACE ReadMeshNodes
  MODULE PROCEDURE ReadMeshNodes
END INTERFACE

INTERFACE StartCommunicateMeshReadin
  MODULE PROCEDURE StartCommunicateMeshReadin
END INTERFACE

INTERFACE FinishCommunicateMeshReadin
  MODULE PROCEDURE FinishCommunicateMeshReadin
END INTERFACE

INTERFACE FinalizeMeshReadin
  MODULE PROCEDURE FinalizeMeshReadin
END INTERFACE

PUBLIC :: ReadMeshBasics
PUBLIC :: ReadMeshElems
PUBLIC :: ReadMeshSides
PUBLIC :: ReadMeshSideNeighbors
PUBLIC :: ReadMeshNodes
PUBLIC :: StartCommunicateMeshReadin
PUBLIC :: FinishCommunicateMeshReadin
PUBLIC :: FinalizeMeshReadin
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadMeshBasics()
!===================================================================================================================================
! Read basic global counters from mesh file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input                ,ONLY: File_ID,ReadAttribute
USE MOD_Particle_Mesh_Vars        ,ONLY: nNonUniqueGlobalSides,nNonUniqueGlobalNodes
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

!CALL ReadAttribute(File_ID,'nUniqueSides',1,IntScalar=nGlobalUniqueSidesFromMesh)
CALL ReadAttribute(File_ID,'nSides'      ,1,IntScalar=nNonUniqueGlobalSides)
CALL ReadAttribute(File_ID,'nNodes'      ,1,IntScalar=nNonUniqueGlobalNodes)

END SUBROUTINE ReadMeshBasics


SUBROUTINE ReadMeshElems()
!===================================================================================================================================
! Create particle mesh arrays for elems
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  ! allocate shared array for ElemInfo
  CALL Allocate_Shared((/ELEMINFOSIZE,nGlobalElems/),ElemInfo_Shared_Win,ElemInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemInfo_Shared_Win,IERROR)

  ElemInfo_Shared(1:ELEMINFOSIZE_H5,offsetElem+1:offsetElem+nElems) = ElemInfo(:,:)
  ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
#if USE_LOADBALANCE
ELSEIF(UseH5IOLoadBalance)THEN
  ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
END IF
#endif /*USE_LOADBALANCE*/
#endif  /*USE_MPI*/

#if USE_MPI
! broadcast elem offset of compute-node root
offsetComputeNodeElem=offsetElem
CALL MPI_BCAST(offsetComputeNodeElem,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)
#else
! allocate local array for ElemInfo
ALLOCATE(ElemInfo_Shared(1:ELEMINFOSIZE,1:nElems))
ElemInfo_Shared(1:ELEMINFOSIZE_H5,1:nElems) = ElemInfo(:,:)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshElems


SUBROUTINE ReadMeshSides()
!===================================================================================================================================
! Create particle mesh arrays for sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
!===================================================================================================================================

FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND ,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)

ALLOCATE(SideInfo_Shared_tmp(offsetSideID+1:offsetSideID+nSideIDs))
SideInfo_Shared_tmp = 0

#if USE_MPI
! all procs on my compute-node communicate the number of non-unique sides
CALL MPI_ALLREDUCE(nSideIDs,nComputeNodeSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

CALL Allocate_Shared((/SIDEINFOSIZE+1,nNonUniqueGlobalSides/),SideInfo_Shared_Win,SideInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideInfo_Shared_Win,IERROR)
SideInfo_Shared(1                :SIDEINFOSIZE  ,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1,offsetSideID+1:offsetSideID+nSideIDs) = 0
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#else
nComputeNodeSides = nSideIDs
ALLOCATE(SideInfo_Shared(1:SIDEINFOSIZE+1         , 1:nSideIDs))
SideInfo_Shared(1                :SIDEINFOSIZE    , 1:nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1  , 1:nSideIDs) = 0
#endif /*USE_MPI*/

END SUBROUTINE ReadMeshSides


SUBROUTINE ReadMeshSideNeighbors(ElemID,SideID)
!===================================================================================================================================
! Fills temporary array to add side neighbors to SideInfo(_Shared)
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
INTEGER,INTENT(IN)             :: SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (ElemID.EQ.0) THEN
  ! no connection
  SideInfo_Shared_tmp(SideID) = 0
ELSE
#if USE_MPI
  IF (.NOT.ElementOnNode(ElemID)) THEN
    ! neighbour element is outside of compute-node
    SideInfo_Shared_tmp(SideID) = 2
  ELSE
#endif /*USE_MPI*/
    SideInfo_Shared_tmp(SideID) = 1
#if USE_MPI
  END IF
#endif /*USE_MPI*/
END IF

END SUBROUTINE ReadMeshSideNeighbors


SUBROUTINE ReadMeshNodes()
!===================================================================================================================================
! Create particle mesh arrays for nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input         ,ONLY: ReadArray,OpenDataFile
USE MOD_IO_HDF5            ,ONLY: CloseDataFile
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iNode
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: FirstNodeInd,LastNodeInd
INTEGER                        :: nNodeIDs,offsetNodeID
INTEGER,ALLOCATABLE            :: NodeInfo(:),NodeInfoTmp(:,:)
REAL,ALLOCATABLE               :: NodeCoords_indx(:,:)
INTEGER                        :: nNodeInfoIDs,NodeID,NodeCounter
INTEGER                        :: CornerNodeIDswitch(8)
!===================================================================================================================================

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nNodeIDs     => INT(nNodeIDs,IK)     ,&
      offsetNodeID => INT(offsetNodeID,IK) )
  ALLOCATE(NodeCoords_indx(3,nNodeIDs))
  ! read all nodes
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadArray('NodeCoords',2,(/3_IK,nNodeIDs/),offsetNodeID,2,RealArray=NodeCoords_indx)
  CALL CloseDataFile()
END ASSOCIATE

! Keep all nodes if elements are curved
IF (useCurveds.OR.NGeo.EQ.1) THEN
  MeshWasCurved = .TRUE.
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nNodeIDs     => INT(nNodeIDs,IK)     ,&
        offsetNodeID => INT(offsetNodeID,IK) )
    ALLOCATE(NodeInfo(FirstNodeInd:LastNodeInd))
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
    CALL ReadArray('GlobalNodeIDs',1,(/nNodeIDs/),offsetNodeID,1,IntegerArray_i4=NodeInfo)
    CALL CloseDataFile()
  END ASSOCIATE

#if USE_MPI
  ! allocate shared array for NodeInfo
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
  CALL Allocate_Shared((/nNonUniqueGlobalNodes/),NodeInfo_Shared_Win,NodeInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
  NodeInfo_Shared(offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeInfo(:)
  CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win,MPI_COMM_SHARED)

  CALL Allocate_Shared((/3,nNonUniqueGlobalNodes/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
  NodeCoords_Shared(:,offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeCoords_indx(:,:)
#else
  nComputeNodeNodes = nNodeIDs
  ALLOCATE(NodeInfo_Shared(1:nNodeIDs))
  NodeInfo_Shared(1:nNodeIDs) = NodeInfo(:)
  ALLOCATE(NodeCoords_Shared(3,nNodeIDs))
  NodeCoords_Shared(:,:) = NodeCoords_indx(:,:)
#endif  /*USE_MPI*/

! Reduce NodeCoords if no curved elements are to be used
ELSE
  ! every proc needs to allocate the array
  ALLOCATE(NodeInfo(1:nNonUniqueGlobalNodes))

#if USE_MPI
  ! root reads NodeInfo for new mapping
  IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (nNonUniqueGlobalNodes     => INT(nNonUniqueGlobalNodes,IK))
      CALL OpenDataFile(MeshFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
      CALL ReadArray('GlobalNodeIDs',1,(/nNonUniqueGlobalNodes/),0_IK,1,IntegerArray_i4=NodeInfo)
      CALL CloseDataFile()
    END ASSOCIATE
#if USE_MPI
  END IF

  ! root broadcasts NodeInfo to all procs on compute node
  CALL MPI_BCAST(NodeInfo,nNonUniqueGlobalNodes,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

  ! Every proc builds new mapping. This step is required for consistency reasons since ElemInfo is not yet communicated
  nNodeInfoIDs = MAXVAL(NodeInfo)
  ALLOCATE(NodeInfoTmp(2,nNodeInfoIDs))
  NodeInfoTmp = 0

  ! Flag unique node IDs we will keep
  DO iNode = 1,nNonUniqueGlobalNodes
    NodeID = NodeInfo(iNode)
    NodeInfoTmp(1,NodeID) = 1
  END DO

  ! Build new NodeInfo IDs
  NodeCounter = 0
  DO iNode = 1,nNodeInfoIDs
    IF (NodeInfoTmp(1,iNode).EQ.0) CYCLE

    NodeCounter = NodeCounter + 1
    NodeInfoTmp(2,iNode) = NodeCounter
  END DO

#if USE_MPI
  CALL Allocate_Shared((/8*nGlobalElems/),NodeInfo_Shared_Win,NodeInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
#else
  ALLOCATE(NodeInfo_Shared(8*nGlobalElems))
#endif /*USE_MPI*/

  ! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
  CornerNodeIDswitch(1)=1
  CornerNodeIDswitch(2)=(Ngeo+1)
  CornerNodeIDswitch(3)=(Ngeo+1)*Ngeo+1
  CornerNodeIDswitch(4)=(Ngeo+1)**2
  CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
  CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
  CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1
  CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2

  ! New crazy corner node switch (philipesque)
  ASSOCIATE(CNS => CornerNodeIDswitch)

  ! Only the 8 corner nodes count for nodes. (NGeo+1)**2 = 8
  nComputeNodeNodes = 8*nComputeNodeElems

#if USE_MPI
  CALL Allocate_Shared((/3,8*nGlobalElems/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
#else
  ALLOCATE(NodeCoords_Shared(3,8*nGlobalElems))
#endif  /*USE_MPI*/

  ! throw away all nodes except the 8 corner nodes of each hexa
  nNonUniqueGlobalNodes = 8*nGlobalElems

  DO iElem = FirstElemInd,LastElemInd
    FirstNodeInd = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) - offsetNodeID
    ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) = 8*(iElem-1)
    ElemInfo_Shared(ELEM_LASTNODEIND ,iElem) = 8* iElem
    DO iNode = 1,8
      NodeCoords_Shared(:,8*(iElem-1) + iNode) = NodeCoords_indx(:,FirstNodeInd+CNS(iNode))
      NodeInfo_Shared  (  8*(iElem-1) + iNode) = NodeInfoTmp(2,NodeInfo(FirstNodeInd+offsetNodeID+CNS(iNode)))
    END DO
  END DO

  DEALLOCATE(NodeInfoTmp)

  END ASSOCIATE

END IF

! Update node counters
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! scale mesh if desired. Mesh deformation currently not supported!
IF (ABS(meshScale-1.).GT.1e-14) THEN
  NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) = NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) * meshScale
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win  ,MPI_COMM_SHARED)
#endif  /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  DEALLOCATE(NodeInfo)
#if USE_MPI
END IF
#endif /*USE_MPI*/
DEALLOCATE(NodeCoords_indx)

END SUBROUTINE ReadMeshNodes


SUBROUTINE StartCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: StartT
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
#if USE_MPI
INTEGER                        :: iProc
INTEGER                        :: offsetNodeID!,nNodeIDs
#endif /*USE_MPI*/
!===================================================================================================================================

! Start timer: finished in FinishCommunicateMeshReadin()
GETTIME(StartT)

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
! nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  IF (myComputeNodeRank.EQ.0) THEN
    LBWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' Updating mesh on shared memory...'

    ! Arrays for the compute node to hold the elem offsets
    ALLOCATE(displsElem(   0:nLeaderGroupProcs-1), recvcountElem(0:nLeaderGroupProcs-1))
    displsElem(myLeaderGroupRank) = offsetComputeNodeElem
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
    END DO
    recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
  END IF

  ! Broadcast compute node side offset on node
  offsetComputeNodeSide=offsetSideID
  CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  IF (myComputeNodeRank.EQ.0) THEN
    ! Arrays for the compute node to hold the side offsets
    ALLOCATE(displsSide(   0:nLeaderGroupProcs-1), recvcountSide(0:nLeaderGroupProcs-1))
    displsSide(myLeaderGroupRank) = offsetComputeNodeSide
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
    END DO
    recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)

    ! Gather mesh information in a non-blocking way
    ALLOCATE(MPI_COMM_LEADERS_REQUEST(1))
    ! ElemInfo_Shared only needs ELEM_RANK updated, performed in loadbalance_tools.f90
    ! CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE       *recvcountElem  &
    !     ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
        ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)

  END IF

  ! Broadcast compute node node offset on node
  offsetComputeNodeNode=offsetNodeID
  CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  RETURN
END IF
#endif /*USE_LOADBALANCE*/

LBWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' Communicating mesh on shared memory...'

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the elem offsets
  ALLOCATE(displsElem(   0:nLeaderGroupProcs-1), recvcountElem(0:nLeaderGroupProcs-1))
  displsElem(myLeaderGroupRank) = offsetComputeNodeElem
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
  END DO
  recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
END IF

! Broadcast compute node side offset on node
offsetComputeNodeSide=offsetSideID
CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the side offsets
  ALLOCATE(displsSide(   0:nLeaderGroupProcs-1), recvcountSide(0:nLeaderGroupProcs-1))
  displsSide(myLeaderGroupRank) = offsetComputeNodeSide
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
  END DO
  recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)
END IF

! Broadcast compute node node offset on node
offsetComputeNodeNode=offsetNodeID
CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsNode(   0:nLeaderGroupProcs-1), recvcountNode(0:nLeaderGroupProcs-1))
  displsNode(myLeaderGroupRank) = offsetComputeNodeNode
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNode,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountNode(iProc-1) = displsNode(iProc)-displsNode(iProc-1)
  END DO
  recvcountNode(nLeaderGroupProcs-1) = nNonUniqueGlobalNodes - displsNode(nLeaderGroupProcs-1)

  ! Gather mesh information in a non-blocking way
  ALLOCATE(MPI_COMM_LEADERS_REQUEST(1:4))
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE    *recvcountElem  &
      ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)*recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(2),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeInfo_Shared,                 recvcountNode  &
      ,displsNode                  ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(3),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeCoords_Shared,3             *recvcountNode  &
      ,3*displsNode                ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(4),IERROR)
END IF
#endif /*USE_MPI*/

END SUBROUTINE StartCommunicateMeshReadin


SUBROUTINE FinishCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: CommMeshReadinWallTime,StartT
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: FirstElemInd,LastElemInd
INTEGER :: iElem,NbElemID
INTEGER :: nSideIDs,offsetSideID
INTEGER :: iSide,sideCount
INTEGER :: iLocSide,jLocSide,nlocSides,nlocSidesNb,NbSideID
REAL    :: EndT
!===================================================================================================================================

#if USE_MPI
! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  ! Finish non-blocking mesh communication
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_WAITALL(1,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
    DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
  END IF

  ! final sync of all mesh shared arrays
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  ! Write compute-node local SIDE_NBELEMTYPE
  IF (myComputeNodeRank.EQ.0) THEN
    SideInfo_Shared(SIDE_NBELEMTYPE,:) = 0
  END IF
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp
  DEALLOCATE(SideInfo_Shared_tmp)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  RETURN
END IF
#endif  /*USE_LOADBALANCE*/

#if USE_MPI
! Finish non-blocking mesh communication
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_WAITALL(4,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
  DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
END IF

! Ensure communication for determination of SIDE_LOCALID
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

DO iElem = FirstElemInd,LastElemInd
  iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  SideInfo_Shared(SIDE_ELEMID,iSide+1:ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)) = iElem
  sideCount = 0
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    ! Big mortar side
    ! 104:bilinear,204:curved
    IF (SideInfo_Shared(SIDE_TYPE,iSide).LE.100) THEN
      sideCount = sideCount + 1
      SideInfo_Shared(SIDE_LOCALID,iSide) = sideCount
    ELSE
      ! Mortar case
      SideInfo_Shared(SIDE_LOCALID,iSide) = -1
    END IF
    ! Check all sides on the small element side to find the small mortar side pointing back
    NbElemID    = SideInfo_Shared(SIDE_NBELEMID,iSide)
    IF(NbElemID.EQ.0) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = 0
    ELSE IF (NbElemID.LE.-1) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = -1
    ELSE
      nlocSidesNb = ElemInfo_Shared(ELEM_LASTSIDEIND,NbElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID)
      DO jLocSide = 1,nlocSidesNb
        NbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID) + jLocSide
        IF (ABS(SideInfo_Shared(SIDE_ID,iSide)).EQ.ABS(SideInfo_Shared(SIDE_ID,NbSideID))) THEN
          SideInfo_Shared(SIDE_NBSIDEID,iSide) = NbSideID
          EXIT
        END IF
      END DO
    END IF
  END DO
END DO

#if USE_MPI
! Perform second communication step to distribute updated SIDE_LOCALID
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)*recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
END IF

! Write compute-node local SIDE_NBELEMTYPE
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  SideInfo_Shared(SIDE_NBELEMTYPE,:) = 0
END IF
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp

! final sync of all mesh shared arrays
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED)
#endif  /*USE_MPI*/

nUniqueGlobalNodes = MAXVAL(NodeInfo_Shared)
SDEALLOCATE(SideInfo_Shared_tmp)

GETTIME(EndT)
CommMeshReadinWallTime=EndT-StartT
CALL DisplayMessageAndTime(CommMeshReadinWallTime, 'DONE!')

END SUBROUTINE FinishCommunicateMeshReadin


SUBROUTINE FinalizeMeshReadin()
!===================================================================================================================================
! Finalizes the shared mesh readin
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! volumes
CALL UNLOCK_AND_FREE(ElemVolume_Shared_Win)
CALL UNLOCK_AND_FREE(ElemCharLength_Shared_Win)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
ADEALLOCATE(ElemVolume_Shared)
ADEALLOCATE(ElemCharLength_Shared)

#if USE_MPI
! Free communication arrays
SDEALLOCATE(displsElem)
SDEALLOCATE(recvcountElem)
SDEALLOCATE(displsSide)
SDEALLOCATE(recvcountSide)

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

! elems
CALL UNLOCK_AND_FREE(ElemInfo_Shared_Win)

! sides
CALL UNLOCK_AND_FREE(SideInfo_Shared_Win)

! nodes
CALL UNLOCK_AND_FREE(NodeInfo_Shared_Win)
CALL UNLOCK_AND_FREE(NodeCoords_Shared_Win)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

ADEALLOCATE(ElemInfo_Shared)
ADEALLOCATE(SideInfo_Shared)
ADEALLOCATE(NodeInfo_Shared)
ADEALLOCATE(NodeCoords_Shared)

! Free communication arrays
#if USE_MPI
SDEALLOCATE(displsNode)
SDEALLOCATE(recvcountNode)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeMeshReadin


END MODULE MOD_Particle_Mesh_Readin
