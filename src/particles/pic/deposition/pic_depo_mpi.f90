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

MODULE MOD_PICDepo_MPI
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE

#if USE_MPI
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
TYPE NodeDepoMapping
  INTEGER                                     :: NodeID
  TYPE (NodeDepoMapping), POINTER             :: next => NULL()
END TYPE NodeDepoMapping
!===================================================================================================================================
PUBLIC :: InitDepoNodesMPI
PUBLIC :: InitDepoSurfNodesMPI
PUBLIC :: ExchangeNodeSource
PUBLIC :: ExchangeNodeSourceExtMPI
PUBLIC :: CollectSurfNodeAreaOnMPIRoot
PUBLIC :: ExchangeSurfNodeSourceMPI
PUBLIC :: LBReverseExchangeSurfNodeSource
PUBLIC :: ReverseExchangeSurfNodeArea
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the MPI communication for the volume node deposition
!===================================================================================================================================
SUBROUTINE InitDepoNodesMPI(DoNodeMapping,SendNode)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_PICDepo_Vars
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: nUniqueGlobalNodes
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_Mesh_Vars              ,ONLY: offsetElem,ELEM_RANK
USE MOD_Particle_Mesh_Vars     ,ONLY: NodeToElemInfo,NodeToElemMapping,ElemNodeID_Shared,NodeInfo_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nProcessors_Global
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(INOUT) :: DoNodeMapping(0:nProcessors_Global-1)
LOGICAL,INTENT(INOUT) :: SendNode(1:nUniqueGlobalNodes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iElem, iNode
INTEGER                   :: UniqueNodeID, testNode
INTEGER                   :: GlobalRankToNodeSendDepoRank(0:nProcessors_Global-1)
INTEGER                   :: jElem,TestElemID
INTEGER                   :: NonUniqueNodeID
INTEGER                   :: SendNodeCount, GlobalElemRank, iProc
INTEGER                   :: GlobalElemRankOrig, iRank
LOGICAL,ALLOCATABLE       :: IsDepoNode(:)
LOGICAL                   :: bordersMyrank
! Non-symmetric particle exchange
TYPE(MPI_Request)         :: SendRequestNonSymDepo(0:nProcessors_Global-1)      , RecvRequestNonSymDepo(0:nProcessors_Global-1)
INTEGER                   :: nSendUniqueNodesNonSymDepo(0:nProcessors_Global-1) , nRecvUniqueNodesNonSymDepo(0:nProcessors_Global-1)
TYPE tElemNodeDepoMap
  TYPE (NodeDepoMapping), POINTER :: first => NULL()
  LOGICAL               :: firstNode
  INTEGER               :: nNodes
END TYPE tElemNodeDepoMap
TYPE(tElemNodeDepoMap), ALLOCATABLE :: ElemNodeDepoMap(:)
TYPE(NodeDepoMapping), POINTER :: node
!===================================================================================================================================
IF(DoDielectricSurfaceCharge)THEN
  ALLOCATE(NodeSourceExtMPI(1:nUniqueGlobalNodes))
  NodeSourceExtMPI = 0.
END IF ! DoDielectricSurfaceCharge

! Loop over the elements of the complete compute-node region (including the halo region)
DO iElem = 1,nComputeNodeTotalElems
  IF (FlagShapeElem(iElem)) THEN
    bordersMyrank = .FALSE.
    ! Loop all local nodes
    TestElemID = GetGlobalElemID(iElem)
    GlobalElemRankOrig = ElemInfo_Shared(ELEM_RANK,TestElemID)
    IF (DoHaloDepo.AND.(GlobalElemRankOrig.NE.myRank)) DoNodeMapping(GlobalElemRankOrig) = .TRUE.

    DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
    ! Loop 1D array [offset + 1 : offset + NbrOfElems]
    ! (all CN elements that are connected to the local nodes)
    DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
      TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
      GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
      IF (DoHaloDepo) THEN
        SendNode(UniqueNodeID) = .TRUE.
        IF (GlobalElemRank.NE.myRank) DoNodeMapping(GlobalElemRank) = .TRUE.
      ELSE
        IF (GlobalElemRank.EQ.myRank) THEN
          bordersMyrank = .TRUE.
          SendNode(UniqueNodeID) = .TRUE.
        END IF
      END IF
    END DO
    IF (.NOT.DoHaloDepo.AND.bordersMyrank) THEN
      DoNodeMapping(GlobalElemRankOrig) = .TRUE.
    END IF
    END DO
  END IF
END DO

! Flag the unique deposition nodes per processor
nDepoNodes = 0
ALLOCATE(IsDepoNode(1:nUniqueGlobalNodes))
IsDepoNode = .FALSE.
DO iElem =1, nElems
  TestElemID = GetCNElemID(iElem + offsetElem)
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,TestElemID)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
    IsDepoNode(UniqueNodeID) = .TRUE.
  END DO
END DO
! Count the number of unique deposition nodes per processor
nDepoNodes = COUNT(IsDepoNode)
! Add number of nodes to be sent
nDepoNodesTotal = nDepoNodes
DO iNode=1, nUniqueGlobalNodes
  IF (.NOT.IsDepoNode(iNode).AND.SendNode(iNode)) THEN
    nDepoNodesTotal = nDepoNodesTotal + 1
  END IF
END DO
! Create mapping from unique deposition node to global unique node
ALLOCATE(DepoNodetoGlobalNode(1:nDepoNodesTotal))
nDepoNodesTotal = 0
DO iNode=1, nUniqueGlobalNodes
  IF (IsDepoNode(iNode)) THEN
    nDepoNodesTotal = nDepoNodesTotal + 1
    DepoNodetoGlobalNode(nDepoNodesTotal) = iNode
  END IF
END DO
DO iNode=1, nUniqueGlobalNodes
  IF (.NOT.IsDepoNode(iNode).AND.SendNode(iNode)) THEN
    nDepoNodesTotal = nDepoNodesTotal + 1
    DepoNodetoGlobalNode(nDepoNodesTotal) = iNode
  END IF
END DO
! Create mapping of exchange processor rank to global rank
GlobalRankToNodeSendDepoRank = -1
nNodeSendExchangeProcs = COUNT(DoNodeMapping)
ALLOCATE(NodeSendDepoRankToGlobalRank(1:nNodeSendExchangeProcs))
NodeSendDepoRankToGlobalRank = 0
nNodeSendExchangeProcs = 0
DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (DoNodeMapping(iRank)) THEN
    nNodeSendExchangeProcs = nNodeSendExchangeProcs + 1
    GlobalRankToNodeSendDepoRank(iRank) = nNodeSendExchangeProcs
    NodeSendDepoRankToGlobalRank(nNodeSendExchangeProcs) = iRank
  END IF
END DO
! ALLOCATE(NodeDepoMapping(1:nNodeSendExchangeProcs, 1:nUniqueGlobalNodes))
! NodeDepoMapping = .FALSE.
ALLOCATE(ElemNodeDepoMap(1:nNodeSendExchangeProcs))
ElemNodeDepoMap(:)%firstNode = .TRUE.
ElemNodeDepoMap(:)%nNodes = 0

DO iNode = 1, nUniqueGlobalNodes
  IF (SendNode(iNode)) THEN
    ElemLoop: DO jElem = NodeToElemMapping(1,iNode) + 1, NodeToElemMapping(1,iNode) + NodeToElemMapping(2,iNode)
      TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
      GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
      IF (GlobalElemRank.NE.myRank) THEN
        iRank = GlobalRankToNodeSendDepoRank(GlobalElemRank)
        IF (iRank.LT.1) CALL ABORT(__STAMP__,'Found not connected Rank!', myRank)
        ! NodeDepoMapping(iRank, iNode) = .TRUE.
        IF (ElemNodeDepoMap(iRank)%firstNode) THEN
          ElemNodeDepoMap(iRank)%firstNode = .FALSE.
          ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
          ALLOCATE(ElemNodeDepoMap(iRank)%first)
          ElemNodeDepoMap(iRank)%first%NodeID = iNode
        ELSE
          ! Check if node already exists
          node => ElemNodeDepoMap(iRank)%first
          DO testNode = 1, ElemNodeDepoMap(iRank)%nNodes
          IF (node%NodeID.EQ.iNode) CYCLE ElemLoop
          IF (.NOT.ASSOCIATED(node%next)) EXIT
          node => node%next
          END DO
          ! Add new node at the end of the list
          ALLOCATE(node%next)
          node%next%NodeID = iNode
          ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
        END IF
      END IF
    END DO ElemLoop
  END IF
END DO
! Get number of send nodes for each proc: Size of each message for each proc for deposition
nSendUniqueNodesNonSymDepo         = 0
nRecvUniqueNodesNonSymDepo(myrank) = 0
ALLOCATE(NodeMappingSend(1:nNodeSendExchangeProcs))
DO iProc = 1, nNodeSendExchangeProcs
  NodeMappingSend(iProc)%nSendUniqueNodes = 0
  ! DO iNode = 1, nUniqueGlobalNodes
  !   IF (NodeDepoMapping(iProc,iNode)) NodeMappingSend(iProc)%nSendUniqueNodes = NodeMappingSend(iProc)%nSendUniqueNodes + 1
  ! END DO
  NodeMappingSend(iProc)%nSendUniqueNodes =  ElemNodeDepoMap(iProc)%nNodes
  ! local to global array
  nSendUniqueNodesNonSymDepo(NodeSendDepoRankToGlobalRank(iProc)) = NodeMappingSend(iProc)%nSendUniqueNodes
END DO

! Open receive buffer for non-symmetric exchange identification
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_IRECV( nRecvUniqueNodesNonSymDepo(iProc)  &
    , 1                                              &
    , MPI_INTEGER                                    &
    , iProc                                          &
    , 2000                                           &
    , MPI_COMM_PICLAS                                &
    , RecvRequestNonSymDepo(iProc)                   &
    , IERROR)
END DO

! Send each proc the number of nodes that can be reached by deposition
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_ISEND( nSendUniqueNodesNonSymDepo(iProc) &
    , 1                                             &
    , MPI_INTEGER                                   &
    , iProc                                         &
    , 2000                                          &
    , MPI_COMM_PICLAS                               &
    , SendRequestNonSymDepo(iProc)                  &
    , IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

nNodeRecvExchangeProcs = COUNT(nRecvUniqueNodesNonSymDepo.GT.0)
ALLOCATE(NodeMappingRecv(1:nNodeRecvExchangeProcs))
ALLOCATE(NodeRecvDepoRankToGlobalRank(1:nNodeRecvExchangeProcs))
NodeRecvDepoRankToGlobalRank = 0
nNodeRecvExchangeProcs = 0
DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (nRecvUniqueNodesNonSymDepo(iRank).GT.0) THEN
    nNodeRecvExchangeProcs = nNodeRecvExchangeProcs + 1
    ! Store global rank of iRecvRank
    NodeRecvDepoRankToGlobalRank(nNodeRecvExchangeProcs) = iRank
    ! Store number of nodes of iRecvRank
    NodeMappingRecv(nNodeRecvExchangeProcs)%nRecvUniqueNodes = nRecvUniqueNodesNonSymDepo(iRank)
  END IF
END DO

! Open receive buffer
ALLOCATE(RecvRequest(1:nNodeRecvExchangeProcs))
DO iProc = 1, nNodeRecvExchangeProcs
  ALLOCATE(NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceCharge(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceCurrent(1:3,1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMappingRecv(iProc)%RecvNodeSourceExt(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID &
    , NodeMappingRecv(iProc)%nRecvUniqueNodes                   &
    , MPI_INTEGER                                               &
    , NodeRecvDepoRankToGlobalRank(iProc)                       &
    , 666                                                       &
    , MPI_COMM_PICLAS                                           &
    , RecvRequest(iProc)                                        &
    , IERROR)
END DO

! Open send buffer
ALLOCATE(SendRequest(1:nNodeSendExchangeProcs))
DO iProc = 1, nNodeSendExchangeProcs
  ALLOCATE(NodeMappingSend(iProc)%SendNodeUniqueGlobalID(1:NodeMappingSend(iProc)%nSendUniqueNodes))
  NodeMappingSend(iProc)%SendNodeUniqueGlobalID=-1
  ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceCharge(1:NodeMappingSend(iProc)%nSendUniqueNodes))
  NodeMappingSend(iProc)%SendNodeSourceCharge=0.
  ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceCurrent(1:3,1:NodeMappingSend(iProc)%nSendUniqueNodes))
  NodeMappingSend(iProc)%SendNodeSourceCurrent=0.
  IF(DoDielectricSurfaceCharge) ALLOCATE(NodeMappingSend(iProc)%SendNodeSourceExt(1:NodeMappingSend(iProc)%nSendUniqueNodes))
  SendNodeCount = 0
  ! DO iNode = 1, nUniqueGlobalNodes
  !   IF (NodeDepoMapping(iProc,iNode)) THEN
  !     SendNodeCount = SendNodeCount + 1
  !     NodeMappingSend(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = iNode
  !   END IF
  ! END DO
  ! ALLOCATE(node)
  ! node => ElemNodeDepoMap(iProc)%first
  ! DO testNode = 1, ElemNodeDepoMap(iProc)%nNodes
  !   SendNodeCount = SendNodeCount + 1
  !   NodeMappingSend(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = node%NodeID
  !   node => node%next
  ! END DO

  ! First loop: Traverse the list and populate NodeMappingSend
  node => ElemNodeDepoMap(iProc)%first
  DO WHILE (ASSOCIATED(node))
    SendNodeCount = SendNodeCount + 1
    NodeMappingSend(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = node%NodeID
    node => node%next
  END DO

  ! node => ElemNodeDepoMap(iProc)%first
  ! DO testNode = 1, ElemNodeDepoMap(iProc)%nNodes
  !   ElemNodeDepoMap(iProc)%first => ElemNodeDepoMap(iProc)%first%next
  !   DEALLOCATE(node)
  !   node => ElemNodeDepoMap(iProc)%first
  ! END DO
  ! IF(ASSOCIATED(ElemNodeDepoMap(iProc)%first)) THEN
  !   DEALLOCATE(ElemNodeDepoMap(iProc)%first)
  ! END IF
  ! IF(ASSOCIATED(node)) THEN
  !   DEALLOCATE(node)
  ! END IF

  ! Deallocate the list
  CALL DeallocateNodeList(ElemNodeDepoMap(iProc)%first)
  NULLIFY(ElemNodeDepoMap(iProc)%first)
  ElemNodeDepoMap(iProc)%nNodes = 0

  CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeUniqueGlobalID                   &
    , NodeMappingSend(iProc)%nSendUniqueNodes                         &
    , MPI_INTEGER                                                 &
    , NodeSendDepoRankToGlobalRank(iProc)                         &
    , 666                                                         &
    , MPI_COMM_PICLAS                                              &
    , SendRequest(iProc)                                          &
    , IERROR)
END DO

! Finish send
DO iProc = 1, nNodeSendExchangeProcs
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

! Finish receive
DO iProc = 1, nNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

END SUBROUTINE InitDepoNodesMPI


!===================================================================================================================================
!> Initialize the MPI communication for the 2D surface node deposition
!>
! 1.) Identify communication partners
! 2.) Create mapping of exchange processor rank to global rank depending on CommunicateWithRank(iRank)
! 3.) Loop over the send FEM vertices and each connected processes and build linked list of node IDs and count them
! 4.) Get the number of send nodes for each communication partner: Size of each message for each process for deposition
! 5.) MPI send/receive the number of deposition nodes
! 6.) From the received messages, determine the message size that is sent from each communication partner.
! 7.) MPI send/receive the vertex IDs of deposition nodes
!===================================================================================================================================
SUBROUTINE InitDepoSurfNodesMPI()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_PICDepo_Vars           ,ONLY: IsDepoSurfNode,SurfNodeSourceMPI,nDepoSurfNodesTotal
USE MOD_PICDepo_Vars           ,ONLY: nSurfNodeRecvExchangeProcs,nSurfNodeSendExchangeProcs
USE MOD_PICDepo_Vars           ,ONLY: SurfRecvRequest,SurfNodeMappingRecv,SurfNodeRecvDepoRankToGlobalRank
USE MOD_PICDepo_Vars           ,ONLY: SurfSendRequest,SurfNodeMappingSend,SurfNodeSendDepoRankToGlobalRank
USE MOD_Mesh_Vars              ,ONLY: nFEMVertices
USE MOD_Particle_Mesh_Vars     ,ONLY: VertexConnectInfo_Shared
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_Mesh_Vars              ,ONLY: ELEM_RANK
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars        ,ONLY: nProcessors_Global
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,VertexInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL :: CommunicateWithRank(0:nProcessors_Global-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iCNElem
INTEGER                   :: testNode
INTEGER                   :: GlobalRankToNodeSendDepoRank(0:nProcessors_Global-1)
INTEGER                   :: SendNodeCount,iProc
INTEGER                   :: iRank
LOGICAL,ALLOCATABLE       :: FEMVertexIDisDone(:)
LOGICAL                   :: NodeAlreadyAssignedToRoot
! Non-symmetric particle exchange
TYPE(MPI_Request)         :: SendRequestNonSymDepo(0:nProcessors_Global-1)      , RecvRequestNonSymDepo(0:nProcessors_Global-1)
INTEGER                   :: nSendUniqueNodesNonSymDepo(0:nProcessors_Global-1) , nRecvUniqueNodesNonSymDepo(0:nProcessors_Global-1)
TYPE tElemNodeDepoMap
  TYPE (NodeDepoMapping), POINTER :: first => NULL()
  LOGICAL               :: firstNode
  INTEGER               :: nNodes
END TYPE tElemNodeDepoMap
TYPE(tElemNodeDepoMap), ALLOCATABLE :: ElemNodeDepoMap(:)
TYPE(NodeDepoMapping), POINTER :: node
INTEGER :: iVertexConnect,GlobalNbElemID
INTEGER :: GlobalNBElemRank
INTEGER :: FirstVertexInd,LastVertexInd,FirstVertexConnectInd,LastVertexConnectInd
INTEGER :: FEMVertexID !< Super unique node ID (folds neighbouring and periodic nodes into a single node index)
INTEGER :: iVertexInd  !< NonUniqueVertexID
INTEGER :: iGlobalElemID
!===================================================================================================================================
! Allocate container for storing the local non-synchronized surface charge, which is always nullified after
! communication/synchronization with other processes
ALLOCATE(SurfNodeSourceMPI(1:nDepoSurfNodesTotal))
SurfNodeSourceMPI = 0.

! Only continue to the communication part when there are multiple processes
IF(nProcessors.LE.1) RETURN

! Nullify container to flag each process if it will receive charge
CommunicateWithRank = .FALSE.
! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
! Step 1 of 2: Force every process to establish a communication with MPIRoot
IF(myrank.NE.0) CommunicateWithRank(0) = .TRUE.

! 1.) Identify communication partners
! Loop over the elements of the complete compute-node region (including the halo region) where the node can deposit charge
! and find the process ranks that have a FEM vertex to which deposition from the myrank might occur
DO iCNElem = 1,nComputeNodeTotalElems
  ! Loop all local nodes
  iGlobalElemID = GetGlobalElemID(iCNElem)
  ! Get local FEMElemInfo of current element
  FirstVertexInd = ElemInfo_Shared(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1 ! this comes from FEMElemInfo() from mesh.h5
  LastVertexInd  = ElemInfo_Shared(ELEM_LASTVERTEXIND,iGlobalElemID)    ! this comes from FEMElemInfo() from mesh.h5
  ! Loop over all non-unique vertices (the total number via iGlobalElemID and iVertexInd corresponds to nVertices in .h5)
  iVertexIndLoop: DO iVertexInd = FirstVertexInd,LastVertexInd
    ! Get topologically unique global vertex ID (via VertexInfo from mesh.h5), includes periodicity (needed for a FEM solver)
    FEMVertexID = VertexInfo_Shared(VERTEX_FEMID,iVertexInd)
    ! Skip vertices without deposition
    IF(.NOT.IsDepoSurfNode(FEMVertexID)) CYCLE iVertexIndLoop ! go to next vertex
    ! Get local vertex connectivity: First and Last connected vertex index
    FirstVertexConnectInd = VertexInfo_Shared(VERTEX_FIRSTCONNECTIND,iVertexInd)+1
    LastVertexConnectInd  = VertexInfo_Shared(VERTEX_LASTCONNECTIND,iVertexInd)
    iVertexConnectLoop: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
      ! Get neighbour infos. Note the ABS() for +/- master/slave notation
      GlobalNbElemID = ABS(VertexConnectInfo_Shared(VERTEXCONNECT_NBELEMID,iVertexConnect))
      ! Do not consider myself
      IF(GlobalNbElemID.EQ.iGlobalElemID) CYCLE iVertexConnectLoop ! go to next connection
      ! Get neighbour rank
      GlobalNBElemRank = ElemInfo_Shared(ELEM_RANK,GlobalNbElemID)
      ! Communicate with this processes
      IF (GlobalNBElemRank.NE.myrank) THEN
        ! Flag the communication partner
        CommunicateWithRank(GlobalNBElemRank) = .TRUE.
      END IF ! GlobalNbElemID.NE.myrank
    END DO iVertexConnectLoop ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
  END DO iVertexIndLoop ! iVertexInd = iFirstVertexInd,LastVertexInd
END DO ! iCNElem = 1,nComputeNodeTotalElems

! 2.) Create mapping of exchange processor rank to global rank depending on CommunicateWithRank(iRank)
!     Initialize mapping global rank to i-th communication partner
GlobalRankToNodeSendDepoRank = -1
! Count the number of processes to communicate with
nSurfNodeSendExchangeProcs = COUNT(CommunicateWithRank)
! Allocate and initialize mapping from i-th communication partner to global rank
ALLOCATE(SurfNodeSendDepoRankToGlobalRank(1:nSurfNodeSendExchangeProcs))
SurfNodeSendDepoRankToGlobalRank = 0
! Nullify the iteration counter
nSurfNodeSendExchangeProcs = 0
! Loop over the global number of procsses
DO iRank= 0, nProcessors_Global-1
  ! Ignore myself
  IF (iRank.EQ.myRank) CYCLE
  ! Only consider communication partners identified in step 1.)
  IF (CommunicateWithRank(iRank)) THEN
    ! Increment the exchange proc index
    nSurfNodeSendExchangeProcs = nSurfNodeSendExchangeProcs + 1
    ! Create mapping from global rank to i-th communication partner
    GlobalRankToNodeSendDepoRank(iRank) = nSurfNodeSendExchangeProcs
    ! Create opposing mapping from i-th communication partner to global rank
    SurfNodeSendDepoRankToGlobalRank(nSurfNodeSendExchangeProcs) = iRank
  END IF
END DO

! 3.) Loop over the send FEM vertices and each connected processes and build linked list of node IDs and their number
! Allocate container with entries for each communication partner
ALLOCATE(ElemNodeDepoMap(1:nSurfNodeSendExchangeProcs))
ElemNodeDepoMap(:)%firstNode = .TRUE.
ElemNodeDepoMap(:)%nNodes = 0
! Initialize with false
ALLOCATE(FEMVertexIDisDone(1:nFEMVertices))
FEMVertexIDisDone = .FALSE.
! Loop over all CN elements, where deposition might occur by the current process
DO iCNElem = 1,nComputeNodeTotalElems
  ! Loop all local nodes
  iGlobalElemID = GetGlobalElemID(iCNElem)
  ! Get local FEMElemInfo of current element
  FirstVertexInd = ElemInfo_Shared(ELEM_FIRSTVERTEXIND,iGlobalElemID)+1 ! this comes from FEMElemInfo() from mesh.h5
  LastVertexInd  = ElemInfo_Shared(ELEM_LASTVERTEXIND,iGlobalElemID)    ! this comes from FEMElemInfo() from mesh.h5
  ! Loop over all non-unique vertices (the total number via iGlobalElemID and iVertexInd corresponds to nVertices in .h5)
  iVertexIndLoop2: DO iVertexInd = FirstVertexInd,LastVertexInd
    ! Get topologically unique global vertex ID (via VertexInfo from mesh.h5), includes periodicity (needed for a FEM solver)
    FEMVertexID = VertexInfo_Shared(VERTEX_FEMID,iVertexInd)
    ! Skip vertices without deposition
    IF(.NOT.IsDepoSurfNode(FEMVertexID)) CYCLE iVertexIndLoop2 ! go to next vertex
    ! Skip vertices that have already been processes
    ! NOTE: The question is still open, if this is required or if this is wrong because some vertices are then not considered and
    ! subsequently connections are not established (might lead to deadlock in MPI_Wait or wrong surface charge)
    IF(FEMVertexIDisDone(FEMVertexID)) CYCLE iVertexIndLoop2 ! go to next vertex
    ! Flag the FEM vertex
    FEMVertexIDisDone(FEMVertexID) = .TRUE.
    ! Get local vertex connectivity: First and Last connected vertex index
    FirstVertexConnectInd = VertexInfo_Shared(VERTEX_FIRSTCONNECTIND,iVertexInd)+1
    LastVertexConnectInd  = VertexInfo_Shared(VERTEX_LASTCONNECTIND,iVertexInd)
    ! Loop over the connections of the vertex
    iVertexConnectLoop2: DO iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
      ! Get neighbour info. Note the ABS() for +/- master/slave notation
      GlobalNbElemID = ABS(VertexConnectInfo_Shared(VERTEXCONNECT_NBELEMID,iVertexConnect))
      ! Do not consider myself
      IF(GlobalNbElemID.EQ.iGlobalElemID) CYCLE iVertexConnectLoop2 ! go to next connection
      ! Get neighbour rank
      GlobalNBElemRank = ElemInfo_Shared(ELEM_RANK,GlobalNbElemID)
      ! Communicate with this processes if it is not myself
      IF (GlobalNBElemRank.NE.myrank) THEN
        ! Get index of the i-th connection partner
        iRank = GlobalRankToNodeSendDepoRank(GlobalNBElemRank)
        ! Sanity check
        IF (iRank.LT.1) CALL ABORT(__STAMP__,'Found not connected Rank!', myRank)
        ! CHeck if the first node for this process is encountered
        IF (ElemNodeDepoMap(iRank)%firstNode) THEN
          ! Flip the first node flag to false
          ElemNodeDepoMap(iRank)%firstNode = .FALSE.
          ! Increment the number of nodes for this communication partner
          ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
          ! Allocate the next link
          ALLOCATE(ElemNodeDepoMap(iRank)%first)
          ! Store the FEMVertexID
          ElemNodeDepoMap(iRank)%first%NodeID = FEMVertexID
        ELSE ! 2nd node encountered
          ! Check if node already exists
          node => ElemNodeDepoMap(iRank)%first
          ! Loop over the stored FEMVertexIDs to not store the same FEMVertexID twice
          DO testNode = 1, ElemNodeDepoMap(iRank)%nNodes
            ! Check for FEMVertexID
            IF (node%NodeID.EQ.FEMVertexID) CYCLE iVertexConnectLoop2 ! Jump to the next connection
            ! Check if the end of the list if encountered
            IF (.NOT.ASSOCIATED(node%next)) EXIT
            ! Next link
            node => node%next
          END DO ! testNode = 1, ElemNodeDepoMap(iRank)%nNodes
          ! Add new node at the end of the list
          ALLOCATE(node%next)
          ! Store the FEMVertexID
          node%next%NodeID = FEMVertexID
          ! Increment the number of nodes for this communication partner
          ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
        END IF ! ElemNodeDepoMap(iRank)%firstNode
      END IF ! GlobalNbElemID.NE.myrank
    END DO iVertexConnectLoop2 ! iVertexConnect = FirstVertexConnectInd, LastVertexConnectInd
    ! ========================================================================================================
    ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
    ! Step 2 of 2: Force every process to send this FEMVertexID to MPIRoot
    ! Remove this link in the furute and replace with a gathered I/O or something different
    GlobalNBElemRank = 0 ! Create artificial link to MPIRoot
    ! Communicate with this processes if it is not myself
    IF (GlobalNBElemRank.NE.myrank) THEN
      ! Get index of the i-th connection partner
      iRank = GlobalRankToNodeSendDepoRank(GlobalNBElemRank)
      ! Sanity check
      IF (iRank.LT.1) CALL ABORT(__STAMP__,'Found not connected Rank!', myRank)
      ! Check if the first node for this process is encountered
      IF (ElemNodeDepoMap(iRank)%firstNode) THEN
        ! Flip the first node flag to false
        ElemNodeDepoMap(iRank)%firstNode = .FALSE.
        ! Increment the number of nodes for this communication partner
        ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
        ! Allocate the next link
        ALLOCATE(ElemNodeDepoMap(iRank)%first)
        ! Store the FEMVertexID
        ElemNodeDepoMap(iRank)%first%NodeID = FEMVertexID
      ELSE ! 2nd node encountered
        ! Check if node already exists
        node => ElemNodeDepoMap(iRank)%first
        ! Loop over the stored FEMVertexIDs to not store the same FEMVertexID twice
        NodeAlreadyAssignedToRoot=.FALSE.
        DO testNode = 1, ElemNodeDepoMap(iRank)%nNodes
          ! Check for FEMVertexID
          IF (node%NodeID.EQ.FEMVertexID) NodeAlreadyAssignedToRoot=.TRUE.
          ! Check if the end of the list if encountered
          IF (.NOT.ASSOCIATED(node%next)) EXIT
          ! Next link
          node => node%next
        END DO ! testNode = 1, ElemNodeDepoMap(iRank)%nNodes
        IF (.NOT.NodeAlreadyAssignedToRoot) THEN
          ! Add new node at the end of the list
          ALLOCATE(node%next)
          ! Store the FEMVertexID
          node%next%NodeID = FEMVertexID
          ! Increment the number of nodes for this communication partner
          ElemNodeDepoMap(iRank)%nNodes = ElemNodeDepoMap(iRank)%nNodes + 1
        END IF ! NodeAlreadyMappedToRoot
      END IF ! ElemNodeDepoMap(iRank)%firstNode
    END IF ! GlobalNbElemID.NE.myrank
    ! ========================================================================================================
  END DO iVertexIndLoop2 ! iVertexInd = iFirstVertexInd,LastVertexInd
END DO ! iCNElem = 1,nComputeNodeTotalElems

! 4.) Get the number of send nodes for each communication partner: Size of each message for each process for deposition
! Initialize
nSendUniqueNodesNonSymDepo         = 0
nRecvUniqueNodesNonSymDepo(myrank) = 0 ! Nullify myself
! Allocate container for sending nodes to each communication partner
ALLOCATE(SurfNodeMappingSend(1:nSurfNodeSendExchangeProcs))
! Loop over each communication partner
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! Store the number of vertices that will be sent to the i-th process
  SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes =  ElemNodeDepoMap(iProc)%nNodes
  ! Store the number of vertices that will be sent to SurfNodeSendDepoRankToGlobalRank(iProc)
  nSendUniqueNodesNonSymDepo(SurfNodeSendDepoRankToGlobalRank(iProc)) = SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes
END DO

! 5.) MPI send/receive the number of deposition nodes
! Open receive buffer for non-symmetric exchange identification
DO iProc = 0,nProcessors_Global-1
  ! Ignore myself
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_IRECV( nRecvUniqueNodesNonSymDepo(iProc)  &
    , 1                                              &
    , MPI_INTEGER                                    &
    , iProc                                          &
    , 20002                                          &
    , MPI_COMM_PICLAS                                &
    , RecvRequestNonSymDepo(iProc)                   &
    , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! Send each communication partner the number of nodes that can be reached by deposition
DO iProc = 0,nProcessors_Global-1
  ! Ignore myself
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_ISEND( nSendUniqueNodesNonSymDepo(iProc) &
    , 1                                             &
    , MPI_INTEGER                                   &
    , iProc                                         &
    , 20002                                         &
    , MPI_COMM_PICLAS                               &
    , SendRequestNonSymDepo(iProc)                  &
    , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_WAIT(RecvRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
  CALL MPI_WAIT(SendRequestNonSymDepo(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! 6.) From the received messages, determine the message size that is sent from each communication partner.
! Count the number of communication partners that have sent a vertex ID count greater than zero
nSurfNodeRecvExchangeProcs = COUNT(nRecvUniqueNodesNonSymDepo.GT.0)
! Allocate container for receiving nodes from each communication partner
ALLOCATE(SurfNodeMappingRecv(1:nSurfNodeRecvExchangeProcs))
! Allocate and initialize mapping from i-th communication partner to global rank
ALLOCATE(SurfNodeRecvDepoRankToGlobalRank(1:nSurfNodeRecvExchangeProcs))
SurfNodeRecvDepoRankToGlobalRank = 0
! Nullify the iteration counter
nSurfNodeRecvExchangeProcs = 0
! Loop over the global number of procsses
DO iRank= 0, nProcessors_Global-1
  ! Ignore myself
  IF (iRank.EQ.myRank) CYCLE
  ! Only consider communication partners that have sent nodes to me
  IF (nRecvUniqueNodesNonSymDepo(iRank).GT.0) THEN
    ! Increment the exchange proc index
    nSurfNodeRecvExchangeProcs = nSurfNodeRecvExchangeProcs + 1
    ! Create mapping from i-th communication partner to global rank: Store global rank of i-th receive rank
    SurfNodeRecvDepoRankToGlobalRank(nSurfNodeRecvExchangeProcs) = iRank
    ! Store number of nodes for the i-th receive rank
    SurfNodeMappingRecv(nSurfNodeRecvExchangeProcs)%nRecvUniqueSurfNodes = nRecvUniqueNodesNonSymDepo(iRank)
  END IF
END DO

! 7.) MPI send/receive the vertex IDs of deposition nodes
! Open receive buffer with the number of nodes received from each process
ALLOCATE(SurfRecvRequest(1:nSurfNodeRecvExchangeProcs))
! Loop over each communication partner
DO iProc = 1, nSurfNodeRecvExchangeProcs
  ! Allocate containers for receiving the FEM vertex IDs and surface charge
  ALLOCATE(SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID(1:SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes))
  ALLOCATE(SurfNodeMappingRecv(iProc)%RecvSurfNodeSource(     1:SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes))
  ALLOCATE(SurfNodeMappingRecv(iProc)%RecvSurfNodeArea(       1:SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes))
  ! Open receive buffer
  CALL MPI_IRECV( SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID &
    , SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes                   &
    , MPI_INTEGER                                                       &
    , SurfNodeRecvDepoRankToGlobalRank(iProc)                           &
    , 6662                                                              &
    , MPI_COMM_PICLAS                                                   &
    , SurfRecvRequest(iProc)                                            &
    , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! Open send buffer
ALLOCATE(SurfSendRequest(1:nSurfNodeSendExchangeProcs))
! Loop over each communication partner
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  ! Allocate containers for sending the FEM vertex IDs and surface charge
  ALLOCATE(SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(1:SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes))
  ALLOCATE(SurfNodeMappingSend(iProc)%SendSurfNodeSource(     1:SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes))
  ALLOCATE(SurfNodeMappingSend(iProc)%SendSurfNodeArea(       1:SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes))
  SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID = -1
  SurfNodeMappingSend(iProc)%SendSurfNodeSource      = 0.
  SurfNodeMappingSend(iProc)%SendSurfNodeArea        = 0.
  ! Nullify iterator
  SendNodeCount = 0

  ! First loop: Traverse the list and populate SurfNodeMappingSend
  node => ElemNodeDepoMap(iProc)%first
  ! Sanity check: This can happen when MPIRoot is forced as communication partner even though no nodes for this process are found
  IF(.NOT.ASSOCIATED(node)) CALL abort(__STAMP__,' Error in InitDepoSurfNodesMPI: node pointer not associated')
  ! Loop until the end of the list is encountered
  DO WHILE (ASSOCIATED(node))
    ! Increment counter
    SendNodeCount = SendNodeCount + 1
    ! Store NodeID for sending
    SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(SendNodeCount) = node%NodeID
    ! Next link
    node => node%next
  END DO

  ! Deallocate the list
  CALL DeallocateNodeList(ElemNodeDepoMap(iProc)%first)
  ! Nullify the pointer
  NULLIFY(ElemNodeDepoMap(iProc)%first)
  ! Nullify the number of nodes
  ElemNodeDepoMap(iProc)%nNodes = 0

  CALL MPI_ISEND( SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID &
    , SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes                   &
    , MPI_INTEGER                                                       &
    , SurfNodeSendDepoRankToGlobalRank(iProc)                           &
    , 6662                                                              &
    , MPI_COMM_PICLAS                                                   &
    , SurfSendRequest(iProc)                                            &
    , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! Finish send
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  CALL MPI_WAIT(SurfSendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! Finish receive
DO iProc = 1, nSurfNodeRecvExchangeProcs
  CALL MPI_WAIT(SurfRecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in InitDepoSurfNodesMPI, IERROR=', IERROR)
END DO

! OPTIMIZE:Check if the received FEMVertexIDs are actually on the receiving process
! Skip MPIRoot, as this process receves all nodes from the other processes for .h5 output
IF(myrank.NE.0) DEALLOCATE(IsDepoSurfNode)
END SUBROUTINE InitDepoSurfNodesMPI


RECURSIVE SUBROUTINE DeallocateNodeList(node)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(NodeDepoMapping), POINTER, INTENT(INOUT) :: node
!===================================================================================================================================

IF (ASSOCIATED(node)) THEN
  CALL DeallocateNodeList(node%next)
  DEALLOCATE(node)
END IF

END SUBROUTINE DeallocateNodeList


!===================================================================================================================================
!> Exchange the node source container between MPI processes
!===================================================================================================================================
SUBROUTINE ExchangeNodeSource(SourceDim,doCalculateCurrentDensity)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars ,ONLY: NodeSource
USE MOD_PICDepo_Vars ,ONLY: nNodeSendExchangeProcs,NodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars ,ONLY: nNodeRecvExchangeProcs
USE MOD_PICDepo_Vars ,ONLY: NodeRecvDepoRankToGlobalRank
USE MOD_PICDepo_Vars ,ONLY: NodeMappingSend,NodeMappingRecv,nNodeSendExchangeProcs,NodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars ,ONLY: nNodeRecvExchangeProcs,NodeRecvDepoRankToGlobalRank
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN) :: SourceDim
LOGICAL, INTENT(IN) :: doCalculateCurrentDensity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iProc
TYPE(MPI_Request) :: RecvRequest(1:nNodeRecvExchangeProcs),SendRequest(1:nNodeSendExchangeProcs)
INTEGER           :: iNode
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)  :: CounterStart,CounterEnd
REAL(KIND=dp)     :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! 1.1) Receive charge density
DO iProc = 1, nNodeRecvExchangeProcs
  ! Open receive buffer
  CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeSourceCharge(:) &
            , NodeMappingRecv(iProc)%nRecvUniqueNodes            &
            , MPI_DOUBLE_PRECISION                               &
            , NodeRecvDepoRankToGlobalRank(iProc)                &
            , 666                                                &
            , MPI_COMM_PICLAS                                    &
            , RecvRequest(iProc)                                 &
            , IERROR)
END DO

! 1.2) Send charge density
DO iProc = 1, nNodeSendExchangeProcs
  ! Send message (non-blocking)
  DO iNode = 1, NodeMappingSend(iProc)%nSendUniqueNodes
    NodeMappingSend(iProc)%SendNodeSourceCharge(iNode) = NodeSource(4,NodeMappingSend(iProc)%SendNodeUniqueGlobalID(iNode))
  END DO
  CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeSourceCharge(:) &
                , NodeMappingSend(iProc)%nSendUniqueNodes        &
                , MPI_DOUBLE_PRECISION                           &
                , NodeSendDepoRankToGlobalRank(iProc)            &
                , 666                                            &
                , MPI_COMM_PICLAS                                &
                , SendRequest(iProc)                             &
                , IERROR)
END DO

! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc = 1, nNodeSendExchangeProcs
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
DO iProc = 1, nNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 2) Send/Receive current density
IF(doCalculateCurrentDensity)THEN
  DO iProc = 1, nNodeRecvExchangeProcs
    ! Open receive buffer
    CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeSourceCurrent(1:3,:) &
        , 3*NodeMappingRecv(iProc)%nRecvUniqueNodes                     &
        , MPI_DOUBLE_PRECISION                                          &
        , NodeRecvDepoRankToGlobalRank(iProc)                           &
        , 666                                                           &
        , MPI_COMM_PICLAS                                               &
        , RecvRequest(iProc)                                            &
        , IERROR)
  END DO

  DO iProc = 1, nNodeSendExchangeProcs
    ! Send message (non-blocking)
    DO iNode = 1, NodeMappingSend(iProc)%nSendUniqueNodes
      NodeMappingSend(iProc)%SendNodeSourceCurrent(1:3,iNode) = NodeSource(1:3,NodeMappingSend(iProc)%SendNodeUniqueGlobalID(iNode))
    END DO
    CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeSourceCurrent(1:3,:) &
        , 3*NodeMappingSend(iProc)%nSendUniqueNodes                     &
        , MPI_DOUBLE_PRECISION                                          &
        , NodeSendDepoRankToGlobalRank(iProc)                           &
        , 666                                                           &
        , MPI_COMM_PICLAS                                               &
        , SendRequest(iProc)                                            &
        , IERROR)
  END DO

  ! Finish communication
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
  DO iProc = 1, nNodeSendExchangeProcs
    CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO
  DO iProc = 1, nNodeRecvExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

  ! 3) Extract messages
  DO iProc = 1, nNodeRecvExchangeProcs
    DO iNode = 1, NodeMappingRecv(iProc)%nRecvUniqueNodes
      ASSOCIATE( NS => NodeSource(SourceDim:4,NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(iNode)))
        NS = NS + (/NodeMappingRecv(iProc)%RecvNodeSourceCurrent(1:3,iNode), NodeMappingRecv(iProc)%RecvNodeSourceCharge(iNode)/)
      END ASSOCIATE
    END DO
  END DO
ELSE ! Only the charge density is communicated
  DO iProc = 1, nNodeRecvExchangeProcs
    DO iNode = 1, NodeMappingRecv(iProc)%nRecvUniqueNodes
      ASSOCIATE( NS => NodeSource(4,NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(iNode)))
        NS = NS + NodeMappingRecv(iProc)%RecvNodeSourceCharge(iNode)
      END ASSOCIATE
    END DO
  END DO
END IF ! doCalculateCurrentDensity
END SUBROUTINE ExchangeNodeSource


!===================================================================================================================================
!> Exchange the node source container between MPI processes (either during load balance or hdf5 output) and nullify the local charge
!> container NodeSourceExtMPI. Updates the node charge container NodeSourceExt at MPI interfaces.
!===================================================================================================================================
SUBROUTINE ExchangeNodeSourceExtMPI()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeMappingRecv,NodeMappingSend,NodeSourceExtMPI
USE MOD_PICDepo_Vars       ,ONLY: nDepoNodesTotal,nNodeSendExchangeProcs,NodeSendDepoRankToGlobalRank,DepoNodetoGlobalNode
USE MOD_PICDepo_Vars       ,ONLY: nNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: NodeRecvDepoRankToGlobalRank
#endif  /*USE_MPI*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                        :: iProc
TYPE(MPI_Request)              :: RecvRequest(1:nNodeRecvExchangeProcs),SendRequest(1:nNodeSendExchangeProcs)
#endif /*USE_MPI*/
INTEGER                        :: globalNode, iNode
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)               :: CounterStart,CounterEnd
REAL(KIND=dp)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! 1) Receive charge density
DO iProc = 1, nNodeRecvExchangeProcs
  ! Open receive buffer
  CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeSourceExt(:) &
      , NodeMappingRecv(iProc)%nRecvUniqueNodes               &
      , MPI_DOUBLE_PRECISION                                  &
      , NodeRecvDepoRankToGlobalRank(iProc)                   &
      , 666                                                   &
      , MPI_COMM_PICLAS                                       &
      , RecvRequest(iProc)                                    &
      , IERROR)
END DO

DO iProc = 1, nNodeSendExchangeProcs
  ! Send message (non-blocking)
  DO iNode = 1, NodeMappingSend(iProc)%nSendUniqueNodes
    NodeMappingSend(iProc)%SendNodeSourceExt(iNode) = NodeSourceExtMPI(NodeMappingSend(iProc)%SendNodeUniqueGlobalID(iNode))
  END DO
  CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeSourceExt(:) &
      , NodeMappingSend(iProc)%nSendUniqueNodes               &
      , MPI_DOUBLE_PRECISION                                  &
      , NodeSendDepoRankToGlobalRank(iProc)                   &
      , 666                                                   &
      , MPI_COMM_PICLAS                                       &
      , SendRequest(iProc)                                    &
      , IERROR)
END DO
! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc = 1, nNodeSendExchangeProcs
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
DO iProc = 1, nNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
DO iProc = 1, nNodeRecvExchangeProcs
  DO iNode = 1, NodeMappingRecv(iProc)%nRecvUniqueNodes
    ASSOCIATE( NS => NodeSourceExtMPI(NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(iNode)))
      NS = NS + NodeMappingRecv(iProc)%RecvNodeSourceExt(iNode)
    END ASSOCIATE
  END DO
END DO

! Add NodeSourceExtMPI values of the last boundary interaction
DO iNode = 1, nDepoNodesTotal
  globalNode = DepoNodetoGlobalNode(iNode)
  NodeSourceExt(globalNode) = NodeSourceExt(globalNode) + NodeSourceExtMPI(globalNode)
END DO
! Reset local surface charge
NodeSourceExtMPI = 0.
END SUBROUTINE ExchangeNodeSourceExtMPI


!===================================================================================================================================
!> Collect the surface area contributions SurfNodeArea(iDepoSurfNodeID) of all processes on the MPIRoot
!> MPIRoot process: only receives data
!> non-MPIRoot processes: only send data
!===================================================================================================================================
SUBROUTINE CollectSurfNodeAreaOnMPIRoot()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeArea
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeMappingRecv,SurfNodeMappingSend
USE MOD_PICDepo_Vars       ,ONLY: nDepoSurfNodesTotal,nSurfNodeSendExchangeProcs,SurfNodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars       ,ONLY: FEMVertexID2DepoSurfNodeID
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeRecvDepoRankToGlobalRank
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
TYPE(MPI_Request)              :: RecvRequest(1:nSurfNodeRecvExchangeProcs),SendRequest(1:nSurfNodeSendExchangeProcs)
INTEGER                        :: iNode, iDepoSurfNodeID, FEMVertexID
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)               :: CounterStart,CounterEnd
REAL(KIND=dp)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! 1) Receive surface node area data
! Only MPIRoot receives data from all other processes
IF (MPIRoot) THEN
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    ! Open receive buffer
    CALL MPI_IRECV( SurfNodeMappingRecv(iProc)%RecvSurfNodeArea(:)      &
        , SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes               &
        , MPI_DOUBLE_PRECISION                                          &
        , SurfNodeRecvDepoRankToGlobalRank(iProc)                       &
        , 16662                                                         &
        , MPI_COMM_PICLAS                                               &
        , RecvRequest(iProc)                                            &
        , IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in CollectSurfNodeAreaOnMPIRoot', IERROR)
  END DO
END IF ! MPIRoot

! Loop over all communication partners and skip all non-MPIRoot processes
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! Skip non-MPIRoot processes
  IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  ! This concerns processes, which do not contribute to the surface charge deposition
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  ! Send message (non-blocking)
  DO iNode = 1, SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes
    ! Get FEMVertexID from mapping
    FEMVertexID = SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(iNode)
    ! Get surface deposition node index
    iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
    ! Store in send array
    SurfNodeMappingSend(iProc)%SendSurfNodeArea(iNode) = SurfNodeArea(iDepoSurfNodeID)
  END DO
  CALL MPI_ISEND( SurfNodeMappingSend(iProc)%SendSurfNodeArea(:)      &
      , SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes               &
      , MPI_DOUBLE_PRECISION                                          &
      , SurfNodeSendDepoRankToGlobalRank(iProc)                       &
      , 16662                                                         &
      , MPI_COMM_PICLAS                                               &
      , SendRequest(iProc)                                            &
      , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in CollectSurfNodeAreaOnMPIRoot', IERROR)
END DO

! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
! Loop over all communication partners and skip all non-MPIRoot processes
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! Skip non-MPIRoot processes
  IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  ! This concerns processes, which do not contribute to the surface charge deposition
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in CollectSurfNodeAreaOnMPIRoot', IERROR)
END DO
! Only MPIRoot receives data from all other processes
IF (MPIRoot) THEN
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in CollectSurfNodeAreaOnMPIRoot', IERROR)
  END DO
END IF ! MPIRoot
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
! Only MPIRoot extracts data
IF (MPIRoot) THEN
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    DO iNode = 1, SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes
      ! Get FEMVertexID from mapping
      FEMVertexID = SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID(iNode)
      IF(FEMVertexID.LE.0) CALL abort(__STAMP__,'ERROR: Invalid FEMVertexID <= 0',FEMVertexID)
      ! Get surface deposition node index
      iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
      IF((iDepoSurfNodeID.LE.0).OR.(iDepoSurfNodeID.GT.nDepoSurfNodesTotal)) CALL abort(__STAMP__,'ERROR: Invalid iDepoSurfNodeID <= 0 or > nDepoSurfNodesTotal',iDepoSurfNodeID)
      ! Unpack in recv array
      SurfNodeArea(iDepoSurfNodeID) = SurfNodeArea(iDepoSurfNodeID) + SurfNodeMappingRecv(iProc)%RecvSurfNodeArea(iNode)
    END DO
  END DO
END IF ! MPIRoot
END SUBROUTINE CollectSurfNodeAreaOnMPIRoot


!===================================================================================================================================
!> Exchange the node source container between MPI processes (either during load balance or hdf5 output) and nullify the local charge
!> container SurfNodeSourceMPI. Updates the node charge container SurfNodeSource at MPI interfaces.
!===================================================================================================================================
SUBROUTINE ExchangeSurfNodeSourceMPI()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeSource
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeMappingRecv,SurfNodeMappingSend,SurfNodeSourceMPI
USE MOD_PICDepo_Vars       ,ONLY: nDepoSurfNodesTotal,nSurfNodeSendExchangeProcs,SurfNodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars       ,ONLY: FEMVertexID2DepoSurfNodeID
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeRecvDepoRankToGlobalRank
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
TYPE(MPI_Request)              :: RecvRequest(1:nSurfNodeRecvExchangeProcs),SendRequest(1:nSurfNodeSendExchangeProcs)
INTEGER                        :: iNode, iDepoSurfNodeID, FEMVertexID
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)               :: CounterStart,CounterEnd
REAL(KIND=dp)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! 1) Receive surface node source data
DO iProc = 1, nSurfNodeRecvExchangeProcs
  ! Open receive buffer
  CALL MPI_IRECV( SurfNodeMappingRecv(iProc)%RecvSurfNodeSource(:)    &
      , SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes               &
      , MPI_DOUBLE_PRECISION                                          &
      , SurfNodeRecvDepoRankToGlobalRank(iProc)                       &
      , 6662                                                          &
      , MPI_COMM_PICLAS                                               &
      , RecvRequest(iProc)                                            &
      , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in ExchangeSurfNodeSourceMPI', IERROR)
END DO

DO iProc = 1, nSurfNodeSendExchangeProcs
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  ! Send message (non-blocking)
  DO iNode = 1, SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes
    ! Get FEMVertexID from mapping
    FEMVertexID = SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(iNode)
    ! Get surface deposition node index
    iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
    ! Store in send array
    SurfNodeMappingSend(iProc)%SendSurfNodeSource(iNode) = SurfNodeSourceMPI(iDepoSurfNodeID)
  END DO
  CALL MPI_ISEND( SurfNodeMappingSend(iProc)%SendSurfNodeSource(:)    &
      , SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes               &
      , MPI_DOUBLE_PRECISION                                          &
      , SurfNodeSendDepoRankToGlobalRank(iProc)                       &
      , 6662                                                          &
      , MPI_COMM_PICLAS                                               &
      , SendRequest(iProc)                                            &
      , IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in ExchangeSurfNodeSourceMPI', IERROR)
END DO

! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
DO iProc = 1, nSurfNodeSendExchangeProcs
  ! OPTIMIZE: All processes communicate with MPIRoot (rank 0) for output to .h5, which is solely done by MPIRoot
  ! Skip MPIRoot, which can happen as it is forced as communication partner even though no nodes for this process are found
  IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in ExchangeSurfNodeSourceMPI', IERROR)
END DO
DO iProc = 1, nSurfNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error in ExchangeSurfNodeSourceMPI', IERROR)
END DO
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
DO iProc = 1, nSurfNodeRecvExchangeProcs
  DO iNode = 1, SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes
    ! Get FEMVertexID from mapping
    FEMVertexID = SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID(iNode)
    IF(FEMVertexID.LE.0) CALL abort(__STAMP__,'ERROR: Invalid FEMVertexID <= 0',FEMVertexID)
    ! Get surface deposition node index
    iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
    IF((iDepoSurfNodeID.LE.0).OR.(iDepoSurfNodeID.GT.nDepoSurfNodesTotal)) CALL abort(__STAMP__,'ERROR: Invalid iDepoSurfNodeID <= 0 or > nDepoSurfNodesTotal',iDepoSurfNodeID)
    ! Unpack in recv array
    ASSOCIATE( NS => SurfNodeSourceMPI(iDepoSurfNodeID) )
      NS = NS + SurfNodeMappingRecv(iProc)%RecvSurfNodeSource(iNode)
    END ASSOCIATE
  END DO
END DO

! Add SurfNodeSourceMPI values of the last boundary interaction
DO iDepoSurfNodeID = 1, nDepoSurfNodesTotal
  ! Add contribution
  SurfNodeSource(iDepoSurfNodeID) = SurfNodeSource(iDepoSurfNodeID) + SurfNodeSourceMPI(iDepoSurfNodeID)
END DO
! Reset local surface charge
SurfNodeSourceMPI = 0.
END SUBROUTINE ExchangeSurfNodeSourceMPI


!===================================================================================================================================
!> Initialize the the SurfNodeSource container on all processes except MPIRoot, which distribtues the data to all others
!> MPIRoot sends all surface charge deposition sending processes the NodeSource data during load balancing and when restarting the
!> simulation because the MPIRoot process has the complete global information.
!> ATTENTION: Do not be confused because Send/Receive containers are used in reverse in this routine
!===================================================================================================================================
SUBROUTINE LBReverseExchangeSurfNodeSource()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeSource
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeMappingRecv,SurfNodeMappingSend
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeSendExchangeProcs,SurfNodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars       ,ONLY: FEMVertexID2DepoSurfNodeID
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeRecvDepoRankToGlobalRank
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
TYPE(MPI_Request)              :: RecvRequest(1:nSurfNodeSendExchangeProcs),SendRequest(1:nSurfNodeRecvExchangeProcs) !> ATTENTION:
! Send/Receive containers are used in reverse here
INTEGER                        :: iNode, iDepoSurfNodeID, FEMVertexID
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)               :: CounterStart,CounterEnd
REAL(KIND=dp)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! OPTIMIZE: check if the SurfNodeSource is sent to processes, which do not have any fem vertices for surface deposition and skip them
! 1) Receive surface charge density
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only open buffer with MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    ! MPIRoot will not send anything, if the process has zero send nodes
    IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
    ! Open receive buffer
    CALL MPI_IRECV( SurfNodeMappingSend(iProc)%SendSurfNodeSource(:)    &
        , SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes               &
        , MPI_DOUBLE_PRECISION                                          &
        , SurfNodeSendDepoRankToGlobalRank(iProc)                       &
        , 666                                                           &
        , MPI_COMM_PICLAS                                               &
        , RecvRequest(iProc)                                            &
        , IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' LBReverseExchangeSurfNodeSource: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! .NOT.MPIRoot

! Only MPIRoot sends data to all other processes that normally send data to MPIRoot
IF (MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    ! Send message (non-blocking)
    DO iNode = 1, SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes
      ! Get FEMVertexID from mapping
      FEMVertexID = SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID(iNode)
      ! Get surface deposition node index
      iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
      ! Store in send array
      SurfNodeMappingRecv(iProc)%RecvSurfNodeSource(iNode) = SurfNodeSource(iDepoSurfNodeID)
    END DO

    CALL MPI_ISEND( SurfNodeMappingRecv(iProc)%RecvSurfNodeSource(:)    &
        , SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes               &
        , MPI_DOUBLE_PRECISION                                          &
        , SurfNodeRecvDepoRankToGlobalRank(iProc)                       &
        , 666                                                           &
        , MPI_COMM_PICLAS                                               &
        , SendRequest(iProc)                                            &
        , IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' LBReverseExchangeSurfNodeSource: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! MPIRoot

! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
! Only MPIRoot sends data to all other processes that normally send data to MPIRoot
IF (MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' LBReverseExchangeSurfNodeSource: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! MPIRoot
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only open buffer with MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    ! MPIRoot will not send anything, if the process has zero send nodes
    IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' LBReverseExchangeSurfNodeSource: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! .NOT.MPIRoot
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only extract buffer from MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    DO iNode = 1, SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes
      ! Get FEMVertexID from mapping
      FEMVertexID = SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(iNode)
      ! Get surface deposition node index
      iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
      if(FEMVertexID.LE.0) CALL abort(__STAMP__,' FEMVertexID <= 0',FEMVertexID)
      ! Unpack in recv array
      SurfNodeSource(iDepoSurfNodeID) = SurfNodeMappingSend(iProc)%SendSurfNodeSource(iNode)
    END DO
  END DO
END IF ! .NOT.MPIRoot

END SUBROUTINE LBReverseExchangeSurfNodeSource


!===================================================================================================================================
!> Initialize the the SurfNodeArea container on all processes except MPIRoot, which distribtues the data to all others
!> MPIRoot sends all surface charge deposition sending processes the SurfNodeArea(iDepoSurfNodeID) data during load balancing and when
!> restarting the simulation because the MPIRoot process has the complete global information.
!> ATTENTION: Do not be confused because Send/Receive containers are used in reverse in this routine
!===================================================================================================================================
SUBROUTINE ReverseExchangeSurfNodeArea()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeArea
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeMappingRecv,SurfNodeMappingSend
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeSendExchangeProcs,SurfNodeSendDepoRankToGlobalRank
USE MOD_PICDepo_Vars       ,ONLY: FEMVertexID2DepoSurfNodeID
USE MOD_PICDepo_Vars       ,ONLY: nSurfNodeRecvExchangeProcs
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeRecvDepoRankToGlobalRank
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars  ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
TYPE(MPI_Request)              :: RecvRequest(1:nSurfNodeSendExchangeProcs),SendRequest(1:nSurfNodeRecvExchangeProcs) !> ATTENTION:
! Send/Receive containers are used in reverse here
INTEGER                        :: iNode, iDepoSurfNodeID, FEMVertexID
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)               :: CounterStart,CounterEnd
REAL(KIND=dp)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
! OPTIMIZE: check if the SurfNodeArea is sent to processes, which do not have any fem vertices for surface deposition and skip them
! 1) Receive surface charge density
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only open buffer with MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    ! MPIRoot will not send anything, if the process has zero send nodes
    IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
    ! Open receive buffer
    CALL MPI_IRECV( SurfNodeMappingSend(iProc)%SendSurfNodeArea(:)      &
        , SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes               &
        , MPI_DOUBLE_PRECISION                                          &
        , SurfNodeSendDepoRankToGlobalRank(iProc)                       &
        , 89666                                                         &
        , MPI_COMM_PICLAS                                               &
        , RecvRequest(iProc)                                            &
        , IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' ReverseExchangeSurfNodeArea: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! .NOT.MPIRoot

! Only MPIRoot sends data to all other processes that normally send data to MPIRoot
IF (MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    ! Send message (non-blocking)
    DO iNode = 1, SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes
      ! Get FEMVertexID from mapping
      FEMVertexID = SurfNodeMappingRecv(iProc)%RecvSurfNodeFEMVertexID(iNode)
      ! Get surface deposition node index
      iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
      ! Store in send array
      SurfNodeMappingRecv(iProc)%RecvSurfNodeArea(iNode) = SurfNodeArea(iDepoSurfNodeID)
    END DO

    CALL MPI_ISEND( SurfNodeMappingRecv(iProc)%RecvSurfNodeArea(:)      &
        , SurfNodeMappingRecv(iProc)%nRecvUniqueSurfNodes               &
        , MPI_DOUBLE_PRECISION                                          &
        , SurfNodeRecvDepoRankToGlobalRank(iProc)                       &
        , 89666                                                         &
        , MPI_COMM_PICLAS                                               &
        , SendRequest(iProc)                                            &
        , IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' ReverseExchangeSurfNodeArea: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! MPIRoot

! Finish communication
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
! Only MPIRoot sends data to all other processes that normally send data to MPIRoot
IF (MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeRecvExchangeProcs
    CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' ReverseExchangeSurfNodeArea: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! MPIRoot
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  ! ATTENTION: Send/Receive containers are used in reverse in this routine
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only open buffer with MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    ! MPIRoot will not send anything, if the process has zero send nodes
    IF(SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes.EQ.0) CYCLE
    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' ReverseExchangeSurfNodeArea: MPI Communication error. IERROR=', IERROR)
  END DO
END IF ! .NOT.MPIRoot
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimePart(6)  = MPIW8TimePart(6) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountPart(6) = MPIW8CountPart(6) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

! 3) Extract messages
! Skip MPIRoot because this process only sends
IF (.NOT.MPIRoot) THEN
  DO iProc = 1, nSurfNodeSendExchangeProcs
    ! Only extract buffer from MPIRoot
    IF(SurfNodeSendDepoRankToGlobalRank(iProc).NE.0) CYCLE
    DO iNode = 1, SurfNodeMappingSend(iProc)%nSendUniqueSurfNodes
      ! Get FEMVertexID from mapping
      FEMVertexID = SurfNodeMappingSend(iProc)%SendSurfNodeFEMVertexID(iNode)
      ! Get surface deposition node index
      iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
      if(FEMVertexID.LE.0) CALL abort(__STAMP__,' FEMVertexID <= 0',FEMVertexID)
      ! Unpack in recv array
      SurfNodeArea(iDepoSurfNodeID) = SurfNodeMappingSend(iProc)%SendSurfNodeArea(iNode)
    END DO
  END DO
END IF ! .NOT.MPIRoot

END SUBROUTINE ReverseExchangeSurfNodeArea
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*USE_MPI*/
END MODULE MOD_PICDepo_MPI