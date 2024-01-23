!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_MPI_Boundary_Sampling
!===================================================================================================================================
! module for MPI communication of particle surface sampling
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE InitSurfCommunication
  MODULE PROCEDURE InitSurfCommunication
END INTERFACE

INTERFACE ExchangeSurfData
  MODULE PROCEDURE ExchangeSurfData
END INTERFACE

INTERFACE FinalizeSurfCommunication
  MODULE PROCEDURE FinalizeSurfCommunication
END INTERFACE

PUBLIC :: InitSurfCommunication
PUBLIC :: ExchangeSurfData
PUBLIC :: FinalizeSurfCommunication
!===================================================================================================================================

CONTAINS


SUBROUTINE InitSurfCommunication()
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeProcessors
USE MOD_MPI_Shared_Vars         ,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_MPI_Shared_Vars         ,ONLY: MPIRankSharedLeader,MPIRankSurfLeader
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank,nSurfLeaders
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides,offsetComputeNodeSurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode,SurfSampSize,nSurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalSurfSides, nGlobalOutputSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfOutputSides,offsetComputeNodeSurfOutputSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared, SideInfo_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeInnerBCs
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                       :: iProc,color
INTEGER                       :: leadersGroup,LeaderID,surfGroup
INTEGER                       :: iSide
INTEGER                       :: sendbuf,recvbuf
INTEGER                       :: nSendSurfSidesTmp(0:nLeaderGroupProcs-1)
INTEGER                       :: nRecvSurfSidesTmp(0:nLeaderGroupProcs-1)
!INTEGER                       :: nSurfSidesLeader(1:2,0:nLeaderGroupProcs-1)
INTEGER                       :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                       :: SendSurfGlobalID(0:nLeaderGroupProcs-1,1:nComputeNodeSurfTotalSides)
INTEGER                       :: SampSizeAllocate
INTEGER                       :: NbGlobalElemID, GlobalSideID, NbGlobalSideID, NbElemRank, NbLeaderID, GlobalElemID, ElemRank
INTEGER                       :: TestCounter(2),iCNinnerBC
INTEGER                       :: SwitchGlobalSideID(1:3,1:SUM(nComputeNodeInnerBCs)),nSideTmp
INTEGER                       :: allocstat
!===================================================================================================================================

nRecvSurfSidesTmp = 0

!--- Open receive buffer (number of sampling surfaces in other node's halo region)
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_IRECV( nRecvSurfSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

!--- count all surf sides per other compute-node which get sampling data from current leader
nSendSurfSidesTmp = 0

TestCounter=0
SwitchGlobalSideID = -1
!--- Loop 1 of 2: Non-HALO sides
! store the GlobalSideID of each sampling side on a different processor (node leader) to send each leader the number and IDs of all
! sampling sides. Special treatment is required for inner BCs as these are only written to .h5 for the smaller GlobalSideID and
! therefor, the processor that samples on the larger side ID must send this information to the leader of the smaller side in order
! to store this information. The receiving leader is sent the smaller GlobalSideID but the larger ID is kept for local sampling and
! is stored here in SwitchGlobalSideID, which is later applied after the message has been sent. In the first part, the local
! (non-HALO) surf sides are considered and in a second step the HALO sides.
DO iSide = 1,nComputeNodeSurfSides
  ! count surf sides per compute node
  LeaderID = SurfSide2GlobalSide(SURF_LEADER,iSide)
  nSendSurfSidesTmp(LeaderID) = nSendSurfSidesTmp(LeaderID) + 1
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  SendSurfGlobalID(LeaderID,nSendSurfSidesTmp(LeaderID)) = GlobalSideID
  ! Check if the side has a neighbour side and is thus an inner BC
  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
    ! Only add sides that are NOT part of the halo region
    NbGlobalSideID = SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)
    ! Skip sides with the smaller global side index as those are on the receiving end
    IF(GlobalSideID.LT.NbGlobalSideID) CYCLE
    NbGlobalElemID = SideInfo_Shared(SIDE_ELEMID,NbGlobalSideID)
    NbElemRank = ElemInfo_Shared(ELEM_RANK,NbGlobalElemID)
    NbLeaderID = INT(NbElemRank/nComputeNodeProcessors)
    IF(NbLeaderID.NE.LeaderID) THEN
      nSendSurfSidesTmp(NbLeaderID) = nSendSurfSidesTmp(NbLeaderID) + 1
      SendSurfGlobalID(NbLeaderID,nSendSurfSidesTmp(NbLeaderID)) = NbGlobalSideID
      ! Store switcheroo information
      TestCounter(1) = TestCounter(1) + 1
      SwitchGlobalSideID(1:3,TestCounter(1) ) = (/NbLeaderID,GlobalSideID,nSendSurfSidesTmp(NbLeaderID)/)
    END IF
  END IF
END DO

!--- Loop 2 of 2: HALO sides
! Note that in this case three nodes might be involved
DO iSide = nComputeNodeSurfSides+1,nComputeNodeSurfTotalSides
  ! Check if the side has a neighbour side and is thus an inner BC
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
    NbGlobalSideID = SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)
    ! Skip sides with the smaller global side index as those are on the receiving end
    IF(GlobalSideID.GT.NbGlobalSideID) THEN
      NbGlobalElemID = SideInfo_Shared(SIDE_ELEMID,NbGlobalSideID)
      NbElemRank = ElemInfo_Shared(ELEM_RANK,NbGlobalElemID)
      NbLeaderID = INT(NbElemRank/nComputeNodeProcessors)
      nSendSurfSidesTmp(NbLeaderID) = nSendSurfSidesTmp(NbLeaderID) + 1
      SendSurfGlobalID(NbLeaderID,nSendSurfSidesTmp(NbLeaderID)) = NbGlobalSideID
      ! Store switcheroo information
      TestCounter(2) = TestCounter(2) + 1
      SwitchGlobalSideID(1:3,TestCounter(1)+TestCounter(2) ) = (/NbLeaderID,GlobalSideID,nSendSurfSidesTmp(NbLeaderID)/)
    ELSE
      GlobalElemID = SideInfo_Shared(SIDE_ELEMID,GlobalSideID)
      ElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
      LeaderID = INT(ElemRank/nComputeNodeProcessors)
      nSendSurfSidesTmp(LeaderID) = nSendSurfSidesTmp(LeaderID) + 1
      SendSurfGlobalID(LeaderID,nSendSurfSidesTmp(LeaderID)) = GlobalSideID
    END IF
  ELSE
    ! Count regular halo sampling on surf sides per compute node
    LeaderID = SurfSide2GlobalSide(SURF_LEADER,iSide)
    nSendSurfSidesTmp(LeaderID) = nSendSurfSidesTmp(LeaderID) + 1
    SendSurfGlobalID(LeaderID,nSendSurfSidesTmp(LeaderID)) = GlobalSideID
  END IF
END DO

! Sanity check for switcheroo:
! Compare the number of counted sides above with the number that was calculated in InitParticleBoundarySampling()
IF(SUM(TestCounter).NE.SUM(nComputeNodeInnerBCs))THEN
  IPWRITE(UNIT_StdOut,*) "Non-HALO: TestCounter(1),nComputeNodeInnerBCs(1) =", TestCounter(1),nComputeNodeInnerBCs(1)
  IPWRITE(UNIT_StdOut,*) "    HALO: TestCounter(2),nComputeNodeInnerBCs(2) =", TestCounter(2),nComputeNodeInnerBCs(2)
  CALL abort(__STAMP__,'InitSurfCommunication: The number of Inner BCs that are redirected to other procs do not match!')
END IF ! TestCounter.ne.nComputeNodeInnerBCs

!--- send all other leaders the number of sampling sides coming from current node
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_ISEND( nSendSurfSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

!--- Split communicator from MPI_COMM_LEADER_SHARED
color = MERGE(1201,MPI_UNDEFINED,SurfTotalSideOnNode)

! create new communicator between node leaders with surfaces. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_LEADERS_SHARED, color, MPI_INFO_NULL, MPI_COMM_LEADERS_SURF, IERROR)

! Do not participate in remainder of communication if no surf sides on node
IF (.NOT.SurfTotalSideOnNode) RETURN

! Find my rank on the shared communicator, comm size and proc name
CALL MPI_COMM_RANK(MPI_COMM_LEADERS_SURF, mySurfRank  , IERROR)
CALL MPI_COMM_SIZE(MPI_COMM_LEADERS_SURF, nSurfLeaders, IERROR)

! Map global rank number into shared rank number. Returns MPI_UNDEFINED if not on the same communicator
ALLOCATE(MPIRankSharedLeader(0:nLeaderGroupProcs-1))
ALLOCATE(MPIRankSurfLeader  (0:nLeaderGroupProcs-1))
DO iProc=0,nLeaderGroupProcs-1
  MPIRankSharedLeader(iProc) = iProc
END DO

! Get handles for each group
CALL MPI_COMM_GROUP(MPI_COMM_LEADERS_SHARED,leadersGroup,IERROR)
CALL MPI_COMM_GROUP(MPI_COMM_LEADERS_SURF  ,surfGroup   ,IERROR)

! Finally translate global rank to local rank
CALL MPI_GROUP_TRANSLATE_RANKS(leadersGroup,nLeaderGroupProcs,MPIRankSharedLeader,surfGroup,MPIRankSurfLeader,IERROR)
IF (mySurfRank.EQ.0) THEN
#if USE_LOADBALANCE
  IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
  WRITE(UNIT_stdOUt,'(A,I0,A)') ' Starting surface communication between ', nSurfLeaders, ' compute nodes...'
END IF

!--- Open receive buffer (mapping from message surface ID to global side ID)
ALLOCATE(SurfMapping(0:nSurfLeaders-1))

SurfMapping(:)%nRecvSurfSides = 0
SurfMapping(:)%nSendSurfSides = 0

DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  ! Save number of send and recv sides
  SurfMapping(MPIRankSurfLeader(iProc))%nRecvSurfSides = nRecvSurfSidesTmp(iProc)
  SurfMapping(MPIRankSurfLeader(iProc))%nSendSurfSides = nSendSurfSidesTmp(iProc)

  ! Only open recv buffer if we are expecting sides from this leader node
  IF (nRecvSurfSidesTmp(iProc).EQ.0) CYCLE

  ALLOCATE(SurfMapping(MPIRankSurfLeader(iProc))%RecvSurfGlobalID(1:nRecvSurfSidesTmp(iProc)))

  CALL MPI_IRECV( SurfMapping(MPIRankSurfLeader(iProc))%RecvSurfGlobalID                         &
                , nRecvSurfSidesTmp(iProc)                 &
                , MPI_INTEGER                                                 &
                , MPIRankSurfLeader(iProc)                                                      &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , RecvRequest(MPIRankSurfLeader(iProc))                                          &
                , IERROR)
END DO

DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  ! Only open send buffer if we are expecting sides from this leader node
  IF (nSendSurfSidesTmp(iProc).EQ.0) CYCLE

  ALLOCATE(SurfMapping(MPIRankSurfLeader(iProc))%SendSurfGlobalID(1:nSendSurfSidesTmp(iProc)))

  SurfMapping(MPIRankSurfLeader(iProc))%SendSurfGlobalID = SendSurfGlobalID(iProc,1:nSendSurfSidesTmp(iProc))

  CALL MPI_ISEND( SurfMapping(MPIRankSurfLeader(iProc))%SendSurfGlobalID                         &
                , nSendSurfSidesTmp(iProc)                 &
                , MPI_INTEGER                                                 &
                , MPIRankSurfLeader(iProc) &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , SendRequest(MPIRankSurfLeader(iProc))                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  IF (nSendSurfSidesTmp(iProc).NE.0) THEN
    CALL MPI_WAIT(SendRequest(MPIRankSurfLeader(iProc)),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF

  IF (nRecvSurfSidesTmp(iProc).NE.0) THEN
    CALL MPI_WAIT(RecvRequest(MPIRankSurfLeader(iProc)),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

!--- Finish switcheroo
DO iCNinnerBC = 1, SUM(nComputeNodeInnerBCs)
  NbLeaderID   = SwitchGlobalSideID(1,iCNinnerBC)
  GlobalSideID = SwitchGlobalSideID(2,iCNinnerBC)
  nSideTmp     = SwitchGlobalSideID(3,iCNinnerBC)
  ! Note that you have to cycle yourself, but that halo sides might also have been changed to a different proc (i.e.  to yourself)
  IF (NbLeaderID .EQ. myLeaderGroupRank) CYCLE
  SurfMapping(MPIRankSurfLeader(NbLeaderID))%SendSurfGlobalID(nSideTmp) = GlobalSideID
END DO ! iSide = 1, nComputeNodeInnerBCs

!--- Allocate send and recv buffer for each surf leader
ALLOCATE(SurfSendBuf(0:nSurfLeaders-1),STAT=allocstat)
IF(allocstat.ne.0) CALL abort(__STAMP__,'Could not allocate SurfSendBuf')
ALLOCATE(SurfRecvBuf(0:nSurfLeaders-1),STAT=allocstat)
IF(allocstat.ne.0) CALL abort(__STAMP__,'Could not allocate SurfRecvBuf')

DO iProc = 0,nSurfLeaders-1
  ! Get message size
  SampSizeAllocate = SurfSampSize
  ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number
  IF(CalcSurfaceImpact) SampSizeAllocate = SampSizeAllocate + 9*nSpecies

  ! Only allocate send buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nSendSurfSides.GT.0) THEN
    ALLOCATE(SurfSendBuf(iProc)%content(SampSizeAllocate*(nSurfSample**2)*SurfMapping(iProc)%nSendSurfSides),STAT=allocstat)
    IF(allocstat.ne.0) CALL abort(__STAMP__,'Could not allocate SurfSendBuf(iProc)%content')
    SurfSendBuf(iProc)%content = 0.
  END IF

  ! Only allocate recv buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nRecvSurfSides.GT.0) THEN
    ALLOCATE(SurfRecvBuf(iProc)%content(SampSizeAllocate*(nSurfSample**2)*SurfMapping(iProc)%nRecvSurfSides),STAT=allocstat)
    IF(allocstat.ne.0) CALL abort(__STAMP__,'Could not allocate SurfRecvBuf(iProc)%content')
    SurfRecvBuf(iProc)%content = 0.
  END IF
END DO ! iProc

!--- Save number of output sides per node (inner BCs are only included once here)
IF (nSurfLeaders.EQ.1) THEN
  offsetComputeNodeSurfOutputSide = 0
  nGlobalOutputSides           = nComputeNodeSurfOutputSides
ELSE
  sendbuf = nComputeNodeSurfOutputSides
  recvbuf = 0
  CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
  offsetComputeNodeSurfOutputSide = recvbuf
  ! last proc knows CN total number of BC elems
  sendbuf = offsetComputeNodeSurfOutputSide + nComputeNodeSurfOutputSides
  CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nSurfLeaders-1,MPI_COMM_LEADERS_SURF,iError)
  nGlobalOutputSides = sendbuf
END IF


!--- Save number of total surf sides
IF (nSurfLeaders.EQ.1) THEN
  offsetComputeNodeSurfSide = 0
  nGlobalSurfSides           = nComputeNodeSurfSides
ELSE
  sendbuf = nComputeNodeSurfSides
  recvbuf = 0
  CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
  offsetComputeNodeSurfSide = recvbuf
  ! last proc knows CN total number of BC elems
  sendbuf = offsetComputeNodeSurfSide + nComputeNodeSurfSides
  CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nSurfLeaders-1,MPI_COMM_LEADERS_SURF,iError)
  nGlobalSurfSides = sendbuf
END IF

IF (mySurfRank.EQ.0) THEN
#if USE_LOADBALANCE
  IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
    WRITE(UNIT_stdOUt,'(A,I0,A)') ' Starting surface communication between ', nSurfLeaders, ' compute nodes... DONE!'
END IF

END SUBROUTINE InitSurfCommunication


SUBROUTINE ExchangeSurfData()
!===================================================================================================================================
! exchange the surface data
!> 1) collect the information on the local compute-node
!> 2) compute-node leaders with sampling sides in their halo region and the original node communicate the sampling information
!> 3) compute-node leaders ensure synchronization of shared arrays on their node
!!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSampSize,nSurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState,SampWallState_Shared,SampWallState_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallPumpCapacity,SampWallPumpCapacity_Shared,SampWallPumpCapacity_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactEnergy,SampWallImpactEnergy_Shared,SampWallImpactEnergy_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactVector,SampWallImpactVector_Shared,SampWallImpactVector_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactAngle ,SampWallImpactAngle_Shared ,SampWallImpactAngle_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactNumber,SampWallImpactNumber_Shared,SampWallImpactNumber_Shared_Win
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iProc,SideID
INTEGER                         :: iPos,p,q
INTEGER                         :: MessageSize,iSurfSide,SurfSideID
INTEGER                         :: nValues
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = SurfSampSize*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(SampWallState,SampWallState_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(SampWallState,0                   ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ENDIF
!
IF(nPorousBC.GT.0) THEN
  MessageSize = nComputeNodeSurfTotalSides
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_REDUCE(SampWallPumpCapacity,SampWallPumpCapacity_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  ELSE
    CALL MPI_REDUCE(SampWallPumpCapacity,0                          ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  END IF
END IF
! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
IF (CalcSurfaceImpact) THEN
  IF (myComputeNodeRank.EQ.0) THEN
    MessageSize = nSpecies*4*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactEnergy,SampWallImpactEnergy_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    MessageSize = nSpecies*3*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactVector,SampWallImpactVector_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    MessageSize = nSpecies*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactAngle ,SampWallImpactAngle_Shared ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    CALL MPI_REDUCE(SampWallImpactNumber,SampWallImpactNumber_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  ELSE
    MessageSize = nSpecies*4*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactEnergy,0                          ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    MessageSize = nSpecies*3*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactVector,0                          ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    MessageSize = nSpecies*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
    CALL MPI_REDUCE(SampWallImpactAngle ,0                          ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
    CALL MPI_REDUCE(SampWallImpactNumber,0                          ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  END IF
END IF

CALL BARRIER_AND_SYNC(SampWallState_Shared_Win         ,MPI_COMM_SHARED)
IF(nPorousBC.GT.0) THEN
  CALL BARRIER_AND_SYNC(SampWallPumpCapacity_Shared_Win,MPI_COMM_SHARED)
END IF
IF (CalcSurfaceImpact) THEN
  CALL BARRIER_AND_SYNC(SampWallImpactEnergy_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactVector_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactAngle_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactNumber_Shared_Win,MPI_COMM_SHARED)
END IF

! prepare buffers for surf leader communication
IF (myComputeNodeRank.EQ.0) THEN
  nValues = SurfSampSize*nSurfSample**2
  ! Sampling of impact energy for each species (trans, rot, vib, elec), impact vector (x,y,z), angle and number: Add 9*nSpecies
  ! to the buffer length
  IF(CalcSurfaceImpact) nValues=nValues+9*nSpecies
  IF(nPorousBC.GT.0) nValues = nValues + 1

  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvSurfSides * nValues
    CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE

    ! Only assemble message if we are expecting sides to send to this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    SurfSendBuf(iProc)%content = 0.

    DO iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      SideID     = SurfMapping(iProc)%SendSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      ! Assemble message
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          SurfSendBuf(iProc)%content(iPos+1:iPos+SurfSampSize) = SampWallState_Shared(:,p,q,SurfSideID)
          iPos = iPos + SurfSampSize
          ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
          IF (CalcSurfaceImpact) THEN
            ! Add average impact energy for each species (trans, rot, vib)
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactEnergy_Shared(:,1,p,q,SurfSideID)
            iPos = iPos + nSpecies
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactEnergy_Shared(:,2,p,q,SurfSideID)
            iPos=iPos + nSpecies
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactEnergy_Shared(:,3,p,q,SurfSideID)
            iPos=iPos + nSpecies
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactEnergy_Shared(:,4,p,q,SurfSideID)
            iPos=iPos + nSpecies

            ! Add average impact vector (x,y,z) for each species
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactVector_Shared(:,1,p,q,SurfSideID)
            iPos = iPos + nSpecies
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactVector_Shared(:,2,p,q,SurfSideID)
            iPos = iPos + nSpecies
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactVector_Shared(:,3,p,q,SurfSideID)
            iPos = iPos + nSpecies

            ! Add average impact angle for each species
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactAngle_Shared(:,p,q,SurfSideID)
            iPos = iPos + nSpecies

            ! Add number of particle impacts
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = SampWallImpactNumber_Shared(:,p,q,SurfSideID)
            iPos = iPos + nSpecies
          END IF ! CalcSurfaceImpact
        END DO ! p=0,nSurfSample
      END DO ! q=0,nSurfSample
      IF(nPorousBC.GT.0) THEN
        SurfSendBuf(iProc)%content(iPos+1:iPos+1) = SampWallPumpCapacity_Shared(SurfSideID)
        iPos = iPos + 1
      END IF

      SampWallState_Shared(:,:,:,SurfSideID)=0.
      ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
      IF (CalcSurfaceImpact) THEN
        SampWallImpactEnergy_Shared(:,:,:,:,SurfSideID) = 0.
        SampWallImpactVector_Shared(:,:,:,:,SurfSideID) = 0.
        SampWallImpactAngle_Shared (:,:,:,SurfSideID)   = 0.
        SampWallImpactNumber_Shared(:,:,:,SurfSideID)   = 0.
      END IF ! CalcSurfaceImpact
      IF(nPorousBC.GT.0) THEN
        SampWallPumpCapacity_Shared(SurfSideID) = 0.
      END IF
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendSurfSides * nValues
    CALL MPI_ISEND( SurfSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
      SideID     = SurfMapping(iProc)%RecvSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      DO q=1,nSurfSample
        DO p=1,nSurfSample
          SampWallState_Shared(:,p,q,SurfSideID) = SampWallState_Shared(:,p,q,SurfSideID) &
                                                 + SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfSampSize)
          iPos = iPos + SurfSampSize
          ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
          IF(CalcSurfaceImpact)THEN
            ! Add average impact energy for each species (trans, rot, vib)
            SampWallImpactEnergy_Shared(:,1,p,q,SurfSideID) = SampWallImpactEnergy_Shared(:,1,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            SampWallImpactEnergy_Shared(:,2,p,q,SurfSideID) = SampWallImpactEnergy_Shared(:,2,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            SampWallImpactEnergy_Shared(:,3,p,q,SurfSideID) = SampWallImpactEnergy_Shared(:,3,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            SampWallImpactEnergy_Shared(:,4,p,q,SurfSideID) = SampWallImpactEnergy_Shared(:,4,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            ! Add average impact vector (x,y,z) for each species
            SampWallImpactVector_Shared(:,1,p,q,SurfSideID) = SampWallImpactVector_Shared(:,1,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            SampWallImpactVector_Shared(:,2,p,q,SurfSideID) = SampWallImpactVector_Shared(:,2,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            SampWallImpactVector_Shared(:,3,p,q,SurfSideID) = SampWallImpactVector_Shared(:,3,p,q,SurfSideID) &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            ! Add average impact angle for each species
            SampWallImpactAngle_Shared(:,p,q,SurfSideID)    = SampWallImpactAngle_Shared(:,p,q,SurfSideID)    &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
            ! Add number of particle impacts
            SampWallImpactNumber_Shared(:,p,q,SurfSideID)   = SampWallImpactNumber_Shared(:,p,q,SurfSideID)   &
                                                            + SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos = iPos + nSpecies
          END IF ! CalcSurfaceImpact
        END DO ! p = 0,nSurfSample
      END DO ! q = 0,nSurfSample
      IF(nPorousBC.GT.0) THEN
        SampWallPumpCapacity_Shared(SurfSideID) = SurfRecvBuf(iProc)%content(iPos+1)
        iPos = iPos + 1
      END IF
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
     ! Nullify buffer
    SurfRecvBuf(iProc)%content = 0.
  END DO ! iProc
END IF

! ensure synchronization on compute node
CALL BARRIER_AND_SYNC(SampWallState_Shared_Win         ,MPI_COMM_SHARED)
IF(nPorousBC.GT.0) THEN
  CALL BARRIER_AND_SYNC(SampWallPumpCapacity_Shared_Win,MPI_COMM_SHARED)
END IF
IF (CalcSurfaceImpact) THEN
  CALL BARRIER_AND_SYNC(SampWallImpactEnergy_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactVector_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactAngle_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactNumber_Shared_Win,MPI_COMM_SHARED)
END IF

END SUBROUTINE ExchangeSurfData


SUBROUTINE FinalizeSurfCommunication()
!----------------------------------------------------------------------------------------------------------------------------------!
! Deallocated arrays used for sampling surface communication
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,mySurfRank
USE MOD_MPI_Shared_Vars         ,ONLY: MPIRankSharedLeader,MPIRankSurfLeader
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc
!===================================================================================================================================

IF (myComputeNodeRank.NE.0) RETURN

! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

SDEALLOCATE(MPIRankSharedLeader)
SDEALLOCATE(MPIRankSurfLeader)

DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
    SDEALLOCATE(SurfMapping(iProc)%RecvSurfGlobalID)
    SDEALLOCATE(SurfMapping(iProc)%RecvPorousGlobalID)
    SDEALLOCATE(SurfRecvBuf(iProc)%content)
  END IF

  IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
    SDEALLOCATE(SurfMapping(iProc)%SendSurfGlobalID)
    SDEALLOCATE(SurfMapping(iProc)%SendPorousGlobalID)
    SDEALLOCATE(SurfSendBuf(iProc)%content)
  END IF
END DO
SDEALLOCATE(SurfMapping)
SDEALLOCATE(SurfSendBuf)
SDEALLOCATE(SurfRecvBuf)

END SUBROUTINE FinalizeSurfCommunication
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Boundary_Sampling
