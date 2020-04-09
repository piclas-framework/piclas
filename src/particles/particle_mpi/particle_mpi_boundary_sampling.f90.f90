!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
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
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE InitSurfCommunication
  MODULE PROCEDURE InitSurfCommunication
END INTERFACE


PUBLIC :: InitSurfCommunication
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
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank,nSurfLeaders,nSurfCommProc
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides,offsetComputeNodeSurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide,SurfSide2GlobalSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                       :: iProc,color
INTEGER                       :: iLeader,LeaderID
INTEGER                       :: leadersGroup,surfGroup
INTEGER                       :: iSide
INTEGER                       :: sendbuf,recvbuf
INTEGER                       :: nSendSurfSidesTmp(0:nLeaderGroupProcs-1)
INTEGER                       :: nRecvSurfSidesTmp(0:nLeaderGroupProcs-1)
!INTEGER                       :: nSurfSidesLeader(1:2,0:nLeaderGroupProcs-1)
INTEGER                       :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                       :: SendSurfGlobalID(0:nLeaderGroupProcs-1,1:nComputeNodeSurfTotalSides)
!===================================================================================================================================
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
nRecvSurfSidesTmp = 0
nSendSurfSidesTmp = 0

DO iSide = 1,nComputeNodeSurfTotalSides
  ! count surf sides per compute node
  LeaderID = SurfSide2GlobalSide(SURF_LEADER,iSide)
  nSendSurfSidesTmp(LeaderID) = nSendSurfSidesTmp(LeaderID) + 1
  SendSurfGlobalID(LeaderID,nSendSurfSidesTmp(LeaderID)) = SurfSide2GlobalSide(SURF_SIDEID,iSide)
END DO

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

  CALL MPI_WAIT(SendRequest(iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

!--- Split communicator from MPI_COMM_LEADER_SHARED
color = MPI_UNDEFINED
IF (SurfOnNode) color = 1201

! create new SurfMesh communicator for SurfMesh communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_LEADERS_SHARED, color, MPI_INFO_NULL, MPI_COMM_LEADERS_SURF, IERROR)

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
SWRITE(UNIT_stdOUt,'(A,I3,A)') ' Starting surface communication between ', nSurfLeaders, ' compute nodes'

!!--- Count all communicated sides and build mapping for other leaders
!ALLOCATE(nSurfSidesLeader(1:2,0:nSurfLeaders-1))
!
!nSurfCommProc = 0
!DO iProc = 0,nLeaderGroupProcs-1
!  ! a leader defines itself as if it has surf sides within its local domain. However, there might be procs which neither send nor
!  ! receive sides from us. We can reduce nSurfLeaders to nSurfCommProc
!  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
!!  IF ((nRecvSurfSidesTmp(iProc).EQ.0) .AND. (nSendSurfSidesTmp(iProc).EQ.0)) CYCLE
!
!  ! MPI ranks, start at 0
!  nSurfSidesLeader(1,nSurfCommProc) = nSendSurfSidesTmp(iProc)
!  nSurfSidesLeader(2,nSurfCommProc) = nRecvSurfSidesTmp(iProc)
!  nSurfCommProc = nSurfCommProc + 1
!END DO

!--- Open receive buffer (mapping from message surface ID to global side ID)
ALLOCATE(SurfMapping(0:nSurfLeaders-1))
DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Save number of send and recv sides
  SurfMapping(iProc)%nRecvSurfSides = nSendSurfSidesTmp(MPIRankSurfLeader(iProc))
  SurfMapping(iProc)%nSendSurfSides = nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))

  ! Only open recv buffer if we are expecting sides from this leader node
  IF (nRecvSurfSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%RecvSurfGlobalID(1:nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))))

  CALL MPI_IRECV( SurfMapping(iProc)%RecvSurfGlobalID                         &
                , nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Only open send buffer if we are expecting sides from this leader node
  IF (nSendSurfSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%SendSurfGlobalID(1:nSendSurfSidesTmp(MPIRankSurfLeader(iProc))))

  SurfMapping(iProc)%SendSurfGlobalID = SendSurfGlobalID(MPIRankSurfLeader(iProc),1:nSendSurfSidesTmp(MPIRankSurfLeader(iProc)))

  CALL MPI_ISEND( SurfMapping(iProc)%SendSurfGlobalID                         &
                , nSendSurfSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  IF (nSendSurfSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(SendRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF

  IF (nRecvSurfSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(RecvRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

!--- Save number of total surf sides
IF (surfOnNode) THEN
  IF (nSurfLeaders.EQ.1) THEN
    offsetComputeNodeSurfSide = 0
    nSurfTotalSides           = nComputeNodeSurfSides
  ELSE
    sendbuf = nComputeNodeSurfSides
    recvbuf = 0
    CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
    offsetComputeNodeSurfSide = recvbuf
    ! last proc knows CN total number of BC elems
    sendbuf = offsetComputeNodeSurfSide + nSurfTotalSides
    CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nSurfLeaders-1,MPI_COMM_LEADERS_SURF,iError)
    nSurfTotalSides = sendbuf
  END IF
END IF


END SUBROUTINE InitSurfCommunication


!SUBROUTINE InitSurfCommunicator()
!!===================================================================================================================================
!! Creates two new subcommunicators.
!! SurfCOMM%COMM contains all MPI-Ranks which have reflective boundary faces in their halo-region and process which have reflective
!! boundary faces in their origin region. This communicator is used to communicate the wall-sampled values of halo-faces to the
!! origin face
!! SurfCOMM%OutputCOMM is another subset. This communicator contains only the processes with origin surfaces. It is used to perform
!! collective writes of the surf-sampled values.
!! Sets also used for communication of adsorption variables
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Particle_Boundary_Vars   ,ONLY:SurfMesh
!USE MOD_Particle_Boundary_Vars   ,ONLY:SurfCOMM
!USE MOD_Particle_MPI_Vars        ,ONLY:PartMPI
!USE MOD_Particle_Boundary_Vars   ,ONLY:OffSetSurfSideMPI,OffSetSurfSide,OffSetInnerSurfSideMPI,OffSetInnerSurfSide
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                   :: color,iProc
!INTEGER                   :: noSurfrank,Surfrank
!LOGICAL                   :: hasSurf
!INTEGER,ALLOCATABLE       :: countSurfSideMPI(:),countInnerSurfSideMPI(:)
!LOGICAL                   :: OutputOnProc, InnerSlaveBCs
!!===================================================================================================================================
!color=MPI_UNDEFINED
!IF(SurfMesh%SurfOnProc) color=1001
!! THEN
!! color=2
!! ELSE
!! color=1
!! END IF
!! create ranks for RP communicator
!IF(PartMPI%MPIRoot) THEN
!  Surfrank=-1
!  noSurfrank=-1
!  SurfCOMM%Myrank=0
!  IF(SurfMesh%SurfOnProc) THEN
!    Surfrank=0
!  ELSE
!    noSurfrank=0
!  END IF
!  DO iProc=1,nProcessors-1
!    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
!    IF(hasSurf) THEN
!      SurfRank=SurfRank+1
!      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
!    ELSE
!      noSurfRank=noSurfRank+1
!      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
!    END IF
!  END DO
!ELSE
!  CALL MPI_SEND(SurfMesh%SurfOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
!  CALL MPI_RECV(SurfCOMM%MyRank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
!END IF
!
!! create new SurfMesh communicator for SurfMesh communication
!CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyRank, SurfCOMM%COMM,iError)
!IF(SurfMesh%SurfOnPRoc) THEN
!  CALL MPI_COMM_SIZE(SurfCOMM%COMM, SurfCOMM%nProcs,iError)
!ELSE
!  SurfCOMM%nProcs = 0
!END IF
!SurfCOMM%MPIRoot=.FALSE.
!IF(SurfCOMM%MyRank.EQ.0 .AND. SurfMesh%SurfOnProc) THEN
!  SurfCOMM%MPIRoot=.TRUE.
!!   WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF COMM:        ',SurfCOMM%nProcs,' procs'
!END IF
!
!! now, create output communicator
!OutputOnProc=.FALSE.
!color=MPI_UNDEFINED
!IF(SurfMesh%nOutputSides.GT.0) THEN
!  OutputOnProc=.TRUE.
!  color=1002
!END IF
!
!IF(PartMPI%MPIRoot) THEN
!  Surfrank=-1
!  noSurfrank=-1
!  SurfCOMM%MyOutputRank=0
!  IF(SurfMesh%nOutputSides.GT.0) THEN
!    Surfrank=0
!  ELSE
!    noSurfrank=0
!  END IF
!  DO iProc=1,nProcessors-1
!    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
!    IF(hasSurf) THEN
!      SurfRank=SurfRank+1
!      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
!    ELSE
!      noSurfRank=noSurfRank+1
!      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
!    END IF
!  END DO
!ELSE
!  CALL MPI_SEND(OutputOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
!  CALL MPI_RECV(SurfCOMM%MyOutputRank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
!END IF
!
!! create new SurfMesh Output-communicator
!CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyOutputRank, SurfCOMM%OutputCOMM,iError)
!IF(OutputOnPRoc)THEN
!  CALL MPI_COMM_SIZE(SurfCOMM%OutputCOMM, SurfCOMM%nOutputProcs,iError)
!ELSE
!  SurfCOMM%nOutputProcs = 0
!END IF
!SurfCOMM%MPIOutputRoot=.FALSE.
!IF(SurfCOMM%MyOutputRank.EQ.0 .AND. OutputOnProc) THEN
!  SurfCOMM%MPIOutputRoot=.TRUE.
!!   WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF OUTPUT-COMM: ',SurfCOMM%nOutputProcs,' procs'
!END IF
!
!IF(SurfMesh%nTotalSides.EQ.0) RETURN
!! check if any proc has innersides and set flag for all proc to do additional communication and slave mapping
!IF((SurfMesh%nSides-SurfMesh%nOutputSides).GT.0) THEN
!  InnerSlaveBCs = .TRUE.
!ELSE
!  InnerSlaveBCs = .FALSE.
!END IF
!CALL MPI_ALLREDUCE(InnerSlaveBCs,SurfCOMM%InnerBCs,1,MPI_LOGICAL,MPI_LOR,SurfCOMM%COMM,iError)
!
!IF(SurfMesh%nOutputSides.EQ.0) RETURN
!! get correct offsets for output of hdf5 file (master sides)
!ALLOCATE(offsetSurfSideMPI(0:SurfCOMM%nOutputProcs))
!offsetSurfSideMPI=0
!ALLOCATE(countSurfSideMPI(0:SurfCOMM%nOutputProcs-1))
!countSurfSideMPI=0
!
!CALL MPI_GATHER(SurfMesh%nOutputSides,1,MPI_INTEGER,countSurfSideMPI,1,MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
!
!! new offsets due to InnerSurfSides
!ALLOCATE(offsetInnerSurfSideMPI(0:SurfCOMM%nOutputProcs))
!offsetInnerSurfSideMPI=0
!ALLOCATE(countInnerSurfSideMPI(0:SurfCOMM%nOutputProcs-1))
!countInnerSurfSideMPI=0
!
!CALL MPI_GATHER(SurfMesh%nInnerSides,1,MPI_INTEGER,countInnerSurfSideMPI,1,MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
!
!IF (SurfCOMM%MPIOutputRoot) THEN
!  DO iProc=1,SurfCOMM%nOutputProcs-1
!    offsetSurfSideMPI(iProc)=SUM(countSurfSideMPI(0:iProc-1))-SUM(countInnerSurfSideMPI(0:iProc-1))
!    offsetInnerSurfSideMPI(iProc)=SUM(countInnerSurfSideMPI(0:iProc-1))
!  END DO
!  offsetSurfSideMPI(SurfCOMM%nOutputProcs)=SUM(countSurfSideMPI(:))-SUM(countInnerSurfSideMPI(:))
!  offsetInnerSurfSideMPI(SurfCOMM%nOutputProcs)=SUM(countInnerSurfSideMPI(:))
!  ! add BC offset to InnerSurfSide offset
!  offsetInnerSurfSideMPI(0:SurfCOMM%nOutputProcs) = &
!  offsetInnerSurfSideMPI(0:SurfCOMM%nOutputProcs) + offsetSurfSideMPI(SurfCOMM%nOutputProcs)
!END IF
!
!CALL MPI_BCAST (offsetSurfSideMPI,size(offsetSurfSideMPI),MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
!offsetSurfSide=offsetSurfSideMPI(SurfCOMM%MyOutputRank)
!CALL MPI_BCAST (offsetInnerSurfSideMPI,size(offsetInnerSurfSideMPI),MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
!offsetInnerSurfSide=offsetInnerSurfSideMPI(SurfCOMM%MyOutputRank)
!
!END SUBROUTINE InitSurfCommunicator


!SUBROUTINE ExchangeSurfData()
!!===================================================================================================================================
!! exchange the surface data
!! only processes with sampling sides in their halo region and the original process participate on the communication
!! structure is similar to particle communication
!! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
!! the receiving process adds the new data to his own sides
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Particle_Vars               ,ONLY:nSpecies
!USE MOD_SurfaceModel_Vars           ,ONLY:Adsorption
!USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample,SampWall,PartBound,CalcSurfaceImpact
!USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                         :: MessageSize,iSurfSide,SurfSideID
!INTEGER                         :: nValues, nReactiveValues
!INTEGER                         :: iPos,p,q,iProc,iReact
!INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!!===================================================================================================================================
!
!nValues = SurfMesh%SampSize*nSurfSample**2
!! additional array entries for Coverage, Accomodation and recombination coefficient
!nReactiveValues=0
!IF(ANY(PartBound%Reactive)) nReactiveValues = SurfMesh%ReactiveSampSize*(nSurfSample)**2
!
!nValues=nValues+nReactiveValues
!
!! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number: Add 8*nSpecies to the
!! buffer length
!IF(CalcSurfaceImpact) nValues=nValues+8*nSpecies
!
!! open receive buffer
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
!  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
!                , MessageSize                                  &
!                , MPI_DOUBLE_PRECISION                         &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
!                , 1009                                         &
!                , SurfCOMM%COMM                                &
!                , SurfExchange%RecvRequest(iProc)              &
!                , IERROR )
!END DO ! iProc
!
!! build message
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  iPos=0
!  SurfSendBuf(iProc)%content = 0.
!  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
!    SurfSideID=SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide))
!    DO q=1,nSurfSample
!      DO p=1,nSurfSample
!        SurfSendBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)= SampWall(SurfSideID)%State(:,p,q)
!        iPos=iPos+SurfMesh%SampSize
!        IF (ANY(PartBound%Reactive)) THEN
!          SurfSendBuf(iProc)%content(iPos+1:iPos+5+nSpecies)= SampWall(SurfSideID)%SurfModelState(:,p,q)
!          iPos=iPos+5+nSpecies
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%Accomodation(:,p,q)
!          iPos=iPos+nSpecies
!          DO iReact=1,2*Adsorption%ReactNum
!            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%SurfModelReactCount(iReact,:,p,q)
!            iPos=iPos+nSpecies
!          END DO
!        END IF
!
!        ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!        IF(CalcSurfaceImpact)THEN
!          ! Add average impact energy for each species (trans, rot, vib)
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactEnergy(:,1,p,q)
!          iPos=iPos+nSpecies
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactEnergy(:,2,p,q)
!          iPos=iPos+nSpecies
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactEnergy(:,3,p,q)
!          iPos=iPos+nSpecies
!
!          ! Add average impact vector (x,y,z) for each species
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactVector(:,1,p,q)
!          iPos=iPos+nSpecies
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactVector(:,2,p,q)
!          iPos=iPos+nSpecies
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactVector(:,3,p,q)
!          iPos=iPos+nSpecies
!
!          ! Add average impact angle for each species
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactAngle(:,p,q)
!          iPos=iPos+nSpecies
!
!          ! Add number of particle impacts
!          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%ImpactNumber(:,p,q)
!          iPos=iPos+nSpecies
!        END IF ! CalcSurfaceImpact
!
!      END DO ! p=0,nSurfSample
!    END DO ! q=0,nSurfSample
!    SampWall(SurfSideID)%State(:,:,:)=0.
!    IF (ANY(PartBound%Reactive)) THEN
!      SampWall(SurfSideID)%SurfModelState(:,:,:)=0.
!      SampWall(SurfSideID)%Accomodation(:,:,:)=0.
!      SampWall(SurfSideID)%SurfModelReactCount(:,:,:,:)=0.
!    END IF
!    ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!    IF(CalcSurfaceImpact)THEN
!      SampWall(SurfSideID)%ImpactEnergy(:,:,:,:)=0.
!      SampWall(SurfSideID)%ImpactVector(:,:,:,:)=0.
!      SampWall(SurfSideID)%ImpactAngle(:,:,:)=0.
!      SampWall(SurfSideID)%ImpactNumber(:,:,:)=0.
!    END IF ! CalcSurfaceImpact
!  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
!END DO
!
!! send message
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
!  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
!                , MessageSize                              &
!                , MPI_DOUBLE_PRECISION                     &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
!                , 1009                                     &
!                , SurfCOMM%COMM                            &
!                , SurfExchange%SendRequest(iProc)          &
!                , IERROR )
!END DO ! iProc
!
!! 4) Finish Received number of particles
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!END DO ! iProc
!
!! add data do my list
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  iPos=0
!  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
!    SurfSideID=SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide))
!    DO q=1,nSurfSample
!      DO p=1,nSurfSample
!        SampWall(SurfSideID)%State(:,p,q)=SampWall(SurfSideID)%State(:,p,q) &
!                                         +SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)
!        iPos=iPos+SurfMesh%SampSize
!        IF (ANY(PartBound%Reactive)) THEN
!          SampWall(SurfSideID)%SurfModelState(:,p,q)=SampWall(SurfSideID)%SurfModelState(:,p,q) &
!                                                +SurfRecvBuf(iProc)%content(iPos+1:iPos+5+nSpecies)
!          iPos=iPos+5+nSpecies
!          SampWall(SurfSideID)%Accomodation(:,p,q)=SampWall(SurfSideID)%Accomodation(:,p,q) &
!                                                  +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!          DO iReact=1,2*Adsorption%ReactNum
!            SampWall(SurfSideID)%SurfModelReactCount(iReact,:,p,q)=SampWall(SurfSideID)%SurfModelReactCount(iReact,:,p,q) &
!                                                       +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!            iPos=iPos+nSpecies
!          END DO
!        END IF
!
!        ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
!        IF(CalcSurfaceImpact)THEN
!          ! Add average impact energy for each species (trans, rot, vib)
!          SampWall(SurfSideID)%ImpactEnergy(:,1,p,q)=SampWall(SurfSideID)%ImpactEnergy(:,1,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!          SampWall(SurfSideID)%ImpactEnergy(:,2,p,q)=SampWall(SurfSideID)%ImpactEnergy(:,2,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!          SampWall(SurfSideID)%ImpactEnergy(:,3,p,q)=SampWall(SurfSideID)%ImpactEnergy(:,3,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!
!          ! Add average impact vector (x,y,z) for each species
!          SampWall(SurfSideID)%ImpactVector(:,1,p,q)=SampWall(SurfSideID)%ImpactVector(:,1,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!          SampWall(SurfSideID)%ImpactVector(:,2,p,q)=SampWall(SurfSideID)%ImpactVector(:,2,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!          SampWall(SurfSideID)%ImpactVector(:,3,p,q)=SampWall(SurfSideID)%ImpactVector(:,3,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!
!          ! Add average impact angle for each species
!          SampWall(SurfSideID)%ImpactAngle(:,p,q)=SampWall(SurfSideID)%ImpactAngle(:,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!
!          ! Add number of particle impacts
!          SampWall(SurfSideID)%ImpactNumber(:,p,q)=SampWall(SurfSideID)%ImpactNumber(:,p,q) &
!                                                    +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
!          iPos=iPos+nSpecies
!        END IF ! CalcSurfaceImpact
!
!      END DO ! p=0,nSurfSample
!    END DO ! q=0,nSurfSample
!  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
!  SurfRecvBuf(iProc)%content = 0.
!END DO ! iProc
!
!END SUBROUTINE ExchangeSurfData
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Boundary_Sampling
