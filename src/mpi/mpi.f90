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

MODULE MOD_MPI
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi_f08
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC::InitMPI

#if USE_MPI
INTERFACE StartReceiveMPIData
  MODULE PROCEDURE StartReceiveMPIData0D
  MODULE PROCEDURE StartReceiveMPIData2D
END INTERFACE

INTERFACE StartSendMPIData
  MODULE PROCEDURE StartSendMPIData0D
  MODULE PROCEDURE StartSendMPIData2D
END INTERFACE

PUBLIC :: InitMPIvars,StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData,FinalizeMPI
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
PUBLIC :: StartReceiveMPIDataType,StartSendMPIDataType,FinishExchangeMPIDataType
#if !(USE_HDG)
PUBLIC :: StartSendMPIDataTypeDielectric,FinishExchangeMPIDataTypeDielectric,StartReceiveMPIDataTypeDielectric
#endif /*!(USE_HDG)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
PUBLIC :: StartExchange_DG_Elems
PUBLIC :: StartReceiveMPIDataInt,StartSendMPIDataInt
#endif /*USE_MPI*/
PUBLIC::DefineParametersMPI
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMPI()
! MODULES
USE MOD_ReadInTools       ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("MPI")
CALL prms%CreateIntOption('GroupSize', "Define size of MPI subgroups, used to e.g. perform grouped IO, where group master\n"//&
                                       "collects and outputs data.",&
                                       '0')
#if defined(PARTICLES)
CALL prms%CreateLogicalOption('CheckExchangeProcs' , 'Check if proc communication of particle info is non-symmetric', '.TRUE.')
CALL prms%CreateLogicalOption('AbortExchangeProcs' , 'Abort if proc communication of particle info is non-symmetric (requires CheckExchangeProcs=T)', '.TRUE.')
CALL prms%CreateLogicalOption('DoParticleLatencyHiding' , 'Determine the elements that require particle communication and use them for latency hiding', '.FALSE.')
#endif /*PARTICLES*/
END SUBROUTINE DefineParametersMPI


!===================================================================================================================================
!> Basic MPI initialization.
!===================================================================================================================================
#if USE_MPI
SUBROUTINE InitMPI(mpi_comm_IN)
#else
SUBROUTINE InitMPI()
#endif /*USE_MPI*/
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_MPI
TYPE(mpi_comm),INTENT(IN),OPTIONAL      :: mpi_comm_IN !< MPI communicator
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
TYPE(mpi_comm) :: MPI_COMM_LOC
LOGICAL :: initDone
!==================================================================================================================================
IF (PRESENT(mpi_comm_IN)) THEN
  MPI_COMM_LOC = mpi_comm_IN
ELSE
  CALL MPI_INIT(iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_INIT',iError)
  CALL MPI_INITIALIZED(initDone,iError)
  IF(.NOT.initDone) CALL MPI_INIT(iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_INITIALIZED',iError)
  ! General communicator
  CALL MPI_COMM_DUP(MPI_COMM_WORLD,MPI_COMM_PICLAS,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_DUP',iError)
  MPI_COMM_LOC = MPI_COMM_PICLAS
END IF

CALL MPI_COMM_RANK(MPI_COMM_LOC, myRank     , iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
CALL MPI_COMM_SIZE(MPI_COMM_LOC, nProcessors, iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Could not get rank and number of processors',iError)
MPIRoot=(myRank .EQ. 0)
#else  /*USE_MPI*/
myRank      = 0
myLocalRank = 0
nProcessors = 1
MPIRoot     =.TRUE.
MPILocalRoot=.TRUE.
#endif  /*USE_MPI*/

END SUBROUTINE InitMPI


#if USE_MPI
!===================================================================================================================================
!> Initialize derived MPI types used for communication and allocate HALO data.
!===================================================================================================================================
SUBROUTINE InitMPIvars()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Readintools        ,ONLY: GETINT
USE MOD_MPI_Shared_Vars    ,ONLY: nProcessors_Global
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_FV
USE MOD_MPI_FV             ,ONLY: InitMPIvarsFV
#endif /*USE_FV*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: color,groupsize
!===================================================================================================================================
IF(.NOT.InterpolationInitIsDone) CALL Abort(__STAMP__,'InitMPITypes called before InitInterpolation')
#if USE_FV
CALL InitMPIvarsFV()
#endif /*USE_FV*/
ALLOCATE(SendRequest_U(nNbProcs)     )
ALLOCATE(SendRequest_U2(nNbProcs)     )
ALLOCATE(SendRequest_GEO(nNbProcs)     )
ALLOCATE(SendRequest_Flux(nNbProcs)  )
ALLOCATE(RecRequest_U(nNbProcs)     )
ALLOCATE(RecRequest_U2(nNbProcs)     )
ALLOCATE(RecRequest_Geo(nNbProcs)     )
ALLOCATE(RecRequest_Flux(nNbProcs)  )
SendRequest_U      = MPI_REQUEST_NULL
SendRequest_U2      = MPI_REQUEST_NULL
SendRequest_Flux   = MPI_REQUEST_NULL
RecRequest_U       = MPI_REQUEST_NULL
RecRequest_U2       = MPI_REQUEST_NULL
RecRequest_Flux    = MPI_REQUEST_NULL
SendRequest_Geo    = MPI_REQUEST_NULL
RecRequest_Geo     = MPI_REQUEST_NULL
DataSizeSide  =(PP_N+1)*(PP_N+1)

! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize','0')
IF(GroupSize.LT.1)THEN ! group procs by node
  ! Split the node communicator (shared memory) from the global communicator on physical processor or node level
#if (CORE_SPLIT==0)
  ! 1.) MPI_COMM_TYPE_SHARED
  ! Note that using SharedMemoryMethod=OMPI_COMM_TYPE_CORE somehow does not work in every case (intel/amd processors)
  ! Also note that OMPI_COMM_TYPE_CORE is undefined when not using OpenMPI
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_PICLAS,SharedMemoryMethod,0,MPI_INFO_NULL,MPI_COMM_NODE,IERROR)
#elif (CORE_SPLIT==1)
  ! 2.) OMPI_COMM_TYPE_CORE
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,myRank,0,MPI_COMM_NODE,iError)
#else
  ! 3.) PICLAS_COMM_TYPE_NODE
  ! Check if more nodes than procs are required or
  ! if the resulting split would create unequal procs per node
  IF((CORE_SPLIT.GE.nProcessors_Global).OR.(MOD(nProcessors_Global,CORE_SPLIT).GT.0))THEN
    LBWRITE (*,'(A,I0,A,I0,A,F0.2,A)') ' WARNING: Either more nodes than cores selected (nodes: ',CORE_SPLIT,', cores: ',&
        nProcessors_Global,') OR unequal number of cores per node (=',REAL(nProcessors_Global)/REAL(CORE_SPLIT),&
        '). Setting 1 core per node for MPI_COMM_NODE!'
    color = myRank
  ELSE
    ! Group procs so that every CORE_SPLIT procs are in the same group
    color = INT(REAL(myrank*CORE_SPLIT)/REAL(nProcessors_Global))+1
  END IF ! (CORE_SPLIT.GE.nProcessors_Global).OR.(MOD().GT.0)
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,color,0,MPI_COMM_NODE,iError)
#endif
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,color,0,MPI_COMM_NODE,iError)
END IF
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
MPILocalRoot=(myLocalRank.EQ.0)

IF (nProcessors_Global.EQ.nLocalProcs) THEN
  LBWRITE(UNIT_stdOUt,'(A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs in ',1,' group'
ELSE
  LBWRITE(UNIT_stdOUt,'(A,I0,A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs each in ',&
                                                        nProcessors_Global/nLocalProcs,' groups for a total number of ',&
                                                        nProcessors_Global,' procs'
END IF

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,0,0,MPI_COMM_LEADERS,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,1,0,MPI_COMM_WORKERS,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
  CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
  CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
  nLeaderProcs=nProcessors-nWorkerProcs
END IF
END SUBROUTINE InitMPIvars


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors.
!===================================================================================================================================
SUBROUTINE StartReceiveMPIData0D(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID                                                 !< defines the send / receive direction -> 1=send MINE
                                                                              !< / receive YOUR, 3=send YOUR / receive MINE
INTEGER,INTENT(IN)  :: firstDim                                               !< size of one entry in array (e.g. one side:
                                                                              !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)  :: LowerBound                                             !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)  :: UpperBound                                             !< upper side index for last dimension of FaceData
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
REAL,INTENT(OUT)    :: FaceData(firstDim,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =firstDim*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIData0D


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors.
!===================================================================================================================================
SUBROUTINE StartReceiveMPIData2D(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID                                                 !< defines the send / receive direction -> 1=send MINE
                                                                              !< / receive YOUR, 3=send YOUR / receive MINE
INTEGER,INTENT(IN)  :: firstDim                                               !< size of one entry in array (e.g. one side:
                                                                              !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)  :: LowerBound                                             !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)  :: UpperBound                                             !< upper side index for last dimension of FaceData
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
REAL,INTENT(OUT)              :: FaceData(firstDim,0:PP_N,0:PP_N,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =firstDim*DataSizeSide*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIData2D


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartReceiveMPIDataType(MPIRequest, SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
#if !(USE_HDG) && !(PP_TimeDiscMethod==700)
USE MOD_PML_Vars ,ONLY: PMLnVar,DoPML
#endif /*!(USE_HDG) && !(PP_TimeDiscMethod==700)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID                                                 !< defines the send / receive direction -> 1=send MINE
                                                                              !< / receive YOUR, 3=send YOUR / receive MINE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    ! Send slave U or Flux when no PML is active
    IF(SendID.EQ.2.OR.(.NOT.DoPML))THEN
      nRecVal = PP_nVar*DataSizeSideRec(iNbProc,SendID)
      ! FaceDataRecvU(1:PP_nVar,1:DataSizeSideRec(iNbProc,SendID))
      CALL MPI_IRECV(DGExchange(iNbProc)%FaceDataRecvU,nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
      nRecVal = (PP_nVar+PMLnVar)*DataSizeSideRec(iNbProc,SendID)
      ! FaceDataRecvFlux(1:PP_nVar+PMLnVar,1:DataSizeSideRec(iNbProc,SendID))
      CALL MPI_IRECV(DGExchange(iNbProc)%FaceDataRecvFlux,nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    END IF ! SendID.EQ.2
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIDataType


#if !(USE_HDG)
!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartReceiveMPIDataTypeDielectric(MPIRequest, SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID !< defines the send / receive direction -> 1=send MINE
                              !< / receive YOUR, 3=send YOUR / receive MINE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs) !< communication handles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    IF(SendID.EQ.2)THEN
      ! Send slave U
      nRecVal = PP_nVar*DataSizeSideRec(iNbProc,SendID)
      CALL MPI_IRECV(DGExchange(iNbProc)%FaceDataRecvU,nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
    ! Send mater U
      nRecVal = PP_nVar*DataSizeSideRecMaster(iNbProc,SendID)
      CALL MPI_IRECV(DGExchange(iNbProc)%FaceDataRecvUMaster,nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    END IF ! SendID.EQ.2
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIDataTypeDielectric
#endif /*!(USE_HDG)*/


!===================================================================================================================================
!> See above, but for for send direction
!===================================================================================================================================
SUBROUTINE StartSendMPIData0D(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
INTEGER, INTENT(IN)          :: firstDim,LowerBound,UpperBound
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(OUT)         :: MPIRequest(nNbProcs)
REAL, INTENT(IN)             :: FaceData(firstDim,LowerBound:UpperBound)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =firstDim*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIData0D


!===================================================================================================================================
!> See above, but for for send direction
!===================================================================================================================================
SUBROUTINE StartSendMPIData2D(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
INTEGER, INTENT(IN)          :: firstDim,LowerBound,UpperBound
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)
REAL, INTENT(IN)              :: FaceData(firstDim,0:PP_N,0:PP_N,LowerBound:UpperBound)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =firstDim*DataSizeSide*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIData2D


#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
!===================================================================================================================================
!> See above, but for for send direction (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartSendMPIDataType(MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_DG_Vars  ,ONLY: U_Surf_N,DG_Elems_slave
#if !(USE_HDG)
USE MOD_PML_Vars ,ONLY: PMLnVar,DoPML
#endif /*!(USE_HDG)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(OUT)         :: MPIRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,p,q,iSide,N_slave
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_send(iNbProc,SendID)

    i = 1
    IF(SendID.EQ.2)THEN
      DO iSide = SideID_start, SideID_end
        N_slave = DG_Elems_slave(iSide)
        DO p = 0, N_slave
          DO q = 0, N_slave
            DGExchange(iNbProc)%FaceDataSendU(1:PP_nVar,i) = U_Surf_N(iSide)%U_Slave(1:PP_nVar,p,q)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
      nSendVal = PP_nVar*DataSizeSideSend(iNbProc,SendID)
      ! FaceDataSendU(1:PP_nVar,1:DataSizeSideSend(iNbProc,SendID))
      CALL MPI_ISEND(DGExchange(iNbProc)%FaceDataSendU,nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
      IF(DoPML)THEN
        DO iSide = SideID_start, SideID_end
          N_slave = DG_Elems_slave(iSide)
          DO p = 0, N_slave
            DO q = 0, N_slave
              DGExchange(iNbProc)%FaceDataSendFlux(1:PP_nVar+PMLnVar,i) = U_Surf_N(iSide)%Flux_Slave(1:PP_nVar+PMLnVar,p,q)
              i = i + 1
            END DO ! q = 0, N_slave
          END DO ! p = 0, N_slave
        END DO ! iSide = SideID_start, SideID_end
        nSendVal = (PP_nVar+PMLnVar)*DataSizeSideSend(iNbProc,SendID)
        ! FaceDataSendFlux(1:PP_nVar+PMLnVar,1:DataSizeSideSend(iNbProc,SendID))
        CALL MPI_ISEND(DGExchange(iNbProc)%FaceDataSendFlux,nSendVal,MPI_DOUBLE_PRECISION,  &
                        nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
      ELSE
        ! not PML
        DO iSide = SideID_start, SideID_end
          N_slave = DG_Elems_slave(iSide)
          DO p = 0, N_slave
            DO q = 0, N_slave
              DGExchange(iNbProc)%FaceDataSendU(1:PP_nVar,i) = U_Surf_N(iSide)%Flux_Slave(1:PP_nVar,p,q)
              i = i + 1
            END DO ! q = 0, N_slave
          END DO ! p = 0, N_slave
        END DO ! iSide = SideID_start, SideID_end
        nSendVal = PP_nVar*DataSizeSideSend(iNbProc,SendID)
        ! FaceDataSendFlux(1:PP_nVar+PMLnVar,1:DataSizeSideSend(iNbProc,SendID))
        CALL MPI_ISEND(DGExchange(iNbProc)%FaceDataSendU,nSendVal,MPI_DOUBLE_PRECISION,  &
                        nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
      END IF ! DoPML
    END IF ! SendID.EQ.2
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)

  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIDataType


#if !(USE_HDG)
!===================================================================================================================================
!> See above, but for for send direction (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartSendMPIDataTypeDielectric(MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_DG_Vars         ,ONLY: DG_Elems_slave,DG_Elems_master
USE MOD_Dielectric_Vars ,ONLY: DielectricSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(OUT)         :: MPIRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,p,q,iSide,N_slave,N_master
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_send(iNbProc,SendID)

    ! SendID: Send either master or slave values
    i = 1
    IF(SendID.EQ.2)THEN
      nSendVal = PP_nVar*DataSizeSideSend(iNbProc,SendID)
      DO iSide = SideID_start, SideID_end
        N_slave = DG_Elems_slave(iSide)
        DO p = 0, N_slave
          DO q = 0, N_slave
            DGExchange(iNbProc)%FaceDataSendU(1:1,i) = DielectricSurf(iSide)%Dielectric_dummy_Slave2(1:1,p,q)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
      ! FaceDataSendU(1:PP_nVar,1:DataSizeSideSend(iNbProc,SendID))
      CALL MPI_ISEND(DGExchange(iNbProc)%FaceDataSendU,nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
      nSendVal = PP_nVar*DataSizeSideSendMaster(iNbProc,SendID)
      DO iSide = SideID_start, SideID_end
        N_master = DG_Elems_master(iSide)
        DO p = 0, N_master
          DO q = 0, N_master
            DGExchange(iNbProc)%FaceDataSendUMaster(1:1,i) = DielectricSurf(iSide)%Dielectric_dummy_Master2(1:1,p,q)
            i = i + 1
          END DO ! q = 0, N_master
        END DO ! p = 0, N_master
      END DO ! iSide = SideID_start, SideID_end

      CALL MPI_ISEND(DGExchange(iNbProc)%FaceDataSendUMaster,nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    END IF ! SendID.EQ.2

    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIDataTypeDielectric
#endif /*not USE_HDG*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/


!==================================================================================================================================
!> Subroutine that performs the send and receive operations for the DG_Elems_slave information at the face
!> that has to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartExchange_DG_Elems(DG_Elems,LowerBound,UpperBound,SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: SendID                          !< defines the send / receive direction -> 1=send MINE/receive YOUR,
                                                         !< 2=send YOUR / receive MINE
INTEGER,INTENT(IN)    :: LowerBound                      !< lower side index for last dimension of DG_Elems
INTEGER,INTENT(IN)    :: UpperBound                      !< upper side index for last dimension of DG_Elems
TYPE(MPI_Request),INTENT(OUT)   :: SendRequest(nNbProcs)           !< communicatio handles for send
TYPE(MPI_Request),INTENT(OUT)   :: RecRequest(nNbProcs)            !< communicatio handles for receive
INTEGER,INTENT(INOUT) :: DG_Elems(LowerBound:UpperBound) !< information about DG_Elems at faces to be communicated
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal     = nMPISides_send(iNbProc,SendID)
    SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(DG_Elems(SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,SendRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    SendRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal      = nMPISides_rec(iNbProc,SendID)
    SideID_start = OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(DG_Elems(SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,RecRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    RecRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartExchange_DG_Elems


!===================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!> SendRequest, RecRequest: communication handles
!> SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
SUBROUTINE FinishExchangeMPIData(SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(INOUT) :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check receive operations first
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(RecRequest(iNbProc) ,MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(2) = MPIW8TimeField(2) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(2) = MPIW8CountField(2) + 1_8
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(SendRequest(iNbProc),MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(1) = MPIW8TimeField(1) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(1) = MPIW8CountField(1) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE FinishExchangeMPIData


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors.
!===================================================================================================================================
SUBROUTINE StartReceiveMPIDataInt(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID                                                 !< defines the send / receive direction -> 1=send MINE
                                                                              !< / receive YOUR, 3=send YOUR / receive MINE
INTEGER,INTENT(IN)  :: firstDim                                               !< size of one entry in array (e.g. one side:
                                                                              !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)  :: LowerBound                                             !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)  :: UpperBound                                             !< upper side index for last dimension of FaceData
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
INTEGER,INTENT(OUT) :: FaceData(firstDim,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =firstDim*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIDataInt


!===================================================================================================================================
!> See above, but for for send direction
!===================================================================================================================================
SUBROUTINE StartSendMPIDataInt(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
INTEGER, INTENT(IN)          :: firstDim,LowerBound,UpperBound
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(OUT)         :: MPIRequest(nNbProcs)
INTEGER, INTENT(IN)          :: FaceData(firstDim,LowerBound:UpperBound)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =firstDim*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIDataInt



#if USE_HDG
#endif /*USE_HDG*/


#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
!===================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!> SendRequest, RecRequest: communication handles
!> SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
SUBROUTINE FinishExchangeMPIDataType(SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_DG_Vars  ,ONLY: U_Surf_N,DG_Elems_slave
#if !(USE_HDG)
USE MOD_PML_Vars ,ONLY: PMLnVar,DoPML
#endif /*!(USE_HDG)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(INOUT)       :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
INTEGER                       :: i,p,q,iSide,N_slave
!===================================================================================================================================
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check receive operations first
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(RecRequest(iNbProc) ,MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(2) = MPIW8TimeField(2) + REAL(CounterEnd-CounterStart,8)/Rate
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(SendRequest(iNbProc),MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(1) = MPIW8TimeField(1) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/

! Unroll data
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    SideID_start = OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_rec(iNbProc,SendID)

    i = 1
    IF(SendID.EQ.2)THEN
      DO iSide = SideID_start, SideID_end
        N_slave = DG_Elems_slave(iSide)
        DO p = 0, N_slave
          DO q = 0, N_slave
            U_Surf_N(iSide)%U_Slave(1:PP_nVar,p,q) = DGExchange(iNbProc)%FaceDataRecvU(1:PP_nVar,i)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
    ELSE
      IF(DoPML)THEN
        DO iSide = SideID_start, SideID_end
          N_slave = DG_Elems_slave(iSide)
          DO p = 0, N_slave
            DO q = 0, N_slave
              U_Surf_N(iSide)%Flux_Slave(1:PP_nVar+PMLnVar,p,q) = DGExchange(iNbProc)%FaceDataRecvFlux(1:PP_nVar+PMLnVar,i)
              i = i + 1
            END DO ! q = 0, N_slave
          END DO ! p = 0, N_slave
        END DO ! iSide = SideID_start, SideID_end
      ELSE
        ! no PML
        DO iSide = SideID_start, SideID_end
          N_slave = DG_Elems_slave(iSide)
          DO p = 0, N_slave
            DO q = 0, N_slave
              U_Surf_N(iSide)%Flux_Slave(1:PP_nVar,p,q) = DGExchange(iNbProc)%FaceDataRecvU(1:PP_nVar,i)
              i = i + 1
            END DO ! q = 0, N_slave
          END DO ! p = 0, N_slave
        END DO ! iSide = SideID_start, SideID_end
      END IF ! DoPML
    END IF ! SendID.EQ.2

  END IF
END DO !iProc=1,nNBProcs

END SUBROUTINE FinishExchangeMPIDataType


#if !(USE_HDG)
!===================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!> SendRequest, RecRequest: communication handles
!> SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
SUBROUTINE FinishExchangeMPIDataTypeDielectric(SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_DG_Vars         ,ONLY: DG_Elems_slave,DG_Elems_master
USE MOD_Dielectric_Vars ,ONLY: DielectricSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(INOUT)       :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
INTEGER                       :: i,p,q,iSide,N_slave,N_master
!===================================================================================================================================
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check receive operations first
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(RecRequest(iNbProc) ,MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(2) = MPIW8TimeField(2) + REAL(CounterEnd-CounterStart,8)/Rate
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(SendRequest(iNbProc),MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(1) = MPIW8TimeField(1) + REAL(CounterEnd-CounterStart,8)/Rate
#endif /*defined(MEASURE_MPI_WAIT)*/

! Unroll data
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    SideID_start = OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_rec(iNbProc,SendID)

    i = 1
    IF(SendID.EQ.2)THEN
      DO iSide = SideID_start, SideID_end
        N_slave = DG_Elems_slave(iSide)
        DO p = 0, N_slave
          DO q = 0, N_slave
            DielectricSurf(iSide)%Dielectric_dummy_Slave2(1:1,p,q) = DGExchange(iNbProc)%FaceDataRecvU(1:1,i)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
    ELSE
      DO iSide = SideID_start, SideID_end
        N_master = DG_Elems_master(iSide)
        DO p = 0, N_master
          DO q = 0, N_master
            DielectricSurf(iSide)%Dielectric_dummy_Master2(1:1,p,q) = DGExchange(iNbProc)%FaceDataRecvUMaster(1:1,i)
            i = i + 1
          END DO ! q = 0, N_master
        END DO ! p = 0, N_master
      END DO ! iSide = SideID_start, SideID_end
    END IF ! SendID.EQ.2

  END IF
END DO !iProc=1,nNBProcs

END SUBROUTINE FinishExchangeMPIDataTypeDielectric
#endif /*not USE_HDG*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/


!----------------------------------------------------------------------------------------------------------------------------------!
!> Finalize DG MPI-Stuff, deallocate arrays with neighbor connections, etc.
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeMPI()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_FV
USE MOD_MPI_FV             ,ONLY: FinalizeMPIFV
#endif /*USE_FV*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(nMPISides_Proc)
SDEALLOCATE(nMPISides_MINE_Proc)
SDEALLOCATE(nMPISides_YOUR_Proc)
SDEALLOCATE(offsetMPISides_MINE)
SDEALLOCATE(offsetMPISides_YOUR)
SDEALLOCATE(NbProc)

! requires knowledge of number of MPI neighbors
#if USE_FV
CALL FinalizeMPIFV()
#endif /*USE_FV*/
SDEALLOCATE(SendRequest_U)
SDEALLOCATE(SendRequest_U2)
SDEALLOCATE(SendRequest_Flux)
SDEALLOCATE(SendRequest_GEO)
SDEALLOCATE(RecRequest_Geo)
SDEALLOCATE(RecRequest_U)
SDEALLOCATE(RecRequest_U2)
SDEALLOCATE(RecRequest_Flux)
SDEALLOCATE(nMPISides_send)
SDEALLOCATE(OffsetMPISides_send)
SDEALLOCATE(nMPISides_rec)
SDEALLOCATE(OffsetMPISides_rec)

! Free the communicators
IERROR=MPI_SUCCESS
IF(MPI_COMM_NODE   .NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_NODE   ,IERROR)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_FREE',iError)
IF(MPI_COMM_LEADERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_LEADERS,IERROR)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_FREE',iError)
IF(MPI_COMM_WORKERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_WORKERS,IERROR)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_FREE',iError)

#if USE_LOADBALANCE
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(offsetElemMPI)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

END SUBROUTINE FinalizeMPI
#endif /*USE_MPI*/

END MODULE MOD_MPI