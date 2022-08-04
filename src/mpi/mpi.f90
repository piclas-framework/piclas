!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_MPI
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitMPI
  MODULE PROCEDURE InitMPI
END INTERFACE

PUBLIC::InitMPI

#if USE_MPI
INTERFACE InitMPIvars
  MODULE PROCEDURE InitMPIvars
END INTERFACE

INTERFACE FinishExchangeMPIData
  MODULE PROCEDURE FinishExchangeMPIData
END INTERFACE

INTERFACE FinalizeMPI
  MODULE PROCEDURE FinalizeMPI
END INTERFACE


PUBLIC::InitMPIvars,StartReceiveMPIData,StartSendMPIData,StartReceiveMPIDataInt,StartSendMPIDataInt,FinishExchangeMPIData,FinalizeMPI
#endif
PUBLIC::DefineParametersMPI
#if defined(MEASURE_MPI_WAIT)
PUBLIC::OutputMPIW8Time
#endif /*defined(MEASURE_MPI_WAIT)*/
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
INTEGER,INTENT(IN),OPTIONAL      :: mpi_comm_IN !< MPI communicator
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER :: MPI_COMM_LOC
LOGICAL :: initDone
!==================================================================================================================================
IF (PRESENT(mpi_comm_IN)) THEN
  MPI_COMM_LOC = mpi_comm_IN
ELSE
  CALL MPI_INIT(iError)
  CALL MPI_INITIALIZED(initDone,iError)
  IF(.NOT.initDone) CALL MPI_INIT(iError)
  IF(iError .NE. 0) &
    CALL Abort(__STAMP__,'Error in MPI_INIT',iError)
  MPI_COMM_LOC = MPI_COMM_WORLD
END IF

CALL MPI_COMM_RANK(MPI_COMM_LOC, myRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_LOC, nProcessors, iError)
IF(iError .NE. 0) &
  CALL Abort(&
  __STAMP__&
  ,'Could not get rank and number of processors',iError)
MPIRoot=(myRank .EQ. 0)
#else  /*USE_MPI*/
myRank      = 0
myLocalRank = 0
nProcessors = 1
MPIRoot     =.TRUE.
MPILocalRoot=.TRUE.
#endif  /*USE_MPI*/

! At this point the initialization is not completed. We first have to create a new MPI communicator. MPIInitIsDone will be set
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
ALLOCATE(SendRequest_U(nNbProcs)     )
ALLOCATE(SendRequest_U2(nNbProcs)     )
ALLOCATE(SendRequest_GEO(nNbProcs)     )
!ALLOCATE(SendRequest_UMinus(nNbProcs)     )
ALLOCATE(SendRequest_Flux(nNbProcs)  )
ALLOCATE(SendRequest_gradUx(nNbProcs))
ALLOCATE(SendRequest_gradUy(nNbProcs))
ALLOCATE(SendRequest_gradUz(nNbProcs))
ALLOCATE(RecRequest_U(nNbProcs)     )
ALLOCATE(RecRequest_U2(nNbProcs)     )
ALLOCATE(RecRequest_Geo(nNbProcs)     )
!ALLOCATE(RecRequest_UMinus(nNbProcs)     )
ALLOCATE(RecRequest_Flux(nNbProcs)  )
ALLOCATE(RecRequest_gradUx(nNbProcs))
ALLOCATE(RecRequest_gradUy(nNbProcs))
ALLOCATE(RecRequest_gradUz(nNbProcs))
SendRequest_U(nNbProcs)      = MPI_REQUEST_NULL
SendRequest_U2(nNbProcs)      = MPI_REQUEST_NULL
!SendRequest_UMinus           = MPI_REQUEST_NULL
SendRequest_Flux(nNbProcs)   = MPI_REQUEST_NULL
SendRequest_gradUx(nNbProcs) = MPI_REQUEST_NULL
SendRequest_gradUy(nNbProcs) = MPI_REQUEST_NULL
SendRequest_gradUz(nNbProcs) = MPI_REQUEST_NULL
RecRequest_U(nNbProcs)       = MPI_REQUEST_NULL
RecRequest_U2(nNbProcs)       = MPI_REQUEST_NULL
!RecRequest_UMinus            = MPI_REQUEST_NULL
RecRequest_Flux(nNbProcs)    = MPI_REQUEST_NULL
RecRequest_gradUx(nNbProcs)  = MPI_REQUEST_NULL
RecRequest_gradUy(nNbProcs)  = MPI_REQUEST_NULL
RecRequest_gradUz(nNbProcs)  = MPI_REQUEST_NULL
SendRequest_Geo(nNbProcs)    = MPI_REQUEST_NULL
RecRequest_Geo(nNbProcs)     = MPI_REQUEST_NULL
DataSizeSide  =(PP_N+1)*(PP_N+1)

! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize','0')
IF(GroupSize.LT.1)THEN ! group procs by node
  ! Split the node communicator (shared memory) from the global communicator on physical processor or node level
#if (CORE_SPLIT==1)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myRank,0,MPI_COMM_NODE,iError)
#elif (CORE_SPLIT==0)
  ! Note that using SharedMemoryMethod=OMPI_COMM_TYPE_CORE somehow does not work in every case (intel/amd processors)
  ! Also note that OMPI_COMM_TYPE_CORE is undefined when not using OpenMPI
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,SharedMemoryMethod,0,MPI_INFO_NULL,MPI_COMM_NODE,IERROR)
#else
  ! Check if more nodes than procs are required or
  ! if the resulting split would create unequal procs per node
  IF((CORE_SPLIT.GE.nProcessors_Global).OR.(MOD(nProcessors_Global,CORE_SPLIT).GT.0))THEN
    SWRITE (*,'(A,I0,A,I0,A,F0.2,A)') ' WARNING: Either more nodes than cores selected (nodes: ',CORE_SPLIT,', cores: ',&
        nProcessors_Global,') OR unequal number of cores per node (=',REAL(nProcessors_Global)/REAL(CORE_SPLIT),&
        '). Setting 1 core per node for MPI_COMM_NODE!'
    color = myRank
  ELSE
    ! Group procs so that every CORE_SPLIT procs are in the same group
    color = INT(REAL(myrank*CORE_SPLIT)/REAL(nProcessors_Global))+1
  END IF ! (CORE_SPLIT.GE.nProcessors_Global).OR.(MOD().GT.0)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,MPI_COMM_NODE,iError)
#endif
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,MPI_COMM_NODE,iError)
END IF
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
MPILocalRoot=(myLocalRank.EQ.0)

IF (nProcessors_Global.EQ.nLocalProcs) THEN
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs in ',1,' group'
ELSE
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs each in ',&
                                                        nProcessors_Global/nLocalProcs,' groups for a total number of ',&
                                                        nProcessors_Global,' procs'
END IF

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,0,0,MPI_COMM_LEADERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,1,0,MPI_COMM_WORKERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
  nLeaderProcs=nProcessors-nWorkerProcs
END IF
END SUBROUTINE InitMPIvars


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors.
!===================================================================================================================================
SUBROUTINE StartReceiveMPIData(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
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
INTEGER,INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
REAL,INTENT(OUT)    :: FaceData(firstDim,0:PP_N,0:PP_N,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =firstDim*DataSizeSide*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIData


!===================================================================================================================================
!> See above, but for for send direction
!===================================================================================================================================
SUBROUTINE StartSendMPIData(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
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
INTEGER, INTENT(OUT)         :: MPIRequest(nNbProcs)
REAL, INTENT(IN)             :: FaceData(firstDim,0:PP_N,0:PP_N,LowerBound:UpperBound)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =firstDim*DataSizeSide*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIData


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
INTEGER, INTENT(INOUT)       :: SendRequest(nNbProcs),RecRequest(nNbProcs)
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
  IF(nMPISides_rec(iNbProc,SendID).GT.0) CALL MPI_WAIT(RecRequest(iNbProc) ,MPIStatus,iError)
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(2) = MPIW8TimeField(2) + REAL(CounterEnd-CounterStart,8)/Rate
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0) CALL MPI_WAIT(SendRequest(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  ! Note: Send and Receive are switched to have the same ordering as for particles (1. Send, 2. Receive)
  MPIW8TimeField(1) = MPIW8TimeField(1) + REAL(CounterEnd-CounterStart,8)/Rate
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
  INTEGER,INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
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
                      nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
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
  INTEGER, INTENT(OUT)         :: MPIRequest(nNbProcs)
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
                      nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
    ELSE
      MPIRequest(iNbProc)=MPI_REQUEST_NULL
    END IF
  END DO !iProc=1,nNBProcs
  END SUBROUTINE StartSendMPIDataInt

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
SDEALLOCATE(SendRequest_U)
SDEALLOCATE(SendRequest_U2)
SDEALLOCATE(SendRequest_Flux)
SDEALLOCATE(SendRequest_GEO)
SDEALLOCATE(RecRequest_Geo)
!ALLOCATE(SendRequest_UMinus(nNbProcs)     )
SDEALLOCATE(SendRequest_gradUx)
SDEALLOCATE(SendRequest_gradUy)
SDEALLOCATE(SendRequest_gradUz)
SDEALLOCATE(RecRequest_U)
SDEALLOCATE(RecRequest_U2)
SDEALLOCATE(RecRequest_Flux)
SDEALLOCATE(RecRequest_gradUx)
SDEALLOCATE(RecRequest_gradUy)
SDEALLOCATE(RecRequest_gradUz)
SDEALLOCATE(nMPISides_send)
SDEALLOCATE(OffsetMPISides_send)
SDEALLOCATE(nMPISides_rec)
SDEALLOCATE(OffsetMPISides_rec)

! Free the communicators
IF(MPI_COMM_NODE   .NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_NODE   ,IERROR)
IF(MPI_COMM_LEADERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_LEADERS,IERROR)

#if USE_LOADBALANCE
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(offsetElemMPI)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

END SUBROUTINE FinalizeMPI
#endif /*USE_MPI*/


#if defined(MEASURE_MPI_WAIT)
!===================================================================================================================================
!> Root writes measured MPI_WAIT() times to disk
!===================================================================================================================================
SUBROUTINE OutputMPIW8Time()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars          ,ONLY: MPIW8TimeGlobal,MPIW8TimeSim,MPIW8TimeProc,MPIW8TimeField,MPIW8Time,MPIW8TimeGlobal,MPIW8TimeBaS
USE MOD_StringTools       ,ONLY: INTTOSTR
#if defined(PARTICLES)
USE MOD_Particle_MPI_Vars ,ONLY: MPIW8TimePart
#endif /*defined(PARTICLES)*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                                :: WriteHeader
INTEGER                                :: ioUnit,i
CHARACTER(LEN=150)                     :: formatStr
CHARACTER(LEN=22),PARAMETER            :: outfile    ='MPIW8Time.csv'
CHARACTER(LEN=22),PARAMETER            :: outfilePerc='MPIW8TimePercent.csv'
CHARACTER(LEN=22),PARAMETER            :: outfileProc='MPIW8TimeProc'
CHARACTER(LEN=30)                      :: outfileProc_loc
CHARACTER(LEN=10)                      :: hilf
INTEGER,PARAMETER                      :: nTotalVars =MPIW8SIZE+2
CHARACTER(LEN=255),DIMENSION(nTotalVars) :: StrVarNames(nTotalVars)=(/ CHARACTER(LEN=255) :: &
    'nProcessors'       , &
    'WallTimeSim'       , &
    'Barrier-and-Sync'    &
#if USE_HDG
   ,'HDG-SendLambda'    , & ! (1)
    'HDG-ReceiveLambda' , & ! (2)
    'HDG-Broadcast'     , & ! (3)
    'HDG-Allreduce'       & ! (4)
#else
   ,'DGSEM-Send'    , &     ! (1)
    'DGSEM-Receive'   &     ! (2)
#endif /*USE_HDG*/
#if defined(PARTICLES)
   ,'SendNbrOfParticles'  , & ! (1)
    'RecvNbrOfParticles'  , & ! (2)
    'SendParticles'       , & ! (3)
    'RecvParticles'       , & ! (4)
    'EmissionParticles'   , & ! (5)
    'PIC-depo-Wait'         & ! (6)
#endif /*defined(PARTICLES)*/
    /)
! CHARACTER(LEN=255),DIMENSION(nTotalVars) :: StrVarNamesProc(nTotalVars)=(/ CHARACTER(LEN=255) :: &
!     'Rank'                  &
! #if defined(PARTICLES)
!    ,'SendNbrOfParticles'  , &
!     'RecvNbrOfParticles'  , &
!     'SendParticles'       , &
!     'RecvParticles'         &
! #endif /*defined(PARTICLES)*/
!     /)
CHARACTER(LEN=255)         :: tmpStr(nTotalVars)
CHARACTER(LEN=1000)        :: tmpStr2
CHARACTER(LEN=1),PARAMETER :: delimiter=","
REAL                       :: TotalSimTime
!===================================================================================================================================
MPIW8Time(               1:1)                              = MPIW8TimeBaS
MPIW8Time(               2:MPIW8SIZEFIELD+1)               = MPIW8TimeField
#if defined(PARTICLES)
MPIW8Time(MPIW8SIZEFIELD+2:MPIW8SIZEFIELD+MPIW8SIZEPART+1) = MPIW8TimePart
#endif /*defined(PARTICLES)*/

! Collect and output measured MPI_WAIT() times
IF(MPIroot)THEN
  ALLOCATE(MPIW8TimeProc(MPIW8SIZE*nProcessors))
  CALL MPI_REDUCE(MPIW8Time , MPIW8TimeGlobal , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_WORLD , iError)
  CALL MPI_GATHER(MPIW8Time , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPIW8TimeProc , MPIW8SIZE , MPI_DOUBLE_PRECISION , 0 , MPI_COMM_WORLD , iError)
ELSE
  CALL MPI_REDUCE(MPIW8Time , 0               , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_WORLD , IError)
  CALL MPI_GATHER(MPIW8Time , MPIW8SIZE , MPI_DOUBLE_PRECISION , 0             , 0         , 0                    , 0 , MPI_COMM_WORLD , iError)
END IF

! --------------------------------------------------
! Only MPI root outputs the data to file
! --------------------------------------------------
IF(.NOT.MPIRoot)RETURN

! --------------------------------------------------
! Output file 1 of 3: MPIW8Time.csv
! --------------------------------------------------
! Either create new file or add info to existing file
WriteHeader = .TRUE.
IF(FILEEXISTS(outfile)) WriteHeader = .FALSE.

IF(WriteHeader)THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=' '
  DO i=1,nTotalVars
    WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO i=1,nTotalVars
    IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF ! WriteHeader

IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
  WRITE(tmpStr2,formatStr)&
      " ",REAL(nProcessors)       ,&
      delimiter,MPIW8TimeSim      ,&
      delimiter,MPIW8TimeGlobal(1),&
      delimiter,MPIW8TimeGlobal(2),&
      delimiter,MPIW8TimeGlobal(3) &
#if USE_HDG
     ,delimiter,MPIW8TimeGlobal(4),&
      delimiter,MPIW8TimeGlobal(5) &
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+1),&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+2),&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+3),&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+4),&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+5),&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+6) &
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  WRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write MPI_WAIT wall time info!"
END IF

! --------------------------------------------------
! Output file 2 of 3: MPIW8Time.csv
! --------------------------------------------------
! Either create new file or add info to existing file
WriteHeader = .TRUE.
IF(FILEEXISTS(outfilePerc)) WriteHeader = .FALSE.

IF(WriteHeader)THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfilePerc),STATUS="UNKNOWN")
  tmpStr=""
  DO i=1,nTotalVars
    WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO i=1,nTotalVars
    IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF ! WriteHeader

IF(FILEEXISTS(outfilePerc))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfilePerc),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
  TotalSimTime = MPIW8TimeSim*REAL(nProcessors)
  WRITE(tmpStr2,formatStr)&
      " ",REAL(nProcessors)       ,&
      delimiter,100.              ,& ! MPIW8TimeSim*nProcessors / TotalSimTime
      delimiter,MPIW8TimeGlobal(1)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(2)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(3)/TotalSimTime &
#if USE_HDG
     ,delimiter,MPIW8TimeGlobal(4)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(5)/TotalSimTime &
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+1)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+2)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+3)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+4)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+5)/TotalSimTime,&
      delimiter,MPIW8TimeGlobal(MPIW8SIZEFIELD+1+6)/TotalSimTime &
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  WRITE(UNIT_StdOut,'(A)')TRIM(outfilePerc)//" does not exist. Cannot write MPI_WAIT wall time info!"
END IF

! --------------------------------------------------
! Output file 3 of 3: MPIW8TimeProc-XXX.csv
! --------------------------------------------------
! Cannot append to proc file, iterate name  -000.csv, -001.csv, -002.csv
WRITE(UNIT=hilf,FMT='(A1,I3.3,A4)') '-',0,'.csv'
outfileProc_loc=TRIM(outfileProc)//'-'//TRIM(ADJUSTL(INTTOSTR(nProcessors)))
DO WHILE(FILEEXISTS(TRIM(outfileProc_loc)//TRIM(hilf)))
  i = i + 1
  WRITE(UNIT=hilf,FMT='(A1,I3.3,A4)') '-',i,'.csv'
END DO
outfileProc_loc=TRIM(outfileProc_loc)//TRIM(hilf)

! Write the file header
OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfileProc_loc),STATUS="UNKNOWN")
tmpStr=' '
DO i=1,nTotalVars
  WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
END DO
WRITE(formatStr,'(A1)')'('
DO i=1,nTotalVars
  IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
    WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
  ELSE
    WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
  END IF
END DO

WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

! Output the processor wait times
WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
DO i = 0,nProcessors-1
  WRITE(tmpStr2,formatStr)&
            " ",REAL(i)                     ,&
      delimiter,MPIW8TimeSim                ,&
      delimiter,MPIW8TimeProc(i*MPIW8SIZE+1),&
      delimiter,MPIW8TimeProc(i*MPIW8SIZE+2),&
      delimiter,MPIW8TimeProc(i*MPIW8SIZE+3) &
#if USE_HDG
     ,delimiter,MPIW8TimeProc(i*MPIW8SIZE+4),&
      delimiter,MPIW8TimeProc(i*MPIW8SIZE+5) &
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+1),&
      delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+2),&
      delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+3),&
      delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+4),&
      delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+5),&
      delimiter,MPIW8TimeProc(MPIW8SIZEFIELD+1+i*MPIW8SIZE+6) &
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
END DO
CLOSE(ioUnit)

DEALLOCATE(MPIW8TimeProc)

END SUBROUTINE OutputMPIW8Time
#endif /*defined(MEASURE_MPI_WAIT)*/

END MODULE MOD_MPI
