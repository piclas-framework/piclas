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

! don't create an interface because some vectors are mapped to arrays
!INTERFACE StartReceiveMPIData
!  MODULE PROCEDURE StartReceiveMPIData
!END INTERFACE

!INTERFACE StartSendMPIData
!  MODULE PROCEDURE StartSendMPIData
!END INTERFACE


INTERFACE FinishExchangeMPIData
  MODULE PROCEDURE FinishExchangeMPIData
END INTERFACE

INTERFACE FinalizeMPI
  MODULE PROCEDURE FinalizeMPI
END INTERFACE


PUBLIC::InitMPIvars,StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData,FinalizeMPI
#endif
PUBLIC::DefineParametersMPI
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMPI()
! MODULES
USE MOD_ReadInTools,              ONLY: prms
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
END SUBROUTINE DefineParametersMPI

SUBROUTINE InitMPI(mpi_comm_IN)
!===================================================================================================================================
! Basic MPI initialization.
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL      :: mpi_comm_IN !< MPI communicator
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
SUBROUTINE InitMPIvars()
!===================================================================================================================================
! Initialize derived MPI types used for communication and allocate HALO data.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
USE MOD_Readintools,       ONLY:GETINT
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
IF(.NOT.InterpolationInitIsDone)THEN
  CALL Abort(&
      __STAMP__&
      ,'InitMPITypes called before InitInterpolation')
END IF
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
! currenlty allocated in prepare_mesh
!ALLOCATE(nMPISides_send(       nNbProcs,2))
!ALLOCATE(OffsetMPISides_send(0:nNbProcs,2))
!ALLOCATE(nMPISides_rec(        nNbProcs,2))
!ALLOCATE(OffsetMPISides_rec( 0:nNbProcs,2))
! Set number of sides and offset for SEND MINE - RECEIVE YOUR case
!nMPISides_send(:,1)     =nMPISides_MINE_Proc
!OffsetMPISides_send(:,1)=OffsetMPISides_MINE
!nMPISides_rec(:,1)      =nMPISides_YOUR_Proc
!OffsetMPISides_rec(:,1) =OffsetMPISides_YOUR
!! Set number of sides and offset for SEND YOUR - RECEIVE MINE case
!nMPISides_send(:,2)     =nMPISides_YOUR_Proc
!OffsetMPISides_send(:,2)=OffsetMPISides_YOUR
!nMPISides_rec(:,2)      =nMPISides_MINE_Proc
!OffsetMPISides_rec(:,2) =OffsetMPISides_MINE


! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize','0')
IF(GroupSize.LT.1)THEN ! group procs by node
#if USE_MPI3
  ! MPI3 directly gives you shared memory groups
  CALL MPI_INFO_CREATE(info,iError)
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,myRank,info,MPI_COMM_NODE,iError)
#else
  ! TODO: Build own node communicator
  !CALL MPI_Get_processor_name(procname,length,iError)
  ! Now generate hash from string
  ! TODO: find good hash function, but beware hash collisions may occur!
  ! Maybe use two different hash functions and check if resulting groups are identical
  !CALL GetHashesFromString(procname,color1,color2) ! beware, hash collisions may occur, check results TODO: find good hash function
  !CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color1,myRank,MPI_COMM_NODE1,iError)
  !CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color2,myRank,MPI_COMM_NODE2,iError)
  ! Compare ...

  ! Fallback:
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myRank,myRank,MPI_COMM_NODE,iError)
#endif
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myRank,MPI_COMM_NODE,iError)
END IF
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
MPILocalRoot=(myLocalRank .EQ. 0)

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,0,myRank,MPI_COMM_LEADERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,1,myRank,MPI_COMM_WORKERS,iError)
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


SUBROUTINE StartSendMPIData(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
!===================================================================================================================================
! See above, but for for send direction
!===================================================================================================================================
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



SUBROUTINE FinishExchangeMPIData(SendRequest,RecRequest,SendID)
!===================================================================================================================================
! We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
! SendRequest, RecRequest: communication handles
! SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
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
!===================================================================================================================================
! Check receive operations first
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0) CALL MPI_WAIT(RecRequest(iNbProc) ,MPIStatus,iError)
END DO !iProc=1,nNBProcs
! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0) CALL MPI_WAIT(SendRequest(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs
END SUBROUTINE FinishExchangeMPIData


SUBROUTINE FinalizeMPI()
!----------------------------------------------------------------------------------------------------------------------------------!
! Finalize DG MPI-Stuff, deallocate arrays with neighbor connections, etc.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(offsetElemMPI)
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



END SUBROUTINE FinalizeMPI
#endif /*USE_MPI*/

END MODULE MOD_MPI
