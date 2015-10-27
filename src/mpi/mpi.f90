#include "boltzplatz.h"

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

#ifdef MPI
INTERFACE InitMPIvars
  MODULE PROCEDURE InitMPIvars
END INTERFACE

INTERFACE StartReceiveMPIData
  MODULE PROCEDURE StartReceiveMPIData
END INTERFACE

INTERFACE StartSendMPIData
  MODULE PROCEDURE StartSendMPIData
END INTERFACE

INTERFACE FinishExchangeMPIData
  MODULE PROCEDURE FinishExchangeMPIData
END INTERFACE

INTERFACE FinalizeMPI
  MODULE PROCEDURE FinalizeMPI
END INTERFACE


PUBLIC::InitMPIvars,StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData, FinalizeMPI
#endif
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMPI()
!===================================================================================================================================
! Basic MPI initialization. 
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
#ifdef MPI
CALL MPI_INIT(iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'Error in MPI_INIT',iError,999.)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProcessors, iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'Could not get rank and number of processors',iError,999.)
MPIRoot=(myRank .EQ. 0)
#else  /*MPI*/
myRank      = 0 
nProcessors = 1 
MPIRoot     =.TRUE.
#endif  /*MPI*/

! At this point the initialization is not completed. We first have to create a new MPI communicator. MPIInitIsDone will be set
END SUBROUTINE InitMPI



#ifdef MPI
SUBROUTINE InitMPIvars()
!===================================================================================================================================
! Initialize derived MPI types used for communication and allocate HALO data. 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.InterpolationInitIsDone)THEN
  CALL Abort(__STAMP__,'InitMPITypes called before InitInterpolation')
END IF
ALLOCATE(SendRequest_U(nNbProcs)     )
!ALLOCATE(SendRequest_UMinus(nNbProcs)     )
ALLOCATE(SendRequest_Flux(nNbProcs)  )
ALLOCATE(SendRequest_gradUx(nNbProcs))
ALLOCATE(SendRequest_gradUy(nNbProcs))
ALLOCATE(SendRequest_gradUz(nNbProcs))
ALLOCATE(RecRequest_U(nNbProcs)     )
!ALLOCATE(RecRequest_UMinus(nNbProcs)     )
ALLOCATE(RecRequest_Flux(nNbProcs)  )
ALLOCATE(RecRequest_gradUx(nNbProcs))
ALLOCATE(RecRequest_gradUy(nNbProcs))
ALLOCATE(RecRequest_gradUz(nNbProcs))
SendRequest_U(nNbProcs)      = MPI_REQUEST_NULL
!SendRequest_UMinus           = MPI_REQUEST_NULL
SendRequest_Flux(nNbProcs)   = MPI_REQUEST_NULL
SendRequest_gradUx(nNbProcs) = MPI_REQUEST_NULL
SendRequest_gradUy(nNbProcs) = MPI_REQUEST_NULL
SendRequest_gradUz(nNbProcs) = MPI_REQUEST_NULL
RecRequest_U(nNbProcs)       = MPI_REQUEST_NULL
!RecRequest_UMinus            = MPI_REQUEST_NULL
RecRequest_Flux(nNbProcs)    = MPI_REQUEST_NULL
RecRequest_gradUx(nNbProcs)  = MPI_REQUEST_NULL
RecRequest_gradUy(nNbProcs)  = MPI_REQUEST_NULL
RecRequest_gradUz(nNbProcs)  = MPI_REQUEST_NULL
DataSizeSide  =(PP_N+1)*(PP_N+1)
ALLOCATE(nMPISides_send(       nNbProcs,2))
ALLOCATE(OffsetMPISides_send(0:nNbProcs,2))
ALLOCATE(nMPISides_rec(        nNbProcs,2))
ALLOCATE(OffsetMPISides_rec( 0:nNbProcs,2))
! Set number of sides and offset for SEND MINE - RECEIVE YOUR case
nMPISides_send(:,1)     =nMPISides_MINE_Proc
OffsetMPISides_send(:,1)=OffsetMPISides_MINE
nMPISides_rec(:,1)      =nMPISides_YOUR_Proc
OffsetMPISides_rec(:,1) =OffsetMPISides_YOUR
! Set number of sides and offset for SEND YOUR - RECEIVE MINE case
nMPISides_send(:,2)     =nMPISides_YOUR_Proc
OffsetMPISides_send(:,2)=OffsetMPISides_YOUR
nMPISides_rec(:,2)      =nMPISides_MINE_Proc
OffsetMPISides_rec(:,2) =OffsetMPISides_MINE
END SUBROUTINE InitMPIvars


SUBROUTINE StartReceiveMPIData(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
!===================================================================================================================================
! Subroutine does the receive operations for the face data that has to be exchanged between processors.
! FaceData: the complete face data (for inner, BC and MPI sides).
! DataSize: size of one entry in array (e.g. one side: nVar*(N+1)*(N+1))
! LowerBound / UpperBound: lower side index and upper side index for last dimension of FaceData
! MPIRequest: communication handles
! SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: SendID
INTEGER,INTENT(IN)          :: firstDim,LowerBound,UpperBound
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)         :: MPIRequest(nNbProcs)
REAL,INTENT(OUT)            :: FaceData(firstDim,0:PP_N,0:PP_N,LowerBound:UpperBound)
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
SDEALLOCATE(SendRequest_Flux)
SDEALLOCATE(SendRequest_gradUx)
SDEALLOCATE(SendRequest_gradUy)
SDEALLOCATE(SendRequest_gradUz)
SDEALLOCATE(RecRequest_U)
SDEALLOCATE(RecRequest_Flux)
SDEALLOCATE(RecRequest_gradUx)
SDEALLOCATE(RecRequest_gradUy)
SDEALLOCATE(RecRequest_gradUz)
SDEALLOCATE(nMPISides_send)
SDEALLOCATE(OffsetMPISides_send)
SDEALLOCATE(nMPISides_rec)
SDEALLOCATE(OffsetMPISides_rec)



END SUBROUTINE FinalizeMPI
#endif /*MPI*/

END MODULE MOD_MPI
