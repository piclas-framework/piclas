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

MODULE MOD_MPI_FV
!===================================================================================================================================
! Module containing subroutines and functions that are required only for Finite Volume (FV) methods
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi_f08
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI
#if USE_FV
PUBLIC::InitMPIvarsFV,StartReceiveMPIDataFV,StartSendMPIDataFV,FinalizeMPIFV
#endif /*USE_FV*/
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI
#if USE_FV
!===================================================================================================================================
!> Allocates the required send/receive buffers for FV
!===================================================================================================================================
SUBROUTINE InitMPIvarsFV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(SendRequest_gradUx(nNbProcs))
ALLOCATE(SendRequest_gradUy(nNbProcs))
ALLOCATE(SendRequest_gradUz(nNbProcs))
ALLOCATE(RecRequest_gradUx(nNbProcs))
ALLOCATE(RecRequest_gradUy(nNbProcs))
ALLOCATE(RecRequest_gradUz(nNbProcs))
SendRequest_gradUx = MPI_REQUEST_NULL
SendRequest_gradUy = MPI_REQUEST_NULL
SendRequest_gradUz = MPI_REQUEST_NULL
RecRequest_gradUx  = MPI_REQUEST_NULL
RecRequest_gradUy  = MPI_REQUEST_NULL
RecRequest_gradUz  = MPI_REQUEST_NULL
END SUBROUTINE InitMPIvarsFV


!===================================================================================================================================
!> Subroutine does the receive operations for the single value face data that has to be exchanged between processors.
!===================================================================================================================================
SUBROUTINE StartReceiveMPIDataFV(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
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
END SUBROUTINE StartReceiveMPIDataFV


!===================================================================================================================================
!> See above, but for for send direction
!===================================================================================================================================
SUBROUTINE StartSendMPIDataFV(firstDim,FaceData,LowerBound,UpperBound,MPIRequest,SendID)
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
TYPE(MPI_Request),INTENT(OUT)         :: MPIRequest(nNbProcs)
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
END SUBROUTINE StartSendMPIDataFV


!===================================================================================================================================
!> Deallocates the required send/receive buffers for FV
!===================================================================================================================================
SUBROUTINE FinalizeMPIFV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(SendRequest_gradUx)
SDEALLOCATE(SendRequest_gradUy)
SDEALLOCATE(SendRequest_gradUz)
SDEALLOCATE(RecRequest_gradUx)
SDEALLOCATE(RecRequest_gradUy)
SDEALLOCATE(RecRequest_gradUz)
END SUBROUTINE FinalizeMPIFV
#endif /*USE_FV*/
#endif /*USE_MPI*/

END MODULE MOD_MPI_FV
