!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel_MPI
!===================================================================================================================================
!> Module for surface model mpi
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

#if USE_MPI
INTERFACE ExchangeSurfaceHaloToOrigin
  MODULE PROCEDURE ExchangeSurfaceHaloToOrigin
END INTERFACE

INTERFACE ExchangeSurfaceOriginToHalo
  MODULE PROCEDURE ExchangeSurfaceHaloToOrigin
END INTERFACE

INTERFACE MapHaloInnerToOriginInnerSurf
  MODULE PROCEDURE MapHaloInnerToOriginInnerSurf
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitSurfModel_MPI
PUBLIC :: InitSMCR_MPI
PUBLIC :: ExchangeSurfaceHaloToOrigin
PUBLIC :: ExchangeSurfaceOriginToHalo
PUBLIC :: ExchangeSurfDistInfo
PUBLIC :: MapHaloInnerToOriginInnerSurf
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI
SUBROUTINE InitSurfModel_MPI()
!===================================================================================================================================
!> Initializing MPI for Surface Model
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars      ,ONLY: SurfExchange
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfModelExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iProc
!===================================================================================================================================
! allocate information for sending from own halo sides to origin sides of neighbor procs
ALLOCATE(SurfModelExchange%nH2OSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%nH2OSidesRecv(1:SurfCOMM%nMPINeighbors) )
ALLOCATE(SurfModelExchange%H2OSendBuf(1:SurfCOMM%nMPINeighbors))
ALLOCATE(SurfModelExchange%H2ORecvBuf(1:SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  ! same number of sides as normal surface exchange      
  SurfModelExchange%nH2OSidesSend(iProc) = SurfExchange%nSidesSend(iProc)
  SurfModelExchange%nH2OSidesRecv(iProc) = SurfExchange%nSidesRecv(iProc)
  IF(SurfModelExchange%nO2HSidesSend(iProc).NE.0) THEN
    ! allocate and initialize mapping of to be send sides
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%H2OSendList(SurfModelExchange%nH2OSidesSend(iProc)))
    SurfCOMM%MPINeighbor(iProc)%H2OSendList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
    ! allocate and initialize buffer of to be send sides
    ALLOCATE(SurfModelExchange%H2OSendBuf(iProc)%content_int((nSurfSample**2)*SurfModelExchange%nH2OSidesSend(iProc)))
    SurfModelExchange%H2OSendBuf(iProc)%content_int=0
    ALLOCATE(SurfModelExchange%H2OSendBuf(iProc)%content((nSurfSample**2)*SurfModelExchange%nH2OSidesSend(iProc)))
    SurfModelExchange%H2OSendBuf(iProc)%content=0.
  END IF
  IF(SurfModelExchange%nH2OSidesRecv(iProc).NE.0) THEN
    ! allocate and initialize mapping of to be received sides
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%H2ORecvList(SurfModelExchange%nH2OSidesRecv(iProc)))
    SurfCOMM%MPINeighbor(iProc)%H2ORecvList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
    ! allocate and initialize buffer of to be received sides
    ALLOCATE(SurfModelExchange%H2ORecvBuf(iProc)%content_int((nSurfSample**2)*SurfModelExchange%nH2OSidesRecv(iProc)))
    SurfModelExchange%H2ORecvBuf(iProc)%content_int=0
    ALLOCATE(SurfModelExchange%H2ORecvBuf(iProc)%content((nSurfSample**2)*SurfModelExchange%nH2OSidesRecv(iProc)))
    SurfModelExchange%H2ORecvBuf(iProc)%content=0
  END IF
END DO

! represents inverse communication where information is send from own origin sides to halo sides of neighbor procs
ALLOCATE(SurfModelExchange%nO2HSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%nO2HSidesRecv(1:SurfCOMM%nMPINeighbors) )
ALLOCATE(SurfModelExchange%O2HSendBuf(1:SurfCOMM%nMPINeighbors))
ALLOCATE(SurfModelExchange%O2HRecvBuf(1:SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  ! inverse mapping of number of sides to communicate
  SurfModelExchange%nO2HSidesSend(iProc) = SurfModelExchange%nH2OSidesRecv(iProc)
  SurfModelExchange%nO2HSidesRecv(iProc) = SurfModelExchange%nH2OSidesSend(iProc)
  IF(SurfModelExchange%nO2HSidesSend(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%O2HSendList(SurfModelExchange%nO2HSidesSend(iProc)))
    ! inverse mapping of list of sides to send
    SurfCOMM%MPINeighbor(iProc)%O2HSendList(:)=SurfCOMM%MPINeighbor(iProc)%H2ORecvList(:)
    ALLOCATE(SurfModelExchange%O2HSendBuf(iProc)%content_int((nSurfSample**2)*SurfModelExchange%nO2HSidesSend(iProc)))
    SurfModelExchange%O2HSendBuf(iProc)%content_int=0
    ALLOCATE(SurfModelExchange%O2HSendBuf(iProc)%content((nSurfSample**2)*SurfModelExchange%nO2HSidesSend(iProc)))
    SurfModelExchange%O2HSendBuf(iProc)%content=0.
  END IF
  IF(SurfModelExchange%nO2HSidesRecv(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%O2HRecvList(SurfModelExchange%nO2HSidesRecv(iProc)))
    ! inverse mapping of list of sides to receive
    SurfCOMM%MPINeighbor(iProc)%O2HRecvList(:)=SurfCOMM%MPINeighbor(iProc)%H2OSendList(:)
    ALLOCATE(SurfModelExchange%O2HRecvBuf(iProc)%content_int((nSurfSample**2)*SurfModelExchange%nO2HSidesRecv(iProc)))
    SurfModelExchange%O2HRecvBuf(iProc)%content_int=0
    ALLOCATE(SurfModelExchange%O2HRecvBuf(iProc)%content((nSurfSample**2)*SurfModelExchange%nO2HSidesRecv(iProc)))
    SurfModelExchange%O2HRecvBuf(iProc)%content=0
  END IF
END DO

ALLOCATE(SurfModelExchange%SendRequest(1:3,SurfCOMM%nMPINeighbors)  &
        ,SurfModelExchange%RecvRequest(1:3,SurfCOMM%nMPINeighbors)  )

END SUBROUTINE InitSurfModel_MPI


SUBROUTINE ExchangeSurfaceHaloToOrigin(IntDataIN,RealDataIn,AddFlag)
!===================================================================================================================================
!> Each process sends his halo-information directly to the halo process by use of a list, containing the surfsideids for sending.
!> The receiving process adds the new data to his origin surfaces
!> Only processes with surfaces to send or receive participate in the communication
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm, nSurfSample, SurfMesh
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfModelExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL :: IntDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
REAL   ,INTENT(INOUT),OPTIONAL :: RealDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
LOGICAL,INTENT(IN)   ,OPTIONAL :: AddFlag
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: MessageSize,iSurfSide,SurfSideID
INTEGER :: iPos,p,q,iProc
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
LOGICAL :: locAddFlag
!===================================================================================================================================
IF (.NOT.ALLOCATED(SurfModelExchange%nH2OSidesSend) .OR. .NOT.ALLOCATED(SurfModelExchange%nH2OSidesRecv)) RETURN
IF (.NOT.PRESENT(IntDataIN) .AND. .NOT.PRESENT(RealDataIN) ) RETURN
IF (PRESENT(AddFlag)) THEN
  locAddFlag = AddFlag
ELSE
  locAddFlag = .FALSE.
END IF

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nH2OSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nH2OSidesRecv(iProc)*(nSurfSample)**2
  IF (PRESENT(IntDataIN)) THEN
    CALL MPI_IRECV( SurfModelExchange%H2ORecvBuf(iProc)%content_int &
                  , MessageSize                                     &
                  , MPI_INT                                         &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1010                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%RecvRequest(1,iProc)          &
                  , IERROR )
  END IF
  IF (PRESENT(RealDataIN)) THEN
    CALL MPI_IRECV( SurfModelExchange%H2ORecvBuf(iProc)%content     &
                  , MessageSize                                     &
                  , MPI_DOUBLE_PRECISION                            &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1011                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%RecvRequest(2,iProc)          &
                  , IERROR )
  END IF
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nH2OSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  IF (PRESENT(IntDataIN))  SurfModelExchange%H2OSendBuf(iProc)%content_int = 0
  IF (PRESENT(RealDataIN)) SurfModelExchange%H2OSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfModelExchange%nH2OSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        iPos=iPos+1
        IF (PRESENT(RealDataIN)) SurfModelExchange%H2OSendBuf(iProc)%content(iPos) = RealDataIN(p,q,SurfSideID)
        IF (PRESENT(IntDataIN)) SurfModelExchange%H2OSendBuf(iProc)%content_int(iPos) = IntDataIN(p,q,SurfSideID)
        IF (locAddFlag) THEN
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = 0.
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = 0
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nH2OSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nH2OSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nH2OSidesSend(iProc)*(nSurfSample)**2
  IF (PRESENT(IntDataIN)) THEN
    CALL MPI_ISEND( SurfModelExchange%H2OSendBuf(iProc)%content_int &
                  , MessageSize                                     &
                  , MPI_INT                                         &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1010                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%SendRequest(1,iProc)          &
                  , IERROR)
  END IF
  IF (PRESENT(RealDataIN)) THEN
    CALL MPI_ISEND( SurfModelExchange%H2OSendBuf(iProc)%content     &
                  , MessageSize                                     &
                  , MPI_DOUBLE_PRECISION                            &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1011                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%SendRequest(2,iProc)          &
                  , IERROR)
  END IF
END DO ! iProc

! 4) Finish communication
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nH2OSidesSend(iProc).NE.0) THEN
    IF (PRESENT(IntDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%SendRequest(1,iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Send MPI Communication error in ExchangeSurfaceHaloToOrigin', IERROR)
    END IF
    IF (PRESENT(RealDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%SendRequest(2,iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Send MPI Communication error in ExchangeSurfaceHaloToOrigin', IERROR)
    END IF
  END IF
  IF(SurfModelExchange%nH2OSidesRecv(iProc).NE.0) THEN
    IF (PRESENT(IntDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Recv MPI Communication error in ExchangeSurfaceHaloToOrigin', IERROR)
    END IF
    IF (PRESENT(RealDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Recv MPI Communication error in ExchangeSurfaceHaloToOrigin', IERROR)
    END IF
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nH2OSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nH2OSidesSend(iProc)*(nSurfSample)**2
  iPos=0
  DO iSurfSide=1,SurfModelExchange%nH2OSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        iPos=iPos+1
        IF (locAddFlag) THEN
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = &
              RealDataIN(p,q,SurfSideID)+SurfModelExchange%H2ORecvBuf(iProc)%content(iPos)
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = &
              IntDataIN(p,q,SurfSideID)+SurfModelExchange%H2ORecvBuf(iProc)%content_int(iPos)
        ELSE
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = &
              SurfModelExchange%H2ORecvBuf(iProc)%content(iPos)
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = &
              SurfModelExchange%H2ORecvBuf(iProc)%content_int(iPos)
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nH2OSidesSend(iProc)
  IF (PRESENT(IntDataIN))  SurfModelExchange%H2ORecvBuf(iProc)%content_int = 0
  IF (PRESENT(RealDataIN)) SurfModelExchange%H2ORecvBuf(iProc)%content = 0.
END DO ! iProc

END SUBROUTINE ExchangeSurfaceHaloToOrigin


SUBROUTINE ExchangeSurfaceOriginToHalo(IntDataIN,RealDataIn,AddFlag)
!===================================================================================================================================
!> Each process sends his origin-information directly to the halo process by use of a list, containing the surfsideids for sending.
!> The receiving process adds the new data to his halo surfaces
!> Only processes with surfaces to send or receive participate in the communication
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm, nSurfSample, SurfMesh
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfModelExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL :: IntDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
REAL   ,INTENT(INOUT),OPTIONAL :: RealDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
LOGICAL,INTENT(IN)   ,OPTIONAL :: AddFlag
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: MessageSize,iSurfSide,SurfSideID
INTEGER :: iPos,p,q,iProc
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
LOGICAL :: locAddFlag
!===================================================================================================================================
IF (.NOT.ALLOCATED(SurfModelExchange%nO2HSidesSend) .OR. .NOT.ALLOCATED(SurfModelExchange%nO2HSidesRecv)) RETURN
IF (.NOT.PRESENT(IntDataIN) .AND. .NOT.PRESENT(RealDataIN) ) RETURN
IF (PRESENT(AddFlag)) THEN
  locAddFlag = AddFlag
ELSE
  locAddFlag = .FALSE.
END IF

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nO2HSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nO2HSidesRecv(iProc)*(nSurfSample)**2
  IF (PRESENT(IntDataIN)) THEN
    CALL MPI_IRECV( SurfModelExchange%O2HRecvBuf(iProc)%content_int &
                  , MessageSize                                     &
                  , MPI_INT                                         &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1013                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%RecvRequest(1,iProc)          &
                  , IERROR )
  END IF
  IF (PRESENT(RealDataIN)) THEN
    CALL MPI_IRECV( SurfModelExchange%O2HRecvBuf(iProc)%content     &
                  , MessageSize                                     &
                  , MPI_DOUBLE_PRECISION                            &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1014                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%RecvRequest(2,iProc)          &
                  , IERROR )
  END IF
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nO2HSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  IF (PRESENT(IntDataIN))  SurfModelExchange%O2HSendBuf(iProc)%content_int = 0
  IF (PRESENT(RealDataIN)) SurfModelExchange%O2HSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfModelExchange%nO2HSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        iPos=iPos+1
        IF (PRESENT(RealDataIN)) SurfModelExchange%O2HSendBuf(iProc)%content(iPos) = RealDataIN(p,q,SurfSideID)
        IF (PRESENT(IntDataIN)) SurfModelExchange%O2HSendBuf(iProc)%content_int(iPos) = IntDataIN(p,q,SurfSideID)
        IF (locAddFlag) THEN
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = 0.
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = 0
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nO2HSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nO2HSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nO2HSidesSend(iProc)*(nSurfSample)**2
  IF (PRESENT(IntDataIN)) THEN
    CALL MPI_ISEND( SurfModelExchange%O2HSendBuf(iProc)%content_int &
                  , MessageSize                                     &
                  , MPI_INT                                         &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1013                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%SendRequest(1,iProc)          &
                  , IERROR)
  END IF
  IF (PRESENT(RealDataIN)) THEN
    CALL MPI_ISEND( SurfModelExchange%O2HSendBuf(iProc)%content     &
                  , MessageSize                                     &
                  , MPI_DOUBLE_PRECISION                            &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID        &
                  , 1014                                            &
                  , SurfCOMM%COMM                                   &
                  , SurfModelExchange%SendRequest(2,iProc)          &
                  , IERROR)
  END IF
END DO ! iProc

! 4) Finish communication
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nO2HSidesSend(iProc).NE.0) THEN
    IF (PRESENT(IntDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%SendRequest(1,iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Send MPI Communication error in ExchangeSurfaceOriginToHalo', IERROR)
    END IF
    IF (PRESENT(RealDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%SendRequest(2,iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Send MPI Communication error in ExchangeSurfaceOriginToHalo', IERROR)
    END IF
  END IF
  IF(SurfModelExchange%nO2HSidesRecv(iProc).NE.0) THEN
    IF (PRESENT(IntDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Recv MPI Communication error in ExchangeSurfaceOriginToHalo', IERROR)
    END IF
    IF (PRESENT(RealDataIN)) THEN
      CALL MPI_WAIT(SurfModelExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
  __STAMP__&
            ,' Recv MPI Communication error in ExchangeSurfaceOriginToHalo', IERROR)
    END IF
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nO2HSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nO2HSidesSend(iProc)*(nSurfSample)**2
  iPos=0
  DO iSurfSide=1,SurfModelExchange%nO2HSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        iPos=iPos+1
        IF (locAddFlag) THEN
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = &
              RealDataIN(p,q,SurfSideID)+SurfModelExchange%O2HRecvBuf(iProc)%content(iPos)
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = &
              IntDataIN(p,q,SurfSideID)+SurfModelExchange%O2HRecvBuf(iProc)%content_int(iPos)
        ELSE
          IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = &
              SurfModelExchange%O2HRecvBuf(iProc)%content(iPos)
          IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = &
              SurfModelExchange%O2HRecvBuf(iProc)%content_int(iPos)
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nO2HSidesSend(iProc)
  IF (PRESENT(IntDataIN))  SurfModelExchange%O2HRecvBuf(iProc)%content_int = 0
  IF (PRESENT(RealDataIN)) SurfModelExchange%O2HRecvBuf(iProc)%content = 0.
END DO ! iProc

END SUBROUTINE ExchangeSurfaceOriginToHalo


SUBROUTINE InitSMCR_MPI()
!===================================================================================================================================
!> Initializing MPI for Kinetic Surface Model (SMCR model specific: communication of distribution data)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfModelExchange, SurfDistSendBuf, SurfDistRecvBuf
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iProc
!===================================================================================================================================
ALLOCATE(SurfModelExchange%nSurfDistSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%nSurfDistSidesRecv(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%SurfDistSendRequest(1:2,1:SurfCOMM%nMPINeighbors)  &
        ,SurfModelExchange%SurfDistRecvRequest(1:2,1:SurfCOMM%nMPINeighbors)  &
        ,SurfModelExchange%NbrOfPos(1:SurfCOMM%nMPINeighbors))
! allocate send and receive buffer
ALLOCATE(SurfDistSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfDistRecvBuf(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  SurfModelExchange%nSurfDistSidesSend(iProc) = SurfModelExchange%nO2HSidesRecv(iProc)
  SurfModelExchange%nSurfDistSidesRecv(iProc) = SurfModelExchange%nO2HSidesSend(iProc)
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(SurfModelExchange%nSurfDistSidesRecv(iProc)))
    SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
    ALLOCATE(SurfModelExchange%NbrOfPos(iProc)%nPosRecv(1:SurfModelExchange%nSurfDistSidesRecv(iProc)))
    SurfModelExchange%NbrOfPos(iProc)%nPosRecv = 0
  END IF
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(SurfModelExchange%nSurfDistSidesSend(iProc)))
    SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
    ALLOCATE(SurfModelExchange%NbrOfPos(iProc)%nPosSend(1:SurfModelExchange%nSurfDistSidesSend(iProc)))
    SurfModelExchange%NbrOfPos(iProc)%nPosSend = 0
  END IF
END DO ! iProc

! communicate the number of surface sites for surfdist communication
CALL ExchangeSurfDistSize()
! fill halo surface distribution through mpi communication
CALL ExchangeSurfDistInfo()

END SUBROUTINE InitSMCR_MPI


SUBROUTINE ExchangeSurfDistSize()
!===================================================================================================================================
!> Calculates and communicates the number of surface distribution sites that need to be send between the participating procs
!>   in subroutine ExchangeSurfDistInfo
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfComm
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfModelExchange
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfDistSendBuf, SurfDistRecvBuf, SurfModelExchange
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: MessageSize,iSurfSide,SurfSideID
INTEGER :: p,q,iProc
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
INTEGER :: iCoord
!===================================================================================================================================
IF (.NOT.ALLOCATED(SurfModelExchange%nSurfDistSidesSend) .OR. .NOT.ALLOCATED(SurfModelExchange%nSurfDistSidesRecv)) RETURN

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nSurfDistSidesRecv(iProc)
  CALL MPI_IRECV( SurfModelExchange%NbrOfPos(iProc)%nPosRecv(:)     &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1011                                         &
                , SurfCOMM%COMM                                &
                , SurfModelExchange%SurfDistRecvRequest(1,iProc)    &
                , IERROR )
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
  DO iSurfSide=1,SurfModelExchange%nSurfDistSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iCoord = 1,3
            SurfModelExchange%NbrOfPos(iProc)%nPosSend(iSurfSide) = SurfModelExchange%NbrOfPos(iProc)%nPosSend(iSurfSide) &
                                                             + 3 + 2*SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
        END DO
      END DO
    END DO
  END DO
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nSurfDistSidesSend(iProc)
  CALL MPI_ISEND( SurfModelExchange%NbrOfPos(iProc)%nPosSend(:)     &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1011                                         &
                , SurfCOMM%COMM                                &
                , SurfModelExchange%SurfDistSendRequest(1,iProc)    &
                , IERROR )
END DO ! iProc

! 4) Finish receiving commsizes
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SurfDistSendRequest(1,iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in number of surface distribution sites (send)', IERROR)
  END IF
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SurfDistRecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in number of surface distribution sites (receive)', IERROR)
  END IF
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    ALLOCATE(SurfDistSendBuf(iProc)%content_int(1:SUM(SurfModelExchange%NbrOfPos(iProc)%nPosSend(:))))
  END IF
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    ALLOCATE(SurfDistRecvBuf(iProc)%content_int(1:SUM(SurfModelExchange%NbrOfPos(iProc)%nPosRecv(:))))
  END IF
END DO

END SUBROUTINE ExchangeSurfDistSize


SUBROUTINE ExchangeSurfDistInfo()
!===================================================================================================================================
!> Exchanges information from origin to the respecitve halosides of corresponding neighbours
!>   The structure is similar to surface sampling/particle communication but has inverse communication path (origin->halo)
!> Exchanged information: Surface distribution
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, SurfComm, nSurfSample, PartBound
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfDistSendBuf, SurfDistRecvBuf, SurfModelExchange
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
USE MOD_Mesh_Vars              ,ONLY: BC,nBCSides,nSides
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: MessageSize, iSurfSide, SurfSideID, iSurf, SideID, PartboundID
INTEGER :: iPos,p,q,iProc
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
INTEGER :: iCoord,nSites,nSitesRemain,iSite,iInteratom,UsedSiteMapPos,iSpec,xpos,ypos
INTEGER :: iSide,TargetHaloSide,SurfSideHaloID
!===================================================================================================================================
IF (.NOT.ALLOCATED(SurfModelExchange%nSurfDistSidesSend) .OR. .NOT.ALLOCATED(SurfModelExchange%nSurfDistSidesRecv)) RETURN

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  Messagesize = 0
  DO iSurf = 1,SurfModelExchange%nSurfDistSidesRecv(iProc)
    MessageSize = MessageSize + SurfModelExchange%NbrOfPos(iProc)%nPosRecv(iSurf)
  END DO
  CALL MPI_IRECV( SurfDistRecvBuf(iProc)%content_int           &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1012                                         &
                , SurfCOMM%COMM                                &
                , SurfModelExchange%SurfDistRecvRequest(2,iProc)    &
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  SurfDistSendBuf(iProc)%content_int = 0
  DO iSurfSide=1,SurfModelExchange%nSurfDistSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iCoord = 1,3
          SurfDistSendBuf(iProc)%content_int(iPos+1) = SurfDistInfo(p,q,SurfSideID)%SitesRemain(iCoord)
          iPos=iPos+1
          SurfDistSendBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)) = &
              SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%Species(:)
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
          SurfDistSendBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)) = &
              SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%UsedSiteMap(:)
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
        END DO
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nSurfDistSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  Messagesize = 0
  DO iSurf = 1,SurfModelExchange%nSurfDistSidesSend(iProc)
    MessageSize = MessageSize + SurfModelExchange%NbrOfPos(iProc)%nPosSend(iSurf)
  END DO
  CALL MPI_ISEND( SurfDistSendBuf(iProc)%content_int        &
                , MessageSize                               &
                , MPI_INT                                   &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID  &
                , 1012                                      &
                , SurfCOMM%COMM                             &
                , SurfModelExchange%SurfDistSendRequest(2,iProc) &
                , IERROR )
END DO ! iProc

! 4) Finish received surface distribution
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SurfDistSendRequest(2,iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in surface distribution (send)', IERROR)
  END IF
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SurfDistRecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in Surface distribution (receive)', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfModelExchange%nSurfDistSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iCoord = 1,3
          SurfDistInfo(p,q,SurfSideID)%SitesRemain(iCoord) = SurfDistRecvBuf(iProc)%content_int(iPos+1)
          iPos=iPos+1
          SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%Species(:) = &
              SurfDistRecvBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord))
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
          SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%UsedSiteMap(:) = &
              SurfDistRecvBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord))
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
        END DO
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample.
  END DO ! iSurfSide=1,nSurfModelExchange%nH2OSidesSend(iProc)
  SurfDistRecvBuf(iProc)%content_int = 0
END DO ! iProc

IF(SurfMesh%nSides.GT.SurfMesh%nMasterSides) THEN ! There are reflective inner BCs on SlaveSide
  DO iSide=nBCSides+1,nSides
    IF(BC(iSide).EQ.0) CYCLE
    IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN
      IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.-1) THEN ! SlaveSide
        DO q=1,nSurfSample
          DO p=1,nSurfSample
            TargetHaloSide = SurfMesh%innerBCSideToHaloMap(iSide)
            SurfSideID=SurfMesh%SideIDToSurfID(iSide)
            SurfSideHaloID=SurfMesh%SideIDToSurfID(TargetHaloSide)
            ! map distribution data (haloinnersurface->innersurface)
            IF (ALLOCATED(SurfDistInfo)) THEN
              SurfDistInfo(p,q,SurfSideID)%SitesRemain(:) = SurfDistInfo(p,q,SurfSideHaloID)%SitesRemain(:)
              DO iCoord = 1,3
                SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%Species(:) = &
                    SurfDistInfo(p,q,SurfSideHaloID)%AdsMap(iCoord)%Species(:)
                SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%UsedSiteMap(:) = &
                    SurfDistInfo(p,q,SurfSideHaloID)%AdsMap(iCoord)%UsedSiteMap(:)
              END DO
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END IF

! assign bond order to surface atoms in the surfacelattice for halo sides
DO iSurfSide = SurfMesh%nMasterSides+1,SurfMesh%nTotalSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%SurfaceModel(PartboundID).NE.3) CYCLE
  DO q=1,nSurfSample
    DO p=1,nSurfSample
      SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(:,:,:) = 0
      DO iCoord = 1,3
        nSitesRemain = SurfDistInfo(p,q,iSurfSide)%SitesRemain(iCoord)
        nSites = SurfDistInfo(p,q,iSurfSide)%nSites(iCoord)
        DO iSite = nSitesRemain+1,nSites
          UsedSiteMapPos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(iSite)
          iSpec = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos)
          DO iInterAtom = 1,SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%nInterAtom
            xpos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
            ypos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
            SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
              SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
          END DO
        END DO
      END DO
    END DO
  END DO
END DO

END SUBROUTINE ExchangeSurfDistInfo


SUBROUTINE MapHaloInnerToOriginInnerSurf(IntDataIN,RealDataIN,AddFlag,Reverse)
!===================================================================================================================================
! Map the surface data from innerBC SlaveSides to corresponding HaloSide.
! All surfacemodel informations of a innerBC SlaveSide is added to corresponding HaloSide and vice versa
! Afterwards, these informations are send to innerBC MasterSide in a second ExchangeData call
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh,nSurfSample,PartBound
USE MOD_Mesh_Vars              ,ONLY: nBCSides,nSides,BC
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL :: IntDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
REAL   ,INTENT(INOUT),OPTIONAL :: RealDataIN(nSurfSample,nSurfSample,SurfMesh%nTotalSides)
LOGICAL,INTENT(IN)   ,OPTIONAL :: AddFlag
LOGICAL,INTENT(IN)   ,OPTIONAL :: Reverse
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iSide,TargetHaloSide,q,p,SurfSideID,SurfSideHaloID
LOGICAL :: locAddFlag, locReverse
!===================================================================================================================================
IF (.NOT.PRESENT(IntDataIN) .AND. .NOT.PRESENT(RealDataIN) ) RETURN
IF (PRESENT(AddFlag)) THEN
  locAddFlag = AddFlag
ELSE
  locAddFlag = .FALSE.
END IF
IF (PRESENT(Reverse)) THEN
  locReverse = Reverse
ELSE
  locReverse = .FALSE.
END IF

IF(SurfMesh%nSides.GT.SurfMesh%nMasterSides) THEN ! There are reflective inner BCs on SlaveSide
  DO iSide=nBCSides+1,nSides
    IF(BC(iSide).EQ.0) CYCLE
    IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN
      IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.-1) THEN ! SlaveSide
        DO q=1,nSurfSample
          DO p=1,nSurfSample
            TargetHaloSide = SurfMesh%innerBCSideToHaloMap(iSide)
            SurfSideID=SurfMesh%SideIDToSurfID(iSide)
            SurfSideHaloID=SurfMesh%SideIDToSurfID(TargetHaloSide)
            IF (locReverse) THEN
              IF (locAddFlag) THEN
                IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = RealDataIN(p,q,SurfSideID)+RealDataIN(p,q,SurfSideHaloID)
                IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = IntDataIN(p,q,SurfSideID)+IntDataIN(p,q,SurfSideHaloID)
              ELSE
                IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = RealDataIN(p,q,SurfSideHaloID)
                IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = IntDataIN(p,q,SurfSideHaloID)
              END IF
            ELSE
              IF (locAddFlag) THEN
                IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = RealDataIN(p,q,SurfSideID)+RealDataIN(p,q,SurfSideHaloID)
                IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = IntDataIN(p,q,SurfSideID)+IntDataIN(p,q,SurfSideHaloID)
              ELSE
                IF (PRESENT(RealDataIN)) RealDataIN(p,q,SurfSideID) = RealDataIN(p,q,SurfSideHaloID)
                IF (PRESENT(IntDataIN)) IntDataIN(p,q,SurfSideID) = IntDataIN(p,q,SurfSideHaloID)
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END IF
END SUBROUTINE MapHaloInnerToOriginInnerSurf
#endif /*USE_MPI*/


END MODULE MOD_SurfaceModel_MPI
