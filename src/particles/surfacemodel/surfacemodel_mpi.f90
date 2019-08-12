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
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI
PUBLIC :: InitSurfModel_MPI
PUBLIC :: InitSMCR_MPI
PUBLIC :: ExchangeAdsorbNum
PUBLIC :: ExchangeSurfDistInfo
PUBLIC :: ExchangeCoverageInfo
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
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars      ,ONLY: SurfExchange
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: AdsorbSendBuf, AdsorbRecvBuf, SurfModelExchange
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfCoverageSendBuf, SurfCoverageRecvBuf
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
ALLOCATE(SurfModelExchange%nSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%nSidesRecv(1:SurfCOMM%nMPINeighbors) &
        ,SurfModelExchange%SendRequest(SurfCOMM%nMPINeighbors)  &
        ,SurfModelExchange%RecvRequest(SurfCOMM%nMPINeighbors)  )
SurfModelExchange%nSidesSend = SurfExchange%nSidesSend
SurfModelExchange%nSidesRecv = SurfExchange%nSidesRecv

! allocate send and receive buffer for sum arrays of adsorbing particles from halo sides to origin sides
ALLOCATE(AdsorbSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(AdsorbRecvBuf(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  ALLOCATE(AdsorbSendBuf(iProc)%content_int(2*nSpecies*(nSurfSample**2)*SurfModelExchange%nSidesSend(iProc)))
  ALLOCATE(AdsorbRecvBuf(iProc)%content_int(2*nSpecies*(nSurfSample**2)*SurfModelExchange%nSidesRecv(iProc)))
  AdsorbSendBuf(iProc)%content_int=0
  AdsorbRecvBuf(iProc)%content_int=0
END DO ! iProc

! currently only needed for SurfaceModel=2
! represents inverse communication where information is send from origin sides to halo sides (local sides of neighbor halo procs)
ALLOCATE(SurfCoverageSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfCoverageRecvBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfModelExchange%nCoverageSidesSend(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfModelExchange%nCoverageSidesRecv(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  SurfModelExchange%nCoverageSidesSend(iProc) = SurfModelExchange%nSidesRecv(iProc)
  SurfModelExchange%nCoverageSidesRecv(iProc) = SurfModelExchange%nSidesSend(iProc)
  IF(SurfModelExchange%nCoverageSidesRecv(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageRecvList(SurfModelExchange%nCoverageSidesRecv(iProc)))
    SurfCOMM%MPINeighbor(iProc)%CoverageRecvList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
    ALLOCATE(SurfCoverageRecvBuf(iProc)%content(1:3*nSpecies*(nSurfSample**2)*SurfModelExchange%nCoverageSidesRecv(iProc)))
    SurfCoverageRecvBuf(iProc)%content = 0.
  END IF
  IF(SurfModelExchange%nCoverageSidesSend(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageSendList(SurfModelExchange%nCoverageSidesSend(iProc)))
    SurfCOMM%MPINeighbor(iProc)%CoverageSendList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
    ALLOCATE(SurfCoverageSendBuf(iProc)%content(1:3*nSpecies*(nSurfSample**2)*SurfModelExchange%nCoverageSidesSend(iProc)))
    SurfCoverageSendBuf(iProc)%content = 0.
  END IF
END DO

END SUBROUTINE InitSurfModel_MPI


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
  SurfModelExchange%nSurfDistSidesSend(iProc) = SurfModelExchange%nSidesRecv(iProc)
  SurfModelExchange%nSurfDistSidesRecv(iProc) = SurfModelExchange%nSidesSend(iProc)
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


SUBROUTINE ExchangeAdsorbNum()
!===================================================================================================================================
!> Exchange the number of particles that adsorbed on a halo surface to origin surface of the corresponding proc
!>   only processes with samling side to send or recieve participate in the communication
!>   structure is similar to surface sampling/particle communication (halo->origin)
!> Each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
!>   the receiving process adds the new data to his own sides
!> Exchanged arrays: SumAdsorbPart, SumERDesorbed
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm, nSurfSample
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: AdsorbSendBuf, AdsorbRecvBuf, SurfModelExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: MessageSize,nValues,iSurfSide,SurfSideID
INTEGER :: iPos,p,q,iProc
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nValues=2*nSpecies*(nSurfSample)**2

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( AdsorbRecvBuf(iProc)%content_int             &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1010                                         &
                , SurfCOMM%COMM                                &
                , SurfModelExchange%RecvRequest(iProc)              &
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  AdsorbSendBuf(iProc)%content_int = 0
  DO iSurfSide=1,SurfModelExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        AdsorbSendBuf(iProc)%content_int(iPos+1:iPos+nSpecies)= SurfModel%SumAdsorbPart(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
        AdsorbSendBuf(iProc)%content_int(iPos+1:iPos+nSpecies)= SurfModel%SumERDesorbed(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( AdsorbSendBuf(iProc)%content_int         &
                , MessageSize                              &
                , MPI_INT                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
                , 1010                                     &
                , SurfCOMM%COMM                            &
                , SurfModelExchange%SendRequest(iProc)          &
                , IERROR)
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfModelExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfModelExchange%nSidesSend(iProc)*nValues
  iPos=0
  DO iSurfSide=1,SurfModelExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SurfModel%SumAdsorbPart(p,q,SurfSideID,:)=SurfModel%SumAdsorbPart(p,q,SurfSideID,:) &
                                         +AdsorbRecvBuf(iProc)%content_int(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
        SurfModel%SumERDesorbed(p,q,SurfSideID,:)=SurfModel%SumERDesorbed(p,q,SurfSideID,:) &
                                         +AdsorbRecvBuf(iProc)%content_int(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfModelExchange%nSidesSend(iProc)
  AdsorbRecvBuf(iProc)%content_int = 0
END DO ! iProc

END SUBROUTINE ExchangeAdsorbNum


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
USE MOD_Mesh_Vars              ,ONLY: BC
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
!===================================================================================================================================

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
  END DO ! iSurfSide=1,nSurfModelExchange%nSidesSend(iProc)
  SurfDistRecvBuf(iProc)%content_int = 0
END DO ! iProc

! assign bond order to surface atoms in the surfacelattice for halo sides
DO iSurfSide = SurfMesh%nSides+1,SurfMesh%nTotalSides
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


SUBROUTINE ExchangeCoverageInfo()
!===================================================================================================================================
!> Exchanges information from origin to the halosides of corresponding neighbours
!>   The structure is similar to surface sampling/particle communication but has inverse communication path (origin->halo)
!> Exchanged information: Surface coverage, constant probabilities
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfComm, nSurfSample
USE MOD_SurfaceModel_MPI_Vars  ,ONLY: SurfCoverageSendBuf, SurfCoverageRecvBuf, SurfModelExchange
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: MessageSize, iSurfSide, SurfSideID
INTEGER :: iPos,p,q,iProc, nValues
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nValues = 3*nSpecies * nSurfSample**2
! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nCoverageSidesRecv(iProc).EQ.0) CYCLE
  Messagesize = SurfModelExchange%nCoverageSidesRecv(iProc)*nValues
  CALL MPI_IRECV( SurfCoverageRecvBuf(iProc)%content           &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1013                                         &
                , SurfCOMM%COMM                                &
                , SurfModelExchange%RecvRequest(iProc)&
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nCoverageSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  SurfCoverageSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfModelExchange%nCoverageSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%CoverageSendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SurfCoverageSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = Adsorption%Coverage(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
        SurfCoverageSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = Adsorption%ProbAds(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
        SurfCoverageSendBuf(iProc)%content(iPos+1:iPos+nSpecies) = Adsorption%ProbDes(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,SurfModelExchange%nSurfDistSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nCoverageSidesSend(iProc).EQ.0) CYCLE
  Messagesize = SurfModelExchange%nCoverageSidesSend(iProc)*nValues
  CALL MPI_ISEND( SurfCoverageSendBuf(iProc)%content       &
                , MessageSize                              &
                , MPI_INT                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
                , 1013                                     &
                , SurfCOMM%COMM                            &
                , SurfModelExchange%SendRequest(iProc)&
                , IERROR )
END DO ! iProc

! 4) Finish received surface distribution
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nCoverageSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in surface distribution (send)', IERROR)
  END IF
  IF(SurfModelExchange%nCoverageSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfModelExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in Surface distribution (receive)', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfModelExchange%nCoverageSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfModelExchange%nCoverageSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%CoverageRecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        Adsorption%Coverage(p,q,SurfSideID,:) = SurfCoverageRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
        Adsorption%ProbAds(p,q,SurfSideID,:) = SurfCoverageRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
        Adsorption%ProbDes(p,q,SurfSideID,:) = SurfCoverageRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample.
  END DO ! iSurfSide=1,SurfModelExchange%nSidesSend(iProc)
  SurfCoverageRecvBuf(iProc)%content = 0.
END DO ! iProc

END SUBROUTINE ExchangeCoverageInfo
#endif /*USE_MPI*/

END MODULE MOD_SurfaceModel_MPI
