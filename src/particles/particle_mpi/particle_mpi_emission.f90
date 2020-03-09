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

MODULE MOD_Particle_MPI_Emission
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE SendEmissionParticlesToProcs
  MODULE PROCEDURE SendEmissionParticlesToProcs
END INTERFACE

!===================================================================================================================================
PUBLIC :: SendEmissionParticlesToProcs
!===================================================================================================================================
CONTAINS



SUBROUTINE SendEmissionParticlesToProcs(chunkSize,DimSend,particle_positions)
!----------------------------------------------------------------------------------------------------------------------------------!
! A particle's host cell in the FIBGM is found and the corresponding procs are notified.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_Particle_Mesh_Tools ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars  ,ONLY: FIBGM_nElems, FIBGM_offsetElem, FIBGM_Element
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: chunkSize
INTEGER,INTENT(IN) :: DimSend
REAL,INTENT(IN)    :: particle_positions(1:chunkSize*DimSend)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iProc,tProc, ijkBGM(3,chunkSize), iDir, ElemID, ProcID
INTEGER          :: iBGMElem,nBGMElems,TotalNbrOfRecvParts,TotalNbrOfSendParts
LOGICAL          :: InsideMyBGM(chunkSize)
REAL,ALLOCATABLE :: recvPartPos(:)
!===================================================================================================================================
ALLOCATE( PartMPIInsert%nPartsSend  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
ALLOCATE( PartMPIInsert%nPartsRecv  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
ALLOCATE( PartMPIInsert%SendRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
ALLOCATE( PartMPIInsert%RecvRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
ALLOCATE( PartMPIInsert%send_message(0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
PartMPIInsert%nPartsSend(:)=0

!--- 1/4 Open receive buffer
DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
  !--- MPI_IRECV lengths of lists of particles entering local mesh
  CALL MPI_IRECV(PartMPIInsert%nPartsRecv(iProc), 1, MPI_INTEGER, iProc, 1011, PartMPI%InitGroup(InitGroup)%COMM, &
      PartMPIInsert%RecvRequest(iProc,1), IERROR)
END DO

! Identify particles that are on the node (or in the halo region of the node) or on other nodes
DO i=1,chunkSize

  ! Set BGM cell index
  ASSOCIATE( xMin   => (/GEO%xminglob  , GEO%yminglob  , GEO%zminglob/)  , &
             BGMMin => (/GEO%FIBGMimin , GEO%FIBGMjmin , GEO%FIBGMkmin/) , &
             BGMMax => (/GEO%FIBGMimax , GEO%FIBGMjmax , GEO%FIBGMkmax/) )
    DO iDir = 1, 3
      ijkBGM(iDir,i) = INT((particle_positions(DimSend*(i-1)+iDir)-xMin(iDir))/GEO%FIBGMdeltas(iDir))+1
    END DO ! iDir = 1, 3
  END ASSOCIATE

  ! Check BGM cell index
  InsideMyBGM(i)=.TRUE.
  DO iDir = 1, 3
    IF(ijkBGM(iDir,i).GT.BGMMin(iDir)) THEN
      InsideMyBGM(i)=.FALSE.
      EXIT
    END IF
    IF(ijkBGM(iDir,i).LT.BGMMax(iDir)) THEN
      InsideMyBGM(i)=.FALSE.
      EXIT
    END IF
  END DO ! iDir = 1, 3

  IF(InsideMyBGM(i))THEN
    !--- check all cells associated with this background mesh cell
    nBGMElems = FIBGM_nElems(ijkBGM(1,i),ijkBGM(2,i),ijkBGM(3,i))

    DO iBGMElem = 1, nBGMElems
      ElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem))

      ! Check if element is on node (or halo region of node)
      IF(ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.1)THEN ! it is 0 or 2
        InsideMyBGM(i) = .FALSE.
      END IF ! ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.1
    END DO ! iBGMElem = 1, nBGMElems
  END IF ! InsideMyBGM(i)
END DO ! i = 1, chunkSize

!--- Find particles for sending to other nodes
DO i = 1, chunkSize
  IF(.NOT.InsideMyBGM(i))THEN
    !--- check all cells associated with this beckground mesh cell
    nBGMElems = FIBGM_nElems(ijkBGM(1,i),ijkBGM(2,i),ijkBGM(3,i))

    ! Loop over all BGM elements and count number of particles per procs for sending
    DO iBGMElem = 1, nBGMElems
      ElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem))
      ProcID = ElemToProcID_Shared(ElemID)

      tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(ProcID)
      IF(tProc.EQ.-1)CYCLE ! Processor is not on emission communicator
      PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
    END DO
  END IF ! .NOT.InsideMyBGM(i)
END DO ! i = 1, chunkSize


!--- 2/4 Send number of particles
DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
  ! send particles
  !--- MPI_ISEND lengths of lists of particles leaving local mesh
  CALL MPI_ISEND(PartMPIInsert%nPartsSend(iProc), 1, MPI_INTEGER, iProc, 1011, PartMPI%InitGroup(InitGroup)%COMM, &
      PartMPIInsert%SendRequest(iProc,1), IERROR)
  IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
    ALLOCATE( PartMPIInsert%send_message(iProc)%content(1:DimSend*PartMPIInsert%nPartsSend(iProc)), STAT=allocStat )
  END IF
END DO


!--- 3/4 Send actual particles
PartMPIInsert%nPartsSend(:)=0
DO i = 1, chunkSize
  IF(.NOT.InsideMyBGM(i))THEN
    !--- check all cells associated with this beckground mesh cell
    nBGMElems = FIBGM_nElems(ijkBGM(1,i),ijkBGM(2,i),ijkBGM(3,i))

    ! Loop over all BGM elements and count number of particles per procs for sending
    DO iBGMElem = 1, nBGMElems
      ElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem))
      ProcID = ElemToProcID_Shared(ElemID)

      tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(ProcID)
      IF(tProc.EQ.-1)CYCLE ! Processor is not on emission communicator
      PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
      TotalNbrOfSendParts=PartMPIInsert%nPartsSend(tProc)
      PartMPIInsert%send_message(tProc)%content(DimSend*(TotalNbrOfSendParts-1)+1:DimSend*TotalNbrOfSendParts)=&
          particle_positions(DimSend*(i-1)+1:DimSend*i)
    END DO
  END IF ! .NOT.InsideMyBGM(i)
END DO ! i = 1, chunkSize


!--- 4/4 Receive actual particles 
DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
  CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,1),msg_status(:),IERROR)
END DO
TotalNbrOfRecvParts=SUM(PartMPIInsert%nPartsRecv)
ALLOCATE(recvPartPos(1:TotalNbrOfRecvParts*DimSend), STAT=allocStat)
TotalNbrOfRecvParts=0
DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
  IF (PartMPIInsert%nPartsRecv(iProc).GT.0) THEN
  !--- MPI_IRECV lengths of lists of particles entering local mesh
    CALL MPI_IRECV(recvPartPos(TotalNbrOfRecvParts*DimSend+1), DimSend*PartMPIInsert%nPartsRecv(iProc),&
                                              MPI_DOUBLE_PRECISION, iProc, 1022,   &
                                              PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%RecvRequest(iProc,2), IERROR)
    TotalNbrOfRecvParts=TotalNbrOfRecvParts+PartMPIInsert%nPartsRecv(iProc)
  END IF
  !--- (non-blocking:) send messages to all procs receiving particles from myself
  IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
    CALL MPI_ISEND(PartMPIInsert%send_message(iProc)%content, DimSend*PartMPIInsert%nPartsSend(iProc),&
     MPI_DOUBLE_PRECISION, iProc, 1022, PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%SendRequest(iProc,2), IERROR)
  END IF
END DO


!--- Locate local (node or halo of node) particles
DO i = 1, chunkSize
  IF(InsideMyBGM(i))THEN
    ! TODO: continue programming philipesque here ...
    CALL LocateParticleInsideELement()
  END IF
END DO ! i = 1, chunkSize


!--- 5/4 Receive actual particles 
DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
  IF (PartMPIInsert%nPartsRecv(iProc).GT.0) THEN
    CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,2),msg_status(:),IERROR)
  END IF
END DO





END SUBROUTINE SendEmissionParticlesToProcs

#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Emission
