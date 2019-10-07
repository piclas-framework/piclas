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

MODULE MOD_PICDepo_MPI
!===================================================================================================================================
! MPI related routines required for deposition
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!===================================================================================================================================

#if USE_MPI
PUBLIC:: MPISourceExchangeBGM,MPIBackgroundMeshInit
#else /*NOT USE_MPI*/
PUBLIC:: PeriodicSourceExchange
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI
SUBROUTINE MPISourceExchangeBGM()
!=================================================================================================================================
! Exchange sources in periodic case for MPI
!==================================================================================================================================
! use MODULES
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI,tMPIMEssage
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)        :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)
TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)
INTEGER                     :: send_request(0:PartMPI%nProcs-1)
INTEGER                     :: recv_request(0:PartMPI%nProcs-1)
INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)
INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)
INTEGER                     :: iProc,k,l,m,n, ppp, Counter, MsgLength(0:PartMPI%nProcs-1), iPer
INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)
INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)
INTEGER                     :: allocStat, Counter2
INTEGER                     :: messageCounterS, messageCounterR
INTEGER                     :: myRealKind, k2,l2,m2
REAL                        :: myRealTestValue
!-----------------------------------------------------------------------------------------------------------------------------------

myRealKind = KIND(myRealTestValue)
IF (myRealKind.EQ.4) THEN
 myRealKind = MPI_REAL
ELSE IF (myRealKind.EQ.8) THEN
 myRealKind = MPI_DOUBLE_PRECISION
ELSE
 myRealKind = MPI_REAL
END IF

!--- Determine which Sources actually need to be sent (<> 0) and build corresponding true/false list
!    One list per process
DO iProc = 0,PartMPI%nProcs-1
   MsgLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      MsgLength(iProc) = (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3) + 1)
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO k = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         MsgLength(iProc) = MsgLength(iProc) + &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,1) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,2) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,3) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
      END DO
   END IF
   IF((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor))THEN
      ALLOCATE(send_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
      ALLOCATE(recv_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
      END IF
   END IF
   !--- check which sources are <> 0
   Counter = 0
   SourceLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
               send_message(iProc)%content_log(Counter) = .TRUE.
               SourceLength(iProc) = SourceLength(iProc) + 1
            ELSE
               send_message(iProc)%content_log(Counter) = .FALSE.
            END IF
          END DO
        END DO
      END DO
   END IF
   !--- same for periodic
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
          DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                 PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
            DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                   PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter = Counter + 1
               IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
                  send_message(iProc)%content_log(Counter) = .TRUE.
                  SourceLength(iProc) = SourceLength(iProc) + 1
               ELSE
                  send_message(iProc)%content_log(Counter) = .FALSE.
               END IF
            END DO
          END DO
         END DO
      END DO
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%Content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF

!--- Assemble actual sources to send
 DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ALLOCATE(send_message(iProc)%content(1:SourceLength(iProc)*4), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
   END IF
   Counter = 0
   Counter2 = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
             Counter = Counter + 1
             IF (send_message(iProc)%content_log(Counter)) THEN
                Counter2 = Counter2 + 1
                DO n = 1,4
                   send_message(iProc)%content((Counter2-1)*4 +n) = BGMSource(k,l,m,n)
                END DO
             END IF
          END DO
        END DO
      END DO
   END IF
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
           DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                Counter = Counter + 1
                IF (send_message(iProc)%content_log(Counter)) THEN
                   Counter2 = Counter2 + 1
                   DO ppp = 1,4
                      send_message(iProc)%content((Counter2-1)*4 +ppp) = BGMSource(k,l,m,ppp)
                   END DO
                END IF
             END DO
           END DO
         END DO
      END DO
   END IF
END DO

!--- allocate actual PartSource receive buffer
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      Counter = 0
      DO k = 1, MsgLength(iProc)
         IF (recv_message(iProc)%Content_log(k)) THEN
            Counter = Counter + 1
         END IF
      END DO
      RecvLength(iProc) = Counter
      IF (RecvLength(iProc).GT.0) THEN
         ALLOCATE(recv_message(iProc)%content(1:Counter*4), STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
         END IF
      END IF
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content,SourceLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(RecvLength(iProc).GT.0)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%content,RecvLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF
!--- Deallocate Send Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      DEALLOCATE(send_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
      END IF
      IF (SourceLength(iProc).GT.0) THEN
         DEALLOCATE(send_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
         END IF
      END IF
   END IF
END DO

!--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
!---           received sources have been added!
IF ((GEO%nPeriodicVectors.GT.0).AND.(GEO%SelfPeriodic)) THEN
   DO iPer = 1, GEO%nPeriodicVectors
      DO k = BGMminX, BGMmaxX
         k2 = k + GEO%PeriodicBGMVectors(1,iPer)
        DO l = BGMminY, BGMmaxY
          l2 = l + GEO%PeriodicBGMVectors(2,iPer)
          DO m = BGMminZ, BGMmaxZ
             m2 = m + GEO%PeriodicBGMVectors(3,iPer)
             IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
             IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
             IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
                BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
                BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
             END IF
             END IF
             END IF
          END DO
        END DO
      END DO
   END DO
END IF

!--- Add Sources and Deallocate Receive Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF (RecvLength(iProc).GT.0) THEN
      Counter = 0
      Counter2 = 0
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
         DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
         DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
         DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF(recv_message(iProc)%content_log(Counter))THEN
               Counter2 = Counter2 + 1
               DO n = 1,4
                 BGMSource(k,l,m,n) = BGMSource(k,l,m,n) + recv_message(iProc)%content((Counter2-1)*4+n)
               END DO
            END IF
         END DO
         END DO
         END DO
      END IF
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
         DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
           DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
               DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                      PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter = Counter + 1
                  IF(recv_message(iProc)%content_log(Counter))THEN
                     Counter2 = Counter2 + 1
                     DO ppp = 1,4
                        BGMSource(k,l,m,ppp) = BGMSource(k,l,m,ppp) + recv_message(iProc)%content((Counter2-1)*4+ppp)
                     END DO
                  END IF
               END DO
             END DO
           END DO
         END DO
      END IF
      IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
         DEALLOCATE(recv_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
         END IF
      END IF
   END IF
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
      DEALLOCATE(recv_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
      END IF
   END IF
END DO

END SUBROUTINE MPISourceExchangeBGM


SUBROUTINE MPIBackgroundMeshInit()
!==================================================================================================================================
! initialize MPI background mesh
!==================================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Globals
! USE MOD_Particle_Vars,       ONLY:
USE MOD_Particle_Mesh_Vars,  ONLY:GEO
USE MOD_Particle_MPI_Vars,   ONLY:PartMPI
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iProc
INTEGER                     :: k,m,n,m0,n0
INTEGER                     :: localminmax(6), maxofmin, minofmax
INTEGER                     :: completeminmax(6*PartMPI%nProcs)
INTEGER                     :: allocStat, NeighCount
INTEGER                     :: TempBorder(1:2,1:3)
INTEGER                     :: Periodicminmax(6), coord, PeriodicVec(1:3)
INTEGER                     :: TempPeriBord(1:26,1:2,1:3)
LOGICAL                     :: CHECKNEIGHBOR
!-----------------------------------------------------------------------------------------------------------------------------------

!Periodic Init stuff
IF(GEO%nPeriodicVectors.GT.0)THEN
  ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
  ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
  IF (allocStat .NE. 0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
  END IF
  DO iProc = 1, GEO%nPeriodicVectors
    GEO%PeriodicBGMVectors(1,iProc) = NINT(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1))
    IF(ABS(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(2,iProc) = NINT(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2))
    IF(ABS(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(3,iProc) = NINT(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3))
    IF(ABS(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
  END DO
  ! Check whether process is periodic with itself
  GEO%SelfPeriodic = .FALSE.
  !--- virtually move myself according to periodic vectors in order to find overlapping areas
  !--- 26 possibilities,
  localminmax(1) = BGMminX
  localminmax(2) = BGMminY
  localminmax(3) = BGMminZ
  localminmax(4) = BGMmaxX
  localminmax(5) = BGMmaxY
  localminmax(6) = BGMmaxZ
  DO k = -1,1
    DO m = -1,1
      DO n = -1,1
        PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
        IF (GEO%nPeriodicVectors.GT.1) THEN
          PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
        END IF
        IF (GEO%nPeriodicVectors.GT.2) THEN
          PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
        END IF
        IF (ALL(PeriodicVec(:).EQ.0)) CYCLE
        periodicminmax(1) = localminmax(1) + PeriodicVec(1)
        periodicminmax(2) = localminmax(2) + PeriodicVec(2)
        periodicminmax(3) = localminmax(3) + PeriodicVec(3)
        periodicminmax(4) = localminmax(4) + PeriodicVec(1)
        periodicminmax(5) = localminmax(5) + PeriodicVec(2)
        periodicminmax(6) = localminmax(6) + PeriodicVec(3)
        !--- find overlap
        DO coord = 1,3  !          x y z direction
          maxofmin = MAX(periodicminmax(coord),localminmax(coord))
          minofmax = MIN(periodicminmax(3+coord),localminmax(3+coord))
          IF (maxofmin.LE.minofmax) GEO%SelfPeriodic = .TRUE.     !  overlapping
        END DO
      END DO
    END DO
  END DO
END IF

! --- send and receive min max indices to and from all processes

!--- enter local min max vector (xmin, ymin, zmin, xmax, ymax, zmax)
localminmax(1) = BGMminX
localminmax(2) = BGMminY
localminmax(3) = BGMminZ
localminmax(4) = BGMmaxX
localminmax(5) = BGMmaxY
localminmax(6) = BGMmaxZ
!--- do allgather into complete min max vector
CALL MPI_ALLGATHER(localminmax,6,MPI_INTEGER,completeminmax,6,MPI_INTEGER,PartMPI%COMM,IERROR)
! Allocate MPIConnect
SDEALLOCATE(PartMPI%DepoBGMConnect)
ALLOCATE(PartMPI%DepoBGMConnect(0:PartMPI%nProcs-1),STAT=allocStat)
IF (allocStat .NE. 0) THEN
  CALL abort(&
  __STAMP__&
  ,' Cannot allocate PartMPI%DepoBGMConnect')
END IF

!--- determine borders indices (=overlapping BGM mesh points) with each process
DO iProc=0,PartMPI%nProcs-1
  PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
   IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
   ELSE
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .TRUE.
      DO k = 1,3          !  x y z direction
         maxofmin = MAX(localminmax(k),completeminmax((iProc*6)+k))
         minofmax = MIN(localminmax(3+k),completeminmax((iProc*6)+3+k))
         IF (maxofmin.LE.minofmax) THEN          !  overlapping
            TempBorder(1,k) = maxofmin
            TempBorder(2,k) = minofmax
         ELSE
            PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
         END IF
      END DO
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)THEN
      SDEALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder)
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,' Cannot allocate PartMPI%DepoMPIConnect%BGMBorder')
      END IF
      PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3) = TempBorder(1:2,1:3)
   END IF
END DO ! iProc=0,PartMPI%nProcs-1

!--- determine border indices for periodic meshes
IF (GEO%nPeriodicVectors.GT.0) THEN
  !--- m-/-n-loops must be executed just once (with 0) if respective PV is not present
  !    (otherwise the unnec. 0 are sent and the 0-0-0-case occurs two 2 more which are not cycled, see below for more info)!
  IF (GEO%nPeriodicVectors.GT.1) THEN
    m0=1
  ELSE
    m0=0
  END IF
  IF (GEO%nPeriodicVectors.GT.2) THEN
    n0=1
  ELSE
    n0=0
  END IF
  DO iProc = 0,PartMPI%nProcs-1
    IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
    ELSE
      !--- virtually move myself according to periodic vectors in order to find overlapping areas
      !--- 26 possibilities, processes need to work through them in opposite direction in order
      !--- to get matching areas.
      !--- Example for 2D:  I am process #3, I compare myself with #7
      !--- Periodic Vectors are p1 and p2.
      !--- I check p1, p2, p1+p2, p1-p2, -p1+p2, -p1-p2, -p2, -p1
      !--- #7 has to check -p1, -p2, -p1-p2, -p1+p2, p1-p2, p1+p1, p2, p1
      !--- This is done by doing 3 loops from -1 to 1 (for the higher process number)
      !--- or 1 to -1 (for the lower process number) and multiplying
      !--- these numbers to the periodic vectors
      NeighCount = 0  ! -- counter: how often is the process my periodic neighbor?
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      DO k = -SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
        DO m = -SIGN(m0,PartMPI%MyRank-iProc),SIGN(m0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
          DO n = -SIGN(n0,PartMPI%MyRank-iProc),SIGN(n0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
            IF ((k.EQ.0).AND.(m.EQ.0).AND.(n.EQ.0)) CYCLE ! this is not periodic and already done above
            CHECKNEIGHBOR = .TRUE.
            PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
            IF (GEO%nPeriodicVectors.GT.1) THEN
              PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
            END IF
            IF (GEO%nPeriodicVectors.GT.2) THEN
              PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
            END IF
            periodicminmax(1) = localminmax(1) + PeriodicVec(1)
            periodicminmax(2) = localminmax(2) + PeriodicVec(2)
            periodicminmax(3) = localminmax(3) + PeriodicVec(3)
            periodicminmax(4) = localminmax(4) + PeriodicVec(1)
            periodicminmax(5) = localminmax(5) + PeriodicVec(2)
            periodicminmax(6) = localminmax(6) + PeriodicVec(3)
            !--- find overlap
            DO coord = 1,3           ! x y z direction
              maxofmin = MAX(periodicminmax(coord),completeminmax((iProc*6)+coord))
              minofmax = MIN(periodicminmax(3+coord),completeminmax((iProc*6)+3+coord))
              IF (maxofmin.LE.minofmax) THEN         !   overlapping
                TempBorder(1,coord) = maxofmin
                TempBorder(2,coord) = minofmax
              ELSE
                CHECKNEIGHBOR = .FALSE.
              END IF
            END DO
            IF(CHECKNEIGHBOR)THEN
              NeighCount = NeighCount + 1
              TempBorder(:,1) = TempBorder(:,1) - PeriodicVec(1)
              TempBorder(:,2) = TempBorder(:,2) - PeriodicVec(2)
              TempBorder(:,3) = TempBorder(:,3) - PeriodicVec(3)
              TempPeriBord(NeighCount,1:2,1:3) = TempBorder(1:2,1:3)
              PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .TRUE.
            END IF
          END DO
        END DO
      END DO
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = NeighCount
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(1:PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
      END IF
      DO k = 1,NeighCount
        ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3),STAT=allocStat)
        IF (allocStat .NE. 0) THEN
          CALL abort(&
          __STAMP__&
          ,'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
        END IF
        PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3) = TempPeriBord(k,1:2,1:3)
      END DO
    END IF
  END DO
ELSE
  !--- initialize to FALSE for completely non-periodic cases
  DO iProc = 0,PartMPI%nProcs-1
    PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  END DO
END IF

END SUBROUTINE MPIBackgroundMeshInit
#else /*NOT USE_MPI*/
SUBROUTINE PeriodicSourceExchange()
!============================================================================================================================
! Exchange sources in periodic case
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,  ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)         :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i,k,l,m,k2,l2,m2
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,GEO%nPeriodicVectors
  DO k = BGMminX, BGMmaxX
    k2 = k + GEO%PeriodicBGMVectors(1,i)
    DO l = BGMminY, BGMmaxY
      l2 = l + GEO%PeriodicBGMVectors(2,i)
      DO m = BGMminZ, BGMmaxZ
        m2 = m + GEO%PeriodicBGMVectors(3,i)
        IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
          IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
            IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
              BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
              BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO

END SUBROUTINE PeriodicSourceExchange
#endif /*USE_MPI*/

END MODULE MOD_PICDepo_MPI
