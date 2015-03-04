#include "boltzplatz.h"

MODULE MOD_Particle_MPI
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitParticleMPI
  MODULE PROCEDURE InitParticleMPI
END INTERFACE

#ifdef MPI
INTERFACE IRecvNbOfParticles
  MODULE PROCEDURE IRecvNbOfParticles
END INTERFACE

INTERFACE FinalizeParticleMPI
  MODULE PROCEDURE FinalizeParticleMPI
END INTERFACE

INTERFACE MPIParticleSend
  MODULE PROCEDURE MPIParticleSend
END INTERFACE

INTERFACE MPIParticleRecv
  MODULE PROCEDURE MPIParticleRecv
END INTERFACE

INTERFACE InitHaloMesh
  MODULE PROCEDURE InitHaloMesh
END INTERFACE

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
END INTERFACE

PUBLIC :: InitParticleMPI,FinalizeParticleMPI,InitHaloMesh, InitParticleCommSize, IRecvNbOfParticles, MPIParticleSend
PUBLIC :: MPIParticleRecv
#else
PUBLIC :: InitParticleMPI
#endif /*MPI*/

!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                             :: myRealTestValue
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) &
  CALL abort(__STAMP__&
  ,' Particle MPI already initialized!')

#ifdef MPI
PartMPI%COMM   = MPI_COMM_WORLD
PartMPI%myrank = myRank
PartMPI%nProcs = nProcessors
PartCommSize   = 0  
IF(PartMPI%MyRank.EQ.0) THEN
  PartMPI%MPIRoot=.TRUE.
ELSE
  PartMPI%MPIRoot=.FALSE.
END IF
#else
PartMPI%myRank = 0 
PartMPI%nProcs = 1 
PartMPI%MPIRoot=.TRUE.
#endif  /*MPI*/
!! determine datatype length for variables to be sent
!myRealKind = KIND(myRealTestValue)
!IF (myRealKind.EQ.4) THEN
!  myRealKind = MPI_REAL
!ELSE IF (myRealKind.EQ.8) THEN
!  myRealKind = MPI_DOUBLE_PRECISION
!ELSE
!  myRealKind = MPI_REAL
!END IF

ParticleMPIInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMPI


#ifdef MPI
SUBROUTINE InitParticleCommSize()
!===================================================================================================================================
! get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
! all required features have to be read from the ini-File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
USE MOD_DSMC_Vars,           ONLY:useDSMC, CollisMode, DSMC
USE MOD_Particle_Vars,       ONLY:usevMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

PartCommSize   = 0  
IF (useDSMC.AND.(CollisMode.NE.1)) THEN
  IF (usevMPF .AND. DSMC%ElectronicState) THEN
    PartCommSize = 18
  ELSE IF (usevMPF ) THEN
    PartCommSize = 17
  ELSE IF ( DSMC%ElectronicState ) THEN
    PartCommSize = 17
  ELSE
    PartCommSize = 16
  END IF
ELSE
  IF (usevMPF) THEN
    PartCommSize = 15
  ELSE
    PartCommSize = 14
  END IF
END IF

#if ((PP_TimeDiscMethod!=1) && (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6))  /* RK3 and RK4 only */
  PartCommSize = PartCommSize - 6
#endif

ALLOCATE( PartMPIExchange%nPartsSend(PartMPI%nMPINeighbors)    & 
        , PartMPIExchange%nPartsRecv(PartMPI%nMPINeighbors)    &
        , PartRecvBuf(1:PartMPI%nMPINeighbors)                 &
        , PartMPIExchange%SendRequest(2,PartMPI%nMPINeighbors) &
        , PartMPIExchange%RecvRequest(2,PartMPI%nMPINeighbors) )
END SUBROUTINE InitParticleCommSize


SUBROUTINE IRecvNbOfParticles()
!===================================================================================================================================
! Open Recv-Buffer for number of received particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartMPIExchange
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

PartMPIExchange%nPartsRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(iProc)                          &
                , 1                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
!  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(iProc)                          &
!                , 1                                                          &
!                , MPI_INTEGER                                                &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1001                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%RecvRequest(1,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE MPIParticleSend()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! 3) Build Message 
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloToProc, PartCommSize,tMPIMessage, PartRecvBuf
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, Pt_temp
USE MOD_DSMC_Vars,                ONLY:useDSMC, CollisMode, DSMC, PartStateIntEn
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID,iPos,iProc
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles
TYPE(tMPIMessage)             :: SendBuf(0:PartMPI%nMPINeighbors)                          
!===================================================================================================================================

! 1) get number of send particles
PartMPIExchange%nPartsSend=0
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  ElemID=PEM%Element(iPart)
  IF(ElemID.GT.PP_nElems)  PartMPIExchange%nPartsSend(PartHaloToProc(LOCAL_PROC_ID,ElemID))=             &
                                        PartMPIExchange%nPartsSend(PartHaloToProc(LOCAL_PROC_ID,ElemID))+1
END DO ! iPart

! 2) send number of send particles
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(iProc)                          &
                , 1                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )

!  CALL MPI_ISEND( PartMPIExchange%nPartsSend(iProc)                          &
!                , 1                                                          &
!                , MPI_INTEGER                                                &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1001                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%SendRequest(1,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

! 3) Build Message
DO iProc=1, PartMPI%nMPINeighbors
  ! allocate SendBuf
  nSendParticles=PartMPIExchange%nPartsSend(iProc)
  iPos=0
  IF(nSendParticles.EQ.0) CYCLE
  MessageSize=nSendParticles*PartCommSize
  ALLOCATE(SendBuf(iProc)%content(MessageSize))
  ! fill message
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    ElemID=PEM%Element(iPart)
    IF(ElemID.GT.PP_nElems) THEN
      !IF(PartHaloToProc(NATIVE_PROC_ID,ElemID).NE.PartMPI%MPINeighbor(iProc))THEN
      IF(PartHaloToProc(LOCAL_PROC_ID,ElemID).NE.iProc)THEN
        WRITE(*,*) " Warning: Target Rank and rank of cell mismatch!!!"
      END IF
      !iPos=iPos+1
      ! fill content
      SendBuf(iProc)%content(1+iPos:6+iPos) = PartState(iPart,1:6)
      SendBuf(iProc)%content(       7+iPos) = REAL(PartSpecies(iPart))
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* only LSERK */
      SendBuf(iProc)%content(8+iPos:13+iPos) = Pt_temp(iPart,1:6)
      SendBuf(iProc)%content(       14+iPos) = REAL(PartHaloToProc(NATIVE_ELEM_ID,ElemID))
      !IF(.NOT.UseLD) THEN   
        IF (useDSMC.AND.(CollisMode.NE.1)) THEN
          IF (usevMPF .AND. DSMC%ElectronicState) THEN
            SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(17+iPos) = PartMPF(iPart)
            SendBuf(iProc)%content(18+iPos) = PartStateIntEn(iPart, 3)
          ELSE IF (usevMPF) THEN
            SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(17+iPos) = PartMPF(iPart)
          ELSE IF ( DSMC%ElectronicState ) THEN
            SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(17+iPos) = PartStateIntEn(iPart, 3)
          ELSE
            SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)
          END IF
        ELSE
          IF (usevMPF) THEN
            SendBuf(iProc)%content(15,iPos) = PartMPF(iPart)
          END IF
        END IF
      !ELSE ! UseLD == true      =>      useDSMC == true
      !  IF (CollisMode.NE.1) THEN
      !    IF (usevMPF .AND. DSMC%ElectronicState) THEN
      !      SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(17+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(18+iPos) = PartStateIntEn(iPart, 3)
      !      SendBuf(iProc)%content(19+iPos:23+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE IF (usevMPF) THEN
      !      SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(17+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(18+iPos:22+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE IF ( DSMC%ElectronicState ) THEN
      !      SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(17+iPos) = PartStateIntEn(iPart, 3)
      !      SendBuf(iProc)%content(18+iPos:22+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE
      !      SendBuf(iProc)%content(15+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(16+iPos) = PartStateIntEn(iPart, 2)
      !      SendBuf(iProc)%content(17+iPos:21+iPos) = PartStateBulkValues(iPart,1:5)
      !    END IF
      !  ELSE
      !    IF (usevMPF) THEN
      !      SendBuf(iProc)%content(15+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(16+iPos:20+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE
      !      SendBuf(iProc)%content(15+iPos:19+iPos) = PartStateBulkValues(iPart,1:5)
      !    END IF
      !  END IF
      !END IF
#else 
      SendBuf(iProc)%content(    8+iPos) = REAL(PartHaloToProc(NATIVE_ELEM_ID,ElemID))
      !IF(.NOT.UseLD) THEN   
        IF (useDSMC.AND.(CollisMode.NE.1)) THEN
          IF (usevMPF .AND. DSMC%ElectronicState) THEN
            SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(11+iPos) = PartMPF(iPart)
            SendBuf(iProc)%content(12+iPos) = PartStateIntEn(iPart, 3)
          ELSE IF (usevMPF) THEN
            SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(11+iPos) = PartMPF(iPart)
          ELSE IF ( DSMC%ElectronicState ) THEN
            SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
            SendBuf(iProc)%content(11+iPos) = PartStateIntEn(iPart, 3)
          ELSE
            SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
            SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)
          END IF
        ELSE
          IF (usevMPF) THEN
            SendBuf(iProc)%content( 9+iPos) = PartMPF(iPart)
          END IF
        END IF
      !ELSE ! UseLD == true      =>      useDSMC == true
      !  IF (CollisMode.NE.1) THEN
      !    IF (usevMPF .AND. DSMC%ElectronicState) THEN
      !      SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(11+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(12+iPos) = PartStateIntEn(iPart, 3)
      !      SendBuf(iProc)%content(13+iPos:17+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE IF (usevMPF) THEN
      !      SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(11+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(12+iPos:16+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE IF ( DSMC%ElectronicState ) THEN
      !      SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)    
      !      SendBuf(iProc)%content(11+iPos) = PartStateIntEn(iPart, 3)
      !      SendBuf(iProc)%content(12+iPos:15+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE
      !      SendBuf(iProc)%content( 9+iPos) = PartStateIntEn(iPart, 1)
      !      SendBuf(iProc)%content(10+iPos) = PartStateIntEn(iPart, 2)
      !      SendBuf(iProc)%content(11+iPos:15+iPos) = PartStateBulkValues(iPart,1:5)
      !    END IF
      !  ELSE
      !    IF (usevMPF) THEN
      !      SendBuf(iProc)%content( 9+iPos) = PartMPF(iPart)
      !      SendBuf(iProc)%content(10+iPos:14+iPos) = PartStateBulkValues(iPart,1:5)
      !    ELSE
      !      SendBuf(iProc)%content( 9+iPos:13+iPos) = PartStateBulkValues(iPart,1:5)
      !    END IF
      !  END IF
      !END IF
#endif 
      iPos=iPos+MessageSize
      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.  
    END IF ! ElemID is HaloElement
  END DO  ! iPart
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
END DO ! iProc
! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(:))

! 5) Allocate received buffer and open MPI_IRECV
DO iProc=1,PartMPI%nMPINeighbors
  nRecvParticles=PartMPIExchange%nPartsRecv(iProc)
  IF(nRecvParticles.EQ.0) CYCLE
  MessageSize=nRecvParticles*PartCommSize
  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize))
  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )

!  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
!                , MessageSize                                                &
!                , MPI_DOUBLE_PRECISION                                       &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1002                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%RecvRequest(2,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

! 6) Send Particles
DO iProc=1,PartMPI%nMPINeighbors
  nSendParticles=PartMPIExchange%nPartsSend(iProc)
  IF(nSendParticles.EQ.0) CYCLE
  MessageSize=nSendParticles*PartCommSize
  !ALLOCATE(SendBuf(iProc)%content(PartCommSize,nSendParticles))
  CALL MPI_ISEND( SendBuf(iProc)%content                                     &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
!  CALL MPI_ISEND( SendBuf(iProc)%content                                     &
!                , MessageSize                                                &
!                , MPI_DOUBLE_PRECISION                                       &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1002                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%SendRequest(2,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

! deallocate send buffer
DO iProc=1,PartMPI%nMPINeighbors
  SDEALLOCATE(SendBuf(iProc)%content)
END DO ! iProc

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! 3) Build Message 
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloToProc, PartCommSize, PartRecvBuf
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, Pt_temp
USE MOD_DSMC_Vars,                ONLY:useDSMC, CollisMode, DSMC, PartStateIntEn
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize
INTEGER                       :: nRecvParticles
!===================================================================================================================================

nRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  nRecvParticles=PartMPIExchange%nPartsRecv(iProc)
  IF(nRecvParticles.EQ.0) CYCLE
  MessageSize=nRecvParticles*PartCommSize
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
  !DO iPos=1,nRecvParticles
  DO iPos=0,MessageSize-1,PartCommSize
    nRecv=nRecv+1
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0)  CALL abort(__STAMP__&
          ,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)
    PartState(PartID,1:6)   = PartRecvBuf(iProc)%content( 1+iPos: 6+iPos)
    PartSpecies(PartID)     = INT(PartRecvBuf(iProc)%content( 7+iPos))
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
    Pt_temp(PartID,1:6)     = PartRecvBuf(iProc)%content( 8+iPos:13+iPos)
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(14+iPos))
    !IF(.NOT.UseLD) THEN
      IF (useDSMC.AND.(CollisMode.NE.1)) THEN
        IF (usevMPF .AND. DSMC%ElectronicState) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(15+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(16+iPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(17+iPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(18+iPos)
        ELSE IF ( usevMPF) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(15+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(16+iPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(17+iPos)
        ELSE IF ( DSMC%ElectronicState ) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(15+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(16+iPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(17+iPos)
        ELSE
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(15+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(16+iPos)
        END IF
      ELSE
        IF (usevMPF) PartMPF(PartID) = PartRecvBuf(iProc)%content(15+iPos)
      END IF
    !ELSE ! UseLD == true      =>      useDSMC == true
      !IF (CollisMode.NE.1) THEN
        !IF (usevMPF .AND. DSMC%ElectronicState) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(15+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(16+iPos)
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content(17+iPos)
          !PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(18+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(19+iPos:23+iPos)
        !ELSE IF ( usevMPF) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(15+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(16+iPos)
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content(17+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(18+iPos:22+iPos)
        !ELSE IF ( DSMC%ElectronicState ) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(15+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(16+iPos)
          !PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(17+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(18+iPos:22+iPos)
        !ELSE
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(15+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(16+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(17+iPos:21+iPos)
        !END IF
      !ELSE
        !IF (usevMPF) THEN
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content(15+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(16+iPos:20+iPos)
        !ELSE
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(15+iPos:19+iPos)
        !END IF
      !END IF
    !END IF
#else 
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(8+iPos))
    !IF(.NOT.UseLD) THEN
      IF (useDSMC.AND.(CollisMode.NE.1)) THEN
        IF (usevMPF .AND. DSMC%ElectronicState) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content( 9+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(10+iPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(11+iPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(12+iPos)
        ELSE IF ( usevMPF ) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content( 9+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(10+iPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(11+iPos)
        ELSE IF ( DSMC%ElectronicState) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content( 9+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(10+iPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(11+iPos)
        ELSE
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content( 9+iPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(10+iPos)
        END IF
      ELSE
        IF (usevMPF) PartMPF(PartID) = PartRecvBuf(iProc)%content( 9+iPos)
      END IF
    !ELSE ! UseLD == true      =>      useDSMC == true
      !IF (CollisMode.NE.1) THEN
        !IF (usevMPF .AND. DSMC%ElectronicState) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content( 9+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(10+iPos)
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content(11+iPos)
          !PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(12+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(13+iPos:17+iPos)
        !ELSE IF ( usevMPF) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content( 9+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(10+iPos)
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content(11+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(12+iPos:16+iPos)
        !ELSE IF ( DSMC%ElectronicState ) THEN
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content( 9+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(10+iPos)
          !PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(11+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(12+iPos:16+iPos)
        !ELSE
          !PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content( 9+iPos)
          !PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(10+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(11+iPos:15+iPos)
        !END IF
      !ELSE
        !IF (usevMPF) THEN
          !PartMPF(PartID)                 = PartRecvBuf(iProc)%content( 9+iPos)
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(10+iPos:14+iPos)
        !ELSE
          !PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(9+iPos:13+iPos)
        !END IF
      !END IF
    !END IF
#endif
    ! Set Flag for received parts in order to localize them later
    PEM%lastElement(PartID) = -888 
    PDM%ParticleInside(PartID) = .TRUE.
  END DO
  ! be nice: deallocate the receive buffer
  ! deallocate non used array
  DEALLOCATE(PartRecvBuf(iProc)%Content)
END DO ! iProc

PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles



END SUBROUTINE MPIParticleRecv


SUBROUTINE FinalizeParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ParticleMPIInitIsDone=.FALSE.
SDEALLOCATE( PartMPI%isMPINeighbor)
SDEALLOCATE( PartMPI%MPINeighbor )
SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
END SUBROUTINE FinalizeParticleMPI

SUBROUTINE InitHaloMesh()
!===================================================================================================================================
! communicate all direct neighbor sides from master to slave
! has to be called after GetSideType and MPI_INIT of DG solver
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_PreProc
USE MOD_Particle_Surfaces_vars
USE MOD_Mesh_Vars,                  ONLY:NGeo,nSides
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI,PartHaloToProc
USE MOD_Particle_MPI_Halo,          ONLY:IdentifyHaloMPINeighborhood,ExchangeHaloGeometry
USE MOD_Particle_Mesh_Vars,         ONLY:nTotalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 ::BezierSideSize,SendID, iElem
INTEGER                 ::iProc,ALLOCSTAT,iMPINeighbor
LOGICAL                 :: TmpNeigh
INTEGER,ALLOCATABLE     ::SideIndex(:)
!===================================================================================================================================

! funny: should not be required, as sides are build for master and slave sides??
! communicate the MPI Master Sides to Slaves
! all processes have now filled sides and can compute the particles inside the proc region
SendID=1
BezierSideSize=3*(NGeo+1)*(NGeo+1)
DO iNbProc=1,nNbProcs
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =BezierSideSize*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest_Flux(iNbProc),iError)
  END IF
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =BezierSideSize*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest_Flux(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs

DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0) CALL MPI_WAIT(RecRequest_Flux(iNbProc) ,MPIStatus,iError)
END DO !iProc=1,nNBProcs
! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0) CALL MPI_WAIT(SendRequest_Flux(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs

ALLOCATE(SideIndex(1:nSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
,'  Cannot allocate SideIndex!')
SideIndex=0
! check epsilondistance
DO iProc=0,PartMPI%nProcs-1
  IF(iProc.EQ.PartMPI%MyRank) CYCLE
  LOGWRITE(*,*)'  - Identify non-immediate MPI-Neighborhood...'
  !--- AS: identifies which of my node have to be sent to iProc w.r.t. to 
  !        eps vicinity region.
  CALL IdentifyHaloMPINeighborhood(iProc,SideIndex)
  LOGWRITE(*,*)'    ...Done'

  LOGWRITE(*,*)'  - Exchange Geometry of MPI-Neighborhood...'
  CALL ExchangeHaloGeometry(iProc,SideIndex)
  LOGWRITE(*,*)'    ...Done'
  SideIndex(:)=0
END DO 

! Make sure PMPIVAR%MPINeighbor is consistent
DO iProc=0,PartMPI%nProcs-1
  IF (PartMPI%MyRank.EQ.iProc) CYCLE
  IF (PartMPI%MyRank.LT.iProc) THEN
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  ELSE IF (PartMPI%MyRank.GT.iProc) THEN
     CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
     CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
  END IF
  IF (TmpNeigh.NEQV.PartMPI%isMPINeighbor(iProc)) THEN
   WRITE(*,*) 'WARNING: MPINeighbor set to TRUE',PartMPI%MyRank,iProc
   PartMPI%isMPINeighbor(iProc) = .TRUE.
   PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
  END IF
END DO


! fill list with neighbor proc id and add local neighbor id to PartHaloToProc
ALLOCATE( PartMPI%MPINeighbor(PartMPI%nMPINeighbors))
iMPINeighbor=0
DO iProc=0,PartMPI%nProcs-1
  IF(PartMPI%isMPINeighbor(iProc))THEN
    iMPINeighbor=iMPINeighbor+1
    PartMPI%MPINeighbor(iMPINeighbor)=iProc
    DO iElem=PP_nElems+1,nTotalElems
      IF(iProc.EQ.PartHaloToProc(NATIVE_PROC_ID,iElem)) PartHaloToProc(LOCAL_PROC_ID,iElem)=iMPINeighbor
    END DO ! iElem
  END IF
END DO
IPWRITE(*,*) ' List Of Neighbor Procs',  PartMPI%nMPINeighbors,PartMPI%MPINeighbor

IF(iMPINeighbor.NE.PartMPI%nMPINeighbors) CALL abort(&
  __STAMP__&
  , ' Found number of mpi neighbors does not match! ', iMPINeighbor,REAL(PartMPI%nMPINeighbors))

!IF(DepositionType.EQ.'shape_function') THEN
!  PMPIVAR%MPINeighbor(PMPIVAR%iProc) = .TRUE.
!ELSE
!  PMPIVAR%MPINeighbor(PMPIVAR%iProc) = .FALSE.
!END IF

END SUBROUTINE InitHaloMesh

#endif /*MPI*/
END MODULE MOD_Particle_MPI
