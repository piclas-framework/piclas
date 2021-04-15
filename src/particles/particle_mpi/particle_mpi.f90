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

MODULE MOD_Particle_MPI
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitParticleMPI
  MODULE PROCEDURE InitParticleMPI
END INTERFACE

#if USE_MPI

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
END INTERFACE

INTERFACE InitEmissionComm
  MODULE PROCEDURE InitEmissionComm
END INTERFACE

INTERFACE IRecvNbOfParticles
  MODULE PROCEDURE IRecvNbOfParticles
END INTERFACE

INTERFACE SendNbOfParticles
  MODULE PROCEDURE SendNbOfParticles
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

!INTERFACE InitHaloMesh
!  MODULE PROCEDURE InitHaloMesh
!END INTERFACE

!INTERFACE ExchangeBezierControlPoints3D
!  MODULE PROCEDURE ExchangeBezierControlPoints3D
!END INTERFACE

PUBLIC :: InitParticleMPI
PUBLIC :: InitEmissionComm
PUBLIC :: InitParticleCommSize
PUBLIC :: SendNbOfParticles
PUBLIC :: IRecvNbOfParticles
PUBLIC :: MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: FinalizeParticleMPI
!PUBLIC :: InitHaloMesh
!PUBLIC :: ExchangeBezierControlPoints3D
#else
PUBLIC :: InitParticleMPI
#endif /*USE_MPI*/

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
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                             :: myRealTestValue
!#if USE_MPI
!INTEGER                         :: color
!#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) &
  CALL ABORT(&
    __STAMP__&
  ,' Particle MPI already initialized!')

#if USE_MPI
!PartMPI%myRank = myRank
!color = 999
!CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,PartMPI%MyRank,PartMPI%COMM,iERROR)
CALL MPI_COMM_DUP (MPI_COMM_WORLD,PartMPI%COMM,iError)
CALL MPI_COMM_RANK(PartMPI%COMM,PartMPI%myRank,iError)
CALL MPI_COMM_SIZE(PartMPI%COMM,PartMPI%nProcs,iError)

IF(PartMPI%nProcs.NE.nProcessors) CALL ABORT(&
    __STAMP__&
    ,' MPI Communicater-size does not match!', IERROR)
PartCommSize   = 0
PartMPI%MPIRoot = .FALSE.
IF(PartMPI%MyRank.EQ.0) PartMPI%MPIRoot=.TRUE.
#else
PartMPI%myRank = 0
PartMPI%nProcs = 1
PartMPI%MPIRoot=.TRUE.
#endif  /*USE_MPI*/

ParticleMPIInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMPI


#if USE_MPI
SUBROUTINE InitParticleCommSize()
!===================================================================================================================================
! get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
! all required features have to be read from the ini-File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
USE MOD_DSMC_Vars,              ONLY:useDSMC, CollisMode, DSMC
USE MOD_Particle_Vars,          ONLY:usevMPF, PDM
USE MOD_Particle_Tracking_vars, ONLY:TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: ALLOCSTAT
!===================================================================================================================================

PartCommSize   = 0
! PartState: position and velocity
PartCommSize   = PartCommSize + 6
! Tracking: Include Reference coordinates
IF(TrackingMethod.EQ.REFMAPPING) PartCommSize=PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! id of element
PartCommSize   = PartCommSize + 1

IF (useDSMC.AND.(CollisMode.GT.1)) THEN
  IF (usevMPF .AND. DSMC%ElectronicModel.GT.0) THEN
    ! vib. , rot and electronic energy and macroparticle factor for each particle
    PartCommSize = PartCommSize + 4
  ELSE IF (usevMPF ) THEN
    ! vib. and rot energy and macroparticle factor for each particle
    PartCommSize = PartCommSize + 3
  ELSE IF (DSMC%ElectronicModel.GT.0) THEN
    ! vib., rot. and electronic energy
    PartCommSize = PartCommSize + 3
  ELSE
    ! vib. and rot. energy
    PartCommSize = PartCommSize + 2
  END IF
ELSE
  ! PIC simulation with MPI
  IF (usevMPF) PartCommSize = PartCommSize+1
END IF

! time integration
#if defined(LSERK)
! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + 6
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1
#endif

! additional stuff for full RK schemes, e.g. implicit and imex RK
#if defined(IMPA) || defined(ROS)
! PartXK
PartCommSize   = PartCommSize + 6
! R_PartXK
PartCommSize   = PartCommSize + 6
! PartQ
PartCommSize   = PartCommSize + 6
! PartDtFrac
PartCommSize   = PartCommSize + 1
! IsNewPart
PartCommSize   = PartCommSize + 1
! FieldAtParticle
PartCommSize   = PartCommSize + 6
! NormVec
PartCommSize   = PartCommSize + 3
! ElementN
PartCommSize   = PartCommSize + 1
! PeriodicMoved
PartCommSize   = PartCommSize + 1
#endif /*IMPA or ROS*/
#if defined(IMPA)
! PartDeltaX
PartCommSize   = PartCommSize + 6
! PartLambdaAccept
PartCommSize   = PartCommSize + 1
! F_PartX0
PartCommSize   = PartCommSize + 6
! F_PartXk
PartCommSize   = PartCommSize + 6
! Norm_F_PartX0, Norm2_F_PartXK, Norm_F_PartXK_old
PartCommSize   = PartCommSize + 3
! DoPartInNewton
PartCommSize   = PartCommSize + 1
! and PartisImplicit
PartCommSize   = PartCommSize + 1
#endif
! if iStage=0, then the PartStateN is not communicated
PartCommSize0  = PartCommSize

ALLOCATE( PartMPIExchange%nPartsSend(4,0:nExchangeProcessors-1)  &
        , PartMPIExchange%nPartsRecv(4,0:nExchangeProcessors-1)  &
        , PartRecvBuf(0:nExchangeProcessors-1)                   &
        , PartSendBuf(0:nExchangeProcessors-1)                   &
        , PartMPIExchange%SendRequest(2,0:nExchangeProcessors-1) &
        , PartMPIExchange%RecvRequest(2,0:nExchangeProcessors-1) &
        , PartTargetProc(1:PDM%MaxParticleNumber)                &
        , STAT=ALLOCSTAT                                         )

IF (ALLOCSTAT.NE.0) &
  CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

PartMPIExchange%nPartsSend=0
PartMPIExchange%nPartsRecv=0

END SUBROUTINE InitParticleCommSize


SUBROUTINE IRecvNbOfParticles()
!===================================================================================================================================
! Open Recv-Buffer for number of received particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartMPIExchange
USE MOD_Particle_MPI_Vars,      ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

PartMPIExchange%nPartsRecv=0
DO iProc=0,nExchangeProcessors-1
  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(:,iProc)                        &
                , 4                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
 ! IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__&
 !         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles(doParticle_In)
!===================================================================================================================================
! This routine sends the number of send particles, for which the following steps are performed:
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! The remaining steps are performed in SendParticles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Part_Tools             ,ONLY: isDepositParticle
USE MOD_DSMC_Vars              ,ONLY: DSMC,SpecDSMC, useDSMC, PolyatomMolDSMC
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI,PartMPIExchange,PartTargetProc
USE MOD_Particle_MPI_Vars,      ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc, halo_eps_velo
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,PEM,PDM,Species
! variables for parallel deposition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: doPartInExists
INTEGER                       :: iPart,ElemID, iPolyatMole
INTEGER                       :: iProc,ProcID
!===================================================================================================================================
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.

! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
PartMPIExchange%nPartsSend=0

PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  !         ! Activate phantom/ghost particles
  !         IF(PartSpecies(iPart).LT.0) PDM%ParticleInside(iPart) = .TRUE.

  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF

  ! This is already the global ElemID
  ElemID = PEM%GlobalElemID(iPart)
  ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)

  ! Particle on local proc, do nothing
  IF (ProcID.EQ.myRank) CYCLE

  ! Sanity check (fails here if halo region is too small)
  IF(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID).LT.0)THEN
    IPWRITE (*,*) "GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) =", GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)
    IPWRITE (*,*) "ProcID                                              =", ProcID
    IPWRITE(UNIT_StdOut,'(I12,A,3(ES25.14E3))') " PartState(4:6,iPart)          =", PartState(4:6,iPart)
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " VECNORM(PartState(4:6,iPart)) =", VECNORM(PartState(4:6,iPart))
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " halo_eps_velo                 =", halo_eps_velo
    CALL abort(&
    __STAMP__&
    ,'Error: GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) is negative. '//&
     'The halo region might be too small. Try increasing Particles-HaloEpsVelo!')
  END IF ! GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID).LT.0

  ! Add particle to target proc count

  PartMPIExchange%nPartsSend(1,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) =  &
    PartMPIExchange%nPartsSend(1,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) + 1
  IF (useDSMC) THEN
    IF ((DSMC%NumPolyatomMolecs.GT.0).OR.(DSMC%ElectronicModel.EQ.2).OR.DSMC%DoAmbipolarDiff) THEN
      IF((DSMC%NumPolyatomMolecs.GT.0).AND.(SpecDSMC(PartSpecies(iPart))%PolyatomicMol)) THEN
        iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
        PartMPIExchange%nPartsSend(2,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) =  &
          PartMPIExchange%nPartsSend(2,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) + PolyatomMolDSMC(iPolyatMole)%VibDOF 
      END IF
      IF ((DSMC%ElectronicModel.EQ.2).AND. & 
          (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized))) THEN
        PartMPIExchange%nPartsSend(3,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) =  &
          PartMPIExchange%nPartsSend(3,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
      END IF
      IF(DSMC%DoAmbipolarDiff.AND.(Species(PartSpecies(iPart))%ChargeIC.GT.0.0)) THEN
        PartMPIExchange%nPartsSend(4,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) =  &
          PartMPIExchange%nPartsSend(4,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) + 3
      END IF
    END IF
  END IF
  PartTargetProc(iPart) = GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)
END DO ! iPart

! 2) send number of send particles
!--- Loop over all neighboring procs. Map local proc ID to global through ExchangeProcToGlobalProc.
!--- Asynchronous communication, just send here and check for success later.
DO iProc=0,nExchangeProcessors-1
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(:,iProc)                        &
                , 4                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE SendNbOfParticles


SUBROUTINE MPIParticleSend()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! first steps are performed in SendNbOfParticles
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! Starting Here:
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
! CAUTION: If particles are sent for deposition, PartTargetProc has the information, if a particle is sent
!          and after the build and wait for number of particles reused to build array with external parts
!          information in PartState,.. can be reused, because they are not overwritten
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DSMC_Vars,               ONLY:useDSMC, CollisMode, DSMC, PartStateIntEn, SpecDSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_DSMC_Vars,               ONLY:ElectronicDistriPart, AmbipolElecVelo
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI,PartMPIExchange,PartCommSize,PartSendBuf,PartRecvBuf,PartTargetProc!,PartHaloElemToProc
USE MOD_Particle_MPI_Vars,       ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Tracking_Vars,  ONLY:TrackingMethod
USE MOD_Particle_Vars,           ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM,PartPosRef,Species
#if defined(LSERK)
USE MOD_Particle_Vars,           ONLY:Pt_temp
#endif
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver_Vars,       ONLY:PartXK,R_PartXK
USE MOD_Particle_Mesh_Vars,      ONLY:ElemToGlobalElemID
USE MOD_Particle_MPI_Vars,       ONLY:PartCommSize0
USE MOD_Particle_Vars,           ONLY:PartStateN,PartStage,PartDtFrac,PartQ
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
USE MOD_Timedisc_Vars,           ONLY:iStage
#endif /*ROS or IMPLICIT*/
#if defined(IMPA)
USE MOD_Particle_Vars,           ONLY:F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton
USE MOD_Particle_Vars,           ONLY:PartDeltaX,PartLambdaAccept
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
#endif /*IMPA*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iPos,iProc,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,0:nExchangeProcessors-1)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles
INTEGER                       :: ALLOCSTAT
#if defined(ROS) || defined(IMPA)
INTEGER                       :: iCounter
#endif /*ROS or IMPA*/
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, MsgRecvLengthPoly, MsgRecvLengthElec, MsgRecvLengthAmbi
INTEGER                       :: MsgLengthPoly(0:nExchangeProcessors-1), pos_poly(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthElec(0:nExchangeProcessors-1), pos_elec(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthAmbi(0:nExchangeProcessors-1), pos_ambi(0:nExchangeProcessors-1)
!===================================================================================================================================

#if defined(ROS)
PartCommSize=PartCommSize0+iStage*6
#endif /*ROS*/
#if defined (IMPA)
PartCommSize=PartCommSize0+iStage*6
#endif /*IMPA*/

!--- Determining the number of additional variables due to VibQuantsPar of polyatomic particles
!--- (size varies depending on the species of particle)
MsgLengthPoly(:) = PartMPIExchange%nPartsSend(2,:) 
MsgLengthElec(:) = PartMPIExchange%nPartsSend(3,:) 
MsgLengthAmbi(:) = PartMPIExchange%nPartsSend(4,:) 

! 3) Build Message
DO iProc=0,nExchangeProcessors-1
  ! allocate SendBuf and prepare to build message
  nSendParticles = PartMPIExchange%nPartsSend(1,iProc)
  iPos           = 0

  ! messageSize is increased for external particles
  IF(nSendParticles.EQ.0) CYCLE

  ! allocate SendBuff of required size
  MessageSize=nSendParticles*PartCommSize

  ! polyatomic molecules follow behind the previous message
  IF (useDSMC) THEN
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      pos_poly(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthPoly(iProc)
    END IF

    IF (DSMC%ElectronicModel.EQ.2) THEN
      pos_elec(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthElec(iProc)
    END IF

    IF (DSMC%DoAmbipolarDiff) THEN
      pos_ambi(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthAmbi(iProc)
    END IF
  END IF

  ! Still no message length, proc can be skipped
  IF (MessageSize.EQ.0) CYCLE

  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) &
    CALL ABORT(__STAMP__,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

  ! build message
  DO iPart=1,PDM%ParticleVecLength

    ! TODO: This seems like a valid check to me, why is it commented out?
    !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE

    ! particle belongs on target proc
    IF (PartTargetProc(iPart).EQ.iProc) THEN
      !>> particle position in physical space
      PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(1:6,iPart)
      jPos=iPos+6
      !>> particle position in reference space
      IF(TrackingMethod.EQ.REFMAPPING) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPosRef(1:3,iPart)
        jPos=jPos+3
      END IF
      !>> particle species
      PartSendBuf(iProc)%content(       1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1

#if defined(LSERK)
      !>> particle acceleration
      PartSendBuf(iProc)%content(1+jPos:6+jPos) = Pt_temp(1:6,iPart)
      !>> flag isNewPart
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(7+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(7+jPos) = 0.
      END IF
      jPos=jPos+7
#endif

#if defined(ROS) || defined(IMPA)
      ! send iStage - 1 messages
      IF (iStage.GT.0) THEN
        IF(iStage.EQ.1) CALL ABORT(__STAMP__,' You should never send particles now!')
        !>>
        PartSendBuf(iProc)%content(1+jpos:6+jpos)        = PartStateN(1:6,iPart)
        DO iCounter=1,iStage-1
          PartSendBuf(iProc)%content(jpos+7+(iCounter-1)*6:jpos+6+iCounter*6) = PartStage(1:6,iCounter,iPart)
        END DO
        jPos=jPos+iStage*6
      ENDIF
      !>>
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartXK(1:6,iPart)
      jPos=jPos+6
      !>>
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = R_PartXK(1:6,iPart)
      jPos=jPos+6
      !>>
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartQ(1:6,iPart)
      jPos=jPos+6
      !>>
      PartSendBuf(iProc)%content(jPos+1) = PartDtFrac(iPart)
      jPos=jPos+1
      !>> flag isNewPart
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.
      END IF
      jPos=jPos+1
      !>> FieldAtParticle
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = FieldAtParticle(1:6,iPart)
      jPos=jPos+6
      !>> particle normVector
      PartSendBuf(iProc)%content(jPos+1:jPos+3) = PEM%NormVec(1:3,iPart)
      PEM%NormVec(1:3,iPart) = 0.
      jPos=jPos+3
      !>> particle elmentN
      PartSendBuf(iProc)%content(jPos+1) = REAL(PEM%GlobalElemID(iPart))
      jPos=jPos+1
      !>> periodic movement
      IF (PEM%PeriodicMoved(iPart)) THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.
      END IF
      PEM%PeriodicMoved(iPart)=.FALSE.
      jPos=jPos+1
#endif /*ROS or IMEX */

#if defined(IMPA)
      ! required for particle newton && closed particle description
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartDeltaX(1:6,iPart)
      IF (PartLambdaAccept(iPart)) THEN
        PartSendBuf(iProc)%content(7+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(7+jPos) = 0.
      END IF
      jPos=jPos+7
      !>>
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = F_PartX0(1:6,iPart)
      jPos=jPos+6
      !>>
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = F_PartXk(1:6,iPart)
      jPos=jPos+6
      !>>
      PartSendBuf(iProc)%content(jPos+1)  = Norm_F_PartX0    (iPart)
      PartSendBuf(iProc)%content(jPos+2)  = Norm_F_PartXK    (iPart)
      PartSendBuf(iProc)%content(jPos+3)  = Norm_F_PartXK_old(iPart)
      !>>
      IF(DoPartInNewton(iPart))THEN
        PartSendBuf(iProc)%content(jPos+4) = 1.0
      ELSE
        PartSendBuf(iProc)%content(jPos+4) = 0.0
      END IF
      jPos=jPos+4
      !>>
      IF(PartIsImplicit(iPart))THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.0
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.0
      END IF
      jPos=jPos+1
#endif /*IMPA*/

      !>> particle element
      PartSendBuf(iProc)%content(    1+jPos) = REAL(PEM%GlobalElemID(iPart),KIND=8)
      jPos=jPos+1

      IF (useDSMC.AND.(CollisMode.GT.1)) THEN
        IF (usevMPF .AND. DSMC%ElectronicModel.GT.0) THEN
          PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn( 1,iPart)
          PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn( 2,iPart)
          PartSendBuf(iProc)%content(3+jPos) = PartMPF(iPart)
          PartSendBuf(iProc)%content(4+jPos) = PartStateIntEn( 3,iPart)
          jPos=jPos+4
        ELSE IF (usevMPF) THEN
          PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn( 1,iPart)
          PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn( 2,iPart)
          PartSendBuf(iProc)%content(3+jPos) = PartMPF(iPart)
          jPos=jPos+3
        ELSE IF (DSMC%ElectronicModel.GT.0) THEN
          PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn( 1,iPart)
          PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn( 2,iPart)
          PartSendBuf(iProc)%content(3+jPos) = PartStateIntEn( 3,iPart)
          jPos=jPos+3
        ELSE
          PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn( 1,iPart)
          PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn( 2,iPart)
          jPos=jPos+2
        END IF
      ELSE
        IF (usevMPF) THEN
          PartSendBuf(iProc)%content(1+jPos) = PartMPF(iPart)
          jPos=jPos+1
        END IF
      END IF

      IF (useDSMC) THEN
        !--- put the polyatomic vibquants per particle at the end of the message
        IF (DSMC%NumPolyatomMolecs.GT.0) THEN
          IF(SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
            PartSendBuf(iProc)%content(pos_poly(iProc)+1:pos_poly(iProc)+PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                                                            = VibQuantsPar(iPart)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
            pos_poly(iProc) = pos_poly(iProc) + PolyatomMolDSMC(iPolyatMole)%VibDOF
          END IF
        END IF

        IF (DSMC%ElectronicModel.EQ.2) THEN
          IF(.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) THEN
            PartSendBuf(iProc)%content(pos_elec(iProc)+1:pos_elec(iProc)+ SpecDSMC(PartSpecies(iPart))%MaxElecQuant) &
                                         = ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(PartSpecies(iPart))%MaxElecQuant)
            pos_elec(iProc) = pos_elec(iProc) + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
          END IF
        END IF

        IF (DSMC%DoAmbipolarDiff) THEN
          IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0)  THEN
            PartSendBuf(iProc)%content(pos_ambi(iProc)+1:pos_ambi(iProc)+ 3) = AmbipolElecVelo(iPart)%ElecVelo(1:3)
            pos_ambi(iProc) = pos_ambi(iProc) + 3
          END IF
        END IF
      END IF

      ! sanity check the message length
      IF(MOD(jPos,PartCommSize).NE.0) THEN
#if defined(ROS) || defined(IMPA)
        IPWRITE(UNIT_stdOut,*)  'iStage',iStage
        IPWRITE(UNIT_stdOut,*)  'PartCommSize0',PartCommSize0,PartCommSize
#else
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
#endif
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL ABORT( __STAMP__,' Particle-wrong sending message size!')
      END IF

      ! increment message position to next element, PartCommSize.EQ.jPos
      iPos=iPos+PartCommSize
      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.

#ifdef IMPA
      DoPartInNewton(  iPart) = .FALSE.
      PartLambdaAccept(iPart) = .TRUE.
      PartIsImplicit(  iPart) = .FALSE.
#endif /*IMPA*/
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart

  IF(iPos.NE.(MessageSize-MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc))) &
      IPWRITE(*,*) ' error message size', iPos,(MessageSize-MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc))

END DO ! iProc

! 4) Finish Received number of particles
DO iProc=0,nExchangeProcessors-1
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))

! nullify data on old particle position for safety
DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(1:6,iPart) = 0.
  PartSpecies(iPart)   = 0
#if defined(LSERK)
  Pt_temp(1:6,iPart)   = 0.
#endif
END DO ! iPart=1,PDM%ParticleVecLength

! 5) Allocate received buffer and open MPI_IRECV
DO iProc=0,nExchangeProcessors-1

  ! skip proc if no particles are to be received
  IF(PartMPIExchange%nPartsRecv(1,iProc).EQ.0) CYCLE

  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize    = nRecvParticles*PartCommSize

  IF (useDSMC) THEN
    ! determine the maximal possible polyatomic addition to the regular recv message
    IF (DSMC%NumPolyatomMolecs.GT.0) THEN
      MsgRecvLengthPoly = PartMPIExchange%nPartsRecv(2,iProc)
      MessageSize       = MessageSize + MsgRecvLengthPoly
    END IF

    IF (DSMC%ElectronicModel.EQ.2) THEN
      MsgRecvLengthElec = PartMPIExchange%nPartsRecv(3,iProc)
      MessageSize       = MessageSize + MsgRecvLengthElec
    END IF

    IF (DSMC%DoAmbipolarDiff) THEN
      MsgRecvLengthAmbi = PartMPIExchange%nPartsRecv(4,iProc)
      MessageSize       = MessageSize + MsgRecvLengthAmbi
    END IF
  END IF

  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received poly particles       ', SUM(PartMPIExchange%nPartsRecv(2,:))
    IPWRITE(*,*) 'sum of total received elec distri particles', SUM(PartMPIExchange%nPartsRecv(3,:))
    IPWRITE(*,*) 'sum of total received ambipolar particles  ', SUM(PartMPIExchange%nPartsRecv(4,:))
    CALL ABORT(__STAMP__,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
  END IF

  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

! 6) Send Particles
DO iProc=0,nExchangeProcessors-1

  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  nSendParticles = PartMPIExchange%nPartsSend(1,iProc)
  MessageSize    = nSendParticles*PartCommSize
  IF (useDSMC) THEN
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      MessageSize = MessageSize + MsgLengthPoly(iProc)
    END IF
    IF (DSMC%ElectronicModel.EQ.2) THEN
      MessageSize = MessageSize + MsgLengthElec(iProc)
    END IF
    IF (DSMC%DoAmbipolarDiff) THEN
      MessageSize = MessageSize + MsgLengthAmbi(iProc)
    END IF
  END IF

  CALL MPI_ISEND( PartSendBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)

  ! Deallocate sendBuffer after send was successful, see MPIParticleRecv
END DO ! iProc

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv()
!===================================================================================================================================
! this routine finishes the communication and places the particle information in the correct arrays. Following steps are performed
! 1) Finish all send requests -> MPI_WAIT
! 2) Finish all recv requests -> MPI_WAIT
! 3) Place particle information in correct arrays
! 4) Deallocate send and recv buffers
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, DSMC, PartStateIntEn, SpecDSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_DSMC_Vars              ,ONLY: ElectronicDistriPart, AmbipolElecVelo
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange,PartCommSize,PartRecvBuf,PartSendBuf!,PartMPI
USE MOD_Particle_MPI_Vars      ,ONLY: nExchangeProcessors
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, PartPosRef, Species
#if defined(LSERK)
USE MOD_Particle_Vars          ,ONLY: Pt_temp
#endif
! variables for parallel deposition
!USE MOD_Mesh_Vars              ,ONLY: nGlobalMortarSides
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemIsMortar
USE MOD_Mesh_Vars              ,ONLY: OffSetElem
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver_Vars      ,ONLY: PartXK,R_PartXK
USE MOD_Particle_Vars          ,ONLY: PartStateN,PartStage,PartDtFrac,PartQ
!USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalElems
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToGlobalElemID
USE MOD_Particle_MPI_Vars      ,ONLY: PartCommSize0
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
USE MOD_Timedisc_Vars          ,ONLY: iStage
#endif /*ROS or IMPA*/
#if defined(IMPA)
USE MOD_Particle_Vars          ,ONLY: F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton
USE MOD_Particle_Vars          ,ONLY: PartDeltaX,PartLambdaAccept
USE MOD_Particle_Vars          ,ONLY: PartIsImplicit
#endif /*IMPA*/
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry          ,ONLY: DSMC_2D_RadialWeighting
!USE MOD_PICDepo_Tools          ,ONLY: DepositParticleOnNodes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID,jPos, iPart, TempNextFreePosition
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,0:nExchangeProcessors-1)
INTEGER                       :: MessageSize, nRecvParticles
#if defined(ROS) || defined(IMPA)
INTEGER                       :: iCounter, LocElemID,iElem
#endif /*ROS or IMPA*/
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, pos_poly, MsgLengthPoly, MsgLengthElec, pos_elec, pos_ambi, MsgLengthAmbi
!===================================================================================================================================

#if defined(ROS)
PartCommSize=PartCommSize0+iStage*6
#endif /*ROS*/
#if defined (IMPA)
PartCommSize=PartCommSize0+iStage*6
#endif /*MPA*/

! wait for all send requests to be successful
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

nRecv=0
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be received
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE

  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize = nRecvParticles*PartCommSize
  IF (useDSMC) THEN
    ! determine the maximal possible polyatomic addition to the regular message
    IF (DSMC%NumPolyatomMolecs.GT.0) THEN
      MsgLengthPoly = PartMPIExchange%nPartsRecv(2,iProc)
    ELSE
      MsgLengthPoly = 0
    END IF

    IF (DSMC%ElectronicModel.EQ.2) THEN
      MsgLengthElec = PartMPIExchange%nPartsRecv(3,iProc)
    ELSE
      MsgLengthElec = 0
    END IF

    IF (DSMC%DoAmbipolarDiff) THEN
      MsgLengthAmbi = PartMPIExchange%nPartsRecv(4,iProc)
    ELSE
      MsgLengthAmbi = 0
    END IF
    pos_poly    = MessageSize
    IF (DSMC%NumPolyatomMolecs.GT.0) MessageSize = MessageSize + MsgLengthPoly
    pos_elec    = MessageSize
    IF (DSMC%ElectronicModel.EQ.2) MessageSize = MessageSize + MsgLengthElec
    pos_ambi    = MessageSize
    IF (DSMC%DoAmbipolarDiff) MessageSize = MessageSize + MsgLengthAmbi
  ELSE
    MsgLengthPoly = 0.
    MsgLengthElec = 0.
    MsgLengthAmbi = 0.
  END IF
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)

  ! place particle information in correct arrays
  !>> correct loop shape
  !>> DO iPart=1,nRecvParticles
  !>> nParts 1 Pos=1..17
  !>> nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1-MsgLengthPoly - MsgLengthElec - MsgLengthAmbi,PartCommSize
    ! find free position in particle array
    nRecv  = nRecv+1
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0) CALL ABORT(__STAMP__,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)

    !>> particle position in physical space
    PartState(1:6,PartID)    = PartRecvBuf(iProc)%content(1+iPos: 6+iPos)
    jPos=iPos+6
    !>> particle position in reference space
    IF(TrackingMethod.EQ.REFMAPPING)THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos: 3+jPos)
      jPos=jPos+3
    END IF
    !>> particle species
    PartSpecies(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1

#if defined(LSERK)
    !>> particle acceleration
    Pt_temp(1:6,PartID)     = PartRecvBuf(iProc)%content( 1+jPos:6+jPos)
    !>> flag isNewPart
    IF (     INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF (INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL ABORT(__STAMP__,'Error with IsNewPart in MPIParticleRecv!')
    END IF
    jPos=jPos+7
#endif

#if defined(ROS) || defined(IMPA)
    IF(iStage.GT.0)THEN
      IF(iStage.EQ.1) CALL ABORT(__STAMP__,' You should never receive particle now!')
      !>>
      PartStateN(1:6,PartID)     = PartRecvBuf(iProc)%content(jpos+1:jpos+6)
      DO iCounter=1,iStage-1
        PartStage(1:6,iCounter,PartID) = PartRecvBuf(iProc)%content(jpos+7+(iCounter-1)*6:jpos+6+iCounter*6)
      END DO
      jPos=jPos+iStage*6
    END IF
    !>>
    PartXK(1:6,PartID)         = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>>
    R_PartXK(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>>
    PartQ(1:6,PartID)          = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>>
    PartDtFrac(PartID) = PartRecvBuf(iProc)%content(jPos+1)
    jPos=jPos+1
    !>> flag isNewPart
    IF (     INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF (INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL ABORT(__STAMP__,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content(1+jPos))
    END IF
    jPos=jPos+1
    !>> FieldAtParticle
    FieldAtParticle(1:6,PartID)  = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>> particle normVector
    PEM%NormVec    (1:3,PartID)  = PartRecvBuf(iProc)%content(jPos+1:jPos+3)
    jPos=jPos+3
    !>> particle elmentN
    LocElemID = INT(PartRecvBuf(iProc)%content(jPos+1),KIND=4)
    IF((LocElemID-OffSetElem.GE.1).AND.(LocElemID-OffSetElem.LE.PP_nElems))THEN
      PEM%GlobalElemID(PartID)=LocElemID-OffSetElem
    ELSE
      ! TODO: This is still broken, halo elems are no longer behind PP_nElems
      CALL ABORT(__STAMP__,'External particles not yet supported with new halo region')

!       PEM%GlobalElemID(PartID)=0
!       DO iElem=PP_nElems+1,nTotalElems
!         IF(ElemToGlobalElemID(iElem).EQ.LocElemID)THEN
!           PEM%GlobalElemID(PartID)=iElem
!           EXIT
!         END IF
!       END DO ! iElem=1,nTotalElems
!       IF(PEM%GlobalElemID(PartID).EQ.0)THEN
!         CALL ABORT(&
!           __STAMP__&
!           ,'Error with IsNewPart in MPIParticleRecv: PEM%GlobalElemID(PartID) = 0!')
!       END IF
    END IF
    ! END TODO
    jPos=jPos+1
    !>> periodic movement
    IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PEM%PeriodicMoved(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PEM%PeriodicMoved(PartID)=.FALSE.
    ELSE
      CALL ABORT(__STAMP__,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content( 1+jPos))
    END IF
    jPos=jPos+1
#endif /*ROS or IMPA*/

#if defined(IMPA)
    !>>
    PartDeltaX(1:6,PartID)     = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    IF (INT(PartRecvBuf(iProc)%content(7+jPos)).EQ.1) THEN
      PartLambdaAccept(PartID)=.TRUE.
    ELSE
      PartLambdaAccept(PartID)=.FALSE.
    END IF
    jPos=jPos+7
    !>>
    F_PartX0(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>>
    F_PartXk(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    !>>
    Norm_F_PartX0    (PartID) = PartRecvBuf(iProc)%content(jPos+1)
    Norm_F_PartXk    (PartID) = PartRecvBuf(iProc)%content(jPos+2)
    Norm_F_PartXk_Old(PartID) = PartRecvBuf(iProc)%content(jPos+3)
    IF(PartRecvBuf(iProc)%content(jPos+4).EQ.1.0)THEN
      DoPartInNewton(PartID) = .TRUE.
    ELSE
      DoPartInNewton(PartID) = .FALSE.
    END IF
    jPos=jPos+4
    !>>
    IF(PartRecvBuf(iProc)%content(jPos+1).EQ.1.0)THEN
        PartIsImplicit(PartID) = .TRUE.
    ELSE
        PartIsImplicit(PartID) = .FALSE.
    END IF
    jPos=jPos+1
#endif /*IMPA*/

    !>> particle element
    PEM%GlobalElemID(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)

!           ! Consider special particles that are communicated with negative species ID (phantom/ghost particles that are deposited on a
!           ! surface with their charge and are then removed here after deposition)
!           IF(PartSpecies(PartID).LT.0)THEN
!             PartSpecies(PartID) = -PartSpecies(PartID) ! make positive species ID again
!             CALL DepositParticleOnNodes(PartID,PartState(1:3,PartID),PEM%GlobalElemID(PartID))
!             PartSpecies(PartID) = 0 ! For safety: nullify the speciesID
!             PDM%ParticleInside(PartID) = .FALSE.
!       #ifdef IMPA
!             PartIsImplicit(PartID) = .FALSE.
!             DoPartInNewton(PartID) = .FALSE.
!       #endif /*IMPA*/
!             CYCLE ! Continue the loop with the next particle
!           END IF ! PartSpecies(PartID).LT.0
    jPos=jPos+1

    IF (useDSMC.AND.(CollisMode.GT.1)) THEN
      IF (usevMPF .AND. DSMC%ElectronicModel.GT.0) THEN
        PartStateIntEn( 1,PartID) = PartRecvBuf(iProc)%content(1+jPos)
        PartStateIntEn( 2,PartID) = PartRecvBuf(iProc)%content(2+jPos)
        PartMPF(PartID)           = PartRecvBuf(iProc)%content(3+jPos)
        PartStateIntEn( 3,PartID) = PartRecvBuf(iProc)%content(4+jPos)
        jPos=jPos+4
      ELSE IF ( usevMPF ) THEN
        PartStateIntEn( 1,PartID) = PartRecvBuf(iProc)%content(1+jPos)
        PartStateIntEn( 2,PartID) = PartRecvBuf(iProc)%content(2+jPos)
        PartMPF(PartID)           = PartRecvBuf(iProc)%content(3+jPos)
        jPos=jPos+3
      ELSE IF ( DSMC%ElectronicModel.GT.0) THEN
        PartStateIntEn( 1,PartID) = PartRecvBuf(iProc)%content(1+jPos)
        PartStateIntEn( 2,PartID) = PartRecvBuf(iProc)%content(2+jPos)
        PartStateIntEn( 3,PartID) = PartRecvBuf(iProc)%content(3+jPos)
        jPos=jPos+3
      ELSE
        PartStateIntEn( 1,PartID) = PartRecvBuf(iProc)%content(1+jPos)
        PartStateIntEn( 2,PartID) = PartRecvBuf(iProc)%content(2+jPos)
        jPos=jPos+2
      END IF
    ELSE
      IF (usevMPF)THEN
        PartMPF(PartID) = PartRecvBuf(iProc)%content(1+jPos)
        jPos=jPos+1
      END IF
    END IF
    IF(MOD(jPos,PartCommSize).NE.0)THEN
#if defined(ROS) || defined(IMPA)
      IPWRITE(UNIT_stdOut,*)  'iStage',iStage
      IPWRITE(UNIT_stdOut,*)  'PartCommSize0',PartCommSize0,PartCommSize
#else
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
#endif
      CALL ABORT(__STAMP__,' Particle-wrong receiving message size!')
    END IF

    IF (useDSMC) THEN
      !--- put the polyatomic vibquants per particle at the end of the message
      IF (DSMC%NumPolyatomMolecs.GT.0) THEN
        IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(PartID))%SpecToPolyArray
          IF(ALLOCATED(VibQuantsPar(PartID)%Quants)) DEALLOCATE(VibQuantsPar(PartID)%Quants)
          ALLOCATE(VibQuantsPar(PartID)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          VibQuantsPar(PartID)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                              = NINT(PartRecvBuf(iProc)%content(pos_poly+1:pos_poly+PolyatomMolDSMC(iPolyatMole)%VibDOF))
          pos_poly = pos_poly + PolyatomMolDSMC(iPolyatMole)%VibDOF
        END IF
      END IF

      IF (DSMC%ElectronicModel.EQ.2) THEN
        IF(.NOT.((SpecDSMC(PartSpecies(PartID))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(PartID))%FullyIonized)) THEN
          IF(ALLOCATED(ElectronicDistriPart(PartID)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(PartID)%DistriFunc)
          ALLOCATE(ElectronicDistriPart(PartID)%DistriFunc(1:SpecDSMC(PartSpecies(PartID))%MaxElecQuant))
          ElectronicDistriPart(PartID)%DistriFunc(1:SpecDSMC(PartSpecies(PartID))%MaxElecQuant) &
                              = PartRecvBuf(iProc)%content(pos_elec+1:pos_elec+SpecDSMC(PartSpecies(PartID))%MaxElecQuant)
          pos_elec = pos_elec + SpecDSMC(PartSpecies(PartID))%MaxElecQuant
        END IF
      END IF

      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
          IF(ALLOCATED(AmbipolElecVelo(PartID)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PartID)%ElecVelo)
          ALLOCATE(AmbipolElecVelo(PartID)%ElecVelo(1:3))
          AmbipolElecVelo(PartID)%ElecVelo(1:3) = PartRecvBuf(iProc)%content(pos_ambi+1:pos_ambi+3)
          pos_ambi = pos_ambi + 3
        END IF
      END IF
    END IF

    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
#if defined(IMPA) || defined(ROS)
    ! only for fully implicit
    !IF(PartIsImplicit(PartID))THEN
    !    ! get LastGlobalElemID in local system
    !    ! element position is verified in armijo if it is needed. This prevents a
    !    ! circular definition.
    !    PEM%LastGlobalElemID(PartID)=-1
    ! END IF
!   ! PEM%StageElement(PartID)= PEM%GlobalElemID(PartID)
!    StagePartPos(PartID,1)  = PartState(1,PartID)
!    StagePartPos(PartID,2)  = PartState(2,PartID)
!    StagePartPos(PartID,3)  = PartState(3,PartID)
#else
    !>> LastGlobalElemID only know to previous proc
    PEM%LastGlobalElemID(PartID) = -888
#endif
  END DO

END DO ! iProc

TempNextFreePosition        = PDM%CurrentNextFreePosition
PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) CALL ABORT(&
    __STAMP__&
    ,' ParticleVecLegnth>MaxParticleNumber due to MPI-communication!')

IF(RadialWeighting%PerformCloning) THEN
  ! Checking whether received particles have to be cloned or deleted
  DO iPart = 1,nrecv
    PartID = PDM%nextFreePosition(iPart+TempNextFreePosition)
    IF ((PEM%GlobalElemID(PartID).GE.1+offSetElem).AND.(PEM%GlobalElemID(PartID).LE.PP_nElems+offSetElem)) &
    CALL DSMC_2D_RadialWeighting(PartID,PEM%GlobalElemID(PartID))
  END DO
END IF

! deallocate send,receive buffer
DO iProc=0,nExchangeProcessors-1
  SDEALLOCATE(PartRecvBuf(iProc)%content)
  SDEALLOCATE(PartSendBuf(iProc)%content)
END DO ! iProc

! last step, nullify number of sent and received particles
PartMPIExchange%nPartsRecv=0
PartMPIExchange%nPartsSend=0
END SUBROUTINE MPIParticleRecv


SUBROUTINE FinalizeParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_MPI_Vars
USE MOD_Particle_Vars,            ONLY:Species,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nInitRegions,iInitRegions,iSpec
!===================================================================================================================================

nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits
END DO ! iSpec
IF(nInitRegions.GT.0) THEN
  DO iInitRegions=1,nInitRegions
    IF(PartMPI%InitGroup(iInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(PartMPI%InitGroup(iInitRegions)%Comm,iERROR)
    END IF
  END DO ! iInitRegions
END IF
IF(PartMPI%COMM.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(PartMPI%COMM,iERROR)

SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPI%InitGroup)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)
SDEALLOCATE( ExchangeProcToGlobalProc)
SDEALLOCATE( GlobalProcToExchangeProc)
SDEALLOCATE( PartShiftVector)
SDEALLOCATE( PartTargetProc )

ParticleMPIInitIsDone=.FALSE.
END SUBROUTINE FinalizeParticleMPI


SUBROUTINE InitEmissionComm()
!===================================================================================================================================
! build emission communicators for particle emission regions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Vars,          ONLY:Species,nSpecies
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
#if !(USE_HDG)
USE MOD_CalcTimeStep,           ONLY:CalcTimeStep
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec,iInit,iNode,iRank
INTEGER                         :: nInitRegions
LOGICAL                         :: RegionOnProc
REAL                            :: xCoords(3,8),lineVector(3),radius,height
REAL                            :: xlen,ylen,zlen
REAL                            :: dt
INTEGER                         :: color,iProc
INTEGER                         :: noInitRank,InitRank
LOGICAL                         :: hasRegion
!===================================================================================================================================

! get number of total init regions
nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits
END DO ! iSpec
IF(nInitRegions.EQ.0) RETURN

! allocate communicators
ALLOCATE( PartMPI%InitGroup(1:nInitRegions))

nInitRegions=0
DO iSpec=1,nSpecies
  RegionOnProc=.FALSE.
  DO iInit=1, Species(iSpec)%NumberOfInits
    nInitRegions=nInitRegions+1
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE ('point')
       xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
       RegionOnProc=PointInProc(xCoords(1:3,1))
    CASE ('line_with_equidistant_distribution')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE ('line')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE('disc')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('photon_SEE_disc')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('2D_landmark','2D_landmark_copy')
       ! Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
       ! for low-temperature partially magnetized plasmas (2019)
       ASSOCIATE( x2 => 1.0e-2       ,& ! m
                  x1 => 0.25e-2      ,& ! m
                  y2 => GEO%ymaxglob ,& ! m
                  y1 => GEO%yminglob ,& ! m
                  z2 => GEO%zmaxglob ,& ! m
                  z1 => GEO%zminglob )
        ! Check all 8 edges
        xCoords(1:3,1) = (/x1,y1,z1/)
        xCoords(1:3,2) = (/x2,y1,z1/)
        xCoords(1:3,3) = (/x1,y2,z1/)
        xCoords(1:3,4) = (/x2,y2,z1/)
        xCoords(1:3,5) = (/x1,y1,z2/)
        xCoords(1:3,6) = (/x2,y1,z2/)
        xCoords(1:3,7) = (/x1,y2,z2/)
        xCoords(1:3,8) = (/x2,y2,z2/)
        RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
    CASE('2D_landmark_neutralization')
      ! Neutralization at const. x-position from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
      ! for low-temperature partially magnetized plasmas (2019)
      ! Check 1st region (emission at fixed x-position x=2.4cm)
      ASSOCIATE( &
                 x2 => 2.4001e-2    ,& ! m
                 x1 => 2.3999e-2    ,& ! m
                 y2 => GEO%ymaxglob ,& ! m
                 y1 => GEO%yminglob ,& ! m
                 z2 => GEO%zmaxglob ,& ! m
                 z1 => GEO%zminglob )
       ! Check all 8 edges
       xCoords(1:3,1) = (/x1,y1,z1/)
       xCoords(1:3,2) = (/x2,y1,z1/)
       xCoords(1:3,3) = (/x1,y2,z1/)
       xCoords(1:3,4) = (/x2,y2,z1/)
       xCoords(1:3,5) = (/x1,y1,z2/)
       xCoords(1:3,6) = (/x2,y1,z2/)
       xCoords(1:3,7) = (/x1,y2,z2/)
       xCoords(1:3,8) = (/x2,y2,z2/)
       RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE

      ! Check 2nd region (left boundary where the exiting particles are counted)
      IF(.NOT.RegionOnProc)THEN
        ASSOCIATE(&
                   x2 => 0.0001e-2    ,& ! m
                   x1 => -0.001e-2    ,& ! m
                   y2 => GEO%ymaxglob ,& ! m
                   y1 => GEO%yminglob ,& ! m
                   z2 => GEO%zmaxglob ,& ! m
                   z1 => GEO%zminglob )
         ! Check all 8 edges
         xCoords(1:3,1) = (/x1,y1,z1/)
         xCoords(1:3,2) = (/x2,y1,z1/)
         xCoords(1:3,3) = (/x1,y2,z1/)
         xCoords(1:3,4) = (/x2,y2,z1/)
         xCoords(1:3,5) = (/x1,y1,z2/)
         xCoords(1:3,6) = (/x2,y1,z2/)
         xCoords(1:3,7) = (/x1,y2,z2/)
         xCoords(1:3,8) = (/x2,y2,z2/)
         RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
      END ASSOCIATE
      END IF ! .NOT.RegionOnProc
    CASE('circle')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('gyrotron_circle')
      Radius=Species(iSpec)%Init(iInit)%RadiusIC+Species(iSpec)%Init(iInit)%RadiusICGyro
      !xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      !ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      xlen=Radius
      ylen=Radius
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      IF(Species(iSpec)%Init(iInit)%ParticleNumber.NE.0)THEN
        lineVector(1:3)=(/0.,0.,Species(iSpec)%Init(iInit)%CylinderHeightIC/)
      ELSE
#if !(USE_HDG)
        dt = CALCTIMESTEP()
#endif /*USE_HDG*/
        lineVector(1:3)= dt* Species(iSpec)%Init(iInit)%VeloIC/Species(iSpec)%Init(iInit)%alpha
        zlen=0.
      END IF
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle_equidistant')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('cuboid')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL ABORT(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC&
                                                           +Species(iSpec)%Init(iInit)%BaseVector2IC

      height= Species(iSpec)%Init(iInit)%CuboidHeightIC
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('sphere')
      ASSOCIATE ( radius => Species(iSpec)%Init(iInit)%RadiusIC        ,&
                  origin => Species(iSpec)%Init(iInit)%BasePointIC(1:3) )
        ! Set the 8 bounding box coordinates depending on the origin and radius
        xCoords(1:3,1)=origin + (/ radius  , -radius , -radius/)
        xCoords(1:3,2)=origin + (/ radius  , radius  , -radius/)
        xCoords(1:3,3)=origin + (/ -radius , radius  , -radius/)
        xCoords(1:3,4)=origin + (/ -radius , -radius , -radius/)
        xCoords(1:3,5)=origin + (/ radius  , -radius , radius /)
        xCoords(1:3,6)=origin + (/ radius  , radius  , radius /)
        xCoords(1:3,7)=origin + (/ -radius , radius  , radius /)
        xCoords(1:3,8)=origin + (/ -radius , -radius , radius /)
      END ASSOCIATE
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder','photon_cylinder')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL ABORT(__STAMP__,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      radius = Species(iSpec)%Init(iInit)%RadiusIC
      ! here no radius, already inclueded
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC-Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                           -Species(iSpec)%Init(iInit)%BaseVector2IC

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC&
                                   +2.0*Species(iSpec)%Init(iInit)%BaseVector2IC

      height= Species(iSpec)%Init(iInit)%CylinderHeightIC
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cell_local')
      RegionOnProc=.TRUE.
    CASE('cuboid_equal')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL ABORT(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)

     !~j CALL ABORT(&
     !~j __STAMP__&
     !~j ,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
    CASE ('cuboid_with_equidistant_distribution')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL ABORT(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE('sin_deviation')
       IF(Species(iSpec)%Init(iInit)%ParticleNumber.NE. &
            (Species(iSpec)%Init(iInit)%maxParticleNumberX * Species(iSpec)%Init(iInit)%maxParticleNumberY &
            * Species(iSpec)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',iSpec,' does not match number of particles in each direction!'
         CALL ABORT(&
         __STAMP__&
         ,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       xCoords(1:3,1) = (/GEO%xminglob,GEO%yminglob,GEO%zminglob/)
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE ('IMD')
       RegionOnProc=.TRUE.
    CASE ('background')
    CASE DEFAULT
      IPWRITE(*,*) 'ERROR: Species ', iSpec, 'of', iInit, 'is using an unknown SpaceIC!'
      CALL ABORT(&
      __STAMP__&
      ,'ERROR: Given SpaceIC is not implemented: '//TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    END SELECT
    ! create new communicator
    color=MPI_UNDEFINED
    !IPWRITE(UNIT_stdOut,*) RegionOnProc
    IF(RegionOnProc) color=nInitRegions!+1
    ! set communicator id
    Species(iSpec)%Init(iInit)%InitComm=nInitRegions

    ! create ranks for RP communicator
    IF(PartMPI%MPIRoot) THEN
      InitRank=-1
      noInitRank=-1
      !myRPRank=0
      iRank=0
      PartMPI%InitGroup(nInitRegions)%MyRank=0
      IF(RegionOnProc) THEN
        InitRank=0
      ELSE
        noInitRank=0
      END IF
      DO iProc=1,PartMPI%nProcs-1
        CALL MPI_RECV(hasRegion,1,MPI_LOGICAL,iProc,0,PartMPI%COMM,MPIstatus,iError)
        IF(hasRegion) THEN
          InitRank=InitRank+1
          CALL MPI_SEND(InitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        ELSE
          noInitRank=noInitRank+1
          CALL MPI_SEND(noInitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        END IF
      END DO
    ELSE
      CALL MPI_SEND(RegionOnProc,1,MPI_LOGICAL,0,0,PartMPI%COMM,iError)
      CALL MPI_RECV(PartMPI%InitGroup(nInitRegions)%MyRank,1,MPI_INTEGER,0,0,PartMPI%COMM,MPIstatus,iError)
    END IF

    ! create new emission communicator
    CALL MPI_COMM_SPLIT(PartMPI%COMM, color, PartMPI%InitGroup(nInitRegions)%MyRank, PartMPI%InitGroup(nInitRegions)%COMM,iError)
    IF(RegionOnProc) CALL MPI_COMM_SIZE(PartMPI%InitGroup(nInitRegions)%COMM,PartMPI%InitGroup(nInitRegions)%nProcs ,iError)
    IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0 .AND. RegionOnProc) &
    WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0,A)') ' Emission-Region,Emission-Communicator:',nInitRegions,' on ',&
    PartMPI%InitGroup(nInitRegions)%nProcs,' procs ('//TRIM(Species(iSpec)%Init(iInit)%SpaceIC)//', iSpec=',iSpec,')'
    IF(PartMPI%InitGroup(nInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0) THEN
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.TRUE.
      ELSE
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.FALSE.
      END IF
      !ALLOCATE(PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1))
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%GroupToComm(PartMPI%InitGroup(nInitRegions)%MyRank)= PartMPI%MyRank
      CALL MPI_ALLGATHER(PartMPI%MyRank,1,MPI_INTEGER&
                        ,PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1)&
                       ,1,MPI_INTEGER,PartMPI%InitGroup(nInitRegions)%COMM,iERROR)
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1)=-1
      DO iRank=0,PartMPI%InitGroup(nInitRegions)%nProcs-1
        PartMPI%InitGroup(nInitRegions)%CommToGroup(PartMPI%InitGroup(nInitRegions)%GroupToComm(iRank))=iRank
      END DO ! iRank
      !IPWRITE(*,*) 'CommToGroup',PartMPI%InitGroup(nInitRegions)%GroupToComm
    END IF
  END DO ! iniT
END DO ! iSpec

END SUBROUTINE InitEmissionComm


FUNCTION BoxInProc(CartNodes,nNodes)
!===================================================================================================================================
! check if bounding box is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: nNodes
REAL,INTENT(IN)   :: CartNodes(1:3,1:nNodes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: BoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

BoxInProc=.FALSE.
! get background of nodes
xmin = HUGE(1)
xmax =-HUGE(1)
ymin = HUGE(1)
ymax =-HUGE(1)
zmin = HUGE(1)
zmax =-HUGE(1)

testval = CEILING((MINVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((MAXVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((MINVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((MAXVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((MINVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((MAXVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) BoxInProc=.TRUE.

END FUNCTION BoxInProc


FUNCTION PointInProc(CartNode)
!===================================================================================================================================
! check if point is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNode(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: PointInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

PointInProc=.FALSE.
! get background of nodes
xmin = HUGE(1)
xmax =-HUGE(1)
ymin = HUGE(1)
ymax =-HUGE(1)
zmin = HUGE(1)
zmax =-HUGE(1)

testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) PointInProc=.TRUE.

END FUNCTION PointInProc

#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI
