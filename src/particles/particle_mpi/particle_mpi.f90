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
USE MOD_ReadInTools       ,ONLY: GETLOGICAL
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) CALL ABORT(__STAMP__,' Particle MPI already initialized!')

! Get flag for ignoring the check and/or abort if the number of global exchange procs is non-symmetric
CheckExchangeProcs = GETLOGICAL('CheckExchangeProcs')
IF(CheckExchangeProcs)THEN
  AbortExchangeProcs = GETLOGICAL('AbortExchangeProcs')
ELSE
  AbortExchangeProcs=.FALSE.
END IF ! .NOT.CheckExchangeProcs

! Get flag for particle latency hiding based on splitting elements in two groups. This first group has particle communication with
! other processors and the second does not.
DoParticleLatencyHiding = GETLOGICAL('DoParticleLatencyHiding')
#if !(PP_TimeDiscMethod==400)
IF(DoParticleLatencyHiding) CALL abort(__STAMP__,'DoParticleLatencyHiding=T not imeplemented for this time disc!')
#endif /*!(PP_TimeDiscMethod==400)*/

ParticleMPIInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

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
USE MOD_Particle_Vars,          ONLY:usevMPF, PDM, UseRotRefFrame
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
! Velocity (rotational reference frame)
IF(UseRotRefFrame) PartCommSize = PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! id of element
PartCommSize   = PartCommSize + 1

IF (useDSMC.AND.(CollisMode.GT.1)) THEN
  PartCommSize = PartCommSize + 2
  IF(DSMC%ElectronicModel.GT.0) PartCommSize   = PartCommSize + 1
END IF
! Simulation with variable particle weights
IF (usevMPF) PartCommSize = PartCommSize+1

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

IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

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
USE MOD_Particle_MPI_Vars,      ONLY:PartMPIExchange
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
                , MPI_COMM_PICLAS                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
 ! IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__&
 !         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


#if defined(IMPA)
SUBROUTINE SendNbOfParticles(doParticle_In)
#else
SUBROUTINE SendNbOfParticles()
#endif /*defined(IMPA)*/
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
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange,PartTargetProc
USE MOD_Particle_MPI_Vars,      ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc, halo_eps_velo
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,PEM,PDM,Species
! variables for parallel deposition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if defined(IMPA)
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
#endif /*defined(IMPA)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(IMPA)
LOGICAL                       :: doPartInExists
#endif /*defined(IMPA)*/
INTEGER                       :: iPart,ElemID, iPolyatMole
INTEGER                       :: iProc,ProcID
!===================================================================================================================================
#if defined(IMPA)
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.
#endif /*defined(IMPA)*/

! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
PartMPIExchange%nPartsSend=0

PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  !         ! Activate phantom/ghost particles
  !         IF(PartSpecies(iPart).LT.0) PDM%ParticleInside(iPart) = .TRUE.

#if defined(IMPA)
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
  ELSE
#endif /*defined(IMPA)*/
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
#if defined(IMPA)
  END IF
#endif /*defined(IMPA)*/

  ! This is already the global ElemID
  ElemID = PEM%GlobalElemID(iPart)
  ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)

  ! Particle on local proc, do nothing
  IF (ProcID.EQ.myRank) CYCLE

  ! Sanity check (fails here if halo region is too small or particle is over speed of light because the time step is too large)
  IF(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID).LT.0)THEN
    IPWRITE (*,*) "GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) =", GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)
    IPWRITE (*,*) "ProcID                                              =", ProcID
    IPWRITE (*,*) "global ElemID                                       =", ElemID
    IPWRITE(UNIT_StdOut,'(I12,A,3(ES25.14E3))') " PartState(1:3,iPart)          =", PartState(1:3,iPart)
    IPWRITE(UNIT_StdOut,'(I12,A,3(ES25.14E3))') " PartState(4:6,iPart)          =", PartState(4:6,iPart)
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " VECNORM(PartState(4:6,iPart)) =", VECNORM(PartState(4:6,iPart))
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " halo_eps_velo                 =", halo_eps_velo
    CALL abort(__STAMP__,'Error: GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) is negative. '//&
                         'The halo region might be too small. Try increasing Particles-HaloEpsVelo! '//&
                         'If this does not help, then maybe the time step is too big. Try reducing ManualTimeStep!')
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
                , MPI_COMM_PICLAS                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE SendNbOfParticles


SUBROUTINE MPIParticleSend(UseOldVecLength)
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
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange,PartCommSize,PartSendBuf,PartRecvBuf,PartTargetProc!,PartHaloElemToProc
USE MOD_Particle_MPI_Vars,       ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Tracking_Vars,  ONLY:TrackingMethod
USE MOD_Particle_Vars,           ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM,PartPosRef,Species
USE MOD_Particle_Vars,           ONLY:UseRotRefFrame,PartVeloRotRef
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
#if defined(LSERK)
USE MOD_Particle_Vars,           ONLY:Pt_temp
#endif
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver_Vars,       ONLY:PartXK,R_PartXK
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
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars,       ONLY:MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL :: UseOldVecLength
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iPos,iProc,jPos, nPartLength
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
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart(2),CounterEnd(2)
REAL(KIND=8)                  :: Rate(2)
#endif /*defined(MEASURE_MPI_WAIT)*/
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

IF (PRESENT(UseOldVecLength)) THEN
  IF (UseOldVecLength) THEN
    nPartLength = PDM%ParticleVecLengthOld
  ELSE
    nPartLength = PDM%ParticleVecLength
  END IF
ELSE
  nPartLength = PDM%ParticleVecLength
END IF

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
  DO iPart=1,nPartLength

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
      !>> particle velocity in rotational reference frame
      IF(UseRotRefFrame) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartVeloRotRef(1:3,iPart)
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
      CALL RemoveParticle(iPart)
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart

  IF(iPos.NE.(MessageSize-MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc))) &
      IPWRITE(*,*) ' error message size', iPos,(MessageSize-MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc))

END DO ! iProc

! 4) Finish Received number of particles
DO iProc=0,nExchangeProcessors-1
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(1))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(1), count_rate=Rate(1))
  CALL SYSTEM_CLOCK(count=CounterStart(2))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(2), count_rate=Rate(2))
  MPIW8TimePart(1)  = MPIW8TimePart(1) + REAL(CounterEnd(1)-CounterStart(1),8)/Rate(1)
  MPIW8CountPart(1) = MPIW8CountPart(1) + 1_8
  MPIW8TimePart(2)  = MPIW8TimePart(2) + REAL(CounterEnd(2)-CounterStart(2),8)/Rate(2)
  MPIW8CountPart(2) = MPIW8CountPart(2) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))

! nullify data on old particle position for safety
DO iPart=1,nPartLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(1:6,iPart) = 0.
  PartSpecies(iPart)   = 0
  IF(UseRotRefFrame) PartVeloRotRef(1:3,iPart) = 0.
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
                , MPI_COMM_PICLAS                                               &
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
                , MPI_COMM_PICLAS                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)

  ! Deallocate sendBuffer after send was successful, see MPIParticleRecv
END DO ! iProc

IF(PDM%UNFPafterMPIPartSend) CALL UpdateNextFreePosition()

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv(DoMPIUpdateNextFreePos)
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
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange,PartCommSize,PartRecvBuf,PartSendBuf
USE MOD_Particle_MPI_Vars      ,ONLY: nExchangeProcessors
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, PartPosRef, Species, LastPartPos
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_Vars          ,ONLY: UseRotRefFrame, PartVeloRotRef
USE MOD_Particle_TimeStep      ,ONLY: GetParticleTimeStep
USE MOD_Particle_Mesh_Vars     ,ONLY: IsExchangeElem
USE MOD_Particle_MPI_Vars      ,ONLY: ExchangeProcToGlobalProc,DoParticleLatencyHiding
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition
#if defined(LSERK)
USE MOD_Particle_Vars          ,ONLY: Pt_temp
#endif
! variables for parallel deposition
!USE MOD_Mesh_Vars              ,ONLY: nGlobalMortarSides
!USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemIsMortar
#if defined(ROS) || defined(IMPA)
USE MOD_Mesh_Vars              ,ONLY: OffSetElem
USE MOD_LinearSolver_Vars      ,ONLY: PartXK,R_PartXK
USE MOD_Particle_Vars          ,ONLY: PartStateN,PartStage,PartDtFrac,PartQ
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
USE MOD_part_tools             ,ONLY: ParticleOnProc, InRotRefFrameCheck
!USE MOD_PICDepo_Tools          ,ONLY: DepositParticleOnNodes
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars,       ONLY:MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL, OPTIONAL             :: DoMPIUpdateNextFreePos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID,jPos, iPart, ElemID
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,0:nExchangeProcessors-1)
INTEGER                       :: MessageSize, nRecvParticles
#if defined(ROS) || defined(IMPA)
INTEGER                       :: iCounter, LocElemID!,iElem
#endif /*ROS or IMPA*/
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, pos_poly, MsgLengthPoly, MsgLengthElec, pos_elec, pos_ambi, MsgLengthAmbi
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart(2),CounterEnd(2)
REAL(KIND=8)                  :: Rate(2)
#endif /*defined(MEASURE_MPI_WAIT)*/
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

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(1))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(1), count_rate=Rate(1))
  MPIW8TimePart(3) = MPIW8TimePart(3) + REAL(CounterEnd(1)-CounterStart(1),8)/Rate(1)
  MPIW8CountPart(3) = MPIW8CountPart(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/
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
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(2))
#endif /*defined(MEASURE_MPI_WAIT)*/
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(2), count_rate=Rate(2))
  MPIW8TimePart(4) = MPIW8TimePart(4) + REAL(CounterEnd(2)-CounterStart(2),8)/Rate(2)
  MPIW8CountPart(4) = MPIW8CountPart(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

  ! place particle information in correct arrays
  !>> correct loop shape
  !>> DO iPart=1,nRecvParticles
  !>> nParts 1 Pos=1..17
  !>> nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1-MsgLengthPoly - MsgLengthElec - MsgLengthAmbi,PartCommSize
    ! find free position in particle array
    nRecv  = nRecv+1
    PartID = GetNextFreePosition(nRecv)

    !>> particle position in physical space
    PartState(1:6,PartID)    = PartRecvBuf(iProc)%content(1+iPos: 6+iPos)
    jPos=iPos+6
    !>> particle position in reference space
    IF(TrackingMethod.EQ.REFMAPPING)THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos: 3+jPos)
      jPos=jPos+3
    END IF
    !>> particle velocity in rotational reference frame
    IF(UseRotRefFrame) THEN
      PDM%InRotRefFrame(PartID) = InRotRefFrameCheck(PartID)
      IF(PDM%InRotRefFrame(PartID)) THEN
        PartVeloRotRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos: 3+jPos)
      ELSE
        PartVeloRotRef(1:3,PartID) = 0.
      END IF
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
    IF(ParticleOnProc(PartID))THEN
      PEM%GlobalElemID(PartID) = LocElemID+OffSetElem
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
    IF (PRESENT(DoMPIUpdateNextFreePos)) THEN
      ElemID = PEM%LocalElemID(PartID)
      IF (ElemID.LT.1) THEN
        CALL abort(__STAMP__,'Particle received in not in proc! Increase halo size! Elem:',PEM%GlobalElemID(PartID))
      END IF
      IF(DoParticleLatencyHiding)THEN
        IF(.NOT.IsExchangeElem(ElemID)) THEN
          IPWRITE(*,*) 'Part Pos + Velo:',PartID,ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc), PartState(1:6,PartID)
          CALL abort(__STAMP__,'Particle received in non exchange elem! Increase halo size! Elem:',PEM%GlobalElemID(PartID))
        END IF
      END IF ! DoParticleLatencyHiding
      IF (useDSMC) THEN
        CALL GetPositionInRefElem(PartState(1:3,PartID),LastPartPos(1:3,PartID),PEM%GlobalElemID(PartID))
      END IF
      IF (useDSMC.OR.usevMPF) THEN
        IF (PEM%pNumber(ElemID).EQ.0) THEN
          PEM%pStart(ElemID) = PartID                    ! Start of Linked List for Particles in Elem
        ELSE
          PEM%pNext(PEM%pEnd(ElemID)) = PartID            ! Next Particle of same Elem (Linked List)
        END IF
        PEM%pEnd(ElemID) = PartID
        ! Number of Particles in Element
        PEM%pNumber(ElemID) = PEM%pNumber(ElemID) + 1
        IF (UseVarTimeStep) THEN
          PartTimeStep(PartID) = GetParticleTimeStep(PartState(1,PartID),PartState(2,PartID),ElemID)
        END IF
      END IF
    END IF
  END DO

END DO ! iProc

IF(PartMPIExchange%nMPIParticles.GT.0) THEN
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
  PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))
END IF
#ifdef CODE_ANALYZE
IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
DO PartID=PDM%ParticleVecLength+1,PDM%maxParticleNumber
  IF (PDM%ParticleInside(PartID)) THEN
    IPWRITE(*,*) PartID,PDM%ParticleVecLength,PDM%maxParticleNumber
    CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=PartID)
  END IF
END DO
#endif

IF(RadialWeighting%PerformCloning) THEN
  ! Checking whether received particles have to be cloned or deleted
  DO iPart = 1,nrecv
    PartID = GetNextFreePosition(iPart-PartMPIExchange%nMPIParticles)
    IF(ParticleOnProc(PartID)) CALL DSMC_2D_RadialWeighting(PartID,PEM%GlobalElemID(PartID))
  END DO
END IF
PartMPIExchange%nMPIParticles = 0
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
    IF(PartMPIInitGroup(iInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(PartMPIInitGroup(iInitRegions)%Comm,iERROR)
    END IF
  END DO ! iInitRegions
END IF

SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPIInitGroup)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)
SDEALLOCATE( ExchangeProcToGlobalProc)
SDEALLOCATE( GlobalProcToExchangeProc)
SDEALLOCATE( PartShiftVector)
SDEALLOCATE( PartTargetProc )

ParticleMPIInitIsDone=.FALSE.
END SUBROUTINE FinalizeParticleMPI
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI
