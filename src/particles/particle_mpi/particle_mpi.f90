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
USE MOD_Globals_Vars, ONLY: i8
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

PUBLIC :: InitParticleMPI
PUBLIC :: InitParticleCommSize
PUBLIC :: SendNbOfParticles
PUBLIC :: IRecvNbOfParticles
PUBLIC :: MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: FinalizeParticleMPI
#else
PUBLIC :: InitParticleMPI
#endif /*USE_MPI*/

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
! Read-in particle communication parameters
!===================================================================================================================================
SUBROUTINE InitParticleMPI()
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
!===================================================================================================================================
!> get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
!> all required features have to be read from the ini-File
!===================================================================================================================================
SUBROUTINE InitParticleCommSize()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
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

! Simulation with variable particle weights
IF (usevMPF) PartCommSize = PartCommSize+1

! time integration
#if defined(LSERK)
! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + 6
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1
#endif

ALLOCATE( PartMPIExchange%nPartsSend(nPartMPIData,0:nExchangeProcessors-1), STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartMPIExchange%nPartsSend. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartMPIExchange%nPartsRecv(nPartMPIData,0:nExchangeProcessors-1), STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartMPIExchange%nPartsRecv. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartRecvBuf(0:nExchangeProcessors-1)                  , STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartRecvBuf. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartSendBuf(0:nExchangeProcessors-1)                  , STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartSendBuf. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartMPIExchange%SendRequest(2,0:nExchangeProcessors-1), STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartMPIExchange%SendRequest. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartMPIExchange%RecvRequest(2,0:nExchangeProcessors-1), STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartMPIExchange%RecvRequest. ALLOCSTAT',ALLOCSTAT)
ALLOCATE( PartTargetProc(1:PDM%MaxParticleNumber)               , STAT=ALLOCSTAT )
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables, PartTargetProc. ALLOCSTAT',ALLOCSTAT)

PartMPIExchange%nPartsSend=0
PartMPIExchange%nPartsRecv=0

END SUBROUTINE InitParticleCommSize


!===================================================================================================================================
!> Open Recv-Buffer for number of received particles
!===================================================================================================================================
SUBROUTINE IRecvNbOfParticles()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPIExchange, nPartMPIData
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
                , nPartMPIData                                               &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , MPI_COMM_PICLAS                                            &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartMPIExchange%nPartsRecv', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


!===================================================================================================================================
!> This routine sends the number of send particles, for which the following steps are performed:
!> 1) Compute number of Send Particles
!> 2) Perform MPI_ISEND with number of particles
!> The remaining steps are performed in SendParticles
!> 3) Build Message
!> 4) MPI_WAIT for number of received particles
!> 5) Open Receive-Buffer for particle message -> MPI_IRECV
!> 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
SUBROUTINE SendNbOfParticles()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Part_Tools             ,ONLY: isDepositParticle
USE MOD_DSMC_Vars              ,ONLY: DSMC,SpecDSMC, useDSMC, PolyatomMolDSMC,CollisMode
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange,PartTargetProc, nPartMPIData
USE MOD_Particle_MPI_Vars,      ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc, halo_eps_velo
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,PEM,PDM,Species, UseGranularSpecies
USE MOD_Mesh_Vars              ,ONLY: ELEM_RANK
#if USE_HDG
USE MOD_Particle_Boundary_Vars ,ONLY: DoVirtualDielectricLayer
USE MOD_Particle_Vars          ,ONLY: ResetVDLSpecID
#endif/*USE_HDG*/
! variables for parallel deposition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID, iPolyatMole
INTEGER                       :: iProc,ProcID, SpecID, ExchangeProc
!===================================================================================================================================
! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
PartMPIExchange%nPartsSend=0

PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  ! This is already the global ElemID
  ElemID = PEM%GlobalElemID(iPart)
  ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)

  ! Particle on local proc, do nothing
  IF (ProcID.EQ.myRank) CYCLE

  ExchangeProc = GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)

  ! Sanity check (fails here if halo region is too small or particle is over speed of light because the time step is too large)
  IF(ExchangeProc.LT.0)THEN
    IPWRITE (*,*) "GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) =", ExchangeProc
    IPWRITE (*,*) "ProcID                                              =", ProcID
    IPWRITE (*,*) "global ElemID                                       =", ElemID
    IPWRITE(UNIT_StdOut,'(I12,A,3(ES25.14E3))') " PartState(1:3,iPart)            =", PartState(1:3,iPart)
    IPWRITE(UNIT_StdOut,'(I12,A,3(ES25.14E3))') " PartState(4:6,iPart)            =", PartState(4:6,iPart)
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " VECNORM3D(PartState(4:6,iPart)) =", VECNORM3D(PartState(4:6,iPart))
    IPWRITE(UNIT_StdOut,'(I12,A,ES25.14E3)')    " halo_eps_velo                   =", halo_eps_velo
    CALL abort(__STAMP__,'Error: GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID) is negative. '//&
                         'The halo region might be too small. Try increasing Particles-HaloEpsVelo! '//&
                         'If this does not help, then maybe the time step is too big. Try reducing ManualTimeStep!')
  END IF ! ExchangeProc.LT.0

  ! Add particle to target proc count
  PartMPIExchange%nPartsSend(1,ExchangeProc) =  PartMPIExchange%nPartsSend(1,ExchangeProc) + 1
  IF (useDSMC) THEN
    SpecID = PartSpecies(iPart)
#if USE_HDG
    ! Check particle index for VDL particles and reset to original species index
    IF(DoVirtualDielectricLayer) SpecID = ResetVDLSpecID(iPart)
#endif/*USE_HDG*/
    IF ((DSMC%NumPolyatomMolecs.GT.0).OR.(DSMC%ElectronicModel.EQ.2).OR.DSMC%DoAmbipolarDiff) THEN
      IF((DSMC%NumPolyatomMolecs.GT.0).AND.(SpecDSMC(SpecID)%PolyatomicMol)) THEN
        iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
        PartMPIExchange%nPartsSend(2,ExchangeProc) =  PartMPIExchange%nPartsSend(2,ExchangeProc) + PolyatomMolDSMC(iPolyatMole)%VibDOF
      END IF
      IF ((DSMC%ElectronicModel.EQ.2).AND.(.NOT.((Species(SpecID)%InterID.EQ.4).OR.SpecDSMC(SpecID)%FullyIonized))) THEN
        PartMPIExchange%nPartsSend(3,ExchangeProc) =  PartMPIExchange%nPartsSend(3,ExchangeProc) + SpecDSMC(SpecID)%MaxElecQuant
      END IF
      IF(DSMC%DoAmbipolarDiff.AND.(Species(SpecID)%ChargeIC.GT.0.0)) THEN
        PartMPIExchange%nPartsSend(4,ExchangeProc) =  PartMPIExchange%nPartsSend(4,ExchangeProc) + 3
      END IF
    END IF
    IF ((CollisMode.GT.1)) THEN
      IF ((Species(SpecID)%InterID.EQ.2).OR.(Species(SpecID)%InterID.EQ.20)) THEN
        PartMPIExchange%nPartsSend(5,ExchangeProc) =  PartMPIExchange%nPartsSend(5,ExchangeProc) + 2
      END IF
      IF ((DSMC%ElectronicModel.GT.0).AND.(Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
        PartMPIExchange%nPartsSend(6,ExchangeProc) =  PartMPIExchange%nPartsSend(6,ExchangeProc) + 1
      END IF
      IF (UseGranularSpecies.AND.(Species(SpecID)%InterID.EQ.100)) THEN
        PartMPIExchange%nPartsSend(7,ExchangeProc) =  PartMPIExchange%nPartsSend(7,ExchangeProc) + 1
      END IF
    END IF
  END IF
  PartTargetProc(iPart) = ExchangeProc
END DO ! iPart

! 2) send number of send particles
!--- Loop over all neighboring procs. Map local proc ID to global through ExchangeProcToGlobalProc.
!--- Asynchronous communication, just send here and check for success later.
DO iProc=0,nExchangeProcessors-1
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(:,iProc)                        &
                , nPartMPIData                                               &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , MPI_COMM_PICLAS                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartMPIExchange%nPartsSend', IERROR)
END DO ! iProc

END SUBROUTINE SendNbOfParticles


!===================================================================================================================================
!> this routine sends the particles. Following steps are performed
!> first steps are performed in SendNbOfParticles
!> 1) Compute number of Send Particles
!> 2) Perform MPI_ISEND with number of particles
!> Starting Here:
!> 3) Build Message
!> 4) MPI_WAIT for number of received particles
!> 5) Open Receive-Buffer for particle message -> MPI_IRECV
!> 6) Send Particles -> MPI_ISEND
!> CAUTION: If particles are sent for deposition, PartTargetProc has the information, if a particle is sent
!>          and after the build and wait for number of particles reused to build array with external parts
!>          information in PartState,.. can be reused, because they are not overwritten
!===================================================================================================================================
SUBROUTINE MPIParticleSend(UseOldVecLength)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DSMC_Vars,               ONLY:useDSMC, CollisMode, DSMC, PartIntEn, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange,PartCommSize,PartSendBuf,PartRecvBuf,PartTargetProc!,PartHaloElemToProc
USE MOD_Particle_MPI_Vars,       ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Tracking_Vars,  ONLY:TrackingMethod
USE MOD_Particle_Vars,           ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM,PartPosRef,Species
USE MOD_Particle_Vars,           ONLY:UseRotRefFrame,PartVeloRotRef,UseGranularSpecies
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
#if defined(LSERK)
USE MOD_Particle_Vars,           ONLY:Pt_temp
#endif
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars,       ONLY:MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
#if USE_HDG
USE MOD_Particle_Boundary_Vars  ,ONLY: DoVirtualDielectricLayer
USE MOD_Particle_Vars           ,ONLY: ResetVDLSpecID
#endif/*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL :: UseOldVecLength
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iPos,iProc,jPos, nPartLength, SpecID
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles
INTEGER                       :: ALLOCSTAT
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, MsgRecvLengthPoly, MsgRecvLengthElec, MsgRecvLengthAmbi
INTEGER                       :: MsgRecvLengthRotVib, MsgRecvLengthElectronic, MsgRecvLengthSolid
INTEGER                       :: MsgLengthPoly(0:nExchangeProcessors-1), pos_poly(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthElec(0:nExchangeProcessors-1), pos_elec(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthAmbi(0:nExchangeProcessors-1), pos_ambi(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthRotVib(0:nExchangeProcessors-1), pos_rotvib(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthElectronic(0:nExchangeProcessors-1), pos_electronic(0:nExchangeProcessors-1)
INTEGER                       :: MsgLengthSolid(0:nExchangeProcessors-1), pos_solid(0:nExchangeProcessors-1)
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)              :: CounterStart(2),CounterEnd(2)
REAL(KIND=dp)                 :: Rate(2)
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

!--- Determining the number of additional variables due to particle internal data
!--- (size varies depending on the species of particle and utilized model)
MsgLengthPoly(:) = PartMPIExchange%nPartsSend(2,:)
MsgLengthElec(:) = PartMPIExchange%nPartsSend(3,:)
MsgLengthAmbi(:) = PartMPIExchange%nPartsSend(4,:)
MsgLengthRotVib(:) = PartMPIExchange%nPartsSend(5,:)
MsgLengthElectronic(:) = PartMPIExchange%nPartsSend(6,:)
MsgLengthSolid(:) = PartMPIExchange%nPartsSend(7,:)

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

  ! Additional particle data follows behind the previous message
  IF (useDSMC) THEN
    ! Quantum numbers of polyatomic vibrational excitation
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      pos_poly(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthPoly(iProc)
    END IF
    ! Distribution function of electron excitation model
    IF (DSMC%ElectronicModel.EQ.2) THEN
      pos_elec(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthElec(iProc)
    END IF
    ! Ambipolar diffusion
    IF (DSMC%DoAmbipolarDiff) THEN
      pos_ambi(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthAmbi(iProc)
    END IF
    IF (CollisMode.GT.1) THEN
      ! Rotational and vibrational excitation
      pos_rotvib(iProc) = MessageSize
      MessageSize = MessageSize + MsgLengthRotVib(iProc)
      ! Electronic excitation
      IF (DSMC%ElectronicModel.GT.0) THEN
        pos_electronic(iProc) = MessageSize
        MessageSize = MessageSize + MsgLengthElectronic(iProc)
      END IF
      ! Granular / solid particle temperature
      IF (UseGranularSpecies) THEN
        pos_solid(iProc) = MessageSize
        MessageSize = MessageSize + MsgLengthSolid(iProc)
      END IF
    END IF
  END IF

  ! Still no message length, proc can be skipped
  IF (MessageSize.EQ.0) CYCLE

  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

  ! build message
  DO iPart=1,nPartLength
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
      PartSendBuf(iProc)%content(       1+jPos) = REAL(PartSpecies(iPart),KIND=i8)
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

      !>> particle element
      PartSendBuf(iProc)%content(    1+jPos) = REAL(PEM%GlobalElemID(iPart),KIND=i8)
      jPos=jPos+1
      SpecID = PartSpecies(iPart)
      IF (usevMPF) THEN
        PartSendBuf(iProc)%content(1+jPos) = PartMPF(iPart)
        jPos=jPos+1
      END IF

      !>> add additional particle information
      IF (useDSMC) THEN
#if USE_HDG
        ! Check particle index for VDL particles and reset to original species index
        IF(DoVirtualDielectricLayer) SpecID = ResetVDLSpecID(iPart)
#endif/*USE_HDG*/
        !--- add the polyatomic vibquants per particle
        IF (DSMC%NumPolyatomMolecs.GT.0) THEN
          IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
            PartSendBuf(iProc)%content(pos_poly(iProc)+1:pos_poly(iProc)+PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                                                            = PartIntEn(iPart)%QVib(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
            pos_poly(iProc) = pos_poly(iProc) + PolyatomMolDSMC(iPolyatMole)%VibDOF
          END IF
        END IF
        !--- add the electronic distribution function per particle
        IF (DSMC%ElectronicModel.EQ.2) THEN
          IF(.NOT.((Species(SpecID)%InterID.EQ.4).OR.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
            PartSendBuf(iProc)%content(pos_elec(iProc)+1:pos_elec(iProc)+ SpecDSMC(SpecID)%MaxElecQuant) &
                                         = PartIntEn(iPart)%DistriFunc(1:SpecDSMC(SpecID)%MaxElecQuant)
            pos_elec(iProc) = pos_elec(iProc) + SpecDSMC(SpecID)%MaxElecQuant
          END IF
        END IF
        !--- add electron velocity for ambipolar diffusion
        IF (DSMC%DoAmbipolarDiff) THEN
          IF(Species(SpecID)%ChargeIC.GT.0.0)  THEN
            PartSendBuf(iProc)%content(pos_ambi(iProc)+1:pos_ambi(iProc)+ 3) = PartIntEn(iPart)%ElecVelo(1:3)
            pos_ambi(iProc) = pos_ambi(iProc) + 3
          END IF
        END IF
        IF (CollisMode.GT.1) THEN
          !--- add rotational and vibrational energy
          IF ((Species(SpecID)%InterID.EQ.2).OR.(Species(SpecID)%InterID.EQ.20)) THEN
            PartSendBuf(iProc)%content(pos_rotvib(iProc)+1) = PartIntEn(iPart)%ERot(1)
            PartSendBuf(iProc)%content(pos_rotvib(iProc)+2) = PartIntEn(iPart)%EVib(1)
            pos_rotvib(iProc) = pos_rotvib(iProc) + 2
          END IF
          !--- add electronic energy
          IF ((DSMC%ElectronicModel.GT.0).AND.(Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
            PartSendBuf(iProc)%content(pos_electronic(iProc)+1) = PartIntEn(iPart)%EElec(1)
            pos_electronic(iProc) = pos_electronic(iProc) + 1
          END IF
          !--- add temperature of the granular / solid particle
          IF (UseGranularSpecies.AND.(Species(SpecID)%InterID.EQ.100)) THEN
            PartSendBuf(iProc)%content(pos_solid(iProc)+1) = PartIntEn(iPart)%TSolid(1)
            pos_solid(iProc) = pos_solid(iProc) + 1
          END IF
        END IF
      END IF

      ! sanity check the message length
      IF(MOD(jPos,PartCommSize).NE.0) THEN
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL ABORT( __STAMP__,' ERROR in MPIParticleSend: wrong sending message size!')
      END IF

      ! increment message position to next element, PartCommSize.EQ.jPos
      iPos=iPos+PartCommSize
      ! particle is ready for send, now it can be deleted
      CALL RemoveParticle(iPart)
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart

  ! Sanity check
  IF(iPos.NE.(MessageSize - MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc) &
      - MsgLengthRotVib(iProc)-MsgLengthElectronic(iProc)-MsgLengthSolid(iProc))) THEN
      IPWRITE(*,*) ' iPos, MessageSize', iPos, (MessageSize-MsgLengthPoly(iProc)-MsgLengthElec(iProc)-MsgLengthAmbi(iProc)-MsgLengthRotVib(iProc)-MsgLengthElectronic(iProc)-MsgLengthSolid(iProc))
      CALL ABORT(__STAMP__,' ERROR in MPIParticleSend: Unexpected message size!')
  END IF
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=0,nExchangeProcessors-1
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(1))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartMPIExchange%SendRequest', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(1), count_rate=Rate(1))
  CALL SYSTEM_CLOCK(count=CounterStart(2))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartMPIExchange%RecvRequest', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(2), count_rate=Rate(2))
  MPIW8TimePart(1)  = MPIW8TimePart(1) + REAL(CounterEnd(1)-CounterStart(1),8)/Rate(1)
  MPIW8CountPart(1) = MPIW8CountPart(1) + 1_i8
  MPIW8TimePart(2)  = MPIW8TimePart(2) + REAL(CounterEnd(2)-CounterStart(2),8)/Rate(2)
  MPIW8CountPart(2) = MPIW8CountPart(2) + 1_i8
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
    IF (CollisMode.GT.1) THEN
      MsgRecvLengthRotVib = PartMPIExchange%nPartsRecv(5,iProc)
      MessageSize = MessageSize + MsgRecvLengthRotVib
      IF (DSMC%ElectronicModel.GT.0) THEN
        MsgRecvLengthElectronic = PartMPIExchange%nPartsRecv(6,iProc)
        MessageSize = MessageSize + MsgRecvLengthElectronic
      END IF
      IF (UseGranularSpecies) THEN
        MsgRecvLengthSolid = PartMPIExchange%nPartsRecv(7,iProc)
        MessageSize = MessageSize + MsgRecvLengthSolid
      END IF
    END IF
  END IF

  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received poly particles       ', SUM(PartMPIExchange%nPartsRecv(2,:))
    IPWRITE(*,*) 'sum of total received elec distri particles', SUM(PartMPIExchange%nPartsRecv(3,:))
    IPWRITE(*,*) 'sum of total received ambipolar particles  ', SUM(PartMPIExchange%nPartsRecv(4,:))
    IPWRITE(*,*) 'sum of total received molecular particles  ', SUM(PartMPIExchange%nPartsRecv(5,:))
    IPWRITE(*,*) 'sum of total received electronic particles ', SUM(PartMPIExchange%nPartsRecv(6,:))
    IPWRITE(*,*) 'sum of total received solid particles      ', SUM(PartMPIExchange%nPartsRecv(7,:))
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
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartRecvBuf', IERROR)
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
    IF (CollisMode.GT.1) THEN
      MessageSize = MessageSize + MsgLengthRotVib(iProc)
      IF (DSMC%ElectronicModel.GT.0) THEN
        MessageSize = MessageSize + MsgLengthElectronic(iProc)
      END IF
      IF (UseGranularSpecies) THEN
        MessageSize = MessageSize + MsgLengthSolid(iProc)
      END IF
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
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartSendBuf', IERROR)

  ! Deallocate sendBuffer after send was successful, see MPIParticleRecv
END DO ! iProc

IF(PDM%UNFPafterMPIPartSend) CALL UpdateNextFreePosition()

END SUBROUTINE MPIParticleSend


!===================================================================================================================================
!> this routine finishes the communication and places the particle information in the correct arrays. Following steps are performed
!> 1) Finish all send requests -> MPI_WAIT
!> 2) Finish all recv requests -> MPI_WAIT
!> 3) Place particle information in correct arrays
!> 4) Deallocate send and recv buffers
!===================================================================================================================================
SUBROUTINE MPIParticleRecv(DoMPIUpdateNextFreePos)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, DSMC, PartIntEn, SpecDSMC, PolyatomMolDSMC
USE MOD_DSMC_Vars              ,ONLY: ParticleWeighting
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange,PartCommSize,PartRecvBuf,PartSendBuf
USE MOD_Particle_MPI_Vars      ,ONLY: nExchangeProcessors
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, PartPosRef, Species, LastPartPos
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_Vars          ,ONLY: UseRotRefFrame, InRotRefFrame, PartVeloRotRef, UseGranularSpecies
USE MOD_Particle_TimeStep      ,ONLY: GetParticleTimeStep
USE MOD_Particle_Mesh_Vars     ,ONLY: IsExchangeElem
USE MOD_Particle_MPI_Vars      ,ONLY: ExchangeProcToGlobalProc,DoParticleLatencyHiding
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition
#if defined(LSERK)
USE MOD_Particle_Vars          ,ONLY: Pt_temp
#endif
#if USE_HDG
USE MOD_Particle_Boundary_Vars ,ONLY: DoVirtualDielectricLayer
USE MOD_Particle_Vars          ,ONLY: ResetVDLSpecID
#endif/*USE_HDG*/
USE MOD_DSMC_Symmetry          ,ONLY: AdjustParticleWeight
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
LOGICAL, INTENT(IN), OPTIONAL :: DoMPIUpdateNextFreePos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, PartID,jPos, ElemID, SpecID
INTEGER                       :: MessageSize, nRecvParticles
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, pos_poly, MsgLengthPoly, MsgLengthElec, pos_elec, pos_ambi, MsgLengthAmbi
INTEGER                       :: MsgLengthRotVib, pos_rotvib, MsgLengthElectronic, pos_electronic, MsgLengthSolid, pos_solid
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=i8)              :: CounterStart(2),CounterEnd(2)
REAL(KIND=dp)                 :: Rate(2)
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

! wait for all send requests to be successful
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(1))
#endif /*defined(MEASURE_MPI_WAIT)*/
  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error for PartMPIExchange%SendRequest', IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(1), count_rate=Rate(1))
  MPIW8TimePart(3) = MPIW8TimePart(3) + REAL(CounterEnd(1)-CounterStart(1),8)/Rate(1)
  MPIW8CountPart(3) = MPIW8CountPart(3) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/
END DO ! iProc

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
    IF (CollisMode.GT.1) THEN
      MsgLengthRotVib = PartMPIExchange%nPartsRecv(5,iProc)
      IF (DSMC%ElectronicModel.GT.0) THEN
        MsgLengthElectronic = PartMPIExchange%nPartsRecv(6,iProc)
      ELSE
        MsgLengthElectronic = 0
      END IF
      IF (UseGranularSpecies) THEN
        MsgLengthSolid = PartMPIExchange%nPartsRecv(7,iProc)
      ELSE
        MsgLengthSolid = 0
      END IF
    ELSE
      MsgLengthRotVib = 0
      MsgLengthElectronic = 0
      MsgLengthSolid = 0
    END IF
    pos_poly    = MessageSize
    IF (DSMC%NumPolyatomMolecs.GT.0) MessageSize = MessageSize + MsgLengthPoly
    pos_elec    = MessageSize
    IF (DSMC%ElectronicModel.EQ.2) MessageSize = MessageSize + MsgLengthElec
    pos_ambi    = MessageSize
    IF (DSMC%DoAmbipolarDiff) MessageSize = MessageSize + MsgLengthAmbi
    pos_rotvib = MessageSize
    IF (CollisMode.GT.1) MessageSize = MessageSize + MsgLengthRotVib
    pos_electronic = MessageSize
    IF ((CollisMode.GT.1).AND.(DSMC%ElectronicModel.GT.0)) MessageSize = MessageSize + MsgLengthElectronic
    pos_solid = MessageSize
    IF ((CollisMode.GT.1).AND.(UseGranularSpecies)) MessageSize = MessageSize + MsgLengthSolid
  ELSE
    MsgLengthPoly = 0.
    MsgLengthElec = 0.
    MsgLengthAmbi = 0.
    MsgLengthRotVib = 0
    MsgLengthElectronic = 0
    MsgLengthSolid = 0
  END IF
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart(2))
#endif /*defined(MEASURE_MPI_WAIT)*/
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),MPI_STATUS_IGNORE,IERROR)
#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd(2), count_rate=Rate(2))
  MPIW8TimePart(4) = MPIW8TimePart(4) + REAL(CounterEnd(2)-CounterStart(2),8)/Rate(2)
  MPIW8CountPart(4) = MPIW8CountPart(4) + 1_i8
#endif /*defined(MEASURE_MPI_WAIT)*/

  ! place particle information in correct arrays
  !>> correct loop shape
  !>> DO iPart=1,nRecvParticles
  !>> nParts 1 Pos=1..17
  !>> nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1-MsgLengthPoly - MsgLengthElec - MsgLengthAmbi-MsgLengthRotVib-MsgLengthElectronic-MsgLengthSolid,PartCommSize
    ! find free position in particle array
    PartID = GetNextFreePosition()

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
      InRotRefFrame(PartID) = InRotRefFrameCheck(PartID)
      IF(InRotRefFrame(PartID)) THEN
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

    !>> particle element
    PEM%GlobalElemID(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1

    SpecID = PartSpecies(PartID)
    IF (usevMPF)THEN
      PartMPF(PartID) = PartRecvBuf(iProc)%content(1+jPos)
      jPos=jPos+1
    END IF
    IF(MOD(jPos,PartCommSize).NE.0)THEN
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
      CALL ABORT(__STAMP__,' ERROR in MPIParticleRecv: wrong receiving message size!')
    END IF

    IF (useDSMC) THEN
#if USE_HDG
      ! Check particle index for VDL particles and reset to original species index
      IF(DoVirtualDielectricLayer) SpecID = ResetVDLSpecID(PartID)
#endif/*USE_HDG*/
      !--- put the polyatomic vibquants per particle at the end of the message
      IF (DSMC%NumPolyatomMolecs.GT.0) THEN
        IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
          IF(ALLOCATED(PartIntEn(PartID)%QVib)) DEALLOCATE(PartIntEn(PartID)%QVib)
          ALLOCATE(PartIntEn(PartID)%QVib(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          PartIntEn(PartID)%QVib(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                              = NINT(PartRecvBuf(iProc)%content(pos_poly+1:pos_poly+PolyatomMolDSMC(iPolyatMole)%VibDOF))
          pos_poly = pos_poly + PolyatomMolDSMC(iPolyatMole)%VibDOF
        END IF
      END IF

      IF (DSMC%ElectronicModel.EQ.2) THEN
        IF(.NOT.((Species(SpecID)%InterID.EQ.4).OR.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
          IF(ALLOCATED(PartIntEn(PartID)%DistriFunc)) DEALLOCATE(PartIntEn(PartID)%DistriFunc)
          ALLOCATE(PartIntEn(PartID)%DistriFunc(1:SpecDSMC(SpecID)%MaxElecQuant))
          PartIntEn(PartID)%DistriFunc(1:SpecDSMC(SpecID)%MaxElecQuant) &
                              = PartRecvBuf(iProc)%content(pos_elec+1:pos_elec+SpecDSMC(SpecID)%MaxElecQuant)
          pos_elec = pos_elec + SpecDSMC(SpecID)%MaxElecQuant
        END IF
      END IF

      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
          IF(ALLOCATED(PartIntEn(PartID)%ElecVelo)) DEALLOCATE(PartIntEn(PartID)%ElecVelo)
          ALLOCATE(PartIntEn(PartID)%ElecVelo(1:3))
          PartIntEn(PartID)%ElecVelo(1:3) = PartRecvBuf(iProc)%content(pos_ambi+1:pos_ambi+3)
          pos_ambi = pos_ambi + 3
        END IF
      END IF

      IF (CollisMode.GT.1) THEN
        IF ((Species(SpecID)%InterID.EQ.2).OR.(Species(SpecID)%InterID.EQ.20)) THEN
          IF(.NOT.ALLOCATED(PartIntEn(PartID)%ERot)) ALLOCATE(PartIntEn(PartID)%ERot(1))
          IF(.NOT.ALLOCATED(PartIntEn(PartID)%EVib)) ALLOCATE(PartIntEn(PartID)%EVib(1))
          PartIntEn(PartID)%ERot(1) = PartRecvBuf(iProc)%content(pos_rotvib+1)
          PartIntEn(PartID)%EVib(1) = PartRecvBuf(iProc)%content(pos_rotvib+2)
          pos_rotvib = pos_rotvib + 2
        END IF
        IF ((DSMC%ElectronicModel.GT.0).AND.(Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
          IF(.NOT.ALLOCATED(PartIntEn(PartID)%EElec)) ALLOCATE(PartIntEn(PartID)%EElec(1))
          PartIntEn(PartID)%EElec(1) = PartRecvBuf(iProc)%content(pos_electronic+1)
          pos_electronic = pos_electronic + 1
        END IF
        IF (UseGranularSpecies.AND.(Species(SpecID)%InterID.EQ.100)) THEN
          IF(.NOT.ALLOCATED(PartIntEn(PartID)%TSolid)) ALLOCATE(PartIntEn(PartID)%TSolid(1))
          PartIntEn(PartID)%TSolid(1) = PartRecvBuf(iProc)%content(pos_solid+1)
          pos_solid = pos_solid + 1
        END IF
      END IF
    END IF
    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
    ! Particle weighting: insert clones either instantly or store for later insertion
    IF(ParticleWeighting%PerformCloning) CALL AdjustParticleWeight(PartID,PEM%GlobalElemID(PartID))

    !>> LastGlobalElemID only know to previous proc
    PEM%LastGlobalElemID(PartID) = -888
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


!===================================================================================================================================
!> Finalize particle MPI communication arrays
!===================================================================================================================================
SUBROUTINE FinalizeParticleMPI()
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
