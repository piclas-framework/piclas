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


#if IMPA
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepByImplicitRK
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStepByImplicitRK()
!===================================================================================================================================
! IMEX time integrator
! ESDIRK for particles
! ERK for field
! from: Kennedy & Carpenter 2003
! Additive Runge-Kutta schemes for convection-diffusion-reaction equations
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
! test timedisc for implicit particle treatment || highly relativistic
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iter,iStage, nRKStages,time
USE MOD_TimeDisc_Vars          ,ONLY: ERK_a,ESDIRK_a,RK_b,RK_c
USE MOD_LinearSolver_Vars      ,ONLY: ImplicitSource, DoPrintConvInfo
USE MOD_DG_Vars                ,ONLY: U
#if USE_HDG
USE MOD_HDG                    ,ONLY: HDG
#else /*pure DG*/
USE MOD_DG_Vars                ,ONLY: Ut,Un
USE MOD_DG                     ,ONLY: DGTimeDerivative_weakForm
USE MOD_Predictor              ,ONLY: Predictor,StorePredictor
USE MOD_LinearSolver_Vars      ,ONLY: LinSolverRHS,FieldStage
USE MOD_Equation               ,ONLY: DivCleaningDamping
USE MOD_Equation               ,ONLY: CalcSource
#ifdef maxwell
USE MOD_Precond                ,ONLY: BuildPrecond
USE MOD_Precond_Vars           ,ONLY: UpdatePrecond
USE MOD_TimeDisc_Vars          ,ONLY: dt_old
#endif /*maxwell*/
#endif /*USE_HDG*/
USE MOD_Newton                 ,ONLY: ImplicitNorm,FullNewton
#ifdef PARTICLES
USE MOD_TimeDisc_Vars          ,ONLY: RK_fillSF
USE MOD_Globals_Vars           ,ONLY: c2_inv
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_PICDepo_Vars           ,ONLY: PartSource,DoDeposition
USE MOD_LinearSolver_Vars      ,ONLY: ExplicitPartSource
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_LinearSolver_Vars      ,ONLY: DoUpdateInStage,PartXk
USE MOD_Predictor              ,ONLY: PartPredictor,PredictorType
USE MOD_Particle_Vars          ,ONLY: PartIsImplicit,PartLorentzType,PartDtFrac &
                                      ,DoForceFreeSurfaceFlux,PartStateN,PartStage,PartQ,DoSurfaceFlux,PEM,PDM  &
                                      , Pt,LastPartPos,DelayTime,PartState,PartMeshHasReflectiveBCs,PartDeltaX
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_part_RHS               ,ONLY: PartVeloToImp
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_Particle_Tracing       ,ONLY: ParticleTracing
USE MOD_Particle_RefTracking   ,ONLY: ParticleRefTracking
USE MOD_Particle_TriaTracking  ,ONLY: ParticleTriaTracking
USE MOD_Particle_Tracking_vars ,ONLY: TrackingMethod
USE MOD_ParticleSolver         ,ONLY: ParticleNewton, SelectImplicitParticles
USE MOD_Part_RHS               ,ONLY: PartRHS
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange
#endif /*USE_MPI*/
USE MOD_PIC_Analyze            ,ONLY: CalcDepositedCharge
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
#ifdef CODE_ANALYZE
USE MOD_Particle_Mesh_Vars     ,ONLY: Geo
USE MOD_Particle_Tracking      ,ONLY: ParticleSanityCheck
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#ifdef maxwell
USE MOD_Precond_Vars           ,ONLY: UpdatePrecondLB
#endif /*maxwell*/
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: tstage
! implicit
REAL               :: alpha
REAL               :: sgamma
! particle surface flux
! RK counter
INTEGER            :: iCounter
#ifdef PARTICLES
INTEGER            :: iElem,i,j,k
REAL               :: dtloc,RandVal
INTEGER            :: iPart,nParts
REAL               :: LorentzFacInv
!LOGICAL            :: NoInterpolation ! fields cannot be interpolated, because particle is "outside", hence, fields and
!                                      ! forces of previous stage are used
REAL               :: n_loc(1:3)
LOGICAL            :: reMap
#ifdef CODE_ANALYZE
REAL               :: IntersectionPoint(1:3)
INTEGER            :: ElemID
LOGICAL            :: ishit
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/
#if USE_LOADBALANCE
REAL               :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
#ifdef maxwell
LOGICAL            :: UpdatePrecondLoc
#endif /*maxwell*/
!===================================================================================================================================

#if !(USE_HDG)
#ifdef maxwell
! caution hard coded
IF (iter==0)THEN
  CALL BuildPrecond(time,time,0,RK_b(nRKStages),dt)
  dt_old=dt
ELSE
  UpdatePrecondLoc=.FALSE.
#if USE_LOADBALANCE
  IF(UpdatePrecondLB) UpdatePrecondLoc=.TRUE.
  UpdatePrecondLb=.FALSE.
#endif /*USE_LOADBALANCE*/
  IF(UpdatePrecond)THEN
    IF(dt.NE.dt_old) UpdatePrecondLoc=.TRUE.
    dt_old=dt
  END IF
  IF(UpdatePrecondLoc) CALL BuildPrecond(time,time,0,RK_b(nRKStages),dt)
END IF
#endif /*maxwell*/
#endif /*DG*/

#ifdef PARTICLES
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
! particle locating
! at the wrong position? depending on how we do it...
IF (time.GE.DelayTime) THEN
  CALL ParticleInserting() ! do not forget to communicate the emitted particles ... for shape function
END IF
#if USE_LOADBALANCE
CALL LBSplitTime(LB_EMISSION,tLBStart)
#endif /*USE_LOADBALANCE*/
! select, if particles are treated implicitly or explicitly
CALL SelectImplicitParticles()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------

tStage=time
#ifdef PARTICLES
! simulation with delay-time, compute the
IF(DelayTime.GT.0.)THEN
  IF((iter.EQ.0).AND.(time.LT.DelayTime))THEN
    ! perform normal deposition
    CALL Deposition()
  END IF
END IF

! compute source of first stage for Maxwell solver
IF (time.GE.DelayTime) THEN
  ! if we call it correctly, we may save here work between different RK-stages
  ! because of emmision and UpdateParticlePosition
  CALL Deposition()
END IF

ImplicitSource=0.
ExplicitPartSource=0.

#if !(USE_HDG)
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL CalcSource(tStage,1.,ImplicitSource)
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#else
! set required for fluid model and HDG, because field may have changed due to different particle distribution
! HDG time measured in HDG
CALL HDG(time,U,iter)
#endif


IF(time.GE.DelayTime)THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! velocity to impulse
  CALL PartVeloToImp(VeloToImp=.TRUE.)
  PartStateN(1:6,1:PDM%ParticleVecLength)=PartState(1:6,1:PDM%ParticleVecLength)
  PEM%GlobalElemID(1:PDM%ParticleVecLength)  =PEM%GlobalElemID(1:PDM%ParticleVecLength)
  IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,1:PDM%ParticleVecLength) =0.
  PEM%PeriodicMoved(1:PDM%ParticleVecLength) = .FALSE.
  IF(iter.EQ.0)THEN ! caution with emission: fields should also be interpolated to new particles, this is missing
                    ! or should be done directly during emission...
    ! should be already be done
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(PartIsImplicit(iPart))THEN
        CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
      END IF ! ParticleIsImplicit
      PDM%IsNewPart(iPart)=.FALSE.
      !PEM%GlobalElemID(iPart) = PEM%GlobalElemID(iPart)
    END DO ! iPart
  END IF
  RKdtFracTotal=0.
  RKdtFrac     =0.
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

  ! surface flux
  IF(DoSurfaceFlux)THEN
    RKdtFrac      = 1.0 ! RK_c(iStage)
    RKdtFracTotal = 1.0 ! RK_c(iStage)
    RK_fillSF     = 1.0
    IF(DoPrintConvInfo) nParts=0
    ! Measurement as well in ParticleSurfaceFlux
    CALL ParticleSurfaceflux()
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    ! compute emission for all particles during dt
    DO iPart=1,PDM%ParticleVecLength
      IF(PDM%ParticleInside(iPart))THEN
        IF(PDM%IsNewPart(iPart))THEN
          IF(DoPrintConvInfo) nParts=nParts+1
          ! nullify
          PartStage(1:6,:,iPart)=0.
          ! f(u^n) for position
          ! CAUTION: position in reference space has to be computed during emission for implicit particles
          ! interpolate field at surface position
          CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
          ! RHS at interface
          CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
          ! f(u^n) for velocity
          IF(.NOT.DoForceFreeSurfaceFlux) PartStage(4:6,1,iPart)=Pt(1:3,iPart)
          ! position NOT known but we backup the state
          PartStateN(1:3,iPart) = PartState(1:3,iPart)
          ! initial velocity equals velocity of surface flux
          PartStateN(4:6,iPart) = PartState(4:6,iPart)
          ! gives entry point into domain
          PEM%GlobalElemID(iPart)      = PEM%GlobalElemID(iPart)
          IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart)   = 0.
          PEM%PeriodicMoved(iPart) = .FALSE.
          CALL RANDOM_NUMBER(RandVal)
          PartDtFrac(iPart)=RandVal
          PartIsImplicit(iPart)=.TRUE.
          ! particle crosses surface at time^n + (1.-RandVal)*dt
          ! for all stages   t_Stage =< time^n + (1.-RandVal)*dt particle is outside of domain
          ! for              t_Stage >  time^n + (1.-RandVal)*dt particle is in domain and can be advanced in time
        ELSE
          ! set DtFrac to unity
          PartDtFrac(iPart)=1.0
        END IF ! IsNewPart
      ELSE
        PartIsImplicit(iPart)=.FALSE.
      END IF ! ParticleInside
    END DO ! iPart
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/
    IF(DoPrintConvInfo)THEN
#if USE_MPI
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,nParts,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
      ELSE
        CALL MPI_REDUCE(nParts       ,iPart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
      END IF
#endif /*USE_MPI*/
      SWRITE(UNIT_StdOut,'(A,I10)') ' SurfaceFlux-Particles: ',nParts
    END IF
  END IF
END IF ! time.GE. DelayTime
#endif /*PARTICLES*/

#if !(USE_HDG)
! LoadBalance Time-Measurement is in DGTimeDerivative_weakForm
IF(iter.EQ.0) CALL DGTimeDerivative_weakForm(time, time, 0,doSource=.FALSE.)
iStage=0
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StorePredictor()
! copy U to U^0
Un = U
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*DG*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to 6
! ----------------------------------------------------------------------------------------------------------------------------------
DO iStage=2,nRKStages
  ! time of current stage
  tStage = time + RK_c(iStage)*dt
  alpha  = ESDIRK_a(iStage,iStage)*dt
  sGamma = 1.0/alpha
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_StdOut,'(A)')    '-----------------------------'
    SWRITE(UNIT_StdOut,'(A,I2)') 'istage:',istage
  END IF
#if !(USE_HDG)
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! store predictor
  CALL StorePredictor()

  ! compute the f(u^s-1)
  ! compute the f(u^s-1)
  ! DG-solver, Maxwell's equations
  FieldStage (:,:,:,:,:,iStage-1) = Ut + ImplicitSource
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*DG*/

  ! and particles
#ifdef PARTICLES
  IF (time.GE.DelayTime) THEN
    IF(DoPrintConvInfo)THEN
      SWRITE(UNIT_StdOut,'(A)') '-----------------------------'
      SWRITE(UNIT_StdOut,'(A)') ' compute last stage value.   '
    END IF
    IF(iStage.EQ.2)THEN
      CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.)
    ELSE
      CALL CountPartsPerElem(ResetNumberOfParticles=.FALSE.)
    END IF
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    ! normal
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        PartIsImplicit(iPart)=.FALSE.
        PDM%IsNewPart(iPart) =.FALSE.
        CYCLE
      END IF
      !IF(PDM%IsNewPart(iPart)) CYCLE ! ignore surface flux particles
      IF(PartIsImplicit(iPart))THEN
        ! caution, implicit is already back-roated for the computation of Pt and Velocity/Momentum
        reMap=.FALSE.
        IF(PartMeshHasReflectiveBCs)THEN
          IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
            PEM%NormVec(1:3,iPart)=0.
            reMap=.TRUE.
          END IF
        END IF
        IF(PEM%PeriodicMoved(iPart))THEN
          PEM%PeriodicMoved(iPart)=.FALSE.
          Remap=.TRUE.
        END IF
        IF(reMap)THEN
          ! recompute the particle position AFTER movement
          ! can be theroretically outside of the mesh, because it is the not-shifted particle position,velocity
          PartState(1:6,iPart)=PartXK(1:6,iPart)+PartDeltaX(1:6,iPart)
        END IF
        IF(PartLorentzType.NE.5)THEN
          PartStage(1:3,iStage-1,iPart) = PartState(4:6,iPart)
          PartStage(4:6,iStage-1,iPart) = Pt(1:3,iPart)
        ELSE
          LorentzFacInv=1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv
          LorentzFacInv=1.0/SQRT(LorentzFacInv)
          PartStage(1  ,iStage-1,iPart) = PartState(4,iPart) * LorentzFacInv
          PartStage(2  ,iStage-1,iPart) = PartState(5,iPart) * LorentzFacInv
          PartStage(3  ,iStage-1,iPart) = PartState(6,iPart) * LorentzFacInv
          PartStage(4:6,iStage-1,iPart) = Pt(1:3,iPart)
        END IF
      ELSE ! PartIsExplicit
        CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
        reMap=.FALSE.
        IF(PartMeshHasReflectiveBCs)THEN
          IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
            n_loc=PEM%NormVec(1:3,iPart)
            ! particle is actually located outside, hence, it moves in the mirror field
            ! mirror electric field, constant B field
            FieldAtParticle(1:3,iPart)=FieldAtParticle(1:3,iPart)-2.*DOT_PRODUCT(FieldAtParticle(1:3,iPart),n_loc)*n_loc
            FieldAtParticle(4:6,iPart)=FieldAtParticle(4:6,iPart)!-2.*DOT_PRODUCT(FieldAtParticle(4:6,iPart),n_loc)*n_loc
            PEM%NormVec(1:3,iPart)=0.
            ! and of coarse, the velocity has to be back-rotated, because the particle has not hit the wall
            reMap=.TRUE.
          END IF
        END IF
        IF(PEM%PeriodicMoved(iPart))THEN
          PEM%PeriodicMoved(iPart)=.FALSE.
          Remap=.TRUE.
        END IF
        IF(reMap)THEN
          ! recompute explicit push, requires only the velocity/momentum
          PartState(4:6,iPart) = ERK_a(iStage-1,iStage-2)*PartStage(4:6,iStage-2,iPart)
          DO iCounter=1,iStage-2
            PartState(4:6,iPart)=PartState(4:6,iPart)+ERK_a(iStage-1,iCounter)*PartStage(4:6,iCounter,iPart)
          END DO ! iCounter=1,iStage-2
          PartState(4:6,iPart)=PartStateN(4:6,iPart)+dt*PartState(4:6,iPart)
          ! luckily - nothing to do
        END IF
        IF(PartLorentzType.EQ.5)THEN
          LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
          CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
        ELSE
          LorentzFacInv = 1.0
          CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
        END IF ! PartLorentzType.EQ.5
        PartStage(1  ,iStage-1,iPart) = PartState(4,iPart)*LorentzFacInv
        PartStage(2  ,iStage-1,iPart) = PartState(5,iPart)*LorentzFacInv
        PartStage(3  ,iStage-1,iPart) = PartState(6,iPart)*LorentzFacInv
        PartStage(4:6,iStage-1,iPart) = Pt(1:3,iPart)
      END IF ! ParticleIsImplicit
    END DO ! iPart
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! new method
  END IF
#endif /*PARTICLES*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! explicit - particle  pusher
  !--------------------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_StdOut,'(A)') '-----------------------------'
    SWRITE(UNIT_StdOut,'(A)') ' explicit particles'
  END IF
  ExplicitPartSource=0.
  ! particle step || only explicit particles
  IF (time.GE.DelayTime) THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(.NOT.PartIsImplicit(iPart))THEN !explicit particles only
        ! particle movement from PartSateN
        LastPartPos(1,iPart)  =PartStateN(1,iPart)
        LastPartPos(2,iPart)  =PartStateN(2,iPart)
        LastPartPos(3,iPart)  =PartStateN(3,iPart)
        PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
        ! delete rotation || periodic movement
        IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart) = 0.
        PEM%PeriodicMoved(iPart) =.FALSE.
        ! compute explicit push
        PartState(1:6,iPart) = ERK_a(iStage,iStage-1)*PartStage(1:6,iStage-1,iPart)
        DO iCounter=1,iStage-2
          PartState(1:6,iPart)=PartState(1:6,iPart)+ERK_a(iStage,iCounter)*PartStage(1:6,iCounter,iPart)
        END DO ! iCounter=1,iStage-2
        PartState(1:6,iPart)=PartStateN(1:6,iPart)+dt*PartState(1:6,iPart)
      END IF ! ParticleIsExplicit
    END DO ! iPart
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
    ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
  ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
    SELECT CASE(TrackingMethod)
    CASE(REFMAPPING)
      CALL ParticleRefTracking(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    CASE(TRACING)
      CALL ParticleTracing(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    CASE(TRIATRACKING)
      CALL ParticleTriaTracking(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    CASE DEFAULT
      CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
    END SELECT
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
    ! send number of particles
    CALL SendNbOfParticles(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! if new
    ! map particle from gamma*v to velocity
    CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! deposit explicit, local particles
    CALL Deposition(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    IF(DoDeposition) THEN
      DO iElem=1,PP_nElems
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
          ExplicitPartSource(1:4,i,j,k,iElem)=PartSource(1:4,i,j,k,iElem)
        END DO; END DO; END DO !i,j,k
      END DO !iElem
    END IF
    ! map particle from v to gamma*v
    CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
#endif /*PARTICLES*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! implicit - particle pusher & Maxwell's field
  !--------------------------------------------------------------------------------------------------------------------------------

#if !(USE_HDG)
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! compute RHS for linear solver
  LinSolverRHS=ESDIRK_a(iStage,iStage-1)*FieldStage(:,:,:,:,:,iStage-1)
  DO iCounter=1,iStage-2
    LinSolverRHS=LinSolverRHS +ESDIRK_a(iStage,iCounter)*FieldStage(:,:,:,:,:,iCounter)
  END DO ! iCoutner=1,iStage-2
  LinSolverRHS=Un+dt*LinSolverRHS
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! get predictor of u^s+1 ! wrong position
  ! CALL Predictor(iStage,dt,Un,FieldStage) ! sets new value for U_DG
#endif /*DG*/

#ifdef PARTICLES
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_StdOut,'(A)') '-----------------------------'
    SWRITE(UNIT_StdOut,'(A)') ' implicit particles '
  END IF
  IF (time.GE.DelayTime) THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        PartIsImplicit(iPart)=.FALSE.
        PDM%IsNewPart(iPart) =.FALSE.
        CYCLE
      END IF
      IF(PartIsImplicit(iPart))THEN
        ! ignore surface flux particle
        ! dirty hack, if particle does not take part in implicit treating, it is removed from this list
        ! surface flux particles
        LastPartPos(1,iPart)=PartStateN(1,iPart)
        LastPartPos(2,iPart)=PartStateN(2,iPart)
        LastPartPos(3,iPart)=PartStateN(3,iPart)
        PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
        IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart) =0.
        PEM%PeriodicMoved(iPart) = .FALSE.
        ! compute Q and U
        PartQ(1:6,iPart) = ESDIRK_a(iStage,iStage-1)*PartStage(1:6,iStage-1,iPart)
        DO iCounter=1,iStage-2
          PartQ(1:6,iPart) = PartQ(1:6,iPart) + ESDIRK_a(iStage,iCounter)*PartStage(1:6,iCounter,iPart)
        END DO ! iCounter=1,iStage-2
        IF(DoSurfaceFlux)THEN
          Dtloc=PartDtFrac(iPart)*dt
        ELSE
          Dtloc=dt
        END IF
        PartQ(1:6,iPart) = PartStateN(1:6,iPart) + dtloc* PartQ(1:6,iPart)
        ! do not use a predictor
        ! position is already safed
        ! caution: has to use 4:6 because particle is not tracked
        !IF(PredictorType.GT.0)THEN
        PartState(1:6,iPart)=PartQ(1:6,iPart)
        !ELSE
        !  PartState(4:6,iPart)=PartQ(4:6,iPart)
        !END IF
        ! store the information before the Newton method || required for init
        PartDeltaX(1:6,iPart) = 0. ! PartState(1:6,iPart)-PartStateN(1:6,iPart)
        PartXk(1:6,iPart)     = PartState(1:6,iPart)
      END IF ! PartIsImplicit
    END DO ! iPart
#if USE_LOADBALANCE
    ! the measurement contains information from SurfaceFlux and implicit particle pushing
    ! simplicity: added to push time
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
    IF(PredictorType.GT.0)THEN
#if USE_LOADBALANCE
      CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
      ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
      ! open receive buffer for number of particles
      CALL IRecvNbofParticles()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
      CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
      SELECT CASE(TrackingMethod)
      CASE(REFMAPPING)
        CALL ParticleRefTracking(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      CASE(TRACING)
        CALL ParticleTracing(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      CASE(TRIATRACKING)
        CALL ParticleTriaTracking(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      CASE DEFAULT
        CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
      END SELECT
#if USE_MPI
#if USE_LOADBALANCE
      CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
      ! send number of particles
      CALL SendNbOfParticles(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      !CALL SendNbOfParticles() ! all particles to get initial deposition right \\ without emmission
      ! finish communication of number of particles and send particles
      CALL MPIParticleSend()
      ! finish communication
      CALL MPIParticleRecv()
#if USE_LOADBALANCE
      CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
      DO iPart=1,PDM%ParticleVecLength
        IF(PartIsImplicit(iPart))THEN
          IF(.NOT.PDM%ParticleInside(iPart)) PartIsImplicit(iPart)=.FALSE.
          IF(.NOT.PDM%ParticleInside(iPart)) PDM%IsNewPart(iPart)=.FALSE.
        END IF ! PartIsImplicit
      END DO ! iPart
#if USE_LOADBALANCE
      CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
    END IF

  END IF
#endif /*PARTICLES*/
  ! full newton for particles and fields
  CALL FullNewton(time,tStage,alpha)
#if !(USE_HDG)
  CALL DivCleaningDamping()
#endif /*DG*/
#ifdef PARTICLES
  IF (time.GE.DelayTime) THEN
    IF(DoUpdateInStage) CALL UpdateNextFreePosition()
  END IF

#ifdef CODE_ANALYZE
  !SWRITE(*,*) 'sanity check'
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
     CALL ParticleSanityCheck(iPart)
  END DO
#endif /*CODE_ANALYZE*/

#endif /*PARTICLES*/
END DO

#ifdef PARTICLES
IF (time.GE.DelayTime) THEN
  IF(DoSurfaceFlux)THEN
    DO iPart=1,PDM%ParticleVecLength
      IF(PDM%ParticleInside(iPart))THEN
        IF(PDM%IsNewPart(iPart)) THEN
          PartDtFrac(iPart)=1.
          PDM%IsNewPart(iPart)=.FALSE.
        END IF
      END IF ! ParticleInside
    END DO ! iPart
  END IF

! particle step || only explicit particles
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! add here new method
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart))CYCLE
    IF(DoSurfaceFlux)THEN
      ! sanity check
      IF(PDM%IsNewPart(iPart)) CALL abort(&
   __STAMP__&
   ,' Particle error with surfaceflux. Part-ID: ',iPart)
    END IF
    LastPartPos(1,iPart)=PartStateN(1,iPart)
    LastPartPos(2,iPart)=PartStateN(2,iPart)
    LastPartPos(3,iPart)=PartStateN(3,iPart)
    PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
    IF(.NOT.PartIsImplicit(iPart))THEN
      ! LastPartPos(1,iPart)=PartState(1,iPart)
      ! LastPartPos(2,iPart)=PartState(2,iPart)
      ! LastPartPos(3,iPart)=PartState(3,iPart)
      ! PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
      reMap=.FALSE.
      IF(PartMeshHasReflectiveBCs)THEN
        IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
          n_loc=PEM%NormVec(1:3,iPart)
          ! particle is actually located outside, hence, it moves in the mirror field
          ! mirror electric field, constant B field
          FieldAtParticle(1:3,iPart)=FieldAtParticle(1:3,iPart)-2.*DOT_PRODUCT(FieldAtParticle(1:3,iPart),n_loc)*n_loc
          FieldAtParticle(4:6,iPart)=FieldAtParticle(4:6,iPart)!-2.*DOT_PRODUCT(FieldAtParticle(4:6,iPart),n_loc)*n_loc
          PEM%NormVec(1:3,iPart)=0.
          ! and of coarse, the velocity has to be back-rotated, because the particle has not hit the wall
          reMap=.TRUE.
        END IF
      END IF
      IF(PEM%PeriodicMoved(iPart))THEN
        PEM%PeriodicMoved(iPart)=.FALSE.
        Remap=.TRUE.
      END IF
      IF(reMap)THEN
        ! recompute explicit push, requires only the velocity/momentum
        iStage=nRKStages
        PartState(4:6,iPart) = ERK_a(iStage,iStage-1)*PartStage(4:6,iStage-1,iPart)
        DO iCounter=1,iStage-2
          PartState(4:6,iPart)=PartState(4:6,iPart)+ERK_a(iStage,iCounter)*PartStage(4:6,iCounter,iPart)
        END DO ! iCounter=1,iStage-2
        PartState(4:6,iPart)=PartStateN(4:6,iPart)+dt*PartState(4:6,iPart)
      END IF
      ! compute acceleration
      IF(PartLorentzType.EQ.5)THEN
        LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
      ELSE
        LorentzFacInv = 1.0
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
      END IF ! PartLorentzType.EQ.5
      PartState(1  ,iPart) = RK_b(nRKStages)*LorentzFacInv*PartState(4,iPart)
      PartState(2  ,iPart) = RK_b(nRKStages)*LorentzFacInv*PartState(5,iPart)
      PartState(3  ,iPart) = RK_b(nRKStages)*LorentzFacInv*PartState(6,iPart)
      PartState(4:6,iPart) = RK_b(nRKSTages)*Pt(1:3,iPart)
      !  stage 1 ,nRKStages-1
      DO iCounter=1,nRKStages-1
        PartState(1:6,iPart) = PartState(1:6,iPart)   &
                             + RK_b(iCounter)*PartStage(1:6,iCounter,iPart)
      END DO ! counter
      PartState(1:6,iPart) = PartStateN(1:6,iPart)+dt*PartState(1:6,iPart)

    END IF ! ParticleIsExplicit
  END DO ! iPart
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

iStage=0
#if USE_MPI
  ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
  CALL PerformTracking()
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
  ! send number of particles
  CALL SendNbOfParticles() !doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
  !CALL SendNbOfParticles() ! all particles to get initial deposition right \\ without emmission
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
! #endif ! old -> new is 9 lines below
#endif
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! map particle from gamma*v to v
  CALL PartVeloToImp(VeloToImp=.FALSE.)
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF
#endif /*PARTICLES*/

!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
IF (useDSMC) THEN ! UpdateNextFreePosition is only required for DSMC
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
END IF

IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL UpdateNextFreePosition()

#endif /*PARTICLES*/

END SUBROUTINE TimeStepByImplicitRK


END MODULE MOD_TimeStep
#endif /*IMPA*/
