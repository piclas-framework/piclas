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


#if ROS
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepByRosenbrock
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStepByRosenbrock()
!===================================================================================================================================
! Rosenbrock-Method
! from: Kaps. 1979
! linear,implicit Runge-Kutta method
! RK-method is 1 Newton-Step
! optimized implementation following E. Hairer in "Solving Ordinary Differential Equations 2" p.111FF
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iter,iStage, nRKStages,dt_inv,dt_old, time
USE MOD_TimeDisc_Vars          ,ONLY: RK_a,RK_c,RK_g,RK_b,RK_gamma
USE MOD_LinearSolver_Vars      ,ONLY: FieldStage,DoPrintConvInfo
USE MOD_DG_Vars                ,ONLY: U,Un
#if USE_HDG
USE MOD_HDG                    ,ONLY: HDG
#else /*pure DG*/
USE MOD_Precond_Vars           ,ONLY: UpdatePrecond
USE MOD_LinearOperator         ,ONLY: MatrixVector
USE MOD_LinearSolver           ,ONLY: LinearSolver
USE MOD_DG_Vars                ,ONLY: Ut
USE MOD_DG                     ,ONLY: DGTimeDerivative_weakForm
USE MOD_LinearSolver_Vars      ,ONLY: LinSolverRHS
USE MOD_Equation               ,ONLY: DivCleaningDamping
USE MOD_Equation               ,ONLY: CalcSource
#ifdef maxwell
USE MOD_Precond                ,ONLY: BuildPrecond
#endif /*maxwell*/
#endif /*USE_HDG*/
#ifdef PARTICLES
USE MOD_Globals_Vars           ,ONLY: c2_inv
USE MOD_LinearOperator         ,ONLY: PartMatrixVector, PartVectorDotProduct
USE MOD_ParticleSolver         ,ONLY: Particle_GMRES
USE MOD_LinearSolver_Vars      ,ONLY: PartXK,R_PartXK,DoFieldUpdate
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Particle_Vars          ,ONLY: PartLorentzType,PartDtFrac,PartStateN,PartStage,PartQ &
    ,DoSurfaceFlux,PEM,PDM,Pt,LastPartPos,DelayTime,PartState,PartMeshHasReflectiveBCs
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_part_RHS               ,ONLY: PartVeloToImp
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_Part_RHS               ,ONLY: PartRHS
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange
#endif /*USE_MPI*/
USE MOD_PIC_Analyze            ,ONLY: CalcDepositedCharge
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
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
REAL               :: coeff,coeff_inv
!INTEGER            :: i,j
REAL               :: tRatio
! particle surface flux
! RK counter
INTEGER            :: iCounter
#ifdef PARTICLES
REAL               :: coeff_loc,dt_inv_loc, LorentzFacInv
REAL               :: PartDeltaX(1:6), PartRHS_loc(1:6), Norm_P2, Pt_tmp(1:6), PartRHS_tild(1:6),FieldAtParticle_loc(1:6)
REAL               :: RandVal!, LorentzFac
REAL               :: AbortCrit
INTEGER            :: iPart,nParts
REAL               :: n_loc(1:3)
LOGICAL            :: reMap
#endif /*PARTICLES*/
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
#ifdef maxwell
LOGICAL            :: UpdatePrecondLoc
#endif /*maxwell*/
!===================================================================================================================================

coeff=dt*RK_gamma
coeff_inv=1./coeff
dt_inv=1.
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
dt_inv=dt_inv/dt
tRatio = 1.

! ! ! sanity check
! ! print*,'RK_gamma',RK_gamma
! DO istage=2,nRKStages
!   DO iCounter=1,nRKStages
!     print*,'a',iStage,iCounter,RK_a(iStage,iCounter)
!   END DO
! END DO
! DO istage=2,nRKStages
!   DO iCounter=1,nRKStages
!     print*,'c',iStage,iCounter,RK_g(iStage,iCounter)
!   END DO
! END DO
! DO iCounter=1,nRKStages
!   print*,'b',iStage,iCounter,RK_b(iCounter)
! END DO
! time level of stage
! DO istage=2,nRKStages
!  print*,'aa',SUM(RK_A(iStage,:))*RK_gamma,RK_c
! END DO

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
CALL LBPauseTime(LB_EMISSION,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization && first linear solver
! ----------------------------------------------------------------------------------------------------------------------------------
tStage=time

#ifdef PARTICLES
! compute number of emitted particles during Rosenbrock-Step
IF(time.GE.DelayTime)THEN
  ! surface flux
  ! major difference to all other routines:
  ! Rosenbrock is one Newton-Step, hence, all particles cross the surface at time^n, however, the time-step size
  ! different for all particles, hence, some move slower, same move faster. this is not 100% correct, however,
  ! typical, RK_c(1) is zero.... and we compute the first stage..
  IF(DoSurfaceFlux)THEN
    IF(DoPrintConvInfo) nParts=0
    ! LB time measurement is performed within ParticleSurfaceflux
    CALL ParticleSurfaceflux()
    ! compute emission for all particles during dt
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    DO iPart=1,PDM%ParticleVecLength
      ! only particle within
      IF(PDM%ParticleInside(iPart))THEN
        ! copy elemnent N
        ! check if new particle
        IF(PDM%IsNewPart(iPart))THEN
          ! initialize of surface-flux particles
          IF(DoPrintConvInfo) nParts=nParts+1
          ! gives entry point into domain
          CALL RANDOM_NUMBER(RandVal)
          PartDtFrac(iPart)=RandVal
          ! particle crosses surface at time^n + (1.-RandVal)*dt
          ! for all stages   t_Stage =< time^n + (1.-RandVal)*dt particle is outside of domain
          ! for              t_Stage >  time^n + (1.-RandVal)*dt particle is in domain and can be advanced in time
          ! particle is no-longer a SF particle
          PDM%IsNewPart(iPart) = .FALSE.
        ELSE
          ! set DtFrac to unity
          PartDtFrac(iPart)=1.0
        END IF ! IsNewPart
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
END IF

! simulation with delay-time, compute the
IF(DelayTime.GT.0.)THEN
  IF((iter.EQ.0).AND.(time.LT.DelayTime))THEN
    ! perform normal deposition
    CALL Deposition()
  END IF
END IF

! compute source of first stage for Maxwell solver
IF (time.GE.DelayTime) THEN
  CALL Deposition()
END IF

#if USE_HDG
! update the fields due to changed particle number: emission or velocity change in DSMC
! LB measurement is performed within HDG
IF(DoFieldUpdate) CALL HDG(time,U,iter)
#endif

IF(time.GE.DelayTime)THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! map velocity space to relativistic momentum
  iStage=1
  CALL PartVeloToImp(VeloToImp=.TRUE.)
  PartStateN(1:PDM%ParticleVecLength,1:6)=PartState(1:6,1:PDM%ParticleVecLength)
  PEM%GlobalElemID(1:PDM%ParticleVecLength)  =PEM%GlobalElemID(1:PDM%ParticleVecLength)
  ! should be already be done
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart))CYCLE
    ! store old particle position
    LastPartPos(1,iPart)  =PartStateN(1,iPart)
    LastPartPos(2,iPart)  =PartStateN(2,iPart)
    LastPartPos(3,iPart)  =PartStateN(3,iPart)
    ! copy date
    PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
    IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart)=0.
    PEM%PeriodicMoved(iPart) = .FALSE.
    ! build RHS of particle with current DG solution and particle position
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
    ! compute particle RHS at time^n
    IF(PartLorentzType.EQ.5)THEN
      LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
      CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
    ELSE
      LorentzFacInv = 1.0
      CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
    END IF ! PartLorentzType.EQ.5
    ! compute current Pt_tmp for the particle
    Pt_tmp(1) =LorentzFacInv*PartState(4,iPart)
    Pt_tmp(2) =LorentzFacInv*PartState(5,iPart)
    Pt_tmp(3) =LorentzFacInv*PartState(6,iPart)
    Pt_tmp(4) =Pt(1,iPart)
    Pt_tmp(5) =Pt(2,iPart)
    Pt_tmp(6) =Pt(3,iPart)
    ! how it works
    ! A x = b
    ! A(x-x0) = b - A x0
    ! set x0 = b
    ! A deltaX = b - A b
    ! xNeu = b  + deltaX
    ! fix matrix during iteration
    PartXK(1:6,iPart)   = PartState(1:6,iPart) ! which is partstateN
    R_PartXK(1:6,iPart) = Pt_tmp(1:6)          ! the delta is not changed
    PartDeltaX=0.
    ! guess for new value is Pt_tmp: remap to reuse old GMRES
    ! OLD
    ! CALL PartMatrixVector(time,Coeff_inv,iPart,Pt_tmp,PartDeltaX)
    ! PartRHS_loc =Pt_tmp - PartDeltaX
    ! CALL PartVectorDotProduct(PartRHS_loc,PartRHS_loc,Norm_P2)
    ! AbortCrit=1e-16
    ! PartDeltaX=0.
    ! CALL Particle_GMRES(time,coeff_inv,iPart,PartRHS_loc,SQRT(Norm_P2),AbortCrit,PartDeltaX)
    ! NEW version is more stable, hence we use it!
    IF(DoSurfaceFlux)THEN
      coeff_loc=PartDtFrac(iPart)*coeff
    ELSE
      coeff_loc=coeff
    END IF
    Pt_tmp=coeff_loc*Pt_Tmp
    CALL PartMatrixVector(time,Coeff_loc,iPart,Pt_tmp,PartDeltaX)
    PartRHS_loc =Pt_tmp - PartDeltaX
    CALL PartVectorDotProduct(PartRHS_loc,PartRHS_loc,Norm_P2)
    AbortCrit=1e-16
    PartDeltaX=0.
    CALL Particle_GMRES(time,coeff_loc,iPart,PartRHS_loc,SQRT(Norm_P2),AbortCrit,PartDeltaX)
    ! update particle
    PartState(1:6,iPart)=Pt_tmp+PartDeltaX(1:6)
    PartStage(1,1,iPart) = PartState(1,iPart)
    PartStage(2,1,iPart) = PartState(2,iPart)
    PartStage(3,1,iPart) = PartState(3,iPart)
    PartStage(4,1,iPart) = PartState(4,iPart)
    PartStage(5,1,iPart) = PartState(5,iPart)
    PartStage(6,1,iPart) = PartState(6,iPart)
  END DO ! iPart
  ! track particle
  iStage=1
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF ! time.GE. DelayTime
IF(DoFieldUpdate)THEN
#endif /*PARTICLES*/

#if !(USE_HDG)
! LB measurement is performed within DGTimeDerivative_weakForm and LinearSolver (again DGTimeDerivative_weakForm)
! the copy time of the arrays is ignored within this measurement
Un = U
! solve linear system for electromagnetic field
! RHS is f(u^n+0) = DG_u^n + source^n
CALL DGTimeDerivative_weakForm(time, time, 0,doSource=.TRUE.) ! source terms are added NOT added in linearsolver
LinSolverRHS =Ut
CALL LinearSolver(tStage,coeff_inv)
FieldStage (:,:,:,:,:,1) = U
#endif /*DG*/
#ifdef PARTICLES
END IF ! DoFieldUpdate
#endif /*PARTICLES*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to nRKStages
! ----------------------------------------------------------------------------------------------------------------------------------
DO iStage=2,nRKStages
  ! time of current stage
  tStage = time + RK_c(iStage)*dt
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_StdOut,'(A)')    '-----------------------------'
    SWRITE(UNIT_StdOut,'(A,I2)') 'istage:',istage
  END IF

  !--------------------------------------------------------------------------------------------------------------------------------
  ! DGSolver: explicit contribution and 1/dt_inv sum_ij RK_g FieldStage
  ! is the state before the linear system is solved
  !--------------------------------------------------------------------------------------------------------------------------------
#if !(USE_HDG)
#ifdef PARTICLES
  IF(DoFieldUpdate) THEN
#endif /*PARTICLES*/
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    ! compute contribution of h Time * sum_j=1^iStage-1
    ! In optimized version sum_j=1^&iStage-1 c_ij/dt*Y_J
    LinSolverRHS = RK_g(iStage,iStage-1)*FieldStage(:,:,:,:,:,iStage-1)
    DO iCounter=1,iStage-2
      LinSolverRHS = LinSolverRHS+RK_g(iStage,iCounter)*FieldStage(:,:,:,:,:,iCounter)
    END DO ! iCounter=1,iStage-2
    ! not required in the optimized implementation, instead
    ! !! matrix vector WITH dt and contribution of Time * sum
    ! !! Jacobian times sum_j^iStage-1 gamma_ij k_j
    ! !U=LinSolverRHS
    ! !CALL DGTimeDerivative_weakForm(time, time, 0,doSource=.FALSE.)
    ! !LinSolverRHS=Ut*dt
    ! OPTIMIZED IMPLEMENTATION
    LinSolverRHS=dt_inv*LinSolverRHS
    ! compute explicit contribution  AGAIN no dt
    U = RK_a(iStage,iStage-1)*FieldStage(:,:,:,:,:,iStage-1)
    DO iCounter=1,iStage-2
      U=U+RK_a(iStage,iCounter)*FieldStage(:,:,:,:,:,iCounter)
    END DO ! iCounter=1,iStage-2
    U=Un+U
    CALL DivCleaningDamping()
    ! this U is required for the particles and the interpolation, hence, we have to continue with the particles
    ! which gives us the source terms, too...
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
  END IF
#endif /*PARTICLES*/
#endif /*NOT HDG->DG*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! particle  pusher: explicit contribution and Time * sum
  ! is the state before the linear system is solved
  !--------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
  IF (time.GE.DelayTime) THEN
    IF(DoPrintConvInfo)THEN
      SWRITE(UNIT_StdOut,'(A)') '-----------------------------'
      SWRITE(UNIT_StdOut,'(A)') ' Implicit step.   '
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
    ! move particle to PartState^n + dt*sum_j=1^(iStage-1)*aij k_j
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
      ! compute contribution of 1/dt* sum_j=1^iStage-1 c(i,j) = diag(gamma)-gamma^inv
      PartQ(1:6,iPart) = RK_g(iStage,iStage-1)*PartStage(1:6,iStage-1,iPart)
      DO iCounter=1,iStage-2
        PartQ(1:6,iPart) = PartQ(1:6,iPart) +RK_g(iStage,iCounter)*PartStage(1:6,iCounter,iPart)
      END DO ! iCounter=1,iStage-2
      IF(DoSurfaceFlux)THEN
        dt_inv_loc=dt_inv/PartDtFrac(iPart)
      ELSE
        dt_inv_loc=dt_inv
      END IF
      PartQ(1:6,iPart) = dt_inv_loc*PartQ(1:6,iPart)
      ! compute explicit contribution which is
      PartState(1:6,iPart) = RK_a(iStage,iStage-1)*PartStage(1:6,iStage-1,iPart)
      DO iCounter=1,iStage-2
        PartState(1:6,iPart)=PartState(1:6,iPart)+RK_a(iStage,iCounter)*PartStage(1:6,iCounter,iPart)
      END DO ! iCounter=1,iStage-2
      PartState(1:6,iPart)=PartStateN(1:6,iPart)+PartState(1:6,iPart)
    END DO ! iPart=1,PDM%ParticleVecLength
    CALL PartVeloToImp(VeloToImp=.FALSE.)
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
    CALL PerformTracking()
#if USE_MPI
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! compute particle source terms on field solver of implicit particles :)
    CALL Deposition()
    ! map particle from v to gamma v
#if USE_HDG
    ! update the fields due to changed particle position and velocity/momentum
    ! LB-TimeMeasurement is performed within HDG
    IF(DoFieldUpdate) CALL HDG(time,U,iter)
#endif
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    CALL PartVeloToImp(VeloToImp=.TRUE.)
    ! should be already be done
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      ! position currently updated, hence store the data for tracking
      LastPartPos(1,iPart)  =PartStateN(1,iPart)
      LastPartPos(2,iPart)  =PartStateN(2,iPart)
      LastPartPos(3,iPart)  =PartStateN(3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      ! build RHS of particle with current DG solution and particle position
      ! CAUTION: we have to use a local variable here. The Jacobian matrix is FIXED for one time step,
      ! hence the fields for the matrix-vector-multiplication shall NOT be updated
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle_loc(1:6))
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
        ! and of coarse, the velocity has to be back-rotated, because the particle has not hit the wall
        ! it is more stable, to recompute the position
        ! compute explicit contribution which is
        PartState(1:6,iPart) = RK_a(iStage,iStage-1)*PartStage(1:6,iStage-1,iPart)
        DO iCounter=1,iStage-2
          PartState(1:6,iPart)=PartState(1:6,iPart)+RK_a(iStage,iCounter)*PartStage(1:6,iCounter,iPart)
        END DO ! iCounter=1,iStage-2
        PartState(1:6,iPart)=PartStateN(1:6,iPart)+PartState(1:6,iPart)
      END IF ! PartMeshHasReflectiveBCs
      ! compute particle RHS at time^n
      IF(PartLorentzType.EQ.5)THEN
        LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
        CALL PartRHS(iPart,FieldAtParticle_loc(1:6),Pt(1:3,iPart),LorentzFacInv)
      ELSE
        LorentzFacInv = 1.0
        CALL PartRHS(iPart,FieldAtParticle_loc(1:6),Pt(1:3,iPart))
      END IF ! PartLorentzType.EQ.5
      ! compute current Pt_tmp for the particle
      Pt_tmp(1) = LorentzFacInv*PartState(4,iPart)
      Pt_tmp(2) = LorentzFacInv*PartState(5,iPart)
      Pt_tmp(3) = LorentzFacInv*PartState(6,iPart)
      Pt_tmp(4) = Pt(1,iPart)
      Pt_tmp(5) = Pt(2,iPart)
      Pt_tmp(6) = Pt(3,iPart)
      IF(DoSurfaceFlux)THEN
        coeff_loc=PartDtFrac(iPart)*coeff
      ELSE
        coeff_loc=coeff
      END IF
      PartRHS_loc =(Pt_tmp + PartQ(1:6,iPart))*coeff_loc
      ! guess for new particleposition is PartState || reuse of OLD GMRES
      CALL PartMatrixVector(time,Coeff_loc,iPart,PartRHS_loc,PartDeltaX)
      PartRHS_tild = PartRHS_loc - PartDeltaX
      PartDeltaX=0.
      CALL PartVectorDotProduct(PartRHS_tild,PartRHS_tild,Norm_P2)
      AbortCrit=1e-16
      CALL Particle_GMRES(time,coeff_loc,iPart,PartRHS_tild,SQRT(Norm_P2),AbortCrit,PartDeltaX)
      ! update particle to k_iStage
      PartState(1:6,iPart)=PartRHS_loc+PartDeltaX(1:6)
      !PartState(1:6,iPart)=PartRHS_loc+PartDeltaX(1:6)
      ! and store value as k_iStage
      IF(iStage.LT.nRKStages)THEN
        PartStage(1,iStage,iPart) = PartState(1,iPart)
        PartStage(2,iStage,iPart) = PartState(2,iPart)
        PartStage(3,iStage,iPart) = PartState(3,iPart)
        PartStage(4,iStage,iPart) = PartState(4,iPart)
        PartStage(5,iStage,iPart) = PartState(5,iPart)
        PartStage(6,iStage,iPart) = PartState(6,iPart)
      END IF
    END DO ! iPart
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
  IF(DoFieldUpdate) THEN
#endif /*PARTICLES*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! DGSolver: now, we can add the contribution of the particles
  ! LB-TimeMeasurement is performed in DGTimeDerivative_weakForm and LinearSolver (DGTimeDerivative_weakForm), hence,
  ! it is neglected here, array copy assumed to be zero
  !--------------------------------------------------------------------------------------------------------------------------------
#if !(USE_HDG)
    ! next DG call is f(u^n + dt sum_j^i-1 a_ij k_j) + source terms
    CALL DGTimeDerivative_weakForm(time, time, 0,doSource=.TRUE.) ! source terms are not-added in linear solver
    ! CAUTION: invert sign of Ut
    LinSolverRHS = LinSolverRHS + Ut
    CALL LinearSolver(tStage,coeff_inv)
    ! and store U in fieldstage
    IF(iStage.LT.nRKStages) FieldStage (:,:,:,:,:,iStage) = U
#endif /*NOT HDG->DG*/
#ifdef PARTICLES
  END IF ! DoFieldUpdate
#endif /*PARTICLES*/
END DO


#if !(USE_HDG)
#ifdef PARTICLES
IF(DoFieldUpdate)THEN
#endif /*PARTICLES*/
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! update field step
  U = RK_b(nRKStages)* U
  DO iCounter=1,nRKStages-1
    U = U +  RK_b(iCounter)* FieldStage(:,:,:,:,:,iCounter)
  END DO ! counter
  U = Un +  U
  CALL DivCleaningDamping()
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
END IF ! DoFieldUpdate
#endif /*PARTICLES*/
#endif /*NOT HDG->DG*/

#ifdef PARTICLES
! particle step || only explicit particles
IF (time.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! add here new method
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart))CYCLE
    !  stage 1 ,nRKStages-1
    PartState(1:6,iPart) = RK_b(nRKStages)*PartState(1:6,iPart)
    DO iCounter=1,nRKStages-1
      PartState(1:6,iPart) = PartState(1:6,iPart) + RK_b(iCounter)*PartStage(1:6,iCounter,iPart)
    END DO ! counter
    PartState(1:6,iPart) = PartStateN(1:6,iPart)+PartState(1:6,iPart)
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
  CALL SendNbOfParticles() ! all particles to get initial deposition right \\ without emmission
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
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
IF (useDSMC) THEN
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
END IF

IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL UpdateNextFreePosition()
#endif /*PARTICLES*/

END SUBROUTINE TimeStepByRosenbrock

END MODULE MOD_TimeStep
#endif /*ROS  */
