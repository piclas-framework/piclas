#include "boltzplatz.h"


MODULE MOD_TimeDisc
!===================================================================================================================================
! Module for the GTS Temporal discretization  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitTimeDisc
  MODULE PROCEDURE InitTimeDisc
END INTERFACE

INTERFACE TimeDisc
  MODULE PROCEDURE TimeDisc
END INTERFACE

INTERFACE FinalizeTimeDisc
  MODULE PROCEDURE FinalizeTimeDisc
END INTERFACE

PUBLIC :: InitTimeDisc,FinalizeTimeDisc
PUBLIC :: TimeDisc
!===================================================================================================================================

CONTAINS

SUBROUTINE InitTimeDisc()
!===================================================================================================================================
! Get information for end time and max time steps from ini file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,          ONLY:GetReal,GetInt, GETLOGICAL
USE MOD_TimeDisc_Vars,        ONLY:CFLScale,dt,TimeDiscInitIsDone
#if (PP_TimeDiscMethod>=100 && PP_TimeDiscMethod<200) 
USE MOD_TimeDisc_Vars,        ONLY:epsTilde_LinearSolver,eps_LinearSolver,eps2_LinearSolver,maxIter_LinearSolver
#endif /*PP_TimeDiscMethod>=100*/
USE MOD_TimeDisc_Vars,        ONLY:IterDisplayStep,DoDisplayIter,IterDisplayStepUser,DoDisplayEmissionWarnings
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(TimeDiscInitIsDone)THEN
   SWRITE(*,*) "InitTimeDisc already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TIMEDISC...'

! Read the normalized CFL number
CFLScale = GETREAL('CFLScale')
CALL fillCFL_DFL()

! read in requested IterDisplayStep (i.e. how often the message "iter: etc" is displayed, might be changed dependent on Particle-dt)
DoDisplayIter=.FALSE.
IterDisplayStepUser = GETINT('IterDisplayStep','1')
DoDisplayEmissionWarnings= GETLOGICAL('DoDisplayEmissionWarning','T')
IterDisplayStep = IterDisplayStepUser
IF(IterDisplayStep.GE.1) DoDisplayIter=.TRUE.


#if (PP_TimeDiscMethod>=100 && PP_TimeDiscMethod<200) 
eps_LinearSolver = GETREAL('eps_LinearSolver')
epsTilde_LinearSolver = eps_LinearSolver
eps2_LinearSolver = eps_LinearSolver *eps_LinearSolver 
maxIter_LinearSolver = GETINT('maxIter_LinearSolver')
#endif

! Set timestep to a large number
dt=HUGE(1.)
#if (PP_TimeDiscMethod==1)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK3-3 '
#elif (PP_TimeDiscMethod==2)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-5 '
#elif (PP_TimeDiscMethod==3)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: TAYLOR'
#elif (PP_TimeDiscMethod==4)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: DSMC-Only'
#elif (PP_TimeDiscMethod==5)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: RK4-Field, Euler-Explicit-Particles'
#elif (PP_TimeDiscMethod==6)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-14 '
#elif (PP_TimeDiscMethod==42)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: DSMC Reservoir and Debug'
#elif (PP_TimeDiscMethod==100)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Implicit'
#elif (PP_TimeDiscMethod==101)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ESDIRK4 Implicit'
#elif (PP_TimeDiscMethod==200)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Static Explicit'
#elif (PP_TimeDiscMethod==201)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Static Explicit with adaptive TimeStep'
#elif (PP_TimeDiscMethod==1000)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LD-Only'
#elif (PP_TimeDiscMethod==1001)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LD-DSMC'
# endif
TimediscInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTimeDisc



SUBROUTINE TimeDisc()
!===================================================================================================================================
! GTS Temporal discretization 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_AnalyzeField,          ONLY: CalcPoyntingIntegral
USE MOD_TimeDisc_Vars,         ONLY: TEnd,dt,tAnalyze,iter,IterDisplayStep,DoDisplayIter,IterDisplayStepUser
USE MOD_Restart_Vars,          ONLY: DoRestart,RestartTime
USE MOD_CalcTimeStep,          ONLY: CalcTimeStep
USE MOD_Analyze,               ONLY: CalcError,PerformAnalyze
USE MOD_Analyze_Vars,          ONLY: Analyze_dt,CalcPoyntingInt
#ifdef PARTICLES
USE MOD_Particle_Analyze,      ONLY: AnalyzeParticles
#else
USE MOD_AnalyzeField,          ONLY: AnalyzeField
#endif /*PARTICLES*/
USE MOD_Output,                ONLY: Visualize
USE MOD_HDF5_output,           ONLY: WriteStateToHDF5
USE MOD_Mesh_Vars,             ONLY: MeshFile,nGlobalElems
USE MOD_PML,                   ONLY: TransformPMLVars,BacktransformPMLVars
USE MOD_PML_Vars,              ONLY: DoPML
USE MOD_Filter,                ONLY: Filter
USE MOD_RecordPoints_Vars,     ONLY: RP_onProc
USE MOD_RecordPoints,          ONLY: RecordPoints,WriteRPToHDF5
#ifdef PARTICLES
USE MOD_PICDepo,               ONLY: Deposition!, DepositionMPF
USE MOD_PICDepo_Vars,          ONLY: DepositionType
USE MOD_Particle_Output,       ONLY: Visualize_Particles
USE MOD_PARTICLE_Vars,         ONLY: ManualTimeStep, Time, useManualTimestep, WriteMacroValues, MacroValSampTime
USE MOD_Particle_Tracking_vars, ONLY: tTracking,tLocalization,nTracks,MeasureTrackTime
#if (PP_TimeDiscMethod==201||PP_TimeDiscMethod==200)
USE MOD_PARTICLE_Vars,         ONLY: dt_maxwell,dt_max_particles,GEO,MaxwellIterNum,NextTimeStepAdjustmentIter
USE MOD_Equation_Vars,         ONLY: c
#endif /*(PP_TimeDiscMethod==201||PP_TimeDiscMethod==200)*/
#if (PP_TimeDiscMethod==201)
USE MOD_PARTICLE_Vars,         ONLY: PDM,Pt,PartState
#endif /*(PP_TimeDiscMethod==201)*/
USE MOD_PARTICLE_Vars,         ONLY : doParticleMerge, enableParticleMerge, vMPFMergeParticleIter
USE MOD_ReadInTools
USE MOD_DSMC_Vars,             ONLY: realtime,nOutput, Iter_macvalout
#ifdef MPI
USE MOD_Particle_MPI,          ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,     ONLY: PartMPIExchange
#endif /*MPI*/
#endif /*PARTICLES*/
#ifdef PP_POIS
USE MOD_Equation,              ONLY: EvalGradient
#endif /*PP_POIS*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: t,tFuture,tZero
INTEGER(KIND=8)              :: nAnalyze
REAL                         :: dt_Min, tEndDiff, tAnalyzeDiff
#if (PP_TimeDiscMethod==201)
REAL                         :: dt_temp,vMax,vMaxx,vMaxy,vMaxz
#endif
INTEGER(KIND=8)              :: iter_loc
REAL                         :: CalcTimeStart,CalcTimeEnd
INTEGER                      :: TimeArray(8)              ! Array for system time
#if (PP_TimeDiscMethod==201)
INTEGER                      :: MaximumIterNum, iPart
#endif
LOGICAL                      :: finalIter
!===================================================================================================================================
! init
SWRITE(UNIT_StdOut,'(132("-"))')
#ifdef PARTICLES
iter_macvalout=0
nOutput = 1
IF(.NOT.DoRestart)THEN
  t=0.
  realtime=0.
  time=0.
ELSE
  t=RestartTime
  realtime=RestartTime
  Time = RestartTime
END IF
IF (WriteMacroValues) MacroValSampTime = Time
#else
IF(.NOT.DoRestart)THEN
  t=0.
  SWRITE(UNIT_StdOut,*)'INITIAL PROJECTION:'
ELSE
  t=RestartTime
  SWRITE(UNIT_StdOut,*)'REWRITING SOLUTION:'
END IF
#endif /*PARTICLES*/
tZero=t
nAnalyze=1
tAnalyze=MIN(tZero+Analyze_dt,tEnd)

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems*(PP_N+1)**3)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)

! Determine next analyze time, since it will be written into output file
tFuture=MIN(t+Analyze_dt,tEnd)
!Evaluate Gradients to get Potential in case of Restart and Poisson Calc
#ifdef PP_POIS
IF(DoRestart) CALL EvalGradient()
#endif /*PP_POIS*/
! Write the state at time=0, i.e. the initial condition
CALL WriteStateToHDF5(TRIM(MeshFile),t,tFuture)

! if measurement of particle tracking time
#ifdef PARTICLES
IF(MeasureTrackTime)THEN
  nTracks=0
  tTracking=0
  tLocalization=0
END IF
#endif /*PARTICLES*/

! No computation needed if tEnd=tStart!
IF(t.EQ.tEnd)RETURN

#ifdef PARTICLES
#ifdef MPI
IF (DepositionType.EQ."shape_function") THEN
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
END IF
#endif /*MPI*/
#endif /*PARTICLES*/

!#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
!  CALL IRecvNbofParticles()
!  CALL MPIParticleSend()
!#endif /*MPI*/
!  CALL Deposition(doInnerParts=.TRUE.)
!#ifdef MPI
!  CALL MPIParticleRecv()
!  ! second buffer
!  CALL Deposition(doInnerParts=.FALSE.)
!#endif /*MPI*/
!#endif
#ifdef PARTICLES
! For tEnd != tStart we have to advance the solution in time
IF(useManualTimeStep)THEN
  ! particle time step is given externally and not calculated through the solver
  dt_Min=ManualTimeStep
#if (PP_TimeDiscMethod==200)
  dt_max_particles = ManualTimeStep
  dt_maxwell = CALCTIMESTEP()
  NextTimeStepAdjustmentIter = 0
  MaxwellIterNum = INT(MAX(GEO%xmaxglob-GEO%xminglob,GEO%ymaxglob-GEO%yminglob,GEO%zmaxglob-GEO%zminglob) &
                 / (c * dt_maxwell))
  IF (MaxwellIterNum*dt_maxwell.GT.dt_max_particles) THEN
    WRITE(*,*) 'WARNING: Time of Maxwell Solver is greater then particle time step!'
  END IF
  IterDisplayStep = MAX(INT(IterDisplayStepUser/MaxwellIterNum),1) !IterDisplayStep refers to dt_maxwell
  IF(MPIroot)THEN
    print*, 'IterNum for MaxwellSolver: ', MaxwellIterNum
    print*, 'Particle TimeStep: ', dt_max_particles  
    print*, 'Maxwell TimeStep: ', dt_maxwell
  END IF
#endif
ELSE ! .NO. ManualTimeStep
#endif /*PARTICLES*/
  ! time step is calculated by the solver
  ! first Maxwell time step for explicit LSRK
  dt_Min=CALCTIMESTEP()
  dt=dt_Min
  ! calculate time step for sub-cycling of divergence correction
  ! automatic particle time step of quasi-stationary time integration is not implemented
#ifdef PARTICLES
#if (PP_TimeDiscMethod==200)
  ! this will not work if particles have velocity of zero
  SWRITE(UNIT_StdOut, '(A)')'ERROR: with Static computations, a maximum delta t (=ManualTimeStep) needs to be given'
  STOP
#endif
END IF ! useManualTimestep

#if (PP_TimeDiscMethod==201)
dt_maxwell = CALCTIMESTEP()
MaximumIterNum = INT(MAX(GEO%xmaxglob-GEO%xminglob,GEO%ymaxglob-GEO%yminglob,GEO%zmaxglob-GEO%zminglob) &
               / (c * dt_maxwell))
IF(MPIroot)THEN
  print*, 'MaxIterNum for MaxwellSolver: ', MaximumIterNum
  print*, 'Maxwell TimeStep: ', dt_maxwell
END IF
dt_temp = 1E-8
#endif
#endif /*PARTICLES*/

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A,ES16.7)')'Initial Timestep  : ', dt_Min
! using sub-cycling requires an addional time step
SWRITE(UNIT_StdOut,*)'CALCULATION RUNNING...'
CalcTimeStart=BOLTZPLATZTIME()
iter=0
iter_loc=0

! fill recordpoints buffer (first iteration)
!IF(RP_onProc) CALL RecordPoints(iter,t,forceSampling=.TRUE.) 

! fill initial analyze stuff
CALL PerformAnalyze(t,iter,0.,forceAnalyze=.TRUE.,OutPut=.FALSE.)

!-----------------------------------------------------------------------------------------------------------------------------------
! iterations starting up from here
!-----------------------------------------------------------------------------------------------------------------------------------
DO !iter_t=0,MaxIter 

#if (PP_TimeDiscMethod==201)
!  IF (vMax.EQ.0) THEN
!    dt_max_particles = dt_maxwell
!  ELSE
!    dt_max_particles = 3.8*dt_maxwell*c/(vMax)!the 3.8 is a factor that lead to a timestep of 2/3*min_celllength
!                      ! this factor should be adjusted
!  END IF
  IF (iter.LE.MaximumIterNum) THEN
      dt_max_particles = dt_maxwell ! initial evolution of field with maxwellts
  ELSE
    DO 
      vMaxx = 0.
      vMaxy = 0.
      vMaxz = 0.
      DO iPart=1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(iPart)) THEN
          vMaxx = MAX( vMaxx , ABS(PartState(iPart, 4) + dt_temp*Pt(iPart,1)) )
          vMaxy = MAX( vMaxy , ABS(PartState(iPart, 5) + dt_temp*Pt(iPart,2)) )
          vMaxz = MAX( vMaxz , ABS(PartState(iPart, 6) + dt_temp*Pt(iPart,3)) )
        END IF
      END DO
!! -- intrinsic logical->real/int-conversion should be avoided!!!
!      vMaxx = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 4) + dt_temp*Pt(1:PDM%ParticleVecLength,1))))
!      vMaxy = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 5) + dt_temp*Pt(1:PDM%ParticleVecLength,2))))
!      vMaxz = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 6) + dt_temp*Pt(1:PDM%ParticleVecLength,3))))
      vMax = MAX(vMaxx,vMaxy,vMaxz,1.0) 
#ifdef MPI
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,vMax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
#endif /*MPI*/

      dt_max_particles = 3.8*dt_maxwell*c/(vMax)
      dt_temp = (dt_max_particles+dt_temp)/2
      IF((dt_temp.GE.dt_max_particles*0.95).AND.(dt_temp.LE.dt_max_particles*1.05)) EXIT
    END DO
  END IF
  dt_Min = dt_max_particles
  MaxwellIterNum = INT(dt_max_particles / dt_maxwell)
  IF (MaxwellIterNum.GT.MaximumIterNum) MaxwellIterNum = MaximumIterNum
  IterDisplayStep = MAX(INT(IterDisplayStepUser/MaxwellIterNum),1) !IterDisplayStep refers to dt_maxwell
  IF ((MPIroot).AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
    SWRITE(UNIT_StdOut,'(132("!"))')
    !print*, 'New IterNum for MaxwellSolver: ', MaxwellIterNum 
    !print*, 'New Particle TimeStep: ', dt_max_particles
  END IF
#endif

#ifdef PARTICLES
  IF(enableParticleMerge) THEN
    IF ((iter.GT.0).AND.(MOD(iter,vMPFMergeParticleIter).EQ.0)) doParticleMerge=.true.
  END IF
#endif /*PARTICLES*/

  tAnalyzeDiff=tAnalyze-t    ! time to next analysis, put in extra variable so number does not change due to numerical errors
  tEndDiff=tEnd-t            ! dito for end time
  dt=MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/))
  IF (tAnalyzeDiff-dt.LT.dt/100.0) dt = tAnalyzeDiff
  IF (tEndDiff-dt.LT.dt/100.0) dt = tEndDiff
  IF ( dt .LT. 0. ) THEN
    SWRITE(UNIT_StdOut,*)'*** ERROR: Is something wrong with the defined tEnd?!? ***'
    CALL abort(__STAMP__,&
      'Error in tEndDiff or tAnalyzeDiff!')
  END IF

! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
#if (PP_TimeDiscMethod==1)
  CALL TimeStepByLSERK(t,iter,tEndDiff)
#elif (PP_TimeDiscMethod==2)
  CALL TimeStepByLSERK(t,iter,tEndDiff)
#elif (PP_TimeDiscMethod==3)
  CALL TimeStepByTAYLOR(t)
#elif (PP_TimeDiscMethod==4)
  CALL TimeStep_DSMC(t)
#elif (PP_TimeDiscMethod==5)
  CALL TimeStepByRK4EulerExpl(t)
#elif (PP_TimeDiscMethod==6)
  CALL TimeStepByLSERK(t,iter,tEndDiff)
#elif (PP_TimeDiscMethod==42)
  CALL TimeStep_DSMC_Debug(t) ! Reservoir and Debug
#elif (PP_TimeDiscMethod==100)
  CALL TimeStepByEulerImplicit(t) ! O1 Euler Implicit
#elif (PP_TimeDiscMethod==101)
  CALL TimeStepByESDIRK4(t) ! O4 Implicit Runge Kutta
#elif (PP_TimeDiscMethod==200)
  CALL TimeStepByEulerStaticExp(t) ! O1 Euler Static Explicit
#elif (PP_TimeDiscMethod==201)
  CALL TimeStepByEulerStaticExpAdapTS(t) ! O1 Euler Static Explicit with adaptive TimeStep
#elif (PP_TimeDiscMethod==1000)
  CALL TimeStep_LD(t)
#elif (PP_TimeDiscMethod==1001)
  CALL TimeStep_LD_DSMC(t)
#endif
  iter=iter+1
  iter_loc=iter_loc+1
  t=t+dt
#ifdef PARTICLES
  realtime=realtime+dt
#endif /*PARTICLES*/
  IF(MPIroot) THEN
    IF(DoDisplayIter)THEN
      IF(MOD(iter,IterDisplayStep).EQ.0) THEN
         SWRITE(*,*) "iter:", iter,"t:",t
      END IF
    END IF
  END IF
#if (PP_TimeDiscMethod!=1) &&  (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6)
  ! calling the analyze routines
  CALL PerformAnalyze(t,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)
#endif
  ! output of state file
  !IF ((dt.EQ.tAnalyzeDiff).OR.(dt.EQ.tEndDiff)) THEN   ! timestep is equal to time to analyze or end
  IF((ALMOSTEQUAL(dt,tAnalyzeDiff)).OR.(ALMOSTEQUAL(dt,tEndDiff)))THEN
    IF(ALMOSTEQUAL(dt,tEndDiff))THEN
      finalIter=.TRUE.
    ELSE
      finalIter=.FALSE.
    END IF
    CalcTimeEnd=BOLTZPLATZTIME()
    IF(MPIroot)THEN
      ! Get calculation time per DOF
      CalcTimeEnd=(CalcTimeEnd-CalcTimeStart)*nProcessors/(nGlobalElems*(PP_N+1)**3*iter_loc)
      CALL DATE_AND_TIME(values=TimeArray) ! get System time
      WRITE(UNIT_StdOut,'(132("-"))')
      WRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
        ' Sys date  :    ',timeArray(3),'.',timeArray(2),'.',timeArray(1),' ',timeArray(5),':',timeArray(6),':',timeArray(7)
      WRITE(UNIT_stdOut,'(A,ES12.5,A)')' CALCULATION TIME PER TSTEP/DOF: [',CalcTimeEnd,' sec ]'
      WRITE(UNIT_StdOut,'(A,ES16.7)')' Timestep  : ',dt_Min
      WRITE(UNIT_stdOut,'(A,ES16.7)')'#Timesteps : ',REAL(iter)
    END IF !MPIroot
    ! Analyze for output
    CALL PerformAnalyze(t,iter,tenddiff,forceAnalyze=.FALSE.,OutPut=.TRUE.,LastIter=finalIter)
    ! future time
    nAnalyze=nAnalyze+1
    tFuture=tZero+REAL(nAnalyze)*Analyze_dt
    ! Write state to file
    IF(DoPML) CALL BacktransformPMLVars()
    CALL WriteStateToHDF5(TRIM(MeshFile),t,tFuture)
    IF(DoPML) CALL TransformPMLVars()
    ! Write recordpoints data to hdf5
    IF(RP_onProc) CALL WriteRPtoHDF5(tAnalyze,.TRUE.)
    iter_loc=0
    tAnalyze=tZero+REAL(nAnalyze)*Analyze_dt
    IF (tAnalyze > tEnd) tAnalyze = tEnd
    SWRITE(UNIT_StdOut,'(132("-"))')
    CalcTimeStart=BOLTZPLATZTIME()
  ENDIF   
  IF(t.GE.tEnd) THEN  ! done, worst case: one additional time step
    EXIT
  END IF
END DO ! iter_t


!CALL FinalizeAnalyze
END SUBROUTINE TimeDisc

#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
SUBROUTINE TimeStepByLSERK(t,iter,tEndDiff)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_Analyze,                 ONLY: PerformAnalyze
USE MOD_TimeDisc_Vars,           ONLY: dt,iStage
USE MOD_TimeDisc_Vars,           ONLY: RK_a,RK_b,RK_c,nRKStages
USE MOD_DG_Vars,                 ONLY: U,Ut,nTotalU
USE MOD_PML_Vars,                ONLY: U2,U2t,nPMLElems,DoPML,nTotalPML
USE MOD_PML,                     ONLY: PMLTimeDerivative,CalcPMLSource
USE MOD_Filter,                  ONLY: Filter
USE MOD_Equation,                ONLY: DivCleaningDamping
USE MOD_Equation,                ONLY: CalcSource
!USE MOD_DG,                      ONLY: DGTimeDerivative_WoSource_weakForm
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
#ifdef PP_POIS
USE MOD_Equation,                ONLY: DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,           ONLY: Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_PICDepo,                 ONLY: Deposition!, DepositionMPF
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_PIC_Vars,                ONLY: PIC
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime, Time, PEM, PDM, usevMPF, & 
                                        doParticleMerge,PartPressureCell
USE MOD_Particle_Vars,           ONLY: nTotalPart,nTotalHalfPart
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS, DSMC
USE MOD_part_MPFtools,           ONLY: StartParticleMerge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#else /*No MPI*/
#endif /*MPI*/
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: tendDiff
REAL,INTENT(INOUT)            :: t
INTEGER(KIND=8),INTENT(INOUT) :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                       :: iPart
REAL                          :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
REAL                          :: U2t_temp(1:6,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems) ! temporal variable for U2t
REAL                          :: tStage,b_dt(1:nRKStages)
#ifdef PARTICLES
REAL                          :: timeStart,timeEnd
#endif /*PARTICLES*/
#ifdef PP_POIS
REAL                          :: Phit_temp(1:4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
!===================================================================================================================================

! RK coefficients
DO iStage=1,nRKStages
  b_dt(iStage)=RK_b(iStage)*dt
END DO
iStage=1

#ifdef PARTICLES
Time=t

IF (t.GE.DelayTime) THEN
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleInserting()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tLocalization=tLocalization+TimeEnd-TimeStart
  END IF
END IF

IF (t.GE.DelayTime) THEN
  ! forces on particle
  ! can be used to hide sending of number of particles
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  CALL CalcPartRHS()
END IF

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
  CALL Deposition(doInnerParts=.FALSE.)
#endif /*MPI*/
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF

#endif /*PARTICLES*/

! field solver
CALL DGTimeDerivative_weakForm(t,t,0)
IF(DoPML) THEN
  CALL CalcPMLSource()
  CALL PMLTimeDerivative()
END IF
CALL DivCleaningDamping()

#ifdef PP_POIS
! Potential
CALL DGTimeDerivative_weakForm_Pois(t,t,0)
CALL DivCleaningDamping_Pois()
#endif

! calling the analyze routines
CALL PerformAnalyze(t,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)

! first RK step
! EM field
Ut_temp = Ut 
U = U + Ut*b_dt(1)
!PML auxiliary variables
IF(DoPML) THEN
  U2t_temp = U2t
  U2 = U2 + U2t*b_dt(1)
END IF

#ifdef PP_POIS
Phit_temp = Phit 
Phi = Phi + Phit*b_dt(1)
CALL EvalGradient()
#endif

#ifdef PARTICLES
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN
  Pt_temp(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,4) 
  Pt_temp(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,5) 
  Pt_temp(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,6) 
  Pt_temp(1:PDM%ParticleVecLength,4) = Pt(1:PDM%ParticleVecLength,1) 
  Pt_temp(1:PDM%ParticleVecLength,5) = Pt(1:PDM%ParticleVecLength,2) 
  Pt_temp(1:PDM%ParticleVecLength,6) = Pt(1:PDM%ParticleVecLength,3)
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                       + PartState(1:PDM%ParticleVecLength,4)*b_dt(1)
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                       + PartState(1:PDM%ParticleVecLength,5)*b_dt(1)
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                       + PartState(1:PDM%ParticleVecLength,6)*b_dt(1)
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                       + Pt(1:PDM%ParticleVecLength,1)*b_dt(1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                       + Pt(1:PDM%ParticleVecLength,2)*b_dt(1)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                       + Pt(1:PDM%ParticleVecLength,3)*b_dt(1)
END IF
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
#ifdef MPI
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTrackingCurved()
  END IF
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tTracking=tTracking+TimeEnd-TimeStart
  END IF
#ifdef MPI
  CALL SendNbOfParticles()
  !CALL MPIParticleSend()
  !CALL MPIParticleRecv()
  !PartMPIExchange%nMPIParticles=0
#endif
  !CALL Filter(U)
END IF
#endif /*PARTICLES*/

DO iStage=2,nRKStages
  tStage=t+dt*RK_c(iStage)
#ifdef PARTICLES
  ! deposition  
  IF (t.GE.DelayTime) THEN 
     CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
#ifdef MPI
     CALL MPIParticleSend()
#endif /*MPI*/
!    ! deposition  
     CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
     CALL MPIParticleRecv()
#endif /*MPI*/
     ! second buffer
     CALL InterpolateFieldToParticle(doInnerParts=.FALSE.)
     CALL CalcPartRHS()
     CALL Deposition(doInnerParts=.FALSE.)
     IF(DoVerifyCharge) CALL VerifyDepositedCharge()
     ! null here, careful
#ifdef MPI
     PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
!    IF (usevMPF) THEN 
!      CALL !DepositionMPF()
!    ELSE 
!      CALL Deposition()
!    END IF
  END IF
#endif /*PARTICLES*/

  ! field solver
  CALL DGTimeDerivative_weakForm(t,tStage,0)
  IF(DoPML) THEN
    CALL CalcPMLSource()
    CALL PMLTimeDerivative()
  END IF
  CALL DivCleaningDamping()

#ifdef PP_POIS
  CALL DGTimeDerivative_weakForm_Pois(t,tStage,0)
  CALL DivCleaningDamping_Pois()
#endif

  ! performe RK steps
  ! field step
  Ut_temp = Ut - Ut_temp*RK_a(iStage)
  U = U + Ut_temp*b_dt(iStage)
  !PML auxiliary variables
  IF(DoPML)THEN
    U2t_temp = U2t - U2t_temp*RK_a(iStage)
    U2 = U2 + U2t_temp*b_dt(iStage)
  END IF

#ifdef PP_POIS
  Phit_temp = Phit - Phit_temp*RK_a(iStage)
  Phi = Phi + Phit_temp*b_dt(iStage)
  CALL EvalGradient()
#endif

#ifdef PARTICLES
  ! particle step
  IF (t.GE.DelayTime) THEN
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    Pt_temp(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,4) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,1)
    Pt_temp(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,5) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,2)
    Pt_temp(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,6) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,3)
    Pt_temp(1:PDM%ParticleVecLength,4) = Pt(1:PDM%ParticleVecLength,1) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,4)
    Pt_temp(1:PDM%ParticleVecLength,5) = Pt(1:PDM%ParticleVecLength,2) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,5)
    Pt_temp(1:PDM%ParticleVecLength,6) = Pt(1:PDM%ParticleVecLength,3) &
                             - RK_a(iStage) * Pt_temp(1:PDM%ParticleVecLength,6)
    PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                       + Pt_temp(1:PDM%ParticleVecLength,1)*b_dt(iStage)
    PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                       + Pt_temp(1:PDM%ParticleVecLength,2)*b_dt(iStage)
    PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                       + Pt_temp(1:PDM%ParticleVecLength,3)*b_dt(iStage)
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                       + Pt_temp(1:PDM%ParticleVecLength,4)*b_dt(iStage)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                       + Pt_temp(1:PDM%ParticleVecLength,5)*b_dt(iStage)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                       + Pt_temp(1:PDM%ParticleVecLength,6)*b_dt(iStage)
    ! particle tracking
#ifdef MPI
    CALL IRecvNbofParticles()
#endif /*MPI*/
    ! actual tracking
!    CALL ParticleTrackingCurved()
    IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
    IF(DoRefMapping)THEN
      CALL ParticleRefTracking()
    ELSE
      CALL ParticleTrackingCurved()
    END IF
    IF(MeasureTrackTime) THEN
      CALL CPU_TIME(TimeEnd)
      tTracking=tTracking+TimeEnd-TimeStart
    END IF
#ifdef MPI
    CALL SendNbOfParticles()
    !CALL MPIParticleSend()
    !CALL MPIParticleRecv()
    !PartMPIExchange%nMPIParticles=0
#endif
  END IF
#endif /*PARTICLES*/

END DO

#ifdef PARTICLES
#ifdef MPI
IF (t.GE.DelayTime) THEN
  CALL MPIParticleSend()
  CALL MPIParticleRecv()
  PartMPIExchange%nMPIParticles=0
END IF
#endif /*MPI*/

IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF

IF (doParticleMerge) THEN
  CALL StartParticleMerge()  
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )
  END IF
  CALL UpdateNextFreePosition()
END IF

IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
#endif /*PARTICLES*/


! vectorisazation is commented out! 
! 
! 
! ! RK coefficients
! DO iStage=1,nRKStages
!   b_dt(iStage)=RK_b(iStage)*dt
! END DO
! iStage=1
! 
! #ifdef PARTICLES
! Time=t
! IF (t.GE.DelayTime) THEN
!   SWRITE(*,*) 'emission'
!   CALL ParticleInserting()
!   SWRITE(*,*) 'emission done'
!   nTotalHalfPart=PDM%ParticleVecLength*3
!   nTotalPart    =PDM%ParticleVecLength*6
! END IF
! 
! IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!   ! because of emmision and UpdateParticlePosition
!   CALL Deposition(doInnerParts=.TRUE.)
! !#ifdef MPI
! !  CALL MPIParticleRecv()
! !  ! second buffer
! !  CALL Deposition(doInnerParts=.FALSE.)
! !#endif /*MPI*/
! END IF
! 
! IF (t.GE.DelayTime) THEN
!   ! forces on particle
!   CALL InterpolateFieldToParticle()
!   CALL CalcPartRHS()
! END IF
! #endif /*PARTICLES*/
! 
! ! field solver
! CALL DGTimeDerivative_weakForm(t,t,0)
! IF(DoPML) THEN
!   CALL CalcPMLSource()
!   CALL PMLTimeDerivative()
! END IF
! CALL DivCleaningDamping()
! 
! #ifdef PP_POIS
! ! Potential
! CALL DGTimeDerivative_weakForm_Pois(t,t,0)
! CALL DivCleaningDamping_Pois()
! #endif
! 
! ! calling the analyze routines
! CALL PerformAnalyze(t,iter,tendDiff,force=.FALSE.)
! 
! ! first RK step
! ! EM field
! !Ut_temp = Ut 
! !U = U + Ut*b_dt(1)
! CALL VCOPY(nTotalU,Ut,Ut_temp)
! CALL VAXPBY(nTotalU,Ut,U,ConstIn=b_dt(1))
! 
! !PML auxiliary variables
! IF(DoPML) THEN
!   CALL VCOPY(nTotalPML,U2t,U2t_temp)
!   CALL VAXPBY(nTotalPML,U2t,U2,ConstIn=b_dt(1))
! END IF
! 
! #ifdef PP_POIS
!   CALL VCOPY(nTotalU,Phit,Phit_temp)
!   CALL VAXPBY(nTotalU,Phit,Phi,ConstIn=b_dt(1))
!   CALL EvalGradient()
! #endif
! 
! #ifdef PARTICLES
! ! particles
! CALL LVCOPY(nTotalHalfPart,PartState(1:PDM%ParticleVecLength,1:3),LastPartPos(1:PDM%ParticleVecLength,1:3))
! DO iPart=1,PDM%ParticleVecLength
!   PEM%LastElement(iPart)=PEM%Element(iPart)
! END DO
! IF (t.GE.DelayTime) THEN
!   CALL LVCOPY(nTotalHalfPart,PartState(1:PDM%ParticleVecLength,4:6),Pt_temp(1:PDM%ParticleVecLength,1:3))
!   CALL LVCOPY(nTotalHalfPart,Pt(1:PDM%ParticleVecLength,:),Pt_temp(1:PDM%ParticleVecLength,4:6))
!   CALL LVAXPBY(nTotalHalfPart,PartState(1:PDM%ParticleVecLength,4:6),PartState(1:PDM%ParticleVecLength,1:3),ConstIn=b_dt(1))
!   CALL LVAXPBY(nTotalHalfPart,Pt(1:PDM%ParticleVecLength,:),PartState(1:PDM%ParticleVecLength,4:6),ConstIn=b_dt(1))
! END IF
! 
! IF (t.GE.DelayTime) THEN
! !CALL ParticleTracking()
! #ifdef MPI
! CALL IRecvNbofParticles()
! #endif /*MPI*/
! CALL ParticleTrackingCurved()
! #ifdef MPI
! CALL MPIParticleSend()
! #endif /*MPI*/
! END IF
! #endif /*PARTICLES*/
! 
! DO iStage=2,nRKStages
!   SWRITE(*,*) " iState", iStage
!   tStage=t+dt*RK_c(iStage)
!   IF (t.GE.DelayTime) THEN 
! #ifdef PARTICLES
!    ! deposition  
!     CALL Deposition(doInnerParts=.TRUE.)
! #ifdef MPI
!     CALL MPIParticleRecv()
!     ! second buffer
!     CALL Deposition(doInnerParts=.FALSE.)
! #endif /*MPI*/
!   ! particle RHS
!     CALL InterpolateFieldToParticle()
!     CALL CalcPartRHS()
!   END IF
! #endif /*PARTICLES*/
! 
! 
!   ! field solver
!   CALL DGTimeDerivative_weakForm(t,tStage,0)
!   IF(DoPML) THEN
!     CALL CalcPMLSource()
!     CALL PMLTimeDerivative()
!   END IF
!   CALL DivCleaningDamping()
! 
! #ifdef PP_POIS
!   CALL DGTimeDerivative_weakForm_Pois(t,tStage,0)
!   CALL DivCleaningDamping_Pois()
! #endif
! 
!   ! performe RK steps
!   ! field step
! 
!   !Ut_temp = Ut - Ut_temp*RK4_a(iStage)
!   !U = U + Ut_temp*b_dt(iStage)
!   CALL VAXPBY(nTotalU,Ut,Ut_temp,ConstOut=-RK_a(iStage))
!   CALL VAXPBY(nTotalU,Ut_temp,U,ConstIn=b_dt(iStage))
! 
!   !PML auxiliary variables
!   IF(DoPML)THEN
!     CALL VAXPBY(nTotalPML,U2t,U2t_temp,ConstOut=-RK_a(iStage))
!     CALL VAXPBY(nTotalPML,U2t_temp,U2,ConstIn=b_dt(iStage))
!   END IF
! 
! #ifdef PP_POIS
!   CALL VAXPBY(nTotalU,Phit,Phit_temp,ConstOut=-RK_a(iStage))
!   CALL VAXPBY(nTotalU,Phit_temp,Phi,ConstIn=b_dt(iStage))
!   CALL EvalGradient()
! #endif
! 
! #ifdef PARTICLES
!   ! particle step
!   IF (t.GE.DelayTime) THEN
!     CALL LVCOPY(nTotalHalfPart,PartState(1:PDM%particleVecLength,1:3),LastPartPos(1:PDM%ParticleVecLength,1:3))
!     DO iPart=1,PDM%ParticleVecLength
!       PEM%LastElement(iPart)=PEM%Element(iPart)
!     END DO
!  CALL LVAXPBY(nTotalHalfPart,PartState(1:PDM%ParticleVecLength,4:6),Pt_temp(1:PDM%ParticleVecLength,1:3),ConstOut=-RK_a(iStage))
!  CALL LVAXPBY(nTotalHalfPart,Pt(1:PDM%ParticleVecLength,1:3),Pt_temp(1:PDM%ParticleVecLength,4:6),ConstOut=-RK_a(iStage))
!  CALL LVAXPBY(nTotalPart,Pt_temp(1:PDM%ParticleVecLength,1:6),PartState(1:PDM%ParticleVecLength,1:6),ConstIn=b_dt(iStage))
! 
!     ! particle tracking
!     !CALL ParticleBoundary()
!    !CALL ParticleTracking()
! #ifdef MPI
!     CALL IRecvNbofParticles()
! #endif
!     CALL ParticleTrackingCurved()
! #ifdef MPI
!     CALL MPIParticleSend()
! #endif
!   END IF
! #endif /*PARTICLES*/
! 
! 
! END DO
! 
! !IF (t.GE.DelayTime) THEN
! !  CALL Deposition(doInnerParts=.TRUE.)
! !#ifdef MPI
! !  CALL MPIParticleRecv()
! !  ! second buffer
! !  CALL Deposition(doInnerParts=.FALSE.)
! !#endif /*MPI*/
! #ifdef PARTICLES
! IF (t.GE.DelayTime) THEN
! #ifdef MPI
!   CALL MPIParticleRecv()
! #endif /*MPI*/
! END IF
! 
! IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!   CALL UpdateNextFreePosition()
!   nTotalHalfPart=PDM%ParticleVecLength*3
!   nTotalPart    =PDM%ParticleVecLength*6
! END IF
! 
! !print*,'time',t1,t2
! 
! IF (useDSMC) THEN
!   IF (t.GE.DelayTime) THEN
!     CALL DSMC_main()
!     CALL LVAXPBY(nTotalHalfPart,DSMC_RHS(1:PDM%ParticleVecLength,1:3),PartState(1:PDM%ParticleVecLength,4:6))
!   END IF
! END IF
! #endif /*PARTICLES*/

END SUBROUTINE TimeStepByLSERK
#endif

#if (PP_TimeDiscMethod==3) 
SUBROUTINE TimeStepByTAYLOR(t)
!===================================================================================================================================
! TAYLOR DG
! time order = N+1
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut,nTotalU
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY:dt
USE MOD_DG,ONLY:DGTimeDerivative_weakForm
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: U_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
REAL             :: Factor
INTEGER          :: tIndex,i
!===================================================================================================================================
Factor = dt
#ifdef OPTIMIZED
DO i=1,nTotalU
  U_temp(i,0,0,0,1) = U(i,0,0,0,1)
END DO
#else
U_temp=U
#endif OPTIMIZED

DO tIndex=0,PP_N
 CALL DGTimeDerivative_weakForm(t,t,tIndex)
#ifdef OPTIMIZED
DO i=1,nTotalU
  U(i,0,0,0,1) = Ut(i,0,0,0,1)
  U_temp(i,0,0,0,1) = U_temp(i,0,0,0,1) + Factor*Ut(i,0,0,0,1)
END DO
#else
 U=Ut
 U_temp = U_temp + Factor*Ut
#endif
 Factor = Factor*dt/REAL(tIndex+2)
END DO ! tIndex

#ifdef OPTIMIZED
DO i=1,nTotalU
  U(i,0,0,0,1) = U_temp(i,0,0,0,1)
END DO
#else
U=U_temp
#endif
END SUBROUTINE TimeStepByTAYLOR
#endif

#if (PP_TimeDiscMethod==4)
SUBROUTINE TimeStep_DSMC(t)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY: dt
USE MOD_Filter,ONLY:Filter
#ifdef PARTICLES
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, Time, PDM,PEM
USE MOD_DSMC_Vars,        ONLY : DSMC_RHS
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_Particle_Tracking_vars, ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,ONLY: ParticleTrackingCurved,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,     ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,ONLY: PartMPIExchange
#endif /*MPI*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: timeEnd, timeStart
!===================================================================================================================================
Time = t

  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)

  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                         + PartState(1:PDM%ParticleVecLength,4) * dt
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                         + PartState(1:PDM%ParticleVecLength,5) * dt
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                         + PartState(1:PDM%ParticleVecLength,6) * dt

#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  ! actual tracking
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTrackingCurved()
  END IF
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tTracking=tTracking+TimeEnd-TimeStart
  END IF
#ifdef MPI
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
#endif /*MPI*/

  CALL ParticleInserting()
  CALL UpdateNextFreePosition()
  CALL DSMC_main()

  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)

END SUBROUTINE TimeStep_DSMC
#endif


#if (PP_TimeDiscMethod==5)
SUBROUTINE TimeStepByRK4EulerExpl(t)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY: U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY: dt
USE MOD_TimeDisc_Vars,ONLY: RK_a,RK_b,RK_c
USE MOD_DG,ONLY:DGTimeDerivative_weakForm
USE MOD_Filter,ONLY:Filter
USE MOD_Equation,ONLY:DivCleaningDamping
#ifdef PP_POIS
USE MOD_Equation,ONLY:DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,ONLY:DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,ONLY:Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos,Time, PEM, PDM, usevMPF, doParticleMerge, DelayTime, PartPressureCell
USE MOD_part_RHS,         ONLY : CalcPartRHS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS
USE MOD_part_MPFtools,    ONLY : StartParticleMerge
!#ifdef MPI
!!USE MOD_part_boundary,    ONLY : ParticleBoundary, Communicate_PIC
!#else
!USE MOD_part_boundary,    ONLY : ParticleBoundary
!#endif
USE MOD_PIC_Analyze,      ONLY: VerifyDepositedCharge
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
#ifdef PP_POIS
REAL                  :: Phit_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
REAL                  :: tStage,b_dt(1:5)
INTEGER               :: rk
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()
!CALL UpdateNextFreePosition()
DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO

!IF(t.EQ.0) CALL Deposition()
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  !IF (usevMPF) THEN 
  !  CALL DepositionMPF()
  !ELSE 
    CALL Deposition()
  !END IF
  !CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles 
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + dt * PartState(1:PDM%ParticleVecLength,4) 
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + dt * PartState(1:PDM%ParticleVecLength,5) 
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + dt * PartState(1:PDM%ParticleVecLength,6) 
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + dt * Pt(1:PDM%ParticleVecLength,1) 
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + dt * Pt(1:PDM%ParticleVecLength,2) 
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + dt * Pt(1:PDM%ParticleVecLength,3) 
END IF
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!  CALL ParticleBoundary()
!#ifdef MPI
!  CALL Communicate_PIC()
!!CALL UpdateNextFreePosition() ! only required for parallel communication
!#endif
END IF

! EM field
CALL DGTimeDerivative_weakForm(t,t,0)
CALL DivCleaningDamping()
Ut_temp = Ut 
U = U + Ut*b_dt(1)

DO rk=2,5
  tStage=t+dt*RK_c(rk)
  ! field RHS
  CALL DGTimeDerivative_weakForm(t,tStage,0)
  CALL DivCleaningDamping()
  ! field step
  Ut_temp = Ut - Ut_temp*RK_a(rk)
  U = U + Ut_temp*b_dt(rk)
END DO

#ifdef PP_POIS
! EM field
CALL DGTimeDerivative_weakForm_Pois(t,t,0)

CALL DivCleaningDamping_Pois()
Phit_temp = Phit 
Phi = Phi + Phit*b_dt(1)

DO rk=2,5
  tStage=t+dt*RK_c(rk)
  ! field RHS
  CALL DGTimeDerivative_weakForm_Pois(t,tStage,0)
  CALL DivCleaningDamping_Pois()
  ! field step
  Phit_temp = Phit - Phit_temp*RK_a(rk)
  Phi = Phi + Phit_temp*b_dt(rk)
END DO

  CALL EvalGradient()
#endif

IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF

IF (doParticleMerge) THEN
  CALL StartParticleMerge()  
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )
  END IF
  CALL UpdateNextFreePosition()
END IF

IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
END SUBROUTINE TimeStepByRK4EulerExpl
#endif

#if (PP_TimeDiscMethod==42)
SUBROUTINE TimeStep_DSMC_Debug(t)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY: dt
USE MOD_Filter,ONLY:Filter
#ifdef PARTICLES
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, Time, PDM,PEM
USE MOD_DSMC_Vars,        ONLY : DSMC_RHS, DSMC
USE MOD_DSMC,             ONLY : DSMC_main
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting
#endif
!#ifdef MPI
!USE MOD_part_boundary,    ONLY : Communicate_PIC
!#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!===================================================================================================================================
Time = t

IF (DSMC%ReservoirSimu) THEN ! fix grid should be defined for reservoir simu
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
ELSE
  CALL ParticleInserting()
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  ! bugfix if more than 2.x mio (2000001) particle per process
  ! tested with 64bit Ubuntu 12.04 backports
  DO i = 1, PDM%ParticleVecLength
    PEM%lastElement(i)=PEM%Element(i)
  END DO
!   PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                         + PartState(1:PDM%ParticleVecLength,4) * dt
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                         + PartState(1:PDM%ParticleVecLength,5) * dt
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                         + PartState(1:PDM%ParticleVecLength,6) * dt
!  CALL ParticleBoundary()
!#ifdef MPI
!  CALL Communicate_PIC()
!#endif
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
END IF

END SUBROUTINE TimeStep_DSMC_Debug
#endif

#if (PP_TimeDiscMethod==100) 
SUBROUTINE TimeStepByEulerImplicit(t)
!===================================================================================================================================
! Euler Implicit method:
! U^n+1 = U^n + dt*R(U^n+1)
! (I -dt*R)*U^n+1 = U^n
! Solve Linear System 
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars,ONLY:dt
USE MOD_DG_Vars                 ,ONLY:U,Ut
USE MOD_DG                      ,ONLY:DGTimeDerivative_weakForm
USE MOD_LinearSolver,     ONLY : LinearSolver
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_PIC_Vars,         ONLY : PIC
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, usevMPF
USE MOD_part_RHS,         ONLY : CalcPartRHS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS, DSMC
!#ifdef MPI
!USE MOD_part_boundary,    ONLY : ParticleBoundary, Communicate_PIC
!#else
!USE MOD_part_boundary,    ONLY : ParticleBoundary
!#endif
USE MOD_PIC_Analyze,ONLY: VerifyDepositedCharge
USE MOD_part_tools,     ONLY : UpdateNextFreePosition
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!  IF (usevMPF) THEN 
!    CALL DepositionMPF()
!  ELSE 
    CALL Deposition()
!  END IF
  !CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles 
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + dt * PartState(1:PDM%ParticleVecLength,4) 
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + dt * PartState(1:PDM%ParticleVecLength,5) 
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + dt * PartState(1:PDM%ParticleVecLength,6) 
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + dt * Pt(1:PDM%ParticleVecLength,1) 
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + dt * Pt(1:PDM%ParticleVecLength,2) 
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + dt * Pt(1:PDM%ParticleVecLength,3) 
END IF
!CALL ParticleBoundary()
!#ifdef MPI
!CALL Communicate_PIC()
!!CALL UpdateNextFreePosition() ! only required for parallel communication
!#endif
! EM field
CALL LinearSolver(t,dt)
CALL DivCleaningDamping()
CALL UpdateNextFreePosition()
IF (useDSMC) THEN
  CALL DSMC_main()
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
END IF

END SUBROUTINE TimeStepByEulerImplicit
#endif /*PP_TimeDiscMethod==100*/

#if (PP_TimeDiscMethod==101) 
SUBROUTINE TimeStepByESDIRK4(t)
!===================================================================================================================================
! RK4 IMPLICIT
! time order = N+1
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
! Butcher Tableau taken from Bijl, Carpenter, Vatsa, Kennedy - doi:10.1006/jcph.2002.7059
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DG_Vars                 ,ONLY:U,Ut
USE MOD_DG                      ,ONLY:DGTimeDerivative_weakForm
USE MOD_PreProc
USE MOD_TimeDisc_Vars
USE MOD_Equation,ONLY:DivCleaningDamping
USE MOD_Equation,      ONLY: CalcSource
USE MOD_LinearSolver,  ONLY : LinearSolver
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_PIC_Vars,         ONLY : PIC
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, usevMPF
USE MOD_part_RHS,         ONLY : CalcPartRHS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS, DSMC
!#ifdef MPI
!USE MOD_part_boundary,    ONLY : ParticleBoundary, Communicate_PIC
!#else
!USE MOD_part_boundary,    ONLY : ParticleBoundary
!#endif
!USE MOD_PIC_Analyze,    ONLY: VerifyDepositedCharge
USE MOD_part_tools,     ONLY : UpdateNextFreePosition
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)

REAL     :: K1(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL     :: K2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL     :: K3(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL     :: K4(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL     :: K5(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)

REAL     :: tStage

Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!    IF (usevMPF) THEN 
!      CALL DepositionMPF()
!    ELSE 
      CALL Deposition()
!    END IF
  !CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles 
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + dt * PartState(1:PDM%ParticleVecLength,4) 
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + dt * PartState(1:PDM%ParticleVecLength,5) 
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + dt * PartState(1:PDM%ParticleVecLength,6) 
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + dt * Pt(1:PDM%ParticleVecLength,1) 
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + dt * Pt(1:PDM%ParticleVecLength,2) 
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + dt * Pt(1:PDM%ParticleVecLength,3) 
END IF
!CALL ParticleBoundary()
!#ifdef MPI
!CALL Communicate_PIC()
!!CALL UpdateNextFreePosition() ! only required for parallel communication
!#endif
! EM field
!===================================================================================================================================
! Save state at t(n)
Un = U
! ----------------------------------------------------------------------------------------------------------------------------------
tStage = t
CALL DGTimeDerivative_weakForm(t, tStage, 0)
CALL DivCleaningDamping()
K1  = Ut
! ----------------------------------------------------------------------------------------------------------------------------------
tStage = t + 0.5 * dt
! Fixed part of the RHS
U = Un + dt*0.25*K1
CALL LinearSolver(tStage,0.25*dt)
CALL DGTimeDerivative_weakForm(t, tStage, 0)
CALL DivCleaningDamping()
K2=Ut
! ----------------------------------------------------------------------------------------------------------------------------------
tStage = t + 83./250. * dt
! Fixed part of the RHS
U = Un + dt*(8611./62500.*K1 - 1743./31250.*K2)
CALL LinearSolver(tStage,0.25*dt)
CALL DGTimeDerivative_weakForm(t, tStage, 0)
CALL DivCleaningDamping()
K3=Ut
! ----------------------------------------------------------------------------------------------------------------------------------
tStage = t + 31./50. * dt
U = Un + dt*( 5012029./34652500.*K1 - 654441./2922500.*K2 + 174375./388108.*K3)
CALL LinearSolver(tStage,0.25*dt)
CALL DGTimeDerivative_weakForm(t, tStage, 0)
CALL DivCleaningDamping()
K4=Ut
! ----------------------------------------------------------------------------------------------------------------------------------
tStage = t + 17./20. * dt
! Fixed part of the RHS
U = Un + dt*(15267082809./155376265600.*K1 - 71443401./120774400.*K2 + 730878875./902184768.*K3 + 2285395./8070912.*K4)
CALL LinearSolver(tStage,0.25*dt)
CALL DGTimeDerivative_weakForm(t, tStage, 0)
CALL DivCleaningDamping()
K5=Ut
!-----------------------------------------------------------------------------------------------------------------------------------
tStage = t + dt
! Fixed part of the RHS
U = Un + dt*(82889./524892.*K1 + 15625./83664.*K3 + 69875./102672.*K4 - 2260./8211.*K5)
CALL LinearSolver(tStage,0.25*dt)
CALL UpdateNextFreePosition()
IF (useDSMC) THEN
  CALL DSMC_main()
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
END IF

CALL DivCleaningDamping()
END SUBROUTINE TimeStepByESDIRK4
#endif /*PP_TimeDiscMethod==101*/

#if (PP_TimeDiscMethod==200)
SUBROUTINE TimeStepByEulerStaticExp(t)
!===================================================================================================================================
! Static (using 4 or 8 variables, depending on compiled equation system (maxwell or electrostatic):
! Field is propagated until steady, then particle is moved
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,          ONLY: U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,    ONLY: dt,IterDisplayStep,iter,IterDisplayStepUser
USE MOD_TimeDisc_Vars,    ONLY: RK_a,RK_b,RK_c
USE MOD_DG,               ONLY:DGTimeDerivative_weakForm
USE MOD_Filter,           ONLY:Filter
USE MOD_Equation,         ONLY:DivCleaningDamping
USE MOD_Globals
#ifdef PP_POIS
USE MOD_Equation,         ONLY:DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,               ONLY:DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,    ONLY:Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, dt_maxwell, MaxwellIterNum, usevMPF
USE MOD_part_RHS,         ONLY : CalcPartRHS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS
!#ifdef MPI
!USE MOD_part_boundary,    ONLY : ParticleBoundary, Communicate_PIC
!#else
!USE MOD_part_boundary,    ONLY : ParticleBoundary
!#endif
USE MOD_PIC_Analyze,      ONLY: VerifyDepositedCharge
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: rk, iLoop
REAL                  :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
#ifdef PP_POIS
REAL                  :: Phit_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
REAL                  :: b_dt(1:5)
REAL                  :: dt_save, tStage, t_rk
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!    IF (usevMPF) THEN 
!      CALL DepositionMPF()
!    ELSE 
      CALL Deposition()
    !END IF
  !CALL VerifyDepositedCharge()
END IF
IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles 
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + dt * PartState(1:PDM%ParticleVecLength,4) 
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + dt * PartState(1:PDM%ParticleVecLength,5) 
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + dt * PartState(1:PDM%ParticleVecLength,6) 
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + dt * Pt(1:PDM%ParticleVecLength,1) 
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + dt * Pt(1:PDM%ParticleVecLength,2) 
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + dt * Pt(1:PDM%ParticleVecLength,3) 
END IF
!CALL ParticleBoundary()
!#ifdef MPI
!CALL Communicate_PIC()
!CALL UpdateNextFreePosition() ! only required for parallel communication
!#endif
! EM field

dt_save = dt  !quick hack
t_rk = t
dt = dt_maxwell

DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO
IF(MOD(iter,IterDisplayStep).EQ.0) THEN
  SWRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
END IF
DO iLoop = 1, MaxwellIterNum
  ! EM field
  IF ((MOD(iter,IterDisplayStep).EQ.0).AND.(MOD(iLoop,IterDisplayStepUser).EQ.0)) THEN
    SWRITE(*,'(A14,I0,A3,I0)') " MaxwellIter: ",iLoop," / ",MaxwellIterNum
  END IF
  CALL DGTimeDerivative_weakForm(t_rk,t_rk,0)
  CALL DivCleaningDamping()
  Ut_temp = Ut 
  U = U + Ut*b_dt(1)
 
  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm(t_rk,tStage,0)
    CALL DivCleaningDamping()
    ! field step
    Ut_temp = Ut - Ut_temp*RK_a(rk)
    U = U + Ut_temp*b_dt(rk)
  END DO
#ifdef PP_POIS
  ! Potential field
  CALL DGTimeDerivative_weakForm_Pois(t_rk,t_rk,0)

  CALL DivCleaningDamping_Pois()
  Phit_temp = Phit
  Phi = Phi + Phit*b_dt(1)

  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm_Pois(t_rk,tStage,0)
    CALL DivCleaningDamping_Pois()
    ! field step
    Phit_temp = Phit - Phit_temp*RK_a(rk)
    Phi = Phi + Phit_temp*b_dt(rk)
  END DO
#endif
  t_rk = t_rk + dt
END DO
dt = dt_save
#ifdef PP_POIS
  CALL EvalGradient()
#endif
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF
IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
END SUBROUTINE TimeStepByEulerStaticExp
#endif

#if (PP_TimeDiscMethod==201)
SUBROUTINE TimeStepByEulerStaticExpAdapTS(t)
!===================================================================================================================================
! Static (using 4 or 8 variables, depending on compiled equation system (maxwell or electrostatic):
! Field is propagated until steady or a particle velo dependent adaptive time, then particle is moved
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY: U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY: dt,IterDisplayStep,iter,IterDisplayStepUser
USE MOD_TimeDisc_Vars,ONLY: RK_a,RK_b,RK_c
USE MOD_DG,ONLY:DGTimeDerivative_weakForm
USE MOD_Filter,ONLY:Filter
USE MOD_Equation,ONLY:DivCleaningDamping
USE MOD_Globals
#ifdef PP_POIS
USE MOD_Equation,ONLY:DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,ONLY:DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,ONLY:Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_PIC_Vars,         ONLY : PIC
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, dt_maxwell, MaxwellIterNum, usevMPF
USE MOD_part_RHS,         ONLY : CalcPartRHS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS, DSMC
!#ifdef MPI
!USE MOD_part_boundary,    ONLY : ParticleBoundary, Communicate_PIC
!#else
!USE MOD_part_boundary,    ONLY : ParticleBoundary
!#endif
USE MOD_PIC_Analyze,    ONLY: VerifyDepositedCharge
USE MOD_part_tools,     ONLY : UpdateNextFreePosition
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i,rk, iLoop
REAL                  :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
#ifdef PP_POIS
REAL                  :: Phit_temp(1:4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
REAL                  :: b_dt(1:5)
REAL                  :: dt_save, tStage, t_rk
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!    IF (usevMPF) THEN 
!      CALL DepositionMPF()
!    ELSE 
      CALL Deposition()
!    END IF
  !CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles 
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + dt * PartState(1:PDM%ParticleVecLength,4) 
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + dt * PartState(1:PDM%ParticleVecLength,5) 
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + dt * PartState(1:PDM%ParticleVecLength,6) 
  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + dt * Pt(1:PDM%ParticleVecLength,1) 
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + dt * Pt(1:PDM%ParticleVecLength,2) 
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + dt * Pt(1:PDM%ParticleVecLength,3) 
END IF
!CALL ParticleBoundary()
!#ifdef MPI
!CALL Communicate_PIC()
!!CALL UpdateNextFreePosition() ! only required for parallel communication
!#endif

! EM field

dt_save = dt  !quick hack
t_rk = t
dt = dt_maxwell

DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO
IF(MOD(iter,IterDisplayStep).EQ.0) THEN
  SWRITE(*,*) "!!!!!!!!!!!!"
END IF
DO iLoop = 1, MaxwellIterNum
  IF ((MOD(iter,IterDisplayStep).EQ.0).AND.(MOD(iLoop,IterDisplayStepUser).EQ.0)) THEN
    SWRITE(*,'(A14,I0,A3,I0)') " MaxwellIter: ",iLoop," / ",MaxwellIterNum
  END IF
  ! EM field

  CALL DGTimeDerivative_weakForm(t_rk,t_rk,0)
  CALL DivCleaningDamping()
  Ut_temp = Ut 
  U = U + Ut*b_dt(1)
 
  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm(t_rk,tStage,0)
    CALL DivCleaningDamping()
    ! field step
    Ut_temp = Ut - Ut_temp*RK_a(rk)
    U = U + Ut_temp*b_dt(rk)
  END DO
#ifdef PP_POIS
  ! Potential field
  CALL DGTimeDerivative_weakForm_Pois(t_rk,t_rk,0)

  CALL DivCleaningDamping_Pois()
  Phit_temp = Phit
  Phi = Phi + Phit*b_dt(1)

  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm_Pois(t_rk,tStage,0)
    CALL DivCleaningDamping_Pois()
    ! field step
    Phit_temp = Phit - Phit_temp*RK_a(rk)
    Phi = Phi + Phit_temp*b_dt(rk)
  END DO
#endif

  t_rk = t_rk + dt
END DO
dt = dt_save

#ifdef PP_POIS
  CALL EvalGradient()
#endif
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF
IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
END SUBROUTINE TimeStepByEulerStaticExpAdapTS
#endif

#if (PP_TimeDiscMethod==1000)
SUBROUTINE TimeStep_LD(t)
!===================================================================================================================================
! Low Diffusion Method (Mirza 2013)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,          ONLY:dt
#ifdef PARTICLES
USE MOD_Particle_Vars,          ONLY:PartState, LastPartPos, Time, PDM,PEM 
USE MOD_LD_Vars,                ONLY:LD_RHS
USE MOD_LD,                     ONLY:LD_main
USE MOD_part_tools,             ONLY:UpdateNextFreePosition
USE MOD_part_emission,          ONLY:ParticleInserting
USE MOD_Particle_Tracking_Vars, ONLY:ntracks,tTracking,tLocalization,MeasureTrackTime
USE MOD_Particle_Tracking,      ONLY:ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,           ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfParticles
#endif /*MPI*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: timeStart, timeEnd
!===================================================================================================================================

  Time = t
  CALL LD_main()
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)

  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + LD_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + LD_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + LD_RHS(1:PDM%ParticleVecLength,3)
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                         + PartState(1:PDM%ParticleVecLength,4) * dt
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                         + PartState(1:PDM%ParticleVecLength,5) * dt
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                         + PartState(1:PDM%ParticleVecLength,6) * dt


#ifdef MPI
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleRefTracking()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tTracking=tTracking+TimeEnd-TimeStart
  END IF
#ifdef MPI
  CALL SendNbOfParticles()
  CALL MPIParticleSend()
  CALL MPIParticleRecv()
#endif

  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleInserting()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tLocalization=tLocalization+TimeEnd-TimeStart
  END IF
  CALL UpdateNextFreePosition()

END SUBROUTINE TimeStep_LD
#endif

#if (PP_TimeDiscMethod==1001)
SUBROUTINE TimeStep_LD_DSMC(t)
!===================================================================================================================================
! Low Diffusion Method (Mirza 2013)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,    ONLY: dt, iter, TEnd
#ifdef PARTICLES
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, Time, PDM,PEM, WriteMacroValues
USE MOD_LD_Vars,          ONLY : LD_DSMC_RHS
USE MOD_LD,               ONLY : LD_main
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_LD_DSMC_TOOLS
!USE MOD_part_boundary,    ONLY : ParticleBoundary
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC_Vars,        ONLY : DSMC
USE MOD_LD_DSMC_DOMAIN_DEC
USE MOD_Particle_Tracking_Vars, ONLY:ntracks,tTracking,tLocalization,MeasureTrackTime
USE MOD_Particle_Tracking,      ONLY:ParticleRefTracking
#endif
#ifdef MPI
USE MOD_Particle_MPI,           ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfParticles
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: nOutput
  REAL              :: timeStart, timeEnd
!===================================================================================================================================

  IF(iter.EQ. 0.0) CALL LD_DSMC_DOMAIN_DECOMPOSITION

  Time = t
  LD_DSMC_RHS(1:PDM%ParticleVecLength,1) = 0
  LD_DSMC_RHS(1:PDM%ParticleVecLength,2) = 0
  LD_DSMC_RHS(1:PDM%ParticleVecLength,3) = 0
  CALL LD_DSMC_Indicate_DSMC_Particles
  CALL DSMC_main() ! first dsmc then ld due to RHS-calculation!
  CALL LD_main()
! ----- Start Analyze Particles
  IF (.NOT.WriteMacroValues) THEN
    IF(Time.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN
      CALL LD_DSMC_data_sampling()  ! Data sampling for output
      IF(DSMC%NumOutput.NE.0) THEN
        nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
        IF(Time.ge.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
          DSMC%NumOutput = DSMC%NumOutput - 1
          CALL LD_DSMC_output_calc()
        END IF
      END IF
    END IF
  END IF
! ----- End Analyze Particles
  CALL LD_DSMC_Clone_Particles
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)

  PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                         + LD_DSMC_RHS(1:PDM%ParticleVecLength,1)
  PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                         + LD_DSMC_RHS(1:PDM%ParticleVecLength,2)
  PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                         + LD_DSMC_RHS(1:PDM%ParticleVecLength,3)
  PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) &
                                         + PartState(1:PDM%ParticleVecLength,4) * dt
  PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) &
                                         + PartState(1:PDM%ParticleVecLength,5) * dt
  PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) &
                                         + PartState(1:PDM%ParticleVecLength,6) * dt
#ifdef MPI
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleRefTracking()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tTracking=tTracking+TimeEnd-TimeStart
  END IF
#ifdef MPI
  CALL SendNbOfParticles()
  CALL MPIParticleSend()
  CALL MPIParticleRecv()
#endif

  CALL LD_DSMC_Clean_Bufferregion

  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleInserting()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tLocalization=tLocalization+TimeEnd-TimeStart
  END IF
  CALL UpdateNextFreePosition()

END SUBROUTINE TimeStep_LD_DSMC
#endif

SUBROUTINE FillCFL_DFL()
!===================================================================================================================================
! scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
! Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" . For N=1-10, 
! input CFLscale can now be 1. and will be scaled adequately depending on PP_N, NodeType and TimeDisc method 
! DFLscale dependance not implemented jet.
!===================================================================================================================================

! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY:CFLScale
#if (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==5) || (PP_TimeDiscMethod==200) || (PP_TimeDiscMethod==201)
USE MOD_TimeDisc_Vars,ONLY:CFLScaleAlpha
#endif
#if (PP_TimeDiscMethod==6) || (PP_TimeDiscMethod==1)
USE MOD_TimeDisc_Vars,ONLY:CFLScaleAlpha
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! CFL in DG depends on the polynomial degree

#if (PP_TimeDiscMethod==3)
#  if (PP_NodeType==1)
  !Gauss  Taylor DG Timeorder=SpaceOrder
  CFLscale=CFLscale*0.55*(2*PP_N+1)/(PP_N+1)
#  elif (PP_NodeType==2)
  !Gauss-Lobatto  Taylor DG Timeorder=SpaceOrder, scales with N+1
  CFLscale=CFLscale*(2*PP_N+1)/(PP_N+1)
#  endif /*PP_NodeType*/
#endif /*PP_TimeDiscMethod*/

#if (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==5) || (PP_TimeDiscMethod==200)||(PP_TimeDiscMethod==201)||(PP_TimeDiscMethod==1)
IF(PP_N.GT.10) CALL abort(&
  __STAMP__,&
  'Polynomial degree is to high!',PP_N,999.)
CFLScale=CFLScale*CFLScaleAlpha(PP_N)
#endif
#if (PP_TimeDiscMethod==6)
IF(PP_N.GT.15) CALL abort(&
  __STAMP__,&
  'Polynomial degree is to high!',PP_N,999.)
CFLScale=CFLScale*CFLScaleAlpha(PP_N)
#endif
!scale with 2N+1
CFLScale = CFLScale/(2.*PP_N+1.)
SWRITE(UNIT_stdOut,'(A,ES16.7)') '   CFL:',CFLScale
END SUBROUTINE fillCFL_DFL



SUBROUTINE FinalizeTimeDisc()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars,ONLY:TimeDiscInitIsDone
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
TimeDiscInitIsDone = .FALSE.
END SUBROUTINE FinalizeTimeDisc


END MODULE MOD_TimeDisc
