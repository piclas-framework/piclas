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

#ifdef IMEX
#ifndef maxwell
CALL abort(__STAMP__,&
  ' Preconditioner is only implemented for Maxwell!')
#endif /*maxwell*/
#endif /*IMEX*/


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
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ERK3-Particles and ESDIRK3-Field'
#elif (PP_TimeDiscMethod==102)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ERK4-Particles and ESDIRK4-Field'
#elif (PP_TimeDiscMethod==103)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ERK3-Particles and ESDIRK3-Field'
  SWRITE(UNIT_stdOut,'(A)') '                             Ascher - 2 - 3 -3               '
#elif (PP_TimeDiscMethod==104)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Implicit with Newton'
#elif (PP_TimeDiscMethod==110)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler-Implicit-Particles and Euler-Explicit-Field'
#elif (PP_TimeDiscMethod==111)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ESDIRK3-Particles and ERK3-Field'
#elif (PP_TimeDiscMethod==112)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ESDIRK4-Particles and ERK4-Field'
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
USE MOD_TimeDisc_Vars,         ONLY: time,TEnd,dt,tAnalyze,iter &
                                     ,IterDisplayStep,DoDisplayIter,sdtCFLOne,CFLtoOne
USE MOD_Restart_Vars,          ONLY: DoRestart,RestartTime
USE MOD_CalcTimeStep,          ONLY: CalcTimeStep
USE MOD_Analyze,               ONLY: CalcError,PerformAnalyze
USE MOD_Analyze_Vars,          ONLY: Analyze_dt
#ifdef PARTICLES
USE MOD_Particle_Analyze,      ONLY: AnalyzeParticles
#else
USE MOD_AnalyzeField,          ONLY: AnalyzeField
#endif /*PARTICLES*/
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
USE MOD_PARTICLE_Vars,         ONLY: ManualTimeStep, useManualTimestep, WriteMacroValues, MacroValSampTime
USE MOD_Particle_Tracking_vars, ONLY: tTracking,tLocalization,nTracks,MeasureTrackTime
#if (PP_TimeDiscMethod==201||PP_TimeDiscMethod==200)
USE MOD_PARTICLE_Vars,         ONLY: dt_maxwell,dt_max_particles,dt_part_ratio,MaxwellIterNum,NextTimeStepAdjustmentIter
USE MOD_Particle_Mesh_Vars,    ONLY: Geo
USE MOD_Equation_Vars,         ONLY: c
#endif /*(PP_TimeDiscMethod==201||PP_TimeDiscMethod==200)*/
#if (PP_TimeDiscMethod==201)
USE MOD_PARTICLE_Vars,         ONLY: PDM,Pt,PartState
#endif /*(PP_TimeDiscMethod==201)*/
USE MOD_PARTICLE_Vars,         ONLY : doParticleMerge, enableParticleMerge, vMPFMergeParticleIter
USE MOD_ReadInTools
USE MOD_DSMC_Vars,             ONLY: Iter_macvalout
#ifdef MPI
USE MOD_Particle_MPI,          ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#endif /*PARTICLES*/
#ifdef PP_POIS
USE MOD_Equation,              ONLY: EvalGradient
#endif /*PP_POIS*/
USE MOD_LoadBalance_Vars,      ONLY: nSkipAnalyze
#ifdef MPI
USE MOD_LoadBalance,           ONLY: LoadBalance,LoadMeasure,ComputeParticleWeightAndLoad,ComputeElemLoad
USE MOD_LoadBalance_Vars,      ONLY: DoLoadBalance
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: tFuture,tZero
INTEGER(KIND=8)              :: nAnalyze
INTEGER                      :: iAnalyze
REAL                         :: dt_Min, tEndDiff, tAnalyzeDiff
#if (PP_TimeDiscMethod==201)
REAL                         :: dt_temp,vMax,vMaxx,vMaxy,vMaxz
#endif
INTEGER(KIND=8)              :: iter_loc
REAL                         :: CalcTimeStart,CalcTimeEnd
INTEGER                      :: TimeArray(8)              ! Array for system time
REAL                         :: CurrentImbalance
LOGICAL                      :: PerformLoadBalance
#if (PP_TimeDiscMethod==201)
INTEGER                      :: MaximumIterNum, iPart
LOGICAL                      :: NoPartInside
#endif
LOGICAL                      :: finalIter
!===================================================================================================================================
! init
SWRITE(UNIT_StdOut,'(132("-"))')
#ifdef PARTICLES
iter_macvalout=0
IF (WriteMacroValues) MacroValSampTime = Time
#endif /*PARTICLES*/
IF(.NOT.DoRestart)THEN
  time=0.
  SWRITE(UNIT_StdOut,*)'INITIAL PROJECTION:'
ELSE
  time=RestartTime
  SWRITE(UNIT_StdOut,*)'REWRITING SOLUTION:'
END IF
tZero=time
nAnalyze=1
iAnalyze=1
tAnalyze=MIN(tZero+Analyze_dt,tEnd)

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems*(PP_N+1)**3)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)

! Determine next analyze time, since it will be written into output file
tFuture=MIN(time+Analyze_dt,tEnd)
!Evaluate Gradients to get Potential in case of Restart and Poisson Calc
#ifdef PP_POIS
IF(DoRestart) CALL EvalGradient()
#endif /*PP_POIS*/
! Write the state at time=0, i.e. the initial condition
CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture)

! if measurement of particle tracking time
#ifdef PARTICLES
IF(MeasureTrackTime)THEN
  nTracks=0
  tTracking=0
  tLocalization=0
END IF
#endif /*PARTICLES*/

! No computation needed if tEnd=tStart!
IF(time.EQ.tEnd)RETURN

#if defined(PARTICLES) && defined(MPI)
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
#endif /*MPI PARTICLES*/

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
  sdtCFLOne  = 1.0/(dt_maxwell*CFLtoOne)

  NextTimeStepAdjustmentIter = 0
  MaxwellIterNum = INT(MAX(GEO%xmaxglob-GEO%xminglob,GEO%ymaxglob-GEO%yminglob,GEO%zmaxglob-GEO%zminglob) &
                 / (c * dt_maxwell))
  IF (MaxwellIterNum*dt_maxwell.GT.dt_max_particles) THEN
    WRITE(*,*) 'WARNING: Time of Maxwell Solver is greater then particle time step!'
  END IF
  IterDisplayStep = MAX(INT(IterDisplayStepUser/(dt_max_particles / dt_maxwell)),1) !IterDisplayStepUser refers to dt_maxwell
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
  sdtCFLOne  = 1.0/(dt_Min*CFLtoOne)
  dt=dt_Min
  ! calculate time step for sub-cycling of divergence correction
  ! automatic particle time step of quasi-stationary time integration is not implemented
#ifdef PARTICLES
#if (PP_TimeDiscMethod==200)
  ! this will not work if particles have velocity of zero
   CALL abort(&
   __STAMP__&
   ,' Error in static computations, a maximum delta t (=ManualTimeStep) needs to be set!')
#endif
END IF ! useManualTimestep

#if (PP_TimeDiscMethod==201)
dt_maxwell = CALCTIMESTEP()
sdtCFLOne  = 1.0/(dt_Maxwell*CFLtoOne)
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
CALL PerformAnalyze(time,iter,0.,forceAnalyze=.TRUE.,OutPut=.FALSE.)

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
    NoPartInside=.TRUE.
    DO 
      vMaxx = 0.
      vMaxy = 0.
      vMaxz = 0.
      DO iPart=1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(iPart)) THEN
          vMaxx = MAX( vMaxx , ABS(PartState(iPart, 4) + dt_temp*Pt(iPart,1)) )
          vMaxy = MAX( vMaxy , ABS(PartState(iPart, 5) + dt_temp*Pt(iPart,2)) )
          vMaxz = MAX( vMaxz , ABS(PartState(iPart, 6) + dt_temp*Pt(iPart,3)) )
          NoPartInside=.FALSE. 
        END IF
      END DO
!! -- intrinsic logical->real/int-conversion should be avoided!!!
!      vMaxx = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 4) + dt_temp*Pt(1:PDM%ParticleVecLength,1))))
!      vMaxy = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 5) + dt_temp*Pt(1:PDM%ParticleVecLength,2))))
!      vMaxz = MAXVAL(ABS(PDM%ParticleInside(1:PDM%ParticleVecLength) &
!            * (PartState(1:PDM%ParticleVecLength, 6) + dt_temp*Pt(1:PDM%ParticleVecLength,3))))
!      vMax = MAX( SQRT(vMaxx*vMaxx + vMaxy*vMaxy + vMaxz*vMaxz) , 1.0 ) 
      vMax = MAX(vMaxx,vMaxy,vMaxz,1.0) 
#ifdef MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vMax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,NoPartInside,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,iError)
#endif /*MPI*/

      IF (NoPartInside) THEN
        dt_max_particles = dt_maxwell
        EXIT
      ELSE
        dt_max_particles =  max(dt_part_ratio*dt_maxwell*c/(vMax),dt_maxwell)
      END IF
      dt_temp = (dt_max_particles+dt_temp)/2
      IF((dt_temp.GE.dt_max_particles*0.95).AND.(dt_temp.LE.dt_max_particles*1.05)) EXIT
    END DO
  END IF
  dt_Min = dt_max_particles
  MaxwellIterNum = INT(dt_max_particles / dt_maxwell)
  IterDisplayStep = MAX(INT(IterDisplayStepUser/(dt_max_particles / dt_maxwell)),1) !IterDisplayStepUser refers to dt_maxwell
  IF (MaxwellIterNum.GT.MaximumIterNum) MaxwellIterNum = MaximumIterNum
  IF ((MPIroot).AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
    SWRITE(UNIT_StdOut,'(132("!"))')
    SWRITE(UNIT_StdOut,*)  'New IterNum for MaxwellSolver: ', MaxwellIterNum 
    SWRITE(UNIT_StdOut,*)  'New Particle TimeStep: ', dt_max_particles
  END IF
#endif

#ifdef PARTICLES
  IF(enableParticleMerge) THEN
    IF ((iter.GT.0).AND.(MOD(iter,vMPFMergeParticleIter).EQ.0)) doParticleMerge=.true.
  END IF
#endif /*PARTICLES*/

  tAnalyzeDiff=tAnalyze-time    ! time to next analysis, put in extra variable so number does not change due to numerical errors
  tEndDiff=tEnd-time            ! dito for end time
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
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==2)
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==3)
  CALL TimeStepByTAYLOR(time)
#elif (PP_TimeDiscMethod==4)
  CALL TimeStep_DSMC(time)
#elif (PP_TimeDiscMethod==5)
  CALL TimeStepByRK4EulerExpl(time)
#elif (PP_TimeDiscMethod==6)
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==42)
  CALL TimeStep_DSMC_Debug(time) ! Reservoir and Debug
#elif (PP_TimeDiscMethod==100)
  CALL TimeStepByEulerImplicit(time) ! O1 Euler Implicit
#elif (PP_TimeDiscMethod==101)
  CALL TimeStepByIMEXRK(time) ! ) O3 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==102)
  CALL TimeStepByIMEXRK(time) ! O4 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==103)
  CALL TimeStepByIMEXRK(time) ! ) O3 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==104)
  CALL TimeStepByEulerNewton(time) ! O1 Euler Implicit
#elif (PP_TimeDiscMethod==110)
  CALL TimeStepByEulerImplicitParticle(time) ! ) O3 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==111)
  CALL TimeStepByIMPA(time) ! ) O3 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==112)
  CALL TimeStepByIMPA(time) ! O4 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==200)
  CALL TimeStepByEulerStaticExp(time) ! O1 Euler Static Explicit
#elif (PP_TimeDiscMethod==201)
  CALL TimeStepByEulerStaticExpAdapTS(time) ! O1 Euler Static Explicit with adaptive TimeStep
#elif (PP_TimeDiscMethod==1000)
  CALL TimeStep_LD(time)
#elif (PP_TimeDiscMethod==1001)
  CALL TimeStep_LD_DSMC(time)
#endif
#ifdef MPI
  CALL LoadMeasure()
#endif /*MPI*/
  iter=iter+1
  iter_loc=iter_loc+1
  time=time+dt
  IF(MPIroot) THEN
    IF(DoDisplayIter)THEN
      IF(MOD(iter,IterDisplayStep).EQ.0) THEN
         SWRITE(*,*) "iter:", iter,"time:",time
      END IF
    END IF
  END IF
#if (PP_TimeDiscMethod!=1) &&  (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6)
  ! calling the analyze routines
  CALL PerformAnalyze(time,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)
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
#ifdef MPI
    !CALL ComputeParticleWeightAndLoad(CurrentImbalance,PerformLoadBalance)
    CALL ComputeElemLoad(CurrentImbalance,PerformLoadBalance)
#endif /*MPI*/
    ! future time
    nAnalyze=nAnalyze+1
    tFuture=tZero+REAL(nAnalyze)*Analyze_dt
#ifdef MPI
    IF(iAnalyze.EQ.nSkipAnalyze .OR. PerformLoadBalance .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#else
    IF( iAnalyze.EQ.nSkipAnalyze .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#endif /*MPI*/
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
      CALL PerformAnalyze(time,iter,tenddiff,forceAnalyze=.TRUE.,OutPut=.TRUE.,LastIter=finalIter)
      ! Write state to file
      IF(DoPML) CALL BacktransformPMLVars()
      CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture)
      IF(DoPML) CALL TransformPMLVars()
      ! Write recordpoints data to hdf5
      IF(RP_onProc) CALL WriteRPtoHDF5(tAnalyze,.TRUE.)
!#ifdef MPI
      IF(iAnalyze.EQ.nSkipAnalyze) iAnalyze=0
!#endif /*MPI*/
      SWRITE(UNIT_StdOut,'(132("-"))')
    END IF
    iAnalyze=iAnalyze+1
    iter_loc=0
    tAnalyze=tZero+REAL(nAnalyze)*Analyze_dt
    IF (tAnalyze > tEnd) tAnalyze = tEnd
#ifdef MPI
    IF(DoLoadBalance .AND. time.LT.tEnd)THEN
      RestartTime=time
      CALL LoadBalance(CurrentImbalance,PerformLoadBalance)
      dt_Min=CALCTIMESTEP()
      dt=dt_Min
      CALL PerformAnalyze(time,iter,tendDiff,forceAnalyze=.TRUE.,OutPut=.FALSE.)
      ! CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture) ! not sure if required
    END IF
#endif /*MPI*/
    CalcTimeStart=BOLTZPLATZTIME()
  ENDIF   
  IF(time.GE.tEnd) THEN  ! done, worst case: one additional time step
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
USE MOD_DG_Vars,                 ONLY: U,Ut!,nTotalU
USE MOD_PML_Vars,                ONLY: U2,U2t,nPMLElems,DoPML
USE MOD_PML,                     ONLY: PMLTimeDerivative,CalcPMLSource
USE MOD_Filter,                  ONLY: Filter
USE MOD_Equation,                ONLY: DivCleaningDamping
USE MOD_Equation,                ONLY: CalcSource
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
#ifdef MPI
USE MOD_LoadBalance_Vars,        ONLY: tCurrent
#endif /*MPI*/
#ifdef PP_POIS
USE MOD_Equation,                ONLY: DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,           ONLY: Phi,Phit,nTotalPhi
#endif /*PP_POIS*/
#ifdef PARTICLES
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_PICDepo,                 ONLY: Deposition!, DepositionMPF
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime, PEM, PDM, & 
                                        doParticleMerge,PartPressureCell
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,           ONLY: StartParticleMerge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#ifdef MPI
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
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
#endif /*PP_POIS*/
#ifdef MPI
! load balance
REAL                          :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

! RK coefficients
DO iStage=1,nRKStages
  b_dt(iStage)=RK_b(iStage)*dt
END DO
iStage=1

#ifdef PARTICLES
#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
IF (t.GE.DelayTime) THEN
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleInserting()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tLocalization=tLocalization+TimeEnd-TimeStart
  END IF
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(4)=tCurrent(4)+tLBEnd-tLBStart
#endif /*MPI*/

#ifdef MPI
CALL CountPartsPerElem()
#endif /*MPI*/


#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
IF (t.GE.DelayTime) THEN
  ! forces on particle
  ! can be used to hide sending of number of particles
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  CALL CalcPartRHS()
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(6)=tCurrent(6)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(7)=tCurrent(7)+tLBEnd-tLBStart
#endif /*MPI*/
#endif /*PARTICLES*/

! field solver
! time measurement in weakForm
CALL DGTimeDerivative_weakForm(t,t,0,doSource=.TRUE.)
CALL DivCleaningDamping()

#ifdef PP_POIS
! Potential
CALL DGTimeDerivative_weakForm_Pois(t,t,0)
CALL DivCleaningDamping_Pois()
#endif /*PP_POIS*/

#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
IF(DoPML) THEN
  CALL CalcPMLSource()
  CALL PMLTimeDerivative()
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(3)=tCurrent(3)+tLBEnd-tLBStart
#endif /*MPI*/

! calling the analyze routines
CALL PerformAnalyze(t,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)

! first RK step
#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/

! EM field
Ut_temp = Ut 
U = U + Ut*b_dt(1)

#ifdef PP_POIS
Phit_temp = Phit 
Phi = Phi + Phit*b_dt(1)
CALL EvalGradient()
#endif /*PP_POIS*/

#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
#endif /*MPI*/

#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
!PML auxiliary variables
IF(DoPML) THEN
  U2t_temp = U2t
  U2 = U2 + U2t*b_dt(1)
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(3)=tCurrent(3)+tLBEnd-tLBStart
#endif /*MPI*/


#ifdef PARTICLES
! particles
#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
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
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(9)=tCurrent(9)+tLBEnd-tLBStart
#endif /*MPI*/
END IF

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
  CALL IRecvNbofParticles()
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
  tLBStart = LOCALTIME() ! LB Time Start
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
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(5)=tCurrent(5)+tLBEnd-tLBStart
  tLBStart = LOCALTIME() ! LB Time Start
  CALL SendNbOfParticles()
  !CALL MPIParticleSend()
  !CALL MPIParticleRecv()
  !PartMPIExchange%nMPIParticles=0
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
#endif /*MPI*/
  !CALL Filter(U)
END IF
#endif /*PARTICLES*/

DO iStage=2,nRKStages
  tStage=t+dt*RK_c(iStage)
#ifdef PARTICLES
  ! deposition  
  IF (t.GE.DelayTime) THEN 
#ifdef MPI
     tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
     CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
#ifdef MPI
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(6)=tCurrent(6)+tLBEnd-tLBStart
     tLBStart = LOCALTIME() ! LB Time Start
     CALL MPIParticleSend()
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
     tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
!    ! deposition  
     CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(7)=tCurrent(7)+tLBEnd-tLBStart
     tLBStart = LOCALTIME() ! LB Time Start
     CALL MPIParticleRecv()
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
     ! second buffer
     tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
     CALL InterpolateFieldToParticle(doInnerParts=.FALSE.)
     CALL CalcPartRHS()
#ifdef MPI
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(6)=tCurrent(6)+tLBEnd-tLBStart
     tLBStart = LOCALTIME() ! LB Time Start
     ! null here, careful
#endif /*MPI*/

     CALL Deposition(doInnerParts=.FALSE.)
     IF(DoVerifyCharge) CALL VerifyDepositedCharge()
#ifdef MPI
     tLBEnd = LOCALTIME() ! LB Time End
     tCurrent(7)=tCurrent(7)+tLBEnd-tLBStart
     ! null here, careful
     PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
!    IF (usevMPF) THEN 
!      CALL !DepositionMPF()
!    ELSE 
!      CALL Deposition()
!    END IF
  END IF
#ifdef MPI
  CALL CountPartsPerElem()
#endif /*MPI*/
#endif /*PARTICLES*/

  ! field solver
  CALL DGTimeDerivative_weakForm(t,tStage,0,doSource=.TRUE.)
  CALL DivCleaningDamping()
#ifdef PP_POIS
  CALL DGTimeDerivative_weakForm_Pois(t,tStage,0)
  CALL DivCleaningDamping_Pois()
#endif

#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  IF(DoPML) THEN
    CALL CalcPMLSource()
    CALL PMLTimeDerivative()
  END IF
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(3)=tCurrent(3)+tLBEnd-tLBStart
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/

  ! performe RK steps
  ! field step
  Ut_temp = Ut - Ut_temp*RK_a(iStage)
  U = U + Ut_temp*b_dt(iStage)
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
#endif /*MPI*/

#ifdef PP_POIS
  Phit_temp = Phit - Phit_temp*RK_a(iStage)
  Phi = Phi + Phit_temp*b_dt(iStage)
  CALL EvalGradient()
#endif

  !PML auxiliary variables
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  IF(DoPML)THEN
    U2t_temp = U2t - U2t_temp*RK_a(iStage)
    U2 = U2 + U2t_temp*b_dt(iStage)
  END IF
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(3)=tCurrent(3)+tLBEnd-tLBStart
#endif /*MPI*/

#ifdef PARTICLES
  ! particle step
  IF (t.GE.DelayTime) THEN
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
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
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(9)=tCurrent(9)+tLBEnd-tLBStart
    ! particle tracking
    tLBStart = LOCALTIME() ! LB Time Start
    CALL IRecvNbofParticles()
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
    ! actual tracking
    tLBStart = LOCALTIME() ! LB Time Start
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
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(5)=tCurrent(5)+tLBEnd-tLBStart
    tLBStart = LOCALTIME() ! LB Time Start
    CALL SendNbOfParticles()
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
    !CALL MPIParticleSend()
    !CALL MPIParticleRecv()
    !PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  END IF
#endif /*PARTICLES*/

END DO

#ifdef PARTICLES
#ifdef MPI
IF (t.GE.DelayTime) THEN
  tLBStart = LOCALTIME() ! LB Time Start
  CALL MPIParticleSend()
  CALL MPIParticleRecv()
  PartMPIExchange%nMPIParticles=0
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
END IF 
#endif /*MPI*/

IF (doParticleMerge) THEN
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(11)=tCurrent(11)+tLBEnd-tLBStart
#endif /*MPI*/
END IF

#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(12)=tCurrent(12)+tLBEnd-tLBStart
#endif /*MPI*/

IF (doParticleMerge) THEN
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  CALL StartParticleMerge()  
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )
  END IF
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(11)=tCurrent(11)+tLBEnd-tLBStart
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  CALL UpdateNextFreePosition()
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(12)=tCurrent(12)+tLBEnd-tLBStart
#endif /*MPI*/
END IF

IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,3)
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(8)=tCurrent(8)+tLBEnd-tLBStart
#endif /*MPI*/
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
  CALL DGTimeDerivative_weakForm(t,t,tIndex,doSource=.TRUE.)
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
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos,  PDM,PEM
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
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS
USE MOD_part_MPFtools,    ONLY : StartParticleMerge
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
CALL DGTimeDerivative_weakForm(t,t,0,doSource=.TRUE.)
CALL DivCleaningDamping()
Ut_temp = Ut 
U = U + Ut*b_dt(1)

DO rk=2,5
  tStage=t+dt*RK_c(rk)
  ! field RHS
  CALL DGTimeDerivative_weakForm(t,tStage,0,doSource=.TRUE.)
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
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, PDM,PEM, Species, PartSpecies
USE MOD_DSMC_Vars,        ONLY : DSMC_RHS, DSMC, Debug_Energy,PartStateIntEn
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_Particle_Tracking_vars, ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,ONLY: ParticleTrackingCurved,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,     ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,ONLY: PartMPIExchange
#endif /*MPI*/
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
REAL                  :: timeStart, timeEnd
!===================================================================================================================================

IF (DSMC%ReservoirSimu) THEN ! fix grid should be defined for reservoir simu
  CALL UpdateNextFreePosition()

!  Debug_Energy=0.0
!  DO i=1,PDM%ParticleVecLength
!    IF (PDM%ParticleInside(i)) THEN
!      Debug_Energy(1)  = Debug_Energy(1) +&
!        0.5* Species(PartSpecies(i))%MassIC*(PartState(i,4)**2+PartState(i,5)**2+PartState(i,6)**2)&
!      + PartStateIntEn(i,3)
!    END IF
!  END DO
  CALL DSMC_main()
!  DO i=1,PDM%ParticleVecLength
!    IF (PDM%ParticleInside(i)) THEN
!      Debug_Energy(2)  = Debug_Energy(2) +&
!        0.5* Species(PartSpecies(i))%MassIC*(PartState(i,4)**2+PartState(i,5)**2+PartState(i,6)**2)&
!      + PartStateIntEn(i,3)
!    END IF
!  END DO



!  IF (Debug_Energy(1)-Debug_Energy(2)>0.0)THEN
!    print*,"energy loss"
!    read*
!  else 
!   print*,Debug_Energy(1),Debug_Energy(2),"   Difference(1-2)=",Debug_Energy(1)-Debug_Energy(2)
!  END IF
!  Debug_Energy=0.0
!  read*

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
USE MOD_DG_Vars,          ONLY:U,Ut
USE MOD_DG,               ONLY:DGTimeDerivative_weakForm
USE MOD_TimeDisc_Vars,    ONLY:dt,time
USE MOD_LinearSolver,     ONLY : LinearSolver
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation, ONLY : InterpolateFieldToParticle
USE MOD_PIC_Vars,         ONLY : PIC
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, DelayTime, PEM, PDM, usevMPF
USE MOD_part_RHS,         ONLY : CalcPartRHS
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS, DSMC
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
REAL               :: tstage,coeff
!===================================================================================================================================

! one Euler implicit step
! time for source is t + dt 
tstage = t + dt
coeff  = dt*1.

#ifdef maxwell
IF(PrecondType.GT.0)THEN
!  IF (iter==0) CALL BuildJacDG()
  IF (iter==0) CALL BuildPrecond(t,t,0,1.,dt)
END IF
#endif /*maxwell*/

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


! EM field
! U predict
!U = U
! b
LinSolverRHS = U
ImplicitSource=0.
CALL LinearSolver(tstage,coeff)
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

#if (PP_TimeDiscMethod==101) || (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==103)
SUBROUTINE TimeStepByIMEXRK(t)
!===================================================================================================================================
! IMEX time integrator
! ERK4 for particles
! ESDIRK4 for field
! from: Kennedy & Carpenter 2003
! Additive Runge-Kutta schemes for convection-diffusion-reaction equations
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DG_Vars,                 ONLY:U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,iter,iStage, nRKStages,time
USE MOD_TimeDisc_Vars,           ONLY: ERK_a,ESDIRK_a,RK_b,RK_c,RK_bs
USE MOD_DG_Vars,                 ONLY: U,Ut
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_LinearSolver,            ONLY: LinearSolver
USE MOD_Predictor,               ONLY: Predictor,StorePredictor
USE MOD_LinearSolver_Vars,       ONLY: ImplicitSource,LinSolverRHS
USE MOD_Equation,                ONLY: DivCleaningDamping
#ifdef maxwell
USE MOD_Precond,                 ONLY: BuildPrecond
USE MOD_Precond_Vars,            ONLY: PrecondType
#endif /*maxwell*/
USE MOD_JacDG,                   ONLY: BuildJacDG
USE MOD_Equation,                ONLY: CalcSource
#ifdef PARTICLES
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_PIC_Vars,                ONLY: PIC
USE MOD_Particle_Vars,           ONLY: PartStateN,PartStage
USE MOD_Particle_Vars,           ONLY: PartState, Pt, LastPartPos, DelayTime, PEM, PDM, usevMPF
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS, DSMC
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_PIC_Analyze,             ONLY: CalcDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
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
REAL               :: tstage
! implicit 
REAL               :: alpha
REAL               :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL               :: FieldStage (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:nRKStages-1)
REAL               :: FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:nRKStages-1)
REAL               :: tRatio, tphi
INTEGER            :: iCounter !, iStage
! explicit
!===================================================================================================================================

#ifdef maxwell
! caution hard coded
!IF (iter==0) CALL BuildPrecond(t,t,0,0.25,dt)
IF (iter==0) CALL BuildPrecond(t,t,0,RK_b(nRKStages),dt)
!IF(PrecondType.GT.0)THEN
  !IF (iter==0) CALL BuildJacDG()
!END IF
#endif /*maxwell*/

Un = U
FieldSource=0.
! prediction with embadded scheme

!tRatio= ! dt^n+1 / dt^n
tRatio = 1. 

#ifdef PARTICLES
! Partilce locating
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL ParticleInserting()
END IF
#endif

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------
tStage=t
! explicit
#ifdef PARTICLES
PartStateN(1:PDM%ParticleVecLength,1:6)=PartState(1:PDM%ParticleVecLength,1:6)
! ----------------------------------------------------------------------------------------------------------------------------------
! implicit - Maxwell solver
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
! u1=un
#endif
ImplicitSource=0.
CALL CalcSource(tStage,1.,ImplicitSource)
! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to 6
! ----------------------------------------------------------------------------------------------------------------------------------
DO iStage=2,nRKStages
  ! compute f(u^s-1)
  CALL DGTimeDerivative_weakForm(tStage, tStage, 0,doSource=.FALSE.)
  ! store old values for use in next stages
  FieldStage (:,:,:,:,:,iStage-1) = Ut
  FieldSource(:,:,:,:,:,iStage-1) = ImplicitSource

  ! time of current stage
  tStage = t + RK_c(iStage)*dt
  ! store predictor
  CALL StorePredictor()

  !--------------------------------------------------------------------------------------------------------------------------------
  ! explicit - Particle pusher
  !--------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
  ! particle RHS
  IF (t.GE.DelayTime) THEN
    CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
    CALL CalcPartRHS()

    PartStage(1:PDM%ParticleVecLength,1:3,iStage-1) = PartState(1:PDM%ParticleVecLength,4:6)
    PartStage(1:PDM%ParticleVecLength,4:6,iStage-1) = Pt       (1:PDM%ParticleVecLength,1:3)
  END IF
  ! particle step
  IF (t.GE.DelayTime) THEN
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    PartState(1:PDM%ParticleVecLength,1:6)=PartStateN(1:PDM%ParticleVecLength,1:6)
    DO iCounter=1,iStage-1
      PartState(1:PDM%ParticleVecLength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)   &
                                                + ERK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,1:6,iCounter)
    END DO ! counter
  END IF

  IF (t.GE.DelayTime) THEN
#ifdef MPI
  ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
#endif /*MPI*/
    IF(DoRefMapping)THEN
      CALL ParticleRefTracking()
    ELSE
      CALL ParticleTrackingCurved()
    END IF
#ifdef MPI
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
#endif
  END IF


  !--------------------------------------------------------------------------------------------------------------------------------
  ! implicit - Maxwell solver
  !--------------------------------------------------------------------------------------------------------------------------------
  IF (t.GE.DelayTime) THEN
    ! because of emmision and UpdateParticlePosition
    CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
    ! here: finish deposition with delta kernal
    !       maps source terms in physical space
    ! ALWAYS require
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    CALL Deposition(doInnerParts=.FALSE.)
    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
  END IF
#endif /*PARTICLES*/
  !--------------------------------------------------------------------------------------------------------------------------------
  ! compute RHS for linear solver
  LinSolverRHS=Un
  DO iCounter = 1,iStage-1
    LinSolverRHS = LinSolverRHS + ESDIRK_a(iStage,iCounter)*dt*(FieldStage(:,:,:,:,:,iCounter)+FieldSource(:,:,:,:,:,iCounter))
  END DO
  ! get predictor of u^s+1
  CALL Predictor(iStage,dt,Un,FieldSource,FieldStage)

  ImplicitSource=0.
  alpha = ESDIRK_a(iStage,iStage)*dt
  ! solve to new stage 
  CALL LinearSolver(tstage,alpha)
    ! damping
  !CALL DivCleaningDamping()
END DO

!----------------------------------------------------------------------------------------------------------------------------------
! update to next time level
!----------------------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
! only required for explicit equations // particle pusher
! implicit coefficients are equal, therefore not required, however, would result in smoother solution
! particle RHS
IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  CALL CalcPartRHS()
!  PartStage(1:PDM%ParticleVecLength,1:3,6) = PartState(1:PDM%ParticleVecLength,4:6)
!  PartStage(1:PDM%ParticleVecLength,4:6,6) = Pt       (1:PDM%ParticleVecLength,1:3)
END IF


IF (t.GE.DelayTime) THEN
  LastPartPos(1:PDM%ParticleVecLength,1)  =PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)  =PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)  =PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  ! compined
  ! stage 6
  PartState(1:PDM%ParticleVeclength,1:3) = PartStateN(1:PDM%ParticleVecLength,1:3)                   &
                                           + RK_b(nRKStages)*dt*PartState(1:PDM%ParticleVecLength,4:6)
  PartState(1:PDM%ParticleVeclength,4:6) = PartStateN(1:PDM%ParticleVecLength,4:6)                   &
                                           + RK_b(nRKSTages)*dt*Pt(1:PDM%ParticleVecLength,1:3)
!  ! stage 1 to 5
!  PartState(1:PDM%ParticleVeclength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)                    &
!                                           + RK_b(1)*dt*PartStage(1:PDM%ParticleVecLength,1:6,1)    &
!                                           + RK_b(2)*dt*PartStage(1:PDM%ParticleVecLength,1:6,2)    &
!                                           + RK_b(3)*dt*PartStage(1:PDM%ParticleVecLength,1:6,3)    &
!                                           + RK_b(4)*dt*PartStage(1:PDM%ParticleVecLength,1:6,4)    &
!                                           + RK_b(5)*dt*PartStage(1:PDM%ParticleVecLength,1:6,5) 
!  PartState(1:PDM%ParticleVecLength,1:6)=PartStateN(1:PDM%ParticleVecLength,1:6)
  DO iCounter=1,nRKStages-1
    PartState(1:PDM%ParticleVecLength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)   &
                                              + RK_b(iCounter)*dt*PartStage(1:PDM%ParticleVecLength,1:6,iCounter)
  END DO ! counter
END IF

! null iStage, because the old f(u^s) have not to be communicated
iStage = 0
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTrackingCurved()
  END IF
#ifdef MPI
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
#endif /*MPI*/
END IF

!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
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
#endif /*PARTICLES*/

! compute source and recompute field

END SUBROUTINE TimeStepByIMEXRK
#endif /*PP_TimeDiscMethod==102 || PP_TimeDiscMethod==101 || PP_TimeDiscMethod==103 */

#if (PP_TimeDiscMethod==104) 
SUBROUTINE TimeStepByEulerNewton(t)
!===================================================================================================================================
! Euler Implicit method:
! U^n+1 = U^n + dt*R(U^n+1)
! (I -dt*R)*U^n+1 = U^n
! Solve Linear System 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_TimeDisc_Vars,    ONLY: dt,iter
USE MOD_DG_Vars,          ONLY: U,Ut
USE MOD_DG,               ONLY: DGTimeDerivative_weakForm
USE MOD_LinearSolver,         ONLY: LinearSolver
USE MOD_LinearSolver_Vars,    ONLY: ImplicitSource,LinSolverRHS
USE MOD_Equation,         ONLY: DivCleaningDamping
#ifdef maxwell
USE MOD_Precond,          ONLY: BuildPrecond
USE MOD_Precond_Vars,     ONLY: PrecondType
#endif /*maxwell*/
!USE MOD_JacDG,            ONLY: BuildJacDG
!USE MOD_Newton,           ONLY: Newton
#ifdef PARTICLES
USE MOD_PICDepo,          ONLY: Deposition, DepositionMPF
USE MOD_PICInterpolation, ONLY: InterpolateFieldToParticle
USE MOD_PIC_Vars,         ONLY: PIC
USE MOD_Particle_Vars,    ONLY: PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, usevMPF
USE MOD_part_RHS,         ONLY: CalcPartRHS
USE MOD_part_emission,    ONLY: ParticleInserting
USE MOD_DSMC,             ONLY: DSMC_main
USE MOD_DSMC_Vars,        ONLY: useDSMC, DSMC_RHS, DSMC
USE MOD_PIC_Analyze,      ONLY: CalcDepositedCharge
USE MOD_part_tools,       ONLY: UpdateNextFreePosition
#endif
USE MOD_LinearSolver_Vars,    ONLY: nNewton,totalIterLinearSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: tstage,coeff
!===================================================================================================================================

Time = t
! one Euler implicit step
! time for source is t + dt 
tstage = t + dt
coeff  = dt*1.

#ifdef maxwell
IF(PrecondType.GT.0)THEN
!  IF (iter==0) CALL BuildJacDG()
  IF (iter==0) CALL BuildPrecond(t,t,0,1.,dt)
END IF
#endif /*maxwell*/

CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  IF (usevMPF) THEN 
    CALL DepositionMPF()
  ELSE 
    CALL Deposition()
  END IF
  !CALL CalcDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
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
CALL ParticleBoundary()
#ifdef MPI
CALL Communicate_PIC()
!CALL UpdateNextFreePosition() ! only required for parallel communication
#endif
! EM field
! U predict
!U = U
! b
LinSolverRHS = U ! equals Q in FLEXI
ImplicitSource=0.
CALL Newton(t,1.,1.)

!CALL LinearSolver(tstage,coeff)
!CALL DivCleaningDamping()

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

!SWRITE(UNIT_stdOut,'(A36,ES16.7,ES16.7)')'Global Iterations of LS and Newton: ' &
!                  ,REAL(totalIterLinearSolver),REAL(nNewton)



END SUBROUTINE TimeStepByEulerNewton
#endif /*PP_TimeDiscMethod==104*/

#if (PP_TimeDiscMethod==110) 
SUBROUTINE TimeStepByEulerImplicitParticle(t)
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
USE MOD_DG_Vars,                 ONLY:U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,iter
USE MOD_DG_Vars,                 ONLY: U,Ut
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_Equation,                ONLY: DivCleaningDamping
USE MOD_Equation,                ONLY: CalcSource
#ifdef PARTICLES
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartQ
USE MOD_Particle_Vars,           ONLY: PartState, LastPartPos, DelayTime, PEM, PDM
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS
USE MOD_ParticleSolver,          ONLY: ParticleNewton
!USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
!USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping
!#ifdef MPI
!USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
!#endif /*MPI*/
USE MOD_PIC_Analyze,             ONLY: CalcDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: tStage,coeff
!===================================================================================================================================

coeff=dt

!#ifdef maxwell
!! caution hard coded
!!IF (iter==0) CALL BuildPrecond(t,t,0,0.25,dt)
!IF (iter==0) CALL BuildPrecond(t,t,0,1.0,dt)
!!IF(PrecondType.GT.0)THEN
!  !IF (iter==0) CALL BuildJacDG()
!!END IF
!#endif /*maxwell*/

#ifdef PARTICLES
! Partilce locating
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL ParticleInserting()
END IF

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------

tStage=t
! source terms of field from U^n
IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
#endif

CALL DGTimeDerivative_weakForm(tStage, tStage, 0,doSource=.TRUE.)
! update to new time-level
U=U+dt*Ut

! implicit step of particles || requires fields of U^n+1
#ifdef PARTICLES
IF (t.GE.DelayTime) THEN
  tstage=t+dt
  PartQ(1,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,1)
  PartQ(2,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,2)
  PartQ(3,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,3)
  PartQ(4,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,4)
  PartQ(5,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,5)
  PartQ(6,1:PDM%ParticleVecLength)=PartState(1:PDM%ParticleVecLength,6)
  ! frozen matrix, e.g. no communication during stage?
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  ! does it make particle movement? depending on scheme?
  ! currently: closed scheme: tracking and pushing in Newton
  CALL ParticleNewton(tstage,coeff)
!  ! if no communication
!#ifdef MPI
!  ! open receive buffer for number of particles
!  CALL IRecvNbofParticles()
!#endif /*MPI*/
!  IF(DoRefMapping)THEN
!    CALL ParticleRefTracking()
!  ELSE
!   CALL ParticleTrackingCurved()
!  END IF
!#ifdef MPI
!  ! send number of particles
!  CALL SendNbOfParticles()
!  ! finish communication of number of particles and send particles
!  CALL MPIParticleSend()
!  ! finish communication
!  CALL MPIParticleRecv()
!#endif
!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
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
END IF
#endif

END SUBROUTINE TimeStepByEulerImplicitParticle
#endif /*PP_TimeDiscMethod==110*/


#if (PP_TimeDiscMethod==111) || (PP_TimeDiscMethod==112) 
SUBROUTINE TimeStepByIMPA(t)
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
USE MOD_DG_Vars,                 ONLY:U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,iter,iStage, nRKStages,time
USE MOD_TimeDisc_Vars,           ONLY: ERK_a,ESDIRK_a,RK_b,RK_c,RK_bs
USE MOD_DG_Vars,                 ONLY: U,Ut
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_LinearSolver,            ONLY: LinearSolver
USE MOD_Predictor,               ONLY: Predictor,StorePredictor
USE MOD_LinearSolver_Vars,       ONLY: ImplicitSource,LinSolverRHS
USE MOD_Equation,                ONLY: DivCleaningDamping
#ifdef maxwell
USE MOD_Precond,                 ONLY: BuildPrecond
USE MOD_Precond_Vars,            ONLY: PrecondType
#endif /*maxwell*/
USE MOD_JacDG,                   ONLY: BuildJacDG
USE MOD_Equation,                ONLY: CalcSource
#ifdef PARTICLES
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_PIC_Vars,                ONLY: PIC
USE MOD_Particle_Vars,           ONLY: PartStateN,PartStage
USE MOD_Particle_Vars,           ONLY: PartState, Pt, LastPartPos, DelayTime, PEM, PDM, usevMPF
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS, DSMC
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_PIC_Analyze,             ONLY: CalcDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: tstage
! implicit 
REAL               :: alpha
REAL               :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL               :: FieldStage (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
REAL               :: FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
REAL               :: tRatio, tphi
INTEGER            :: iCounter !, iStage
! explicit
!===================================================================================================================================

!#ifdef maxwell
!! caution hard coded
!!IF (iter==0) CALL BuildPrecond(t,t,0,0.25,dt)
!IF (iter==0) CALL BuildPrecond(t,t,0,RK_b(nRKStages),dt)
!!IF(PrecondType.GT.0)THEN
!  !IF (iter==0) CALL BuildJacDG()
!!END IF
!#endif /*maxwell*/
!
!Un = U
!FieldSource=0.
!! prediction with embadded scheme
!
!!tRatio= ! dt^n+1 / dt^n
!tRatio = 1. 
!
!#ifdef PARTICLES
!! Partilce locating
!IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!  CALL ParticleInserting()
!END IF
!#endif
!
! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------
!tStage=t
!! implicit
!#ifdef PARTICLES
!PartStateN(1:PDM%ParticleVecLength,1:6)=PartState(1:PDM%ParticleVecLength,1:6)
!IF (t.GE.DelayTime) THEN
!  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
!  CALL CalcPartRHS()
!END IF
!
! ----------------------------------------------------------------------------------------------------------------------------------
!!! explicit - Maxwell solver
!!IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!!  ! because of emmision and UpdateParticlePosition
!!  CALL Deposition(doInnerParts=.TRUE.)
!!#ifdef MPI
!!  ! here: finish deposition with delta kernal
!!  !       maps source terms in physical space
!!  ! ALWAYS require !! not so correct, implement particle send mode 1,2,3
!!  PartMPIExchange%nMPIParticles=0
!!#endif /*MPI*/
!!  CALL Deposition(doInnerParts=.FALSE.)
!!  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
!!END IF
!!! u1=un
!!#endif
!!ImplicitSource=0.
!!CALL CalcSource(tStage,1.,ImplicitSource)
!
! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to 6
! ----------------------------------------------------------------------------------------------------------------------------------
!DO iStage=2,nRKStages
!  ! compute f(u^s-1)
!  CALL DGTimeDerivative_weakForm(tStage, tStage, 0,doSource=.FALSE.)
!  ! store old values for use in next stages
!  FieldStage (:,:,:,:,:,iStage-1) = Ut
!  FieldSource(:,:,:,:,:,iStage-1) = ImplicitSource
!
!  ! time of current stage
!  tStage = t + RK_c(iStage)*dt
!  ! store predictor
! !--------------------------------------------------------------------------------------------------------------------------------
! ! explicit - field pusher
! !--------------------------------------------------------------------------------------------------------------------------------
!  ! requires source terms of particle at stage level iStage-1
!  IF ((t.GE.DelayTime).OR.((t.EQ.0).AND.(iStage.EQ.2)) THEN
!    ! because of emmision and UpdateParticlePosition
!    CALL Deposition(doInnerParts=.TRUE.)
!#ifdef MPI
!    ! here: finish deposition with delta kernal
!    !       maps source terms in physical space
!    ! ALWAYS require
!    PartMPIExchange%nMPIParticles=0
!#endif /*MPI*/
!    CALL Deposition(doInnerParts=.FALSE.)
!    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
!  END IF
!#endif /*PARTICLES*/
!  
!  U=Un
!  DO iCounter=1,iStage-1
!    U=U+ERK_a(iStage,iCounter)*dt*(FieldStage(:,:,:,:,:,iCounter)+FieldSource(:,:,:,:,:,iCounter))
!  END DO ! counter
!
!
!
!
!
!
!
!  !--------------------------------------------------------------------------------------------------------------------------------
!  ! implicit - Particle pusher
!  !--------------------------------------------------------------------------------------------------------------------------------
!#ifdef PARTICLES
!  ! particle RHS
!  IF (t.GE.DelayTime) THEN
!    CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
!    CALL CalcPartRHS()
!
!    PartStage(1:PDM%ParticleVecLength,1:3,iStage-1) = PartState(1:PDM%ParticleVecLength,4:6)
!    PartStage(1:PDM%ParticleVecLength,4:6,iStage-1) = Pt       (1:PDM%ParticleVecLength,1:3)
!  END IF
!  ! particle step
!  IF (t.GE.DelayTime) THEN
!    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
!    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
!    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
!    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
!    PartState(1:PDM%ParticleVecLength,1:6)=PartStateN(1:PDM%ParticleVecLength,1:6)
!    DO iCounter=1,iStage-1
!      PartState(1:PDM%ParticleVecLength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)   &
!                                                + ERK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,1:6,iCounter)
!    END DO ! counter
!  END IF
!
!#ifdef MPI
!  ! open receive buffer for number of particles
!  CALL IRecvNbofParticles()
!#endif /*MPI*/
!  IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!    IF(DoRefMapping)THEN
!      CALL ParticleRefTracking()
!    ELSE
!      CALL ParticleTrackingCurved()
!    END IF
!#ifdef MPI
!    ! send number of particles
!    CALL SendNbOfParticles()
!    ! finish communication of number of particles and send particles
!    CALL MPIParticleSend()
!    ! finish communication
!    CALL MPIParticleRecv()
!#endif
!  END IF
!
!
!  !--------------------------------------------------------------------------------------------------------------------------------
!  ! implicit - Maxwell solver
!  !--------------------------------------------------------------------------------------------------------------------------------
!  !--------------------------------------------------------------------------------------------------------------------------------
!  ! compute RHS for linear solver
!  LinSolverRHS=Un
!  DO iCounter = 1,iStage-1
!    LinSolverRHS = LinSolverRHS + ESDIRK_a(iStage,iCounter)*dt*(FieldStage(:,:,:,:,:,iCounter)+FieldSource(:,:,:,:,:,iCounter))
!  END DO
!  ! get predictor of u^s+1
!  CALL Predictor(iStage,dt,Un,FieldSource,FieldStage)
!
!  ImplicitSource=0.
!  alpha = ESDIRK_a(iStage,iStage)*dt
!  ! solve to new stage 
!  CALL LinearSolver(tstage,alpha)
!    ! damping
!  !CALL DivCleaningDamping()
!END DO
!
!!----------------------------------------------------------------------------------------------------------------------------------
!! update to next time level
!!----------------------------------------------------------------------------------------------------------------------------------
!
!#ifdef PARTICLES
!! only required for explicit equations // particle pusher
!! implicit coefficients are equal, therefore not required, however, would result in smoother solution
!! particle RHS
!IF (t.GE.DelayTime) THEN
!  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
!  CALL CalcPartRHS()
!!  PartStage(1:PDM%ParticleVecLength,1:3,6) = PartState(1:PDM%ParticleVecLength,4:6)
!!  PartStage(1:PDM%ParticleVecLength,4:6,6) = Pt       (1:PDM%ParticleVecLength,1:3)
!END IF
!
!
!IF (t.GE.DelayTime) THEN
!  LastPartPos(1:PDM%ParticleVecLength,1)  =PartState(1:PDM%ParticleVecLength,1)
!  LastPartPos(1:PDM%ParticleVecLength,2)  =PartState(1:PDM%ParticleVecLength,2)
!  LastPartPos(1:PDM%ParticleVecLength,3)  =PartState(1:PDM%ParticleVecLength,3)
!  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
!  ! compined
!  ! stage 6
!  PartState(1:PDM%ParticleVeclength,1:3) = PartStateN(1:PDM%ParticleVecLength,1:3)                   &
!                                           + RK_b(nRKStages)*dt*PartState(1:PDM%ParticleVecLength,4:6)
!  PartState(1:PDM%ParticleVeclength,4:6) = PartStateN(1:PDM%ParticleVecLength,4:6)                   &
!                                           + RK_b(nRKSTages)*dt*Pt(1:PDM%ParticleVecLength,1:3)
!!  ! stage 1 to 5
!!  PartState(1:PDM%ParticleVeclength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)                    &
!!                                           + RK_b(1)*dt*PartStage(1:PDM%ParticleVecLength,1:6,1)    &
!!                                           + RK_b(2)*dt*PartStage(1:PDM%ParticleVecLength,1:6,2)    &
!!                                           + RK_b(3)*dt*PartStage(1:PDM%ParticleVecLength,1:6,3)    &
!!                                           + RK_b(4)*dt*PartStage(1:PDM%ParticleVecLength,1:6,4)    &
!!                                           + RK_b(5)*dt*PartStage(1:PDM%ParticleVecLength,1:6,5) 
!!  PartState(1:PDM%ParticleVecLength,1:6)=PartStateN(1:PDM%ParticleVecLength,1:6)
!  DO iCounter=1,nRKStages-1
!    PartState(1:PDM%ParticleVecLength,1:6) = PartState(1:PDM%ParticleVecLength,1:6)   &
!                                              + RK_b(iCounter)*dt*PartStage(1:PDM%ParticleVecLength,1:6,iCounter)
!  END DO ! counter
!END IF
!
!! null iStage, because the old f(u^s) have not to be communicated
!iStage = 0
!IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
!#ifdef MPI
!  ! open receive buffer for number of particles
!  CALL IRecvNbofParticles()
!#endif /*MPI*/
!  IF(DoRefMapping)THEN
!    CALL ParticleRefTracking()
!  ELSE
!    CALL ParticleTrackingCurved()
!  END IF
!#ifdef MPI
!  ! send number of particles
!  CALL SendNbOfParticles()
!  ! finish communication of number of particles and send particles
!  CALL MPIParticleSend()
!  ! finish communication
!  CALL MPIParticleRecv()
!#endif /*MPI*/
!END IF
!
!!----------------------------------------------------------------------------------------------------------------------------------
!! DSMC
!!----------------------------------------------------------------------------------------------------------------------------------
!CALL UpdateNextFreePosition()
!IF (useDSMC) THEN
! CALL DSMC_main()
! PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
!                                        + DSMC_RHS(1:PDM%ParticleVecLength,1)
! PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
!                                        + DSMC_RHS(1:PDM%ParticleVecLength,2)
! PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
!                                        + DSMC_RHS(1:PDM%ParticleVecLength,3)
!END IF
!#endif /*PARTICLES*/
!
!! compute source and recompute field

END SUBROUTINE TimeStepByIMPA
#endif /*PP_TimeDiscMethod==111 || PP_TimeDiscMethod==112  */


#if (PP_TimeDiscMethod==200)
SUBROUTINE TimeStepByEulerStaticExp(t)
!===================================================================================================================================
! Static (using 4 or 8 variables, depending on compiled equation system (maxwell or electrostatic):
! Field is propagated until steady, then particle is moved
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,                 ONLY: U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,IterDisplayStep,iter,IterDisplayStepUser
USE MOD_TimeDisc_Vars,           ONLY: RK_a,RK_b,RK_c
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm
USE MOD_Filter,                  ONLY:Filter
USE MOD_Equation,                ONLY:DivCleaningDamping
USE MOD_Globals
#ifdef PP_POIS
USE MOD_Equation,                ONLY:DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,           ONLY:Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_PICDepo,                 ONLY : Deposition!, DepositionMPF
USE MOD_PICInterpolation,        ONLY : InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY : PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, dt_maxwell, MaxwellIterNum, usevMPF
USE MOD_part_RHS,                ONLY : CalcPartRHS
USE MOD_part_emission,           ONLY : ParticleInserting
USE MOD_DSMC,                    ONLY : DSMC_main
USE MOD_DSMC_Vars,               ONLY : useDSMC, DSMC_RHS
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_part_tools,              ONLY : UpdateNextFreePosition
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
#endif /*Particles*/
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
REAL                  :: Phit_temp(1:4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
REAL                  :: b_dt(1:5)
REAL                  :: dt_save, tStage, t_rk
#ifdef PARTICLES
REAL                  :: timeStart,timeEnd
#endif /*PARTICLES*/
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
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

dt_save = dt  !quick hack
t_rk = t
dt = dt_maxwell

DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO
IF(MOD(iter,IterDisplayStep).EQ.0) THEN
  SWRITE(*,'(A44,I0,A30)') " !!! EM field is being propagated... (every ",IterDisplayStepUser &
    ,". dt_maxwell is displayed) !!!"
END IF
DO iLoop = 1, MaxwellIterNum
  ! EM field
  IF ((MOD(iter,IterDisplayStep).EQ.0).AND.(MOD(iLoop,IterDisplayStepUser).EQ.0)) THEN
    SWRITE(*,'(A14,I0,A3,I0)') " MaxwellIter: ",iLoop," / ",MaxwellIterNum
  END IF
  CALL DGTimeDerivative_weakForm(t_rk,t_rk,0,doSource=.TRUE.)
  CALL DivCleaningDamping()
  Ut_temp = Ut 
  U = U + Ut*b_dt(1)
 
  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm(t_rk,tStage,0,doSource=.TRUE.)
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
USE MOD_DG_Vars,                 ONLY:U,Ut
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY:dt,IterDisplayStep,iter,IterDisplayStepUser
USE MOD_TimeDisc_Vars,           ONLY:RK_a,RK_b,RK_c
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm
USE MOD_Filter,                  ONLY:Filter
USE MOD_Equation,                ONLY:DivCleaningDamping
USE MOD_Globals
#ifdef PP_POIS
USE MOD_Equation,                ONLY:DivCleaningDamping_Pois,EvalGradient
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm_Pois
USE MOD_Equation_Vars,           ONLY:Phi,Phit,nTotalPhi
#endif
#ifdef PARTICLES
USE MOD_PICDepo,                 ONLY:Deposition!, DepositionMPF
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToParticle
USE MOD_PIC_Vars,                ONLY:PIC
USE MOD_Particle_Vars,           ONLY:PartState, Pt, LastPartPos, DelayTime, Time, PEM, PDM, dt_maxwell, MaxwellIterNum, usevMPF
USE MOD_part_RHS,                ONLY:CalcPartRHS
USE MOD_part_emission,           ONLY:ParticleInserting
USE MOD_DSMC,                    ONLY:DSMC_main
USE MOD_DSMC_Vars,               ONLY:useDSMC, DSMC_RHS, DSMC
USE MOD_PIC_Analyze,             ONLY:VerifyDepositedCharge
USE MOD_part_tools,              ONLY:UpdateNextFreePosition
USE MOD_Particle_Analyze_Vars,   ONLY:DoVerifyCharge
USE MOD_Particle_Tracking_vars,  ONLY:tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,       ONLY:ParticleTrackingCurved,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange
#endif /*MPI*/
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
#ifdef PARTICLES
REAL                          :: timeStart,timeEnd
#endif /*PARTICLES*/
!===================================================================================================================================
Time = t
IF (t.GE.DelayTime) CALL ParticleInserting()

IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
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

! EM field
dt_save = dt  !quick hack
t_rk = t
dt = dt_maxwell

DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO
IF(MOD(iter,IterDisplayStep).EQ.0) THEN
  SWRITE(*,'(A44,I0,A30)') " !!! EM field is being propagated... (every ",IterDisplayStepUser &
    ,". dt_maxwell is displayed) !!!"
END IF
DO iLoop = 1, MaxwellIterNum
  ! EM field
  IF ((MOD(iter,IterDisplayStep).EQ.0).AND.(MOD(iLoop,IterDisplayStepUser).EQ.0)) THEN
    SWRITE(*,'(A14,I0,A3,I0)') " MaxwellIter: ",iLoop," / ",MaxwellIterNum
  END IF

  CALL DGTimeDerivative_weakForm(t_rk,t_rk,0,doSource=.TRUE.)
  CALL DivCleaningDamping()
  Ut_temp = Ut 
  U = U + Ut*b_dt(1)
 
  DO rk=2,5
    tStage=t_rk+dt*RK_c(rk)
    ! field RHS
    CALL DGTimeDerivative_weakForm(t_rk,tStage,0,doSource=.TRUE.)
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
USE MOD_TimeDisc_Vars,ONLY:CFLScale,CFLtoOne
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

CFLToOne=1.0/CFLScale
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
IF(PP_N.GT.15) CALL abort(&
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
