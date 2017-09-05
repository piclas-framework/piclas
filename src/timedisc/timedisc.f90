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

#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTERFACE TimeStepPoissonByLSERK
  MODULE PROCEDURE TimeStepPoissonByLSERK
END INTERFACE
#endif

#if (PP_TimeDiscMethod==500)
INTERFACE TimeStepPoisson
  MODULE PROCEDURE TimeStepPoisson
END INTERFACE
#endif

PUBLIC :: InitTimeDisc,FinalizeTimeDisc
PUBLIC :: TimeDisc
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
PUBLIC :: TimeStepPoissonByLSERK
#endif
#if (PP_TimeDiscMethod==500)
PUBLIC :: TimeStepPoisson
#endif
!===================================================================================================================================

CONTAINS

SUBROUTINE InitTimeDisc()
!===================================================================================================================================
! Get information for end time and max time steps from ini file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,          ONLY:GetReal,GetInt, GETLOGICAL
USE MOD_TimeDisc_Vars,        ONLY:CFLScale,dt,TimeDiscInitIsDone,RKdtFrac,RKdtFracTotal
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
CALL abort(&
__STAMP__&
,' Preconditioner is only implemented for Maxwell!')
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
#elif (PP_TimeDiscMethod==110)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler-Implicit-Particles and Euler-Explicit-Field'
#elif (PP_TimeDiscMethod==111)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ESDIRK3-Particles and ERK3-Field'
#elif (PP_TimeDiscMethod==112)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ESDIRK4-Particles and ERK4-Field'
#elif (PP_TimeDiscMethod==120)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Heun/Crank-Nicolson1-2-2' 
#elif (PP_TimeDiscMethod==121)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ERK3/ESDIRK3-Particles and ESDIRK3-Field'
#elif (PP_TimeDiscMethod==122)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ERK4/ESDIRK4-Particles and ESDIRK4-Field'
#elif (PP_TimeDiscMethod==200)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Static Explicit'
#elif (PP_TimeDiscMethod==201)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler Static Explicit with adaptive TimeStep'
#elif (PP_TimeDiscMethod==300)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: FP Flow'
#elif (PP_TimeDiscMethod==500)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler, Poisson'
#elif (PP_TimeDiscMethod==501)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK3-3, Poisson'
#elif (PP_TimeDiscMethod==502)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-5, Poisson'
#elif (PP_TimeDiscMethod==506)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-14, Poisson'
#elif (PP_TimeDiscMethod==1000)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LD-Only'
#elif (PP_TimeDiscMethod==1001)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LD-DSMC'
# endif
RKdtFrac=1.
RKdtFracTotal=1.
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
USE MOD_TimeDisc_Vars,         ONLY: time,TEnd,dt,tAnalyze,iter,IterDisplayStep,DoDisplayIter,dt_Min
USE MOD_TimeAverage_vars,      ONLY: doCalcTimeAverage
USE MOD_TimeAverage,           ONLY: CalcTimeAverage
#if (PP_TimeDiscMethod==201)                                                                                                         
USE MOD_TimeDisc_Vars,         ONLY: dt_temp, MaximumIterNum 
#endif
USE MOD_Restart_Vars,          ONLY: DoRestart,RestartTime
#ifndef PP_HDG
USE MOD_CalcTimeStep,          ONLY: CalcTimeStep
USE MOD_PML_Vars,              ONLY: DoPML,DoPMLTimeRamp,PMLTimeRamp
USE MOD_PML,                   ONLY: PMLTimeRamping
#endif /*PP_HDG*/
USE MOD_Analyze,               ONLY: CalcError,PerformAnalyze
USE MOD_Analyze_Vars,          ONLY: Analyze_dt
#ifdef PARTICLES
USE MOD_Particle_Analyze,      ONLY: AnalyzeParticles
USE MOD_HDF5_output,           ONLY: WriteIMDStateToHDF5
#else
USE MOD_AnalyzeField,          ONLY: AnalyzeField
#endif /*PARTICLES*/
USE MOD_HDF5_output,           ONLY: WriteStateToHDF5
USE MOD_Mesh_Vars,             ONLY: MeshFile,nGlobalElems,DoWriteStateToHDF5
USE MOD_Mesh,                  ONLY: SwapMesh
USE MOD_Filter,                ONLY: Filter
USE MOD_RecordPoints_Vars,     ONLY: RP_onProc
USE MOD_RecordPoints,          ONLY: RecordPoints,WriteRPToHDF5
#ifdef PARTICLES
USE MOD_PICDepo,               ONLY: Deposition
USE MOD_PICDepo_Vars,          ONLY: DepositionType
USE MOD_Particle_Output,       ONLY: Visualize_Particles
USE MOD_PARTICLE_Vars,         ONLY: WriteMacroVolumeValues, MacroValSampTime,DoImportIMDFile
USE MOD_Particle_Tracking_vars,ONLY: tTracking,tLocalization,nTracks,MeasureTrackTime,CountNbOfLostParts,nLostParts
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
USE MOD_DSMC_Vars,             ONLY: Iter_macvalout,Iter_macsurfvalout
#ifdef MPI
USE MOD_Particle_MPI,          ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#endif /*PARTICLES*/
#ifdef PP_POIS
USE MOD_Equation,              ONLY: EvalGradient
#endif /*PP_POIS*/
USE MOD_LoadBalance_Vars,      ONLY: nSkipAnalyze
#ifdef MPI
USE MOD_LoadBalance,           ONLY: LoadBalance,LoadMeasure,ComputeElemLoad
USE MOD_LoadBalance_Vars,      ONLY: DoLoadBalance
!USE MOD_Particle_MPI_Vars,    ONLY: PartMPI
#endif /*MPI*/
#if defined(IMEX) || (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
USE MOD_LinearSolver_Vars,     ONLY:totalIterLinearSolver
#endif /*IMEX*/
#ifdef IMPA
USE MOD_LinearSolver_vars,     ONLY:nPartNewton
USE MOD_LinearSolver_Vars,     ONLY:TotalPartIterLinearSolver
#endif /*IMPA*/
#if (PP_TimeDiscMethod==120)||(PP_TimeDiscMethod==121||PP_TimeDiscMethod==122)
USE MOD_LinearSolver_Vars,    ONLY: totalFullNewtonIter
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: tFuture,tZero
INTEGER(KIND=8)              :: nAnalyze
INTEGER                      :: iAnalyze
REAL                         :: tEndDiff, tAnalyzeDiff
#if (PP_TimeDiscMethod==201)
REAL                         :: vMax,vMaxx,vMaxy,vMaxz
#endif
INTEGER(KIND=8)              :: iter_loc
REAL                         :: CalcTimeStart,CalcTimeEnd
INTEGER                      :: TimeArray(8)              ! Array for system time
REAL                         :: CurrentImbalance
LOGICAL                      :: PerformLoadBalance
#if (PP_TimeDiscMethod==201)
INTEGER                      :: iPart
LOGICAL                      :: NoPartInside
#endif
LOGICAL                      :: finalIter
#ifdef PARTICLES
INTEGER                      :: nLostPartsTot
#endif /*PARTICLES*/
!===================================================================================================================================
! init
SWRITE(UNIT_StdOut,'(132("-"))')
IF(.NOT.DoRestart)THEN
  time=0.
  SWRITE(UNIT_StdOut,*)'INITIAL PROJECTION:'
ELSE
  time=RestartTime
  SWRITE(UNIT_StdOut,*)'REWRITING SOLUTION:'
END IF
#ifdef PARTICLES
iter_macvalout=0
iter_macsurfvalout=0
IF (WriteMacroVolumeValues) MacroValSampTime = Time
#endif /*PARTICLES*/
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
#ifdef PARTICLES
IF(DoImportIMDFile) CALL WriteIMDStateToHDF5(time) ! write IMD particles to state file (and TTM if it exists)
#endif /*PARTICLES*/
IF(DoWriteStateToHDF5) CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture)

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
IF ((TRIM(DepositionType).EQ."shape_function")             &
.OR.(TRIM(DepositionType).EQ."shape_function_1d")          &
.OR.(TRIM(DepositionType).EQ."shape_function_spherical")   &
.OR.(TRIM(DepositionType).EQ."shape_function_simple")      &
.OR.(TRIM(DepositionType).EQ."shape_function_cylindrical"))THEN
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

!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
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
CALL InitTimeStep() ! Initial time step calculation
CalcTimeStart=BOLTZPLATZTIME()
iter=0
iter_loc=0

! fill recordpoints buffer (first iteration)
!IF(RP_onProc) CALL RecordPoints(iter,t,forceSampling=.TRUE.) 

! fill initial analyze stuff
tAnalyzeDiff=tAnalyze-time    ! time to next analysis, put in extra variable so number does not change due to numerical errors
tEndDiff=tEnd-time            ! dito for end time
dt=MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/)) ! quick fix: set dt for initial write DSMCHOState (WriteMacroVolumeValues=T)
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
#endif /*(PP_TimeDiscMethod==201)*/

#ifdef PARTICLES
  IF(enableParticleMerge) THEN
    IF ((iter.GT.0).AND.(MOD(iter,vMPFMergeParticleIter).EQ.0)) doParticleMerge=.true.
  END IF
#endif /*PARTICLES*/

  tAnalyzeDiff=tAnalyze-time    ! time to next analysis, put in extra variable so number does not change due to numerical errors
  tEndDiff=tEnd-time            ! dito for end time

  !IF(time.LT.3e-8)THEN
  !    !RETURN
  !ELSE
  !  IF(time.GT.4e-8) dt_Min=MIN(dt_Min*1.2,2e-8)
  !END IF

  dt=MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/))
  IF (tAnalyzeDiff-dt.LT.dt/100.0) dt = tAnalyzeDiff
  IF (tEndDiff-dt.LT.dt/100.0) dt = tEndDiff
  IF ( dt .LT. 0. ) THEN
    SWRITE(UNIT_StdOut,*)'*** ERROR: Is something wrong with the defined tEnd?!? ***'
    CALL abort(&
    __STAMP__&
    ,'Error in tEndDiff or tAnalyzeDiff!')
  END IF

  IF(doCalcTimeAverage) CALL CalcTimeAverage(.FALSE.,dt,time,tFuture)

#ifndef PP_HDG
  IF(DoPML)THEN
    IF(DoPMLTimeRamp)THEN
      CALL PMLTimeRamping(time,PMLTimeRamp)
    END IF
  END IF
#endif /*NOT PP_HDG*/

! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
#if (PP_TimeDiscMethod==1)
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==2)
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==3)
  CALL TimeStepByTAYLOR(time)
#elif (PP_TimeDiscMethod==4)
  CALL TimeStep_DSMC()
#elif (PP_TimeDiscMethod==5)
  CALL TimeStepByRK4EulerExpl(time)
#elif (PP_TimeDiscMethod==6)
  CALL TimeStepByLSERK(time,iter,tEndDiff)
#elif (PP_TimeDiscMethod==42)
  CALL TimeStep_DSMC_Debug() ! Reservoir and Debug
#elif (PP_TimeDiscMethod==100)
  CALL TimeStepByEulerImplicit(time) ! O1 Euler Implicit
#elif (PP_TimeDiscMethod==101)
  CALL TimeStepByIMEXRK(time) ! ) O3 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==102)
  CALL TimeStepByIMEXRK(time) ! O4 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==103)
  CALL TimeStepByIMEXRK(time) ! ) O3 ERK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==110)
  CALL TimeStepByEulerImplicitParticle(time) ! ) O3 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==111)
  CALL TimeStepByIMPA(time) ! ) O3 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==112)
  CALL TimeStepByIMPA(time) ! O4 ESDIRK Particles + ERK Field 
#elif (PP_TimeDiscMethod==120)
  CALL TimeStepByImplicitRK(time) !  O3 ERK/ESDIRK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==121)
  CALL TimeStepByImplicitRK(time) !  O3 ERK/ESDIRK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==122)
  CALL TimeStepByImplicitRK(time) ! O4 ERK/ESDIRK Particles + ESDIRK Field 
#elif (PP_TimeDiscMethod==200)
  CALL TimeStepByEulerStaticExp(time) ! O1 Euler Static Explicit
#elif (PP_TimeDiscMethod==201)
  CALL TimeStepByEulerStaticExpAdapTS(time) ! O1 Euler Static Explicit with adaptive TimeStep
#elif (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=506)
#ifdef PP_HDG
#if (PP_TimeDiscMethod==500)
  CALL TimeStepPoisson(time) ! Euler Explicit, Poisson
#else
  CALL TimeStepPoissonByLSERK(time,iter,tEndDiff) ! Runge Kutta Explicit, Poisson
#endif
#else
  CALL abort(&
  __STAMP__&
  ,'Timedisc 50x only available for EQNSYS Poisson!',PP_N,999.)
#endif /*PP_HDG*/
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
#if (PP_TimeDiscMethod!=1)&&(PP_TimeDiscMethod!=2)&&(PP_TimeDiscMethod!=6)&&(PP_TimeDiscMethod<501||PP_TimeDiscMethod>506)
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
    CALL ComputeElemLoad(CurrentImbalance,PerformLoadBalance,time)
#endif /*MPI*/
    ! future time
    nAnalyze=nAnalyze+1
    tFuture=tZero+REAL(nAnalyze)*Analyze_dt
#ifdef MPI
    IF(iAnalyze.EQ.nSkipAnalyze .OR. PerformLoadBalance .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#else
    IF( iAnalyze.EQ.nSkipAnalyze .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#endif /*MPI*/
#ifdef PARTICLES
    IF(CountNbOfLostParts)THEN
#ifdef MPI
      IF(MPIRoot) THEN
        CALL MPI_REDUCE(nLostParts,nLostPartsTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
      ELSE ! no Root
        CALL MPI_REDUCE(nLostParts,nLostPartsTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
      END IF
#else
      nLostPartsTot=nLostParts
#endif /*MPI*/
    END IF
#endif /*PARICLES*/
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
#ifdef PARTICLES
        IF(CountNbOfLostParts)THEN
          WRITE(UNIT_stdOut,'(A,I12)')' NbOfLostParticle : ',nLostPartsTot
        END IF
#endif /*PARICLES*/
      END IF !MPIroot
#if defined(IMEX) || (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
      SWRITE(UNIT_stdOut,'(132("="))')
      SWRITE(UNIT_stdOut,'(A32,I12)') ' Total iteration Linear Solver    ',totalIterLinearSolver
      TotalIterLinearSolver=0
#endif /*IMEX*/
#ifdef IMPA
      SWRITE(UNIT_stdOut,'(A32,I12)')  ' IMPLICIT PARTICLE TREATMENT    '
      SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration Newton         ',nPartNewton
      SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration GMRES          ',TotalPartIterLinearSolver
      IF(nPartNewton.GT.0)THEN
        SWRITE(UNIT_stdOut,'(A35,F12.2)')' Average GMRES steps per Newton    ',REAL(TotalPartIterLinearSolver)&
                                                                              /REAL(nPartNewton)
      END IF
#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
      SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration outer-Newton    ',TotalFullNewtonIter
      totalFullNewtonIter=0
#endif 
      nPartNewTon=0
      TotalPartIterLinearSolver=0
      SWRITE(UNIT_stdOut,'(132("="))')
#endif /*IMPA*/
      ! Analyze for output
      CALL PerformAnalyze(tAnalyze,iter,tenddiff,forceAnalyze=.FALSE.,OutPut=.TRUE.,LastIter_In=finalIter)
#ifndef PP_HDG
#endif /*PP_HDG*/
      ! Write state to file
      IF(DoWriteStateToHDF5) CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture)
      IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,time,tFuture)
#ifndef PP_HDG
      ! Write state to file
#endif /*PP_HDG*/
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
!#ifndef PP_HDG
!      dt_Min=CALCTIMESTEP()
!#endif /*PP_HDG*/
      IF(PerformLoadBalance .AND. iAnalyze.NE.nSkipAnalyze) &
        CALL PerformAnalyze(time,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.TRUE.)
      IF(PerformLoadBalance) CALL InitTimeStep() ! re-calculate time step after load balance is performed
!      dt=dt_Min !not sure if nec., was here before InitTimtStep was created, overwritten in next iter anyway
      ! CALL WriteStateToHDF5(TRIM(MeshFile),time,tFuture) ! not sure if required
    END IF
#endif /*MPI*/
    CalcTimeStart=BOLTZPLATZTIME()
  ENDIF   
  IF(time.GE.tEnd)EXIT ! done, worst case: one additional time step
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
USE MOD_PML_Vars,                ONLY: U2,U2t,nPMLElems,DoPML,PMLnVar
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
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,           ONLY: StartParticleMerge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#ifdef MPI
USE MOD_Particle_MPI_Vars,       ONLY: DoExternalParts
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars,       ONLY: ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM
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
REAL                          :: U2t_temp(1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems) ! temporal variable for U2t
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
CALL CountPartsPerElem()
#endif /*MPI*/

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
  END IF
#endif /*MPI*/
END IF

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

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
  END IF 
#endif /*MPI*/
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  IF(DoExternalParts)THEN
    ! finish communication
    CALL MPIParticleRecv()
  END IF
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
CALL DivCleaningDamping()

#ifdef PP_POIS
! Potential
CALL DGTimeDerivative_weakForm_Pois(t,t,0)
CALL DivCleaningDamping_Pois()
#endif /*PP_POIS*/


! calling the analyze routines
! Analysis is called in first RK-stage of NEXT iteration, however, the iteration count is performed AFTER the time step,
! hence, this is the correct iteration for calling the analysis routines.
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

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
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
    CALL ParticleTracing()
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
  CALL DivCleaningDamping()
#ifdef PP_POIS
  CALL DGTimeDerivative_weakForm_Pois(t,tStage,0)
  CALL DivCleaningDamping_Pois()
#endif



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
      CALL ParticleTracing()
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
IF (t.GE.DelayTime) THEN
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  CALL ParticleInserting()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tLocalization=tLocalization+TimeEnd-TimeStart
  END IF
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(4)=tCurrent(4)+tLBEnd-tLBStart
  tLBStart = LOCALTIME() ! LB Time Start
  CALL MPIParticleSend()
  CALL MPIParticleRecv()
  PartMPIExchange%nMPIParticles=0
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(10)=tCurrent(10)+tLBEnd-tLBStart
#endif /*MPI*/
END IF 

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
SUBROUTINE TimeStep_DSMC()
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,    ONLY: dt, IterDisplayStep, iter, TEnd, Time
#ifdef PARTICLES
USE MOD_Globals,          ONLY : abort
USE MOD_Particle_Vars,    ONLY : KeepWallParticles
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, PDM, PEM, DoSurfaceFlux, WriteMacroVolumeValues
USE MOD_DSMC_Vars,        ONLY : DSMC_RHS, DSMC, CollisMode
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting, ParticleSurfaceflux
USE MOD_Particle_Tracking_vars, ONLY: tTracking,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,ONLY: ParticleTracing,ParticleRefTracking
USE MOD_DSMC_SurfModel_Tools,   ONLY: Calc_PartNum_Wall_Desorb, DSMC_Update_Wall_Vars
#ifdef MPI
USE MOD_Particle_MPI,     ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: timeEnd, timeStart
INTEGER :: iPart
REAL    :: RandVal, dtFrac
!===================================================================================================================================

  IF (DoSurfaceFlux) THEN
    ! Calculate desobing particles for Surfaceflux
    IF ((.NOT.KeepWallParticles) .AND. (DSMC%WallModel.GT.0)) THEN
     CALL Calc_PartNum_Wall_Desorb()
    END IF
    CALL ParticleSurfaceflux()
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%dtFracPush(iPart)) THEN
          LastPartPos(iPart,1)=PartState(iPart,1)
          LastPartPos(iPart,2)=PartState(iPart,2)
          LastPartPos(iPart,3)=PartState(iPart,3)
          PEM%lastElement(iPart)=PEM%Element(iPart)
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dt
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dt
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dt
        ELSE !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RandVal
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
    END DO
  ELSE
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + PartState(1:PDM%ParticleVecLength,4) * dt
    PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + PartState(1:PDM%ParticleVecLength,5) * dt
    PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + PartState(1:PDM%ParticleVecLength,6) * dt
  END IF
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*MPI*/
  CALL DSMC_Update_Wall_Vars()
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  ! actual tracking
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTracing()
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
  IF (CollisMode.NE.0) THEN
    CALL UpdateNextFreePosition()
  ELSE IF ( (MOD(iter,IterDisplayStep).EQ.0) .OR. &
            (Time.ge.(1-DSMC%TimeFracSamp)*TEnd) .OR. &
            WriteMacroVolumeValues ) THEN
    CALL UpdateNextFreePosition() !postpone UNFP for CollisMode=0 to next IterDisplayStep or when needed for DSMC-Sampling
  ELSE IF (PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).GT.PDM%maxParticleNumber .OR. &
           PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).EQ.0) THEN
    CALL abort(&
    __STAMP__&
    ,'maximum nbr of particles reached!')  !gaps in PartState are not filled until next UNFP and array might overflow more easily!
  END IF
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
USE MOD_TimeDisc_Vars,ONLY: dt,time
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
USE MOD_Particle_Vars,    ONLY : PartState, Pt, LastPartPos, PEM, PDM, usevMPF, doParticleMerge, DelayTime, PartPressureCell
USE MOD_part_RHS,         ONLY : CalcPartRHS
USE MOD_part_emission,    ONLY : ParticleInserting
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_DSMC_Vars,        ONLY : useDSMC, DSMC_RHS
USE MOD_part_MPFtools,    ONLY : StartParticleMerge
USE MOD_PIC_Analyze,      ONLY: VerifyDepositedCharge
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_Particle_Tracking_vars, ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,ONLY: ParticleTracing,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
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
REAL                  :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) ! temporal variable for Ut
#ifdef PP_POIS
REAL                  :: Phit_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#endif
REAL                  :: tStage,b_dt(1:5)
INTEGER               :: rk
REAL                  :: timeEnd, timeStart
!===================================================================================================================================
IF (t.GE.DelayTime) CALL ParticleInserting()
!CALL UpdateNextFreePosition()
DO rk=1,5
  b_dt(rk)=RK_b(rk)*dt   ! TBD: put in initiation (with maxwell we are linear!!!)
END DO

!IF(t.EQ.0) CALL Deposition()
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)

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
IF ((t.GE.DelayTime)) THEN
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  ! actual tracking
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTracing()
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
SUBROUTINE TimeStep_DSMC_Debug()
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
USE MOD_Particle_Vars,    ONLY : DoSurfaceFlux, KeepWallParticles
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos, PDM,PEM!, Species, PartSpecies
USE MOD_DSMC_Vars,        ONLY : DSMC_RHS, DSMC!, Debug_Energy,PartStateIntEn
USE MOD_DSMC,             ONLY : DSMC_main
USE MOD_part_tools,       ONLY : UpdateNextFreePosition
USE MOD_part_emission,    ONLY : ParticleInserting, ParticleSurfaceflux
USE MOD_Particle_Tracking_vars, ONLY: tTracking,DoRefMapping,MeasureTrackTime
USE MOD_Particle_Tracking,ONLY: ParticleTracing,ParticleRefTracking
USE MOD_DSMC_SurfModel_Tools,   ONLY: Calc_PartNum_Wall_Desorb, DSMC_Update_Wall_Vars
#ifdef MPI
USE MOD_Particle_MPI,     ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPart
REAL                  :: timeStart, timeEnd, RandVal, dtFrac
!===================================================================================================================================

IF (DSMC%ReservoirSimu) THEN ! fix grid should be defined for reservoir simu
  ! Calculate desobing particles for Surfaceflux
  IF ((.NOT.KeepWallParticles) .AND. (DSMC%WallModel.GT.0)) THEN
    CALL Calc_PartNum_Wall_Desorb()
  END IF
  CALL DSMC_Update_Wall_Vars()
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
  IF (DoSurfaceFlux) THEN
    ! Calculate desobing particles for Surfaceflux
    IF ((.NOT.KeepWallParticles) .AND. (DSMC%WallModel.GT.0)) THEN
     CALL Calc_PartNum_Wall_Desorb()
    END IF
    
    CALL ParticleSurfaceflux()
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%dtFracPush(iPart)) THEN
          LastPartPos(iPart,1)=PartState(iPart,1)
          LastPartPos(iPart,2)=PartState(iPart,2)
          LastPartPos(iPart,3)=PartState(iPart,3)
          PEM%lastElement(iPart)=PEM%Element(iPart)
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dt
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dt
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dt
        ELSE !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RandVal
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
    END DO
  ELSE
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    ! bugfix if more than 2.x mio (2000001) particle per process
    ! tested with 64bit Ubuntu 12.04 backports
    DO iPart = 1, PDM%ParticleVecLength
      PEM%lastElement(iPart)=PEM%Element(iPart)
    END DO
    !PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + PartState(1:PDM%ParticleVecLength,4) * dt
    PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + PartState(1:PDM%ParticleVecLength,5) * dt
    PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + PartState(1:PDM%ParticleVecLength,6) * dt
  END IF
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbOfParticles()
#endif /*MPI*/
  CALL DSMC_Update_Wall_Vars()
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  ! actual tracking
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTracing()
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
    CALL ParticleTracing()
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
USE MOD_TimeDisc_Vars,           ONLY: ERK_a,ESDIRK_a,RK_b,RK_c
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
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking
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
IF (iter==0) CALL BuildPrecond(t,t,0,RK_b(nRKStages),dt)
#endif /*maxwell*/

Un = U
FieldStage =0.
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
      CALL ParticleTracing()
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
    CALL ParticleTracing()
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
!USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking
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
!   CALL ParticleTracing()
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
USE MOD_TimeDisc_Vars,           ONLY: ERK_a,ESDIRK_a,RK_b,RK_c
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
USE MOD_Particle_Vars,           ONLY: PartStateN,PartStage, PartQ
USE MOD_Particle_Vars,           ONLY: PartState, Pt, LastPartPos, DelayTime, PEM, PDM, usevMPF
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS, DSMC
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping,MeasureTrackTime
USE MOD_ParticleSolver,          ONLY: ParticleNewton
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
!REAL               :: FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
REAL               :: tRatio, tphi
INTEGER            :: iCounter !, iStage
! explicit
!===================================================================================================================================

tRatio = 1.0
Un          = U
!FieldSource = 0.

#ifdef PARTICLES
! Partilce locating
IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
  CALL ParticleInserting()
END IF
PartStateN(1:PDM%ParticleVecLength,1:6)=PartState(1:PDM%ParticleVecLength,1:6)
#endif

!! stage 1, explicit
!IF (nNewtonIterGlobal==0) THEN
!  CALL DGTimeDerivative_weakForm(t,t,0)
!  R_Stage(:,:,:,:,:,1)=Ut
!ELSE
!  ! Ut is stored in R_Xk from the last Newton iteration
!  R_Stage(:,:,:,:,:,1)=R_Xk
!END IF

 
! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------

!tStage=t
! explicit - Maxwell solver
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
!  ImplicitSource=0.
!  CALL CalcSource(tStage,1.,ImplicitSource)
END IF

!#ifdef PARTICLES
!IF (t.GE.DelayTime) THEN
!  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
!  CALL CalcPartRHS()
!END IF
!#endif /*PARTICLES*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to 6
! ----------------------------------------------------------------------------------------------------------------------------------
DO iStage=2,nRKStages
  ! compute f(u^s-1)
  ! it's for the implicit part, here: PARTICLES
#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
    CALL CalcPartRHS()
    PartStage(1:PDM%ParticleVecLength,1:3,iStage-1) = PartState(1:PDM%ParticleVecLength,4:6)
    PartStage(1:PDM%ParticleVecLength,4:6,iStage-1) = Pt       (1:PDM%ParticleVecLength,1:3)
  END IF
#endif /*PARTICLES*/

  ! time of current stage
  tStage = t + RK_c(iStage)*dt
  !--------------------------------------------------------------------------------------------------------------------------------
  ! explicit - field pusher
  !--------------------------------------------------------------------------------------------------------------------------------
  CALL DGTimeDerivative_weakForm(tStage, tStage, 0,doSource=.TRUE.)
  ! store old values for use in next stages
  FieldStage (:,:,:,:,:,iStage-1) = Ut
  !FieldSource(:,:,:,:,:,iStage-1) = ImplicitSource
  U=Un
  DO iCounter = 1,iStage-1
    U = U + ERK_a(iStage,iCounter)*dt*(FieldStage(:,:,:,:,:,iCounter)) !+FieldSource(:,:,:,:,:,iCounter))
  END DO

  !--------------------------------------------------------------------------------------------------------------------------------
  ! implicit - Particle pusher
  !--------------------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    ! old state
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    ! compute Q and U
    PartQ(1,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,1)
    PartQ(2,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,2)
    PartQ(3,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,3)
    PartQ(4,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,4)
    PartQ(5,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,5)
    PartQ(6,1:PDM%ParticleVecLength)=PartStateN(1:PDM%ParticleVecLength,6)
    DO iCounter = 1, iStage-1
      PartQ(1,1:PDM%ParticleVecLength)=PartQ(1,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,1,iCounter)
      PartQ(2,1:PDM%ParticleVecLength)=PartQ(2,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,2,iCounter)
      PartQ(3,1:PDM%ParticleVecLength)=PartQ(3,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,3,iCounter)
      PartQ(4,1:PDM%ParticleVecLength)=PartQ(4,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,4,iCounter)
      PartQ(5,1:PDM%ParticleVecLength)=PartQ(5,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,5,iCounter)
      PartQ(6,1:PDM%ParticleVecLength)=PartQ(6,1:PDM%ParticleVecLength) &
                                      + ESDIRK_a(iStage,iCounter)*dt*PartStage(1:PDM%ParticleVecLength,6,iCounter)
    END DO 
    ! predictor
    tphi = RK_c(iStage)
    ! do not change PartState .... working??
    !PartState(1:PDM%ParticleVecLength,1:6)=PartStateN(1:PDM%ParticleVecLength,1:6)
    !DO iCounter = 1,iStage-1
    !  PartState(1:PDM%ParticleVecLength,1:6)=PartState(1:PDM%ParticleVecLength,1:6) &
    !                +         (RK_bs(iCounter,1)*tphi+RK_bs(iCounter,2)*tphi**2)*dt &
    !                +          PartStage(1:PDM%ParticleVecLength,1:6,iCounter)
    !END DO
    alpha = ESDIRK_a(iStage,iStage)*dt
    CALL ParticleNewton(tstage,alpha)
    ! move particle
#ifdef MPI
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
#endif /*MPI*/
    IF(DoRefMapping)THEN
      CALL ParticleRefTracking()
    ELSE
      CALL ParticleTracing()
    END IF
#ifdef MPI
    ! here: could use deposition as hiding, not done yet
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
#endif /*MPI*/
    
    ! funny, can use interpolation and deposition as hiding :)
    ! recompute particle source terms on field solver :)
    CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
    ! here: finish deposition with delta kernal
    !       maps source terms in physical space
    ! ALWAYS require
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    CALL Deposition(doInnerParts=.FALSE.)
  END IF
#endif /*PARTICLES*/
END DO
! DEBUGGGGGG
! missssisssssing !
!  FieldStage (:,:,:,:,:,iStage-1) = Ut
!  FieldSource(:,:,:,:,:,iStage-1) = ImplicitSource


! particles implicit done
!----------------------------------------------------------------------------------------------------------------------------------
! update explicit to next time level
!----------------------------------------------------------------------------------------------------------------------------------

CALL DGTimeDerivative_weakForm(tStage, tStage, 0,doSource=.TRUE.)
! store old values for use in next stages
U=Un+RK_b(nRKStages)*dt*Ut
DO iCounter = 1,nRKStages-1
  U = U + RK_b(iCounter)*dt*(FieldStage(:,:,:,:,:,iCounter)) !+FieldSource(:,:,:,:,:,iCounter))
END DO


!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
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


END SUBROUTINE TimeStepByIMPA
#endif /*PP_TimeDiscMethod==111 || PP_TimeDiscMethod==112  */


#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
SUBROUTINE TimeStepByImplicitRK(t)
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
USE MOD_TimeDisc_Vars,           ONLY:dt,iter,iStage, nRKStages
USE MOD_TimeDisc_Vars,           ONLY:ERK_a,ESDIRK_a,RK_b,RK_c,RKdtFrac
USE MOD_LinearSolver_Vars,       ONLY:ImplicitSource, ExplicitSource,DoPrintConvInfo
#ifdef PP_HDG
USE MOD_Equation,                ONLY:CalcSourceHDG
#else /*pure DG*/
USE MOD_DG_Vars,                 ONLY:U
USE MOD_DG_Vars,                 ONLY:Ut
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm
USE MOD_Predictor,               ONLY:Predictor,StorePredictor
USE MOD_LinearSolver_Vars,       ONLY:LinSolverRHS
USE MOD_Equation,                ONLY:DivCleaningDamping
USE MOD_Equation,                ONLY:CalcSource
USE MOD_TimeDisc_Vars,      ONLY: dt_Min
#ifdef maxwell
USE MOD_Precond,                 ONLY:BuildPrecond
#endif /*maxwell*/
#endif /*PP_HDG*/
USE MOD_Newton,                  ONLY:ImplicitNorm,FullNewton
USE MOD_Equation_Vars,           ONLY:c2_inv
#ifdef PARTICLES
USE MOD_Timedisc_Vars,           ONLY:RKdtFrac,RKdtFracTotal
USE MOD_LinearSolver_Vars,       ONLY:DoUpdateInStage
USE MOD_Predictor,               ONLY:PartPredictor,PredictorType
USE MOD_Particle_Vars,           ONLY:PartIsImplicit,PartLorentzType, doParticleMerge,PartPressureCell
USE MOD_Particle_Analyze_Vars,   ONLY:DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY:VerifyDepositedCharge
USE MOD_PICDepo,                 ONLY:Deposition
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY:PartStateN,PartStage, PartQ,Species,nSpecies
USE MOD_Particle_Vars,           ONLY:PartState, Pt, LastPartPos, DelayTime, PEM, PDM,  DoSurfaceFlux!, StagePartPos
USE MOD_part_RHS,                ONLY:CalcPartRHS,PartVeloToImp
USE MOD_part_emission,           ONLY:ParticleInserting, ParticleSurfaceflux
USE MOD_DSMC,                    ONLY:DSMC_main
USE MOD_DSMC_Vars,               ONLY:useDSMC, DSMC_RHS
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_ParticleSolver,          ONLY:ParticleNewton, SelectImplicitParticles
USE MOD_Part_RHS,                ONLY:SLOW_RELATIVISTIC_PUSH,FAST_RELATIVISTIC_PUSH&
                                     ,RELATIVISTIC_PUSH,NON_RELATIVISTIC_PUSH
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
USE MOD_part_MPFtools,           ONLY:StartParticleMerge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange
USE MOD_Particle_MPI_Vars,       ONLY:DoExternalParts
USE MOD_Particle_MPI_Vars,       ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM
#endif /*MPI*/
USE MOD_PIC_Analyze,             ONLY:CalcDepositedCharge
USE MOD_part_tools,              ONLY:UpdateNextFreePosition
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
REAL               :: sgamma
#ifndef PP_HDG
REAL               :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL               :: FieldStage (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
REAL               :: FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
#else 
INTEGER            :: iElem,i,j,k
#endif /*DG*/
REAL               :: tRatio, LorentzFacInv
! particle surface flux
REAL               :: dtFrac,RandVal, LorentzFac
! RK counter
INTEGER            :: iCounter !, iStage
#ifdef PARTICLES
INTEGER            :: iPart,iSpec,iSF
#endif /*PARTICLES*/
!===================================================================================================================================

#ifndef PP_HDG
#ifdef maxwell
! caution hard coded
IF (iter==0) CALL BuildPrecond(t,t,0,RK_b(nRKStages),dt)
#endif /*maxwell*/
Un          = U
FieldSource = 0.
#endif /*DG*/
tRatio = 1.

#ifdef PARTICLES
! Partilce locating
IF (t.GE.DelayTime) THEN
  CALL ParticleInserting() ! do not forget to communicate the emmitted particles ... for shape function
END IF
! select, if particles are treated implicitly or explicitly
CALL SelectImplicitParticles()
!RKSumC_inv=1.0/SUM(RK_C)
! partstaten is set AFTER VeloToImp
#endif

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 1 - initialization
! ----------------------------------------------------------------------------------------------------------------------------------

tStage=t

#ifdef PARTICLES
IF((t.GE.DelayTime).OR.(iter.EQ.0))THEN
! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
    ! set exchanged number of particles to zero
    PartMPIExchange%nMPIParticles=0
  END IF
#endif /*MPI*/
END IF

! simulation with delay-time, compute the
IF(DelayTime.GT.0.)THEN
  IF((iter.EQ.0).AND.(t.LT.DelayTime))THEN
    ! perform normal deposition
    CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
    ! here: finish deposition with delta kernal
    !       maps source terms in physical space
    ! ALWAYS require
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    CALL Deposition(doInnerParts=.FALSE.)
  END IF
END IF

! compute source of first stage for Maxwell solver
IF (t.GE.DelayTime) THEN
  ! if we call it correctly, we may save here work between different RK-stages
  ! because of emmision and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.)
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.)
END IF

ImplicitSource=0.
#ifndef PP_HDG
CALL CalcSource(tStage,1.,ImplicitSource)
#else
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL CalcSourceHDG(i,j,k,iElem,ExplicitSource(1:PP_nVar,i,j,k,iElem))
  END DO; END DO; END DO !i,j,k    
END DO !iElem 
#endif
IF(DoVerifyCharge) CALL VerifyDepositedCharge()

IF(t.GE.DelayTime)THEN
  ! velocity to impulse
  CALL PartVeloToImp(VeloToImp=.TRUE.) 
  PartStateN(1:PDM%ParticleVecLength,1:6)=PartState(1:PDM%ParticleVecLength,1:6)
  IF(iter.EQ.0)THEN ! caution with emission: fields should also be interpolated to new particles, this is missing
                    ! or should be done directly during emission...
    ! should be already be done
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(PartIsImplicit(iPart))THEN
        CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
        SELECT CASE(PartLorentzType)
        CASE(0)
          Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(1)
          Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(3)
          Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(5)
          Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE DEFAULT
        END SELECT
      END IF ! ParticleIsImplicit
      PDM%IsNewPart(iPart)=.FALSE.
    END DO ! iPart
  END IF
  IF(DoSurfaceFlux)THEN
    DO iSpec=1,nSpecies
      IF (Species(iSpec)%nSurfacefluxBCs .EQ. 0) CYCLE
      DO iSF=1,Species(iSpec)%nSurfacefluxBCs
        Species(iSpec)%Surfaceflux(iSF)%tmpInsertedParticle        = Species(iSpec)%Surfaceflux(iSF)%InsertedParticle 
        Species(iSpec)%Surfaceflux(iSF)%tmpInsertedParticleSurplus = Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus
      END DO
    END DO ! iSpec
  END IF
END IF
RKdtFracTotal=0.
RKdtFrac     =0.
#endif /*PARTICLES*/

#ifndef PP_HDG
IF(iter.EQ.0) CALL DGTimeDerivative_weakForm(t, t, 0,doSource=.FALSE.)
#endif /*DG*/

! ----------------------------------------------------------------------------------------------------------------------------------
! stage 2 to 6
! ----------------------------------------------------------------------------------------------------------------------------------
DO iStage=2,nRKStages
  ! time of current stage
  tStage = t + RK_c(iStage)*dt
  alpha  = ESDIRK_a(iStage,iStage)*dt
  sGamma = 1.0/alpha
  IF(MPIRoot)THEN
    IF(DoPrintConvInfo)THEN
      WRITE(*,*) '-----------------------------'
      WRITE(*,*) 'istage:',istage
    END IF
  END IF
#ifndef PP_HDG
  ! store predictor
  CALL StorePredictor()

  ! compute the f(u^s-1)
  ! compute the f(u^s-1)
  ! DG-solver, Maxwell's equations
  FieldStage (:,:,:,:,:,iStage-1) = Ut
  FieldSource(:,:,:,:,:,iStage-1) = ImplicitSource
#endif /*DG*/

  ! and particles
#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    ! normal
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(PartIsImplicit(iPart))THEN
        IF(PartLorentzType.NE.5)THEN
          PartStage(iPart,1:3,iStage-1) = PartState(iPart,4:6) 
          PartStage(iPart,4:6,iStage-1) = Pt       (iPart,1:3)
        ELSE
          LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
          LorentzFacInv=1.0/SQRT(LorentzFacInv)
          PartStage(iPart,1  ,iStage-1) = PartState(iPart,4  ) * LorentzFacInv
          PartStage(iPart,2  ,iStage-1) = PartState(iPart,5  ) * LorentzFacInv
          PartStage(iPart,3  ,iStage-1) = PartState(iPart,6  ) * LorentzFacInv
          PartStage(iPart,4:6,iStage-1) = Pt       (iPart,1:3)
        END IF
      ELSE
        CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
        SELECT CASE(PartLorentzType)
        CASE(0)
          Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(1)
          Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
          LorentzFacInv = 1.0
        CASE(3)
          Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
          LorentzFacInv = 1.0
        CASE(5)
          LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
          LorentzFacInv=1.0/SQRT(LorentzFacInv)
          Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
        CASE DEFAULT
        END SELECT
        PartStage(iPart,1,iStage-1) = PartState(iPart,4)*LorentzFacInv
        PartStage(iPart,2,iStage-1) = PartState(iPart,5)*LorentzFacInv
        PartStage(iPart,3,iStage-1) = PartState(iPart,6)*LorentzFacInv
        PartStage(iPart,4:6,iStage-1) = Pt       (iPart,1:3)
      END IF ! ParticleIsImplicit
    END DO ! iPart
    ! new method
  END IF
#endif /*PARTICLES*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! explicit - particle  pusher
  !--------------------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
  ExplicitSource=0.
  ! particle step || only explicit particles
  IF (t.GE.DelayTime) THEN
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(.NOT.PartIsImplicit(iPart))THEN
        LastPartPos(iPart,1)=PartState(iPart,1)
        LastPartPos(iPart,2)=PartState(iPart,2)
        LastPartPos(iPart,3)=PartState(iPart,3)
        PEM%lastElement(iPart)=PEM%Element(iPart)
        ! compute explicit push
        PartState(iPart,1:6) = ERK_a(iStage,iStage-1)*PartStage(iPart,1:6,iStage-1)
        DO iCounter=1,iStage-2
          PartState(iPart,1:6)=PartState(iPart,1:6)+ERK_a(iStage,iCounter)*PartStage(iPart,1:6,iCounter)
        END DO ! iCounter=1,iStage-2
        PartState(iPart,1:6)=PartStateN(iPart,1:6)+dt*PartState(iPart,1:6)
      END IF ! ParticleIsExplicit
    END DO ! iPart

    !------------------------------------------------------------------------------------------------------------------------------
    ! particle surface flux || inflow boundary condition
    !------------------------------------------------------------------------------------------------------------------------------
    IF(DoSurfaceFlux)THEN
      ! insert the paritcles during each RK stage, later removed 
      RKdtFrac      = RK_c(iStage)
      RKdtFracTotal = RK_c(iStage)

      CALL ParticleSurfaceflux() 
      ! PO: normal RK: the timefrac is part of the total RKdtFrac because during this step, a certain number of particles is 
      !     emitted
      dtFrac= dt*RKdtFrac ! added total instead of normal....?
      ! now push the particles by their 
      !SF, new in current RKStage (no forces assumed in this stage)
      DO iPart=1,PDM%ParticleVecLength
        IF (PDM%IsNewPart(iPart)) THEN
          CALL RANDOM_NUMBER(RandVal)
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac*RandVal
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac*RandVal
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac*RandVal
          PartIsImplicit(iPart)=.FALSE.
          IF(PartLorentzType.EQ.5)THEN
            ! lorentztype
            LorentzFac=1.0-DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv
            LorentzFac=1.0/SQRT(LorentzFac)
            PartState(iPart,4) = LorentzFac*PartState(iPart,4)
            PartState(iPart,5) = LorentzFac*PartState(iPart,5)
            PartState(iPart,6) = LorentzFac*PartState(iPart,6)
          END IF
        END IF ! ParticleIsNew
      END DO ! iPart
    END IF ! DoSurfaceFlux
      
#ifdef MPI
    ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
  ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
#endif /*MPI*/
    IF(DoRefMapping)THEN
      ! tracking routines has to be extended for optional flag, like deposition
      CALL ParticleRefTracking(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    ELSE
      CALL ParticleTracing(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    END IF
#ifdef MPI
    ! send number of particles
    CALL SendNbOfParticles(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
#endif /*MPI*/
#ifdef MPI
    ! finish communication
    CALL MPIParticleRecv()
    ! set exchanged number of particles to zero
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    ! if new
    ! map particle from gamma*v to velocity
    CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    ! deposit explicit, local particles
    CALL Deposition(doInnerParts=.TRUE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
    CALL Deposition(doInnerParts=.FALSE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength)) ! external particles arg
    !PartMPIExchange%nMPIParticles=0
#ifndef PP_HDG
    CALL CalcSource(tStage,1.,ExplicitSource)
#else
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        CALL CalcSourceHDG(i,j,k,iElem,ExplicitSource(1:PP_nVar,i,j,k,iElem))
      END DO; END DO; END DO !i,j,k    
    END DO !iElem 
#endif
    ! map particle from v to gamma*v
    CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))

    ! delete surface flux particle
    IF(DoSurfaceFlux)THEN
      IF(iStage.NE.nRKStages)THEN
        DO iPart=1,PDM%ParticleVecLength
          IF (PDM%IsNewPart(iPart)) THEN
            PDM%IsNewPart(iPart)=.FALSE.
            PDM%ParticleInside(iPart)=.FALSE.
          END IF ! ParticleIsNew
        END DO ! iPart
        DO iSpec=1,nSpecies
          IF (Species(iSpec)%nSurfacefluxBCs .EQ. 0) CYCLE
          DO iSF=1,Species(iSpec)%nSurfacefluxBCs
            Species(iSpec)%Surfaceflux(iSF)%InsertedParticle        = Species(iSpec)%Surfaceflux(iSF)%tmpInsertedParticle 
            Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = Species(iSpec)%Surfaceflux(iSF)%tmpInsertedParticleSurplus
          END DO
        END DO ! iSpec
      END IF
    END IF ! DoSurfaceFlux
  END IF
#endif /*PARTICLES*/

  !--------------------------------------------------------------------------------------------------------------------------------
  ! implicit - particle pusher & Maxwell's field
  !--------------------------------------------------------------------------------------------------------------------------------

#ifndef PP_HDG
  ! compute RHS for linear solver
  LinSolverRHS=ESDIRK_a(iStage,iStage-1)*(FieldStage(:,:,:,:,:,iStage-1)+FieldSource(:,:,:,:,:,iStage-1))
  DO iCounter=1,iStage-2
    LinSolverRHS=LinSolverRHS +ESDIRK_a(iStage,iCounter)*(FieldStage(:,:,:,:,:,iCounter)+FieldSource(:,:,:,:,:,iCounter))
  END DO ! iCoutner=1,iStage-2
  LinSolverRHS=Un+dt*LinSolverRHS

  ! get predictor of u^s+1
  CALL Predictor(iStage,dt,Un,FieldSource,FieldStage) ! sets new value for U_DG
#endif /*DG*/

#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    DO iPart=1,PDM%ParticleVecLength
      IF(PartIsImplicit(iPart))THEN
        ! old position of stage
        ! StagePartPos(iPart,1)=PartState(iPart,1)
        ! StagePartPos(iPart,2)=PartState(iPart,2)
        ! StagePartPos(iPart,3)=PartState(iPart,3)
        ! PEM%StageElement(iPart)=PEM%Element(iPart)
        LastPartPos(iPart,1)=PartState(iPart,1)
        LastPartPos(iPart,2)=PartState(iPart,2)
        LastPartPos(iPart,3)=PartState(iPart,3)
        PEM%lastElement(iPart)=PEM%Element(iPart)
        ! compute Q and U
        PartQ(1:6,iPart) = ESDIRK_a(iStage,iStage-1)*PartStage(iPart,1:6,iStage-1)
        DO iCounter=1,iStage-2
          PartQ(1:6,iPart) = PartQ(1:6,iPart) + ESDIRK_a(iStage,iCounter)*PartStage(iPart,1:6,iCounter)
        END DO ! iCounter=1,iStage-2
        PartQ(1:6,iPart) = PartStateN(iPart,1:6) + dt* PartQ(1:6,iPart)
        ! predictor
        CALL PartPredictor(iStage,dt,iPart) ! sets new value for U_DG
      END IF ! PartIsImplicit
    END DO ! iPart
    IF(PredictorType.GT.0)THEN
#ifdef MPI
      ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
      ! open receive buffer for number of particles
      CALL IRecvNbofParticles()
#endif /*MPI*/
      IF(DoRefMapping)THEN
        ! tracking routines has to be extended for optional flag, like deposition
        !CALL ParticleRefTracking()
        CALL ParticleRefTracking(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      ELSE
        !CALL ParticleTracing()
        CALL ParticleTracing(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      END IF
#ifdef MPI
      ! send number of particles
      CALL SendNbOfParticles(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
      !CALL SendNbOfParticles() ! all particles to get initial deposition right \\ without emmission
      ! finish communication of number of particles and send particles
      CALL MPIParticleSend()
      ! finish communication
      CALL MPIParticleRecv()
      PartMPIExchange%nMPIParticles=0
      IF(DoExternalParts)THEN
        ! as we do not have the shape function here, we have to deallocate something
        SDEALLOCATE(ExtPartState)
        SDEALLOCATE(ExtPartSpecies)
        SDEALLOCATE(ExtPartToFIBGM)
        SDEALLOCATE(ExtPartMPF)
      END IF
#endif
      DO iPart=1,PDM%ParticleVecLength
        IF(PartIsImplicit(iPart))THEN
          IF(.NOT.PDM%ParticleInside(iPart)) PartIsImplicit(iPart)=.FALSE.
        END IF ! PartIsImplicit
      END DO ! iPart
    END IF
  END IF
#endif /*PARTICLES*/
  ! full newton for particles and fields
  CALL FullNewton(t,tStage,alpha)
#ifndef PP_HDG
  CALL DivCleaningDamping()
#endif /*DG*/
#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    IF(DoUpdateInStage)THEN
      CALL UpdateNextFreePosition()
    END IF
    ! surface flux, select if implicit or explicit particle treatment
  END IF
#endif /*PARTICLES*/
END DO

#ifdef PARTICLES
! particle step || only explicit particles
IF (t.GE.DelayTime) THEN
  ! add here new method
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart))CYCLE
    IF(DoSurfaceFlux)THEN
      IF(PDM%IsNewPart(iPart))THEN
        PDM%IsNewPart(iPart)=.FALSE.
        PEM%lastElement(iPart)=PEM%Element(iPart)
        LastPartPos(iPart,1)=PartState(iPart,1)
        LastPartPos(iPart,2)=PartState(iPart,2)
        LastPartPos(iPart,3)=PartState(iPart,3)
        LorentzFac=1.0-DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv
        LorentzFac=1.0/SQRT(LorentzFac)
        PartState(iPart,4) = LorentzFac*PartState(iPart,4)
        PartState(iPart,5) = LorentzFac*PartState(iPart,5)
        PartState(iPart,6) = LorentzFac*PartState(iPart,6)
        CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
        SELECT CASE(PartLorentzType)
        CASE(0)
          Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(1)
          Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(3)
          Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE(5)
          Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        CASE DEFAULT
        END SELECT
        CYCLE
      END IF
    END IF
    IF(.NOT.PartIsImplicit(iPart))THEN
      LastPartPos(iPart,1)=PartState(iPart,1)
      LastPartPos(iPart,2)=PartState(iPart,2)
      LastPartPos(iPart,3)=PartState(iPart,3)
      PEM%lastElement(iPart)=PEM%Element(iPart)
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
      SELECT CASE(PartLorentzType)
      CASE(0)
        Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      CASE(1)
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(3)
        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(5)
        LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
        LorentzFacInv=1.0/SQRT(LorentzFacInv)
        Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
      CASE DEFAULT
      END SELECT
      PartState(iPart,1  ) = RK_b(nRKStages)*LorentzFacInv*PartState(iPart,4)
      PartState(iPart,2  ) = RK_b(nRKStages)*LorentzFacInv*PartState(iPart,5)
      PartState(iPart,3  ) = RK_b(nRKStages)*LorentzFacInv*PartState(iPart,6)
      PartState(iPart,4:6) = RK_b(nRKSTages)*Pt(iPart,1:3)
      !  stage 1 ,nRKStages-1
      DO iCounter=1,nRKStages-1
        PartState(iPart,1:6) = PartState(iPart,1:6)   &
                             + RK_b(iCounter)*PartStage(iPart,1:6,iCounter)
      END DO ! counter
      PartState(iPart,1:6) = PartStateN(iPart,1:6)+dt*PartState(iPart,1:6)

    END IF ! ParticleIsExplicit
  END DO ! iPart

!  ! do the same stuff with the implicit pushed particles
!  DO iPart=1,PDM%ParticleVecLength
!    IF(.NOT.PDM%ParticleInside(iPart))CYCLE
!    IF(PartIsImplicit(iPart))THEN
!      LastPartPos(iPart,1)=PartState(iPart,1)
!      LastPartPos(iPart,2)=PartState(iPart,2)
!      LastPartPos(iPart,3)=PartState(iPart,3)
!      PEM%lastElement(iPart)=PEM%Element(iPart)
!      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
!      SELECT CASE(PartLorentzType)
!      CASE(1)
!        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
!      CASE(3)
!        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
!      CASE DEFAULT
!      END SELECT
!      PartState(iPart,1:3) = RK_b(nRKStages)*PartState(iPart,4:6)
!      PartState(iPart,4:6) = RK_b(nRKSTages)*Pt(iPart,1:3)
!      !  stage 1 ,nRKStages-1
!      DO iCounter=1,nRKStages-1
!        PartState(iPart,1:6) = PartState(iPart,1:6)   &
!                             + RK_b(iCounter)*PartStage(iPart,1:6,iCounter)
!      END DO ! counter
!      PartState(iPart,1:6) = PartStateN(iPart,1:6)+dt*PartState(iPart,1:6)
!
!    END IF ! ParticleIsImplicit
!  END DO ! iPart

  iStage=0  
#ifdef MPI
  ! mpi-routines should be extended by additional input: PartisImplicit, better criterion, saves computational time
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#endif /*MPI*/
  IF(DoRefMapping)THEN
    ! tracking routines has to be extended for optional flag, like deposition
    CALL ParticleRefTracking()
    !CALL ParticleRefTracking(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
  ELSE
    CALL ParticleTracing()
    !CALL ParticleTracing(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
  END IF
#ifdef MPI
  ! send number of particles
  CALL SendNbOfParticles(doParticle_In=.NOT.PartIsImplicit(1:PDM%ParticleVecLength))
  !CALL SendNbOfParticles() ! all particles to get initial deposition right \\ without emmission
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
! #endif ! old -> new is 9 lines below
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
  END IF
#endif

  ! map particle from gamma*v to v
  CALL PartVeloToImp(VeloToImp=.FALSE.) 

END IF
#endif /*PARTICLES*/

!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
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

!----------------------------------------------------------------------------------------------------------------------------------
! split and merge
!----------------------------------------------------------------------------------------------------------------------------------

IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
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
#endif /*PARTICLES*/

END SUBROUTINE TimeStepByImplicitRK
#endif /*PP_TimeDiscMethod==120 || PP_TimeDiscMethod==121 || PP_TimeDiscMethod==122  */


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
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking
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
    CALL ParticleTracing()
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
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
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
    CALL ParticleTracing()
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

#ifdef PP_HDG
#if (PP_TimeDiscMethod==500)
SUBROUTINE TimeStepPoisson(t)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,                 ONLY: U
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,iter
USE MOD_HDG,                     ONLY: HDG
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping!,MeasureTrackTime
#ifdef PARTICLES
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartState, Pt, LastPartPos,PEM, PDM, usevMPF, doParticleMerge, DelayTime, PartPressureCell
!USE MOD_Particle_Vars,           ONLY : Time
USE MOD_Particle_Vars,           ONLY: DoSurfaceFlux
USE MOD_part_RHS,                ONLY: CalcPartRHS
!USE MOD_part_boundary,           ONLY : ParticleBoundary
USE MOD_part_emission,           ONLY: ParticleInserting, ParticleSurfaceflux
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,           ONLY: StartParticleMerge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars,      ONLY:  DoExternalParts
USE MOD_Particle_MPI_Vars,       ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM
#endif
!USE MOD_PIC_Analyze,      ONLY: CalcDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping!,MeasureTrackTime
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleCollectCharges
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart
REAL    :: RandVal, dtFrac
!===================================================================================================================================
!Time = t

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
  END IF
#endif /*MPI*/
  CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
  IF(DoExternalParts)THEN
    ! finish communication
    CALL MPIParticleRecv()
  END IF
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF

! EM field
CALL HDG(t,U,iter)

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
  CALL CalcPartRHS()
END IF

IF (DoSurfaceFlux) THEN
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  IF (t.GE.DelayTime) THEN
    CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
  END IF
  IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%dtFracPush(iPart)) THEN
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dt
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dt
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dt
          PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1) * dt
          PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2) * dt
          PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3) * dt
        ELSE
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RandVal
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
    END DO
  END IF
ELSE
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles
    PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + PartState(1:PDM%ParticleVecLength,4) * dt
    PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + PartState(1:PDM%ParticleVecLength,5) * dt
    PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + PartState(1:PDM%ParticleVecLength,6) * dt
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + Pt(1:PDM%ParticleVecLength,1) * dt
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + Pt(1:PDM%ParticleVecLength,2) * dt
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + Pt(1:PDM%ParticleVecLength,3) * dt
  END IF
END IF

IF (t.GE.DelayTime) THEN
#ifdef MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTracing()
  END IF
  CALL ParticleInserting()
#ifdef MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
#endif
  CALL ParticleCollectCharges()
END IF


IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
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
END SUBROUTINE TimeStepPoisson
#endif /*(PP_TimeDiscMethod==500)*/

#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
SUBROUTINE TimeStepPoissonByLSERK(t,iter,tEndDiff)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY: Abort
USE MOD_PreProc
USE MOD_Analyze,               ONLY: PerformAnalyze
USE MOD_TimeDisc_Vars,         ONLY: dt,iStage,RKdtFrac,RKdtFracTotal
USE MOD_TimeDisc_Vars,         ONLY: RK_a,RK_b,RK_c,nRKStages
USE MOD_DG_Vars,               ONLY: U
USE MOD_Particle_Tracking_vars,ONLY: DoRefMapping
#ifdef PARTICLES
USE MOD_PICDepo,               ONLY: Deposition
USE MOD_PICInterpolation,      ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,         ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime,  PEM, PDM, & 
                                     doParticleMerge,PartPressureCell,DoSurfaceFlux
USE MOD_part_RHS,              ONLY: CalcPartRHS
USE MOD_part_emission,         ONLY: ParticleInserting, ParticleSurfaceflux
USE MOD_DSMC,                  ONLY: DSMC_main
USE MOD_DSMC_Vars,             ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,         ONLY: StartParticleMerge
USE MOD_PIC_Analyze,           ONLY: VerifyDepositedCharge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars,      ONLY:  DoExternalParts
USE MOD_Particle_MPI_Vars,       ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM
#endif
USE MOD_Particle_Tracking_vars, ONLY: DoRefMapping
USE MOD_part_tools,             ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking,      ONLY: ParticleTracing,ParticleRefTracking,ParticleCollectCharges
#endif /*PARTICLES*/
USE MOD_HDG           ,ONLY: HDG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: tendDiff
REAL,INTENT(INOUT)            :: t
INTEGER(KIND=8),INTENT(INOUT) :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL           :: tStage,b_dt(1:nRKStages)
REAL           :: Pa_rebuilt_coeff(1:nRKStages),Pa_rebuilt(1:3,1:nRKStages),Pv_rebuilt(1:3,1:nRKStages),v_rebuilt(1:3,0:nRKStages-1)
INTEGER        :: iPart, iStage_loc
REAL           :: RandVal
!===================================================================================================================================

DO iStage_loc=1,nRKStages
  ! RK coefficients
  b_dt(iStage_loc)=RK_b(iStage_loc)*dt
  ! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RK_a(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO
iStage=1

! first RK step
#ifdef PARTICLES
!Time=t
RKdtFrac = RK_c(2)
RKdtFracTotal=RKdtFrac

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! communicate shape function particles
#ifdef MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
  END IF
#endif /*MPI*/
  CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
  IF(DoExternalParts)THEN
    ! finish communication
    CALL MPIParticleRecv()
  END IF
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
#endif /*PARTICLES*/

CALL HDG(t,U,iter)

! calling the analyze routines
CALL PerformAnalyze(t,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)

#ifdef PARTICLES
! set last data already here, since surfaceflux moved before interpolation
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN
  IF (DoSurfaceFlux) THEN
    CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
  END IF
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
  !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
  CALL CalcPartRHS()
END IF

! particles
IF (t.GE.DelayTime) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Pt is not known only for new Surfaceflux-Parts -> change IsNewPart back to F for other Parts
      IF (.NOT.DoSurfaceFlux) THEN
        PDM%IsNewPart(iPart)=.FALSE.
      ELSE
        IF (.NOT.PDM%dtFracPush(iPart)) PDM%IsNewPart(iPart)=.FALSE.
      END IF
      !-- Particle Push
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(iPart,1) = PartState(iPart,4)
        Pt_temp(iPart,2) = PartState(iPart,5)
        Pt_temp(iPart,3) = PartState(iPart,6)
        Pt_temp(iPart,4) = Pt(iPart,1)
        Pt_temp(iPart,5) = Pt(iPart,2)
        Pt_temp(iPart,6) = Pt(iPart,3)
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4)*b_dt(1)
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5)*b_dt(1)
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6)*b_dt(1)
        PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1)*b_dt(1)
        PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2)*b_dt(1)
        PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3)*b_dt(1)
      ELSE !IsNewPart: no Pt_temp history available!
        IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage
          CALL RANDOM_NUMBER(RandVal)
        ELSE
          CALL abort(&
__STAMP__&
,'Error in LSERK-HDG-Timedisc: This case should be impossible...')
        END IF
        Pa_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(iPart,1:3)
        END DO
        v_rebuilt(:,:)=0.
        DO iStage_loc=iStage-1,0,-1
          IF (iStage_loc.EQ.iStage-1) THEN
            v_rebuilt(1:3,iStage_loc) = PartState(iPart,4:6) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          ELSE
            v_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc+1) - b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          END IF
        END DO
        Pv_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          IF (iStage_loc.EQ.1) THEN
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,0)
          ELSE
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc-1) - RK_a(iStage_loc)*Pv_rebuilt(1:3,iStage_loc-1)
          END IF
        END DO
        Pt_temp(iPart,1:3) = Pv_rebuilt(1:3,iStage)
        Pt_temp(iPart,4:6) = Pa_rebuilt(1:3,iStage)
        PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)*RandVal
        PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)*RandVal
        PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)*RandVal
        PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)*RandVal
        PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)*RandVal
        PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)*RandVal
        PDM%dtFracPush(iPart) = .FALSE.
        PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
      END IF !IsNewPart
    END IF
  END DO
END IF

IF (t.GE.DelayTime) THEN ! removed .OR.(iter.EQ.0) because particles have not moved
#ifdef MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTracing()
  END IF
  CALL ParticleInserting()
#ifdef MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
  CALL ParticleCollectCharges()
END IF

#endif /*PARTICLES*/

! perform RK steps
DO iStage=2,nRKStages
  tStage=t+dt*RK_c(iStage)
#ifdef PARTICLES
  IF (iStage.NE.nRKStages) THEN
    RKdtFrac = RK_c(iStage+1)-RK_c(iStage)
    RKdtFracTotal=RKdtFracTotal+RKdtFrac
  ELSE
    RKdtFrac = 1.-RK_c(nRKStages)
    RKdtFracTotal=1.
  END IF

  ! deposition 
  IF (t.GE.DelayTime) THEN
    CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
    ! here: finish deposition with delta kernal
    !       maps source terms in physical space
    ! ALWAYS require
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
  END IF
#endif /*PARTICLES*/

  CALL HDG(tStage,U,iter)

#ifdef PARTICLES
  ! set last data already here, since surfaceflux moved before interpolation
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  IF (t.GE.DelayTime) THEN
    IF (DoSurfaceFlux) CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
    ! forces on particle
    CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
    !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
    CALL CalcPartRHS()

    ! particle step
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%IsNewPart(iPart)) THEN
          Pt_temp(iPart,1) = PartState(iPart,4) - RK_a(iStage) * Pt_temp(iPart,1)
          Pt_temp(iPart,2) = PartState(iPart,5) - RK_a(iStage) * Pt_temp(iPart,2)
          Pt_temp(iPart,3) = PartState(iPart,6) - RK_a(iStage) * Pt_temp(iPart,3)
          Pt_temp(iPart,4) = Pt(iPart,1) - RK_a(iStage) * Pt_temp(iPart,4)
          Pt_temp(iPart,5) = Pt(iPart,2) - RK_a(iStage) * Pt_temp(iPart,5)
          Pt_temp(iPart,6) = Pt(iPart,3) - RK_a(iStage) * Pt_temp(iPart,6)
          PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)
          PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)
          PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)
          PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)
          PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)
          PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)
        ELSE !IsNewPart: no Pt_temp history available!
          IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage
            CALL RANDOM_NUMBER(RandVal)
            PDM%dtFracPush(iPart) = .FALSE.
          ELSE !new but without SF in current RKStage (i.e., from ParticleInserting or diffusive wall reflection)
               ! -> rebuild Pt_tmp-coefficients assuming F=const. (value at last Pos) in previous stages
            RandVal=1. !"normal" particles (i.e. not from SurfFlux) are pushed with whole timestep!
          END IF
          Pa_rebuilt(:,:)=0.
          DO iStage_loc=1,iStage
            Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(iPart,1:3)
          END DO
          v_rebuilt(:,:)=0.
          DO iStage_loc=iStage-1,0,-1
            IF (iStage_loc.EQ.iStage-1) THEN
              v_rebuilt(1:3,iStage_loc) = PartState(iPart,4:6) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
            ELSE
              v_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc+1) - b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
            END IF
          END DO
          Pv_rebuilt(:,:)=0.
          DO iStage_loc=1,iStage
            IF (iStage_loc.EQ.1) THEN
              Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,0)
            ELSE
              Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc-1) - RK_a(iStage_loc)*Pv_rebuilt(1:3,iStage_loc-1)
            END IF
          END DO
          Pt_temp(iPart,1:3) = Pv_rebuilt(1:3,iStage)
          Pt_temp(iPart,4:6) = Pa_rebuilt(1:3,iStage)
          PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)*RandVal
          PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)*RandVal
          PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)*RandVal
          PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)*RandVal
          PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)*RandVal
          PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)*RandVal
          PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
        END IF !IsNewPart
      END IF
    END DO

    ! particle tracking
#ifdef MPI
    CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
    IF(DoRefMapping)THEN
      CALL ParticleRefTracking()
    ELSE
      CALL ParticleTracing()
    END IF
    CALL ParticleInserting()
#ifdef MPI
    CALL SendNbOfParticles() ! send number of particles
    CALL MPIParticleSend()   ! finish communication of number of particles and send particles
    CALL MPIParticleRecv()   ! finish communication
#endif
    CALL ParticleCollectCharges()
  END IF
#endif /*PARTICLES*/
END DO

#ifdef PARTICLES
#ifdef MPI
PartMPIExchange%nMPIParticles=0 ! and set number of received particles to zero for deposition
#endif
IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
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

END SUBROUTINE TimeStepPoissonByLSERK
#endif /*(PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)*/
#endif /*PP_HDG*/

#if (PP_TimeDiscMethod==1000)
SUBROUTINE TimeStep_LD(t)
!===================================================================================================================================
! Low Diffusion Method (Mirza 2013)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars,          ONLY:dt
#ifdef PARTICLES
USE MOD_Particle_Vars,          ONLY:PartState, LastPartPos, PDM,PEM !, Time
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

  !Time = t
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
USE MOD_Particle_Vars,    ONLY : PartState, LastPartPos,  PDM,PEM, WriteMacroVolumeValues!Time
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

  !Time = t
  LD_DSMC_RHS(1:PDM%ParticleVecLength,1) = 0
  LD_DSMC_RHS(1:PDM%ParticleVecLength,2) = 0
  LD_DSMC_RHS(1:PDM%ParticleVecLength,3) = 0
  CALL LD_DSMC_Indicate_DSMC_Particles
  CALL DSMC_main() ! first dsmc then ld due to RHS-calculation!
  CALL LD_main()
! ----- Start Analyze Particles
  IF (.NOT.WriteMacroVolumeValues) THEN
    IF(t.ge.(1-DSMC%TimeFracSamp)*TEnd) THEN
      CALL LD_DSMC_data_sampling()  ! Data sampling for output
      IF(DSMC%NumOutput.NE.0) THEN
        nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
        IF(t.ge.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
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
  __STAMP__&
  ,'Polynomial degree is to high!',PP_N,999.)
CFLScale=CFLScale*CFLScaleAlpha(PP_N)
#endif
#if (PP_TimeDiscMethod==6)
IF(PP_N.GT.15) CALL abort(&
  __STAMP__&
  ,'Polynomial degree is to high!',PP_N,999.)
CFLScale=CFLScale*CFLScaleAlpha(PP_N)
#endif
!scale with 2N+1
CFLScale = CFLScale/(2.*PP_N+1.)
SWRITE(UNIT_stdOut,'(A,ES16.7)') '   CFL:',CFLScale
END SUBROUTINE fillCFL_DFL


SUBROUTINE InitTimeStep()
!===================================================================================================================================
!initial time step calculations for new timedisc-loop
!===================================================================================================================================
! MODULES
USE MOD_Globals
#ifdef PARTICLES
USE MOD_PARTICLE_Vars,      ONLY: ManualTimeStep, useManualTimestep
#if (PP_TimeDiscMethod==200)
USE MOD_PARTICLE_Vars,      ONLY: dt_maxwell,dt_max_particles,MaxwellIterNum,NextTimeStepAdjustmentIter
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_TimeDisc_Vars,      ONLY: IterDisplayStep
#endif
#endif /*PARTICLES*/
#ifndef PP_HDG
USE MOD_CalcTimeStep,       ONLY:CalcTimeStep
USE MOD_TimeDisc_Vars,      ONLY: CFLtoOne
#endif
USE MOD_TimeDisc_Vars,      ONLY: dt, dt_Min
USE MOD_TimeDisc_Vars,      ONLY: sdtCFLOne
#if (PP_TimeDiscMethod==201)                                                                                                         
USE MOD_TimeDisc_Vars,      ONLY: dt_temp, MaximumIterNum 
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#ifdef PARTICLES
! For tEnd != tStart we have to advance the solution in time
IF(useManualTimeStep)THEN
  ! particle time step is given externally and not calculated through the solver
  dt_Min=ManualTimeStep
#if (PP_TimeDiscMethod==200)
  dt_max_particles = ManualTimeStep
  dt_maxwell = CalcTimeStep()
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
#ifndef PP_HDG
  dt_Min=CalcTimeStep()
  sdtCFLOne  = 1.0/(dt_Min*CFLtoOne)
#else
  dt_Min=0
  sdtCFLOne  = -1.0 !dummy for HDG!!!
#endif /*PP_HDG*/

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
END SUBROUTINE InitTimeStep



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
