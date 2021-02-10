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


MODULE MOD_TimeDisc
!===================================================================================================================================
! Module for the GTS Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeDisc
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeDisc()
!===================================================================================================================================
! GTS Temporal discretization
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: SimulationEfficiency,PID,WallTime
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: time,TEnd,dt,iter,IterDisplayStep,DoDisplayIter,dt_Min,tAnalyze,dtWeight,tEndDiff,tAnalyzeDiff
#if (PP_TimeDiscMethod==509)
USE MOD_TimeDisc_Vars          ,ONLY: dt_old
#endif /*(PP_TimeDiscMethod==509)*/
USE MOD_TimeAverage_vars       ,ONLY: doCalcTimeAverage
USE MOD_TimeAverage            ,ONLY: CalcTimeAverage
USE MOD_Analyze                ,ONLY: PerformAnalyze
USE MOD_Analyze_Vars           ,ONLY: Analyze_dt,iAnalyze
USE MOD_Restart_Vars           ,ONLY: RestartTime,RestartWallTime
USE MOD_HDF5_output            ,ONLY: WriteStateToHDF5
USE MOD_Mesh_Vars              ,ONLY: MeshFile,nGlobalElems
USE MOD_RecordPoints_Vars      ,ONLY: RP_onProc
USE MOD_RecordPoints           ,ONLY: WriteRPToHDF5!,RecordPoints
USE MOD_LoadBalance_Vars       ,ONLY: nSkipAnalyze
#if !(USE_HDG)
USE MOD_PML_Vars               ,ONLY: DoPML,DoPMLTimeRamp,PMLTimeRamp
USE MOD_PML                    ,ONLY: PMLTimeRamping
#if USE_LOADBALANCE
#ifdef maxwell
#if defined(ROS) || defined(IMPA)
USE MOD_Precond_Vars           ,ONLY:UpdatePrecondLB
#endif /*ROS or IMPA*/
#endif /*maxwell*/
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
#ifdef PP_POIS
USE MOD_Restart_Vars           ,ONLY: DoRestart
USE MOD_Equation               ,ONLY: EvalGradient
#endif /*PP_POIS*/
#if USE_MPI
#if USE_LOADBALANCE
USE MOD_LoadBalance            ,ONLY: LoadBalance,ComputeElemLoad
USE MOD_LoadBalance_Vars       ,ONLY: DoLoadBalance,ElemTime
USE MOD_LoadBalance_Vars       ,ONLY: LoadBalanceSample,PerformLBSample,PerformLoadBalance,LoadBalanceMaxSteps,nLoadBalanceSteps
USE MOD_Restart_Vars           ,ONLY: DoInitialAutoRestart,InitialAutoRestartSample,IAR_PerformPartWeightLB
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimeField
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#ifdef PARTICLES
USE MOD_Particle_Vars          ,ONLY: WriteMacroVolumeValues, WriteMacroSurfaceValues, MacroValSampTime
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_HDF5_Output_Tools      ,ONLY: WriteIMDStateToHDF5
#endif /*PARTICLES*/
#ifdef PARTICLES
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_Particle_Vars          ,ONLY: DoImportIMDFile
#if USE_MPI
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#endif /*USE_MPI*/
USE MOD_Particle_Vars          ,ONLY: doParticleMerge, enableParticleMerge, vMPFMergeParticleIter, UseAdaptive
USE MOD_Particle_Tracking_vars ,ONLY: tTracking,tLocalization,nTracks,MeasureTrackTime
#if (USE_MPI) && (USE_LOADBALANCE) && defined(PARTICLES)
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimePart
#endif /* USE_LOADBALANCE && PARTICLES*/
USE MOD_Part_Emission          ,ONLY: AdaptiveBCAnalyze
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation       ,ONLY: InitAnalyticalParticleState
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/
USE MOD_Output                 ,ONLY: PrintStatusLine
USE MOD_TimeStep
USE MOD_TimeDiscInit           ,ONLY: InitTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: tStart                   !> simulation time at the beginning of the simulation
REAL                         :: tPreviousAnalyze         !> time of previous analyze.
                                                         !> Used for Nextfile info written into previous file if greater tAnalyze
REAL                         :: tPreviousAverageAnalyze  !> time of previous Average analyze.
REAL                         :: tZero
INTEGER(KIND=8)              :: iter_PID                 !> iteration counter since last InitPiclas call for PID calculation
REAL                         :: WallTimeStart            !> wall time of simulation start
REAL                         :: WallTimeEnd              !> wall time of simulation end
#if USE_LOADBALANCE
INTEGER                      :: tmp_LoadBalanceSample    !> loadbalance sample saved until initial autorestart ist finished
LOGICAL                      :: tmp_DoLoadBalance        !> loadbalance flag saved until initial autorestart ist finished
#endif /*USE_LOADBALANCE*/
LOGICAL                      :: finalIter
!===================================================================================================================================
tPreviousAnalyze=RestartTime
! first average analyze is not written at start but at first tAnalyze
tPreviousAverageAnalyze=tAnalyze
! saving the start of the simulation as restart time is overwritten during load balance step
!   In case the overwritten one is used, state write out is performed only next Nth analze-dt after restart instead after analyze-dt
!   w/o  tZero: nSkipAnalyze=5 , restart after iAnalyze=2 , next write out witout any restarts after iAnalyze=7
!   with tZero: nSkipAnalyze=5 , restart after iAnalyze=2 , next write out witout any restarts after iAnalyze=5
tZero = RestartTime

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems*(PP_N+1)**3)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs/Proc : ',REAL(nGlobalElems*(PP_N+1)**3/nProcessors)

!Evaluate Gradients to get Potential in case of Restart and Poisson Calc
#ifdef PP_POIS
IF(DoRestart) CALL EvalGradient()
#endif /*PP_POIS*/
! Write the state at time=0, i.e. the initial condition

#if defined(PARTICLES) && (USE_MPI)
! e.g. 'shape_function'
IF(StringBeginsWith(DepositionType,'shape_function'))THEN
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
!#endif /*USE_MPI*/
!  CALL Deposition()
!#if USE_MPI
!  CALL MPIParticleRecv()
!  ! second buffer
!  CALL Deposition()
!#endif /*USE_MPI*/
!#endif

tStart = time
CALL PrintStatusLine(time,dt,tStart,tEnd)
CALL InitTimeStep() ! Initial time step calculation
WallTimeStart=PICLASTIME()
iter=0
iter_PID=0

! fill recordpoints buffer (first iteration)
!IF(RP_onProc) CALL RecordPoints(iter,t,forceSampling=.TRUE.)

dt=MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/)) ! quick fix: set dt for initial write DSMCHOState (WriteMacroVolumeValues=T)

#if USE_LOADBALANCE
IF (DoInitialAutoRestart) THEN
  tmp_DoLoadBalance     = DoLoadBalance
  DoLoadBalance         = .TRUE.
  tmp_LoadbalanceSample = LoadBalanceSample
  LoadBalanceSample     = InitialAutoRestartSample
  ! correct initialautrestartSample if partweight_initialautorestart is enabled so tAnalyze is calculated correctly
  ! LoadBalanceSample still needs to be zero
  IF (IAR_PerformPartWeightLB) InitialAutoRestartSample=1
  ! correction for first analyzetime due to auto initial restart
  IF (MIN(RestartTime+iAnalyze*Analyze_dt,tEnd,RestartTime+InitialAutoRestartSample*dt).LT.tAnalyze) THEN
    tAnalyze     = MIN(RestartTime+iAnalyze*Analyze_dt,tEnd,RestartTime+InitialAutoRestartSample*dt)
    tAnalyzeDiff = tAnalyze-time
    dt           = MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/))
  END IF
END IF
#endif /*USE_LOADBALANCE*/

#ifdef PARTICLES
#ifdef CODE_ANALYZE
! Set specific particle position and velocity (calculated from an analytical expression)
CALL InitAnalyticalParticleState()
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/

CALL PerformAnalyze(time,FirstOrLastIter=.TRUE.,OutPutHDF5=.FALSE.)

#ifdef PARTICLES
IF(DoImportIMDFile) CALL WriteIMDStateToHDF5() ! write IMD particles to state file (and TTM if it exists)
#endif /*PARTICLES*/
! Write initial state to file
CALL WriteStateToHDF5(TRIM(MeshFile),time,tPreviousAnalyze)

! if measurement of particle tracking time (used for analyze, load balancing uses own time measurement for tracking)
#ifdef PARTICLES
IF(MeasureTrackTime)THEN
  nTracks=0
  tTracking=0
  tLocalization=0
END IF
#endif /*PARTICLES*/

! No computation needed if tEnd=tStart!
IF(time.EQ.tEnd)RETURN

!-----------------------------------------------------------------------------------------------------------------------------------
! iterations starting up from here
!-----------------------------------------------------------------------------------------------------------------------------------
DO !iter_t=0,MaxIter

#ifdef PARTICLES
  IF(enableParticleMerge) THEN
    IF ((iter.GT.0).AND.(MOD(iter,INT(vMPFMergeParticleIter,8)).EQ.0)) doParticleMerge=.true.
  END IF
#endif /*PARTICLES*/

  tAnalyzeDiff=tAnalyze-time    ! time to next analysis, put in extra variable so number does not change due to numerical errors
  tEndDiff=tEnd-time            ! dito for end time

  !IF(time.LT.3e-8)THEN
  !    !RETURN
  !ELSE
  !  IF(time.GT.4e-8) dt_Min=MIN(dt_Min*1.2,2e-8)
  !END IF
#if (PP_TimeDiscMethod==509)
  IF (iter.GT.0) THEN
    dt_old=dt
  END IF
#endif /*(PP_TimeDiscMethod==509)*/
  dt=MINVAL((/dt_Min,tAnalyzeDiff,tEndDiff/))
  dtWeight=dt/dt_Min !might be further descreased by rk-stages
#if (PP_TimeDiscMethod==509)
  IF (iter.EQ.0) THEN
    dt_old=dt
  ELSE IF (ABS(dt-dt_old).GT.1.0E-6*dt_old) THEN
    SWRITE(UNIT_StdOut,'(A,G0)')'WARNING: dt changed from last iter by a relative difference of ',(dt-dt_old)/dt_old
  END IF
#endif /*(PP_TimeDiscMethod==509)*/
#if USE_LOADBALANCE
  ! check if loadbalancing is enabled with elemtime calculation and only LoadBalanceSample number of iteration left until analyze
  ! --> set PerformLBSample true
  IF ((tAnalyzeDiff.LE.LoadBalanceSample*dt &                                 ! all iterations in LoadbalanceSample interval
      .OR. (ALMOSTEQUALRELATIVE(tAnalyzeDiff,LoadBalanceSample*dt,1e-5))) &   ! make sure to get the first iteration in interval
      .AND. .NOT.PerformLBSample .AND. DoLoadBalance) PerformLBSample=.TRUE.  ! make sure Loadbalancing is enabled
#endif /*USE_LOADBALANCE*/
  IF (tAnalyzeDiff-dt.LT.dt/100.0) dt = tAnalyzeDiff
  IF (tEndDiff-dt.LT.dt/100.0) dt = tEndDiff
  IF ( dt .LT. 0. ) THEN
    SWRITE(UNIT_StdOut,*)'*** ERROR: Is something wrong with the defined tEnd?!? ***'
    CALL abort(&
    __STAMP__&
    ,'Error in tEndDiff or tAnalyzeDiff!')

  END IF

  IF(doCalcTimeAverage) CALL CalcTimeAverage(.FALSE.,dt,time,tPreviousAverageAnalyze) ! tPreviousAnalyze not used if finalize_flag=false

#if !(USE_HDG)
  IF(DoPML)THEN
    IF(DoPMLTimeRamp)THEN
      CALL PMLTimeRamping(time,PMLTimeRamp)
    END IF
  END IF
#endif /*NOT USE_HDG*/

  CALL PrintStatusLine(time,dt,tStart,tEnd)

! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
#if (PP_TimeDiscMethod==1)
  CALL TimeStepByLSERK()
#elif (PP_TimeDiscMethod==2)
  CALL TimeStepByLSERK()
#elif (PP_TimeDiscMethod==3)
  CALL TimeStepByTAYLOR()
#elif (PP_TimeDiscMethod==4)
  CALL TimeStep_DSMC()
#elif (PP_TimeDiscMethod==6)
  CALL TimeStepByLSERK()
#elif (PP_TimeDiscMethod==42)
  CALL TimeStep_DSMC_Debug() ! Reservoir and Debug
#elif (PP_TimeDiscMethod==100)
  CALL TimeStepByEulerImplicit() ! O1 Euler Implicit
#elif (PP_TimeDiscMethod==120)
  CALL TimeStepByImplicitRK() ! O3 ERK/ESDIRK Particles + ESDIRK Field
#elif (PP_TimeDiscMethod==121)
  CALL TimeStepByImplicitRK() ! O3 ERK/ESDIRK Particles + ESDIRK Field
#elif (PP_TimeDiscMethod==122)
  CALL TimeStepByImplicitRK() ! O4 ERK/ESDIRK Particles + ESDIRK Field
#elif (PP_TimeDiscMethod==123)
  CALL TimeStepByImplicitRK() ! O3 ERK/ESDIRK Particles + ESDIRK Field
#elif (PP_TimeDiscMethod==130)
  CALL TimeStepByRosenbrock() ! linear Rosenbrock implicit
#elif (PP_TimeDiscMethod==131)
  CALL TimeStepByRosenbrock() ! linear Rosenbrock implicit
#elif (PP_TimeDiscMethod==132)
  CALL TimeStepByRosenbrock() ! linear Rosenbrock implicit
#elif (PP_TimeDiscMethod==133)
  CALL TimeStepByRosenbrock() ! linear Rosenbrock implicit
#elif (PP_TimeDiscMethod==134)
  CALL TimeStepByRosenbrock() ! linear Rosenbrock implicit
#elif (PP_TimeDiscMethod==300)
  CALL TimeStep_FPFlow()
#elif (PP_TimeDiscMethod==400)
  CALL TimeStep_BGK()
#elif (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
#if USE_HDG
#if (PP_TimeDiscMethod==500) || (PP_TimeDiscMethod==509)
  CALL TimeStepPoisson() ! Euler Explicit or leapfrog, Poisson
#elif (PP_TimeDiscMethod==508)
  CALL TimeStepPoissonByBorisLeapfrog() ! Boris-Leapfrog, Poisson
#else
  CALL TimeStepPoissonByLSERK() ! Runge Kutta Explicit, Poisson
#endif
#else
  CALL abort(&
  __STAMP__&
  ,'Timedisc 50x only available for EQNSYS Poisson!',PP_N,999.)
#endif /*USE_HDG*/
#endif
  ! calling the analyze routines
  iter=iter+1
  iter_PID=iter_PID+1
  time=time+dt
  IF(MPIroot) THEN
    IF(DoDisplayIter)THEN
      IF(MOD(iter,IterDisplayStep).EQ.0) THEN
         SWRITE(UNIT_stdOut,'(A,I21,A6,ES26.16E3,25X)')" iter:", iter,"time:",time ! new format for analyze time output
      END IF
    END IF
  END IF
  ! calling the analyze routines
  IF(ALMOSTEQUAL(dt,tEndDiff))THEN
    finalIter=.TRUE.
  ELSE
    finalIter=.FALSE.
  END IF
  CALL PerformAnalyze(time,FirstOrLastIter=finalIter,OutPutHDF5=.FALSE.)
#ifdef PARTICLES
  ! sampling of near adaptive boundary element values
  IF(UseAdaptive.OR.(nPorousBC.GT.0)) CALL AdaptiveBCAnalyze()
#endif /*PARICLES*/
  ! output of state file
  !IF ((dt.EQ.tAnalyzeDiff).OR.(dt.EQ.tEndDiff)) THEN   ! timestep is equal to time to analyze or end
  IF((ALMOSTEQUAL(dt,tAnalyzeDiff)).OR.(ALMOSTEQUAL(dt,tEndDiff)))THEN
    WallTimeEnd=PICLASTIME()
    IF(MPIroot)THEN ! determine the SimulationEfficiency and PID here,
                    ! because it is used in ComputeElemLoad -> WriteElemTimeStatistics
      WallTime = WallTimeEnd-StartTime
      SimulationEfficiency = (time-RestartTime)/((WallTimeEnd-RestartWallTime)*nProcessors/3600.) ! in [s] / [CPUh]
      PID=(WallTimeEnd-WallTimeStart)*nProcessors/(nGlobalElems*(PP_N+1)**3*iter_PID)
    END IF

#if USE_MPI
#ifdef PARTICLES
#if !defined(LSERK) && !defined(IMPA) && !defined(ROS)
    CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
#endif
#endif /*PARICLES*/
#if USE_LOADBALANCE
#ifdef PARTICLES
    ! Check if loadbalancing is enabled with partweight and set PerformLBSample true to calculate elemtimes with partweight
    ! LoadBalanceSample is 0 if partweightLB or IAR_partweighlb are enabled. If only one of them is set Loadbalancesample switches
    ! during time loop
    IF (LoadBalanceSample.EQ.0 .AND. DoLoadBalance .AND. .NOT.PerformLBSample) PerformLBSample=.TRUE.
#endif /*PARICLES*/
    ! routine calculates imbalance and if greater than threshold PerformLoadBalance=.TRUE.
    CALL ComputeElemLoad()
#ifdef maxwell
#if defined(ROS) || defined(IMPA)
    UpdatePrecondLB=PerformLoadBalance
#endif /*ROS or IMPA*/
#endif /*maxwell*/
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

#if USE_LOADBALANCE
    IF(MOD(iAnalyze,nSkipAnalyze).EQ.0 .OR. PerformLoadBalance .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#else
    IF( MOD(iAnalyze,nSkipAnalyze).EQ.0 .OR. ALMOSTEQUAL(dt,tEndDiff))THEN
#endif /*USE_LOADBALANCE*/
      ! Analyze for output
      CALL PerformAnalyze(tAnalyze,FirstOrLastIter=finalIter,OutPutHDF5=.TRUE.)
      ! write information out to std-out of console
      CALL WriteInfoStdOut()
      ! Write state to file
      CALL WriteStateToHDF5(TRIM(MeshFile),time,tPreviousAnalyze)
      IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,time,tPreviousAverageAnalyze)
      ! Write recordpoints data to hdf5
      IF(RP_onProc) CALL WriteRPtoHDF5(tAnalyze,.TRUE.)
      tPreviousAnalyze=tAnalyze
      tPreviousAverageAnalyze=tAnalyze
      SWRITE(UNIT_StdOut,'(132("-"))')
    END IF ! actual analyze is done
    iter_PID=0
#if USE_LOADBALANCE
    ! Check if load balancing must be performed
    IF(DoLoadBalance.AND.PerformLBSample.AND.(LoadBalanceMaxSteps.GT.nLoadBalanceSteps))THEN
      IF(time.LT.tEnd)THEN ! do not perform a load balance restart when the last timestep is performed
        IF(PerformLoadBalance) THEN
          ! DO NOT DELETE THIS: ONLY recalculate the timestep when the mesh is changed!
          !CALL InitTimeStep() ! re-calculate time step after load balance is performed
          RestartTime=time ! Set restart simulation time to current simulation time because the time is not read from the state file
          RestartWallTime=PICLASTIME() ! Set restart wall time if a load balance step is performed
          dtWeight=1. ! is intialized in InitTimeDisc which is not called in LoadBalance, but needed for restart (RestartHDG)
        END IF
        CALL LoadBalance()
        IF(PerformLoadBalance .AND. MOD(iAnalyze,nSkipAnalyze).NE.0) &
          CALL PerformAnalyze(time,FirstOrLastIter=.FALSE.,OutPutHDF5=.TRUE.)
        !      dt=dt_Min !not sure if nec., was here before InitTimeStep was created, overwritten in next iter anyway
        ! CALL WriteStateToHDF5(TRIM(MeshFile),time,tPreviousAnalyze) ! not sure if required
      END IF
    ELSE
      ElemTime=0. ! nullify ElemTime before measuring the time in the next cycle
#ifdef PARTICLES
      ElemTimePart    = 0.
#endif /*PARTICLES*/
      ElemTimeField    = 0.
    END IF
    PerformLBSample=.FALSE.
    IF (DoInitialAutoRestart) THEN
      DoInitialAutoRestart = .FALSE.
      DoLoadBalance = tmp_DoLoadBalance
      LoadBalanceSample = tmp_LoadBalanceSample
      iAnalyze=0
#ifdef PARTICLES
      DSMC%SampNum=0
      IF (WriteMacroVolumeValues .OR. WriteMacroSurfaceValues) MacroValSampTime = Time
#endif /*PARTICLES*/
    END IF
#endif /*USE_LOADBALANCE*/
    ! count analyze dts passed
    iAnalyze=iAnalyze+1
    tAnalyze=MIN(tZero+REAL(iAnalyze)*Analyze_dt,tEnd)
    WallTimeStart=PICLASTIME()
  END IF !dt_analyze
  IF(time.GE.tEnd)EXIT ! done, worst case: one additional time step
#ifdef PARTICLES
  ! Switch flag to false after the number of particles has been written to std out and before the time next step is started
  GlobalNbrOfParticlesUpdated = .FALSE.    
#endif /*PARTICLES*/
END DO ! iter_t
END SUBROUTINE TimeDisc


SUBROUTINE WriteInfoStdOut()
!===================================================================================================================================
!> writes information to console std_out
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: SimulationEfficiency,PID
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: iter,dt_Min
#if defined(IMPA) || defined(ROS)
USE MOD_LinearSolver_Vars      ,ONLY: totalIterLinearSolver
#endif /*IMPA || ROS*/
#ifdef PARTICLES
USE MOD_Particle_Tracking_vars ,ONLY: CountNbrOfLostParts,NbrOfLostParticles,NbrOfLostParticlesTotal,NbrOfLostParticlesTotal_old
#ifdef IMPA
USE MOD_LinearSolver_vars      ,ONLY: nPartNewton
USE MOD_LinearSolver_Vars      ,ONLY: totalFullNewtonIter
#endif /*IMPA*/
#if defined(IMPA) || defined(ROS)
USE MOD_LinearSolver_Vars      ,ONLY: TotalPartIterLinearSolver
#endif /*IMPA || ROS*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: TimeArray(8)              ! Array for system time
INTEGER             :: NbrOfLostParticlesTotal_old_tmp
!===================================================================================================================================
#ifdef PARTICLES
IF(CountNbrOfLostParts)THEN
#if USE_MPI
  NbrOfLostParticlesTotal_old_tmp = NbrOfLostParticlesTotal ! keep old value
  ! Allreduce is required because of the particle output to .h5 in which all processors must take place
  CALL MPI_ALLREDUCE(NbrOfLostParticles , NbrOfLostParticlesTotal , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_WORLD , IERROR)
  NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + NbrOfLostParticlesTotal_old_tmp ! add old value
#else
  NbrOfLostParticlesTotal=NbrOfLostParticles
#endif /*USE_MPI*/
END IF
#endif /*PARICLES*/
IF(MPIroot)THEN
  ! simulation time per CPUh efficiency in [s]/[CPUh]
  !SimulationEfficiency = (time-RestartTime)/((WallTimeEnd-StartTime)*nProcessors/3600.) ! in [s] / [CPUh]
  ! Get calculation time per DOF
  !PID=(WallTimeEnd-WallTimeStart)*nProcessors/(nGlobalElems*(PP_N+1)**3*iter_PID)
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
    ' Sys date  :    ',TimeArray(3),'.',TimeArray(2),'.',TimeArray(1),' ',TimeArray(5),':',TimeArray(6),':',TimeArray(7)
  WRITE(UNIT_stdOut,'(A,ES12.5,A)')' PID: CALCULATION TIME PER TSTEP/DOF: [',PID,' sec ]'
  WRITE(UNIT_stdOut,'(A,ES12.5,A)')' EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [',SimulationEfficiency,&
                                                                                        ' sec/h ]'
  WRITE(UNIT_StdOut,'(A,ES16.7)')' Timestep  : ',dt_Min
  WRITE(UNIT_stdOut,'(A,ES16.7)')'#Timesteps : ',REAL(iter)
#ifdef PARTICLES
  IF(CountNbrOfLostParts.AND.(NbrOfLostParticlesTotal.GT.0))THEN
    WRITE(UNIT_stdOut,'(A,I0,A,I0,A1)')' Total number of lost particles : ',NbrOfLostParticlesTotal,' (difference to last output: +',&
        NbrOfLostParticlesTotal-NbrOfLostParticlesTotal_old,')'
    NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal
  END IF
#endif /*PARICLES*/
END IF !MPIroot
#if defined(IMPA) || defined(ROS)
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A32,I12)') ' Total iteration Linear Solver    ',totalIterLinearSolver
TotalIterLinearSolver=0
#endif /*IMPA && ROS*/
#if defined(IMPA) && defined(PARTICLES)
SWRITE(UNIT_stdOut,'(A32,I12)')  ' IMPLICIT PARTICLE TREATMENT    '
SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration Newton         ',nPartNewton
SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration GMRES          ',TotalPartIterLinearSolver
IF(nPartNewton.GT.0)THEN
  SWRITE(UNIT_stdOut,'(A35,F12.2)')' Average GMRES steps per Newton    ',REAL(TotalPartIterLinearSolver)&
                                                                        /REAL(nPartNewton)
END IF
nPartNewTon=0
TotalPartIterLinearSolver=0
#if IMPA
SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration outer-Newton    ',TotalFullNewtonIter
totalFullNewtonIter=0
#endif
SWRITE(UNIT_stdOut,'(132("="))')
#endif /*IMPA && PARTICLE*/
#if defined(ROS) && defined(PARTICLES)
SWRITE(UNIT_stdOut,'(A32,I12)')  ' IMPLICIT PARTICLE TREATMENT    '
SWRITE(UNIT_stdOut,'(A32,I12)')  ' Total iteration GMRES          ',TotalPartIterLinearSolver
TotalPartIterLinearSolver=0
SWRITE(UNIT_stdOut,'(132("="))')
#endif /*IMPA && PARTICLE*/
END SUBROUTINE WriteInfoStdOut


END MODULE MOD_TimeDisc
