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
USE MOD_Globals_Vars           ,ONLY: SimulationEfficiency,PID,WallTime,ProjectName
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: time,TEnd,dt,iter,IterDisplayStep,DoDisplayIter,dt_Min,tAnalyze
USE MOD_TimeDisc_Vars          ,ONLY: time_start
#if USE_LOADBALANCE
USE MOD_TimeDisc_Vars          ,ONLY: dtWeight
#if defined(PARTICLES)
USE MOD_Particle_Vars          ,ONLY: WriteMacroVolumeValues, WriteMacroSurfaceValues, MacroValSampTime
#endif /*defined(PARTICLES)*/
#endif /*USE_LOADBALANCE*/
USE MOD_TimeAverage_vars       ,ONLY: doCalcTimeAverage
USE MOD_TimeAverage            ,ONLY: CalcTimeAverage
USE MOD_Analyze                ,ONLY: PerformAnalyze
USE MOD_Analyze_Vars           ,ONLY: Analyze_dt,iAnalyze,nSkipAnalyze,SkipAnalyzeWindow,SkipAnalyzeSwitchTime,nSkipAnalyzeSwitch
USE MOD_Restart_Vars           ,ONLY: RestartTime,RestartWallTime
USE MOD_HDF5_Output_State      ,ONLY: WriteStateToHDF5
USE MOD_Mesh_Vars              ,ONLY: MeshFile,nGlobalElems
USE MOD_RecordPoints_Vars      ,ONLY: RP_onProc
USE MOD_RecordPoints           ,ONLY: WriteRPToHDF5!,RecordPoints
USE MOD_Restart_Vars           ,ONLY: DoRestart,FlushInitialState
#if !(USE_HDG)
USE MOD_PML_Vars               ,ONLY: DoPML,PMLTimeRamp
USE MOD_PML                    ,ONLY: PMLTimeRamping
#if USE_LOADBALANCE
#ifdef maxwell
#if defined(ROS) || defined(IMPA)
USE MOD_Precond_Vars           ,ONLY:UpdatePrecondLB
#endif /*ROS or IMPA*/
#endif /*maxwell*/
#endif /*USE_LOADBALANCE*/
#else
USE MOD_HDG_Vars               ,ONLY: iterationTotal,RunTimeTotal
#endif /*USE_HDG*/
#ifdef PP_POIS
USE MOD_Equation               ,ONLY: EvalGradient
#endif /*PP_POIS*/
#if USE_MPI
#if USE_LOADBALANCE
USE MOD_LoadBalance            ,ONLY: LoadBalance,ComputeElemLoad
USE MOD_LoadBalance_Vars       ,ONLY: DoLoadBalance,ElemTime,DoLoadBalanceBackup,LoadBalanceSampleBackup,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: LoadBalanceSample,PerformLBSample,PerformLoadBalance,LoadBalanceMaxSteps,nLoadBalanceSteps
USE MOD_Restart_Vars           ,ONLY: DoInitialAutoRestart
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimeField
USE MOD_Restart_Vars           ,ONLY: RestartFile
USE MOD_HDF5_output            ,ONLY: RemoveHDF5
#endif /*USE_LOADBALANCE*/
#else
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
#endif /*USE_MPI*/
#ifdef PARTICLES
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_HDF5_Output_Particles  ,ONLY: WriteElectroMagneticPICFieldToHDF5
USE MOD_HDF5_Output_State      ,ONLY: WriteIMDStateToHDF5
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcEMFieldOutput
USE MOD_HDF5_Output_Particles  ,ONLY: FillParticleData
#endif /*PARTICLES*/
#ifdef PARTICLES
USE MOD_RayTracing             ,ONLY: RayTracing
!USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_Particle_Vars          ,ONLY: DoImportIMDFile
#if USE_MPI
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#endif /*USE_MPI*/
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
USE MOD_Particle_Tracking_vars ,ONLY: tTracking,tLocalization,nTracks,MeasureTrackTime
#if (USE_MPI) && (USE_LOADBALANCE) && defined(PARTICLES)
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimePart
#endif /* USE_LOADBALANCE && PARTICLES*/
USE MOD_Particle_Sampling_Adapt,ONLY: AdaptiveBCSampling
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation       ,ONLY: InitAnalyticalParticleState
#endif /*CODE_ANALYZ*/
#endif /*PARTICLES*/
USE MOD_Output                 ,ONLY: PrintStatusLine
USE MOD_TimeStep
USE MOD_TimeDiscInit           ,ONLY: InitTimeStep,UpdateTimeStep
#if defined(PARTICLES) && USE_HDG
USE MOD_Part_BR_Elecron_Fluid  ,ONLY: SwitchBRElectronModel,UpdateVariableRefElectronTemp
USE MOD_HDG_Vars               ,ONLY: BRConvertMode,BRTimeStepBackup,BRTimeStepMultiplier,UseBRElectronFluid
USE MOD_HDG_Vars               ,ONLY: CalcBRVariableElectronTemp
#endif /*defined(PARTICLES) && USE_HDG*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars               ,ONLY: MPIW8TimeSim
#endif /*defined(MEASURE_MPI_WAIT)*/
#if defined(PARTICLES)
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPointsPerDebyeLength,CalcPICTimeStep
USE MOD_Part_Tools             ,ONLY: ReduceMaxParticleNumber
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: tStart                   !> simulation time at the beginning of the simulation
REAL            :: tPreviousAnalyze         !> time of previous analyze.
                                            !> Used for Nextfile info written into previous file if greater tAnalyze
REAL            :: tPreviousAverageAnalyze  !> time of previous Average analyze.
REAL            :: tZero
INTEGER(KIND=8) :: iter_PID                 !> iteration counter since last InitPiclas call for PID calculation
REAL            :: WallTimeStart            !> wall time of simulation start
REAL            :: WallTimeEnd              !> wall time of simulation end
LOGICAL         :: finalIter
#if USE_LOADBALANCE
REAL            :: RestartTimeBackup
LOGICAL         :: ForceInitialLoadBalance  !> Set true when initial load balance steps are completed and force the load balance
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
ForceInitialLoadBalance = .FALSE. ! Initialize
#endif /*USE_LOADBALANCE*/
tPreviousAnalyze=RestartTime
! first average analyze is not written at start but at first tAnalyze
tPreviousAverageAnalyze=tAnalyze
! saving the start of the simulation as restart time is overwritten during load balance step
! In case the overwritten one is used, state write out is performed only next Nth Analyze_dt after restart instead after Analyze_dt
!   w/o  tZero: nSkipAnalyze=5 , restart after iAnalyze=2 , next write out without any restarts after iAnalyze=7
!   with tZero: nSkipAnalyze=5 , restart after iAnalyze=2 , next write out without any restarts after iAnalyze=5
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
CALL InitTimeStep() ! Initial time step calculation for dt_Min
iter     = 0
iter_PID = 0

! fill recordpoints buffer (first iteration)
!IF(RP_onProc) CALL RecordPoints(iter,t,forceSampling=.TRUE.)

! Ray tracing
#if defined(PARTICLES)
IF(.NOT.DoRestart) CALL RayTracing()
#endif /*defined(PARTICLES)*/

CALL PrintStatusLine(time,dt,tStart,tEnd,1)

#if defined(PARTICLES) && defined(CODE_ANALYZE)
! Set specific particle position and velocity (calculated from an analytical expression)
CALL InitAnalyticalParticleState() ! Requires dt
#endif /*defined(PARTICLES) && defined(CODE_ANALYZE)*/

#if defined(PARTICLES)
IF(CalcPointsPerDebyeLength.OR.CalcPICTimeStep)THEN
  CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
END IF ! CalcPointsPerDebyeLength.OR.CalcPICTimeStep
#endif
CALL PerformAnalyze(time,FirstOrLastIter=.TRUE.,OutPutHDF5=.FALSE.)

#ifdef PARTICLES
IF(DoImportIMDFile)THEN
  CALL WriteIMDStateToHDF5() ! Write IMD particles to state file (and TTM if it exists)
  IF(.NOT.DoRestart) RETURN
END IF ! DoImportIMDFile
#endif /*PARTICLES*/
IF((.NOT.DoRestart).OR.FlushInitialState.OR.(.NOT.FILEEXISTS(TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',time))//'.h5'))) THEN
#if defined(PARTICLES)
  CALL FillParticleData() ! Fill the SFC-ordered particle arrays
#endif /*defined(PARTICLES)*/
  CALL WriteStateToHDF5(TRIM(MeshFile),time,tPreviousAnalyze) ! Write initial state to file
END IF

! if measurement of particle tracking time (used for analyze, load balancing uses own time measurement for tracking)
#ifdef PARTICLES
IF(MeasureTrackTime)THEN
  nTracks=0
  tTracking=0
  tLocalization=0
END IF
IF(CalcEMFieldOutput) CALL WriteElectroMagneticPICFieldToHDF5() ! Write magnetic field to file
#endif /*PARTICLES*/

! No computation needed if tEnd=tStart!
IF(ALMOSTEQUALRELATIVE(time,tEnd,1e-10))RETURN

!-----------------------------------------------------------------------------------------------------------------------------------
! iterations starting up from here
!-----------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_StdOut,*)'CALCULATION RUNNING...'
WallTimeStart=PICLASTIME()
CALL CPU_TIME(time_start)

DO !iter_t=0,MaxIter

#if defined(PARTICLES) && USE_HDG
  ! Check if BR<->kin switch is active
  IF(BRConvertMode.NE.0) CALL SwitchBRElectronModel()
  ! Restore the initial dt_Min from InitTimeStep() as it might have been changed in the previous time step due to BR<->kin switch
  dt_Min(DT_MIN) = BRTimeStepBackup
  ! Adjust the time step when BR electron fluid is active. Usually BR electron time step is XX times larger than the fully kinetic
  IF(UseBRElectronFluid) dt_Min(DT_MIN) = BRTimeStepMultiplier*dt_Min(DT_MIN)
#endif /*defined(PARTICLES) && USE_HDG*/

  CALL UpdateTimeStep()

  IF(doCalcTimeAverage) CALL CalcTimeAverage(.FALSE.,dt,time,tPreviousAverageAnalyze) ! tPreviousAnalyze not used if finalize_flag=false

#if !(USE_HDG)
  IF(DoPML) CALL PMLTimeRamping(time,PMLTimeRamp)
#else
  IF(MPIroot)THEN
    iterationTotal = 0
    RunTimeTotal   = 0.
  END IF ! MPIroot
#endif /*NOT USE_HDG*/

  CALL PrintStatusLine(time,dt,tStart,tEnd,1)

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
#elif (PP_TimeDiscMethod==507)
  CALL TimeStepPoissonByHigueraCary() ! Higuera-Cary, Poisson
#elif (PP_TimeDiscMethod==508)
  CALL TimeStepPoissonByBorisLeapfrog() ! Boris-Leapfrog, Poisson
#else
  CALL TimeStepPoissonByLSERK() ! Runge Kutta Explicit, Poisson
#endif
#else
  CALL abort(__STAMP__,'Timedisc 50x only available for EQNSYS Poisson! PP_N=',IntInfoOpt=PP_N)
#endif /*USE_HDG*/
#elif (PP_TimeDiscMethod==600)
  CALL TimeStep_Radiation()
#endif
  ! calling the analyze routines
  iter     = iter+1
  iter_PID = iter_PID+1
  time     = time+dt

  IF(MPIroot) THEN
    IF(DoDisplayIter.AND.(MOD(iter,IterDisplayStep).EQ.0)) WRITE(UNIT_stdOut,'(A,I21,A6,ES26.16E3,25X)')" iter:", iter,"time:",time
  END IF

  ! Calling the analyze routines in last iteration
  finalIter = ALMOSTEQUAL(dt,dt_Min(DT_END))

#if defined(PARTICLES) && USE_HDG
  ! Depending on kinetic/BR model, set the reference electron temperature for t^n+1, therefore "add" -dt to the calculation
  IF(CalcBRVariableElectronTemp) CALL UpdateVariableRefElectronTemp(-dt)
#endif /*defined(PARTICLES) && USE_HDG*/
  CALL PerformAnalyze(time,FirstOrLastIter=finalIter,OutPutHDF5=.FALSE.) ! analyze routines are not called here in last iter
#ifdef PARTICLES
  ! sampling of near adaptive boundary element values
  IF(UseAdaptiveBC.OR.(nPorousBC.GT.0)) CALL AdaptiveBCSampling()
#endif /*PARICLES*/

  ! Analysis (possible PerformAnalyze+WriteStateToHDF5 and/or LoadBalance)
  !IF ((dt.EQ.dt_Min(DT_ANALYZE)).OR.(dt.EQ.dt_Min(DT_END))) THEN   ! timestep is equal to time to analyze or end
#if USE_LOADBALANCE
  ! For automatic initial restart, check if the number of sampling steps has been achieved and force a load balance step, but skip
  ! this procedure in the final iteration after which the simulation is finished
  !      DoInitialAutoRestart: user-activated load balance restart in first time step (could also be during a normal restart)
  ! iter.GE.LoadBalanceSample: as soon as the number of time steps for sampling is reached, perform the load balance restart
  !                 finalIter: prevent removal of last state file even though no load balance restart was performed
  IF(DoInitialAutoRestart.AND.(iter.GE.LoadBalanceSample).AND.(.NOT.finalIter)) ForceInitialLoadBalance=.TRUE.

  IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.finalIter.OR.ForceInitialLoadBalance)THEN
#else
  IF(ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)).OR.finalIter)THEN
#endif /*USE_LOADBALANCE*/
    WallTimeEnd=PICLASTIME()
    IF(MPIroot)THEN ! determine the SimulationEfficiency and PID here,
                    ! because it is used in ComputeElemLoad -> WriteElemTimeStatistics
      WallTime             = WallTimeEnd-StartTime
      SimulationEfficiency = (time-RestartTime)/((WallTimeEnd-RestartWallTime)*nProcessors/3600.) ! in [s] / [CPUh]
      PID                  = (WallTimeEnd-WallTimeStart)*nProcessors/(nGlobalElems*(PP_N+1)**3*iter_PID)
    END IF
#if defined(MEASURE_MPI_WAIT)
    MPIW8TimeSim = MPIW8TimeSim + (WallTimeEnd-WallTimeStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if defined(PARTICLES) && !defined(LSERK) && !defined(IMPA) && !defined(ROS)
    CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
#endif
#if USE_MPI

#if USE_LOADBALANCE
#ifdef PARTICLES
    ! Check if loadbalancing is enabled with partweight and set PerformLBSample true to calculate elemtimes with partweight
    ! LoadBalanceSample is 0 if partweightLB or IAR_partweighlb are enabled. If only one of them is set Loadbalancesample switches
    ! during time loop
    IF (LoadBalanceSample.EQ.0 .AND. DoLoadBalance) PerformLBSample=.TRUE.
#endif /*PARICLES*/
    ! Routine calculates imbalance and if greater than threshold sets PerformLoadBalance=.TRUE.
    CALL ComputeElemLoad()
    ! Force load balance step after elem time has been calculated when doing an initial load balance step at iter=0
    IF(ForceInitialLoadBalance) PerformLoadBalance=.TRUE.
    ! Do not perform a load balance restart when the last timestep is performed
    IF(finalIter) PerformLoadBalance=.FALSE.
#if defined(maxwell) && (defined(ROS) || defined(IMPA))
    UpdatePrecondLB=PerformLoadBalance
#endif /*MAXWELL AND (ROS or IMPA)*/
#endif /*USE_LOADBALANCE*/

#else /*NOT USE_MPI*/
    CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time_opt=time)
#endif /*USE_MPI*/

    ! Adjust nSkipAnalyze when, e.g., also adjusting the time step (during load balance, this value is over-written)
    IF(MOD(time,SkipAnalyzeWindow).GT.SkipAnalyzeSwitchTime) nSkipAnalyze = nSkipAnalyzeSwitch

    !--- Perform analysis and write state file .h5
    ! MOD(iAnalyze,nSkipAnalyze).EQ.0: Use nSkipAnalyze to skip analyze steps
    ! finalIter=T: last iteration of the simulation is reached, hence, always perform analysis and output to hdf5
#if USE_LOADBALANCE
    ! PerformLoadBalance.AND.UseH5IOLoadBalance: Load balance step will be performed and load balance restart via hdf5 IO active
    ! .NOT.(DoInitialAutoRestart.AND.(.NOT.UseH5IOLoadBalance)): Skip I/O for initial LB because MOD might give 0
    IF( ((.NOT.(DoInitialAutoRestart.AND.(.NOT.UseH5IOLoadBalance))).AND.MOD(iAnalyze,nSkipAnalyze).EQ.0)&
        .OR. (PerformLoadBalance.AND.UseH5IOLoadBalance) .OR. finalIter )THEN
#else
    IF(MOD(iAnalyze,nSkipAnalyze).EQ.0 .OR. finalIter)THEN
#endif /*USE_LOADBALANCE*/
      ! Analyze for output
#if defined(PARTICLES)
      IF(CalcPointsPerDebyeLength.OR.CalcPICTimeStep)THEN
        CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
      END IF ! CalcPointsPerDebyeLength.OR.CalcPICTimeStep
#endif
      CALL PerformAnalyze(time, FirstOrLastIter=finalIter, OutPutHDF5=.TRUE.) ! analyze routines are called here in last iter
      ! write information out to std-out of console
      CALL PrintStatusLine(time,dt,tStart,tEnd,2)
      CALL WriteInfoStdOut()
#if defined(PARTICLES)
      CALL FillParticleData() ! Fill the SFC-ordered particle arrays for LB or I/O
#endif /*defined(PARTICLES)*/
      ! Write state to file
      CALL WriteStateToHDF5(TRIM(MeshFile),time,tPreviousAnalyze)
      IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,time,tPreviousAverageAnalyze)
      ! Write recordpoints data to hdf5
      IF(RP_onProc) CALL WriteRPtoHDF5(tAnalyze,.TRUE.)
      tPreviousAnalyze        = tAnalyze
      tPreviousAverageAnalyze = tAnalyze
      SWRITE(UNIT_StdOut,'(132("-"))')
#if USE_LOADBALANCE && defined(PARTICLES)
    ELSEIF(PerformLoadBalance) THEN
      CALL FillParticleData() ! Fill the SFC-ordered particle arrays for LB
#endif /*USE_LOADBALANCE && defined(PARTICLES)*/
    END IF ! actual analyze is done

    iter_PID=0

    !--- Check if load balancing must be performed
#if USE_LOADBALANCE
    IF((DoLoadBalance.AND.PerformLBSample.AND.(LoadBalanceMaxSteps.GT.nLoadBalanceSteps)).OR.ForceInitialLoadBalance)THEN
      IF(PerformLoadBalance) THEN
        ! DO NOT DELETE THIS: ONLY recalculate the timestep when the mesh is changed!
        !CALL InitTimeStep() ! re-calculate time step after load balance is performed
        RestartTimeBackup = RestartTime! make backup of original restart time
        RestartTime       = time       ! Set restart simulation time to current simulation time because the time is not read from
                                       ! the state file
        RestartWallTime = PICLASTIME() ! Set restart wall time if a load balance step is performed
        dtWeight        = 1.           ! is initialized in InitTimeDisc which is not called in LoadBalance, but needed for restart
                                       ! (RestartHDG)
      END IF
      CALL LoadBalance()
    ELSE
      ElemTime      = 0. ! nullify ElemTime before measuring the time in the next cycle
#ifdef PARTICLES
      ElemTimePart  = 0.
#endif /*PARTICLES*/
      ElemTimeField = 0.
    END IF
    PerformLBSample=.FALSE. ! Deactivate load balance sampling

    ! Switch off Initial Auto Restart (initial load balance) after the restart was performed
    IF (DoInitialAutoRestart) THEN
      ! Remove the extra state file written for load balance (only when load balance restart via hdf5 IO was performed)
      IF(PerformLoadBalance.AND.UseH5IOLoadBalance) CALL RemoveHDF5(RestartFile)
      ! Get original settings from backup variables
      DoInitialAutoRestart = .FALSE.
      ForceInitialLoadBalance = .FALSE.
      DoLoadBalance        = DoLoadBalanceBackup
      LoadBalanceSample    = LoadBalanceSampleBackup
      ! Set to iAnalyze zero so that this first analysis is not counted and the next analysis is the first one,
      ! but only if the initial load balance restart and dt_Analyze did not coincide
      IF(.NOT.ALMOSTEQUALRELATIVE(dt, dt_Min(DT_ANALYZE), 1E-5)) iAnalyze=0
      ! Set time of the state file that was created before automatic initial restart (to be written in the next state file)
      tPreviousAnalyze = RestartTimeBackup
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
  CALL ReduceMaxParticleNumber()
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
USE MOD_Particle_Tracking_vars ,ONLY: NbrOfNewLostParticlesTotal
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
#if USE_MPI && defined(PARTICLES)
INTEGER             :: NbrOfLostParticlesTotal_old_tmp
#endif /*USE_MPI && defined(PARTICLES)*/
!===================================================================================================================================
#ifdef PARTICLES
IF(CountNbrOfLostParts)THEN
#if USE_MPI
  NbrOfLostParticlesTotal_old_tmp = NbrOfLostParticlesTotal ! keep old value
  ! Allreduce is required because of the particle output to .h5 in which all processors must take place
  CALL MPI_ALLREDUCE(NbrOfLostParticles , NbrOfLostParticlesTotal , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_PICLAS , IERROR)
  NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + NbrOfLostParticlesTotal_old_tmp ! add old value
#else
  NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + NbrOfLostParticles
#endif /*USE_MPI*/
  NbrOfNewLostParticlesTotal  = NbrOfLostParticlesTotal-NbrOfLostParticlesTotal_old
  NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal
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
  WRITE(UNIT_StdOut,'(A,ES16.7)')' Timestep  : ',dt_Min(DT_MIN)
  WRITE(UNIT_stdOut,'(A,ES16.7)')'#Timesteps : ',REAL(iter)
#ifdef PARTICLES
  IF(CountNbrOfLostParts.AND.(NbrOfLostParticlesTotal.GT.0))THEN
    WRITE(UNIT_stdOut,'(A,I0,A,I0,A1)')' Total number of lost particles : ',NbrOfLostParticlesTotal,' (difference to last output: +',&
        NbrOfNewLostParticlesTotal,')'
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
