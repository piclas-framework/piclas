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


!===================================================================================================================================
!> Module for the Temporal discretization
!===================================================================================================================================
MODULE MOD_TimeDiscInit
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DefineParametersTimeDisc,InitTime,InitTimeDisc,InitTimeStep,UpdateTimeStep,FinalizeTimeDisc
!===================================================================================================================================
CONTAINS

!=================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTimeDisc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TimeDisc")
CALL prms%CreateRealOption('ManualTimeStep'  , 'Manual timestep [sec].'                                                , '-1.0')
CALL prms%CreateRealOption('TEnd'            , "End time of the simulation (mandatory).")
CALL prms%CreateRealOption('CFLScale'        , "Scaling factor for the theoretical CFL number; typical range 0.1..1.0" , '1.0')
CALL prms%CreateIntOption( 'maxIter'         , "Stop simulation when specified number of timesteps has been performed.", '-1')
CALL prms%CreateIntOption( 'NCalcTimeStepMax', "Compute dt at least after every Nth timestep."                         , '1')
CALL prms%CreateIntOption( 'IterDisplayStep' , "Step size of iteration that are displayed."                            , '1')
END SUBROUTINE DefineParametersTimeDisc


!===================================================================================================================================
!> The temporal genesis
!===================================================================================================================================
SUBROUTINE InitTime()
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars          ,ONLY: Time,TEnd,tAnalyze,dt_Min
USE MOD_Restart_Vars           ,ONLY: RestartTime
#ifdef PARTICLES
USE MOD_DSMC_Vars              ,ONLY: Iter_macvalout,Iter_macsurfvalout
USE MOD_Particle_Vars          ,ONLY: WriteMacroVolumeValues,WriteMacroSurfaceValues,MacroValSampTime
#endif /*PARTICLES*/
USE MOD_Analyze_Vars           ,ONLY: Analyze_dt,iAnalyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Setting the time variable, RestartTime is either zero or the actual restart time
Time=RestartTime
#ifdef PARTICLES
iter_macvalout=0
iter_macsurfvalout=0
IF (WriteMacroVolumeValues.OR.WriteMacroSurfaceValues) MacroValSampTime = Time
#endif /*PARTICLES*/
iAnalyze=1
! Determine the first analyze time
tAnalyze=MIN(RestartTime+Analyze_dt,tEnd)

! fill initial analyze stuff
dt_Min(DT_ANALYZE) = tAnalyze-Time ! Time to next analysis, put in extra variable so number does not change due to numerical errors
dt_Min(DT_END)     = tEnd    -Time ! dito for end Time
#if defined(PARTICLES) && USE_HDG
dt_Min(DT_BR_SWITCH) = HUGE(1.)
#endif/*defined(PARTICLES) && USE_HDG*/

END SUBROUTINE InitTime


SUBROUTINE InitTimeDisc()
!===================================================================================================================================
! Get information for end time and max time steps from ini file
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools   ,ONLY: GetReal,GetInt, GETLOGICAL
USE MOD_TimeDisc_Vars ,ONLY: IterDisplayStepUser
USE MOD_TimeDisc_Vars ,ONLY: CFLScale,dt,TimeDiscInitIsDone,RKdtFrac,RKdtFracTotal,dtWeight
USE MOD_TimeDisc_Vars ,ONLY: IterDisplayStep,DoDisplayIter
USE MOD_TimeDisc_Vars ,ONLY: ManualTimeStep,useManualTimestep
USE MOD_TimeDisc_Vars ,ONLY: TEnd
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
!USE MOD_TimeDisc_Vars          ,ONLY: U2t_temp
USE MOD_PML_Vars               ,ONLY: nPMLElems
USE MOD_PML_Vars               ,ONLY: PMLnVar
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
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
CALL FillCFL_DFL()

!--- Read Manual Time Step
useManualTimeStep = .FALSE.
ManualTimeStep = GETREAL('ManualTimeStep')
IF (ManualTimeStep.GT.0.0) useManualTimeStep=.True.

! Read the maximum number of time steps MaxIter and the end time TEnd from ini file
TEnd=GetReal('TEnd') ! must be read in here due to DSMC_init

! read in requested IterDisplayStep (i.e. how often the message "iter: etc" is displayed, might be changed dependent on Particle-dt)
DoDisplayIter=.FALSE.
IterDisplayStepUser = GETINT('IterDisplayStep','1')
IterDisplayStep = IterDisplayStepUser
IF(IterDisplayStep.GE.1) DoDisplayIter=.TRUE.

! Set timestep to a large number
dt=HUGE(1.)
#if (PP_TimeDiscMethod==1)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK3-3 '
#elif (PP_TimeDiscMethod==2)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-5 '
#elif (PP_TimeDiscMethod==4)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Direct Simulation Monte Carlo (DSMC)'
#elif (PP_TimeDiscMethod==6)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-14 '
#elif (PP_TimeDiscMethod==300)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Fokker-Planck (FP) Collision Operator'
#elif (PP_TimeDiscMethod==400)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Bhatnagar-Gross-Krook (BGK) Collision Operator'
#elif (PP_TimeDiscMethod==500)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Euler, Poisson'
#elif (PP_TimeDiscMethod==501)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK3-3, Poisson'
#elif (PP_TimeDiscMethod==502)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-5, Poisson'
#elif (PP_TimeDiscMethod==506)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: LSERK4-14, Poisson'
#elif (PP_TimeDiscMethod==507)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Higuera-Cary, Poisson'
#elif (PP_TimeDiscMethod==508)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Boris-Leapfrog, Poisson'
#elif (PP_TimeDiscMethod==509)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Leapfrog, Poisson'
#elif (PP_TimeDiscMethod==600)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Radiation'
#elif (PP_TimeDiscMethod==700)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: ED-DVM or DUGKS'
#elif (PP_TimeDiscMethod==701)
  SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: Explicit finite volumes'
#endif

RKdtFrac      = 1.
RKdtFracTotal = 1.
dtWeight      = 1.

TimediscInitIsDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitTimeDisc


!===================================================================================================================================
!> Initial time step calculation for new timedisc-loop
!===================================================================================================================================
SUBROUTINE InitTimeStep()
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars         ,ONLY: dt, dt_Min,ManualTimeStep,useManualTimestep,sdtCFLOne
#if ! (USE_HDG)
USE MOD_CalcTimeStep          ,ONLY: CalcTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: CFLtoOne
#endif
USE MOD_ReadInTools           ,ONLY: GETREAL
#if defined(PARTICLES) && USE_HDG
USE MOD_HDG_Vars              ,ONLY: BRTimeStepMultiplier,UseBRElectronFluid,BRTimeStepBackup,BRConvertMode
USE MOD_Part_BR_Elecron_Fluid ,ONLY: GetNextBRSwitchTime
#endif /*defined(PARTICLES) && USE_HDG*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: DoLoadBalanceBackup,LoadBalanceSampleBackup,DoLoadBalance
USE MOD_LoadBalance_Vars      ,ONLY: LoadBalanceSample,PerformLBSample
USE MOD_Restart_Vars          ,ONLY: DoInitialAutoRestart,InitialAutoRestartSample
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! For tEnd != tStart we have to advance the solution in time
IF(useManualTimeStep)THEN
  ! particle time step is given externally and not calculated through the solver
  dt_Min(DT_MIN)=ManualTimeStep
ELSE ! .NO. ManualTimeStep
  ! time step is calculated by the solver
  ! first Maxwell time step for explicit LSRK
#if !(USE_HDG)
  dt_Min(DT_MIN) = CalcTimeStep()
  sdtCFLOne      = 1.0/(dt_Min(DT_MIN)*CFLtoOne)
#else
  dt_Min(DT_MIN) = 0
  sdtCFLOne      = -1.0 !dummy for HDG!!!
#endif /*USE_HDG*/

END IF ! useManualTimestep

#if defined(PARTICLES) && USE_HDG
! Check if BR<->kin switch is active
IF(BRConvertMode.NE.0) CALL GetNextBRSwitchTime()

! Adjust the time step when BR electron fluid is active
BRTimeStepBackup = dt_Min(DT_MIN)
IF(UseBRElectronFluid) dt_Min(DT_MIN) = BRTimeStepMultiplier*dt_Min(DT_MIN)
#endif /*defined(PARTICLES) && USE_HDG*/

! Select the smallest time delta
dt = MINVAL(dt_Min)
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A,ES16.7)')'Initial Timestep  : ', dt

! Sanity check: dt must be greater zero, this can happen when using a fixed timestep TD and not supplying a value for dt
IF(dt.LE.0.)THEN
  IPWRITE(UNIT_StdOut,*) "Initial timestep dt_Min=[",dt_Min,"]"
  !CALL abort(__STAMP__,'Time step is less/equal zero: dt = ',RealInfoOpt=dt)
END IF

#if USE_LOADBALANCE
IF (DoInitialAutoRestart) THEN

  ! Set general load balance flag ON
  DoLoadBalanceBackup   = DoLoadBalance ! Backup
  DoLoadBalance         = .TRUE.        ! Force TRUE (during automatic restart, original variable might be false)

  ! Backup number of samples required for each load balance
  LoadBalanceSampleBackup = LoadBalanceSample        ! Backup: this is zero when PerformPartWeightLB=.TRUE.
  LoadBalanceSample       = InitialAutoRestartSample ! this is zero when InitialAutoRestartPartWeight=.TRUE.

  ! Activate sampling in first time step
  PerformLBSample=.TRUE.

  ! Sanity check: initial automatic restart must happen before tAnalyze is reached (tAnalyze < LoadBalanceSample*dt not implemented)
  DO WHILE(LoadBalanceSample*dt.GT.dt_Min(DT_ANALYZE).AND.(LoadBalanceSample.GT.1))
   LoadBalanceSample = LoadBalanceSample-1
  END DO

END IF
#endif /*USE_LOADBALANCE*/
END SUBROUTINE InitTimeStep


!===================================================================================================================================
!> Update time step at the beginning of each timedisc loop
!===================================================================================================================================
SUBROUTINE UpdateTimeStep()
! MODULES
USE MOD_Globals          ,ONLY: abort,UNIT_StdOut,LESSEQUALTOLERANCE
#if USE_MPI
USE MOD_Globals          ,ONLY: myrank
#endif /*USE_MPI*/
USE MOD_TimeDisc_Vars    ,ONLY: dt,time,tEnd,tAnalyze,dt_Min,dtWeight
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: DoLoadBalance,LoadBalanceSample,PerformLBSample
#endif /*USE_LOADBALANCE*/
#if (PP_TimeDiscMethod==509)
USE MOD_TimeDisc_Vars    ,ONLY: iter,dt_old
#if USE_MPI
USE MOD_Globals          ,ONLY: MPIRoot
#endif /*USE_MPI*/
#endif /*(PP_TimeDiscMethod==509)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
dt_Min(DT_ANALYZE) = tAnalyze-time ! Time to next analysis, put in extra variable so number does not change due to numerical errors
dt_Min(DT_END)     = tEnd    -time ! Do the same for end time

#if (PP_TimeDiscMethod==509)
IF (iter.GT.0) dt_old=dt
#endif /*(PP_TimeDiscMethod==509)*/

dt       = MINVAL(dt_Min)
dtWeight = dt/dt_Min(DT_MIN) ! Might be further decreased by RK-stages

#if (PP_TimeDiscMethod==509)
IF (iter.EQ.0) THEN
  dt_old=dt
ELSE IF (ABS(dt-dt_old).GT.1.0E-6*dt_old) THEN
  SWRITE(UNIT_StdOut,'(A,G0)')'WARNING: dt changed from last iter by a relative difference of ',(dt-dt_old)/dt_old
END IF
#endif /*(PP_TimeDiscMethod==509)*/

#if USE_LOADBALANCE
! Activate normal load balancing (NOT initial restart load balancing)
! 1.) Catch all iterations within sampling interval (make sure to get the first iteration in interval): LESSEQUALTOLERANCE(a,b,tol)
! 2.)             Load balancing is activated: DoLoadBalance=T
IF( LESSEQUALTOLERANCE(dt_Min(DT_ANALYZE), LoadBalanceSample*dt, 1e-5) &
    .AND. DoLoadBalance) PerformLBSample=.TRUE. ! Activate load balancing in this time step
#endif /*USE_LOADBALANCE*/

IF(dt_Min(DT_ANALYZE)-dt.LT.dt/100.0) dt = dt_Min(DT_ANALYZE) ! Increase time step if the NEXT time step would be smaller than dt/100
IF(    dt_Min(DT_END)-dt.LT.dt/100.0) dt = dt_Min(DT_END)     ! Increase time step if the LAST time step would be smaller than dt/100

! Sanity check: dt must be greater zero
IF(dt.LE.0.)THEN
  IPWRITE(UNIT_StdOut,*) "time=[", time,"], dt_Min=[",dt_Min,"], tAnalyze=[",tAnalyze,"]"
  IF(dt.LT.0.)THEN
    CALL abort(__STAMP__,'dt < 0: Is something wrong with tEnd? Error in dt_Min(DT_END) or dt_Min(DT_ANALYZE)! dt=',RealInfoOpt=dt)
  ELSE
    CALL abort(__STAMP__,'Time step is less/equal zero: dt=',RealInfoOpt=dt)
  END IF
END IF

END SUBROUTINE UpdateTimeStep


!===================================================================================================================================
!> scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
!> Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" . For N=1-10,
!> input CFLscale can now be 1. and will be scaled adequately depending on PP_N, NodeType and TimeDisc method
!> DFLscale dependance not implemented jet.
!===================================================================================================================================
SUBROUTINE FillCFL_DFL()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY:CFLScale,CFLtoOne
#if (PP_TimeDiscMethod==2)
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

#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
IF(PP_N.GT.15) CALL abort(__STAMP__,'Polynomial degree is to high!',PP_N,999.)
CFLScale=CFLScale*CFLScaleAlpha(PP_N)
#endif
!scale with 2N+1
CFLScale = CFLScale/(2.*PP_N+1.)
SWRITE(UNIT_stdOut,'(A,ES16.7)') '   CFL:',CFLScale
END SUBROUTINE FillCFL_DFL


!===================================================================================================================================
!> Deallocate timedisc variables
!===================================================================================================================================
SUBROUTINE FinalizeTimeDisc()
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars, ONLY:TimeDiscInitIsDone
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
USE MOD_TimeDisc_Vars          ,ONLY: Ut_N
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
SDEALLOCATE(Ut_N)
#endif

TimeDiscInitIsDone = .FALSE.

END SUBROUTINE FinalizeTimeDisc

END MODULE MOD_TimeDiscInit