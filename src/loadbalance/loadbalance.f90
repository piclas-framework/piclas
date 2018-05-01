#include "boltzplatz.h"

!===================================================================================================================================
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#ifdef MPI
INTERFACE InitLoadBalance
  MODULE PROCEDURE InitLoadBalance
END INTERFACE

INTERFACE FinalizeLoadBalance
  MODULE PROCEDURE FinalizeLoadBalance
END INTERFACE

INTERFACE ComputeElemLoad
  MODULE PROCEDURE ComputeElemLoad
END INTERFACE

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

!INTERFACE LoadMeasure 
  !MODULE PROCEDURE LoadMeasure
!END INTERFACE

PUBLIC::InitLoadBalance,FinalizeLoadBalance,LoadBalance
!PUBLIC::LoadMeasure
PUBLIC::ComputeElemLoad
#endif /*MPI*/

PUBLIC::DefineParametersLoadBalance
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersLoadBalance()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("LoadBalance")

CALL prms%CreateLogicalOption( 'DoLoadBalance'                 ,  "Set flag for doing LoadBalance.", '.FALSE.')
CALL prms%CreateRealOption(    'Load-DeviationThreshold'       ,  "Define threshold for dynamic load-balancing.\n"//&
  "Restart performed if (Maxweight-Targetweight)/Targetweight > defined value." , value='0.10')
CALL prms%CreateRealOption(    'Particles-MPIWeight'           ,  "Define weight of particles for elem loads.\n"//&
  "(only used if ElemTime does not exist or DoLoadBalance=F).", value='0.02')
CALL prms%CreateIntOption(     'WeightDistributionMethod'      ,  "Method for distributing the elem to procs.\n"//&
  "DEFAULT: 1 if Elemtime exits else -1\n"//&
  "-1: elements are equally distributed\n"//&
  " 0: distribute to procs using elemloads\n"//&
  " 1: distribute to procs using elemloads, last proc recieves least\n"//&
  " 2: NOT WORKING\n"//&
  " 3: TODO DEFINE\n"//&
  " 4: TODO DEFINE\n"//&
  " 5/6: iterative smoothing of loads towards last proc\n")
CALL prms%CreateIntOption(     'LoadBalanceSample'           ,  "Define number of iterations (before analyze_dt)"//&
  " that are used for calculation of elemtime information" , value='1')

END SUBROUTINE DefineParametersLoadBalance

#ifdef MPI
SUBROUTINE InitLoadBalance()
!===================================================================================================================================
! init load balancing, new initialization of variables for load balancing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars
USE MOD_ReadInTools,          ONLY:GETLOGICAL, GETREAL, GETINT
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE ...'

#ifdef PARTICLES
! Read particle MPI weight in order to determine the ElemTime when no time measurement is performed
! Must be read in init (only once) and before the first load balance is determined because if no ElemTimes are used they are
! calculated with ParticleMPIWeight 
ParticleMPIWeight = GETREAL('Particles-MPIWeight','0.02')
IF (ParticleMPIWeight.LE.0.0) THEN
  CALL abort(&
      __STAMP__&
      ,' ERROR: Particle weight cannot be negative!')
END IF
#endif /*PARTICLES*/

IF(nProcessors.EQ.1)THEN
  DoLoadBalance=.FALSE. ! deactivate loadbalance for single computations
  SWRITE(UNIT_StdOut,'(a3,a45,a3,L33,a3,a7,a3)')' | ',TRIM("No LoadBalance (nProcessors=1): DoLoadBalance")       ,' | ',&
      DoLoadBalance   ,' | ',TRIM("INFO"),' | '
  DeviationThreshold=HUGE(1.0)
ELSE
  DoLoadBalance = GETLOGICAL('DoLoadBalance','F')
  DeviationThreshold  = GETREAL('Load-DeviationThreshold','0.10')
END IF
nLoadBalance = 0
nLoadBalanceSteps = 0
LoadBalanceSample  = GETINT('LoadBalanceSample')
PerformLBSample = .FALSE.

ALLOCATE( tTotal(1:LB_NTIMES)   )
ALLOCATE( tCurrent(1:LB_NTIMES) )
! Allocation length (1:number of loadbalance times)
! look into boltzplatz.h for more info about time names

tTotal=0.
tCurrent=0.

InitLoadBalanceIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance



SUBROUTINE ComputeElemLoad(PerformLoadbalance,time)
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the element load
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,nLoadBalance,tTotal,tCurrent
#ifndef PP_HDG
USE MOD_PML_Vars               ,ONLY: DoPML,nPMLElems,ElemToPML
#endif /*PP_HDG*/
USE MOD_LoadBalance_Vars       ,ONLY: DeviationThreshold
#ifdef PARTICLES
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nDeposPerElem,nTracksPerElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacefluxPerElem,nPartsPerBCElem
USE MOD_Particle_Tracking_vars ,ONLY: DoRefMapping,TriaTracking
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight
#endif /*PARTICLES*/
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance,PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
REAL ,INTENT(IN)             :: time
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: PerformLoadBalance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
REAL                  :: tDG, tPML
#ifdef PARTICLES
INTEGER(KIND=8)       :: HelpSum
REAL                  :: stotalDepos,stotalParts,sTotalTracks,stotalSurfacefluxes, sTotalBCParts
REAL                  :: tParts
#endif /*PARTICLES*/
!====================================================================================================================================

! If elem times are calculated by time measurement (PerformLBSample)
IF(PerformLBSample) THEN
  ! time measurement over LBSample before dt_analyze for all elements of one proc
  tTotal   = tTotal+tCurrent
  tCurrent = 0.

  ! number of load balance calls to Compute Elem Load
  nLoadBalance=nLoadBalance+1

  ! time used in dg routines
  tDG=(tTotal(LB_DG)+tTotal(LB_DGANALYZE))/REAL(PP_nElems)
  ! time used in pml routines
  tPML=0.
#ifndef PP_HDG
  IF(DoPML)THEN
    IF(nPMLElems.GT.0) tPML=tTotal(LB_PML)/REAL(nPMLElems)
  END IF
#endif /*PP_HDG*/

#ifdef PARTICLES
  stotalDepos=1.0
  sTotalTracks=1.0
  stotalSurfacefluxes=1.0
  ! calculate and weight particle number per element
  helpSum=SUM(nPartsPerElem)
  IF(helpSum.GT.0) THEN
    stotalParts=1.0/REAL(helpSum)
  ELSE
    stotalParts=1.0/REAL(PP_nElems)
    nPartsPerElem=1
  END IF
  ! sum of particle dependant load balance time 
  tParts = tTotal(LB_INTERPOLATION)+tTotal(LB_PUSH)+tTotal(LB_UNFP)+tTotal(LB_DSMC) !interpolation+push+unfp+dsmc(analyze)
  ! set and weight tracks per element
  IF (DoRefMapping .OR. TriaTracking) THEN
    helpSum=SUM(nTracksPerElem)
    IF(SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks=1.0/REAL(helpSum)
    ELSE
      stotalTracks=1.0/REAL(PP_nElems)
      nTracksPerElem=1
    END IF
  END IF
  ! set and weight depositions per element
  helpSum=SUM(nDeposPerElem)
  IF(helpSum.GT.0) THEN
    stotalDepos=1.0/REAL(helpSum)
  END IF
  ! calculate and weight number of surface fluxes (during setvelo,setinterenergies...)
  helpSum=SUM(nSurfacefluxPerElem)
  IF(helpSum.GT.0) THEN
    stotalSurfacefluxes=1.0/REAL(helpSum)
  END IF
  ! calculate weight number of particles used in sampling of bc element for adaptive particle bcs
  helpSum=SUM(nPartsPerBCElem)
  IF(helpSum.GT.0) THEN
    stotalBCParts=1.0/REAL(helpSum)
  ELSE
    stotalBCParts=1.0/REAL(PP_nElems)
    nPartsPerBCElem=1
  END IF
#endif /*PARTICLES*/

  ! distribute times of different routines on elements with respective weightings
  DO iElem=1,PP_nElems
    ElemTime(iElem) = ElemTime(iElem) + tDG
#ifndef PP_HDG
    IF(DoPML)THEN
      IF(ElemToPML(iElem).GT.0 ) ElemTime(iElem) = ElemTime(iElem) + tPML
    END IF
#endif /*PP_HDG*/

#ifdef PARTICLES
    ! add particle LB times to elements with respective weightings
    ElemTime(iElem) = ElemTime(iElem)                              &
        + tTotal(LB_ADAPTIVE) * nPartsPerBCElem(iElem)*sTotalBCParts &
        + tParts * nPartsPerElem(iElem)*sTotalParts    &
        + tTotal(LB_CARTMESHDEPO) * nPartsPerElem(iElem)*sTotalParts &
        + tTotal(LB_TRACK) * nTracksPerElem(iElem)*sTotalTracks &
        + tTotal(LB_SURFFLUX) * nSurfacefluxPerElem(iElem)*stotalSurfacefluxes
    ! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical'
    IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
      ElemTime(iElem) = ElemTime(iElem)                              &
          + tTotal(LB_DEPOSITION) * nDeposPerElem(iElem)*stotalDepos
    END IF
#endif /*PARTICLES*/
  END DO ! iElem=1,PP_nElems
#ifdef PARTICLES
ELSE ! no time measurement and particles are present: simply add the ParticleMPIWeight times the number of particles present
  DO iElem=1,PP_nElems
    ElemTime(iElem) = ElemTime(iElem) &
        + nPartsPerElem(iElem)*ParticleMPIWeight + 1.0
  END DO ! iElem=1,PP_nElems
  IF((MAXVAL(nPartsPerElem).GT.0).AND.(MAXVAL(ElemTime).LE.1.0))THEN
    IPWRITE (*,*) "parts, time =", MAXVAL(nPartsPerElem),MAXVAL(ElemTime)
    CALL abort(&
        __STAMP__&
        ,' ERROR: MAXVAL(nPartsPerElem).GT.0 but MAXVAL(ElemTime).LE.1.0 with ParticleMPIWeight=',RealInfoOpt=ParticleMPIWeight)
  END IF
#endif /*PARTICLES*/
END IF

! Determine sum of balance and calculate target balanced weight and communicate via MPI_ALLREDUCE
CALL ComputeImbalance()

! Fill .csv file for performance analysis and load balance: write data line
CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time=time)

! only check if imbalance is > a given threshold
PerformLoadBalance=.FALSE.
IF(CurrentImbalance.GT.DeviationThreshold) PerformLoadBalance=.TRUE.

#ifdef PARTICLES
nTracksPerElem=0
nDeposPerElem=0
nPartsPerElem=0
nSurfacefluxPerElem=0
nPartsPerBCElem=0
#endif /*PARTICLES*/

END SUBROUTINE ComputeElemLoad


SUBROUTINE LoadBalance(PerformLoadBalance)
!===================================================================================================================================
! routine perfoming the load balancing stuff
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Globals_vars     ,ONLY: InitializationWallTime
USE MOD_Preproc
USE MOD_Restart          ,ONLY: Restart
USE MOD_Boltzplatz_Init  ,ONLY: InitBoltzplatz,FinalizeBoltzplatz
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,nLoadBalanceSteps,NewImbalance,MinWeight,MaxWeight
#ifdef PARTICLES
USE MOD_PICDepo_Vars     ,ONLY: DepositionType
USE MOD_Particle_MPI     ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*PARTICLES*/
USE MOD_LoadBalance_Vars ,ONLY: CurrentImbalance, MaxWeight, MinWeight
USE MOD_LoadBalance_Vars ,ONLY: Currentimbalance
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)   :: PerformLoadBalance
!REAL,INTENT(IN)      :: Currentimbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: LB_Time,LB_StartTime
!===================================================================================================================================
! only do load-balance if necessary
IF(.NOT.PerformLoadBalance) THEN
  ElemTime=0.
  InitializationWallTime=0.
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'
! Measure init duration
LB_StartTime=BOLTZPLATZTIME()

nLoadBalanceSteps=nLoadBalanceSteps+1
! finialize all arrays
CALL FinalizeBoltzplatz(IsLoadBalance=.TRUE.)
! reallocate
CALL InitBoltzplatz(IsLoadBalance=.TRUE.) ! determines new imbalance in InitMesh() -> ReadMesh()

! restart
CALL Restart()

! compute imbalance 
!CALL ComputeImbalance(NewImbalance,MaxWeight,MinWeight,ElemTime)  ! --> new imbalance has already been calculated in mesh_readin

! zero ElemTime, the measurement starts again
ElemTime=0.

IF( NewImbalance.GT.CurrentImbalance ) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' LoadBalance successful!'
END IF
SWRITE(UNIT_stdOut,'(A25,E15.7)') ' OldImbalance: ', CurrentImbalance
SWRITE(UNIT_stdOut,'(A25,E15.7)') ' NewImbalance: ', NewImbalance
SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MaxWeight:    ', MaxWeight
SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MinWeight:    ', MinWeight

#ifdef PARTICLES
! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical'
IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
END IF
#endif /*PARTICLES*/

! Measure init duration
LB_Time=BOLTZPLATZTIME()
InitializationWallTime=LB_Time-LB_StartTime
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',InitializationWallTime,' sec ]'
SWRITE(UNIT_stdOut,'(A)')' LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE LoadBalance


SUBROUTINE ComputeImbalance()
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
! Maxweight:        maximum weight of all processes
! Minweight:        minimum weight of all processes
! targetweight:     current target weight
! CurrentImbalance: rel. deviation of current imbalance
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_LoadBalance_Vars,    ONLY:WeightSum, TargetWeight,CurrentImbalance, MaxWeight, MinWeight
USE MOD_LoadBalance_Vars,    ONLY:ElemTime, PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(.NOT. PerformLBSample) THEN
  WeightSum        = 0.
  TargetWeight     = 0.
  CurrentImbalance = -1.0
ELSE
  WeightSum=SUM(ElemTime)

  IF(ALMOSTZERO(WeightSum))THEN
    IPWRITE(*,*) 'Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)=.TRUE., WeightSum=',WeightSum
  END IF

  CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum,MaxWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum,MinWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)

  WeightSum    = TargetWeight ! Set total weight for writing to file
  TargetWeight = TargetWeight/nProcessors ! Calculate the average value that is supposed to be the optimally distributed weight

  !IF(ALMOSTZERO(WeightSum))CALL abort( __STAMP__&
      !,' ERROR: after ALLREDUCE, weight sum cannot be zero! WeightSum=',RealInfoOpt=WeightSum)

  ! new computation of current imbalance
  IF(ABS(TargetWeight).EQ.0.)THEN
    CurrentImbalance = 0.
  ELSE IF(ABS(TargetWeight).LT.0.0)THEN
    SWRITE(UNIT_stdOut,'(A,F8.2,A1)')&
        ' ERROR: after ALLREDUCE, WeightSum/TargetWeight cannot be zero! TargetWeight=[',TargetWeight,']'
    CurrentImbalance = HUGE(1.0)
  ELSE
    CurrentImbalance =  (MaxWeight-TargetWeight ) / TargetWeight
  END IF
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MaxWeight:        ', MaxWeight
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MinWeight:        ', MinWeight
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' TargetWeight:     ', TargetWeight
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' CurrentImbalance: ', CurrentImbalance

END IF

END SUBROUTINE ComputeImbalance


SUBROUTINE FinalizeLoadBalance()
!===================================================================================================================================
! Deallocate arrays
!===================================================================================================================================
! MODULES
USE MOD_LoadBalance_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( tTotal  )
SDEALLOCATE( tCurrent  )
InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
