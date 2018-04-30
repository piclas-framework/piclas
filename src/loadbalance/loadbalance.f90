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

CALL prms%CreateLogicalOption( 'DoLoadBalance'                 ,  "Flag for doing LoadBalance.", '.FALSE.')
CALL prms%CreateRealOption(    'Load-DeviationThreshold'       ,  "TODO-DEFINE-PARAMETER\n"//&
                                                                  "Threshold for load-balancing" , value='0.10')
CALL prms%CreateRealOption(    'Particles-MPIWeight'           ,  "Weight of particles for elem loads\n"//&
                                                                  "(only used if ElemTime does not exist or DoLoadBalance=F)."&
                                                               , value='0.02')
CALL prms%CreateIntOption(     'WeightDistributionMethod'      ,  "Method for distributing the elem loads. (def.: 1 or -1)")
CALL prms%CreateIntOption(     'LoadBalanceSample'           ,  "Number of iterations used for calculation of elemtime information"&
                                                               , value='10')

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
!DeviationThreshold  = 1.0+DeviationThreshold
nLoadBalance = 0
nLoadBalanceSteps = 0
LoadBalanceSample  = GETINT('LoadBalanceSample')
PerformLBSample = .FALSE.

ALLOCATE( tTotal(1:LB_NTIMES)   )
ALLOCATE( tCurrent(1:LB_NTIMES) )
!ALLOCATE( LoadSum(1:14)  )
!  1 -tDG
!     2 -tDGComm ! (not used for ElemTime!)
!  3 -tPML/QDS
!     4 -tEmission ! (not used for ElemTime!)
!     5 -tTrack ! (not used for ElemTime!)
!  6 -tInterpolation/RHS
!  7 -tDeposition (only used with shapefunctions for ElemTime!)
!     8 -tDSMC ! (not used for ElemTime!)
!  9 -tPush
!    10 -tPartComm ! (not used for ElemTime!)
!    11 -tSplit&Merge ! (not used for ElemTime!)
! 12 -UNFP
! 13 -DGAnalyze
! 14 -PartAnalyze


tCartMesh=0.
tTracking=0.

tTotal=0.
!LoadSum=0.
tCurrent=0.

nTotalParts=0 
!nLoadIter  =0

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
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nDeposPerElem,nTracksPerElem,tTracking,tCartMesh
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
  ! time measurement over whole dt_analyze 
  tTotal   = tTotal+tCurrent ! Moved from LoadMeasure
  tCurrent = 0.

  ! number of load balance calls to Compute Elem Load
  nLoadBalance=nLoadBalance+1

  ! time per dg elem
  tDG=(tTotal(LB_DG)+tTotal(LB_DGANALYZE))/REAL(PP_nElems)
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

  helpSum=SUM(nPartsPerElem)
  IF(helpSum.GT.0) THEN
    stotalParts=1.0/REAL(helpSum)
  ELSE
    stotalParts=1.0/REAL(PP_nElems)
    nPartsPerElem=1
  END IF
  tParts = tTotal(LB_INTERPOLATION)+tTotal(LB_PUSH)+tTotal(LB_UNFP)+tTotal(LB_DSMC)!+tTotal(LB_PARTANALYZE) ! interpolation+unfp+analyze
  IF(DoRefMapping)THEN
    helpSum=SUM(nTracksPerElem)
    IF(SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks=1.0/REAL(helpSum)
    ELSE
      stotalTracks=1.0/REAL(PP_nElems)
      nTracksPerElem=1
    END IF
  ELSE IF(TriaTracking)THEN
    tTracking = tTotal(LB_TRACK)
    helpSum=SUM(nTracksPerElem)
    IF(SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks=1.0/REAL(helpSum)
    ELSE
      stotalTracks=1.0/REAL(PP_nElems)
      nTracksPerElem=1
    END IF
  ELSE
    stotalTracks=1.0/REAL(PP_nElems)
    nTracksPerElem=1
  END IF
  helpSum=SUM(nDeposPerElem)
  IF(helpSum.GT.0) THEN
    stotalDepos=1.0/REAL(helpSum)
  END IF
  helpSum=SUM(nSurfacefluxPerElem)
  IF(helpSum.GT.0) THEN
    stotalSurfacefluxes=1.0/REAL(helpSum)
  END IF
  helpSum=SUM(nPartsPerBCElem)
  IF(helpSum.GT.0) THEN
    stotalBCParts=1.0/REAL(helpSum)
  ELSE
    stotalBCParts=1.0/REAL(PP_nElems)
    nPartsPerBCElem=1
  END IF
#endif /*PARTICLES*/

  DO iElem=1,PP_nElems
    ElemTime(iElem) = ElemTime(iElem) + tDG
    !IF(ElemTime(iElem).GT.1000) THEN
    !  IPWRITE(*,*) 'ElemTime already above 1000'
    !END IF
#ifndef PP_HDG
    IF(DoPML)THEN
      IF(ElemToPML(iElem).GT.0 ) ElemTime(iElem) = ElemTime(iElem) + tPML
    END IF
#endif /*PP_HDG*/

#ifdef PARTICLES
    !IF(tParts * nPartsPerElem(iElem)*sTotalParts.GT.1000)THEN
    !  IPWRITE(*,*) 'tParts above 1000',tParts * nPartsPerElem(iElem)*sTotalParts,nPartsPerElem(iElem),sTotalParts,1.0/sTotalParts
    !END IF
    !IF(tTracking * nTracksPerElem(iElem)*sTotalTracks.GT.1000)THEN
    !  IPWRITE(*,*) 'tTracking above 1000',tTracking * nTracksPerElem(iElem)*sTotalTracks,&
    !                                     nTracksPerElem(iElem),sTotalTracks,1.0/sTotalTracks
    !END IF
    ElemTime(iElem) = ElemTime(iElem)                              &
        + tTotal(LB_ADAPTIVE) * nPartsPerBCElem(iElem)*sTotalBCParts &
        + tParts * nPartsPerElem(iElem)*sTotalParts    &
        + tCartMesh * nPartsPerElem(iElem)*sTotalParts &
        + tTracking * nTracksPerElem(iElem)*sTotalTracks &
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

tCartMesh  =0.
tTracking  =0.
#endif /*PARTICLES*/

!tTotal     =0.
!LoadSum    =0.
!tCurrent   =0.
!nTotalParts=0.
!nLoadIter  =0

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


!   SUBROUTINE LoadMeasure() 
!   !----------------------------------------------------------------------------------------------------------------------------------!
!   ! Performs the load measure stuff
!   !----------------------------------------------------------------------------------------------------------------------------------!
!   ! MODULES                                                                                                                          !
!   !----------------------------------------------------------------------------------------------------------------------------------!
!   USE MOD_Globals
!   USE MOD_Preproc
!   USE MOD_LoadBalance_Vars,       ONLY:tCurrent,LoadSum,tTotal,nLoadIter
!   #ifdef PARTICLES
!   USE MOD_LoadBalance_Vars,       ONLY:nTotalParts
!   USE MOD_Particle_Tracking_Vars, ONLY:nCurrentParts
!   #endif /*PARTICLES*/
!   #ifndef PP_HDG
!   USE MOD_PML_Vars,               ONLY:DoPML,nPMLElems
!   #endif /*PP_HDG*/
!   #if defined(LSERK) || defined(IMPA) || defined(ROS)
!   #if (PP_TimeDiscMethod!=110)
!   USE MOD_TimeDisc_Vars,          ONLY:nRKStages
!   #endif
!   #endif
!   !----------------------------------------------------------------------------------------------------------------------------------!
!   IMPLICIT NONE
!   ! INPUT VARIABLES 
!   !----------------------------------------------------------------------------------------------------------------------------------!
!   ! OUTPUT VARIABLES
!   !-----------------------------------------------------------------------------------------------------------------------------------
!   ! LOCAL VARIABLES
!   REAL                :: nLocalParts
!   !===================================================================================================================================
!   
!   nloadIter=nloaditer+1
!   tTotal=tTotal+tCurrent
!   
!   #ifdef PARTICLES
!   #if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
!       nLocalParts=REAL(nCurrentParts)/REAL(nRKStages) ! parts per stage
!   #else
!       nLocalParts=REAL(nCurrentParts)
!   #endif /*TimeDiscMethod*/
!     nTotalParts=nTotalParts+nLocalParts
!   #else
!     nLocalParts=0
!   #endif /*PARTICLES*/
!   
!   ! compute load sum
!   ! dg
!   LoadSum(LB_DG)    =LoadSum(LB_DG)    +tCurrent(LB_DG)    /REAL(PP_nElems)
!   LoadSum(LB_DGCOMM)=LoadSum(LB_DGCOMM)+tCurrent(LB_DGCOMM)/REAL(PP_nElems)
!   #ifndef PP_HDG
!   IF(DoPML)THEN
!     IF(nPMLElems.GT.0) LoadSum(LB_PML)=LoadSum(LB_PML)+tCurrent(LB_PML)/REAL(nPMLElems)
!   END IF
!   #endif /*PP_HDG*/
!   
!   LoadSum(LB_DGANALYZE)=LoadSum(LB_DGANALYZE)+tCurrent(LB_DGANALYZE)/REAL(PP_nElems)
!   
!   #ifdef PARTICLES
!   ! particles
!   IF(nLocalParts.GT.0)THEN
!     nLocalParts=1.0/nLocalParts
!     LoadSum(LB_EMISSION:LB_UNFP)=LoadSum(LB_EMISSION:LB_UNFP)+tCurrent(LB_EMISSION:LB_UNFP)*nLocalParts
!     LoadSum(LB_PARTANALYZE)=LoadSum(LB_PARTANALYZE)+tCurrent(LB_PARTANALYZE)*nLocalParts
!   END IF
!   #endif /*PARTICLES*/
!   
!   ! last operation
!   tCurrent=0.
!   #ifdef PARTICLES
!   nCurrentParts=0
!   #endif /*PARTICLES*/
!   
!   END SUBROUTINE LoadMeasure


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
!SDEALLOCATE( LoadSum )
InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
