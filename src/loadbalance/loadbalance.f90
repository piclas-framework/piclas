#include "boltzplatz.h"

MODULE MOD_LoadBalance
!===================================================================================================================================
! Module contains the routines for load balancing
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#ifdef MPI
INTERFACE InitLoadBalance
  MODULE PROCEDURE InitLoadBalance
END INTERFACE

INTERFACE FinalizeLoadBalance
  MODULE PROCEDURE FinalizeLoadBalance
END INTERFACE

!INTERFACE SingleStepOptimalPartition
!  MODULE PROCEDURE SingleStepOptimalPartition
!END INTERFACE

INTERFACE ComputeElemLoad
  MODULE PROCEDURE ComputeElemLoad
END INTERFACE

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

INTERFACE LoadMeasure 
  MODULE PROCEDURE LoadMeasure
END INTERFACE

PUBLIC::InitLoadBalance,FinalizeLoadBalance,LoadBalance,LoadMeasure
PUBLIC::ComputeElemLoad
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

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

DoLoadBalance= GETLOGICAL('DoLoadBalance','F')
DeviationThreshold  = GETREAL('Load-DeviationThreshold','0.10')
!DeviationThreshold  = 1.0+DeviationThreshold
nLoadBalance = 0
nLoadBalanceSteps = 0

ParticleMPIWeight = GETREAL('Particles-MPIWeight','0.02')
IF (ParticleMPIWeight.LT.0) THEN
  CALL abort(&
    __STAMP__&
    ,' ERROR: Particle weight cannot be negative!')
END IF

PartWeightMethod  = GETINT('Particles-WeightMethod','1')
WeightAverageMethod = GETINT('Particles-WeightAverageMethod','1')
IF ( (WeightAverageMethod.NE.1) .AND. (WeightAverageMethod.NE.2) ) THEN
  CALL abort(&
    __STAMP__&
    ,' ERROR: WeightAverageMethod must be 1 (per iter) or 2 (per dt_analyze)!')
END IF

ALLOCATE( tTotal(1:14)    &
        , tCurrent(1:14)  &
        , LoadSum(1:14)   )
!  1 -tDG
!  2 -tDGComm
!  3 -tPML
!  4 -tEmission
!  5 -tTrack
!  6 -tInterpolation
!  7 -tDeposition
!  8 -tDSMC
!  9 -tPush
! 10 -tPartComm
! 11 -tSplit&Merge
! 12 -UNFP
! 13 -DGAnalyze
! 14 -PartAnalyze


tCartMesh=0.
tTracking=0.

tTotal=0.
LoadSum=0.
tCurrent=0.

nTotalParts=0 
nLoadIter  =0

InitLoadBalanceIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance



SUBROUTINE ComputeElemLoad(CurrentImbalance,PerformLoadbalance,time) 
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the element load
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,        ONLY:ElemTime,nLoadBalance,tTotal,tCurrent
#ifndef PP_HDG
USE MOD_PML_Vars,                ONLY:DoPML,nPMLElems,ElemToPML
#endif /*PP_HDG*/
USE MOD_LoadBalance_Vars,        ONLY:DeviationThreshold,LoadSum,nLoadIter
#ifdef PARTICLES
USE MOD_LoadBalance_Vars,        ONLY:nPartsPerElem,nDeposPerElem,nTracksPerElem,tTracking,tCartMesh
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_PICDepo_Vars,            ONLY:DepositionType
#endif /*PARTICLES*/
USE MOD_LoadBalance_Vars,        ONLY:TargetWeight

!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
REAL ,INTENT(IN)             :: time
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: PerformLoadBalance
REAL ,INTENT(OUT)             :: Currentimbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem,ioUnit
REAL                  :: tDG, tPML
#ifdef PARTICLES
INTEGER(KIND=8)       :: HelpSum
REAL                  :: stotalDepos,stotalParts,sTotalTracks
REAL                  :: tParts
#endif /*PARTICLES*/
REAL                  :: MaxWeight, MinWeight
CHARACTER(LEN=255)    :: outfile
!===================================================================================================================================

nLoadBalance=nLoadBalance+1

! time per dg elem
tDG=(tTotal(1)+tTotal(13))/REAL(PP_nElems)
tPML=0.
#ifndef PP_HDG
IF(DoPML)THEN
  IF(nPMLElems.GT.0) tPML=tTotal(3)/REAL(nPMLElems)
END IF
#endif /*PP_HDG*/

#ifdef PARTICLES
stotalDepos=1.0
sTotalTracks=1.0

helpSum=SUM(nPartsPerElem)
IF(helpSum.GT.0) THEN
  stotalParts=1.0/REAL(helpSum)
ELSE
  stotalParts=1.0/REAL(PP_nElems)
  nPartsPerElem=1
END IF
tParts = tTotal(6)+tTotal(9)+tTotal(12)+tTotal(14) ! interpolation+unfp+analyze
IF(DoRefMapping)THEN
  helpSum=SUM(nTracksPerElem)
  IF(SUM(nTracksPerEleM).GT.0) THEN
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
                  + tParts * nPartsPerElem(iElem)*sTotalParts    &
                  + tCartMesh * nPartsPerElem(iElem)*sTotalParts &
                  + tTracking * nTracksPerElem(iElem)*sTotalTracks
  IF(   (TRIM(DepositionType).EQ.'shape_function')             &
   .OR. (TRIM(DepositionType).EQ.'shape_function_1d')          &    
   .OR. (TRIM(DepositionType).EQ.'shape_function_cylindrical') &    
   .OR. (TRIM(DepositionType).EQ.'shape_function_simple')      &    
   .OR. (TRIM(DepositionType).EQ.'shape_function_spherical') )THEN
    !IF(tTotal(7) * nDeposPerElem(iElem)*sTotalDepos.GT.1000)THEN
    !  IPWRITE(*,*) 'deposition above 1000',tTotal(7) * nDeposPerElem(iElem)*sTotalDepos,nDeposPerElem(iElem)& 
    !                                      ,sTotalDepos,1.0/sTotalDepos
    !END IF
    ElemTime(iElem) = ElemTime(iElem)                              &
                    + tTotal(7) * nDeposPerElem(iElem)*stotalDepos 
  END IF
#endif /*PARTICLES*/
END DO ! iElem=1,PP_nElems

CALL ComputeImbalance(CurrentImbalance,MaxWeight,MinWeight,ElemTime)

! Fill .csv file for parformance analysis and load blaaaance
IF(MPIRoot)THEN
  outfile='ElemTimeStatistics.csv'
  IF(FILEEXISTS(outfile))THEN
    ioUnit=GETFREEUNIT()
    OPEN(UNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') time
    WRITE(ioUnit,'(A1,ES25.10)',ADVANCE='NO') ',',MinWeight 
    WRITE(ioUnit,'(A1,ES25.10)',ADVANCE='NO') ',',MaxWeight
    WRITE(ioUnit,'(A1,ES25.10)',ADVANCE='NO') ',',CurrentImbalance
    WRITE(ioUnit,'(A1,ES25.10)',ADVANCE='NO') ',',TargetWeight  
    WRITE(ioUnit,'(A1)') ' '
    CLOSE(ioUnit) 
  END IF
END IF !IF(MPIRoot)


PerformLoadBalance=.FALSE.
! only check if imbalance is > a given threshold
IF(CurrentImbalance.GT.DeviationThreshold) PerformLoadBalance=.TRUE.

#ifdef PARTICLES
nTracksPerElem=0
nDeposPerElem=0
nPartsPerElem=0


tCartMesh  =0.
tTracking  =0.
#endif /*PARTICLES*/
tTotal     =0.
LoadSum    =0.
tCurrent   =0.
!nTotalParts=0.
nLoadIter  =0

END SUBROUTINE ComputeElemLoad


SUBROUTINE LoadBalance(CurrentImbalance,PerformLoadBalance)
!===================================================================================================================================
! routine perfoming the load balancing stuff
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Restart,               ONLY:Restart
USE MOD_Boltzplatz_Tools,      ONLY:InitBoltzplatz,FinalizeBoltzplatz
USE MOD_LoadBalance_Vars,      ONLY:ElemTime,nLoadBalanceSteps
#ifdef PARTICLES
USE MOD_PICDepo_Vars,          ONLY:DepositionType
USE MOD_Particle_MPI,          ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)   :: PerformLoadBalance
REAL,INTENT(IN)      :: Currentimbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: Newimbalance, MaxWeight, MinWeight
!===================================================================================================================================


! only do load-balance if necessary
IF(.NOT.PerformLoadBalance) THEN
  ElemTime=0.
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'

nLoadBalanceSteps=nLoadBalanceSteps+1
! finialize all arrays
CALL FinalizeBoltzplatz(IsLoadBalance=.TRUE.)
! reallocate
CALL InitBoltzplatz(IsLoadBalance=.TRUE.)

! restart
CALL Restart()

! compute imbalance 
CALL ComputeImbalance(NewImbalance,MaxWeight,MinWeight,ElemTime)
! zero ElemTime, the measurement starts again
ElemTime=0.

IF( NewImbalance.GT.CurrentImbalance ) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' OldImbalance: ', CurrentImbalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' NewImbalance: ', NewImbalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MaxWeight:    ', MaxWeight
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MinWeight: '   , MinWeight
END IF
 
!CurrentImbalance=NewImbalance

#ifdef PARTICLES
IF(   (TRIM(DepositionType).EQ.'shape_function')             &
 .OR. (TRIM(DepositionType).EQ.'shape_function_1d')          &    
 .OR. (TRIM(DepositionType).EQ.'shape_function_cylindrical') &    
 .OR. (TRIM(DepositionType).EQ.'shape_function_simple')      &    
 .OR. (TRIM(DepositionType).EQ.'shape_function_spherical') )THEN
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

SWRITE(UNIT_stdOut,'(A)')' LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE LoadBalance


SUBROUTINE LoadMeasure() 
!----------------------------------------------------------------------------------------------------------------------------------!
! Performs the load measure stuff
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,       ONLY:tCurrent,LoadSum,tTotal,nLoadIter
#ifdef PARTICLES
USE MOD_LoadBalance_Vars,       ONLY:nTotalParts
USE MOD_Particle_Tracking_Vars, ONLY:nCurrentParts
#endif /*PARTICLES*/
#ifndef PP_HDG
USE MOD_PML_Vars,               ONLY:DoPML,nPMLElems
#endif /*PP_HDG*/
#if defined(LSERK) || defined(IMPA) || defined(IMEX)
#if (PP_TimeDiscMethod!=110)
USE MOD_TimeDisc_Vars,          ONLY:nRKStages
#endif
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: nLocalParts
!===================================================================================================================================

nloadIter=nloaditer+1
tTotal=tTotal+tCurrent

#ifdef PARTICLES
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
  nLocalParts=REAL(nCurrentParts)/REAL(nRKStages) ! parts per stage
#else
  nLocalParts=REAL(nCurrentParts)
#endif /*TimeDiscMethod*/
nTotalParts=nTotalParts+nLocalParts
#else
  nlocalParts=0
#endif /*PARTICLES*/

! compute load sum
! dg
LoadSum(1:2)=LoadSum(1:2)+tCurrent(1:2)/REAL(PP_nElems)
#ifndef PP_HDG
IF(DoPML)THEN
  IF(nPMLElems.GT.0) LoadSum(3)=LoadSum(3)+tCurrent(3)/REAL(nPMLElems)
END IF
#endif /*PP_HDG*/

LoadSum(13)=LoadSum(13)+tCurrent(13)/REAL(PP_nElems)

! particles
IF(nLocalParts.GT.0)THEN
  nLocalParts=1.0/nLocalParts
  LoadSum(4:12)=LoadSum(4:12)+tCurrent(4:12)*nLocalParts
  LoadSum(14)=LoadSum(14)+tCurrent(14)*nLocalParts
END IF

! last operation
tCurrent=0.
#ifdef PARTICLES
nCurrentParts=0
#endif /*PARTICLES*/

END SUBROUTINE LoadMeasure


SUBROUTINE ComputeImbalance(CurrentImbalance,MaxWeight,MinWeight,ElemTime)
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
! Maxweight: maximum weight of all processes
! Minweight: minimum weight of all processes
! targetweight: current target weight
! CurrentImbalance: rel. deviation of current imbalance
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc,             ONLY:PP_nElems
USE MOD_LoadBalance_Vars,    ONLY:WeightSum, TargetWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)            :: ElemTime(1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: CurrentImbalance, MaxWeight, MinWeight
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                       :: MaxWeight
!===================================================================================================================================

WeightSum=SUM(ElemTime)

IF(ALMOSTZERO(WeightSum))THEN
  IPWRITE(*,*) 'Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)'
END IF
CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(WeightSum,MaxWeight,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(WeightSum,MinWeight,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)

TargetWeight=TargetWeight/nProcessors

! new computation of current imbalance
!CurrentImbalance = MAX( MaxWeight-TargetWeight, ABS(MinWeight-TargetWeight)  )/ TargetWeight
CurrentImbalance =  (MaxWeight-TargetWeight ) / TargetWeight
!CurrentImbalance=MaxWeight/TargetWeight

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
SDEALLOCATE( LoadSum )
InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
