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

INTERFACE ComputeParticleWeightAndLoad
  MODULE PROCEDURE ComputeParticleWeightAndLoad
END INTERFACE

INTERFACE ComputeElemLoad
  MODULE PROCEDURE ComputeElemLoad
END INTERFACE

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

INTERFACE CalculateProcWeights
  MODULE PROCEDURE CalculateProcWeights
END INTERFACE

INTERFACE LoadMeasure 
  MODULE PROCEDURE LoadMeasure
END INTERFACE

PUBLIC::InitLoadBalance,FinalizeLoadBalance,LoadBalance,LoadMeasure,CalculateProcWeights
PUBLIC::ComputeParticleWeightAndLoad,ComputeElemLoad
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
DeviationThreshold  = 1.0+DeviationThreshold
nLoadBalance = 0

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

LastImbalance=0.
tTotal=0.
LoadSum=0.
tCurrent=0.

nTotalParts=0 
nLoadIter  =0

InitLoadBalanceIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance


SUBROUTINE ComputeParticleWeightAndLoad(CurrentImbalance,PerformLoadbalance) 
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the current particle weight
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,      ONLY:tCurrent,LoadSum,tTotal,nloaditer,nTotalParts,nLoadBalance,PartWeightMethod,DeviationThreshold&
                                   ,WeightAverageMethod,ParticleMPIWeight,LastImbalance,ElemWeight
#ifndef PP_HDG
USE MOD_PML_Vars,              ONLY:DoPML,nPMLElems
#endif /*PP_HDG*/
USE MOD_Utils,                 ONLY:InsertionSort
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: PerformLoadBalance
REAL ,INTENT(OUT)             :: Currentimbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: TotalLoad(1:14)
REAL,ALLOCATABLE  :: GlobalLoad(:,:,:), GlobalWeights(:,:,:)
REAL              :: totalRatio(2,0:nProcessors-1),WeightList(1:nProcessors)
CHARACTER(LEN=64) :: filename, firstLine(2)
CHARACTER(LEN=4)  :: hilf
INTEGER           :: iounit,iProc,iOut,iWeigth,iList,listZeros
INTEGER           :: nWeights(2,3)
REAL              :: PartWeight(3,2)
REAL              :: MaxWeight, MinWeight !dummy so far, used in timedisc-call

!===================================================================================================================================

nLoadBalance=nLoadBalance+1

! per iter
! finish load measure
LoadSum=LoadSum/REAL(nLoadIter)

TotalLoad=0.
! per dt_analyze
! dg
TotalLoad(1:2)=tTotal(1:2)/(REAL(nloaditer)*REAL(PP_nElems))
#ifndef PP_HDG
IF(DoPML)THEN
  IF(nPMLElems.GT.0) TotalLoad(3)=tTotal(3)/(REAL(nPMLElems)*REAL(nloaditer))
END IF
#endif /*PP_HDG*/
TotalLoad(13)=tTotal(13)/(REAL(nloaditer)*REAL(PP_nElems))

! particles
IF(nTotalParts.GT.0)THEN
  TotalLoad(4:12)=tTotal(4:12)/nTotalParts
  TotalLoad(14)=tTotal(14)/nTotalParts
ELSE
  TotalLoad(4:12)=0.
  TotalLoad(14)=0.
END IF

  ! communication to root
  IF(MPIRoot)THEN
    ALLOCATE(GlobalLoad(1:2,1:14,0:nProcessors-1))
    GlobalLoad=0.
  ELSE
    ALLOCATE(GlobalLoad(2,1,0))
  END IF
  CALL MPI_GATHER(LoadSum  ,14,MPI_DOUBLE_PRECISION,GlobalLoad(1,:,:),14,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
  CALL MPI_GATHER(TotalLoad,14,MPI_DOUBLE_PRECISION,GlobalLoad(2,:,:),14,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)

! no mpi
!  ALLOCATE(GlobalLoad(1:2,1:13,0))
!  GlobalLoad(1,:,:) = LoadSum
!  GlobalLoad(2,:,:) = TotalLoad

! compute weights
IF(MPIRoot)THEN
  PartWeight=0.
  nWeights=0
  ALLOCATE(GlobalWeights(1:2,1:5,0:nProcessors-1))
  ! 1 - pure DG time
  ! 2 - pure particle time
  ! 3 - ratio t_particle / t_DG
  ! 4 - ratio communication time
  ! 5 - total ratio (incl. comm+analyze)
  GlobalWeights=0.
  DO iOut=1,2 !WeightAverageMethod
    DO iProc=0,nProcessors-1
      GlobalWeights(iOut,1,iProc) = GlobalLoad(iOut,1 ,iProc)+GlobalLoad(iOut,3 ,iProc) 
      GlobalWeights(iOut,2,iProc) = SUM(GlobalLoad(iOut,4:9,iProc))+SUM(GlobalLoad(iOut,11:12,iProc))
      IF(.NOT.ALMOSTZERO(GlobalWeights(iOut,1,iProc)))THEN
        GlobalWeights(iOut,3,iProc) =  GlobalWeights(iOut,2,iProc) /  GlobalWeights(iOut,1,iProc)
        PartWeight(1,iOut)=PartWeight(1,iOut)+GlobalWeights(iOut,3,iProc)
        nWeights(iOut,1)=nWeights(iOut,1)+1
      END IF
      IF(.NOT.ALMOSTZERO(GlobalLoad(iOut,2,iProc))) THEN
        GlobalWeights(iOut,4,iProc) = GlobalLoad(iOut,10,iProc) / GlobalLoad(iOut,2,iProc)
        PartWeight(2,iOut)=PartWeight(2,iOut)+GlobalWeights(iOut,4,iProc)
        nWeights(iOut,2)=nWeights(iOut,2)+1
      END IF
      IF(.NOT.ALMOSTZERO(SUM(GlobalLoad(iOut,1:3,iProc))+GlobalLoad(iOut,13,iProc))) THEN
        GlobalWeights(iOut,5,iProc) = (SUM(GlobalLoad(iOut,4:12,iProc))+GlobalLoad(iOut,14,iProc)) &
                                     / (SUM(GlobalLoad(iOut,1:3,iProc))+GlobalLoad(iOut,13,iProc))
        PartWeight(3,iOut)=PartWeight(3,iOut)+GlobalWeights(iOut,5,iProc)
        nWeights(iOut,3)=nWeights(iOut,3)+1
      END IF
    END DO ! iProc
    DO iProc=0,nProcessors-1
      IF(ALMOSTZERO(GlobalWeights(iOut,5,iProc)))THEN
        totalRatio(iOut,iProc)=HUGE(0.)
      ELSE
        totalRatio(iOut,iProc)=GlobalWeights(iOut,5,iProc)
      END IF
    END DO ! iProc
    DO iWeigth=1,3
      IF (nWeights(iOut,iWeigth).NE.0) THEN
        PartWeight(iWeigth,iOut)=PartWeight(iWeigth,iOut)/nWeights(iOut,iWeigth)
      ELSE !iWeigth=0, for safety...
        PartWeight(iWeigth,iOut)=0.
      END IF
    END DO
  END DO ! iOut
  !ParticleMPIWeight=MINVAL(totalRatio(WeightAverageMethod,:)) ! currently per analyze_dt

  SELECT CASE(PartWeightMethod)
  CASE(1)
!    ParticleMPIWeight=PartWeight(3,WeightAverageMethod)
    WeightList(1:nProcessors)=GlobalWeights(WeightAverageMethod,5,0:nProcessors-1)
    listZeros=0
    ParticleMPIWeight=0.
    DO iList=1,nProcessors
      IF(ALMOSTZERO(WeightList(iList)))THEN
        listZeros=listZeros+1
      ELSE
        ParticleMPIWeight=ParticleMPIWeight+WeightList(iList)
      END IF
    END DO
    ParticleMPIWeight=ParticleMPIWeight/(nProcessors-listZeros)
  CASE(2)
    WeightList(1:nProcessors)=GlobalWeights(WeightAverageMethod,5,0:nProcessors-1)
    CALL InsertionSort(a=WeightList,len=nProcessors)
    listZeros=0
    DO iList=1,nProcessors
      IF(ALMOSTZERO(WeightList(iList)))THEN
        listZeros=iList
      ELSE
        EXIT
      END IF
    END DO
!SWRITE(*,*) 'listZeros: ',listZeros
    IF (listZeros.EQ.nProcessors) THEN
      ParticleMPIWeight=0.
    ELSE IF (listZeros.EQ.nProcessors-1) THEN
      ParticleMPIWeight=WeightList(nProcessors)
    ELSE IF(MOD(nProcessors-listZeros,2).EQ.0)THEN
      ParticleMPIWeight=0.5*( WeightList((nProcessors-listZeros)/2 + listZeros) &
                             +WeightList((nProcessors-listZeros)/2 + 1 + listZeros) )
    ELSE
      ParticleMPIWeight=WeightList((nProcessors-listZeros)/2 + 1 + listZeros)
    END IF
  !CASE(3)
  CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,' No valid PartWeightMethod defined!')
  END SELECT
END IF

  CALL MPI_BCAST (ParticleMPIWeight,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
 
! write load info to file
IF(MPIRoot)THEN
  firstLine(1)='Averaged Load per Iter Information:'
  firstLine(2)='Load per analyze-dt Information:'
  WRITE( hilf,'(I4.4)') nLoadBalance
  DO iOut=1,2
    IF(iOut.EQ.1)THEN
      filename='LoadAveragedPerIter_'//TRIM(hilf)//'.csv'
    ELSE
      filename='LoadAveragedPerAnalyze'//TRIM(hilf)//'.csv'
    END IF
    ioUnit=GETFREEUNIT()
    OPEN(UNIT=ioUnit,FILE=filename,STATUS='REPLACE')
    WRITE(ioUnit,'(A)') firstLine(iOut)
    WRITE(ioUnit,'(A,I10)')'total number of Procs,',nProcessors
    WRITE(ioUnit,'(A20,ES17.5)') 'PartMPIWeight:  ', ParticleMPIWeight
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-PartWeight: ', PartWeight(1,iOut)
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-Comm-Weight:', PartWeight(2,iOut)
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-TotalWeight:', PartWeight(3,iOut)
    WRITE(ioUnit,'(6(A20))')'Rank','DG','Part','Ratio-Part-DG','Ratio-Comm','Ratio-total'
    WRITE(ioUnit,'(A80,A40)')&
          '================================================================================',&
          '========================================'
    DO iProc=0,nProcessors-1
      WRITE(ioUnit,'(5X,I10,5x,5(3x,ES17.5))') iProc,GlobalWeights(iOut,1:5,iProc)
      WRITE(ioUnit,'(A80,A40)')&
            '--------------------------------------------------------------------------------',&
            '----------------------------------------'
    END DO ! iProc
    WRITE(ioUnit,'(A)') ''
    WRITE(ioUnit,'(A)') ''
    WRITE(ioUnit,'(15(A20))')'Rank','DG','DGComm','PML','Part-Emission','Part-Tracking','PIC-Inter','PIC-Depo',&
                             'DSMC','Push','PartComm' ,'SplitMerge','UNFP','DG-Analyze','Part-Analyze'
    WRITE(ioUnit,'(A120,A120,A60)')&
      '========================================================================================================================',&
      '========================================================================================================================',&
      '============================================================='
    DO iProc=0,nProcessors-1
      WRITE(ioUnit,'(5X,I10,5x,14(3x,ES17.5))') iProc,GlobalLoad(iOut,:,iProc)
      WRITE(ioUnit,'(A120,A120,A60)')&
      '------------------------------------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------------------------------------------------------------------',& 
      '-------------------------------------------------------------'
    END DO ! iProc
    CLOSE(ioUnit) 
  END DO ! iOut

END IF

DEALLOCATE(GlobalLoad)
SDEALLOCATE(GlobalWeights)
TotalLoad  =0.
tTotal     =0.
LoadSum    =0.
tCurrent   =0.
nTotalParts=0.
nLoadIter  =0

! compute current load
CALL CalculateProcWeights()
! compute impalance 
CALL ComputeImbalance(CurrentImbalance,MaxWeight,MinWeight,ElemWeight)

PerformLoadBalance=.FALSE.
IF(CurrentImbalance.GT.DeviationThreshold*LastImbalance) PerformLoadBalance=.TRUE.


END SUBROUTINE ComputeParticleWeightAndLoad


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
USE MOD_LoadBalance_Vars,        ONLY:DeviationThreshold,LastImbalance,LoadSum,nLoadIter
#ifdef PARTICLES
USE MOD_LoadBalance_Vars,        ONLY:nPartsPerElem,nDeposPerElem,nTracksPerElem,tTracking,tCartMesh
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_PICDepo_Vars,            ONLY:DepositionType
#endif /*PARTICLES*/
USE MOD_Particle_Analyze_Vars,   ONLY:IsRestart
USE MOD_LoadBalance_Vars,        ONLY:TargetWeight,WeightOutput
USE MOD_TimeDisc_Vars,           ONLY:iter

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
INTEGER(KIND=8)       :: HelpSum
#ifdef PARTICLES
REAL                  :: stotalDepos,stotalParts,sTotalTracks
REAL                  :: tParts
#endif /*PARTICLES*/
REAL                  :: MaxWeight, MinWeight
LOGICAL               :: isOpen, FileExists
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
  WeightOutput(3) = (MaxWeight - WeightOutput(2))/(TargetWeight - WeightOutput(4)) ! real current imbalance (no average)
  outfile='ElemTimeStatistics.csv'
  INQUIRE(FILE=TRIM(outfile),EXIST=FileExists)
  !IF (isRestart .and. FileExists) THEN
  IF(FileExists)THEN
    ioUnit=GETFREEUNIT()
    OPEN(UNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') time
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') MinWeight        - WeightOutput(1)
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') MaxWeight        - WeightOutput(2)
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') WeightOutput(3)
    WRITE(ioUnit,'(ES25.10)',ADVANCE='NO') TargetWeight     - WeightOutput(4)
    WRITE(ioUnit,'(A1)') ' '
    CLOSE(ioUnit) 
  END IF
  WeightOutput(1) = MinWeight
  WeightOutput(2) = MaxWeight
  WeightOutput(4) = TargetWeight
END IF !IF(MPIRoot)


PerformLoadBalance=.FALSE.
IF(CurrentImbalance.GT.DeviationThreshold*LastImbalance) PerformLoadBalance=.TRUE.

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
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Restart,               ONLY:Restart
!USE MOD_Boltzplatz_Tools,      ONLY:InitBoltzplatz,FinalizeBoltzplatz
USE MOD_Boltzplatz_Tools,      ONLY:InitBoltzplatz,FinalizeBoltzplatz
USE MOD_LoadBalance_Vars,      ONLY:DeviationThreshold,LastImbalance,ElemWeight
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
IF(.NOT.PerformLoadBalance) RETURN

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'

! finialize all arrays
CALL FinalizeBoltzplatz(IsLoadBalance=.TRUE.)
SDEALLOCATE(ElemWeight)
! reallocate
CALL InitBoltzplatz(IsLoadBalance=.TRUE.)

! restart
CALL Restart()

! check new distribution
! compute current load
! read in elemweight
! new method uses read in elemweigfht
!CALL CalculateProcWeights()

!! compute impalance 
CALL ComputeImbalance(NewImbalance,MaxWeight,MinWeight,ElemWeight)

IF( NewImbalance.GT.DeviationThreshold*LastImbalance   &
  .OR. NewImbalance.GT.CurrentImbalance ) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' LastImbalance:    ', LastImBalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' CurrentImbalance: ', CurrentImbalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' NewImbalance: '    , NewImbalance
  !SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MaxWeight:    '    , MaxWeight
  !SWRITE(UNIT_stdOut,'(A25,E15.7)') ' MinWeight: '       , MinWeight
END IF
 
LastImbalance=NewImBalance

#ifdef PARTICLES
IF(   (TRIM(DepositionType).EQ.'shape_function')             &
 .OR. (TRIM(DepositionType).EQ.'shape_function_1d')          &    
 .OR. (TRIM(DepositionType).EQ.'shape_function_cylindrical') &    
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
USE MOD_LoadBalance_Vars,       ONLY:tCurrent,LoadSum,tTotal,nLoadIter,nTotalParts
#ifdef PARTICLES
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


SUBROUTINE CalculateProcWeights() 
!----------------------------------------------------------------------------------------------------------------------------------!
! calculation of weights of each processor
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,  ONLY: ElemWeight,ParticleMPIWeight
#ifdef PARTICLES
USE MOD_Particle_Vars, ONLY: PDM,PEM
#endif /*PARTICLES*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iElem, iPart,locnParts
!===================================================================================================================================

IF(.NOT.ALLOCATED(ElemWeight)) ALLOCATE(ElemWeight(1:PP_nElems))

ElemWeight=1.0

#ifdef PARTICLES
DO iElem=1,PP_nElems ! loop only over internal elems, if particle is already in HALO, it shall not be found here
  locnParts=0
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    IF(PEM%LastElement(iPart).NE.iElem) CYCLE
    locnParts=locnParts+1
  END DO
  ElemWeight(iElem)=ElemWeight(iElem)+ParticleMPIWeight*locnParts
END DO ! iElem
#endif /*PARTICLES*/

END SUBROUTINE CalculateProcWeights


SUBROUTINE ComputeImbalance(CurrentImbalance,MaxWeight,MinWeight,ElemWeight)
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
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
REAL,INTENT(IN)            :: ElemWeight(1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: CurrentImbalance, MaxWeight, MinWeight
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                       :: MaxWeight
!===================================================================================================================================

WeightSum=SUM(ElemWeight)

IF(ALMOSTZERO(WeightSum))THEN
  IPWRITE(*,*) 'Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)'
END IF
CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(WeightSum,MaxWeight,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(WeightSum,MinWeight,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)

TargetWeight=TargetWeight/nProcessors

CurrentImbalance=MaxWeight/TargetWeight

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
SDEALLOCATE( ElemWeight )
InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
