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
PUBLIC::ComputeParticleWeightAndLoad
!===================================================================================================================================

CONTAINS

SUBROUTINE InitLoadBalance()
!===================================================================================================================================
! init load balancing, new initialization of variables for load balancing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars
USE MOD_ReadInTools,          ONLY:GETLOGICAL, GETREAL, GETINT
USE MOD_Analyze_Vars,         ONLY:Analyze_dt
#ifdef MPI
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

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE ...'

DoLoadBalance= GETLOGICAL('Static-LoadBalance','F')
DeviationThreshold  = GETREAL('Load-DeviationTreshold','0.10')
DeviationThreshold  = 1.0+DeviationThreshold
OutputRank= GETLOGICAL('OutputRank','F')
nLoadBalance = 0

ParticleMPIWeight = GETREAL('Particles-MPIWeight','0.02')
IF (ParticleMPIWeight.LT.0) THEN
  CALL abort(&
    __STAMP__,' ERROR: Particle weight cannot be negative!')
END IF

PartWeightMethod  = GETINT('Particles-WeightMethod','1')
WeightAverageMethod = GETINT('Particles-WeightAverageMethod','2')
IF ( (WeightAverageMethod.NE.1) .AND. (WeightAverageMethod.NE.2) ) THEN
  CALL abort(&
    __STAMP__&
    ,' ERROR: WeightAverageMethod must be 1 (per iter) or 2 (per dt_analyze)!')
END IF

ALLOCATE( tTotal(1:13)    &
        , tCurrent(1:13)  &
        , LoadSum(1:13)   )
!  1 -tDG
!  2 -tDGComm
!  3 -tPML
!  4 -tEmission
!  5 -tTrack
!  6 -tPIC
!  7 -tDSMC
!  8 -tPush
!  9 -tPartComm
! 10 -tSplit&Merge
! 11 -UNFP
! 12 -DGAnalyze
! 13 -PartAnalyze


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
                                   ,WeightAverageMethod,ParticleMPIWeight,LastImbalance
USE MOD_PML_Vars,              ONLY:DoPML,nPMLElems
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
REAL              :: TotalLoad(1:13)
REAL,ALLOCATABLE  :: GlobalLoad(:,:,:), GlobalWeights(:,:,:)
REAL              :: totalRatio(2,0:nProcessors-1),WeightList(1:nProcessors)
CHARACTER(LEN=64) :: filename, firstLine(2)
CHARACTER(LEN=4)  :: hilf
INTEGER           :: iounit,iProc,iOut,iWeigth,iList,listZeros
INTEGER           :: nWeights(2,3)
REAL              :: PartWeight(3,2)

!===================================================================================================================================

nLoadBalance=nLoadBalance+1

! per iter
! finish load measure
LoadSum=LoadSum/REAL(nLoadIter)

TotalLoad=0.
! per dt_analyze
! dg
TotalLoad(1:2)=tTotal(1:2)/(REAL(nloaditer)*REAL(PP_nElems))
IF(DoPML)THEN
  IF(nPMLElems.GT.0) TotalLoad(3)=tTotal(3)/(REAL(nPMLElems)*REAL(nloaditer))
END IF
TotalLoad(12)=tTotal(12)/(REAL(nloaditer)*REAL(PP_nElems))

! particles
IF(nTotalParts.GT.0)THEN
  TotalLoad(4:11)=tTotal(4:11)/nTotalParts
  TotalLoad(13)=tTotal(13)/nTotalParts
ELSE
  TotalLoad(4:11)=0.
  TotalLoad(13)=0.
END IF

#ifdef MPI
  ! communication to root
  IF(MPIRoot)THEN
    ALLOCATE(GlobalLoad(1:2,1:13,0:nProcessors-1))
    GlobalLoad=0.
  ELSE
    ALLOCATE(GlobalLoad(2,1,0))
  END IF
  CALL MPI_GATHER(LoadSum  ,13,MPI_DOUBLE_PRECISION,GlobalLoad(1,:,:),13,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
  CALL MPI_GATHER(TotalLoad,13,MPI_DOUBLE_PRECISION,GlobalLoad(2,:,:),13,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
#else
  ALLOCATE(GlobalLoad(1:2,1:13,0))
  GlobalLoad(1,:,:) = LoadSum
  GlobalLoad(2,:,:) = TotalLoad
#endif /*MPI*/

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
      GlobalWeights(iOut,2,iProc) = SUM(GlobalLoad(iOut,4:8,iProc))+SUM(GlobalLoad(iOut,10:11,iProc))
      IF(.NOT.ALMOSTZERO(GlobalWeights(iOut,1,iProc)))THEN
        GlobalWeights(iOut,3,iProc) =  GlobalWeights(iOut,2,iProc) /  GlobalWeights(iOut,1,iProc)
        PartWeight(1,iOut)=PartWeight(1,iOut)+GlobalWeights(iOut,3,iProc)
        nWeights(iOut,1)=nWeights(iOut,1)+1
      END IF
      IF(.NOT.ALMOSTZERO(GlobalLoad(iOut,2,iProc))) THEN
        GlobalWeights(iOut,4,iProc) = GlobalLoad(iOut,9,iProc) / GlobalLoad(iOut,2,iProc)
        PartWeight(2,iOut)=PartWeight(2,iOut)+GlobalWeights(iOut,4,iProc)
        nWeights(iOut,2)=nWeights(iOut,2)+1
      END IF
      IF(.NOT.ALMOSTZERO(SUM(GlobalLoad(iOut,1:3,iProc))+GlobalLoad(iOut,12,iProc))) THEN
        GlobalWeights(iOut,5,iProc) = (SUM(GlobalLoad(iOut,4:11,iProc))+GlobalLoad(iOut,13,iProc)) &
                                     / (SUM(GlobalLoad(iOut,1:3,iProc))+GlobalLoad(iOut,12,iProc))
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
      __STAMP__,' No valid PartWeightMethod defined!')
  END SELECT
END IF
#ifdef MPI
  CALL MPI_BCAST (ParticleMPIWeight,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
 
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
    WRITE(ioUnit,'(14(A20))')'Rank','DG','DGComm','PML','Part-Emission','Part-Tracking','PIC','DSMC','Push','PartComm' &
                                 ,'SplitMerge','UNFP','DG-Analyze','Part-Analyze'
    WRITE(ioUnit,'(A120,A120,A40)')&
      '========================================================================================================================',&
      '========================================================================================================================',&
      '========================================='
    DO iProc=0,nProcessors-1
      WRITE(ioUnit,'(5X,I10,5x,13(3x,ES17.5))') iProc,GlobalLoad(iOut,:,iProc)
      WRITE(ioUnit,'(A120,A120,A40)')&
      '------------------------------------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------------------------------------------------------------------',& 
      '-----------------------------------------'
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
CALL ComputeImbalance(CurrentImbalance)

PerformLoadBalance=.FALSE.
IF(CurrentImbalance.GT.DeviationThreshold*LastImbalance) PerformLoadBalance=.TRUE.


END SUBROUTINE ComputeParticleWeightAndLoad


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
USE MOD_LoadBalance_Vars,      ONLY:DeviationThreshold,LastImbalance
#ifdef PARTICLES
USE MOD_PICDepo_Vars,          ONLY:DepositionType
#ifdef MPI
USE MOD_Particle_MPI,          ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,     ONLY:PartMPIExchange
#endif
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
REAL              :: Newimbalance
!===================================================================================================================================


! only do load-balance if necessary
IF(.NOT.PerformLoadBalance) RETURN

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'

! finialize all arrays
CALL FinalizeBoltzplatz(IsLoadBalance=.TRUE.)
! reallocate
CALL InitBoltzplatz(IsLoadBalance=.TRUE.)

! restart
CALL Restart()

! check new distribution
! compute current load
CALL CalculateProcWeights()
! compute impalance 
CALL ComputeImbalance(NewImbalance)

IF( NewImbalance.GT.DeviationThreshold*LastImbalance   &
  .OR. NewImbalance.GT.CurrentImbalance ) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' LastImbalance:    ', LastImBalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' CurrentImbalance: ', CurrentImbalance
  SWRITE(UNIT_stdOut,'(A25,E15.7)') ' NewImbalance: '    , NewImbalance
END IF

LastImbalance=NewImBalance

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
USE MOD_Particle_Tracking_Vars, ONLY:nCurrentParts
USE MOD_PML_Vars,               ONLY:DoPML,nPMLElems
USE MOD_TimeDisc_Vars,          ONLY:nRKStages
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

#if ((PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6))  /* RK3 and RK4 only */
  nLocalParts=REAL(nCurrentParts)/REAL(nRKStages) ! parts per stage
#else
  nLocalParts=REAL(nCurrentParts)
#endif

nTotalParts=nTotalParts+nLocalParts

! compute load sum
! dg
LoadSum(1:2)=LoadSum(1:2)+tCurrent(1:2)/REAL(PP_nElems)
IF(DoPML)THEN
  IF(nPMLElems.GT.0) LoadSum(3)=LoadSum(3)+tCurrent(3)/REAL(nPMLElems)
END IF
LoadSum(12)=LoadSum(12)+tCurrent(12)/REAL(PP_nElems)

! particles
IF(nLocalParts.GT.0)THEN
  nLocalParts=1.0/nLocalParts
  LoadSum(4:11)=LoadSum(4:11)+tCurrent(4:11)*nLocalParts
  LoadSum(13)=LoadSum(13)+tCurrent(13)*nLocalParts
END IF

! last operation
tCurrent=0.
nCurrentParts=0

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


SUBROUTINE ComputeImbalance(CurrentImbalance)
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_LoadBalance_Vars,    ONLY:ElemWeight, WeightSum, TargetWeight
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: CurrentImbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                       :: MaxWeight
!===================================================================================================================================

WeightSum=SUM(ElemWeight)

CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(WeightSum,MaxWeight,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)

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

END MODULE MOD_LoadBalance
