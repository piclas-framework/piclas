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

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

INTERFACE LoadMeasure 
  MODULE PROCEDURE LoadMeasure
END INTERFACE

PUBLIC::InitLoadBalance,FinalizeLoadBalance,LoadBalance,LoadMeasure
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
USE MOD_ReadInTools,          ONLY:GETLOGICAL
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
OutputRank= GETLOGICAL('OutputRank','F')
nLoadBalance = 0

ALLOCATE( tTotal(1:11)    &
        , tCurrent(1:11)  &
        , LoadSum(1:11)   )
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


tTotal=0.
LoadSum=0.
tCurrent=0.

nTotalParts=0 
nLoadIter  =0

InitLoadBalanceIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance


SUBROUTINE LoadBalance()
!===================================================================================================================================
! routine perfoming the load balancing stuff
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Restart,               ONLY:Restart
USE MOD_Boltzplatz_Tools,      ONLY:InitBoltzplatz,FinalizeBoltzplatz
USE MOD_PICDepo_Vars,          ONLY:DepositionType
USE MOD_LoadBalance_Vars,      ONLY:tCurrent,LoadSum,tTotal,nloaditer,nTotalParts,nLoadBalance
USE MOD_PML_Vars,              ONLY:DoPML,nPMLElems
USE MOD_Mesh_Vars,             ONLY:ParticleMPIWeight
#ifdef MPI
USE MOD_Particle_MPI,          ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,     ONLY:PartMPIExchange
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: TotalLoad(1:11)
REAL,ALLOCATABLE  :: GlobalLoad(:,:,:), GlobalWeights(:,:,:)
CHARACTER(LEN=64) :: filename, firstLine(2)
CHARACTER(LEN=4)  :: hilf
INTEGER           :: iounit,iProc,iOut
INTEGER           :: nWeights(2)
REAL              :: PartWeight(3,2)
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'

nLoadBalance=nLoadBalance+1

! per iter
! finish load measure
LoadSum=LoadSum/REAL(nloaditer)


TotalLoad=0.
! per dt_analyze
! dg
TotalLoad(1:2)=tTotal(1:2)/(REAL(nloaditer)*REAL(PP_nElems))

IF(DoPML)THEN
  IF(nPMLElems.GT.0) TotalLoad(3)=tTotal(3)/(REAL(nPMLElems)*REAL(nloaditer))
END IF

! particles
IF(nTotalParts.GT.0)THEN
  TotalLoad(4:11)=tTotal(4:11)/REAL(nTotalParts)
ELSE
  TotalLoad(4:11)=0.
END IF

#ifdef MPI
  ! communication to root
  IF(MPIRoot)THEN
    ALLOCATE(GlobalLoad(1:2,1:11,0:nProcessors-1))
    GlobalLoad=0.
  ELSE
    ALLOCATE(GlobalLoad(2,1,0))
  END IF
  CALL MPI_GATHER(LoadSum  ,11,MPI_DOUBLE_PRECISION,GlobalLoad(1,:,:),11,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
  CALL MPI_GATHER(TotalLoad,11,MPI_DOUBLE_PRECISION,GlobalLoad(2,:,:),11,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
#else
  ALLOCATE(GlobalLoad(1:2,1:11,0))
  GlobalLoad(1,:,:) = LoadSum
  GlobalLoad(2,:,:) = TotalLoad
#endif /*MPI*/

! compute weights
IF(MPIRoot)THEN
  PartWeight=0.
  nWeights=0
  ALLOCATE(GlobalWeights(1:2,1:5,0:nProcessors-1))
  GlobalWeights=0.
  DO iOut=1,2
    DO iProc=0,nProcessors-1
      GlobalWeights(iOut,1,iProc) = GlobalLoad(iOut,1 ,iProc)+GlobalLoad(iOut,3 ,iProc)
      GlobalWeights(iOut,2,iProc) = SUM(GlobalLoad(iOut,4:8,iProc))+SUM(GlobalLoad(iOut,10:11,iProc))
      IF(.NOT.ALMOSTZERO(GlobalWeights(iOut,2,iProc)))THEN
        GlobalWeights(iOut,3,iProc) =  GlobalWeights(iOut,2,iProc) /  GlobalWeights(iOut,1,iProc)
        nWeights(iOut)=nWeights(iOut)+1
        IF(.NOT.ALMOSTZERO(GlobalLoad(iOut,2,iProc)))&
          GlobalWeights(iOut,4,iProc) = GlobalLoad(iOut,9,iProc) / GlobalLoad(iOut,2,iProc)
        GlobalWeights(iOut,5,iProc) = SUM(GlobalLoad(iOut,4:11,iProc)) / SUM(GlobalLoad(iOut,1:3,iProc))
        PartWeight(1:3,iOut)=PartWeight(1:3,iOut)+GlobalWeights(iOut,3:5,iProc)
      END IF
    END DO ! iProc
  END DO ! iOut
  PartWeight(1,:)=PartWeight(1,:)/nWeights
  PartWeight(2,:)=PartWeight(2,:)/nWeights
  PartWeight(3,:)=PartWeight(3,:)/nWeights
  ParticleMPIWeight=MIN(PartWeight(3,1),PartWeight(1,1))
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
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-PartWeight: ', PartWeight(1,iOut)
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-Comm-Weight:', PartWeight(2,iOut)
    WRITE(ioUnit,'(A20,ES17.5)') 'Mean-TotalWeight:', PartWeight(3,iOut)
    WRITE(ioUnit,'(4(A20))')'Rank','DG','Part','Ratio-Part-DG'
    WRITE(ioUnit,'(A80)')&
          '================================================================================'
    DO iProc=0,nProcessors-1
      WRITE(ioUnit,'(5X,I10,5x,3(3x,ES17.5))') iProc,GlobalWeights(iOut,1:3,iProc)
      WRITE(ioUnit,'(A80)')&
            '--------------------------------------------------------------------------------'
    END DO ! iProc
    WRITE(ioUnit,'(A)') ''
    WRITE(ioUnit,'(12(A20))')'Rank','DG','DGComm','PML','Part-Emission','Part-Tracking','PIC','DSMC','Push','PartComm' &
                                 ,'SplitMerge','UNFP'
    WRITE(ioUnit,'(A120,A120)')&
      '========================================================================================================================',&
      '========================================================================================================================'
    DO iProc=0,nProcessors-1
      WRITE(ioUnit,'(5X,I10,5x,11(3x,ES17.5))') iProc,GlobalLoad(iOut,:,iProc)
      WRITE(ioUnit,'(A120,A120)')&
      '------------------------------------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------------------------------------------------------------------'
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
nTotalParts=0 
nLoadIter  =0


! finialize all arrays
CALL FinalizeBoltzplatz(IsLoadBalance=.TRUE.)
! reallocate
CALL InitBoltzplatz(IsLoadBalance=.TRUE.)

! restart
CALL Restart()

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
  nLocalParts=REAL(nCurrentParts)*0.2 ! parts per stage
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

! particles
IF(nCurrentParts.GT.0)THEN
  LoadSum(4:11)=LoadSum(4:11)+tCurrent(4:11)/nLocalParts
END IF

! last operation
tCurrent=0.
nCurrentParts=0

END SUBROUTINE LoadMeasure


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

END MODULE MOD_LoadBalance
