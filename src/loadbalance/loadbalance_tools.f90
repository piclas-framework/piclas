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
!> Module contains the tools for load_balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_Tools
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI
INTERFACE DomainDecomposition
  MODULE PROCEDURE DomainDecomposition
END INTERFACE

PUBLIC::DomainDecomposition
#endif /*USE_MPI*/

CONTAINS


#if USE_MPI
SUBROUTINE DomainDecomposition()
!===================================================================================================================================
!> Read ElemTime from .h5 container either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: DomainDecompositionWallTime
USE MOD_Restart_Vars         ,ONLY: DoRestart
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_LoadBalance_Vars     ,ONLY: ElemTimeField
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: ElemTime,PerformLoadBalance
#if USE_HDG
USE MOD_LoadBalance_Vars     ,ONLY: ElemHDGSides,TotalHDGSides
USE MOD_Analyze_Vars         ,ONLY: CalcMeshInfo
#endif /*USE_HDG*/
#endif /*USE_LOADBALANCE*/
USE MOD_MPI_Vars             ,ONLY: offsetElemMPI
USE MOD_LoadDistribution     ,ONLY: ApplyWeightDistributionMethod
#ifdef PARTICLES
USE MOD_Particle_TimeStep    ,ONLY: VarTimeStep_InitDistribution
USE MOD_Particle_Vars        ,ONLY: VarTimeStep
USE MOD_LoadBalance_Vars     ,ONLY: ElemTimePart
#endif /*PARTICLES*/
USE MOD_LoadBalance_Vars     ,ONLY: NewImbalance,MaxWeight,MinWeight,ElemGlobalTime,LoadDistri,TargetWeight
USE MOD_IO_HDF5
USE MOD_LoadBalance_Vars  ,ONLY: MPInElemSend,MPIoffsetElemSend,MPInElemRecv,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars  ,ONLY: MPInSideSend,MPIoffsetSideSend,MPInSideRecv,MPIoffsetSideRecv
USE MOD_LoadBalance_Vars  ,ONLY: nElemsOld,offsetElemOld
USE MOD_LoadBalance_Vars  ,ONLY: ElemInfoRank_Shared,ElemInfoRank_Shared_Win
USE MOD_MPI_Shared_Vars   ,ONLY: myComputeNodeRank,MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars,ONLY: ElemInfo_Shared
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_Particle_Mesh_Vars,ONLY: ElemInfo_Shared_Win
USE MOD_LoadDistribution  ,ONLY: WeightDistribution_Equal
USE MOD_MPI_Shared        ,ONLY: BARRIER_AND_SYNC
USE MOD_LoadBalance_Vars  ,ONLY: PartDistri
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!INTEGER,INTENT(IN)  :: single !< read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL          :: ElemTimeExists
REAL,ALLOCATABLE :: WeightSum_proc(:)
REAL             :: StartT,EndT
INTEGER          :: offsetElemSend,offsetElemRecv
INTEGER          :: iProc,iElem,ElemRank,nElemsProc
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("."))')
LBWRITE(UNIT_stdOut,'(A)', ADVANCE='NO')' DOMAIN DECOMPOSITION ...'

GETTIME(StartT)

IF (PerformLoadBalance) THEN
  nElemsOld     = nElems
  offsetElemOld = offsetElem
#if ! defined(PARTICLES)
  CALL CollectiveStop(__STAMP__,'Load balance not implemented for PARTICLES=OFF')
#endif /*defined(PARTICLES)*/
  IF (myComputeNodeRank.EQ.0) &
    ElemInfoRank_Shared  = ElemInfo_Shared(ELEM_RANK,:)
  CALL BARRIER_AND_SYNC(ElemInfoRank_Shared_Win,MPI_COMM_SHARED)
END IF

!simple partition: nGlobalelems/nprocs, do this on proc 0
ALLOCATE(LoadDistri(0:nProcessors-1))
LoadDistri(:)=0.
SDEALLOCATE(PartDistri)
ALLOCATE(PartDistri(0:nProcessors-1))
PartDistri(:)=0

ElemTimeExists=.FALSE.

#ifdef PARTICLES
IF(VarTimeStep%UseDistribution) THEN
! Initialize variable time step distribution (done before domain decomposition to include time step as weighting for load balance)
! Get the time step factor distribution or calculate it from quality factors from the DSMC state (from the MacroRestartFileName)
  CALL VarTimeStep_InitDistribution()
END IF
#endif

IF (DoRestart.OR.PerformLoadBalance) THEN
  !--------------------------------------------------------------------------------------------------------------------------------!
  ! Readin of ElemTime: Read in only by MPIRoot in single mode, only communicate logical ElemTimeExists
  ! because the root performs the distribution of elements (domain decomposition) due to the load distribution scheme
  ALLOCATE(ElemGlobalTime(1:nGlobalElems)) ! Allocate ElemGlobalTime for all MPI ranks
  ElemGlobalTime = 0.

  IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
    CALL ReadElemTime(single=.TRUE.)
  ELSEIF (MPIRoot) THEN
    ! 1) Only MPIRoot does readin of ElemTime during restart
    CALL ReadElemTime(single=.TRUE.)
  END IF

  IF(MPIRoot)THEN
    ! if the elemtime is 0.0, the value must be changed in order to prevent a division by zero
    IF(MAXVAL(ElemGlobalTime).LE.0.0) THEN
      ElemGlobalTime = 1.0
      ElemTimeExists = .FALSE.
    ELSE
      ElemTimeExists = .TRUE.
    END IF
  END IF

  ! 2) Distribute logical information ElemTimeExists
  CALL MPI_BCAST(ElemTimeExists,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)

  ! Distribute the elements according to the selected distribution method
  CALL ApplyWeightDistributionMethod(ElemTimeExists)
ELSE
  ! Simple partition: nGlobalelems/nProcessors
  CALL WeightDistribution_Equal(nProcessors,nGlobalElems,offsetElemMPI)

  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offsetElemMPI,nProcessors+1,MPI_INTEGER,0,MPI_COMM_PICLAS,iERROR)

END IF ! IF(DoRestart.OR.PerformLoadBalance)

! Set local number of elements
nElems=offsetElemMPI(myRank+1)-offsetElemMPI(myRank)
! Set element offset for every processor and write info to log file
offsetElem=offsetElemMPI(myRank)
LOGWRITE(*,'(4(A,I8))')'offsetElem = ',offsetElem,' ,nElems = ', nElems, &
             ' , firstGlobalElemID= ',offsetElem+1,', lastGlobalElemID= ',offsetElem+nElems

IF(PerformLoadBalance)THEN
  ! Only update the mapping of element to rank
  IF (myComputeNodeRank.EQ.0) THEN
    ! Shared array is allocated on compute-node level, compute-node root must update the mapping
    DO iProc = 0,nProcessors-1
      nElemsProc = offsetElemMPI(iProc+1) - offsetElemMPI(iProc)
      ElemInfo_Shared(ELEM_RANK,offsetElemMPI(iProc)+1:offsetElemMPI(iProc)+nElemsProc) = iProc
    END DO ! iProc = 0,nProcessors-1
  END IF ! myComputeNodeRank.EQ.0
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)

  IF(.NOT.UseH5IOLoadBalance)THEN
    ! ------------------------------------------------
    ! Element- and Side-based LB
    ! ------------------------------------------------
    ! Calculate the elements to send
    MPInElemSend      = 0
    MPIoffsetElemSend = 0
    MPInSideSend      = 0
    MPIoffsetSideSend = 0
    ! Loop with the old element over the new elem distribution
    DO iElem = 1,nElemsOld
      ElemRank               = ElemInfo_Shared(ELEM_RANK,offsetElemOld+iElem)+1
      MPInElemSend(ElemRank) = MPInElemSend(ElemRank) + 1
      MPInSideSend(ElemRank) = MPInSideSend(ElemRank) + &
          ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElemOld+iElem) - &
          ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElemOld+iElem)
    END DO

    offsetElemSend = 0
    DO iProc = 2,nProcessors
      MPIoffsetElemSend(iProc) = SUM(MPInElemSend(1:iProc-1))
      MPIoffsetSideSend(iProc) = SUM(MPInSideSend(1:iProc-1))
    END DO

    ! Calculate the elements to send
    MPInElemRecv      = 0
    MPIoffsetElemRecv = 0
    MPInSideRecv      = 0
    MPIoffsetSideRecv = 0
    ! Loop with the new element over the old elem distribution
    DO iElem = 1,nElems
      ElemRank               = ElemInfoRank_Shared(offsetElem+iElem)+1
      MPInElemRecv(ElemRank) = MPInElemRecv(ElemRank) + 1
      MPInSideRecv(ElemRank) = MPInSideRecv(ElemRank) + &
          ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem+iElem) - &
          ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+iElem)
    END DO

    offsetElemRecv = 0
    DO iProc = 2,nProcessors
      MPIoffsetElemRecv(iProc) = SUM(MPInElemRecv(1:iProc-1))
      MPIoffsetSideRecv(iProc) = SUM(MPInSideRecv(1:iProc-1))
    END DO
  END IF ! .NOT.UseH5IOLoadBalance
END IF ! PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)

! Sanity check: local nElems and offset
IF(nElems.LE.0) CALL abort(__STAMP__,' Process did not receive any elements/load! ')

! Read the ElemTime again, but this time with every proc, depending on the domain decomposition in order to write the data
! to the state file (keep ElemTime on restart, if no new ElemTime is calculated during the run or replace with newly measured values
! if LoadBalance is on)
#if USE_LOADBALANCE
IF(ElemTimeExists) CALL ReadElemTime(single=.FALSE.)

#if USE_HDG
! Allocate container for number of master sides for the HDG solver for each element
ALLOCATE(ElemHDGSides(1:nElems))
ElemHDGSides=0
IF(CalcMeshInfo)THEN
  CALL AddToElemData(ElementOut,'ElemHDGSides',IntArray=ElemHDGSides(1:nElems))
END IF ! CalcMeshInfo
TotalHDGSides=0
#endif /*USE_HDG*/

! Set new ElemTime depending on new load distribution
SDEALLOCATE(ElemTime)
ALLOCATE(ElemTime(1:nElems))
ElemTime=0.
CALL AddToElemData(ElementOut,'ElemTime',RealArray=ElemTime(1:nElems))
#endif /*USE_LOADBALANCE*/

#ifdef PARTICLES
ElemTimePart    = 0.
#endif /*PARTICLES*/
ElemTimeField    = 0.

! Calculate new (theoretical) imbalance with offsetElemMPI information
IF(ElemTimeExists.AND.MPIRoot)THEN
  ALLOCATE(WeightSum_proc(0:nProcessors-1))
  DO iProc=0,nProcessors-1
    WeightSum_proc(iProc) = SUM(ElemGlobalTime(1+offsetElemMPI(iProc):offsetElemMPI(iProc+1)))
  END DO
  MaxWeight = MAXVAL(WeightSum_proc)
  MinWeight = MINVAL(WeightSum_proc)
  ! WeightSum (Mesh global value) is already set in BalanceMethod scheme

  ! new computation of current imbalance
  TargetWeight = SUM(WeightSum_proc)/nProcessors
  NewImbalance = (MaxWeight-TargetWeight)/TargetWeight

  IF(TargetWeight.LE.0.0) CALL abort(__STAMP__, &
      ' LoadBalance: TargetWeight = ',RealInfoOpt=TargetWeight)
  LBWRITE(UNIT_stdOut,'(A)') ' Calculated new (theoretical) imbalance with offsetElemMPI information'
  SWRITE(UNIT_stdOut,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A)')&
      ' MinWeight: ', MinWeight, '    MaxWeight: ', MaxWeight, '    TargetWeight: ', TargetWeight,'    NewImbalance: ',&
        NewImbalance,' (theoretical)'
  DEALLOCATE(WeightSum_proc)
ELSE
  LBWRITE(UNIT_stdOut,'(A)', ADVANCE='NO') ' No ElemTime found in restart file ... '
  NewImbalance = -1.
  MaxWeight = -1.
  MinWeight = -1.
END IF
GETTIME(EndT)
DomainDecompositionWallTime=EndT-StartT
CALL DisplayMessageAndTime(DomainDecompositionWallTime, 'DONE!')
END SUBROUTINE DomainDecomposition


SUBROUTINE ReadElemTime(single)
!===================================================================================================================================
!> Read ElemTime from .h5 container either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input             ,ONLY: ReadArray,DatasetExists
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,ElemGlobalTime
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime_tmp
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars       ,ONLY: offsetElemMPIOld
USE MOD_MPI_Vars               ,ONLY: offsetElemMPI
USE MOD_Restart_Vars           ,ONLY: FlushInitialState
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Restart_Vars           ,ONLY: RestartFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)  :: single !< read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL             :: ElemTimeExists
INTEGER             :: iProc
REAL                :: StartT,EndT
REAL,ALLOCATABLE    :: ElemTimeTmp(:)
INTEGER             :: ElemPerProc(0:nProcessors-1)
INTEGER,ALLOCATABLE :: ElemProc(:)
CHARACTER(LEN=32)   :: hilf
!===================================================================================================================================

IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  IF (single) THEN
    DO iProc = 0,nProcessors-1
      ElemPerProc(iProc) = offsetElemMPIOld(iProc+1) - offsetElemMPIOld(iProc)
    END DO
    CALL MPI_GATHERV(ElemTime,nElems,MPI_DOUBLE_PRECISION,ElemGlobalTime,ElemPerProc,offsetElemMPIOld(0:nProcessors-1),MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
  ELSE
    ALLOCATE(ElemTimeTmp(1:nElems))

    ASSOCIATE (&
            counts_send  => INT(MPInElemSend     ) ,&
            disp_send    => INT(MPIoffsetElemSend) ,&
            counts_recv  => INT(MPInElemRecv     ) ,&
            disp_recv    => INT(MPIoffsetElemRecv))
      ! Communicate PartInt over MPI
      CALL MPI_ALLTOALLV(ElemTime,counts_send,disp_send,MPI_DOUBLE_PRECISION,ElemTimeTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_PICLAS,iError)
    END ASSOCIATE

    DEALLOCATE(ElemTime)
    ALLOCATE(ElemTime(1:nElems))
    ElemTime = ElemTimeTmp
    DEALLOCATE(ElemTimeTmp)
  END IF
ELSE
  IF(single)THEN
    hilf='(only MPI root)'
  ELSE
    hilf='with all processors'
  END IF ! single
  IF(MPIRoot)THEN
    WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO') ' | Reading ElemTime '//TRIM(hilf)//' from restart file: ',TRIM(RestartFile),' ...'
    GETTIME(StartT)
  END IF

  ! Read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
  IF (single) THEN
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    CALL DatasetExists(File_ID,'ElemTime',ElemTimeExists)
    IF (ElemTimeExists) &
      CALL ReadArray('ElemTime',2,(/1_IK,INT(nGlobalElems,IK)/),0_IK,2,RealArray=ElemGlobalTime)
    CALL CloseDataFile()
  ELSE
    IF(UseH5IOLoadBalance)THEN
      ! Sanity check: some processors will return ElemTimeExists=F even though it is actually present on the disk
      ! When this happens, the root process and its processors that are on the same node always return ElemTimeExists=T
      ! This points to a corrupt state file (accompanied by SpecID=0 particles within the file)
      ! If the load balance step is performed without h5 I/O in the future, this check can be removed
      CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
      CALL DatasetExists(File_ID,'ElemTime',ElemTimeExists)
      IF(.NOT.ElemTimeExists) CALL abort(__STAMP__,'ElemTime does not exit for some processors in .h5 which indicates a corrupt state file')
      CALL CloseDataFile()

      ! Check if the original ElemTime needs to be communicated to all procs for output to a state file
      IF(FlushInitialState)THEN
        ! Check is allocated and re-allocate
        SDEALLOCATE(ElemTime_tmp)
        ALLOCATE(ElemTime_tmp(1:nElems))
        ElemTime_tmp=0.

        ! Because on some file systems the data of the new state file which was read here might not be completed yet before it accessed
        ! here, it is instead synchronized via mpi scatterv from the root process to all other processes.
        ! This is because HDF5 only guarantees that the data is in the kernel buffer, but not on the disk itself.
        IF(MPIRoot)THEN
          ALLOCATE(ElemProc(0:nProcessors-1))
          DO iProc=0,nProcessors-1
            ElemProc(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
          END DO ! iPRoc
          ! Is this necessary for the root process?
          ElemTime_tmp(1:nElems) = ElemGlobalTime(1:nElems)
        END IF ! MPIRoot

        ! Send from root to all other processes
        CALL MPI_SCATTERV(ElemGlobalTime, ElemProc, offsetElemMPI, MPI_DOUBLE_PRECISION, ElemTime_tmp, nElems, MPI_DOUBLE_PRECISION, 0, MPI_COMM_PICLAS, IERROR)

        ! Deallocate temporary array
        IF(MPIRoot) DEALLOCATE(ElemProc)
      END IF ! FlushInitialState
     ELSE

      ! Normal restart
      SDEALLOCATE(ElemTime_tmp)
      ALLOCATE(ElemTime_tmp(1:nElems))
      ElemTime_tmp  = 0.
      CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
      CALL ReadArray('ElemTime',2,(/1_IK,INT(nElems,IK)/),INT(offsetElem,IK),2,RealArray=ElemTime_tmp)
      CALL CloseDataFile()

     END IF ! UseH5IOLoadBalance
  END IF ! single

  IF(MPIRoot)THEN
    GETTIME(EndT)
    WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE! [',EndT-StartT,'s]'
  END IF
END IF

END SUBROUTINE ReadElemTime

#endif /*USE_MPI*/


END MODULE MOD_LoadBalance_Tools
