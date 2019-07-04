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

MODULE MOD_LoadDistribution
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
INTERFACE ApplyWeightDistributionMethod
  MODULE PROCEDURE ApplyWeightDistributionMethod
END INTERFACE
PUBLIC::ApplyWeightDistributionMethod
#endif /*MPI*/

INTERFACE WriteElemTimeStatistics
  MODULE PROCEDURE WriteElemTimeStatistics
END INTERFACE
PUBLIC::WriteElemTimeStatistics
!===================================================================================================================================

CONTAINS

#ifdef MPI
SUBROUTINE SingleStepOptimalPartition(OldElems,NewElems,ElemTime) 
!----------------------------------------------------------------------------------------------------------------------------------!
! Calculate the optimal load partition, subroutine taken from sparta.f90 of HALO
! Modification for performance on root
!
! Algorithm can lead to zero elements per proc!!!!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES
!
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,   ONLY:TargetWeight
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)                :: OldElems
REAL,INTENT(IN)                   :: ElemTime(1:OldElems)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)               :: NewElems
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: preSum(:)
REAL                           :: LowerBoundary
REAL                           :: UpperBoundary
INTEGER                        :: iElem, iRank
INTEGER                        :: minRank, maxRank, leftOff, lb, ub,mid
! MPI-Stuff
REAL                           :: LoadSend, opt_split, WeightSplit
INTEGER, ALLOCATABLE           ::  split(:), sEND_count(:), recv_count(:)
!===================================================================================================================================

ALLOCATE(PreSum(1:OldElems)           &
        ,send_count(0:nProcessors-1)  &
        ,split(0:nProcessors-1)       &
        ,recv_count(0:nProcessors-1)  )

PreSum(1)=ElemTime(1)
DO iElem=2,OldElems
  PreSum(iElem)=Presum(iElem-1)+ElemTime(iElem)
END DO ! iElem

LoadSend=PreSum(OldElems)

CALL MPI_EXSCAN(LoadSend, LowerBoundary, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, iERROR)
IF(MPIRoot) LowerBoundary=0.

UpperBoundary=LowerBoundary+PreSum(OldElems)

send_count = 0
split = 0
recv_count = 0
minRank = MAX(FLOOR(lowerBoundary/TargetWeight),0)
maxRank = MIN(CEILING(upperboundary / TargetWeight),nProcessors-1)

leftOff=1
! retest algorithm with 1 element per proc!!
DO iRank=minRank,maxRank
  lb = leftoff
  ub = OldElems
  opt_split = (iRank+1)*TargetWeight
  IF(iRank*TargetWeight.LT.UpperBoundary)THEN
    DO
      mid = (lb+ub)/2
      WeightSplit= PreSum(mid)+LowerBoundary
      IF(WeightSplit .EQ. Opt_split) EXIT
      IF(WeightSplit .LT. Opt_split) THEN
        lb = mid
      ELSE
        ub = mid
      END IF
      IF(lb.GE.ub-1) EXIT
      ! EXIT IF a single element was found, need to DO this
      !                    here, to have mid and wsplit set.
    END DO
    IF(ABS(WeightSplit - Opt_Split) .GT. ABS(WeightSplit-Opt_Split-ElemTime(mid)))THEN
      ! return 0 IF the splitter is left of the lower boundary
      mid = mid - 1 
    ELSE
      IF (mid+1 .LE. OldElems) THEN
        IF (ABS(WeightSplit - opt_split) .GT. ABS(WeightSplit - opt_split + ElemTime(mid+1))) THEN
          ! return myElems at most
          mid = mid + 1 
        END IF
      ELSE
        IF (opt_split .GT. UpperBoundary) mid = OldElems
      END IF
    END IF
    split(iRank) = mid
    !IF (r .GT. 0) THEN
    !    IF (split(r) .NE. split(r-1)) THEN
    send_count(iRank) = mid - leftoff + 1
    leftoff = mid + 1
    !   END IF
    ! ELSE
    !   send_count(r) = mid - left_off + 1
    !   left_off = mid + 1
    ! END IF
  END IF
END DO ! iRank=minRank,maxRank

CALL MPI_ALLTOALL(send_count, 1, MPI_INTEGER, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, iERROR)
newElems = SUM(recv_count)

DEALLOCATE(PreSum, send_count, recv_count, split)


END SUBROUTINE SingleStepOptimalPartition


SUBROUTINE ApplyWeightDistributionMethod(ElemTimeExists)
!----------------------------------------------------------------------------------------------------------------------------------!
! Description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_MPI_Vars         ,ONLY: offsetElemMPI
USE MOD_LoadBalance_Vars ,ONLY: ElemGlobalTime
USE MOD_Mesh_Vars        ,ONLY: nGlobalElems,nElems
USE MOD_LoadBalance_Vars ,ONLY: LoadDistri,ParticleMPIWeight,WeightSum
#ifdef PARTICLES
USE MOD_LoadBalance_Vars ,ONLY: PartDistri
USE MOD_HDF5_Input       ,ONLY: File_ID,ReadArray,DatasetExists,OpenDataFile,CloseDataFile
USE MOD_Restart_Vars     ,ONLY: RestartFile
USE MOD_Particle_Vars     ,ONLY: VarTimeStep
#endif /*PARTICLES*/
USE MOD_ReadInTools      ,ONLY: GETINT,GETREAL
!----------------------------------------------------------------------------------------------------------------------------------!
! Insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
LOGICAL,INTENT(IN)             :: ElemTimeExists
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: WeightDistributionMethod
INTEGER                        :: iProc,curiElem,MyElems,jProc,NewElems
REAL                           :: MaxLoadDiff,LastLoadDiff,LoadDiff(0:nProcessors-1)
REAL                           :: LastProcDiff
REAL                           :: MinLoadVal,MaxLoadVal,MaxLoadVal_opt,MaxLoadVal_opt0
INTEGER                        :: ElemDistri(0:nProcessors-1),getElem,MinLoadIdx,MaxLoadIdx,MinLoadIdx_glob,lastopt
INTEGER                        :: itershiftmax, iDistriItermax
INTEGER                        :: iElem,optIter
REAL                           :: diffLower,diffUpper
INTEGER                        :: ErrorCode
INTEGER,ALLOCATABLE            :: PartsInElem(:)
LOGICAL                        :: FoundDistribution,exitoptimization,globalshift,identical
REAL                           :: curWeight
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: offsetElemMPI_opt(0:nProcessors),offsetElemMPI_opt0(0:nProcessors)
INTEGER                        :: offsetElemMPI_tmp(0:nProcessors)
INTEGER                        :: iDistriIter,itershift,imax,numOfCalls,nthMinLoad_Idx,startIdx,iShiftLocal,currentRight
#ifdef PARTICLES
INTEGER(KIND=IK),ALLOCATABLE   :: PartInt(:,:)
INTEGER(KIND=IK)               :: locnPart
LOGICAL                        :: PartIntExists
INTEGER,PARAMETER              :: ELEM_FirstPartInd=1
INTEGER,PARAMETER              :: ELEM_LastPartInd=2
REAL                           :: timeWeight(1:nGlobalElems)
#endif /*PARTICLES*/
REAL                           :: TargetWeight_loc
!===================================================================================================================================
WeightSum = 0.0
CurWeight = 0.0

! Load balancing for particles: read in particle data
#ifdef PARTICLES
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'PartInt',PartIntExists)
IF(PartIntExists)THEN
  ALLOCATE(PartInt(1:nGlobalElems,2))
  PartInt(:,:)=0
  ! Check integer KIND=8 possibility
  CALL ReadArray('PartInt',2,(/INT(nGlobalElems,IK),2_IK/),0_IK,1,IntegerArray=PartInt)
END IF
CALL CloseDataFile() 
#endif /*PARTICLES*/
ALLOCATE(PartsInElem(1:nGlobalElems))
PartsInElem=0

#ifdef PARTICLES
timeWeight = 1.0
IF(VarTimeStep%UseDistribution) THEN
  ! If the time step distribution was adapted, the elements should be weighted with the new time step factor
  ! If the distribution is only read-in and not changed, the particle numbers should already fit the time step distribution
  IF(VarTimeStep%AdaptDistribution) THEN
    timeWeight(1:nGlobalelems) = VarTimeStep%ElemWeight(1:nGlobalelems)
  END IF
END IF

IF (PartIntExists) THEN
  DO iElem = 1, nGlobalElems
    locnPart=PartInt(iElem,ELEM_LastPartInd)-PartInt(iElem,ELEM_FirstPartInd)
    PartsInElem(iElem)=INT(locnPart,4) ! switch to KIND=4
    IF(.NOT.ElemTimeExists) ElemGlobalTime(iElem) = locnPart*ParticleMPIWeight*timeWeight(iElem) + 1.0
  END DO
END IF
#endif /*PARTICLES*/
WeightSum = SUM(ElemGlobalTime(:))
!IF(WeightSum.LE.0.0) CALL abort(&
!__STAMP__, &
!' LoadBalance: WeightSum = ',RealInfoOpt=WeightSum)

IF (.NOT.ElemTimeExists .AND. ALL(PartsInElem(:).EQ.0)) THEN
  WeightDistributionMethod = GETINT('WeightDistributionMethod','-1') !-1 is optimum distri for const. elem-weight
  IF (WeightDistributionMethod.NE.-1) THEN
    SWRITE(*,*) 'WARNING: WeightDistributionMethod.NE.-1 with neither particles nor ElemTimes!'
  END IF
ELSE
  WeightDistributionMethod = GETINT('WeightDistributionMethod','1')
END IF

SELECT CASE(WeightDistributionMethod)
CASE(-1) ! same as in no-restart: the elements are equally distributed
  IF(MPIRoot)THEN
    nElems=nGlobalElems/nProcessors
    iElem=nGlobalElems-nElems*nProcessors
    DO iProc=0,nProcessors-1
      offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
    END DO
    offsetElemMPI(nProcessors)=nGlobalElems
  END IF
  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(0) ! old scheme
  IF(MPIRoot)THEN
    IF(nGlobalElems.EQ.nProcessors) THEN
      DO iProc=0, nProcessors-1
        offsetElemMPI(iProc) = iProc
      END DO
    ELSE
      curiElem = 1
      WeightSum=WeightSum/REAL(nProcessors)
      DO iProc=0, nProcessors-1
        offsetElemMPI(iProc)=curiElem - 1
        DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc
          CurWeight=CurWeight+ElemGlobalTime(iElem)
          IF (CurWeight.GE.WeightSum*(iProc+1)) THEN
            curiElem = iElem + 1
            EXIT
          END IF
        END DO
      END DO
    END IF
    offsetElemMPI(nProcessors)=nGlobalElems
  END IF
  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(1)
  ! 1: last Proc receives the least load
  IF(MPIRoot)THEN
    FoundDistribution=.FALSE.
    TargetWeight_loc=WeightSum/REAL(nProcessors)
    LastProcDiff=0.
    iDistriIter=0
    DO WHILE(.NOT.FoundDistribution)
      iDistriIter=iDistriIter+1
      SWRITE(*,'(A19,I4,A19,G0)') '... LoadDistriIter ',iDistriIter,' with TargetWeight=',TargetWeight_loc
      SWRITE(*,*)  lastprocdiff
      TargetWeight_loc=TargetWeight_loc+LastProcDiff/REAL(nProcessors)
      curiElem=1
      offSetElemMPI=0
      offsetElemMPI(nProcessors)=nGlobalElems
      LoadDistri=0.
      LoadDiff=0.
      DO iProc=0,nProcessors-1
        offSetElemMPI(iProc)=curiElem-1
        CurWeight=0.
        getElem=0
        DO iElem=curiElem, nGlobalElems - nProcessors +1 + iProc
          CurWeight=CurWeight+ElemGlobalTime(iElem)
          getElem=getElem+1
          IF((CurWeight.GT.TargetWeight_loc) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
            diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight_loc
            diffUpper=Curweight-TargetWeight_loc
            IF(getElem.GT.1)THEN
              IF(iProc.EQ.nProcessors-1)THEN
                LoadDistri(iProc)=CurWeight
                LoadDiff(iProc)=diffUpper
                curiElem=iElem+1
                EXIT
              ELSE
                IF(ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcessors+1+iProc)THEN
                  LoadDiff(iProc)=diffLower
                  curiElem=iElem
                  LoadDistri(iProc)=CurWeight-ElemGlobalTime(iElem)
                  EXIT
                ELSE
                  LoadDiff(iProc)=diffUpper
                  curiElem=iElem+1
                  LoadDistri(iProc)=CurWeight
                  EXIT
                END IF
              END IF
            ELSE
              LoadDiff(iProc)=diffUpper
              curiElem=iElem+1
              LoadDistri(iProc)=CurWeight
              EXIT
            END IF
          END IF
        END DO ! iElem
      END DO ! iProc
      ElemDistri=0
      DO iProc=0,nProcessors-1
        ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
        ! sanity check
        IF(ElemDistri(iProc).LE.0) CALL abort(&
            __STAMP__&
            ,' Process received zero elements during load distribution',iProc)
      END DO ! iPRoc
      ! Determine the remaining load on the last proc
      IF(ElemDistri(nProcessors-1).EQ.1)THEN
        LoadDistri(nProcessors-1)=ElemGlobalTime(nGlobalElems)
      ELSE
        LoadDistri(nProcessors-1)=SUM(ElemGlobalTime(offSetElemMPI(nProcessors-1)+1:nGlobalElems))
      END IF
      LastLoadDiff = LoadDistri(nProcessors-1)-TargetWeight_loc
      LoadDiff(nProcessors-1)=LastLoadDiff
      MaxLoadDiff=MAXVAL(LoadDiff(0:nProcessors-2))
      LastProcDiff=LastLoadDiff-MaxLoadDiff
      IF(LastProcDiff.LT.0.01*TargetWeight_loc)THEN
        FoundDistribution=.TRUE.
      END IF
      IF(iDistriIter.GT.nProcessors) THEN
        SWRITE(UNIT_StdOut,'(A)') &
            'No valid load distribution throughout the processes found! Alter ParticleMPIWeight!'
        FoundDistribution=.TRUE.
      END IF
      IF(ABS(WeightSum-SUM(LoadDistri)).GT.0.5) THEN
        CALL abort(&
            __STAMP__&
            ,' Lost Elements and/or Particles during load distribution!')
      END IF
    END DO  ! .NOT.FoundDistribution
  END IF    ! MPIRoot
  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(2)
  CALL abort(__STAMP__,'WeightDistributionMethod=2 is not working!')
  ! 1: last Proc receives the least load
  ! 2: Root receives the least load
  IF(MPIRoot)THEN
    FoundDistribution=.FALSE.
    TargetWeight_loc=WeightSum/REAL(nProcessors)
    LastProcDiff=0.
    iDistriIter=0
    DO WHILE(.NOT.FoundDistribution)
      iDistriIter=iDistriIter+1
      SWRITE(*,*) 'LoadDistriIter', iDistriIter
      TargetWeight_loc=TargetWeight_loc+LastProcDiff/REAL(nProcessors)
      curiElem=nGlobalElems
      offSetElemMPI=0
      LoadDistri=0.
      LoadDiff=0.
      DO iProc=nProcessors-1,0,-1
        offSetElemMPI(iProc+1)=curiElem
        CurWeight=0.
        getElem=0
        DO iElem=curiElem,  iProc+1,-1
          CurWeight=CurWeight+ElemGlobalTime(iElem)
          getElem=getElem+1
          IF((CurWeight.GT.TargetWeight_loc) .OR. (iElem .EQ. iProc+1))THEN ! take lower and upper is special case
            diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight_loc
            diffUpper=Curweight-TargetWeight_loc
            IF(getElem.GT.1)THEN
              IF(iProc.EQ.0)THEN
                LoadDistri(iProc)=CurWeight
                LoadDiff(iProc)=diffLower
                curiElem=iElem
                EXIT
              ELSE
                IF(ABS(diffLower).GT.ABS(diffUpper) .AND. iElem.GT.iProc+1)THEN
                  ! take upper
                  LoadDiff(iProc)=diffUpper
                  curiElem=iElem+1
                  LoadDistri(iProc)=CurWeight-ElemGlobalTime(iElem)
                  EXIT
                ELSE
                  LoadDiff(iProc)=diffLower
                  curiElem=iElem
                  LoadDistri(iProc)=CurWeight!+ElemGlobalTime(iElem)
                  EXIT
                END IF
              END IF
            ELSE
              LoadDiff(iProc)=diffLower
              curiElem=iElem
              LoadDistri(iProc)=CurWeight
              EXIT
            END IF
          END IF
        END DO ! iElem
      END DO ! iProc
      ElemDistri=0
      DO iProc=0,nProcessors-1
        ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
        ! sanity check
        IF(ElemDistri(iProc).LE.0) CALL abort(&
            __STAMP__&
            ,' Process received zero elements during load distribution',iProc)
      END DO ! iPRoc
      IF(ElemTimeExists)THEN
        LoadDistri(0)=SUM(ElemGlobalTime(1:offSetElemMPI(1)))
      ELSE
        LoadDistri(0)=ElemDistri(0) +&
            SUM(PartsInElem(1:offSetElemMPI(1)))*ParticleMPIWeight
      END IF
      LastLoadDiff = LoadDistri(0)-TargetWeight_loc
      LoadDiff(0)=LastLoadDiff
      MaxLoadDiff=MAXVAL(LoadDiff(1:nProcessors-1))
      LastProcDiff=LastLoadDiff-MaxLoadDiff
      IF(LastProcDiff.GT.0)THEN
        FoundDistribution=.TRUE.
      END IF
      IF(iDistriIter.EQ.nProcessors) CALL abort(&
          __STAMP__&
          ,'No valid load distribution throughout the processes found! Alter ParticleMPIWeight!')
      IF(ABS(WeightSum-SUM(LoadDistri)).GT.0.5) THEN
        WRITE(*,*) WeightSum-SUM(LoadDistri)
        WRITE(*,*) OffSetElemMPI(1)
        WRITE(*,*) ElemDistri
        WRITE(*,*) LoadDistri
        CALL abort(&
            __STAMP__&
            ,' Lost Elements and/or Particles during load distribution!')
      END IF
    END DO
  END IF
  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(3)
  ! 1: last Proc receives the least load
  ! Distribute ElemGlobalTime to all procs
  CALL MPI_BCAST(ElemGlobalTime,nGlobalElems,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)

  ! Do Rebalance
  WeightSum = 0.
  DO iElem = 1, nGlobalElems
    WeightSum = WeightSum + ElemGlobalTime(iElem)
  END DO
  TargetWeight_loc=WeightSum/REAL(nProcessors)
  LastProcDiff=0.
  curiElem=1
  offSetElemMPI=0
  offsetElemMPI(nProcessors)=nGlobalElems
  LoadDistri=0.
  LoadDiff=0.
  DO iProc=0,nProcessors-1
    offSetElemMPI(iProc)=curiElem-1
    CurWeight=0.
    getElem=0
    DO iElem=curiElem, nGlobalElems - nProcessors +1 + iProc
      CurWeight=CurWeight+ElemGlobalTime(iElem)
      getElem=getElem+1
      IF((CurWeight.GT.TargetWeight_loc) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
        diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight_loc
        diffUpper=Curweight-TargetWeight_loc
        IF(getElem.GT.1)THEN
          IF(iProc.EQ.nProcessors-1)THEN
            LoadDistri(iProc)=CurWeight
            LoadDiff(iProc)=diffUpper
            curiElem=iElem+1
            EXIT
          ELSE
            IF(ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcessors+1+iProc)THEN
              LoadDiff(iProc)=diffLower
              curiElem=iElem
              LoadDistri(iProc)=CurWeight-ElemGlobalTime(iElem)
              EXIT
            ELSE
              LoadDiff(iProc)=diffUpper
              curiElem=iElem+1
              LoadDistri(iProc)=CurWeight
              EXIT
            END IF
          END IF
        ELSE
          LoadDiff(iProc)=diffUpper
          curiElem=iElem+1
          LoadDistri(iProc)=CurWeight
          EXIT
        END IF
      END IF
    END DO ! iElem
  END DO ! iProc
  ElemDistri=0
  DO iProc=0,nProcessors-1
    ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
    ! sanity check
    IF(ElemDistri(iProc).LE.0) CALL abort(&
        __STAMP__&
        ,' Process received zero elements during load distribution',iProc)
  END DO ! iPRoc
  ! redistribute element weight
  DO iProc=1,nProcessors
    FirstElemInd=OffSetElemMPI(MyRank)+1
    LastElemInd =OffSetElemMPI(MyRank+1)
    MyElems=ElemDistri(MyRank)
    CALL SingleStepOptimalPartition(MyElems,NewElems,ElemGlobalTime(FirstElemInd:LastElemInd))
    ElemDistri=0
    CALL MPI_ALLGATHER(NewElems,1,MPI_INTEGER,ElemDistri(:),1,MPI_INTEGER,MPI_COMM_WORLD,iERROR)
    ! calculate proc offset
    OffSetElemMPI(0)=0
    DO jProc=0,nProcessors-1
      OffSetElemMPI(jProc+1) = OffsetElemMPI(jProc) + ElemDistri(jProc)
    END DO ! jProc=0,nProcessors-1
  END DO ! iProc=1,nProcessors
  ! compute load distri
  DO iProc=0,nProcessors-1
    FirstElemInd=OffSetElemMPI(iProc)+1
    LastElemInd =OffSetElemMPI(iProc+1)
    IF(ElemTimeExists)THEN
      LoadDistri(iProc) = SUM(ElemGlobalTime(FirstElemInd:LastElemInd))
    ELSE
      LoadDistri(iProc) = LastElemInd-OffSetElemMPI(iProc) &
          + SUM(PartsInElem(FirstElemInd:LastElemInd))*ParticleMPIWeight
    END IF
    !  SWRITE(*,*) FirstElemInd,LastElemInd,LoadDistri(iProc),SUM(PartsInElem(FirstElemInd:LastElemInd))
  END DO ! iPRoc
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(4)
  ! Distribute ElemGlobalTime to all procs
  CALL MPI_BCAST(ElemGlobalTime,nGlobalElems,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
  ! Do Rebalance
  WeightSum = 0.
  DO iElem = 1, nGlobalElems
    WeightSum = WeightSum + ElemGlobalTime(iElem)
  END DO
  ! predistribute elements
  curiElem = 1
  TargetWeight_loc=WeightSum/REAL(nProcessors)
  SWRITE(*,*) 'TargetWeight', TargetWeight_loc !,ParticleMPIWeight
  offsetElemMPI(nProcessors)=nGlobalElems
  DO iProc=0, nProcessors-1
    offsetElemMPI(iProc)=curiElem - 1 
    DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc  
      CurWeight=CurWeight+ElemGlobalTime(iElem)  
      IF (CurWeight.GE.TargetWeight_loc*(iProc+1)) THEN
        curiElem = iElem + 1 
        EXIT
      END IF
    END DO   
  END DO
  ElemDistri=0
  DO iProc=0,nProcessors-1
    ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
    ! sanity check
    IF(ElemDistri(iProc).LE.0) CALL abort(&
        __STAMP__&
        ,' Process received zero elements during load distribution',iProc)
  END DO ! iPRoc
  ! redistribute element weight
  DO iProc=1,nProcessors
    !SWRITE(*,*) 'distri',iProc
    ErrorCode=0
    FirstElemInd=OffSetElemMPI(MyRank)+1
    LastElemInd =OffSetElemMPI(MyRank+1)
    MyElems=ElemDistri(MyRank)
    CALL SingleStepOptimalPartition(MyElems,NewElems,ElemGlobalTime(FirstElemInd:LastElemInd))
    ElemDistri=0
    IF(NewElems.LE.0) ErrorCode=ErrorCode+100
    CALL MPI_ALLGATHER(NewElems,1,MPI_INTEGER,ElemDistri(:),1,MPI_INTEGER,MPI_COMM_WORLD,iERROR)
    ! calculate proc offset
    OffSetElemMPI(0)=0
    DO jProc=0,nProcessors-1
      OffSetElemMPI(jProc+1) = OffsetElemMPI(jProc) + ElemDistri(jProc)
    END DO ! jProc=1,nProcessors
    IF(OffSetElemMPI(nProcessors).NE.nGlobalElems) ErrorCode=ErrorCode+10
    IF(SUM(ElemDistri).NE.nGlobalElems) ErrorCode=ErrorCode+1
    IF(ErrorCode.NE.0) CALL abort(&
        __STAMP__&
        ,' Error during re-distribution! ErrorCode:', ErrorCode)
  END DO ! jProc=0,nProcessors
  ! compute load distri
  LoadDistri=0.
  DO iProc=0,nProcessors-1
    FirstElemInd=OffSetElemMPI(iProc)+1
    LastElemInd =OffSetElemMPI(iProc+1)
    IF(ElemTimeExists)THEN
      LoadDistri(iProc) = SUM(PartsInElem(FirstElemInd:LastElemInd))
    ELSE
      LoadDistri(iProc) = LastElemInd-OffSetElemMPI(iProc) &
          + SUM(PartsInElem(FirstElemInd:LastElemInd))*ParticleMPIWeight
    END IF
    !  SWRITE(*,*) FirstElemInd,LastElemInd,LoadDistri(iProc),SUM(PartsInElem(FirstElemInd:LastElemInd))
  END DO ! iPRoc
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(5,6)
  ! 5,6: trying to minimize max load of all procs based on CASE(-1,0) with iterative smoothing towards last proc
  IF(MPIRoot)THEN
    itershiftmax=nGlobalElems*nProcessors*2 !estimation, might be set to lower value...
    exitoptimization=.FALSE.
    FoundDistribution=.FALSE.
    nElems=nGlobalElems/nProcessors
    iElem=nGlobalElems-nElems*nProcessors
    itershift=0
    IF (WeightDistributionMethod.EQ.5) THEN !-- init as for CASE(-1)
      DO iProc=0,nProcessors-1
        offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
      END DO
    ELSE ! WeightDistributionMethod.EQ.6    !-- init as for CASE(0)
      IF(nGlobalElems.EQ.nProcessors) THEN
        DO iProc=0, nProcessors-1
          offsetElemMPI(iProc) = iProc
        END DO
      ELSE
        curiElem = 1
        WeightSum=WeightSum/REAL(nProcessors)
        DO iProc=0, nProcessors-1
          offsetElemMPI(iProc)=curiElem - 1
          DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc
            CurWeight=CurWeight+ElemGlobalTime(iElem)
            IF (CurWeight.GE.WeightSum*(iProc+1)) THEN
              curiElem = iElem + 1
              EXIT
            END IF
          END DO
        END DO
      END IF
    END IF !WeightDistributionMethod 5 or 6
    offsetElemMPI(nProcessors)=nGlobalElems
    !-- calc inital distri
    CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
#ifdef CODE_ANALYZE
    SWRITE(*,*)'******** initial distri:',MaxLoadIdx,MinLoadIdx,MinLoadIdx_glob
    DO iProc=0,nProcessors-1
      SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
    END DO
#endif /*CODE_ANALYZE*/
    !-- check for special cases that cannot be further optimized
    IF (ElemDistri(MaxLoadIdx).EQ.1) THEN !proc with maxload has only one element -> no further optimization possible
      FoundDistribution=.TRUE.
      exitoptimization=.TRUE.
      SWRITE(*,*)'WARNING: Max. load is defined by single element!'
    ELSE IF (nProcessors.EQ.1 .OR. nGlobalElems.EQ.nProcessors) THEN !trivial, non-optimizable distri
      FoundDistribution=.TRUE.
      exitoptimization=.TRUE.
      SWRITE(*,*)'WARNING: trivial, non-optimizable elem-distribution!'
    ELSE IF (MinLoadIdx_glob.LT.MaxLoadIdx) THEN
      !-- global minimum is left of maximum, so all need to be shifted (and afterwards smoothed to the right)
      !-- --> shift all to left and calc resulting distri until "left" has more than "right"
      !-- (must be at somepoint for >1 procs when last elem is not more exp. than all other):
      SWRITE(*,*)'shifting all to the left (init)'
      currentRight=nProcessors-1
      DO WHILE (MinLoadIdx_glob.LT.MaxLoadIdx)
        !-- shift and calc new distri
        offSetElemMPI(1:currentRight)=offSetElemMPI(1:currentRight)+1
        CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
        !-- again, check for single element as max load
        IF (ElemDistri(MaxLoadIdx).EQ.1) THEN
          FoundDistribution=.TRUE.
          exitoptimization=.TRUE.
          SWRITE(*,*)'WARNING: Max. load is defined by single element!'
          EXIT
        END IF
        !-- check if last proc has now only one elem left, then it has to be shifted from second last proc and so on...
        IF (ElemDistri(currentRight).EQ.1) THEN
          currentRight=currentRight-1
          IF (currentRight.LT.1) THEN
            SWRITE(*,*)'WARNING: already all elements shifted to left!'
            EXIT
          END IF
        END IF
      END DO
    END IF
#ifdef CODE_ANALYZE
    SWRITE(*,*)'******** adapted distri:',MaxLoadIdx,MinLoadIdx,MinLoadIdx_glob
    DO iProc=0,nProcessors-1
      SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
    END DO
#endif /*CODE_ANALYZE*/
    MaxLoadVal_opt=MaxLoadVal
    offsetElemMPI_opt=offsetElemMPI
    optIter=0
    lastopt=-1
    iDistriIter=0
    globalshift=.TRUE.
    numOfCalls=0
    !-- loop for "shifts" (possibly multiple) shift(s) to left, i.e. elem(s) from last proc are moved to left)
    DO WHILE (globalshift)
      iDistriItermax=(nGlobalElems-nProcessors+1)*(nProcessors-1)*(itershift+2)  !certain maximum, might be set to lower value...
      iShiftLocal=0
      !-- loop for "smoothing" (moving elems from heavy intervals to light intervals)
      DO WHILE(.NOT.FoundDistribution .AND. MaxLoadIdx.NE.nProcessors-1)
        IF (MaxLoadIdx.GE.MinLoadIdx) CALL abort(& !should be catched before...
            __STAMP__, &
            'MaxLoadIdx.GE.MinLoadIdx! ')
        iShiftLocal=iShiftLocal+1
        iDistriIter=iDistriIter+1
        offsetElemMPI_tmp=offsetElemMPI
        MaxLoadVal_opt0=HUGE(MaxLoadVal_opt)
        ! when shifting the same config, the shifted loads are initially smoothed towards nth minimum index only
        IF (numOfCalls.EQ.0 .OR. iShiftLocal.GT.1) nthMinLoad_Idx=MinLoadIdx
        startIdx=MaxLoadIdx+1
#ifdef CODE_ANALYZE
        SWRITE(*,*)'numOfCalls,startIdx,nthMinLoad_Idx,MinLoadIdx',numOfCalls,startIdx,nthMinLoad_Idx,MinLoadIdx
#endif /*CODE_ANALYZE*/
        IF (startIdx.LE.nthMinLoad_Idx) THEN
          !-- smooth towards respective minimum (but also allow further "left" minima as target and ultimatively take best one)
          DO imax=startIdx,nthMinLoad_Idx
            offSetElemMPI(startIdx:imax)=offSetElemMPI_tmp(startIdx:imax)-1
            CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
                ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
            IF (MaxLoadVal.LT.MaxLoadVal_opt0) THEN
              MaxLoadVal_opt0=MaxLoadVal
              offsetElemMPI_opt0=offsetElemMPI
            END IF
          END DO
          offsetElemMPI=offsetElemMPI_opt0
        END IF
        !-- calc best smoothing distri (or same as before when no smooth applicable)
        CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
#ifdef CODE_ANALYZE
        SWRITE(*,*)'******** iDistriIter',iDistriIter,MaxLoadIdx,MinLoadIdx
        DO iProc=0,nProcessors-1
          SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
        END DO
#endif /*CODE_ANALYZE*/
        !-- check for special cases that cannot be further optimized
        IF (ElemDistri(MaxLoadIdx).EQ.1) THEN
          FoundDistribution=.TRUE.
          exitoptimization=.TRUE.
          SWRITE(*,*)'WARNING: Max. load is defined by single element!'
        ELSE IF(iDistriIter.GE.iDistriItermax) THEN
          FoundDistribution=.TRUE.
          exitoptimization=.TRUE.
          SWRITE(*,*)'WARNING: max iternum reached: iDistriIter'
        ELSE IF(MaxLoadIdx.GE.MinLoadIdx) THEN !go to next shift...
          FoundDistribution=.TRUE.
#ifdef CODE_ANALYZE
          SWRITE(*,*)'MaxLoadIdx.GE.MinLoadIdx...'
#endif /*CODE_ANALYZE*/
        END IF    
        !-- save optimal distri
        IF (MaxLoadVal.LT.MaxLoadVal_opt) THEN
          MaxLoadVal_opt=MaxLoadVal
          offsetElemMPI_opt=offsetElemMPI
          optIter=iDistriIter
        END IF
      END DO
      FoundDistribution=.FALSE. !caution: this flag is now used for identifying to-be shiftes config
#ifdef CODE_ANALYZE
      SWRITE(*,*)'******** optIter',optIter
#endif /*CODE_ANALYZE*/
      !-- check if shifts are applicable for further optimization (first try to shift optimum, then last iter of smoothing-loop)
      offsetElemMPI_tmp=offsetElemMPI
      !check if last iter of previous iteration needs shift
      offsetElemMPI=offsetElemMPI_opt
      CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
          ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
      IF (.NOT.exitoptimization .AND. &
          MaxLoadIdx.EQ.nProcessors-1 .AND. &
          itershift.LT.itershiftmax   .AND. &
          lastopt.NE.optIter) THEN
        FoundDistribution=.TRUE.
        lastopt=optIter
#ifdef CODE_ANALYZE
        SWRITE(*,*)'shifting all to the left for opt...'
#endif /*CODE_ANALYZE*/
        !check if last iter of previous iteration needs shift
      ELSE
        offsetElemMPI=offsetElemMPI_tmp
        CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
        IF (.NOT.exitoptimization .AND. &
            MaxLoadIdx.EQ.nProcessors-1 .AND. &
            itershift.LT.itershiftmax) THEN
          FoundDistribution=.TRUE.
#ifdef CODE_ANALYZE
          SWRITE(*,*)'shifting all to the left for last iter...'
        ELSE IF (itershift.EQ.itershiftmax) THEN
          SWRITE(*,*)'WARNING: max iternum reached: itershift'
        ELSE IF (.NOT.exitoptimization) THEN
          SWRITE(*,*)'exiting shift-iteration, since neither opt-distri nor last iter ended with max load at last proc...'
#endif /*CODE_ANALYZE*/
        END IF
      END IF
      numOfCalls=0
      IF (FoundDistribution) THEN
        !-- opt-distri or last iter ended with max load at last proc -> "shift to left" (and afterwards smoothed to the right)
        !-- --> shift all to left and calc resulting distri until last proc has not maxload anymore (or only one elem left):
        itershift=itershift+1
        currentRight=nProcessors-1
        DO WHILE (MaxLoadIdx.EQ.nProcessors-1)
          offSetElemMPI(1:currentRight)=offSetElemMPI(1:currentRight)+1
          CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
              ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
          !-- check if last proc has now only one elem left, then it has to be shifted from second last proc and so on...
          IF (ElemDistri(currentRight).EQ.1) THEN
            currentRight=currentRight-1
            IF (currentRight.LT.1) THEN
              SWRITE(*,*)'WARNING: already all elements shifted to left!'
              EXIT
            END IF
          END IF
        END DO
        !-- check if resulting config was already present after shift and if or many times (give numOfCalls=0 if no valid found)
        CALL checkList(offSetElemMPI,identical,numOfCalls)
        !-- calc again distri and give position of nths(=numOfCalls) minimum-index for ensuring different shift-iteration
        CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob,numOfCalls,nthMinLoad_Idx)
#ifdef CODE_ANALYZE
        SWRITE(*,*)'numOfCalls:',numOfCalls
#endif /*CODE_ANALYZE*/
        IF (numOfCalls.GT.0) THEN
          globalshift=.TRUE.
        ELSE !set to opt and done...
          offsetElemMPI=offsetElemMPI_opt
          CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
              ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
          globalshift=.FALSE.
#ifdef CODE_ANALYZE
          SWRITE(*,*)'no valid shift left...'
#endif /*CODE_ANALYZE*/
        END IF
        FoundDistribution=.FALSE.
      ELSE !set to opt and done... (corresponding reason already printed above)
        offsetElemMPI=offsetElemMPI_opt
        CALL CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
        globalshift=.FALSE.
      END IF
#ifdef CODE_ANALYZE
      SWRITE(*,*)'******** final, itershift',itershift,MaxLoadIdx,MinLoadIdx
      DO iProc=0,nProcessors-1
        SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
      END DO
#endif /*CODE_ANALYZE*/
    END DO
    CALL freeList()
  END IF
  ! at the end communicate the found distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
#ifdef CODE_ANALYZE
  SWRITE(*,*) 'done'
#endif /*CODE_ANALYZE*/
  !------------------------------------------------------------------------------------------------------------------------------!
CASE DEFAULT
  CALL abort(&
      __STAMP__&
      , ' Error in mesh-readin: Invalid load balance distribution for WeightDistributionMethod = ',IntInfoOpt=WeightDistributionMethod)
END SELECT ! WeightDistributionMethod

! Set element offset for last processor
offsetElemMPI(nProcessors)=nGlobalElems

#ifdef PARTICLES
! Set PartDistri for every processor
DO iProc=0,nProcessors-1
  DO iElem=offSetElemMPI(iProc)+1,offSetElemMPI(iProc+1)
    PartDistri(iProc)=PartDistri(iProc)+PartsInElem(iElem)
  END DO
END DO ! iProc
#endif /*PARTICLES*/

SDEALLOCATE(PartsInElem)

END SUBROUTINE ApplyWeightDistributionMethod


SUBROUTINE CalcDistriFromOffsets(nProcessors,nGlobalElems,ElemGlobalTime,offSetElemMPI &
  ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob,nth_opt,nthMinLoad_Idx)
!===================================================================================================================================
! Calculate Distribution from offSetElemMPI
!===================================================================================================================================
! MODULES
USE MOD_Globals ,ONLY: abort
USE MOD_Utils   ,ONLY: InsertionSort
#ifdef CODE_ANALYZE
USE MOD_Globals ,ONLY: mpiroot
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: nProcessors,nGlobalElems
REAL,INTENT(IN)              :: ElemGlobalTime(1:nGlobalElems)
INTEGER,INTENT(IN)           :: offsetElemMPI(0:nProcessors)
INTEGER,INTENT(INOUT),OPTIONAL :: nth_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: ElemDistri(0:nProcessors-1)
REAL,INTENT(OUT)             :: LoadDistri(0:nProcessors-1)
REAL,INTENT(OUT)             :: MinLoadVal,MaxLoadVal
INTEGER,INTENT(OUT)          :: MinLoadIdx,MaxLoadIdx,MinLoadIdx_glob
INTEGER,INTENT(OUT),OPTIONAL :: nthMinLoad_Idx
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iProc,procList(0:nProcessors-1),counter,nth
REAL                         :: MinLoadVal_glob,LoadDistri_sort(0:nProcessors-1)
!===================================================================================================================================

IF (PRESENT(nth_opt)) THEN
  nth=nth_opt
ELSE
  nth=0
END IF

DO iProc=0,nProcessors-1
  procList(iProc)=iProc
  ElemDistri(iProc)=offSetElemMPI(iProc+1)-offSetElemMPI(iProc)
  LoadDistri(iProc)=SUM(ElemGlobalTime(offSetElemMPI(iProc)+1:offSetElemMPI(iProc+1)))
END DO
LoadDistri_sort=LoadDistri
CALL InsertionSort(LoadDistri_sort,procList,nProcessors)

MaxLoadIdx=procList(nProcessors-1)
MaxLoadVal=LoadDistri_sort(nProcessors-1)

MinLoadIdx_glob=procList(0)
MinLoadVal_glob=LoadDistri_sort(0)

DO iProc=0,nProcessors-1
  IF(procList(iProc).GE.MaxLoadIdx) THEN
    MinLoadIdx=procList(iProc)
    MinLoadVal=LoadDistri_sort(iProc)
    EXIT
  END IF
END DO
#ifdef CODE_ANALYZE
SWRITE(*,*)'maxIdx:',MaxLoadIdx,'minIdx:',MinLoadIdx
#endif /*CODE_ANALYZE*/

IF (nth.GT.1) THEN
  counter=1
  nthMinLoad_Idx=-1
  IF (.NOT.PRESENT(nthMinLoad_Idx)) CALL abort(&
__STAMP__&
,'nthMinLoad_Idx not present!')
  DO iProc=MinLoadIdx+1,nProcessors-1
    IF(procList(iProc).GT.MaxLoadIdx) THEN
      counter=counter+1
      IF (counter.EQ.nth) THEN
        nthMinLoad_Idx=procList(iProc)
        EXIT
      END IF
    END IF
  END DO
  IF (counter.EQ.1 .OR. nthMinLoad_Idx.EQ.-1) THEN
    nth=0
#ifdef CODE_ANALYZE
    SWRITE(*,*)'counter set to 0!'
#endif /*CODE_ANALYZE*/
  END IF
ELSE IF (nth.EQ.1) THEN
  IF (.NOT.PRESENT(nthMinLoad_Idx)) CALL abort(&
__STAMP__&
,'nthMinLoad_Idx not present!')
  nthMinLoad_Idx=MinLoadIdx
END IF
IF (PRESENT(nth_opt)) nth_opt = nth

END SUBROUTINE CalcDistriFromOffsets


SUBROUTINE checkList(offSetElemMPI,identical,numOfCalls)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_LoadBalance_Vars
USE MOD_Globals          ,ONLY: nProcessors
#ifdef CODE_ANALYZE
USE MOD_Globals          ,ONLY: mpiroot
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)::offsetElemMPI(0:nProcessors)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)::identical
INTEGER,INTENT(OUT)::numOfCalls
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tdata), POINTER :: newData, tmpData
!===================================================================================================================================
tmpData => firstData
identical=.FALSE.
DO WHILE(ASSOCIATED(tmpData))
#ifdef CODE_ANALYZE
  SWRITE(*,*)'stored:',tmpData%offSetElemMPI
#endif /*CODE_ANALYZE*/
  IF (ALL(tmpData%offSetElemMPI.EQ.offSetElemMPI)) THEN !first access
    identical=.TRUE.
    tmpData%numOfCalls=tmpData%numOfCalls+1
    !EXIT !no exit just for printing all...
    newData => tmpData !point to tmpData for output-Var
  END IF
  tmpData => tmpData%nextData
END DO
#ifdef CODE_ANALYZE
SWRITE(*,*)'current:',offSetElemMPI
SWRITE(*,*)'identical:',identical
#endif /*CODE_ANALYZE*/
!read*
IF(.NOT.identical) THEN
  ALLOCATE(newData)
  ALLOCATE(newData%offSetElemMPI(0:nProcessors))
  newData%offSetElemMPI=offSetElemMPI
  newData%numOfCalls=1
  !insert at beginning of list
  IF (.NOT. ASSOCIATED(firstData)) then
    firstData => newData
  ELSE
    tmpData => firstData
    firstData => newData
    firstData%nextData => tmpData
  END IF
END IF
numOfCalls=newData%numOfCalls

END SUBROUTINE checkList

SUBROUTINE freeList()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_LoadBalance_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tdata), POINTER :: tmpData
!===================================================================================================================================

DO
  tmpData => firstData
  IF (.NOT.ASSOCIATED(tmpData)) EXIT
  firstData => firstData%nextData
  DEALLOCATE(tmpData)
END DO

END SUBROUTINE freeList


#endif /*MPI*/


!----------------------------------------------------------------------------------------------------------------------------------!
!> Write load balance info to ElemTimeStatistics.csv file
!>
!>   WriteHeader = T: only write the header line to the file removing old data
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteElemTimeStatistics(WriteHeader,time,iter)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: TargetWeight,nLoadBalanceSteps,CurrentImbalance,MinWeight,MaxWeight,WeightSum
USE MOD_Globals          ,ONLY: MPIRoot,FILEEXISTS,unit_stdout
USE MOD_Globals_Vars     ,ONLY: SimulationEfficiency,PID,WallTime,InitializationWallTime
USE MOD_Restart_Vars     ,ONLY: DoRestart
USE MOD_Globals          ,ONLY: abort
USE MOD_Globals          ,ONLY: nProcessors
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
LOGICAL,INTENT(IN)                  :: WriteHeader
REAL,INTENT(IN),OPTIONAL            :: time
INTEGER(KIND=8),INTENT(IN),OPTIONAL :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: time_loc
CHARACTER(LEN=22),PARAMETER              :: outfile='ElemTimeStatistics.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=50)                        :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=12
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    'time', &
    'Procs', &
    'MinWeight', &
    'MaxWeight', &
    'CurrentImbalance', &
    'TargetWeight (mean)', &
    'nLoadBalanceSteps', &
    'WeightSum', &
    'SimulationEfficiency',&
    'PID', &
    'SimulationWallTime',&
    'InitializationWallTime'/)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called mutiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2 
CHARACTER(LEN=1),PARAMETER               :: delimiter="," 
!===================================================================================================================================
IF(.NOT.MPIRoot)RETURN

! Either create new file or add info to existing file
IF(WriteHeader)THEN ! create new file
  IF(.NOT.PRESENT(iter))CALL abort(&
      __STAMP__, &
      ' WriteElemTimeStatistics: When creating ElemTimeStatistics.csv (WriteHeader=T) then supply [iter] variable')
  IF(iter.GT.0)                             RETURN ! don't create new file if this is not the first iteration
  IF((DoRestart).AND.(FILEEXISTS(outfile))) RETURN ! don't create new file if this is a restart and the file already exists;
  !                                                ! assume continued simulation and old load balance data is still needed

  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible relimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit) 
ELSE ! 
  IF(.NOT.PRESENT(time))THEN
    time_loc=-1.
  ELSE
    time_loc=time
  END IF
  IF(FILEEXISTS(outfile))THEN
    OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
    WRITE(formatStr,'(A2,I2,A14)')'(',nOutputVar,CSVFORMAT
    WRITE(tmpStr2,formatStr)&
              " ",time_loc, &
        delimiter,REAL(nProcessors), &
        delimiter,MinWeight, &
        delimiter,MaxWeight, &
        delimiter,CurrentImbalance, &
        delimiter,TargetWeight, &
        delimiter,REAL(nLoadBalanceSteps), &
        delimiter,WeightSum, &
        delimiter,SimulationEfficiency,&
        delimiter,PID,&
        delimiter,WallTime,&
        delimiter,InitializationWallTime
    WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
    CLOSE(ioUnit) 
  ELSE
    SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write load balance info!"
  END IF
END IF
END SUBROUTINE WriteElemTimeStatistics


END MODULE MOD_LoadDistribution
