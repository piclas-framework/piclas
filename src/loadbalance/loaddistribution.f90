#include "boltzplatz.h"

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
USE MOD_LoadBalance_Vars ,ONLY: LoadDistri,ParticleMPIWeight,WeightSum,TargetWeight
#ifdef PARTICLES
USE MOD_LoadBalance_Vars ,ONLY: PartDistri
#endif /*PARTICLES*/
USE MOD_ReadInTools      ,ONLY: GETINT
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
INTEGER                        :: iDistriIter
INTEGER                        :: ElemDistri(0:nProcessors-1),getElem
INTEGER                        :: iElem
LOGICAL                        :: FoundDistribution
REAL                           :: diffLower,diffUpper
INTEGER                        :: ErrorCode
INTEGER,ALLOCATABLE            :: PartsInElem(:)
#ifdef PARTICLES
INTEGER,ALLOCATABLE            :: PartInt(:,:)
INTEGER                        :: locnPart
LOGICAL                        :: PartIntExists
#endif /*PARTICLES*/
REAL                           :: curWeight
INTEGER                        :: FirstElemInd,LastElemInd
!===================================================================================================================================
WeightDistributionMethod = GETINT('WeightDistributionMethod','1')

WeightSum = 0.0
CurWeight = 0.0

! Load balancing for particles: read in particle data
#ifdef PARTICLES
CALL H5LEXISTS_F(File_ID,'PartInt',PartIntExists,iERROR)
IF(PartIntExists)THEN
  ALLOCATE(PartInt(1:nGlobalElems,2))
  PartInt(:,:)=0
  CALL ReadArray('PartInt',2,(/nGlobalElems,2/),0,1,IntegerArray=PartInt)
  CALL CloseDataFile() 
END IF
#endif /*PARTICLES*/
ALLOCATE(PartsInElem(1:nGlobalElems))

DO iElem = 1, nGlobalElems
#ifdef PARTICLES
  locnPart=PartInt(iElem,ELEM_LastPartInd)-PartInt(iElem,ELEM_FirstPartInd)
  PartsInElem(iElem)=locnPart
  IF(.NOT.ElemTimeExists) ElemGlobalTime(iElem) = locnPart*ParticleMPIWeight + 1.0
#else
  PartsInElem(iElem)=0
  IF(.NOT.ElemTimeExists) ElemGlobalTime(iElem) = 1.0
#endif /*PARTICLES*/
  WeightSum = WeightSum + ElemGlobalTime(iElem)
END DO
!IF(WeightSum.LE.0.0) CALL abort(&
!__STAMP__, &
!' LoadBalance: WeightSum = ',RealInfoOpt=WeightSum)









SELECT CASE(WeightDistributionMethod)
CASE(-1) ! same as in no-restart: the elements are equally distributed
  nElems=nGlobalElems/nProcessors
  iElem=nGlobalElems-nElems*nProcessors
  DO iProc=0,nProcessors-1
    offsetElemMPI(iProc)=nElems*iProc+MIN(iProc,iElem)
  END DO
  offsetElemMPI(nProcessors)=nGlobalElems
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(0) ! old scheme
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
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(1)
  ! 1: last Proc receives the least load
  ! 2: Root receives the least load
  IF(MPIRoot)THEN
    FoundDistribution=.FALSE.
    TargetWeight=WeightSum/REAL(nProcessors)
    LastProcDiff=0.
    iDistriIter=0
    DO WHILE(.NOT.FoundDistribution)
      iDistriIter=iDistriIter+1
      SWRITE(*,'(A19,I4,A19,G0)') '... LoadDistriIter ',iDistriIter,' with TargetWeight=',TargetWeight
      SWRITE(*,*)  lastprocdiff
      TargetWeight=TargetWeight+LastProcDiff/REAL(nProcessors)
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
          IF((CurWeight.GT.TargetWeight) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
            diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight
            diffUpper=Curweight-TargetWeight
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
      IF(ElemTimeExists)THEN
        IF(ElemDistri(nProcessors-1).EQ.1)THEN
          LoadDistri(nProcessors-1)=ElemGlobalTime(nGlobalElems)
          LastLoadDiff = LoadDistri(nProcessors-1)-TargetWeight
        ELSE
          LoadDistri(nProcessors-1)=SUM(ElemGlobalTime(offSetElemMPI(nProcessors-1)+1:nGlobalElems))
          LastLoadDiff = LoadDistri(nProcessors-1)-TargetWeight
        END IF
      ELSE
        LoadDistri(nProcessors-1)=ElemDistri(nProcessors-1) +&
            SUM(PartsInElem(offSetElemMPI(nProcessors-1)+1:nGlobalElems))*ParticleMPIWeight
        LastLoadDiff = LoadDistri(nProcessors-1)-TargetWeight
      END IF
      LoadDiff(nProcessors-1)=LastLoadDiff
      MaxLoadDiff=MAXVAL(LoadDiff(0:nProcessors-2))
      LastProcDiff=LastLoadDiff-MaxLoadDiff
      IF(LastProcDiff.LT.0.01*TargetWeight)THEN
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
    END DO
  END IF
  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offSetElemMPI,nProcessors+1, MPI_INTEGER,0,MPI_COMM_WORLD,iERROR)
  !------------------------------------------------------------------------------------------------------------------------------!
CASE(2)
  CALL abort(&
      __STAMP__&
      ,' error in load distritubion. please fix me!')
  ! 1: last Proc receives the least load
  ! 2: Root receives the least load
  IF(MPIRoot)THEN
    FoundDistribution=.FALSE.
    TargetWeight=WeightSum/REAL(nProcessors)
    LastProcDiff=0.
    iDistriIter=0
    DO WHILE(.NOT.FoundDistribution)
      iDistriIter=iDistriIter+1
      SWRITE(*,*) 'LoadDistriIter', iDistriIter
      TargetWeight=TargetWeight+LastProcDiff/REAL(nProcessors)
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
          IF((CurWeight.GT.TargetWeight) .OR. (iElem .EQ. iProc+1))THEN ! take lower and upper is special case
            diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight
            diffUpper=Curweight-TargetWeight
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
      LastLoadDiff = LoadDistri(0)-TargetWeight
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
  TargetWeight=WeightSum/REAL(nProcessors)

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
      IF((CurWeight.GT.TargetWeight) .OR. (iElem .EQ. nGlobalElems - nProcessors +1 + iProc))THEN
        diffLower=CurWeight-ElemGlobalTime(iElem)-TargetWeight
        diffUpper=Curweight-TargetWeight
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
  TargetWeight=WeightSum/REAL(nProcessors)
  SWRITE(*,*) 'TargetWeight', TargetWeight,ParticleMPIWeight
  offsetElemMPI(nProcessors)=nGlobalElems
  DO iProc=0, nProcessors-1
    offsetElemMPI(iProc)=curiElem - 1 
    DO iElem = curiElem, nGlobalElems - nProcessors + 1 + iProc  
      CurWeight=CurWeight+ElemGlobalTime(iElem)  
      IF (CurWeight.GE.TargetWeight*(iProc+1)) THEN
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
USE MOD_Globals_Vars     ,ONLY: SimulationEfficiency,PID,SimulationTime
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
LOGICAL,INTENT(IN)                  :: WriteHeader
REAL,INTENT(IN),OPTIONAL            :: time
INTEGER(KIND=8),INTENT(IN),OPTIONAL :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=22),PARAMETER :: outfile='ElemTimeStatistics.csv'
INTEGER                     :: ioUnit,I
CHARACTER(LEN=50)           :: formatStr
INTEGER,PARAMETER           :: nOutputVar=10
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: 'time', &
                                                                           'MinWeight', &
                                                                           'MaxWeight', &
                                                                           'CurrentImbalance', &
                                                                           'TargetWeight (mean)', &
                                                                           'nLoadBalanceSteps', &
                                                                           'WeightSum', &
                                                                           'SimulationEfficiency',&
                                                                           'PID', &
                                                                           'SimulationTime'/)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called mutiple times at the beginning
!===================================================================================================================================
IF(.NOT.MPIRoot)RETURN

IF(WriteHeader)THEN
  IF(iter.EQ.0)THEN
    OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
    tmpStr=""
    DO I=1,nOutputVar
      WRITE(tmpStr(I),'(A)')' "'//TRIM(StrVarNames(I))//'" '
    END DO
    WRITE(formatStr,'(A1)')'('
    DO I=1,nOutputVar
      IF(I.EQ.nOutputVar)THEN
        WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
      ELSE
        WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
      END IF
    END DO
    WRITE(formatStr,'(A,A1)')TRIM(formatStr),')'
    write(ioUnit,formatStr)tmpStr
    CLOSE(ioUnit) 
  END IF
ELSE
  IF(FILEEXISTS(outfile))THEN
    OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
    WRITE(formatStr,'(A1,I1,A14)')'(',nOutputVar,'(1X,E21.14E3))'
    WRITE(ioUnit,formatStr)(/time, MinWeight, MaxWeight, CurrentImbalance, TargetWeight, REAL(nLoadBalanceSteps), WeightSum, &
        SimulationEfficiency,PID,SimulationTime/)
    CLOSE(ioUnit) 
  ELSE
    SWRITE(UNIT_StdOut,'(A)')"ElemTimeStatistics.csv does not exist. cannot write load balance info!"
  END IF
END IF
END SUBROUTINE WriteElemTimeStatistics


END MODULE MOD_LoadDistribution
