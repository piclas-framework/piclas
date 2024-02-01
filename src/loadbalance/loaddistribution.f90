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
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_LoadDistribution
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE ApplyWeightDistributionMethod
  MODULE PROCEDURE ApplyWeightDistributionMethod
END INTERFACE
#endif /*USE_MPI*/

INTERFACE WriteElemTimeStatistics
  MODULE PROCEDURE WriteElemTimeStatistics
END INTERFACE

#if USE_MPI
PUBLIC::ApplyWeightDistributionMethod
PUBLIC::WeightDistribution_Equal
#endif /*USE_MPI*/
PUBLIC::WriteElemTimeStatistics
!===================================================================================================================================

CONTAINS

#if USE_MPI
!===================================================================================================================================
!> Calculate the optimal load partition, subroutine taken from sparta.f90 of HALO
!> Modification for performance on root
!>
!> Algorithm can lead to zero elements per proc!!!!
!===================================================================================================================================
SUBROUTINE SingleStepOptimalPartition(nProcs,OldElems,NewElems,ElemTime)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars,   ONLY:TargetWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: OldElems
REAL,INTENT(IN)                :: ElemTime(1:OldElems)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: NewElems
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: preSum(:)
REAL                           :: LowerBoundary
REAL                           :: UpperBoundary
INTEGER                        :: iElem,iRank
INTEGER                        :: minRank,maxRank,leftOff,lb,ub,mid
! MPI
REAL                           :: LoadSend,opt_split,WeightSplit
INTEGER, ALLOCATABLE           :: split(:),sEND_count(:),recv_count(:)
!===================================================================================================================================

ALLOCATE(PreSum(    1:OldElems) &
        ,send_count(0:nProcs-1) &
        ,split(     0:nProcs-1) &
        ,recv_count(0:nProcs-1))

PreSum(1) = ElemTime(1)
DO iElem = 2,OldElems
  PreSum(iElem) = Presum(iElem-1) + ElemTime(iElem)
END DO ! iElem

LoadSend = PreSum(OldElems)

CALL MPI_EXSCAN(LoadSend,LowerBoundary,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iERROR)
IF(MPIRoot) LowerBoundary = 0.

UpperBoundary = LowerBoundary + PreSum(OldElems)

send_count = 0
split      = 0
recv_count = 0
minRank    = MAX(FLOOR  (lowerBoundary / TargetWeight),0)
maxRank    = MIN(CEILING(upperboundary / TargetWeight),nProcs-1)

leftOff = 1
! retest algorithm with 1 element per proc!!
DO iRank = minRank,maxRank
  lb = leftoff
  ub = OldElems
  opt_split = (iRank+1)*TargetWeight
  IF (iRank*TargetWeight.LT.UpperBoundary) THEN
    DO
      mid = (lb+ub) / 2
      WeightSplit = PreSum(mid) + LowerBoundary
      IF(WeightSplit.EQ.Opt_split) EXIT
      IF(WeightSplit.LT.Opt_split) THEN
        lb = mid
      ELSE
        ub = mid
      END IF
      IF (lb.GE.ub-1) EXIT
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
        IF (opt_split.GT.UpperBoundary) mid = OldElems
      END IF
    END IF

    split(iRank)      = mid
    send_count(iRank) = mid - leftoff + 1
    leftoff           = mid + 1
  END IF
END DO ! iRank=minRank,maxRank

CALL MPI_ALLTOALL(send_count,1,MPI_INTEGER,recv_count,1, MPI_INTEGER,MPI_COMM_PICLAS,iERROR)
NewElems = SUM(recv_count)

DEALLOCATE(PreSum,send_count,recv_count,split)

END SUBROUTINE SingleStepOptimalPartition


!===================================================================================================================================
!> Calculate elem distribution according to selected method
!===================================================================================================================================
SUBROUTINE ApplyWeightDistributionMethod(ElemTimeExists)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_LoadBalance_Vars ,ONLY: ElemGlobalTime,WeightDistributionMethod
USE MOD_Mesh_Vars        ,ONLY: nGlobalElems
USE MOD_MPI_Vars         ,ONLY: offsetElemMPI
USE MOD_ReadInTools      ,ONLY: GETINT,GETREAL
USE MOD_StringTools      ,ONLY: set_formatting,clear_formatting
#ifdef PARTICLES
USE MOD_HDF5_Input       ,ONLY: File_ID,ReadArray,DatasetExists,OpenDataFile,CloseDataFile,ReadAttribute
USE MOD_LoadBalance_Vars ,ONLY: PartDistri,ParticleMPIWeight,nElemsOld
USE MOD_Particle_Vars    ,ONLY: VarTimeStep,PartIntSize
USE MOD_Restart_Vars     ,ONLY: RestartFile
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars ,ONLY: offsetElemMPIOld
USE MOD_Particle_Vars    ,ONLY: PartInt
USE MOD_HDF5_Output_Particles  ,ONLY: FillParticleData
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)             :: ElemTimeExists
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: WeightDistributionMethod_loc
INTEGER                        :: nProcs
INTEGER,ALLOCATABLE            :: PartsInElem(:)
#ifdef PARTICLES
INTEGER,PARAMETER              :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER              :: ELEM_LastPartInd  = 2
INTEGER(KIND=IK),ALLOCATABLE   :: PartIntTmp(:,:)
INTEGER(KIND=IK)               :: locnPart
INTEGER                        :: iElem,iProc
LOGICAL                        :: PartIntExists
REAL                           :: timeWeight(1:nGlobalElems)
LOGICAL                        :: FileVersionExists
REAL                           :: FileVersionHDF5Real
INTEGER                        :: FileVersionHDF5Int
INTEGER(KIND=IK),ALLOCATABLE   :: PartIntGlob(:,:)
INTEGER                        :: ElemPerProc(0:nProcessors-1)
#endif /*PARTICLES*/
!===================================================================================================================================
ALLOCATE(PartsInElem(1:nGlobalElems))
PartsInElem = 0
nProcs      = nProcessors

! Readin of PartInt: Read in only by MPIRoot in single mode because the root performs the distribution of elements (domain decomposition)
! due to the load distribution scheme
#ifdef PARTICLES
IF (MPIRoot) ALLOCATE(PartIntGlob(PartIntSize,1:nGlobalElems))

! Redistribute/read PartInt array
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  DO iProc = 0,nProcessors-1
    ElemPerProc(iProc) = offsetElemMPIOld(iProc+1) - offsetElemMPIOld(iProc)
  END DO
  ! Sanity check
  IF(.NOT.ALLOCATED(PartInt)) CALL abort(__STAMP__,'PartInt is not allocated') ! Missing call to FillParticleData()
  CALL MPI_GATHERV(PartInt,nElemsOld,MPI_DOUBLE_PRECISION,PartIntGlob,ElemPerProc,offsetElemMPIOld(0:nProcessors-1),MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
  PartIntExists = .TRUE.
ELSE
  ! Readin of PartInt: Read in only by MPIRoot in single mode because the root performs the distribution of elements (domain decomposition)
  ! due to the load distribution scheme
  IF (MPIRoot) THEN

    ! Load balancing for particles: read in particle data
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    CALL DatasetExists(File_ID,'PartInt',PartIntExists)
    IF (PartIntExists) THEN

      !check file version
      CALL DatasetExists(File_ID,'File_Version',FileVersionExists,attrib=.TRUE.)
      IF(FileVersionExists)THEN
        CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionHDF5Real)

        ! Depending on the file version, PartInt may have switched dimensions
        IF(FileVersionHDF5Real.LT.2.8)THEN
          ALLOCATE(PartIntTmp(1:nGlobalElems,PartIntSize))
          ! Check integer KIND=8 possibility
          CALL ReadArray('PartInt',2,(/INT(nGlobalElems,IK),INT(PartIntSize,IK)/),0_IK,1,IntegerArray=PartIntTmp)
          ! Switch dimensions
          DO iElem = 1, nGlobalElems
            PartIntGlob(:,iElem) = PartIntTmp(iElem,:)
          END DO ! iElem = FirstElemInd, LastElemInd
          DEALLOCATE(PartIntTmp)
        ELSE
          ! Check integer KIND=8 possibility
          CALL ReadArray('PartInt',2,(/INT(PartIntSize,IK),INT(nGlobalElems,IK)/),0_IK,2,IntegerArray=PartIntGlob)
        END IF ! FileVersionHDF5Real.LT.2.7
      ELSE
        CALL DatasetExists(File_ID,'Piclas_VersionInt',FileVersionExists,attrib=.TRUE.)
        IF (FileVersionExists) THEN
          CALL ReadAttribute(File_ID,'Piclas_VersionInt',1,IntScalar=FileVersionHDF5Int)
        ELSE
          CALL abort(__STAMP__,'Error in ApplyWeightDistributionMethod(): Attribute "Piclas_VersionInt" does not exist!')
        END IF
      ENDIF

    END IF
    CALL CloseDataFile()

  END IF ! MPIRoot
END IF ! PerformLoadBalance

IF(MPIRoot)THEN
  ! Account for variable time stepping
  timeWeight = 1.0
  IF(VarTimeStep%UseDistribution) THEN
    ! If the time step distribution was adapted, the elements should be weighted with the new time step factor
    ! If the distribution is only read-in and not changed, the particle numbers should already fit the time step distribution
    IF(VarTimeStep%AdaptDistribution) THEN
      timeWeight(1:nGlobalelems) = VarTimeStep%ElemWeight(1:nGlobalelems)
    END IF
  END IF

  IF (PartIntExists) THEN
    DO iElem = 1,nGlobalElems
      locnPart           = PartIntGlob(ELEM_LastPartInd,iElem)-PartIntGlob(ELEM_FirstPartInd,iElem)
      PartsInElem(iElem) = INT(locnPart,4) ! switch to KIND=4

      ! Calculate ElemTime according to number of particles in elem if we have no historical information
      IF(.NOT.ElemTimeExists) ElemGlobalTime(iElem) = locnPart*ParticleMPIWeight*timeWeight(iElem) + 1.0
    END DO
  END IF ! PartIntExists

  SDEALLOCATE(PartIntGlob)
END IF ! MPIRoot
#endif /*PARTICLES*/

! Distribute PartsInElem to all procs (Every proc needs to get the information to arrive at the same timedisc)
CALL MPI_BCAST(PartsInElem,nGlobalElems,MPI_INTEGER,0,MPI_COMM_PICLAS,iError)

! Every proc needs to get the information to arrive at the same timedisc
! No historical data and no particles in restart file
IF (.NOT.ElemTimeExists .AND. ALL(PartsInElem(:).EQ.0)) THEN
  ! Check if distribution method is ideal for pure DG
  !> -1 is optimum distri for const. elem-weight
  IF (WeightDistributionMethod.NE.-1) THEN
    CALL set_formatting("red")
    SWRITE(UNIT_stdOut,'(A)') 'WARNING: WeightDistributionMethod other than -1 with neither particles nor ElemTimes!'
    CALL clear_formatting()
  END IF
  WeightDistributionMethod_loc = -1
ELSE
  ! ElemTime always exists during load balance, attempt to respect the user choice
  WeightDistributionMethod_loc = WeightDistributionMethod
END IF ! .NOT.ElemTimeExists .AND. ALL(PartsInElem(:).EQ.0)

! Some methods require all procs, check within the specific routines, if non-MPIRoots immediately return
SELECT CASE(WeightDistributionMethod_loc)

  ! Same as in no-restart: the elements are equally distributed
  CASE(-1)
    CALL WeightDistribution_Equal(nProcs,nGlobalElems,offsetElemMPI)

  ! Distribute to procs using ElemLoads
  CASE(0)
    CALL WeightDistribution_ElemTime(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)

  ! Distribute to procs using ElemLoads, last proc receives the least load
  CASE(1)
    CALL WeightDistribution_ElemTimeLeast(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)

  !------------------------------------------------------------------------------------------------------------------------------!
  CASE(3,4)
    CALL WeightDistribution_SingleStepOptimal(nProcs,nGlobalElems,ElemGlobalTime,ElemTimeExists,PartsInElem,offsetElemMPI)

  ! Minimize max load of all procs based on CASE(-1,0) with iterative smoothing towards last proc
  CASE(5,6)
    CALL WeightDistribution_Iterative(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)

  CASE DEFAULT
    CALL abort(__STAMP__, ' Error in mesh-readin: Invalid load balance distribution for WeightDistributionMethod = ',&
        IntInfoOpt=WeightDistributionMethod_loc)
END SELECT ! WeightDistributionMethod_loc

! Send the load distribution to all other procs
CALL MPI_BCAST(offsetElemMPI,nProcs+1,MPI_INTEGER,0,MPI_COMM_PICLAS,iERROR)

#ifdef PARTICLES
! Set PartDistri for every processor
DO iProc = 0 ,nProcs-1
  DO iElem = offsetElemMPI(iProc)+1,offsetElemMPI(iProc+1)
    PartDistri(iProc) = PartDistri(iProc)+PartsInElem(iElem)
  END DO
END DO ! iProc
#endif /*PARTICLES*/

SDEALLOCATE(PartsInElem)

END SUBROUTINE ApplyWeightDistributionMethod


!===================================================================================================================================
!> Elements are equally distributed
!===================================================================================================================================
SUBROUTINE WeightDistribution_Equal(nProcs,nGlobalElems,offsetElemMPI)
! MODULES
USE MOD_Globals          ,ONLY: MPIRoot
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: nGlobalElems
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nElems,iElem
INTEGER                        :: iProc
!===================================================================================================================================
IF (.NOT.MPIRoot) RETURN

nElems = nGlobalElems/nProcs
iElem  = nGlobalElems-nElems*nProcs
DO iProc = 0,nProcs-1
  offsetElemMPI(iProc) = nElems*iProc+MIN(iProc,iElem)
END DO

! Set element offset for last processor
offsetElemMPI(nProcs) = nGlobalElems

END SUBROUTINE WeightDistribution_Equal


!===================================================================================================================================
!> Distribute to procs using ElemLoads
!===================================================================================================================================
SUBROUTINE WeightDistribution_ElemTime(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)
! MODULES
USE MOD_Globals          ,ONLY: MPIRoot
USE MOD_LoadBalance_Vars ,ONLY: WeightSum
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: nGlobalElems
REAL,INTENT(IN)                :: ElemGlobalTime(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,curiElem
INTEGER                        :: iProc
REAL                           :: CurWeight
!===================================================================================================================================
IF (.NOT.MPIRoot) RETURN

WeightSum = SUM(ElemGlobalTime(:))
CurWeight = 0.

IF (nGlobalElems.EQ.nProcs) THEN
  DO iProc = 0,nProcs-1
    offsetElemMPI(iProc) = iProc
  END DO
ELSE
  curiElem  = 1
  WeightSum = WeightSum/REAL(nProcs)
  DO iProc = 0,nProcs-1
    offsetElemMPI(iProc) = curiElem - 1

    DO iElem = curiElem,nGlobalElems-nProcs+iProc+1
      CurWeight = CurWeight+ElemGlobalTime(iElem)

      IF (CurWeight.GE.WeightSum*(iProc+1)) THEN
        curiElem = iElem + 1
        EXIT
      END IF
    END DO
  END DO
END IF

! Set element offset for last processor
offsetElemMPI(nProcs) = nGlobalElems

END SUBROUTINE WeightDistribution_ElemTime


!===================================================================================================================================
!> Distribute to procs using ElemLoads, last proc recieves least load
!===================================================================================================================================
SUBROUTINE WeightDistribution_ElemTimeLeast(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars ,ONLY: LoadDistri,WeightSum,PerformLoadBalance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: nGlobalElems
REAL,INTENT(IN)                :: ElemGlobalTime(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: FoundDistribution
INTEGER                        :: iDistriIter
INTEGER                        :: iProc
INTEGER                        :: iElem,CuriElem,getElem
REAL                           :: TargetWeight_loc,CurWeight
REAL                           :: MaxLoadDiff,LastLoadDiff,LastProcDiff
REAL                           :: diffLower,diffUpper
INTEGER                        :: ElemDistri(0:nProcs-1)
REAL                           :: LoadDiff(  0:nProcs-1)
!===================================================================================================================================
IF (.NOT.MPIRoot) RETURN

FoundDistribution = .FALSE.
WeightSum         = SUM(ElemGlobalTime(:))
TargetWeight_loc  = WeightSum/REAL(nProcs)
LastProcDiff      = 0.
iDistriIter       = 0
CurWeight         = 0.

LBWRITE(UNIT_stdOut,'(A)') ' Performing iterative search for new load distribution...'

DO WHILE(.NOT.FoundDistribution)
  iDistriIter     = iDistriIter+1
  LBWRITE(UNIT_stdOut,'(A,I4,A,ES15.7)') ' | Iteration ',iDistriIter,' with TargetWeight ',TargetWeight_loc

  TargetWeight_loc = TargetWeight_loc+LastProcDiff/REAL(nProcs)
  curiElem         = 1
  offsetElemMPI    = 0
  offsetElemMPI(nProcs) = nGlobalElems
  LoadDistri       = 0.
  LoadDiff         = 0.

  DO iProc = 0,nProcs-1
    offsetElemMPI(iProc) = curiElem-1
    CurWeight = 0.
    getElem   = 0

    DO iElem = curiElem, nGlobalElems - nProcs +1 + iProc
      CurWeight = CurWeight+ElemGlobalTime(iElem)
      getElem   = getElem+1

      IF ((CurWeight.GT.TargetWeight_loc) .OR. (iElem .EQ. nGlobalElems - nProcs +1 + iProc)) THEN
        diffLower = CurWeight-ElemGlobalTime(iElem)-TargetWeight_loc
        diffUpper = Curweight-TargetWeight_loc

        IF (getElem.GT.1) THEN
          IF (iProc.EQ.nProcs-1) THEN
            LoadDistri(iProc) = CurWeight
            LoadDiff(iProc)   = diffUpper
            curiElem = iElem+1
            EXIT
          ELSE
            IF (ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcs+1+iProc) THEN
              LoadDiff(iProc)   = diffLower
              curiElem = iElem
              LoadDistri(iProc)=CurWeight-ElemGlobalTime(iElem)
              EXIT
            ELSE
              LoadDiff(iProc)   = diffUpper
              curiElem = iElem+1
              LoadDistri(iProc)=CurWeight
              EXIT
            END IF
          END IF
        ELSE ! getElem.GT.1
          LoadDiff(iProc)   = diffUpper
          curiElem = iElem+1
          LoadDistri(iProc)=CurWeight
          EXIT
        END IF ! getElem.GT.1
      END IF
    END DO ! iElem
  END DO ! iProc
  ElemDistri = 0

  DO iProc = 0,nProcs-1
    ElemDistri(iProc) = offsetElemMPI(iProc+1)-offsetElemMPI(iProc)
    ! sanity check
    IF(ElemDistri(iProc).LE.0) CALL abort(__STAMP__,' Process received zero elements during load distribution',iProc)
  END DO ! iPRoc

  ! Determine the remaining load on the last proc
  IF (ElemDistri(nProcs-1).EQ.1) THEN
    LoadDistri(nProcs-1) =     ElemGlobalTime(nGlobalElems)
  ELSE
    LoadDistri(nProcs-1) = SUM(ElemGlobalTime(offsetElemMPI(nProcs-1)+1:nGlobalElems))
  END IF

  LastLoadDiff = LoadDistri(nProcs-1)-TargetWeight_loc
  LoadDiff(nProcs-1) = LastLoadDiff
  MaxLoadDiff  = MAXVAL(LoadDiff(0:nProcs-2))
  LastProcDiff = LastLoadDiff-MaxLoadDiff

  IF (LastProcDiff.LT.0.01*TargetWeight_loc) THEN
    FoundDistribution = .TRUE.
  END IF

  IF(iDistriIter.GT.nProcs) THEN
    SWRITE(UNIT_stdOut,'(A)') ' No valid load distribution throughout the processes found! Alter ParticleMPIWeight!'
    FoundDistribution = .TRUE.
  END IF

  IF (ABS(WeightSum-SUM(LoadDistri)).GT.0.5) CALL Abort(__STAMP__,' Lost Elements and/or Particles during load distribution!')

END DO ! .NOT.FoundDistribution

LBWRITE(UNIT_stdOut,'(A,A17,ES11.4,A,ES11.4,A)') ' Accepted distribution','    TargetWeight: ',TargetWeight_loc,'    (last proc: ',LastProcDiff,')'

END SUBROUTINE WeightDistribution_ElemTimeLeast


!===================================================================================================================================
!> Elements are equally distributed
!===================================================================================================================================
SUBROUTINE WeightDistribution_SingleStepOptimal(nProcs,nGlobalElems,ElemGlobalTime,ElemTimeExists,PartsInElem,offsetElemMPI)
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars ,ONLY: WeightDistributionMethod
USE MOD_LoadBalance_Vars ,ONLY: LoadDistri,ParticleMPIWeight,WeightSum
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: nGlobalElems
REAL,INTENT(IN)                :: ElemGlobalTime(:)
LOGICAL,INTENT(IN)             :: ElemTimeExists
INTEGER,INTENT(IN)             :: PartsInElem(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,MyElems,NewElems,CuriElem,getElem
INTEGER                        :: iProc,jProc
REAL                           :: LastProcDiff
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: ErrorCode
REAL                           :: TargetWeight_loc,CurWeight
REAL                           :: diffLower,diffUpper
INTEGER                        :: ElemDistri(0:nProcs-1)
REAL                           :: LoadDiff(  0:nProcs-1)
!===================================================================================================================================
! Distribute ElemGlobalTime to all procs
CALL MPI_BCAST(ElemGlobalTime,nGlobalElems,MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)

! Do Rebalance
WeightSum = 0.
CurWeight = 0.

DO iElem = 1, nGlobalElems
  WeightSum = WeightSum + ElemGlobalTime(iElem)
END DO

curiElem      = 1
offsetElemMPI = 0
offsetElemMPI(nProcs) = nGlobalElems

TargetWeight_loc = WeightSum/REAL(nProcs)
SWRITE(UNIT_stdOut,'(A,ES16.7)') 'TargetWeight', TargetWeight_loc

IF (WeightDistributionMethod.EQ.3) THEN
  LastProcDiff = 0.
  LoadDistri   = 0.
  LoadDiff     = 0.
END IF

DO iProc = 0,nProcs-1
  offsetElemMPI(iProc) = curiElem-1

  IF (WeightDistributionMethod.EQ.3) THEN
    CurWeight = 0.
    getElem   = 0
  END IF

ElemLoop: DO iElem = curiElem,nGlobalElems-nProcs+iProc+1
    CurWeight = CurWeight+ElemGlobalTime(iElem)

    SELECT CASE(WeightDistributionMethod)
      CASE(3)
        getElem = getElem+1

        IF (CurWeight.GT.TargetWeight_loc .OR. iElem.EQ.nGlobalElems-nProcs +iProc+1) THEN
          diffLower = CurWeight-ElemGlobalTime(iElem)-TargetWeight_loc
          diffUpper = Curweight-TargetWeight_loc

          IF (getElem.GT.1) THEN
            IF (iProc.EQ.nProcs-1) THEN
              LoadDistri(iProc) = CurWeight
              LoadDiff(iProc)   = diffUpper
              curiElem          = iElem+1
              EXIT ElemLoop
            ELSE
              IF (ABS(diffLower).LT.ABS(diffUpper) .AND. iElem.LT.nGlobalElems-nProcs+iProc+1) THEN
                LoadDiff(iProc)   = diffLower
                curiElem          = iElem
                LoadDistri(iProc) = CurWeight-ElemGlobalTime(iElem)
                EXIT
              ELSE
                LoadDiff(iProc)   = diffUpper
                curiElem          = iElem+1
                LoadDistri(iProc) = CurWeight
                EXIT ElemLoop
              END IF
            END IF
          ELSE
            LoadDiff(iProc)   = diffUpper
            curiElem          = iElem+1
            LoadDistri(iProc) = CurWeight
            EXIT ElemLoop
          END IF
        END IF

      CASE(4)
        IF (CurWeight.GE.TargetWeight_loc*(iProc+1)) THEN
          curiElem = iElem + 1
          EXIT ElemLoop
        END IF

    END SELECT
  END DO ElemLoop
END DO ! iProc

ElemDistri = 0

DO iProc = 0,nProcs-1
  ElemDistri(iProc) = offsetElemMPI(iProc+1)-offsetElemMPI(iProc)
  ! sanity check
  IF (ElemDistri(iProc).LE.0) CALL Abort(__STAMP__,'Process received zero elements during load distribution',iProc)
END DO ! iPRoc

! redistribute element weight
DO iProc = 1,nProcs
  FirstElemInd = offsetElemMPI(MyRank)+1
  LastElemInd  = offsetElemMPI(MyRank+1)
  MyElems      = ElemDistri(MyRank)
  CALL SingleStepOptimalPartition(nProcs,MyElems,NewElems,ElemGlobalTime(FirstElemInd:LastElemInd))
  ElemDistri   = 0

  ErrorCode    = 0
  IF(NewElems.LE.0) ErrorCode = ErrorCode+100

  CALL MPI_ALLGATHER(NewElems,1,MPI_INTEGER,ElemDistri(:),1,MPI_INTEGER,MPI_COMM_PICLAS,iERROR)

  ! calculate proc offset
  offsetElemMPI(0) = 0

  DO jProc = 0,nProcs-1
    offsetElemMPI(jProc+1) = offsetElemMPI(jProc) + ElemDistri(jProc)
  END DO ! jProc=1,nProcs

  IF(offsetElemMPI(nProcs).NE.nGlobalElems) ErrorCode = ErrorCode+10
  IF(SUM(ElemDistri)           .NE.nGlobalElems) ErrorCode = ErrorCode+1
  IF(ErrorCode.NE.0) CALL Abort(__STAMP__,' Error during re-distribution! ErrorCode:', ErrorCode)
END DO ! iProc

! compute load distri
LoadDistri = 0.
DO iProc = 0,nProcs-1
  FirstElemInd = offsetElemMPI(iProc)+1
  LastElemInd  = offsetElemMPI(iProc+1)
  IF (ElemTimeExists) THEN
    SELECT CASE(WeightDistributionMethod)
      CASE(3)
        LoadDistri(iProc) = SUM(ElemGlobalTime(FirstElemInd:LastElemInd))
      CASE(4)
        LoadDistri(iProc) = SUM(PartsInElem(FirstElemInd:LastElemInd))
      END SELECT
  ELSE
    LoadDistri(iProc) = LastElemInd-offsetElemMPI(iProc) + SUM(PartsInElem(FirstElemInd:LastElemInd))*ParticleMPIWeight
  END IF
END DO ! iProc

END SUBROUTINE WeightDistribution_SingleStepOptimal


!===================================================================================================================================
!> Minimize max load of all procs based on CASE(-1,0) with iterative smoothing towards last proc
!===================================================================================================================================
SUBROUTINE WeightDistribution_Iterative(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars ,ONLY: WeightDistributionMethod
USE MOD_LoadBalance_Vars ,ONLY: LoadDistri
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: nGlobalElems
REAL,INTENT(IN)                :: ElemGlobalTime(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: FoundDistribution,exitoptimization,globalshift,identical
INTEGER                        :: optIter,itershiftmax,iDistriItermax
INTEGER                        :: iDistriIter,itershift,imax,numOfCalls,nthMinLoad_Idx,startIdx,iShiftLocal,currentRight
INTEGER                        :: MinLoadIdx,MaxLoadIdx,MinLoadIdx_glob,lastopt
REAL                           :: MinLoadVal,MaxLoadVal,MaxLoadVal_opt,MaxLoadVal_opt0
INTEGER                        :: ElemDistri(        0:nProcs-1)
INTEGER                        :: offsetElemMPI_opt( 0:nProcs)
INTEGER                        :: offsetElemMPI_opt0(0:nProcs)
INTEGER                        :: offsetElemMPI_tmp( 0:nProcs)
#if CODE_ANALYZE
INTEGER                        :: iProc
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
IF (.NOT.MPIRoot) RETURN

! estimation, might be set to lower value...
itershiftmax      = nGlobalElems*nProcs*2
exitoptimization  = .FALSE.
FoundDistribution = .FALSE.
itershift         = 0

! Obtain starting distribution
SELECT CASE(WeightDistributionMethod)
  CASE(5)
    CALL WeightDistribution_Equal(nProcs,nGlobalElems,offsetElemMPI)
  CASE(6)
    CALL WeightDistribution_ElemTime(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI)
END SELECT

!-- calc inital distri
CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
    ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)

#if CODE_ANALYZE
SWRITE(*,*)'******** initial distri:',MaxLoadIdx,MinLoadIdx,MinLoadIdx_glob
DO iProc = 0,nProcs-1
  SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
END DO
#endif /*CODE_ANALYZE*/

!-- check for special cases that cannot be further optimized
!>> proc with maxload has only one element -> no further optimization possible
IF (ElemDistri(MaxLoadIdx).EQ.1) THEN
  FoundDistribution = .TRUE.
  exitoptimization  = .TRUE.
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: Max. load is defined by single element!'
!>> trivial, non-optimizable distri
ELSE IF (nProcs.EQ.1 .OR. nGlobalElems.EQ.nProcs) THEN
  FoundDistribution = .TRUE.
  exitoptimization  = .TRUE.
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: trivial, non-optimizable elem-distribution!'
!>> possible optimization
ELSE IF (MinLoadIdx_glob.LT.MaxLoadIdx) THEN
  !-- global minimum is left of maximum, so all need to be shifted (and afterwards smoothed to the right)
  !-- --> shift all to left and calc resulting distri until "left" has more than "right"
  !-- (must be at somepoint for >1 procs when last elem is not more exp. than all other):
  SWRITE(UNIT_stdOut,'(A)') ' | Shifting all to the left (init)'
  currentRight = nProcs-1
  DO WHILE (MinLoadIdx_glob.LT.MaxLoadIdx)
    !-- shift and calc new distri
    offsetElemMPI(1:currentRight) = offsetElemMPI(1:currentRight)+1
    CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
    !-- again, check for single element as max load
    IF (ElemDistri(MaxLoadIdx).EQ.1) THEN
      FoundDistribution = .TRUE.
      exitoptimization  = .TRUE.
      SWRITE(UNIT_stdOut,'(A)') ' WARNING: Max. load is defined by single element!'
      EXIT
    END IF
    !-- check if last proc has now only one elem left, then it has to be shifted from second last proc and so on...
    IF (ElemDistri(currentRight).EQ.1) THEN
      currentRight = currentRight-1
      IF (currentRight.LT.1) THEN
        SWRITE(UNIT_stdOut,'(A)') ' WARNING: already all elements shifted to left!'
        EXIT
      END IF
    END IF
  END DO
END IF

#if CODE_ANALYZE
SWRITE(*,*) '******** adapted distri:',MaxLoadIdx,MinLoadIdx,MinLoadIdx_glob
DO iProc = 0,nProcs-1
  SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
END DO
#endif /*CODE_ANALYZE*/

MaxLoadVal_opt    = MaxLoadVal
offsetElemMPI_opt = offsetElemMPI
optIter     = 0
lastopt     = -1
iDistriIter = 0
globalshift = .TRUE.
numOfCalls  = 0

!-- loop for "shifts" (possibly multiple) shift(s) to left, i.e. elem(s) from last proc are moved to left)
DO WHILE (globalshift)
  iDistriItermax = (nGlobalElems-nProcs+1)*(nProcs-1)*(itershift+2)  ! certain maximum, might be set to lower value...
  iShiftLocal    = 0

  !-- loop for "smoothing" (moving elems from heavy intervals to light intervals)
  DO WHILE (.NOT.FoundDistribution .AND. MaxLoadIdx.NE.nProcs-1)
    ! Should be check and terminated before arriving here
    IF (MaxLoadIdx.GE.MinLoadIdx) CALL Abort(__STAMP__,'MaxLoadIdx.GE.MinLoadIdx! ')

    iShiftLocal       = iShiftLocal+1
    iDistriIter       = iDistriIter+1
    offsetElemMPI_tmp = offsetElemMPI
    MaxLoadVal_opt0   = HUGE(MaxLoadVal_opt)

    ! when shifting the same config, the shifted loads are initially smoothed towards nth minimum index only
    IF (numOfCalls.EQ.0 .OR. iShiftLocal.GT.1) nthMinLoad_Idx = MinLoadIdx
    startIdx = MaxLoadIdx+1

#if CODE_ANALYZE
    SWRITE(*,*)'numOfCalls,startIdx,nthMinLoad_Idx,MinLoadIdx',numOfCalls,startIdx,nthMinLoad_Idx,MinLoadIdx
#endif /*CODE_ANALYZE*/

    !-- smooth towards respective minimum (but also allow further "left" minima as target and ultimatively take best one)
    IF (startIdx.LE.nthMinLoad_Idx) THEN
      DO imax = startIdx,nthMinLoad_Idx
        offsetElemMPI(startIdx:imax) = offsetElemMPI_tmp(startIdx:imax)-1
        CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
            ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
        IF (MaxLoadVal.LT.MaxLoadVal_opt0) THEN
          MaxLoadVal_opt0    = MaxLoadVal
          offsetElemMPI_opt0 = offsetElemMPI
        END IF
      END DO
      offsetElemMPI = offsetElemMPI_opt0
    END IF

    !-- calc best smoothing distri (or same as before when no smooth applicable)
    CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)

#if CODE_ANALYZE
    SWRITE(*,*) '******** iDistriIter',iDistriIter,MaxLoadIdx,MinLoadIdx
    DO iProc = 0,nProcs-1
      SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
    END DO
#endif /*CODE_ANALYZE*/

    !-- check for special cases that cannot be further optimized
    IF (ElemDistri(MaxLoadIdx).EQ.1) THEN
      FoundDistribution = .TRUE.
      exitoptimization  = .TRUE.
      SWRITE(UNIT_stdOut,'(A)') ' WARNING: Max. load is defined by single element!'
    ELSE IF (iDistriIter.GE.iDistriItermax) THEN
      FoundDistribution = .TRUE.
      exitoptimization  = .TRUE.
      SWRITE(UNIT_stdOut,'(A)') ' WARNING: max iternum reached: iDistriIter'
    ! go to next shift...
    ELSE IF (MaxLoadIdx.GE.MinLoadIdx) THEN
      FoundDistribution = .TRUE.
#if CODE_ANALYZE
      SWRITE(*,*) 'MaxLoadIdx.GE.MinLoadIdx...'
#endif /*CODE_ANALYZE*/
    END IF
    !-- save optimal distri
    IF (MaxLoadVal.LT.MaxLoadVal_opt) THEN
      MaxLoadVal_opt    = MaxLoadVal
      offsetElemMPI_opt = offsetElemMPI
      optIter = iDistriIter
    END IF
  END DO

  FoundDistribution = .FALSE. !caution: this flag is now used for identifying to-be shiftes config

#if CODE_ANALYZE
  SWRITE(*,*) '******** optIter',optIter
#endif /*CODE_ANALYZE*/

  !-- check if shifts are applicable for further optimization (first try to shift optimum, then last iter of smoothing-loop)
  offsetElemMPI_tmp = offsetElemMPI
  !check if last iter of previous iteration needs shift
  offsetElemMPI     = offsetElemMPI_opt
  CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
      ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)

  IF (.NOT.exitoptimization         .AND. &
        MaxLoadIdx.EQ.nProcs-1 .AND. &
        itershift .LT.itershiftmax  .AND. &
        lastopt   .NE.optIter) THEN
    FoundDistribution = .TRUE.
    lastopt = optIter

#if CODE_ANALYZE
    SWRITE(*,*) 'shifting all to the left for opt...'
#endif /*CODE_ANALYZE*/

  !-- check if last iter of previous iteration needs shift
  ELSE
    offsetElemMPI = offsetElemMPI_tmp
    CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
    IF (.NOT.exitoptimization       .AND. &
        MaxLoadIdx.EQ.nProcs-1 .AND. &
        itershift .LT.itershiftmax) THEN
      FoundDistribution = .TRUE.

#if CODE_ANALYZE
      SWRITE(*,*)'shifting all to the left for last iter...'
    ELSE IF (itershift.EQ.itershiftmax) THEN
      SWRITE(*,*)'WARNING: max iternum reached: itershift'
    ELSE IF (.NOT.exitoptimization) THEN
      SWRITE(*,*)'exiting shift-iteration, since neither opt-distri nor last iter ended with max load at last proc...'
#endif /*CODE_ANALYZE*/

    END IF
  END IF

  numOfCalls = 0

  IF (FoundDistribution) THEN
    !-- opt-distri or last iter ended with max load at last proc -> "shift to left" (and afterwards smoothed to the right)
    !-- --> shift all to left and calc resulting distri until last proc has not maxload anymore (or only one elem left):
    itershift    = itershift+1
    currentRight = nProcs-1

    DO WHILE (MaxLoadIdx.EQ.nProcs-1)
      offsetElemMPI(1:currentRight) = offsetElemMPI(1:currentRight)+1
      CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
          ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
      !-- check if last proc has now only one elem left, then it has to be shifted from second last proc and so on...
      IF (ElemDistri(currentRight).EQ.1) THEN
        currentRight = currentRight-1
        IF (currentRight.LT.1) THEN
          SWRITE(UNIT_stdOut,'(A)') ' WARNING: already all elements shifted to left!'
          EXIT
        END IF
      END IF
    END DO

    !-- check if resulting config was already present after shift and if or many times (give numOfCalls=0 if no valid found)
    CALL checkList(nProcs,offsetElemMPI,identical,numOfCalls)
    !-- calc again distri and give position of nths(=numOfCalls) minimum-index for ensuring different shift-iteration
    CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob,numOfCalls,nthMinLoad_Idx)

#if CODE_ANALYZE
    SWRITE(*,*)'numOfCalls:',numOfCalls
#endif /*CODE_ANALYZE*/

    IF (numOfCalls.GT.0) THEN
      globalshift = .TRUE.
    ! set to opt and done...
    ELSE
      offsetElemMPI = offsetElemMPI_opt
      CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
          ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
      globalshift   = .FALSE.

#if CODE_ANALYZE
      SWRITE(*,*) 'no valid shift left...'
#endif /*CODE_ANALYZE*/

    END IF
    FoundDistribution = .FALSE.
  ! set to opt and done... (corresponding reason already printed above)
  ELSE
    offsetElemMPI = offsetElemMPI_opt
    CALL CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
        ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob)
    globalshift = .FALSE.
  END IF

#if CODE_ANALYZE
  SWRITE(*,*) '******** final, itershift',itershift,MaxLoadIdx,MinLoadIdx
  DO iProc = 0,nProcs-1
    SWRITE(*,*)'iProc',iProc,LoadDistri(iProc),ElemDistri(iProc)
  END DO
#endif /*CODE_ANALYZE*/

END DO

! Set element offset for last processor
offsetElemMPI(nProcs) = nGlobalElems

CALL freeList()

END SUBROUTINE WeightDistribution_Iterative


!===================================================================================================================================
!> Calculate Distribution from offsetElemMPI
!===================================================================================================================================
SUBROUTINE CalcDistriFromOffsets(nProcs,nGlobalElems,ElemGlobalTime,offsetElemMPI &
  ,ElemDistri,LoadDistri,MaxLoadIdx,MaxLoadVal,MinLoadIdx,MinLoadVal,MinLoadIdx_glob,nth_opt,nthMinLoad_Idx)
! MODULES
USE MOD_Globals           ,ONLY: abort
USE MOD_Utils             ,ONLY: InsertionSort
#ifdef CODE_ANALYZE
USE MOD_Globals           ,ONLY: MPIRoot
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs,nGlobalElems
REAL,INTENT(IN)                :: ElemGlobalTime(1:nGlobalElems)
INTEGER,INTENT(IN)             :: offsetElemMPI(0:nProcs)
INTEGER,INTENT(INOUT),OPTIONAL :: nth_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: ElemDistri(0:nProcs-1)
REAL,INTENT(OUT)               :: LoadDistri(0:nProcs-1)
REAL,INTENT(OUT)               :: MinLoadVal,MaxLoadVal
INTEGER,INTENT(OUT)            :: MinLoadIdx,MaxLoadIdx,MinLoadIdx_glob
INTEGER,INTENT(OUT),OPTIONAL   :: nthMinLoad_Idx
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc,procList(0:nProcs-1),counter,nth
REAL                           :: MinLoadVal_glob,LoadDistri_sort(0:nProcs-1)
!===================================================================================================================================

IF (PRESENT(nth_opt)) THEN
  nth = nth_opt
ELSE
  nth = 0
END IF

DO iProc = 0,nProcs-1
  procList(iProc)   = iProc
  ElemDistri(iProc) = offsetElemMPI(iProc+1)-offsetElemMPI(iProc)
  LoadDistri(iProc) = SUM(ElemGlobalTime(offsetElemMPI(iProc)+1:offsetElemMPI(iProc+1)))
END DO
LoadDistri_sort = LoadDistri
CALL InsertionSort(LoadDistri_sort,procList,nProcs)

MaxLoadIdx = procList(nProcs-1)
MaxLoadVal = LoadDistri_sort(nProcs-1)

MinLoadIdx_glob = procList(0)
MinLoadVal_glob = LoadDistri_sort(0)

DO iProc = 0,nProcs-1
  IF (procList(iProc).GE.MaxLoadIdx) THEN
    MinLoadIdx = procList(iProc)
    MinLoadVal = LoadDistri_sort(iProc)
    EXIT
  END IF
END DO

#ifdef CODE_ANALYZE
SWRITE(*,*)'maxIdx:',MaxLoadIdx,'minIdx:',MinLoadIdx
#endif /*CODE_ANALYZE*/

IF (nth.GT.1) THEN
  counter        = 1
  nthMinLoad_Idx =-1
  IF (.NOT.PRESENT(nthMinLoad_Idx)) CALL Abort(__STAMP__,'nthMinLoad_Idx not present!')

  DO iProc = MinLoadIdx+1,nProcs-1
    IF (procList(iProc).GT.MaxLoadIdx) THEN
      counter = counter+1
      IF (counter.EQ.nth) THEN
        nthMinLoad_Idx=procList(iProc)
        EXIT
      END IF
    END IF
  END DO

  IF (counter.EQ.1 .OR. nthMinLoad_Idx.EQ.-1) THEN
    nth = 0
#ifdef CODE_ANALYZE
    SWRITE(*,*) 'counter set to 0!'
#endif /*CODE_ANALYZE*/
  END IF
ELSE IF (nth.EQ.1) THEN
  IF (.NOT.PRESENT(nthMinLoad_Idx)) CALL abort(__STAMP__,'nthMinLoad_Idx not present!')

  nthMinLoad_Idx = MinLoadIdx
END IF

IF (PRESENT(nth_opt)) nth_opt = nth

END SUBROUTINE CalcDistriFromOffsets


!===================================================================================================================================
!> Check if resulting config was already present after shift and if or many times (give numOfCalls=0 if no valid found)
!===================================================================================================================================
SUBROUTINE checkList(nProcs,offsetElemMPI,identical,numOfCalls)
! MODULES
USE MOD_LoadBalance_Vars
#ifdef CODE_ANALYZE
USE MOD_Globals          ,ONLY: MPIRoot,UNIT_stdOut
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nProcs
INTEGER,INTENT(IN)             :: offsetElemMPI(0:nProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)            :: identical
INTEGER,INTENT(OUT)            :: numOfCalls
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tdata),POINTER            :: newData,tmpData
!===================================================================================================================================
tmpData   => firstData
newData   => tmpData
identical = .FALSE.

DO WHILE(ASSOCIATED(tmpData))
#ifdef CODE_ANALYZE
  SWRITE(UNIT_stdOut,'(A,I0)') 'stored:',tmpData%offsetElemMPI
#endif /*CODE_ANALYZE*/

  ! first access
  IF (ALL(tmpData%offsetElemMPI.EQ.offsetElemMPI)) THEN
    identical = .TRUE.
    tmpData%numOfCalls = tmpData%numOfCalls+1
    ! point to tmpData for output-Var
    newData   => tmpData
  END IF
  tmpData => tmpData%nextData
END DO

#if CODE_ANALYZE
SWRITE(UNIT_stdOut,'(A,I0)') 'current:  ',offsetElemMPI
SWRITE(UNIT_stdOut,'(A,I0)') 'identical:',identical
#endif /*CODE_ANALYZE*/

! read*
IF (.NOT.identical) THEN
  ALLOCATE(newData)
  ALLOCATE(newData%offsetElemMPI(0:nProcs))
  newData%offsetElemMPI = offsetElemMPI
  newData%numOfCalls    = 1
  ! insert at beginning of list
  IF (.NOT.ASSOCIATED(firstData)) THEN
    firstData => newData
  ELSE
    tmpData   => firstData
    firstData => newData
    firstData%nextData => tmpData
  END IF
END IF

numOfCalls = newData%numOfCalls

END SUBROUTINE checkList


!===================================================================================================================================
!> Deallocate list
!===================================================================================================================================
SUBROUTINE freeList()
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
#endif /*USE_MPI*/


!===================================================================================================================================
!> Write load balance info to ElemTimeStatistics.csv file
!> WriteHeader = T: only write the header line to the file removing old data
!===================================================================================================================================
SUBROUTINE WriteElemTimeStatistics(WriteHeader,time_opt,iter_opt)
! MODULES
USE MOD_LoadBalance_Vars ,ONLY: TargetWeight,nLoadBalanceSteps,CurrentImbalance,MinWeight,MaxWeight,WeightSum
USE MOD_Globals          ,ONLY: MPIRoot,FILEEXISTS,unit_stdout,abort,nProcessors,ProcessMemUsage,nProcessors
USE MOD_Globals_Vars     ,ONLY: SimulationEfficiency,PID,WallTime,InitializationWallTime,ReadMeshWallTime,memory
USE MOD_Globals_Vars     ,ONLY: DomainDecompositionWallTime,CommMeshReadinWallTime
USE MOD_Restart_Vars     ,ONLY: DoRestart
#ifdef PARTICLES
USE MOD_LoadBalance_Vars ,ONLY: ElemTimeField
USE MOD_LoadBalance_Vars ,ONLY: ElemTimePart
USE MOD_Globals          ,ONLY: nGlobalNbrOfParticles
#endif /*PARTICLES*/
#if USE_MPI
USE MOD_MPI_Shared_Vars  ,ONLY: myComputeNodeRank,myLeaderGroupRank
USE MOD_Globals
USE MOD_MPI_Shared_Vars  ,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_SHARED
#if ! (CORE_SPLIT==0)
USE MOD_MPI_Shared_Vars  ,ONLY: NbrOfPhysicalNodes,nLeaderGroupProcs
#endif /*! (CORE_SPLIT==0)*/
#endif /*USE_MPI*/
USE MOD_StringTools      ,ONLY: set_formatting,clear_formatting
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars          ,ONLY: MPIW8TimeMM,MPIW8CountMM
#endif /*defined(MEASURE_MPI_WAIT)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)                      :: WriteHeader
REAL,INTENT(IN),OPTIONAL                :: time_opt
INTEGER(KIND=8),INTENT(IN),OPTIONAL     :: iter_opt
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                                     :: time_loc
CHARACTER(LEN=22),PARAMETER              :: outfile='ElemTimeStatistics.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
#ifdef PARTICLES
REAL                                     :: SumElemTime,ElemTimeFieldPercent,ElemTimePartPercent
INTEGER,PARAMETER                        :: nOutputVar=23
#else
INTEGER,PARAMETER                        :: nOutputVar=18
#endif /*PARTICLES*/
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    'time'                   , &
    'Procs'                  , &
    'MinWeight'              , &
    'MaxWeight'              , &
    'CurrentImbalance'       , &
    'TargetWeight (mean)'    , &
    'nLoadBalanceSteps'      , &
    'WeightSum'              , &
    'SimulationEfficiency'   , &
    'PID'                    , &
    'SimulationWallTime'     , &
    'DomainDecompositionWallTime', &
    'CommMeshReadinWallTime' , &
    'ReadMeshWallTime'       , &
    'InitializationWallTime' , &
    'MemoryUsed'             , &
    'MemoryAvailable'        , &
    'MemoryTotal'              &
#ifdef PARTICLES
  , '#Particles'             , &
    'FieldTime'              , &
    'PartTime'               , &
    'FieldTimePercent'       , &
    'PartTimePercent'          &
#endif /*PARTICLES*/
    /)
CHARACTER(LEN=255)         :: tmpStr(nOutputVar) ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)        :: tmpStr2
CHARACTER(LEN=1),PARAMETER :: delimiter=","
#if USE_MPI
REAL                       :: ProcMemoryUsed    ! Used memory on a single proc
REAL                       :: NodeMemoryUsed    ! Sum of used memory across one compute node
#endif /*USE_MPI*/
REAL                       :: MemUsagePercent
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

! Get process memory info
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal

! Only CN roots communicate available and total memory info (count once per node)
#if USE_MPI
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
IF(nProcessors.GT.1)THEN
  ! Collect data on node roots
  ProcMemoryUsed = memory(1)
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_REDUCE(ProcMemoryUsed , NodeMemoryUsed , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
    memory(1) = NodeMemoryUsed
  ELSE
    CALL MPI_REDUCE(ProcMemoryUsed , 0              , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
  END IF

  ! collect data from node roots on first root node
  IF (myComputeNodeRank.EQ.0) THEN ! only leaders
    IF (myLeaderGroupRank.EQ.0) THEN ! first node leader MUST be MPIRoot
      CALL MPI_REDUCE(MPI_IN_PLACE , memory(1:3) , 3 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_LEADERS_SHARED , IERROR)
    ELSE
      CALL MPI_REDUCE(memory(1:3)       , 0      , 3 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_LEADERS_SHARED , IERROR)
    END IF ! myLeaderGroupRank.EQ.0
  END IF ! myComputeNodeRank.EQ.0
END IF ! nProcessors.EQ.1
#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimeMM  = MPIW8TimeMM + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountMM = MPIW8CountMM + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/
#endif /*USE_MPI*/

! --------------------------------------------------
! Only MPI root outputs the data to file
! --------------------------------------------------
IF(.NOT.MPIRoot)RETURN

! Convert kB to GB
memory(1:3)=memory(1:3)/1048576.
#if ! (CORE_SPLIT==0)
! When core-level splitting is used, it is not clear how many cores are on the same physical compute node.
! Therefore, the values are set to -1.
#if USE_MPI
IF(NbrOfPhysicalNodes.GT.0)THEN
  memory(2:3) = memory(2:3) * REAL(NbrOfPhysicalNodes) / REAL(nLeaderGroupProcs)
ELSE
  memory(2:3) = -1.
END IF ! NbrOfPhysicalNodes.GT.0
#endif /*USE_MPI*/
#endif /*! (CORE_SPLIT==0)*/

! Check for memory leaks
#if USE_MPI
IF(MemoryMonitor)THEN
#endif /*USE_MPI*/
  IF(memory(4).LE.0.)THEN
    memory(4) = memory(1) ! Store total initially used memory
  ELSE
    ! Check if more memory is used than at the beginning AND 95% of the total memory available is reached
    MemUsagePercent = (memory(1)/memory(3))*100.0
    !MemUsagePercent = 99.32
    IF((memory(1).GT.memory(4)).AND.(MemUsagePercent.GT.95.0))THEN
      CALL set_formatting("red")
      SWRITE(UNIT_stdOut,'(A,F6.2,A)') ' WARNING: Memory reaching maximum, RAM is at ',MemUsagePercent,'%'
      CALL clear_formatting()
    END IF
  END IF ! WriteHeader
#if USE_MPI
END IF ! MemoryMonitor
#endif /*USE_MPI*/

! Either create new file or add info to existing file
!> create new file
IF(WriteHeader)THEN
  IF(.NOT.PRESENT(iter_opt))CALL abort(__STAMP__, &
      ' WriteElemTimeStatistics: When creating ElemTimeStatistics.csv (WriteHeader=T) then supply [iter_opt] variable')
  IF(iter_opt.GT.0)                         RETURN ! don't create new file if this is not the first iteration
  IF((DoRestart).AND.(FILEEXISTS(outfile))) RETURN ! don't create new file if this is a restart and the file already exists;
  !                                                ! assume continued simulation and old load balance data is still needed

  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr = ' '

  DO I = 1,nOutputVar
    WRITE(tmpStr(I),'(A)') delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO

  WRITE(formatStr,'(A1)') '('

  DO I = 1,nOutputVar
    ! skip writing "," and the end of the line
    IF (I.EQ.nOutputVar) THEN
      WRITE(formatStr,'(A,A1,I2)')    TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)') TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = ' '                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF

IF(PRESENT(time_opt))THEN
  time_loc = time_opt
ELSE
  time_loc = -1.
END IF
#ifdef PARTICLES
! Calculate elem time proportions for field and particle routines
SumElemTime=ElemTimeField+ElemTimePart
IF(SumElemTime.LE.0.)THEN
  ElemTimeFieldPercent = 0.
  ElemTimePartPercent  = 0.
ELSE
  ElemTimeFieldPercent = 100. * ElemTimeField / SumElemTime
  ElemTimePartPercent  = 100. * ElemTimePart / SumElemTime
END IF ! ElemTimeField+ElemTimePart.LE.0.
#endif /*PARTICLES*/

IF (FILEEXISTS(outfile)) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nOutputVar,CSVFORMAT,')'
  WRITE(tmpStr2,formatStr)&
            " ",time_loc                ,&
      delimiter,REAL(nProcessors)       ,&
      delimiter,MinWeight               ,&
      delimiter,MaxWeight               ,&
      delimiter,CurrentImbalance        ,&
      delimiter,TargetWeight            ,&
      delimiter,REAL(nLoadBalanceSteps) ,&
      delimiter,WeightSum               ,&
      delimiter,SimulationEfficiency    ,&
      delimiter,PID                     ,&
      delimiter,WallTime                ,&
      delimiter,DomainDecompositionWallTime,&
      delimiter,CommMeshReadinWallTime  ,&
      delimiter,ReadMeshWallTime        ,&
      delimiter,InitializationWallTime  ,&
      delimiter,memory(1)               ,&
      delimiter,memory(2)               ,&
      delimiter,memory(3)                &
#ifdef PARTICLES
     ,delimiter,REAL(nGlobalNbrOfParticles(3)),&
      delimiter,ElemTimeField              ,&
      delimiter,ElemTimePart               ,&
      delimiter,ElemTimeFieldPercent       ,&
      delimiter,ElemTimePartPercent
#endif /*PARTICLES*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  WRITE(UNIT_stdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write load balance info!"
END IF

END SUBROUTINE WriteElemTimeStatistics

END MODULE MOD_LoadDistribution
