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
INTERFACE SingleStepOptimalPartition
  MODULE PROCEDURE SingleStepOptimalPartition
END INTERFACE

PUBLIC::SingleStepOptimalPartition
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

#ifdef MPI
SUBROUTINE SingleStepOptimalPartition(OldElems,NewElems,ElemWeight) 
!----------------------------------------------------------------------------------------------------------------------------------!
! calculate the optimal load partiton, subroutine taken from sparta.f90 of HALO
! modification for performance on root
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES
!
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,          ONLY:nGlobalElems
USE MOD_LoadBalance_Vars,   ONLY:WeightSum,TargetWeight
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)                :: OldElems
REAL,INTENT(IN)                   :: ElemWeight(1:OldElems)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)               :: NewElems
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: preSum(:)
REAL                           :: LowerBoundary
REAL                           :: UpperBoundary
INTEGER                        :: iElem, iProc,lowerBound,UpperBound, iRank
INTEGER                        :: minRank, maxRank, leftOff, lb, ub,mid
! MPI-Stuff
REAL                           :: LoadSend, opt_split, WeightSplit
INTEGER, ALLOCATABLE           ::  split(:), sEND_count(:), recv_count(:)
!===================================================================================================================================

ALLOCATE(PreSum(1:OldElems)           &
        ,send_count(0:nProcessors-1)  &
        ,split(0:nProcessors-1)       &
        ,recv_count(0:nProcessors-1)  )

PreSum(1)=ElemWeight(1)
DO iElem=2,OldElems
  PreSum(iElem)=Presum(iElem-1)+ElemWeight(iElem)
END DO ! iElem

LoadSend=PreSum(OldElems)

CALL MPI_EXSCAN(LoadSend, LowerBoundary, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, iERROR)
IF(MPIRoot) LowerBoundary=0.

UpperBoundary=LowerBoundary+PreSum(OldElems)

send_count = 0
split = 0
recv_count = 0
minRank = MAX(FLOOR(lowerBoundary/targetWeight),0)
maxRank = MIN(CEILING(upperboundary / targetWeight),nProcessors-1)

leftOff=1

DO iRank=minRank,maxRank
  lb = leftoff
  ub = OldElems
  opt_split = (iRank+1)*targetWeight
  IF(iRank*targetweight.LT.UpperBoundary)THEN
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
    IF(ABS(WeightSplit - Opt_Split) .GT. ABS(WeightSplit-Opt_Split-ElemWeight(mid)))THEN
      ! return 0 IF the splitter is left of the lower boundary
      mid = mid - 1 
    ELSE
      IF (mid+1 .LE. OldElems) THEN
        IF (ABS(WeightSplit - opt_split) .GT. ABS(WeightSplit - opt_split + ElemWeight(mid+1))) THEN
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
#endif /*MPI*/

END MODULE MOD_LoadDistribution
