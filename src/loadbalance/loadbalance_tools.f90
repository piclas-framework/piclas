#include "boltzplatz.h"

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
INTERFACE LBStartTime
  MODULE PROCEDURE LBStartTime
END INTERFACE

INTERFACE LBSplitTime
  MODULE PROCEDURE LBSplitTime
END INTERFACE

INTERFACE LBPauseTime
  MODULE PROCEDURE LBPauseTime
END INTERFACE

INTERFACE LBElemSplitTime
  MODULE PROCEDURE LBElemSplitTime
END INTERFACE

INTERFACE LBElemPauseTime
  MODULE PROCEDURE LBElemPauseTime
END INTERFACE

PUBLIC::LBStartTime
PUBLIC::LBSplitTime
PUBLIC::LBPauseTime
PUBLIC::LBElemSplitTime
PUBLIC::LBElemPauseTime

CONTAINS

SUBROUTINE LBStartTime(tLBStart)
!===================================================================================================================================
!> calculates and sets start time for Loadbalance.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: LOCALTIME
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(INOUT)  :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBStart = LOCALTIME() ! LB Time Start
END SUBROUTINE LBStartTime

SUBROUTINE LBSplitTime(LB_index,tLBStart)
!===================================================================================================================================
!> Splits the time and resets LB_start. Adds time to tcurrent(LB_index) for current proc
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: LOCALTIME
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample,tCurrent
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: LB_index
REAL,INTENT(INOUT)  :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(LB_index)=tCurrent(LB_index)+tLBEnd-tLBStart
tLBStart = tLBEnd !LOCALTIME() ! LB Time Start
END SUBROUTINE LBSplitTime

SUBROUTINE LBPauseTime(LB_index,tLBStart)
!===================================================================================================================================
!> calculates end time and adds time to tcurrent(LB_index) for current proc
!> does not reset tLBstart
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: LOCALTIME
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample,tCurrent
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: LB_index
REAL,INTENT(IN)     :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(LB_index)=tCurrent(LB_index)+tLBEnd-tLBStart
END SUBROUTINE LBPauseTime


SUBROUTINE LBElemSplitTime(ElemID,tLBStart)
!===================================================================================================================================
!> Splits the time and resets LB_start. Adds time to Elemtime(ElemID)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: LOCALTIME
USE MOD_LoadBalance_Vars ,ONLY: ElemTime, PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: ElemID
REAL,INTENT(INOUT)  :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
ElemTime(ELemID)=ElemTime(ElemID)+tLBEnd-tLBStart
tLBStart = tLBEnd !LOCALTIME() ! LB Time Start
END SUBROUTINE LBElemSplitTime

SUBROUTINE LBElemPauseTime(ElemID,tLBStart)
!===================================================================================================================================
!> calculates end time and adds time to Elemtime(ElemID)
!> does not reset tLBstart
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals          ,ONLY: LOCALTIME
USE MOD_LoadBalance_Vars ,ONLY: ElemTime, PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: ElemID
REAL,INTENT(IN)     :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
ElemTime(ELemID)=ElemTime(ElemID)+tLBEnd-tLBStart
END SUBROUTINE LBElemPauseTime


END MODULE MOD_LoadBalance_Tools
