#include "boltzplatz.h"
MODULE MOD_TimeDisc_Vars
!===================================================================================================================================
! Contains global variables used by the Timedisc modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL      :: TEnd
REAL      :: TAnalyze
REAL      :: dt
REAL      :: CFLScale 
REAL      :: eps_LinearSolver,eps2_LinearSolver,maxIter_LinearSolver,epsTilde_LinearSolver
INTEGER(KIND=8) :: iter
LOGICAL   :: ViscousTimeStep=.FALSE.
LOGICAl   :: TimediscInitIsDone = .FALSE.
#if (PP_TimeDiscMethod==1)
REAL,PARAMETER  :: RK3_a2= 5./9. 
REAL,PARAMETER  :: RK3_a3= 153./128.
REAL,PARAMETER  :: RK3_a(2:3) = (/RK3_a2,RK3_a3/)

REAL,PARAMETER  :: RK3_b1= 1./3. 
REAL,PARAMETER  :: RK3_b2= 15./16. 
REAL,PARAMETER  :: RK3_b3= 8./15.
REAL,PARAMETER  :: RK3_b(1:3) = (/RK3_b1,RK3_b2,RK3_b3/)

REAL,PARAMETER  :: RK3_c2= 1./3. 
REAL,PARAMETER  :: RK3_c3= 0.75
REAL,PARAMETER  :: RK3_c(2:3) = (/RK3_c2,RK3_c3/)
#endif
#if ((PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==5) || (PP_TimeDiscMethod==200) || (PP_TimeDiscMethod==201))
REAL,PARAMETER  :: RK4_a2=  567301805773.0/  1357537059087.0
REAL,PARAMETER  :: RK4_a3= 2404267990393.0/  2016746695238.0
REAL,PARAMETER  :: RK4_a4= 3550918686646.0/  2091501179385.0
REAL,PARAMETER  :: RK4_a5= 1275806237668.0/   842570457699.0
REAL,PARAMETER  :: RK4_a(2:5) = (/RK4_a2,RK4_a3,RK4_a4,RK4_a5/)

REAL,PARAMETER  :: RK4_b1= 1432997174477.0/  9575080441755.0
REAL,PARAMETER  :: RK4_b2= 5161836677717.0/ 13612068292357.0
REAL,PARAMETER  :: RK4_b3= 1720146321549.0/  2090206949498.0
REAL,PARAMETER  :: RK4_b4= 3134564353537.0/  4481467310338.0
REAL,PARAMETER  :: RK4_b5= 2277821191437.0/ 14882151754819.0
REAL,PARAMETER  :: RK4_b(1:5) = (/RK4_b1,RK4_b2,RK4_b3,RK4_b4,RK4_b5 /)

REAL,PARAMETER  :: RK4_c2= 1432997174477.0/  9575080441755.0
REAL,PARAMETER  :: RK4_c3= 2526269341429.0/  6820363962896.0
REAL,PARAMETER  :: RK4_c4= 2006345519317.0/  3224310063776.0
REAL,PARAMETER  :: RK4_c5= 2802321613138.0/  2924317926251.0
REAL,PARAMETER  :: RK4_c(2:5) = (/RK4_c2,RK4_c3,RK4_c4,RK4_c5/)
#endif
!===================================================================================================================================
END MODULE MOD_TimeDisc_Vars
