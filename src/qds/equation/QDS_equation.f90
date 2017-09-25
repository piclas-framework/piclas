#include "boltzplatz.h"

MODULE MOD_QDS_Equation
!===================================================================================================================================
!> Contains the routines to
!> - determine exact functions for QDS-Ddetermine exact functions for QDS-DG
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ExactFuncQDS
  MODULE PROCEDURE ExactFuncQDS
END INTERFACE

PUBLIC::ExactFuncQDS
!===================================================================================================================================
CONTAINS



SUBROUTINE ExactFuncQDS(ExactFunction,t,tDeriv,x,resu) 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_QDS_DG_Vars,        ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)              
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:QDSnVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Resu_t(1:QDSnVar),Resu_tt(1:QDSnVar) ! state in conservative variables
REAL                            :: Cent(3)
!===================================================================================================================================
Cent=x
SELECT CASE (ExactFunction)
#ifdef PARTICLES
CASE(0) ! Particles
  Resu=0.
  !resu(1:3)= x(1:3)!*x(1) 
#endif
CASE(1) ! Constant 
  Resu=1.
  Resu_t=0.
  Resu_tt=0.

CASE DEFAULT
  SWRITE(*,*)'Exact function not specified'
END SELECT ! ExactFunction



# if (PP_TimeDiscMethod==1)
! For O3 RK, the boundary condition has to be adjusted
! Works only for O3 RK!!
SELECT CASE(tDeriv)
CASE(0)
  ! resu = g(t)
CASE(1)
  ! resu = g(t) + dt/3*g'(t)
  Resu=Resu + dt/3.*Resu_t
CASE(2)
  ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
  Resu=Resu + 0.75*dt*Resu_t+5./16.*dt*dt*Resu_tt
CASE DEFAULT
  ! Stop, works only for 3 Stage O3 LS RK
  CALL abort(&
      __STAMP__&
      ,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE ExactFuncQDS


END MODULE MOD_QDS_Equation
