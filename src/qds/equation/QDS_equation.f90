#if USE_QDS_DG
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
INTERFACE QDS_ExactFunc
  MODULE PROCEDURE QDS_ExactFunc
END INTERFACE

INTERFACE QDS_InitEquation
  MODULE PROCEDURE QDS_InitEquation
END INTERFACE

PUBLIC::QDS_ExactFunc
PUBLIC::QDS_InitEquation
!===================================================================================================================================
CONTAINS


SUBROUTINE QDS_InitEquation() 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_ReadInTools,     ONLY: GETINT
USE MOD_QDS_Equation_vars, ONLY:QDSIniExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
QDSIniExactFunc = GETINT('QDSIniExactFunc','0')
END SUBROUTINE QDS_InitEquation


SUBROUTINE QDS_ExactFunc(QDSExactFunction,t,tDeriv,x,resu) 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_QDS_DG_Vars,        ONLY:QDSnVar_macro
USE MOD_Particle_Vars,      ONLY:BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)              
INTEGER,INTENT(IN)              :: QDSExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:QDSnVar_macro)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Resu_t(1:QDSnVar_macro),Resu_tt(1:QDSnVar_macro) ! state in conservative variables
REAL                            :: Cent(3)
REAL                            :: EOS_R
!===================================================================================================================================
Cent=x
SELECT CASE (QDSExactFunction)
#ifdef PARTICLES
CASE(0) ! Particles
  Resu=0.
  !resu(1:3)= x(1:3)!*x(1) 
#endif
CASE(1) ! Constant 
  Resu=1.
  Resu_t=0.
  Resu_tt=0.
CASE(100) ! SOD test case
  EOS_R=208.1101 ! [J/(kg K)] for Argon
  IF(x(1).GT.0.0)THEN ! right side of interface
    Resu(1)   = 0.125!*EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 0.1E5/(0.125*EOS_R)
  ELSE ! left side of interface
    Resu(1)   = 1.!EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 1E5/(EOS_R)
  END IF
  !print*,Resu
  !read*

CASE DEFAULT
  SWRITE(*,*)'Exact function not specified'
END SELECT ! QDSExactFunction



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
END SUBROUTINE QDS_ExactFunc


END MODULE MOD_QDS_Equation
#endif /*USE_QDS_DG*/
