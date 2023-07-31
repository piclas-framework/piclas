#include "piclas.h"

MODULE MOD_Riemann
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

PUBLIC::Riemann
!===================================================================================================================================

CONTAINS

SUBROUTINE Riemann(F,U_L,U_R,nv,GradSide,E)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars
USE MOD_TimeDisc_Vars, ONLY : dt
USE MOD_Globals,  ONLY :abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                                  :: GradSide
REAL,INTENT(IN)                                  :: E(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3), driftVelo, mu, D
INTEGER                                          :: Count_1,Count_2
!===================================================================================================================================
! Gauss point i,j
  DO Count_2=0,PP_N
    DO Count_1=0,PP_N
      n_loc(:)=nv(:,Count_1,Count_2)
      mu = 1.
      driftVelo=-mu*DOT_PRODUCT(n_loc,E)
      D = 1.
      F(1,Count_1,Count_2)=0.5*((driftVelo+abs(driftVelo))*U_L(1,Count_1,Count_2)+(driftVelo-abs(driftVelo))*U_R(1,Count_1,Count_2)) &
                          +D*GradSide
    END DO
  END DO
END SUBROUTINE Riemann


END MODULE MOD_Riemann
