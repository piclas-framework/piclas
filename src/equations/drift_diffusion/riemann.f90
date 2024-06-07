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

SUBROUTINE Riemann(F,U_L,U_R,nv,GradSide,E_L,E_R)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! 0
USE MOD_TimeDisc_Vars, ONLY : dt
USE MOD_Globals,  ONLY :abort, vecnorm
USE MOD_Transport_Data ,ONLY: CalcDriftDiffusionCoeff
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar_FV,0:0,0:0),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:0,0:0)
REAL,INTENT(IN)                                  :: GradSide
REAL,INTENT(IN)                                  :: E_L(3),E_R(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar_FV,0:0,0:0)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: driftVelo, mu_L, mu_R, D_L, D_R

!===================================================================================================================================
CALL CalcDriftDiffusionCoeff(VECNORM(E_L),1.E25,mu_L,D_L)
CALL CalcDriftDiffusionCoeff(VECNORM(E_R),1.E25,mu_R,D_R)
driftVelo=-DOT_PRODUCT(nv(:,0,0),(E_L+E_R)/2.)

! print*, driftVelo
! print*, nv
! read*
! IF (mu_L.GT.0.)print*, 'E', 1e21*E/U_L(1,0,0), 1e21*E/U_R(1,0,0)
! IF (mu_L.GT.0.)print*, 'mu', mu_L*U_L(1,0,0), mu_R*U_R(1,0,0)
! IF (mu_L.GT.0.) print*, 'D', D_L/mu_L, D_R/mu_R
! IF (mu_L.GT.0.)read*
! F(1,0,0)=0.5*((driftVelo+abs(driftVelo))*mu_L*U_L(1,0,0)+(driftVelo-abs(driftVelo))*mu_R*U_R(1,0,0)) &
!                     +D_L*GradSide
IF (driftVelo.GT.0.) THEN
  F(1,0,0)=driftVelo*mu_L*U_L(1,0,0) + D_L*GradSide
ELSE
  F(1,0,0)=driftVelo*mu_R*U_R(1,0,0) + D_R*GradSide
END IF

END SUBROUTINE Riemann


END MODULE MOD_Riemann
