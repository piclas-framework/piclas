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

SUBROUTINE Riemann(F,U_L,U_R,nv,GradSide)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! 0
USE MOD_Globals        ,ONLY :abort, vecnorm
USE MOD_Transport_Data ,ONLY: CalcDriftDiffusionCoeff
USE MOD_DSMC_Vars      ,ONLY: BGGas
USE MOD_Particle_Vars  ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar_FV+3,0:0,0:0),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:0,0:0)
REAL,INTENT(IN)                                  :: GradSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar_FV,0:0,0:0)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: normalE, mu_L, mu_R, D_L, D_R, E_L(1:3), E_R(1:3)
INTEGER                                          :: iSpec
!===================================================================================================================================
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    EXIT
  END IF
END DO

E_L = U_L(PP_nVar_FV+1:PP_nVar_FV+3,0,0)
E_R = U_R(PP_nVar_FV+1:PP_nVar_FV+3,0,0)


CALL CalcDriftDiffusionCoeff(VECNORM(E_L),BGGas%NumberDensity(iSpec),mu_L,D_L)
CALL CalcDriftDiffusionCoeff(VECNORM(E_R),BGGas%NumberDensity(iSpec),mu_R,D_R)

normalE = -DOT_PRODUCT(nv(:,0,0),(E_L+E_R)/2.)

IF (normalE.GT.0.) THEN
  F(1,0,0)=normalE*mu_L*U_L(1,0,0) + D_L*GradSide
ELSE
  F(1,0,0)=normalE*mu_R*U_R(1,0,0) + D_R*GradSide
END IF

END SUBROUTINE Riemann


END MODULE MOD_Riemann