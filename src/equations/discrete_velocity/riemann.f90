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

SUBROUTINE Riemann(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_DistFunc, ONLY  : MacroValuesFromDistribution
USE MOD_DistFunc, ONLY  : MaxwellDistribution, MaxwellDistributionCons, ShakhovDistribution, ESBGKDistribution
USE MOD_Equation_Vars_FV,ONLY: DVMDim, DVMnVelos, DVMVelos, DVMMethod, DVMBGKModel
USE MOD_TimeDisc_Vars, ONLY : dt
USE MOD_Globals,  ONLY :abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar_FV,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar_FV,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),MacroVal_L(14), MacroVal_R(14), tau_L, tau_R
REAL                                             :: Velo
REAL,DIMENSION(PP_nVar_FV)                       :: fTarget_L, fTarget_R, UTemp_L, UTemp_R
REAL                                             :: gamma_R, gamma_L
INTEGER                                          :: Count_1,Count_2, iVel, jVel, kVel, upos
!===================================================================================================================================
! Gauss point i,j
  DO Count_2=0,PP_N
    DO Count_1=0,PP_N
      n_loc(:)=nv(:,Count_1,Count_2)
      CALL MacroValuesFromDistribution(MacroVal_L,U_L(:,Count_1,Count_2),dt/2.,tau_L,1)
      CALL MacroValuesFromDistribution(MacroVal_R,U_R(:,Count_1,Count_2),dt/2.,tau_R,1)
      SELECT CASE (DVMBGKModel)
        CASE(1)
          CALL ESBGKDistribution(MacroVal_L,fTarget_L)
          CALL ESBGKDistribution(MacroVal_R,fTarget_R)
        CASE(2)
          CALL ShakhovDistribution(MacroVal_L,fTarget_L)
          CALL ShakhovDistribution(MacroVal_R,fTarget_R)
        CASE(3)
          CALL MaxwellDistribution(MacroVal_L,fTarget_L)
          CALL MaxwellDistribution(MacroVal_R,fTarget_R)
        CASE(4)
          CALL MaxwellDistributionCons(MacroVal_L,fTarget_L)
          CALL MaxwellDistributionCons(MacroVal_R,fTarget_R)
        CASE DEFAULT
          CALL abort(__STAMP__,'DVM BGK Model not implemented.',999,999.)
      END SELECT
      IF (dt.EQ.0.) THEN
        UTemp_L = 0.
        UTemp_R = 0.
      ELSE
        SELECT CASE (DVMMethod)
        CASE(1)
          gamma_L = 2.*tau_L*(1.-EXP(-dt/2./tau_L))/dt
          gamma_R = 2.*tau_R*(1.-EXP(-dt/2./tau_R))/dt
        CASE(2)
          gamma_L = 2.*tau_L/(2.*tau_L+dt)
          gamma_R = 2.*tau_R/(2.*tau_R+dt)
        END SELECT
        UTemp_L = gamma_L*U_L(:,Count_1,Count_2) + (1.-gamma_L)*fTarget_L
        UTemp_R = gamma_R*U_R(:,Count_1,Count_2) + (1.-gamma_R)*fTarget_R
      END IF

      DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
        upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
        Velo= n_loc(1)*DVMVelos(iVel,1) + n_loc(2)*DVMVelos(jVel,2) + n_loc(3)*DVMVelos(kVel,3)
        F(upos,Count_1,Count_2)=0.5*((Velo+abs(Velo))*Utemp_L(upos)+(Velo-abs(Velo))*Utemp_R(upos))
        IF (DVMDim.LT.3) THEN
          F(PP_nVar_FV/2+upos,Count_1,Count_2) = &
            0.5*((Velo+abs(Velo))*Utemp_L(PP_nVar_FV/2+upos)+(Velo-abs(Velo))*Utemp_R(PP_nVar_FV/2+upos))
        END IF
      END DO; END DO; END DO;
    END DO
  END DO
END SUBROUTINE Riemann


END MODULE MOD_Riemann
