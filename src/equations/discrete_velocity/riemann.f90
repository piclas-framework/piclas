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
USE MOD_PreProc
USE MOD_DistFunc,        ONLY: MacroValuesFromDistribution, TargetDistribution, MoleculeRelaxEnergy
USE MOD_Equation_Vars_FV,ONLY: DVMDim, DVMSpecData, DVMnSpecies, DVMMethod, DVMColl, DVMnMacro
USE MOD_TimeDisc_Vars,   ONLY: dt
USE MOD_Globals,         ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar_FV),INTENT(IN)            :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar_FV)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: MacroVal_L(DVMnMacro,DVMnSpecies+1), MacroVal_R(DVMnMacro,DVMnSpecies+1)
REAL                                             :: tau_L, tau_R, Velo, rho_L, rho_R, Pr_L, Pr_R
REAL,ALLOCATABLE                                 :: fTarget_L(:), fTarget_R(:), UTemp_L(:), UTemp_R(:)
REAL                                             :: gamma_R, gamma_L, relaxFac
INTEGER                                          :: iVel, jVel, kVel, upos, iSpec, vFirstID, vLastID
REAL                                             :: Erot_L(DVMnSpecies+1),Erot_R(DVMnSpecies+1)
REAL                                             :: Evib_L(DVMnSpecies+1),Evib_R(DVMnSpecies+1)
REAL                                             :: ErelaxTrans_L,ErelaxTrans_R
REAL                                             :: ErelaxRot_L(DVMnSpecies),ErelaxRot_R(DVMnSpecies)
REAL                                             :: ErelaxVib_L(DVMnSpecies),ErelaxVib_R(DVMnSpecies)
!===================================================================================================================================
IF (DVMColl.AND.DVMMethod.GT.0) THEN
  CALL MacroValuesFromDistribution(MacroVal_L,U_L,dt/2.,tau_L,1,MassDensity=rho_L,PrandtlNumber=Pr_L,Erot=Erot_L,Evib=Evib_L)
  CALL MacroValuesFromDistribution(MacroVal_R,U_R,dt/2.,tau_R,1,MassDensity=rho_R,PrandtlNumber=Pr_R,Erot=Erot_R,Evib=Evib_R)
  CALL MoleculeRelaxEnergy(ErelaxTrans_L,ErelaxRot_L,ErelaxVib_L,MacroVal_L(5,DVMnSpecies+1),&
                            ERot_L(1:DVMnSpecies),Evib_L(1:DVMnSpecies),Pr_L)
  CALL MoleculeRelaxEnergy(ErelaxTrans_R,ErelaxRot_R,ErelaxVib_R,MacroVal_R(5,DVMnSpecies+1),&
                            ERot_R(1:DVMnSpecies),Evib_R(1:DVMnSpecies),Pr_R)
  SELECT CASE (DVMMethod)
  CASE(1)
    gamma_L = 0.
    IF (tau_L.GT.0.) THEN
      relaxFac = dt/2./tau_L
      IF(CHECKEXP(relaxFac)) THEN
        gamma_L = 2.*tau_L*(1.-EXP(-relaxFac))/dt
      END IF
    END IF
    gamma_R = 0.
    IF (tau_R.GT.0.) THEN
      relaxFac = dt/2./tau_R
      IF(CHECKEXP(relaxFac)) THEN
        gamma_R = 2.*tau_R*(1.-EXP(-relaxFac))/dt
      END IF
    END IF
  CASE(2)
    gamma_L = 2.*tau_L/(2.*tau_L+dt/2.)
    gamma_R = 2.*tau_R/(2.*tau_R+dt/2.)
  END SELECT
END IF
vFirstID=1
vLastID=0
DO iSpec=1,DVMnSpecies
  vLastID = vLastID + DVMSpecData(iSpec)%nVar
  ALLOCATE(fTarget_L(DVMSpecData(iSpec)%nVar))
  ALLOCATE(fTarget_R(DVMSpecData(iSpec)%nVar))
  ALLOCATE(UTemp_L(DVMSpecData(iSpec)%nVar))
  ALLOCATE(UTemp_R(DVMSpecData(iSpec)%nVar))
  IF (DVMColl.AND.DVMMethod.GT.0) THEN
    CALL TargetDistribution(MacroVal_L(:,DVMnSpecies+1),fTarget_L,iSpec,MacroVal_L(1,iSpec),rho_L,Pr_L,ErelaxTrans_L,ErelaxRot_L(iSpec),ErelaxVib_L(iSpec))
    CALL TargetDistribution(MacroVal_R(:,DVMnSpecies+1),fTarget_R,iSpec,MacroVal_R(1,iSpec),rho_R,Pr_R,ErelaxTrans_R,ErelaxRot_R(iSpec),ErelaxVib_R(iSpec))
    IF (dt.EQ.0.) THEN
      UTemp_L = 0.
      UTemp_R = 0.
    ELSE
      UTemp_L = gamma_L*U_L(vFirstID:vLastID) + (1.-gamma_L)*fTarget_L
      UTemp_R = gamma_R*U_R(vFirstID:vLastID) + (1.-gamma_R)*fTarget_R
    END IF
  ELSE ! first order method or no collisions
    UTemp_L = U_L(vFirstID:vLastID)
    UTemp_R = U_R(vFirstID:vLastID)
  END IF

  DO kVel=1, DVMSpecData(iSpec)%nVelos(3); DO jVel=1, DVMSpecData(iSpec)%nVelos(2); DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
    upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
    Velo= nv(1)*DVMSpecData(iSpec)%Velos(iVel,1) &
        + nv(2)*DVMSpecData(iSpec)%Velos(jVel,2) &
        + nv(3)*DVMSpecData(iSpec)%Velos(kVel,3)
    F(vFirstID+upos-1)=0.5*((Velo+abs(Velo))*Utemp_L(upos)+(Velo-abs(Velo))*Utemp_R(upos))
    IF (DVMDim.LT.3) THEN
      F(vFirstID+DVMSpecData(iSpec)%nVarReduced+upos-1) = &
                                              0.5*((Velo+abs(Velo))*Utemp_L(DVMSpecData(iSpec)%nVarReduced+upos) &
                                                  + (Velo-abs(Velo))*Utemp_R(DVMSpecData(iSpec)%nVarReduced+upos))
    END IF
    IF (DVMSpecData(iSpec)%Xi_Rot.GT.0.) THEN
      ! rotational energy reduced distribution
      F(vFirstID+DVMSpecData(iSpec)%nVarErotStart+upos-1) = &
                                              0.5*((Velo+abs(Velo))*Utemp_L(DVMSpecData(iSpec)%nVarErotStart+upos) &
                                                  + (Velo-abs(Velo))*Utemp_R(DVMSpecData(iSpec)%nVarErotStart+upos))
    END IF
    IF (DVMSpecData(iSpec)%T_Vib.GT.0.) THEN
      ! vibrational energy reduced distribution
      F(vFirstID+DVMSpecData(iSpec)%nVarEvibStart+upos-1) = &
                                              0.5*((Velo+abs(Velo))*Utemp_L(DVMSpecData(iSpec)%nVarEvibStart+upos) &
                                                  + (Velo-abs(Velo))*Utemp_R(DVMSpecData(iSpec)%nVarEvibStart+upos))
    END IF
  END DO; END DO; END DO;
  DEALLOCATE(fTarget_L)
  DEALLOCATE(fTarget_R)
  DEALLOCATE(UTemp_L)
  DEALLOCATE(UTemp_R)
  vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
END DO ! iSpec
END SUBROUTINE Riemann


END MODULE MOD_Riemann
