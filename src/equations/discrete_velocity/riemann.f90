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
USE MOD_PreProc ! 0
USE MOD_DistFunc,        ONLY: MacroValuesFromDistribution, MaxwellDistribution, MaxwellDistributionCons
USE MOD_DistFunc,        ONLY: ShakhovDistribution, ESBGKDistribution, GradDistributionPrandtl
USE MOD_DistFunc,        ONLY: SkewNormalDistribution, SkewtDistribution
USE MOD_Equation_Vars_FV,ONLY: DVMDim, DVMSpecData, DVMnSpecies, DVMMethod, DVMBGKModel, DVMColl
USE MOD_TimeDisc_Vars,   ONLY: dt
USE MOD_Globals,         ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar_FV,0:0,0:0),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:0,0:0)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar_FV,0:0,0:0)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),MacroVal_L(14,DVMnSpecies+1), MacroVal_R(14,DVMnSpecies+1), tau_L, tau_R
REAL                                             :: Velo
REAL,ALLOCATABLE                                 :: fTarget_L(:), fTarget_R(:), UTemp_L(:), UTemp_R(:)
REAL                                             :: gamma_R, gamma_L
INTEGER                                          :: Count_1,Count_2, iVel, jVel, kVel, upos, iSpec, vFirstID, vLastID
!===================================================================================================================================
! Gauss point i,j
  DO Count_2=0,0
    DO Count_1=0,0
      n_loc(:)=nv(:,Count_1,Count_2)
      IF (DVMColl) THEN
        CALL MacroValuesFromDistribution(MacroVal_L,U_L(:,Count_1,Count_2),dt/2.,tau_L,1)
        CALL MacroValuesFromDistribution(MacroVal_R,U_R(:,Count_1,Count_2),dt/2.,tau_R,1)
        SELECT CASE (DVMMethod)
        CASE(1)
          gamma_L = 2.*tau_L*(1.-EXP(-dt/2./tau_L))/dt
          gamma_R = 2.*tau_R*(1.-EXP(-dt/2./tau_R))/dt
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
        MacroVal_L(2:14,iSpec) = MacroVal_L(2:14,DVMnSpecies+1)
        MacroVal_R(2:14,iSpec) = MacroVal_R(2:14,DVMnSpecies+1)
        IF (DVMColl) THEN
          SELECT CASE (DVMBGKModel)
            CASE(1)
              CALL ESBGKDistribution(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL ESBGKDistribution(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(2)
              CALL ShakhovDistribution(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL ShakhovDistribution(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(3)
              CALL MaxwellDistribution(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL MaxwellDistribution(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(4)
              CALL MaxwellDistributionCons(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL MaxwellDistributionCons(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(5)
              CALL SkewNormalDistribution(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL SkewNormalDistribution(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(6)
              CALL SkewtDistribution(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL SkewtDistribution(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE(7)
              CALL GradDistributionPrandtl(MacroVal_L(:,iSpec),fTarget_L,iSpec)
              CALL GradDistributionPrandtl(MacroVal_R(:,iSpec),fTarget_R,iSpec)
            CASE DEFAULT
              CALL abort(__STAMP__,'DVM BGK Model not implemented.',999,999.)
          END SELECT
          IF (dt.EQ.0.) THEN
            UTemp_L = 0.
            UTemp_R = 0.
          ELSE
            UTemp_L = gamma_L*U_L(vFirstID:vLastID,Count_1,Count_2) + (1.-gamma_L)*fTarget_L
            UTemp_R = gamma_R*U_R(vFirstID:vLastID,Count_1,Count_2) + (1.-gamma_R)*fTarget_R
          END IF
        ELSE ! no collisions
          UTemp_L = U_L(vFirstID:vLastID,Count_1,Count_2)
          UTemp_R = U_R(vFirstID:vLastID,Count_1,Count_2)
        END IF

        DO kVel=1, DVMSpecData(iSpec)%nVelos(3); DO jVel=1, DVMSpecData(iSpec)%nVelos(2); DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
          upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
          Velo= n_loc(1)*DVMSpecData(iSpec)%Velos(iVel,1) &
              + n_loc(2)*DVMSpecData(iSpec)%Velos(jVel,2) &
              + n_loc(3)*DVMSpecData(iSpec)%Velos(kVel,3)
          F(vFirstID+upos-1,Count_1,Count_2)=0.5*((Velo+abs(Velo))*Utemp_L(upos)+(Velo-abs(Velo))*Utemp_R(upos))
          IF (DVMDim.LT.3) THEN
            F(vFirstID+DVMSpecData(iSpec)%nVar/2+upos-1,Count_1,Count_2) = &
                                                    0.5*((Velo+abs(Velo))*Utemp_L(DVMSpecData(iSpec)%nVar/2+upos) &
                                                       + (Velo-abs(Velo))*Utemp_R(DVMSpecData(iSpec)%nVar/2+upos))
          END IF
        END DO; END DO; END DO;
        DEALLOCATE(fTarget_L)
        DEALLOCATE(fTarget_R)
        DEALLOCATE(UTemp_L)
        DEALLOCATE(UTemp_R)
        vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
      END DO ! iSpec
    END DO
  END DO
END SUBROUTINE Riemann


END MODULE MOD_Riemann
