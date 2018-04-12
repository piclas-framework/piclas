#include "boltzplatz.h"

MODULE MOD_ESBGK_Phase
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ComputePhasePotential
  MODULE PROCEDURE ComputePhasePotential
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ComputePhasePotential,EvalPhaseForce, InterpolatePhaseForceToParticle
!===================================================================================================================================

CONTAINS

SUBROUTINE ComputePhasePotential()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_PreProc      ,ONLY: PP_N
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_ESBGK_Vars   ,ONLY: BGKPhasePreFactor, BGKCriticalDens, BGKPhiPhase , BGKPhiPhase2
USE MOD_PICDepo_Vars ,ONLY: PartSource
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
REAL                  :: m,l,k
!===================================================================================================================================
DO iElem = 1, nElems
  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    BGKPhiPhase(k,l,m,iElem) = BGKPhasePreFactor*(1. - EXP(-PartSource( 4 ,k,l,m,iElem)/BGKCriticalDens)) !**2.0
!      BGKPhiPhase(k,l,m,iElem) = BGKPhasePreFactor*(1. - EXP(-PartSource( 4 ,k,l,m,iElem)/BGKCriticalDens)) &
!          *(1. - EXP(-PartSource( 1 ,k,l,m,iElem)/BGKCriticalDens))       
!    BGKPhiPhase(k,l,m,iElem) = BGKPhasePreFactor*PartSource( 4 ,k,l,m,iElem)/BGKCriticalDens*PartSource( 1 ,k,l,m,iElem)/BGKCriticalDens
!    BGKPhiPhase2(k,l,m,iElem) = BGKPhasePreFactor*PartSource( 1 ,k,l,m,iElem)/BGKCriticalDens
  END DO; END DO; END DO
END DO
END SUBROUTINE ComputePhasePotential


SUBROUTINE EvalPhaseForce()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_PreProc
USE MOD_ESBGK_Vars ,ONLY: BGKPhiPhase, BGKForcePhase, D_BGK,  BGKPhiPhase2, BGKForcePhase2
USE MOD_Mesh_Vars  ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde, sJ
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)   :: gradPhi_xi,gradPhi_eta,gradPhi_zeta
INTEGER                                :: i,j,k,l,iElem
INTEGER,SAVE                           :: N_old=0
!===================================================================================================================================

DO iElem = 1, PP_nElems
  ! Compute the gradient in the reference system
  gradPhi_xi  = 0.
  gradPhi_eta = 0.
  gradPhi_zeta= 0.
  DO l=0,PP_N
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          gradPhi_xi(i,j,k)  = gradPhi_xi(i,j,k)   + D_BGK(i,l) * BGKPhiPhase(l,j,k,iElem)
          gradPhi_eta(i,j,k) = gradPhi_eta(i,j,k)  + D_BGK(j,l) * BGKPhiPhase(i,l,k,iElem)
          gradPhi_zeta(i,j,k)= gradPhi_zeta(i,j,k) + D_BGK(k,l) * BGKPhiPhase(i,j,l,iElem)
       END DO ! i 
     END DO ! j 
    END DO ! k 
  END DO ! l 
  ! Transform the gradients from the reference system to the xyz-System. Only exact for cartesian mesh!
  DO k=0,N
    DO j=0,N
      DO i=0,N
        BGKForcePhase(1,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(1,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(1,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(1,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
        BGKForcePhase(2,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(2,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(2,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(2,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
        BGKForcePhase(3,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(3,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(3,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(3,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
      END DO ! i 
    END DO ! j 
  END DO ! k 
END DO

!DO iElem = 1, PP_nElems
!  ! Compute the gradient in the reference system
!  gradPhi_xi  = 0.
!  gradPhi_eta = 0.
!  gradPhi_zeta= 0.
!  DO l=0,PP_N
!    DO k=0,PP_N
!      DO j=0,PP_N
!        DO i=0,PP_N
!          gradPhi_xi(i,j,k)  = gradPhi_xi(i,j,k)   + D_BGK(i,l) * BGKPhiPhase2(l,j,k,iElem)
!          gradPhi_eta(i,j,k) = gradPhi_eta(i,j,k)  + D_BGK(j,l) * BGKPhiPhase2(i,l,k,iElem)
!          gradPhi_zeta(i,j,k)= gradPhi_zeta(i,j,k) + D_BGK(k,l) * BGKPhiPhase2(i,j,l,iElem)
!       END DO ! i 
!     END DO ! j 
!    END DO ! k 
!  END DO ! l 
!  ! Transform the gradients from the reference system to the xyz-System. Only exact for cartesian mesh!
!  DO k=0,N
!    DO j=0,N
!      DO i=0,N
!        BGKForcePhase2(1,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &   
!                          Metrics_fTilde(1,i,j,k,iElem) * gradPhi_xi(i,j,k)   + & 
!                          Metrics_gTilde(1,i,j,k,iElem) * gradPhi_eta(i,j,k)  + & 
!                          Metrics_hTilde(1,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )   
!        BGKForcePhase2(2,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &   
!                          Metrics_fTilde(2,i,j,k,iElem) * gradPhi_xi(i,j,k)   + & 
!                          Metrics_gTilde(2,i,j,k,iElem) * gradPhi_eta(i,j,k)  + & 
!                          Metrics_hTilde(2,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )   
!        BGKForcePhase2(3,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &   
!                          Metrics_fTilde(3,i,j,k,iElem) * gradPhi_xi(i,j,k)   + & 
!                          Metrics_gTilde(3,i,j,k,iElem) * gradPhi_eta(i,j,k)  + & 
!                          Metrics_hTilde(3,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )   
!      END DO ! i 
!    END DO ! j 
!  END DO ! k 
!END DO
END SUBROUTINE EvalPhaseForce


SUBROUTINE InterpolatePhaseForceToParticle()
!===================================================================================================================================
!> interpolates phase forces to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DG_Vars               ,ONLY: U
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_PICInterpolation_Vars
USE MOD_PreProc
USE MOD_Eval_xyz
USE MOD_ESBGK_Vars            ,ONLY: BGKForcePhase, BGKForcePhase2
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                    
REAL                             :: Pos(3), alpha1, alpha2, alpha3
REAL                             :: field(3)
INTEGER                          :: m,iPart,iElem, kk, ll, mm, N_FieldDOF
!===================================================================================================================================
           ! skip if no self fields are calculated
SELECT CASE(TRIM(InterpolationType))
CASE('nearest_blurrycenter')
  ! add fields to fields at particle position
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      m = INT(PP_N/2)+1
      Pt(iPart,1:3) = BGKForcePhase(1:3,m,m,m,iElem)
    END IF
  END DO
CASE('particle_position')
  DO iPart = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      Pos = PartState(iPart,1:3)
        !--- evaluate at Particle position
!        CALL eval_xyz_curved(Pos,3,PP_N,BGKForcePhase(1:3,:,:,:,iElem),field,iElem)
!        IF (PartSpecies(iPart).EQ.1) THEN
!          CALL eval_xyz_curved(Pos,3,PP_N,BGKForcePhase2(1:3,:,:,:,iElem),field,iElem)
!        ELSE
        CALL eval_xyz_curved(Pos,3,PP_N,BGKForcePhase(1:3,:,:,:,iElem),field,iElem)
!        END IF
      Pt(iPart,1:3) = field(1:3)
    END IF
  END DO
CASE DEFAULT
  CALL abort(__STAMP__, &
      'ERROR: Unknown InterpolationType!')
END SELECT

RETURN
END SUBROUTINE InterpolatePhaseForceToParticle


END MODULE  MOD_ESBGK_Phase
