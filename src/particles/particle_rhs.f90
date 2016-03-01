#include "boltzplatz.h"

MODULE MOD_part_RHS                                                                                
!===================================================================================================================================
! Subroutine to compute the particle right hand side, therefore the acceleration due to the Lorentz-force with 
! respect to the Lorentz factor
!===================================================================================================================================
IMPLICIT NONE                                                                                    
PRIVATE                                                                                          

INTERFACE CalcPartRHS
  MODULE PROCEDURE CalcPartRHS
END INTERFACE

INTERFACE SLOW_RELATIVISTIC_PUSH
  MODULE PROCEDURE SLOW_RELATIVISTIC_PUSH
END INTERFACE

INTERFACE FAST_RELATIVISTIC_PUSH
  MODULE PROCEDURE FAST_RELATIVISTIC_PUSH
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC            :: CalcPartRHS                                                                 
PUBLIC            :: SLOW_RELATIVISTIC_PUSH
PUBLIC            :: FAST_RELATIVISTIC_PUSH
!----------------------------------------------------------------------------------------------------------------------------------
                                                                                                   
CONTAINS                                                                                           
                                                                                                   
SUBROUTINE CalcPartRHS()                                                                           
!===================================================================================================================================
! Computes the acceleration from the Lorentz force with respect to the species data and velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY : PDM, PartState, Pt, Species, PartSpecies, PartLorentzType
USE MOD_Equation_Vars,          ONLY : c2_inv, c ,c2
USE MOD_TimeDisc_Vars,          ONLY : dt
USE MOD_PICInterpolation_Vars,  ONLY : FieldAtParticle                                           
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE 
INTEGER                          :: iPart                                               
REAL                             :: E(1:3),B(1:3)                                                
REAL                             :: qmt, LorentzFac, velosq
REAL                             :: ax, ay, az, bx, by, bz, dx, dy, snx, sny, snz
!===================================================================================================================================

! Lorentzforce
Pt(1:PDM%ParticleVecLength,:)=0.
SELECT CASE(PartLorentzType)
  CASE(1) ! default
    ! constant Lorentz factor over time step
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE(2)
  ! Lorentz-Pusher, wrong
  ! prevent particles from acceleration above speed of light
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
    ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        velosq = PartState(iPart,4) * PartState(iPart,4) &
               + PartState(iPart,5) * PartState(iPart,5) &
               + PartState(iPart,6) * PartState(iPart,6)  
        IF(velosq.GT.c2) CALL abort(&
            __STAMP__,&
          'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',iPart,velosq/c2)
          ! MPF in ChargeIC and MassIC cancels out.
        qmt = Species(PartSpecies(iPart))%ChargeIC/Species(PartSpecies(iPart))%MassIC
        E(1:3) = FieldAtParticle(iPart,1:3) * qmt
#if (PP_nVar==8)
        B(1:3) = FieldAtParticle(iPart,4:6) * qmt
#endif
      ! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
        Pt(iPart,1) = E(1) + PartState(iPart,5) * B(3) - PartState(iPart,6) * B(2)
        Pt(iPart,2) = E(2) + PartState(iPart,6) * B(1) - PartState(iPart,4) * B(3)
        Pt(iPart,3) = E(3) + PartState(iPart,4) * B(2) - PartState(iPart,5) * B(1)
#else
        Pt(iPart,1) = E(1) 
        Pt(iPart,2) = E(2) 
        Pt(iPart,3) = E(3) 
#endif

        LorentzFac = 1/sqrt(1.0 - velosq * c2_inv)
        bx = Pt(iPart,1) *dt + LorentzFac * PartState(iPart,4)
        snx = sign(1.0,bx)
        bx = bx*bx*c2_inv
        bx = bx/(1+bx)

        by = Pt(iPart,2) *dt + LorentzFac * PartState(iPart,5)
        sny = sign(1.0,by)
        by = by*by*c2_inv
        by = by/(1+by)

        bz = Pt(iPart,3) *dt + LorentzFac * PartState(iPart,6)
        snz = sign(1.0,bz)
        bz = bz*bz*c2_inv
        bz = bz/(1+bz)

        dx = (bx-bx*bz)/(1-bx*bz)
        dy = (by-by*bz)/(1-by*bz)
          
        ax = (dx-dx*dy)/(1-dx*dy)
        ay = (dy-dy*ax)
        az = (bz*(1-ax-ay))

        Pt(iPart,1) = (snx * sqrt(ax)*c - PartState(iPart,4)) / dt
        Pt(iPart,2) = (sny * sqrt(ay)*c - PartState(iPart,5)) / dt
        Pt(iPart,3) = (snz * sqrt(az)*c - PartState(iPart,6)) / dt

      END IF
    END DO

  CASE(3)
    ! derivation of relativistic equation of motion
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        Pt(iPart,1:3)=FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE DEFAULT
    CALL abort(&
        __STAMP__,&
      'This Type of Lorentz-force calculation is not implemented:.',PartLorentzType,999.)
END SELECT

END SUBROUTINE CalcPartRHS


FUNCTION SLOW_RELATIVISTIC_PUSH(PartID,FieldAtParticle)
!===================================================================================================================================
! Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : abort
USE MOD_Particle_Vars,     ONLY : PartState, Species, PartSpecies
USE MOD_Equation_Vars,     ONLY : c2_inv, c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: SLOW_RELATIVISTIC_PUSH(1:3) ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq, LorentzFac,qmt
REAL                :: E(1:3),B(1:3)                                                
!===================================================================================================================================

velosq = PartState(PartID,4) * PartState(PartID,4) &
       + PartState(PartID,5) * PartState(PartID,5) &
       + PartState(PartID,6) * PartState(PartID,6)  

IF(velosq.GT.c2) CALL abort(&
    __STAMP__,&
  'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)

! MPF in ChargeIC and MassIC cancels out.
LorentzFac = (SQRT(1.0 - velosq * c2_inv))
qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC * LorentzFac
E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
SLOW_RELATIVISTIC_PUSH(1) = E(1) + PartState(PartID,5) * B(3) - PartState(PartID,6) * B(2)
SLOW_RELATIVISTIC_PUSH(2) = E(2) + PartState(PartID,6) * B(1) - PartState(PartID,4) * B(3)
SLOW_RELATIVISTIC_PUSH(3) = E(3) + PartState(PartID,4) * B(2) - PartState(PartID,5) * B(1)
#else
SLOW_RELATIVISTIC_PUSH(1) = E(1) 
SLOW_RELATIVISTIC_PUSH(2) = E(2) 
SLOW_RELATIVISTIC_PUSH(3) = E(3) 
#endif

END FUNCTION SLOW_RELATIVISTIC_PUSH


FUNCTION FAST_RELATIVISTIC_PUSH(PartID,FieldAtParticle)
!===================================================================================================================================
! Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : abort
USE MOD_Particle_Vars,     ONLY : PartState, Species, PartSpecies
USE MOD_Equation_Vars,     ONLY : c2_inv, c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: FAST_RELATIVISTIC_PUSH(1:3) ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq, LorentzFac,qmt
REAL                :: E(1:3),B(1:3),Pt(1:3)
REAL                :: LorentzFac2,LorentzFac3, v1s,v2s,v3s, Vinv(3,3), v1,v2,v3, normfac
!===================================================================================================================================

! required helps
v1  = PartState(PartID,4)
v2  = PartState(PartID,5)
v3  = PartState(PartID,6)

v1s = v1*v1
v2s = v2*v2
v3s = v3*v3
velosq = v1s+v2s+v3s
IF(velosq.GT.c2) THEN
 CALL abort(&
     __STAMP__,&
 'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
END IF

LorentzFac=SQRT(1.0 - velosq*c2_inv)
LorentzFac2=LorentzFac*LorentzFac
LorentzFac3=LorentzFac2*LorentzFac
normfac=1.0/(c2*LorentzFac2 + velosq)
! define inverted matrix
Vinv(1,1) = (c2*LorentzFac3 + (v2s+v3s)*LorentzFac)*normfac
Vinv(1,2) =-(v1*v2*LorentzFac)*normfac
Vinv(1,3) =-(v1*v3*LorentzFac)*normfac
Vinv(2,1) = Vinv(1,2)
Vinv(2,2) = (c2*LorentzFac3 + (v1s+v3s)*LorentzFac)*normfac
Vinv(2,3) =-(v2*v3*LorentzFac)*normfac
Vinv(3,1) = Vinv(1,3)
Vinv(3,2) = Vinv(2,3)
Vinv(3,3) = (c2*LorentzFac3 + (v1s+v2s)*LorentzFac)*normfac

qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC

E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
Pt(1) = E(1) + PartState(PartID,5) * B(3) - PartState(PartID,6) * B(2)
Pt(2) = E(2) + PartState(PartID,6) * B(1) - PartState(PartID,4) * B(3)
Pt(3) = E(3) + PartState(PartID,4) * B(2) - PartState(PartID,5) * B(1)
#else
Pt(1) = E(1) 
Pt(2) = E(2) 
Pt(3) = E(3) 
#endif

FAST_RELATIVISTIC_PUSH = MATMUL(Vinv,Pt)

END FUNCTION FAST_RELATIVISTIC_PUSH

END MODULE MOD_part_RHS
