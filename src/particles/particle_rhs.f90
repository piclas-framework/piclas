!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

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

INTERFACE RELATIVISTIC_PUSH
  MODULE PROCEDURE RELATIVISTIC_PUSH
END INTERFACE

INTERFACE NON_RELATIVISTIC_PUSH
  MODULE PROCEDURE NON_RELATIVISTIC_PUSH
END INTERFACE

INTERFACE PartVeloToImp
  MODULE PROCEDURE PartVeloToImp
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC            :: CalcPartRHS
PUBLIC            :: SLOW_RELATIVISTIC_PUSH
PUBLIC            :: FAST_RELATIVISTIC_PUSH
PUBLIC            :: RELATIVISTIC_PUSH
PUBLIC            :: NON_RELATIVISTIC_PUSH
PUBLIC            :: PartVeloToImp
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
REAL                             :: E(1:3)
#if (PP_nVar==8)
REAL                             :: B(1:3)
#endif
REAL                             :: qmt, LorentzFac, velosq
REAL                             :: ax, ay, az, bx, by, bz, dx, dy, snx, sny, snz
!===================================================================================================================================

! Lorentzforce
Pt(1:PDM%ParticleVecLength,:)=0.
SELECT CASE(PartLorentzType)
  CASE(0) ! default
    ! non-relativistic
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Don't push neutral particles!
        IF(.NOT.PUSHPARTICLE(iPart)) CYCLE
        Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE(1) ! default
    ! constant Lorentz factor over time step
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Don't push neutral particles!
        IF(.NOT.PUSHPARTICLE(iPart)) CYCLE
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE(2)
  ! Lorentz-Pusher, wrong
  ! prevent particles from acceleration above speed of light
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Don't push neutral particles!
        IF(.NOT.PUSHPARTICLE(iPart)) CYCLE
        ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        velosq = PartState(iPart,4) * PartState(iPart,4) &
               + PartState(iPart,5) * PartState(iPart,5) &
               + PartState(iPart,6) * PartState(iPart,6)
        IF(velosq.GT.c2) THEN
          IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
          CALL abort(&
          __STAMP__&
          ,'Particle is faster than the speed of light. Maybe reducing the time step would help. Particle-Nr., velosq/c2:'&
          ,iPart,velosq/c2)
        END IF
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
        ! Don't push neutral particles!
        IF(.NOT.PUSHPARTICLE(iPart)) CYCLE
        ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        Pt(iPart,1:3)=FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE(31)
    ! derivation of relativistic equation of motion
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Don't push neutral particles!
        IF(.NOT.PUSHPARTICLE(iPart)) CYCLE
        ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        Pt(iPart,1:3)=ACCELERATION_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      END IF
    END DO
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'This Type of Lorentz-force calculation is not implemented:.',PartLorentzType,999.)
END SELECT

END SUBROUTINE CalcPartRHS


FUNCTION SLOW_RELATIVISTIC_PUSH(PartID,FieldAtParticle)
!===================================================================================================================================
! Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : abort
#ifdef MPI
USE MOD_Globals,           ONLY : MyRank
#endif
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
REAL                :: E(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3)
#endif
!===================================================================================================================================

velosq = PartState(PartID,4) * PartState(PartID,4) &
       + PartState(PartID,5) * PartState(PartID,5) &
       + PartState(PartID,6) * PartState(PartID,6)

IF(velosq.GT.c2) THEN
 IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
 IPWRITE(*,*) ' Species-ID',PartSpecies(PartID)
  CALL abort(&
  __STAMP__&
  ,'Particle is faster than the speed of light. Maybe reducing the time step would help. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
END IF

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
#ifdef MPI
USE MOD_Globals,           ONLY : MyRank
#endif
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
REAL                :: E(1:3),Pt(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3)
#endif
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
 IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
 IPWRITE(*,*) ' Species-ID',PartSpecies(PartID)
 IPWRITE(*,*) ' x=',PartState(PartID,1),' y=',PartState(PartID,2),' z=',PartState(PartID,3)
 CALL abort(&
  __STAMP__&
  ,'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
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


FUNCTION ACCELERATION_RELATIVISTIC_PUSH(PartID,FieldAtParticle)
!===================================================================================================================================
! Returns the relativistic acceleration a = dv/dt
! see W. Rindler, Relativity: Special, General, and Cosmological, 2006, Oxford University Press, New York, p.125
!
! CAUTION: This routines is used for HDG in combination with magnetic (external) fields
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : abort
#ifdef MPI
USE MOD_Globals,           ONLY : MyRank
#endif
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
REAL                :: ACCELERATION_RELATIVISTIC_PUSH(1:3) ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq,F(1:3)
!===================================================================================================================================
ASSOCIATE (&
      qmt => Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC ,& ! charge/m_0
      v1  => PartState(PartID,4) ,& ! Velocity in x
      v2  => PartState(PartID,5) ,& ! Velocity in y
      v3  => PartState(PartID,6) ,& ! Velocity in z
      E1  => FieldAtParticle(1)  ,& ! Electric field in x
      E2  => FieldAtParticle(2)  ,& ! Electric field in y
      E3  => FieldAtParticle(3)  ,& ! Electric field in z
      B1  => FieldAtParticle(4)  ,& ! Magnetic field in x
      B2  => FieldAtParticle(5)  ,& ! Magnetic field in y
      B3  => FieldAtParticle(6)   & ! Magnetic field in z
      )

  ! Check squared velocity with c^2
  velosq = v1*v1 + v2*v2 + v3*v3
  IF(velosq.GT.c2) THEN
    IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
    IPWRITE(*,*) ' Species-ID',PartSpecies(PartID)
    IPWRITE(*,*) ' x=',PartState(PartID,1),' y=',PartState(PartID,2),' z=',PartState(PartID,3)
    CALL abort(&
        __STAMP__&
        ,'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
  END IF
  ASSOCIATE ( gammas => SQRT(1.0-velosq*c2_inv) ) ! Inverse of Lorentz factor

! The following ifdef is commented out in order to push particles for the HDG solver in combination with a magnetic field
!#if (PP_nVar==8)
    F(1) = E1 + v2 * B3 - v3 * B2
    F(2) = E2 + v3 * B1 - v1 * B3
    F(3) = E3 + v1 * B2 - v2 * B1
!#else
!    F(1) = E1
!    F(2) = E2
!    F(3) = E3
!#endif

    ! Calculate the acceleration
    ACCELERATION_RELATIVISTIC_PUSH = gammas * qmt * ( F - DOT_PRODUCT(F,PartState(PartID,4:6))*PartState(PartID,4:6)*c2_inv )
  END ASSOCIATE
END ASSOCIATE

END FUNCTION ACCELERATION_RELATIVISTIC_PUSH



FUNCTION RELATIVISTIC_PUSH(PartID,FieldAtParticle,LorentzFacInvIn)
!===================================================================================================================================
! full relativistic push in case that the particle velocity*gamma is updated in time
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : cross
USE MOD_Particle_Vars,     ONLY : PartState, Species, PartSpecies
USE MOD_Equation_Vars,     ONLY : c2_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL::LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: RELATIVISTIC_PUSH(1:3) ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: LorentzFacInv,qmt
REAL                :: E(1:3),Pt(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3),Velo(3)
#endif
!===================================================================================================================================

IF(PRESENT(LorentzFacInvIn))THEN
  LorentzFacInv=LorentzFacInvIn
ELSE
  LorentzFacInv=1.0+DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))*c2_inv
  LorentzFacInv=1.0/SQRT(LorentzFacInv)
END IF

qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC

E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
Velo(1) = LorentzFacInv*PartState(PartID,4)
Velo(2) = LorentzFacInv*PartState(PartID,5)
Velo(3) = LorentzFacInv*PartState(PartID,6)
Pt = E + CROSS(Velo,B)
!Pt(1) = E(1) + LorentzFacInv*(PartState(PartID,5) * B(3) - PartState(PartID,6) * B(2))
!Pt(2) = E(2) + LorentzFacInv*(PartState(PartID,6) * B(1) - PartState(PartID,4) * B(3))
!Pt(3) = E(3) + LorentzFacInv*(PartState(PartID,4) * B(2) - PartState(PartID,5) * B(1))
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

RELATIVISTIC_PUSH = Pt

END FUNCTION RELATIVISTIC_PUSH


FUNCTION NON_RELATIVISTIC_PUSH(PartID,FieldAtParticle)
!===================================================================================================================================
! NON relativistic push
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : cross
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies
#if (PP_nVar==8)
USE MOD_Particle_Vars,     ONLY : PartState
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: NON_RELATIVISTIC_PUSH(1:3) ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: E(1:3),Pt(1:3),qmt
#if (PP_nVar==8)
REAL                :: B(1:3),Velo(3)
#endif
!===================================================================================================================================

qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC

E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
Velo(1) = PartState(PartID,4)
Velo(2) = PartState(PartID,5)
Velo(3) = PartState(PartID,6)
Pt = E + CROSS(Velo,B)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

NON_RELATIVISTIC_PUSH = Pt

END FUNCTION NON_RELATIVISTIC_PUSH


SUBROUTINE PartVeloToImp(VeloToImp,doParticle_In)
!===================================================================================================================================
! map the particle velocity to gamma*velocity
! or
! gamma*velocity to velocity
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY : PDM, PartState, PartLorentzType
USE MOD_Equation_Vars,          ONLY : c2_inv
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: VeloToImp
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: doParticle(1:PDM%ParticleVecLength)
INTEGER                          :: iPart
REAL                             :: LorentzFac,LorentzFacInv
!===================================================================================================================================

IF(PartLorentzType.NE.5) RETURN

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

IF(VeloToImp)THEN
  DO iPart=1,PDM%ParticleVecLength
    IF(DoParticle(iPart))THEN
      LorentzFac=1.0-DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv
      LorentzFac=1.0/SQRT(LorentzFac)
      PartState(iPart,4) = LorentzFac*PartState(iPart,4)
      PartState(iPart,5) = LorentzFac*PartState(iPart,5)
      PartState(iPart,6) = LorentzFac*PartState(iPart,6)
    END IF ! DoParticle
  END DO ! iPart
ELSE
  DO iPart=1,PDM%ParticleVecLength
    IF(DoParticle(iPart))THEN
      LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv
      LorentzFacInv=1.0/SQRT(LorentzFacInv)
      PartState(iPart,4) = LorentzFacInv*PartState(iPart,4)
      PartState(iPart,5) = LorentzFacInv*PartState(iPart,5)
      PartState(iPart,6) = LorentzFacInv*PartState(iPart,6)
    END IF ! DoParticle
  END DO ! iPart
END IF

END SUBROUTINE PartVeloToImp

END MODULE MOD_part_RHS
