!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
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

INTERFACE CalcPartRHSSingleParticle
  MODULE PROCEDURE CalcPartRHSSingleParticle
END INTERFACE

INTERFACE PartVeloToImp
  MODULE PROCEDURE PartVeloToImp
END INTERFACE

INTERFACE PartRHS
  PROCEDURE PartRHS
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalcPartRHS
PUBLIC :: PartVeloToImp
PUBLIC :: PartRHS
PUBLIC :: CalcPartRHSSingleParticle
PUBLIC :: CalcPartRHSRotRefFrame
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE PartRHSInterface(PartID,FieldAtParticle,Push,LorentzFacInvIn)
    INTEGER,INTENT(IN)              :: PartID          ! Particle ID
    REAL,DIMENSION(1:6),INTENT(IN)  :: FieldAtParticle ! Electric and magnetic fields E, B at the position of the particle
    REAL,DIMENSION(1:3),INTENT(OUT) :: Push            ! dP/dt: acceleration of the particle (change of momentum)
    REAL,INTENT(IN),OPTIONAL        :: LorentzFacInvIn ! Reciprocal Lorentz factor
  END SUBROUTINE
END INTERFACE

PROCEDURE(PartRHSInterface),POINTER :: PartRHS    !< pointer defining the standard inner Riemann solver

INTEGER,PARAMETER      :: PRM_PART_RHS_NR  = 0   ! non-relativistic
INTEGER,PARAMETER      :: PRM_PART_RHS_D   = 1   ! default
INTEGER,PARAMETER      :: PRM_PART_RHS_W   = 2   ! wrong
INTEGER,PARAMETER      :: PRM_PART_RHS_RN  = 3   ! relativistic-new
INTEGER,PARAMETER      :: PRM_PART_RHS_REM = 31  ! relativistic-EM (electromagnetic)
INTEGER,PARAMETER      :: PRM_PART_RHS_RM  = 5   ! relativistic, momentum-based
INTEGER,PARAMETER      :: PRM_PART_RHS_CEM = 9   ! constant-EM (acceleration due to an electro-magnetic field that is constant)


INTERFACE InitPartRHS
  MODULE PROCEDURE InitPartRHS
END INTERFACE

PUBLIC :: InitPartRHS
!==================================================================================================================================

PUBLIC :: DefineParametersParticleRHS
CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersParticleRHS()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Particle RHS")
CALL prms%CreateIntFromStringOption('Part-LorentzType', "Lorentz force calculation for charged particles: "//&
                                                        "non-relativistic ("//TRIM(int2strf(PRM_PART_RHS_NR))//"), "//&
                                                        "default ("//TRIM(int2strf(PRM_PART_RHS_D))//"), "//&
                                                        "wrong ("//TRIM(int2strf(PRM_PART_RHS_W))//"), "//&
                                                        "relativistic-new ("//TRIM(int2strf(PRM_PART_RHS_RN))//"), "//&
                                                        "relativistic-EM ("//TRIM(int2strf(PRM_PART_RHS_REM))//"), "//&
                                                        "relativistic-momentum ("//TRIM(int2strf(PRM_PART_RHS_RM))//"), "//&
                                                        "constant-EM ("//TRIM(int2strf(PRM_PART_RHS_CEM))//")", "non-relativistic")
CALL addStrListEntry('Part-LorentzType' , 'non-relativistic'      , PRM_PART_RHS_NR)
CALL addStrListEntry('Part-LorentzType' , 'default'               , PRM_PART_RHS_D)
CALL addStrListEntry('Part-LorentzType' , 'wrong'                 , PRM_PART_RHS_W)
CALL addStrListEntry('Part-LorentzType' , 'relativistic-new'      , PRM_PART_RHS_RN)
CALL addStrListEntry('Part-LorentzType' , 'relativistic-EM'       , PRM_PART_RHS_REM)
CALL addStrListEntry('Part-LorentzType' , 'relativistic-momentum' , PRM_PART_RHS_RM)
CALL addStrListEntry('Part-LorentzType' , 'constant-EM'           , PRM_PART_RHS_CEM)
END SUBROUTINE DefineParametersParticleRHS


!==================================================================================================================================!
!> Initialize particle RHS functions
!==================================================================================================================================!
SUBROUTINE InitPartRHS()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools   ,ONLY: GETINTFROMSTR
USE MOD_Particle_Vars ,ONLY: PartLorentzType
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: dummy(1:3)
!==================================================================================================================================
PartLorentzType = GETINTFROMSTR('Part-LorentzType')
SELECT CASE(PartLorentzType)
CASE(PRM_PART_RHS_NR) ! 0
  PartRHS => PartRHS_NR
CASE(PRM_PART_RHS_D) ! 1
  PartRHS => PartRHS_D
CASE(PRM_PART_RHS_W) ! 2
  PartRHS => PartRHS_W
CASE(PRM_PART_RHS_RN) ! 3
  PartRHS => PartRHS_RN
CASE(PRM_PART_RHS_REM) ! 31
  PartRHS => PartRHS_REM
CASE(PRM_PART_RHS_RM) ! 5
  PartRHS => PartRHS_RM
CASE(PRM_PART_RHS_CEM) ! 9
  PartRHS => PartRHS_CEM
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Part-LorentzType-new not defined!')
END SELECT

! Suppress compiler warning
RETURN
CALL PartRHS_NR(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_D(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_W(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_RN(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_REM(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_RM(0,(/0.,0.,0.,0.,0.,0./),dummy)
CALL PartRHS_CEM(0,(/0.,0.,0.,0.,0.,0./),dummy)
END SUBROUTINE InitPartRHS


SUBROUTINE CalcPartRHS()
!===================================================================================================================================
! Computes the acceleration from the Lorentz force with respect to the species data and velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PDM,Pt
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Part_Tools            ,ONLY: isPushParticle
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
INTEGER                          :: iPart
!===================================================================================================================================
! Loop all particles and call particle right-hand-side calculation
DO iPart = 1,PDM%ParticleVecLength
  ! Particle is inside and not a neutral particle
  IF(PDM%ParticleInside(iPart))THEN
     IF(isPushParticle(iPart))THEN
       CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
       CYCLE
     END IF ! isPushParticle(iPart)
  END IF ! PDM%ParticleInside(iPart)
  ! Pt(:,iPart)=0.
END DO
END SUBROUTINE CalcPartRHS


SUBROUTINE CalcPartRHSSingleParticle(iPart)
!===================================================================================================================================
! Computes the acceleration from the Lorentz force with respect to the species data and velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PDM,Pt
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
INTEGER,INTENT(IN)            :: iPart
!===================================================================================================================================
! Particle is inside and not a neutral particle
IF(PDM%ParticleInside(iPart))THEN
  CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
  RETURN
END IF ! PDM%ParticleInside(iPart)
Pt(:,iPart)=0.
END SUBROUTINE CalcPartRHSSingleParticle


PPURE SUBROUTINE PartRHS_NR(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'non-relativistic'
! Particle Right-Hand-Side: Non-relativistic push
! Former FUNCTION NON_RELATIVISTIC_PUSH
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
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: E(1:3),qmt
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
Velo(1) = PartState(4,PartID)
Velo(2) = PartState(5,PartID)
Velo(3) = PartState(6,PartID)
Pt = E + CROSS(Velo,B)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

! Suppress compiler warning
RETURN
qmt=LorentzFacInvIn ! dummy statement
END SUBROUTINE PartRHS_NR


SUBROUTINE PartRHS_D(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'default'
! Particle Right-Hand-Side: relativistic push (old slow function)
! Former FUNCTION SLOW_RELATIVISTIC_PUSH
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort,DOTPRODUCT
#if USE_MPI
USE MOD_Globals       ,ONLY: MyRank
#endif
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars  ,ONLY: c2_inv, c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq, LorentzFac,qmt
REAL                :: E(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3)
#endif
!===================================================================================================================================
velosq = DOTPRODUCT(PartState(4:6,PartID))

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
Pt(1) = E(1) + PartState(5,PartID) * B(3) - PartState(6,PartID) * B(2)
Pt(2) = E(2) + PartState(6,PartID) * B(1) - PartState(4,PartID) * B(3)
Pt(3) = E(3) + PartState(4,PartID) * B(2) - PartState(5,PartID) * B(1)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

! Suppress compiler warning
RETURN
qmt=LorentzFacInvIn ! dummy statement
END SUBROUTINE PartRHS_D


SUBROUTINE PartRHS_W(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'wrong'
! Particle Right-Hand-Side: relativistic push (old wrong function)
! Lorentz-Pusher, wrong
! prevent particles from acceleration above speed of light
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars  ,ONLY: c2_inv, c ,c2
USE MOD_TimeDisc_Vars ,ONLY: dt
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Pt(1:3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
REAL                :: E(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3)
#endif
REAL                :: qmt, LorentzFac, velosq
REAL                :: ax, ay, az, bx, by, bz, dx, dy, snx, sny, snz
!===================================================================================================================================
! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
velosq = DOTPRODUCT(PartState(4:6,PartID))
IF(velosq.GT.c2) THEN
  IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
  CALL abort(&
      __STAMP__&
      ,'Particle is faster than the speed of light. Maybe reducing the time step would help. Particle-Nr., velosq/c2:'&
      ,PartID,velosq/c2)
END IF
! MPF in ChargeIC and MassIC cancels out.
qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC
E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
Pt(1) = E(1) + PartState(5,PartID) * B(3) - PartState(6,PartID) * B(2)
Pt(2) = E(2) + PartState(6,PartID) * B(1) - PartState(4,PartID) * B(3)
Pt(3) = E(3) + PartState(4,PartID) * B(2) - PartState(5,PartID) * B(1)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

LorentzFac = 1/sqrt(1.0 - velosq * c2_inv)
bx = Pt(1) *dt + LorentzFac * PartState(4,PartID)
snx = sign(1.0,bx)
bx = bx*bx*c2_inv
bx = bx/(1+bx)

by = Pt(2) *dt + LorentzFac * PartState(5,PartID)
sny = sign(1.0,by)
by = by*by*c2_inv
by = by/(1+by)

bz = Pt(3) *dt + LorentzFac * PartState(6,PartID)
snz = sign(1.0,bz)
bz = bz*bz*c2_inv
bz = bz/(1+bz)

dx = (bx-bx*bz)/(1-bx*bz)
dy = (by-by*bz)/(1-by*bz)

ax = (dx-dx*dy)/(1-dx*dy)
ay = (dy-dy*ax)
az = (bz*(1-ax-ay))

Pt(1) = (snx * sqrt(ax)*c - PartState(4,PartID)) / dt
Pt(2) = (sny * sqrt(ay)*c - PartState(5,PartID)) / dt
Pt(3) = (snz * sqrt(az)*c - PartState(6,PartID)) / dt

! Suppress compiler warning
RETURN
qmt=LorentzFacInvIn ! dummy statement
END SUBROUTINE PartRHS_W


SUBROUTINE PartRHS_RN(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'relativistic-new'
! Particle Right-Hand-Side: relativistic push (old fast function)
! Former FUNCTION FAST_RELATIVISTIC_PUSH
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
#if USE_MPI
USE MOD_Globals       ,ONLY: MyRank
#endif
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars  ,ONLY: c2_inv, c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq, LorentzFac,qmt
REAL                :: E(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3)
#endif
REAL                :: LorentzFac2,LorentzFac3, v1s,v2s,v3s, Vinv(3,3), v1,v2,v3, normfac
!===================================================================================================================================
! required helps
v1  = PartState(4,PartID)
v2  = PartState(5,PartID)
v3  = PartState(6,PartID)

v1s = v1*v1
v2s = v2*v2
v3s = v3*v3
velosq = v1s+v2s+v3s
IF(velosq.GT.c2) THEN
 IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
 IPWRITE(*,*) ' Species-ID',PartSpecies(PartID)
 IPWRITE(*,*) ' x=',PartState(1,PartID),' y=',PartState(2,PartID),' z=',PartState(3,PartID)
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
Pt(1) = E(1) + PartState(5,PartID) * B(3) - PartState(6,PartID) * B(2)
Pt(2) = E(2) + PartState(6,PartID) * B(1) - PartState(4,PartID) * B(3)
Pt(3) = E(3) + PartState(4,PartID) * B(2) - PartState(5,PartID) * B(1)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

Pt = MATMUL(Vinv,Pt)

! Suppress compiler warning
RETURN
qmt=LorentzFacInvIn ! dummy statement
END SUBROUTINE PartRHS_RN


SUBROUTINE PartRHS_REM(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'relativistic-EM' (electromagnetic)
! Particle Right-Hand-Side: relativistic push (old fast function)
! Former FUNCTION ACCELERATION_RELATIVISTIC_PUSH
!
! Returns the relativistic acceleration a = dv/dt
! see W. Rindler, Relativity: Special, General, and Cosmological, 2006, Oxford University Press, New York, p.125
!
! CAUTION: This routines is used for HDG in combination with magnetic (external) fields
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
#if USE_MPI
USE MOD_Globals       ,ONLY: MyRank
#endif
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars  ,ONLY: c2_inv, c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq,F(1:3)
!===================================================================================================================================
ASSOCIATE (&
      qmt => Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC ,& ! charge/m_0
      v1  => PartState(4,PartID) ,& ! Velocity in x
      v2  => PartState(5,PartID) ,& ! Velocity in y
      v3  => PartState(6,PartID) ,& ! Velocity in z
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
    IPWRITE(*,*) ' x=',PartState(1,PartID),' y=',PartState(2,PartID),' z=',PartState(3,PartID)
    CALL abort(__STAMP__,'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
  END IF
  ASSOCIATE ( gammas => SQRT(1.0 - velosq*c2_inv) ) ! Inverse of Lorentz factor

    F(1) = E1 + v2 * B3 - v3 * B2
    F(2) = E2 + v3 * B1 - v1 * B3
    F(3) = E3 + v1 * B2 - v2 * B1

    ! Calculate the acceleration
    Pt = gammas * qmt * ( F - DOT_PRODUCT(F,PartState(4:6,PartID))*PartState(4:6,PartID)*c2_inv )
  END ASSOCIATE
END ASSOCIATE

! Suppress compiler warning
RETURN
velosq=LorentzFacInvIn ! dummy statement
END SUBROUTINE PartRHS_REM


SUBROUTINE PartRHS_RM(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'relativistic-momentum' (electromagnetic)
! Particle Right-Hand-Side: relativistic push (old fast function)
! Former FUNCTION RELATIVISTIC_PUSH
!
! full relativistic push in case that the particle velocity*gamma is updated in time
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: cross,DOTPRODUCT
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars  ,ONLY: c2_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)         :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: LorentzFacInv,qmt
REAL                :: E(1:3)
#if (PP_nVar==8)
REAL                :: B(1:3),Velo(3)
#endif
!===================================================================================================================================
! Calculate Lorentz factor gamma -> 1/gamma
IF(PRESENT(LorentzFacInvIn))THEN
  LorentzFacInv=LorentzFacInvIn
ELSE
  LorentzFacInv=1.0+DOTPRODUCT(PartState(4:6,PartID))*c2_inv
  LorentzFacInv=1.0/SQRT(LorentzFacInv)
END IF

qmt = Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC

E(1:3) = FieldAtParticle(1:3) * qmt
#if (PP_nVar==8)
B(1:3) = FieldAtParticle(4:6) * qmt
#endif
! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
Velo(1) = LorentzFacInv*PartState(4,PartID)
Velo(2) = LorentzFacInv*PartState(5,PartID)
Velo(3) = LorentzFacInv*PartState(6,PartID)
Pt = E + CROSS(Velo,B)
#else
Pt(1) = E(1)
Pt(2) = E(2)
Pt(3) = E(3)
#endif

END SUBROUTINE PartRHS_RM


SUBROUTINE PartRHS_CEM(PartID,FieldAtParticle,Pt,LorentzFacInvIn)
!===================================================================================================================================
! 'constant-EM'
!
! A constant electromagnetic field E = (/Ex, Ey, Ez/) = const. and B = (/Bx, By, Bz/) = const. is used for all charged
! particles (simply uses the "externalField" variable)
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: abort
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies
USE MOD_Globals_Vars          ,ONLY: c2_inv,c2
#if USE_MPI
USE MOD_Globals               ,ONLY: MyRank
#endif
USE MOD_PICInterpolation_Vars ,ONLY: externalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: PartID
REAL,INTENT(IN)          :: FieldAtParticle(1:6)
REAL,INTENT(IN),OPTIONAL :: LorentzFacInvIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)         :: Pt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: velosq,F(1:3)
!===================================================================================================================================
ASSOCIATE (&
      qmt => Species(PartSpecies(PartID))%ChargeIC/Species(PartSpecies(PartID))%MassIC ,& ! charge/m_0
      vx  => PartState(4,PartID) ,& ! Velocity in x
      vy  => PartState(5,PartID) ,& ! Velocity in y
      vz  => PartState(6,PartID) ,& ! Velocity in z
      E1  => externalField(1)    ,& ! Electric field in x
      E2  => externalField(2)    ,& ! Electric field in y
      E3  => externalField(3)    ,& ! Electric field in z
      B1  => externalField(4)    ,& ! Magnetic field in x
      B2  => externalField(5)    ,& ! Magnetic field in y
      B3  => externalField(6)     & ! Magnetic field in z
      )

  ! Check squared velocity with c^2
  velosq = vx*vx + vy*vy + vz*vz
  IF(velosq.GT.c2) THEN
    IPWRITE(*,*) ' Particle is faster than the speed of light (v_x^2 + v_y^2 + v_z^2 > c^2)'
    IPWRITE(*,*) ' Species-ID',PartSpecies(PartID)
    IPWRITE(*,*) ' x=',PartState(1,PartID),' y=',PartState(2,PartID),' z=',PartState(3,PartID)
    CALL abort(&
        __STAMP__&
        ,'Particle is faster than the speed of light. Particle-Nr., velosq/c2:',PartID,velosq*c2_inv)
  END IF
  ASSOCIATE ( gammas => SQRT(1.0-velosq*c2_inv) ) ! Inverse of Lorentz factor

    F(1) = E1 + vy * B3 - vz * B2
    F(2) = E2 + vz * B1 - vx * B3
    F(3) = E3 + vx * B2 - vy * B1

    ! Calculate the acceleration
    Pt = gammas * qmt * ( F - DOT_PRODUCT(F,PartState(4:6,PartID))*PartState(4:6,PartID)*c2_inv )
  END ASSOCIATE
END ASSOCIATE

! Suppress compiler warning
RETURN
velosq=LorentzFacInvIn ! dummy statement
velosq=FieldAtParticle(1) ! dummy statement

END SUBROUTINE PartRHS_CEM


PPURE FUNCTION CalcPartRHSRotRefFrame(PosRotRef,VeloRotRef)
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: CROSS
USE MOD_Particle_Vars ,ONLY: RotRefFrameOmega
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: PosRotRef(1:3), VeloRotRef(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: CalcPartRHSRotRefFrame(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(ALL(ALMOSTZERO(VeloRotRef(1:3)))) THEN
  CalcPartRHSRotRefFrame(1:3) = - CROSS(RotRefFrameOmega(1:3),CROSS(RotRefFrameOmega(1:3),PosRotRef(1:3)))
ELSE
  CalcPartRHSRotRefFrame(1:3) = - CROSS(RotRefFrameOmega(1:3),CROSS(RotRefFrameOmega(1:3),PosRotRef(1:3))) &
                                - 2.*CROSS(RotRefFrameOmega(1:3),VeloRotRef(1:3))
END IF

END FUNCTION CalcPartRHSRotRefFrame


SUBROUTINE PartVeloToImp(VeloToImp,doParticle_In)
!===================================================================================================================================
! map the particle velocity to gamma*velocity
! or
! gamma*velocity to velocity
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: PDM, PartState, PartLorentzType
USE MOD_Globals_Vars  ,ONLY: c2_inv
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
      LorentzFac=1.0-DOTPRODUCT(PartState(4:6,iPart))*c2_inv
      LorentzFac=1.0/SQRT(LorentzFac)
      PartState(4:6,iPart) = LorentzFac*PartState(4:6,iPart)
    END IF ! DoParticle
  END DO ! iPart
ELSE
  DO iPart=1,PDM%ParticleVecLength
    IF(DoParticle(iPart))THEN
      LorentzFacInv=1.0+DOTPRODUCT(PartState(4:6,iPart))*c2_inv
      LorentzFacInv=1.0/SQRT(LorentzFacInv)
      PartState(4:6,iPart) = LorentzFacInv*PartState(4:6,iPart)
    END IF ! DoParticle
  END DO ! iPart
END IF

END SUBROUTINE PartVeloToImp

END MODULE MOD_part_RHS
