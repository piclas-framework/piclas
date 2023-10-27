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

MODULE MOD_Particle_Analyze_Code
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#if defined(PARTICLES) && defined(CODE_ANALYZE)
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: WriteParticleTrackingDataAnalytic
PUBLIC :: CalcAnalyticalParticleState
PUBLIC :: AnalyticParticleMovement
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Calculate the analytical position and velocity depending on the pre-defined function
!===================================================================================================================================
SUBROUTINE CalcAnalyticalParticleState(t,PartStateAnalytic,alpha_out,theta_out)
! MODULES
USE MOD_Preproc
USE MOD_Globals               ,ONLY: DOTPRODUCT,abort,CROSS
USE MOD_Globals_Vars          ,ONLY: PI,c,c2_inv,c2
USE MOD_PICInterpolation_Vars ,ONLY: AnalyticInterpolationType,AnalyticInterpolationSubType,AnalyticInterpolationP
USE MOD_PICInterpolation_Vars ,ONLY: AnalyticInterpolationPhase,AnalyticInterpolationGamma,AnalyticInterpolationE,AnalyticPartDim
USE MOD_PICInterpolation_Vars ,ONLY: TimeReset,r_WallVec,v_WallVec
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_PARTICLE_Vars         ,ONLY: PartSpecies,Species,RotRefFrameOmega,RotRefFrameFreq
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t                        !< simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: PartStateAnalytic(1:AnalyticPartDim)   !< analytic position and velocity
REAL,INTENT(OUT),OPTIONAL     :: alpha_out                              !< dimensionless parameter: alpha_out = q*B_0*l / (m*v_perpendicular)
REAL,INTENT(OUT),OPTIONAL     :: theta_out                              !< angle
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL    :: p
REAL    :: gamma_0
REAL    :: phi_0
REAL    :: Theta
REAL    :: beta
REAL    :: v_perp !< Perpendicular velocity
REAL    :: B_0    !< Magnetic flux density
REAL    :: c1
REAL    :: gamma1
REAL    :: TempArrayCross1(3),TempArrayCross2(3),TempArrayCross3(3),TempArrayCross4(3),TempArrayCross5(3),TempArrayCross6(3)
REAL    :: r_0Vec(3),v_0Vec(3),r_tVec(3),OmegaNormVec(3),omega
REAL    :: New_t
INTEGER :: iPart, iSpec
!===================================================================================================================================
PartStateAnalytic=0. ! default

iPart = 1
iSpec = PartSpecies(iPart)
! Select analytical solution depending on the type of the selected (analytic) interpolation
SELECT CASE(AnalyticInterpolationType)

! 0: const. magnetostatic field: B = B_z = (/ 0 , 0 , 1 T /) = const.
CASE(0)
  ! 0: non-relativistic, 1: relativistic
  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0) ! 0: non-relativistic
    ASSOCIATE(  B_0    => 1.0                                    ,& ! [T] cons. magnetic field
                v_perp => 1.0                                    ,& ! [m/s] perpendicular velocity (to guiding center)
                m      => Species(iSpec)%MassIC     ,& ! [kg] particle mass
                q      => Species(iSpec)%ChargeIC   ,& ! [C] particle charge
                phi    => AnalyticInterpolationPhase             )  ! [rad] phase shift
      ASSOCIATE( omega_c => ABS(q)*B_0/m )
        ASSOCIATE( r_c => v_perp/omega_c )
          PartStateAnalytic(1) = COS(omega_c*t + phi)*r_c
          PartStateAnalytic(2) = SIN(omega_c*t + phi)*r_c
          PartStateAnalytic(3) = 0.
          PartStateAnalytic(4) = -SIN(omega_c*t + phi)*v_perp
          PartStateAnalytic(5) =  COS(omega_c*t + phi)*v_perp
          PartStateAnalytic(6) = 0.
        END ASSOCIATE
      END ASSOCIATE
    END ASSOCIATE
  CASE(1) ! 1: relativistic
    ASSOCIATE(  gamma1 => AnalyticInterpolationGamma ,& ! Lorentz factor
                m      => 1.0                        ,& ! [kg] particle mass
                q      => 1.0                        ,& ! [C] particle charge
                phi    => AnalyticInterpolationPhase )  ! [rad] phase shift
      !-- get Lorentz factor gamma1(n)
      v_perp = c*SQRT(1.0 - 1/(gamma1**2))
      !v_perp = 2.99792457999850103770999962525942749e8
      IF(v_perp.GE.c) CALL abort(__STAMP__,'Velocity is geater than c',RealInfoOpt=v_perp)
      !-- Set const. magnetic field [T]
      B_0 = gamma1*v_perp
      !IPWRITE(UNIT_StdOut,*) "v_perp,v_perp/c,B_0 =", v_perp,v_perp/c,B_0

      !write(*,*) gamma1,v_perp,v_perp/c,B_0,B_0/c,1.0/SQRT(1.0-v_perp**2*c2_inv)
      ASSOCIATE( omega_c => ABS(q)*B_0/(gamma1*m) )
        !WRITE (*,*) "omega_c =", omega_c
        ASSOCIATE( r_c => v_perp/omega_c )
          PartStateAnalytic(1) = COS(omega_c*t + phi)*r_c
          PartStateAnalytic(2) = SIN(omega_c*t + phi)*r_c
          PartStateAnalytic(3) = 0.
          PartStateAnalytic(4) = -SIN(omega_c*t + phi)*v_perp
          PartStateAnalytic(5) =  COS(omega_c*t + phi)*v_perp
          PartStateAnalytic(6) = 0.
        END ASSOCIATE
      END ASSOCIATE
    END ASSOCIATE
  END SELECT

! 1: magnetostatic field: B = B_z = (/ 0 , 0 , B_0 * EXP(x/l) /) = const.
CASE(1)
  SELECT CASE(AnalyticInterpolationSubType)
  CASE(1,2)
    ASSOCIATE(  p       => AnalyticInterpolationP , &
                Theta_0 => -PI/2.0                     , &
                t       => t - TEnd/2. )
                !t       => t )
      ! gamma
      gamma_0 = SQRT(ABS(p*p - 1.))

      ! angle
      Theta   = -2.*ATAN( SQRT((1.+p)/(1.-p)) * TANH(0.5*gamma_0*t) ) + Theta_0

      ! x-pos
      PartStateAnalytic(1) = LOG(-SIN(Theta) + p )

      ! y-pos
      PartStateAnalytic(2) = p*t + Theta - Theta_0
    END ASSOCIATE
  CASE(3)
    ASSOCIATE(  p       => AnalyticInterpolationP , &
                Theta_0 => -PI/2.0                      &
                )
      ! gamma
      gamma_0 = SQRT(ABS(p*p - 1.))

      ! angle
      Theta   = -2.*ATAN( SQRT((p+1.)/(p-1.)) * TAN(0.5*gamma_0*t) ) -2.*PI*REAL(NINT((gamma_0*t)/(2.*PI))) + Theta_0

      ! x-pos
      PartStateAnalytic(1) = LOG(-SIN(Theta) + p )

      ! y-pos
      PartStateAnalytic(2) = p*t + Theta - Theta_0
    END ASSOCIATE
  CASE(11,21) ! old version of CASE(1,2)
    ASSOCIATE( p       => AnalyticInterpolationP , &
          Theta_0 => 0.d0 ) !0.785398163397448d0    )
      beta = ACOS(p)
      !beta = ASIN(-p)
      ! phase shift
      phi_0   = ATANH( (1./TAN(beta/2.)) * TAN(Theta_0/2.) )
      ! angle
      Theta   = -2.*ATANH( TAN(beta/2.) * TANH(0.5*t*SIN(beta)-phi_0) )
      Theta   = -2.*ATANH( TAN(beta/2.) * TANH(0.5*SIN(beta*t)-phi_0) )
      ! x-pos
      PartStateAnalytic(1) = LOG((COS(Theta)-p)/(COS(Theta_0)-p))
      ! y-pos
      PartStateAnalytic(2) = p*t - (Theta-Theta_0)
    END ASSOCIATE
  CASE(31) ! old version of CASE(3)
    ASSOCIATE( p       => AnalyticInterpolationP , &
                Theta_0 => 0.d0                   )
      gamma_0 = SQRT(p*p-1.)
      ! phase shift
      phi_0   = ATAN( (gamma_0/(p-1.)) * TAN(Theta_0/2.) )
      ! angle
      Theta   = 2.*ATAN( SQRT((p-1)/(p+1)) * TAN(0.5*gamma_0*t - phi_0) ) + 2*Pi*REAL(NINT((t*gamma_0)/(2*Pi) - phi_0/Pi))
      ! x-pos
      PartStateAnalytic(1) = LOG((COS(Theta)-p)/(COS(Theta_0)-p))
      ! y-pos
      PartStateAnalytic(2) = p*t - (Theta-Theta_0)
    END ASSOCIATE
  END SELECT

  SELECT CASE(AnalyticInterpolationSubType)
  CASE(1,2,3)
    ! Set analytic velocity
    PartStateAnalytic(4) = COS(Theta)
    PartStateAnalytic(5) = SIN(Theta)
    PartStateAnalytic(6) = 0.
  CASE(11,21,31)
    ! Set analytic velocity
    PartStateAnalytic(4) = SIN(Theta)
    PartStateAnalytic(5) = COS(Theta)
    PartStateAnalytic(6) = 0.
  END SELECT

  ! Optional output variables
  IF(PRESENT(alpha_out))THEN
    ASSOCIATE( dot_theta => SIN(Theta) - AnalyticInterpolationP )
      ASSOCIATE( alpha_0 => -dot_theta / EXP(PartStateAnalytic(1)) )
        alpha_out = alpha_0
        WRITE (*,*) "alpha_out =", alpha_out
      END ASSOCIATE
    END ASSOCIATE
  END IF
  IF(PRESENT(theta_out))THEN
    theta_out = Theta
    WRITE (*,*) "theta_out =", theta_out
  END IF

! 2: const. electromagnetic field: B = B_z = (/ 0 , 0 , (x^2+y^2)^0.5 /) = const.
!                                  E = 1e-2/(x^2+y^2)^(3/2) * (/ x , y , 0. /)
CASE(2)
  ! missing ...

! 3: const. electric field: E = E_x = (/ 1 V/m , 0 , 0 /) = const.
CASE(3)

  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0) ! 0: non-relativistic
    ASSOCIATE(  m      => Species(iSpec)%MassIC     ,& ! [kg] particle mass
                q      => Species(iSpec)%ChargeIC   ,& ! [C] particle charge
                E      => AnalyticInterpolationE                 )  ! [V/m] Electric field strength in x-direction

      ASSOCIATE( a => q*E/m )
        PartStateAnalytic(1) = 0.5*a*t*t
        PartStateAnalytic(2) = 0.
        PartStateAnalytic(3) = 0.
        PartStateAnalytic(4) = a*t
        PartStateAnalytic(5) = 0.
        PartStateAnalytic(6) = 0.
      END ASSOCIATE
    END ASSOCIATE
  CASE(1) ! 1: relativistic
    ASSOCIATE(  m      => Species(iSpec)%MassIC     ,& ! [kg] particle mass
                q      => Species(iSpec)%ChargeIC   ,& ! [C] particle charge
                E      => AnalyticInterpolationE                 )  ! [V/m] Electric field strength in x-direction
      ASSOCIATE( aStar => q*E/m )
        ASSOCIATE( tStar => (c/aStar)*ASINH(aStar*t/c) ,&
                    b     => c2/aStar )
        ASSOCIATE( c1 => aStar*tStar/c )
          ASSOCIATE( gamma1 => COSH(c1)&
                      )
            PartStateAnalytic(1) = b*gamma1 -b
            !WRITE (*,*) "q,E,m,t =", q,E,m,t
        !PartStateAnalytic(1) = 0.5*q*E/m*t*t
            PartStateAnalytic(2) = 0.
            PartStateAnalytic(3) = 0.
            PartStateAnalytic(4) = c*TANH(c1)
            PartStateAnalytic(5) = 0.
            PartStateAnalytic(6) = 0.
            !WRITE (*,*) "t,PartStateAnalytic =", t,PartStateAnalytic,gamma1
          END ASSOCIATE
          END ASSOCIATE
        END ASSOCIATE
      END ASSOCIATE
    END ASSOCIATE
  END SELECT

! 4: const. electric field: E = E_x = (/ X V/m , 0 , 0 /) = const.
CASE(4)

  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0) ! 0: non-relativistic
  CASE(1) ! 1: relativistic
    ASSOCIATE(  m      => Species(iSpec)%MassIC     ,& ! [kg] particle mass
                q      => Species(iSpec)%ChargeIC   ,& ! [C] particle charge
                E      => AnalyticInterpolationE                 )  ! [V/m] Electric field strength in x-direction
            c1 = q*E/m
            gamma_0 = SQRT(1.0+(c1*t/c)**2)
            PartStateAnalytic(1) = c**2/c1*(gamma_0-1.0)
            PartStateAnalytic(2) = 0.
            PartStateAnalytic(3) = 0.
            PartStateAnalytic(4) = c1*t/gamma_0
            PartStateAnalytic(5) = 0.
            PartStateAnalytic(6) = 0.
            !WRITE (*,*) "x,v,gamma_0,c =", PartStateAnalytic(1),PartStateAnalytic(4),gamma_0,c
    END ASSOCIATE
  END SELECT
! 5: motion of the particle in RotRefFrame without Collisions
CASE(5)
  r_0Vec          = Species(iSpec)%Init(1)%BasePointIC(1:3)
  v_0Vec          = Species(iSpec)%Init(1)%VeloVecIC(1:3) * Species(iSpec)%Init(1)%VeloIC
  v_0Vec          = v_0Vec - CROSS(RotRefFrameOmega(1:3),r_0Vec)
  r_tVec          = r_0Vec + v_0Vec * t
  omega           = 2.*PI*RotRefFrameFreq
  OmegaNormVec    = RotRefFrameOmega/omega
  TempArrayCross1 = CROSS(r_tVec,OmegaNormVec)
  TempArrayCross2 = CROSS(OmegaNormVec,TempArrayCross1)
  TempArrayCross3 = CROSS(r_0Vec,OmegaNormVec)
  TempArrayCross4 = CROSS(OmegaNormVec,TempArrayCross3)
  TempArrayCross5 = CROSS(v_0Vec,OmegaNormVec)
  TempArrayCross6 = CROSS(OmegaNormVec,TempArrayCross5)

  PartStateAnalytic(1) = r_0Vec(1) + v_0Vec(1) * t                                                          &
                        - TempArrayCross2(1)                                                                 &
                        + SIN(omega * t) * TempArrayCross3(1) + COS(omega * t) * TempArrayCross4(1)          &
                        + omega * t * SIN(omega * t) * ( TempArrayCross4(1) + 1/omega * TempArrayCross5(1) ) &
                        - omega * t * COS(omega * t) * ( TempArrayCross3(1) - 1/omega * TempArrayCross6(1) )

  PartStateAnalytic(2) = r_0Vec(2) + v_0Vec(2) * t                                                          &
                        - TempArrayCross2(2)                                                                 &
                        + SIN(omega * t) * TempArrayCross3(2) + COS(omega * t) * TempArrayCross4(2)          &
                        + omega * t * SIN(omega * t) * ( TempArrayCross4(2) + 1/omega * TempArrayCross5(2) ) &
                        - omega * t * COS(omega * t) * ( TempArrayCross3(2) - 1/omega * TempArrayCross6(2) )

  PartStateAnalytic(3) = r_0Vec(3) + v_0Vec(3) * t                                                          &
                        - TempArrayCross2(3)                                                                 &
                        + SIN(omega * t) * TempArrayCross3(3) + COS(omega * t) * TempArrayCross4(3)          &
                        + omega * t * SIN(omega * t) * ( TempArrayCross4(3) + 1/omega * TempArrayCross5(3) ) &
                        - omega * t * COS(omega * t) * ( TempArrayCross3(3) - 1/omega * TempArrayCross6(3) )

  PartStateAnalytic(4) = v_0Vec(1)           &
                        - TempArrayCross6(1)  &
                        + omega * ( COS(omega * t) * TempArrayCross3(1) - SIN(omega * t) * TempArrayCross4(1) ) &
                        + ( omega * (SIN(omega * t) + omega * t * COS(omega * t) ) ) &
                                                        * ( TempArrayCross4(1) + 1/omega * TempArrayCross5(1) ) &
                        - ( omega * (COS(omega * t) - omega * t * SIN(omega * t) ) ) &
                                                        * ( TempArrayCross3(1) - 1/omega * TempArrayCross6(1) )

  PartStateAnalytic(5) = v_0Vec(2)           &
                        - TempArrayCross6(2)  &
                        + omega * ( COS(omega * t) * TempArrayCross3(2) - SIN(omega * t) * TempArrayCross4(2) ) &
                        + ( omega * (SIN(omega * t) + omega * t * COS(omega * t) ) ) &
                                                        * ( TempArrayCross4(2) + 1/omega * TempArrayCross5(2) ) &
                        - ( omega * (COS(omega * t) - omega * t * SIN(omega * t) ) ) &
                                                        * ( TempArrayCross3(2) - 1/omega * TempArrayCross6(2) )

  PartStateAnalytic(6) = v_0Vec(3)           &
                        - TempArrayCross6(3)  &
                        + omega * ( COS(omega * t) * TempArrayCross3(3) - SIN(omega * t) * TempArrayCross4(3) ) &
                        + ( omega * (SIN(omega * t) + omega * t * COS(omega * t) ) ) &
                                                        * ( TempArrayCross4(3) + 1/omega * TempArrayCross5(3) ) &
                        - ( omega * (COS(omega * t) - omega * t * SIN(omega * t) ) ) &
                                                        * ( TempArrayCross3(3) - 1/omega * TempArrayCross6(3) )
! 51: motion of the particle in RotRefFrame with wall collisions at x=0
CASE(51)
  IF(TimeReset.GT.0.0) THEN
    r_0Vec          = r_WallVec
    v_0Vec          = v_WallVec
    New_t           = t - TimeReset
  ELSE
    r_0Vec          = Species(iSpec)%Init(1)%BasePointIC(1:3)
    v_0Vec          = Species(iSpec)%Init(1)%VeloVecIC(1:3) * Species(iSpec)%Init(1)%VeloIC
    New_t           = t
    v_0Vec = v_0Vec - CROSS(RotRefFrameOmega(1:3),r_0Vec)
  END IF
  r_tVec          = r_0Vec + v_0Vec * New_t
  omega = 2.*PI*RotRefFrameFreq
  OmegaNormVec    = RotRefFrameOmega/omega
  TempArrayCross1 = CROSS(r_tVec,OmegaNormVec)
  TempArrayCross2 = CROSS(OmegaNormVec,TempArrayCross1)
  TempArrayCross3 = CROSS(r_0Vec,OmegaNormVec)
  TempArrayCross4 = CROSS(OmegaNormVec,TempArrayCross3)
  TempArrayCross5 = CROSS(v_0Vec,OmegaNormVec)
  TempArrayCross6 = CROSS(OmegaNormVec,TempArrayCross5)

  PartStateAnalytic(1) = r_0Vec(1) + v_0Vec(1) * New_t                                                          &
                        - TempArrayCross2(1)                                                                 &
                        + SIN(omega * New_t) * TempArrayCross3(1) + COS(omega * New_t) * TempArrayCross4(1)          &
                        + omega * New_t * SIN(omega * New_t) * ( TempArrayCross4(1) + 1/omega * TempArrayCross5(1) ) &
                        - omega * New_t * COS(omega * New_t) * ( TempArrayCross3(1) - 1/omega * TempArrayCross6(1) )

  PartStateAnalytic(2) = r_0Vec(2) + v_0Vec(2) * New_t                                                          &
                        - TempArrayCross2(2)                                                                 &
                        + SIN(omega * New_t) * TempArrayCross3(2) + COS(omega * New_t) * TempArrayCross4(2)          &
                        + omega * New_t * SIN(omega * New_t) * ( TempArrayCross4(2) + 1/omega * TempArrayCross5(2) ) &
                        - omega * New_t * COS(omega * New_t) * ( TempArrayCross3(2) - 1/omega * TempArrayCross6(2) )

  PartStateAnalytic(3) = r_0Vec(3) + v_0Vec(3) * New_t                                                          &
                        - TempArrayCross2(3)                                                                 &
                        + SIN(omega * New_t) * TempArrayCross3(3) + COS(omega * New_t) * TempArrayCross4(3)          &
                        + omega * New_t * SIN(omega * New_t) * ( TempArrayCross4(3) + 1/omega * TempArrayCross5(3) ) &
                        - omega * New_t * COS(omega * New_t) * ( TempArrayCross3(3) - 1/omega * TempArrayCross6(3) )


  PartStateAnalytic(4) = v_0Vec(1)           &
                        - TempArrayCross6(1)  &
                        + omega * ( COS(omega * New_t) * TempArrayCross3(1) - SIN(omega * New_t) * TempArrayCross4(1) ) &
                        + ( omega * (SIN(omega * New_t) + omega * New_t * COS(omega * New_t) ) ) &
                                                        * ( TempArrayCross4(1) + 1/omega * TempArrayCross5(1) ) &
                        - ( omega * (COS(omega * New_t) - omega * New_t * SIN(omega * New_t) ) ) &
                                                        * ( TempArrayCross3(1) - 1/omega * TempArrayCross6(1) )

  PartStateAnalytic(5) = v_0Vec(2)           &
                        - TempArrayCross6(2)  &
                        + omega * ( COS(omega * New_t) * TempArrayCross3(2) - SIN(omega * New_t) * TempArrayCross4(2) ) &
                        + ( omega * (SIN(omega * New_t) + omega * New_t * COS(omega * New_t) ) ) &
                                                        * ( TempArrayCross4(2) + 1/omega * TempArrayCross5(2) ) &
                        - ( omega * (COS(omega * New_t) - omega * New_t * SIN(omega * New_t) ) ) &
                                                        * ( TempArrayCross3(2) - 1/omega * TempArrayCross6(2) )

  PartStateAnalytic(6) = v_0Vec(3)           &
                        - TempArrayCross6(3)  &
                        + omega * ( COS(omega * New_t) * TempArrayCross3(3) - SIN(omega * New_t) * TempArrayCross4(3) ) &
                        + ( omega * (SIN(omega * New_t) + omega * New_t * COS(omega * New_t) ) ) &
                                                        * ( TempArrayCross4(3) + 1/omega * TempArrayCross5(3) ) &
                        - ( omega * (COS(omega * New_t) - omega * New_t * SIN(omega * New_t) ) ) &
                                                        * ( TempArrayCross3(3) - 1/omega * TempArrayCross6(3) )

!              IF((PartStateAnalytic(1).GE.0.0).AND.(TimeReset.LE.0.0)) THEN
  IF(ABS(PartStateAnalytic(1)).LT.1E-8) THEN
    TimeReset    = t
    r_WallVec    = PartStateAnalytic(1:3)
    v_WallVec(1) = -PartStateAnalytic(4)
    v_WallVec(2) = PartStateAnalytic(5)
    v_WallVec(3) = PartStateAnalytic(6)
  END IF
!              PartStateAnalytic(4:6) = 0.0
END SELECT

! Calculate analytical Lorentz factor
gamma1 = DOTPRODUCT(PartStateAnalytic(4:6))*c2_inv
! Sanity check: Lorentz factor must be below 1.0
IF(gamma1.GE.1.0)THEN
  PartStateAnalytic(7)=-1.0
ELSE
  PartStateAnalytic(7)=1.0/SQRT(1.-gamma1)
END IF

END SUBROUTINE CalcAnalyticalParticleState


!===================================================================================================================================
!> Calculates "running" L_2 norms
!> running means: use the old L_2 error from the previous iteration in order to determine the L_2 error over time (simulation time)
!>
!> -------------------------------------------------------------------------
!> OLD METHOD: assuming constant timestep (ignoring the total time tEnd -> Delta t = tEnd / Niter)
!> L_2(t) = SQRT( ( L_2(t-1)^2 * (iter-1) + delta(t)^2 ) / iter )
!>
!> -------------------------------------------------------------------------
!> NEW METHOD: assuming variable timestep
!> L_2(t) = SQRT(  L_2(t-1)^2   +   (t - t_old) * delta(t)^2  )
!>
!> t     : simulation time
!> t_old : simulation time of the last iteration
!> L_2   : error norm
!> delta : difference numerical to analytical solution
!> iter  : simulation iteration counter
!===================================================================================================================================
SUBROUTINE CalcErrorParticle(t,iter,PartStateAnalytic)
! MODULES
USE MOD_PICInterpolation_Vars ,ONLY: L_2_Error_Part,AnalyticPartDim,L_2_Error_Part_time
USE MOD_Particle_Vars         ,ONLY: PartState, PDM
USE MOD_Globals_Vars          ,ONLY: c2_inv
USE MOD_globals               ,ONLY: DOTPRODUCT
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Particle_Vars         ,ONLY: Pt
USE MOD_part_RHS              ,ONLY: CalcPartRHSSingleParticle
!USE MOD_PICInterpolation      ,ONLY: InterpolateFieldToParticle ! already known from previous call (causes circular definition)
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
#if (PP_TimeDiscMethod==4)
USE MOD_Particle_Vars         ,ONLY: PartVeloRotRef, UseRotRefFrame
#endif /*(PP_TimeDiscMethod==4)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)    :: iter                     !< simulation iteration counter
REAL,INTENT(IN)               :: t                        !< simulation time
REAL,INTENT(INOUT)            :: PartStateAnalytic(1:AnalyticPartDim)   !< analytic position and velocity
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iPartState
REAL                          :: PartStateLoc(1:AnalyticPartDim)
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
REAL                          :: PartStateLocAnalytic(1:AnalyticPartDim)
INTEGER,PARAMETER             :: Method=2 !< 1: shift numerical solution from v(n-1/2) to v(n), gives O(2) for x and v for Boris-LF
                                          !<    for const. magnetic field
                                          !< 2: shift analytical solution from v(n) to v(n-1/2), gives no order of convergence for
                                          !<    the velocity components Boris-LF and const. magnetic field due to conservation of
                                          !<    energy
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
REAL                          :: gamma1
!===================================================================================================================================
! Get analytic particle position
CALL CalcAnalyticalParticleState(t,PartStateAnalytic)

! Depending on the iteration counter, set the L_2 error (re-use the value in the next loop)
IF(iter.LT.1)THEN ! first iteration
  L_2_Error_Part(1:AnalyticPartDim) = 0.
  L_2_Error_Part_time = 0.
ELSE
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN

      ! Store particle state in temp. variable
      PartStateLoc(1:6) = PartState(1:6,iPart)

      !-- Only for DSMC timedisk method for testing rotational reference
      ! Use PartVeloRotRef for comparison with analytical solution
#if (PP_TimeDiscMethod==4)
      IF(UseRotRefFrame) THEN
        PartStateLoc(4:6) = PartVeloRotRef(1:3,iPart)
      END IF
#endif /*(PP_TimeDiscMethod==4)*/
      !-- Only for time-staggered methods (Leapfrog and Boris-Leapfrog):
      ! Set analytic velocity at v(n-0.5) from analytic particle solution
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
      IF(Method.EQ.1)THEN
        !CALL InterpolateFieldToParticle()   ! forces on particles, already known from previous call (causes circular definition)
        CALL CalcPartRHSSingleParticle(iPart)

        !-- v(n+0.5) => v(n+1) by a(n+1):
        PartStateLoc(4:6) = PartState(4:6,iPart) + Pt(1:3,iPart) * dt*0.5
      ELSE
        CALL CalcAnalyticalParticleState(t-dt*0.5,PartStateLocAnalytic)
        !-- get analytical v(n-0.5)
        PartStateAnalytic(4:7) = PartStateLocAnalytic(4:7)
      END IF ! Method.EQ.1
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/

        ! Calculate new Lorentz factor
        gamma1 = DOTPRODUCT(PartStateLoc(4:6))*c2_inv
        ! Sanity check: Lorentz factor must be below 1.0
        IF(gamma1.GE.1.0)THEN
          PartStateLoc(7)=-1.0
        ELSE
          PartStateLoc(7)=1.0/SQRT(1.-gamma1)
        END IF

      DO iPartState = 1, AnalyticPartDim
        ! OLD METHOD: original
        ! L_2_Error_Part(iPartState) = SQRT( ( (L_2_Error_Part(iPartState))**2*REAL(iter-1) + &
        !                               (PartStateAnalytic(iPartState)-PartStateLoc(iPartState))**2 )/ REAL(iter))

        ! OLD METHOD: considering TEnd
        ! L_2_Error_Part(iPartState) = SQRT( Tend * ( (L_2_Error_Part(iPartState))**2*REAL(iter-1) + &
        !                               (PartStateAnalytic(iPartState)-PartStateLoc(iPartState))**2 ) &
        !                      / REAL(iter))

        ! NEW METHOD: considering variable time step
        L_2_Error_Part(iPartState) = SQRT(  (L_2_Error_Part(iPartState))**2 + &
                                   (t-L_2_Error_Part_time)*(PartStateAnalytic(iPartState)-PartStateLoc(iPartState))**2 )

        ! Additional method: Consider only the last difference
        !L_2_Error_Part(iPartState) = PartStateAnalytic(iPartState)-PartStateLoc(iPartState)
      END DO ! iPartState = 1, 6
      !WRITE (*,*) "PartStateAnalytic =", PartStateAnalytic
      !WRITE (*,*) "PartStateLoc      =", PartStateLoc
      !WRITE (*,*) "L_2_Error_Part    =", L_2_Error_Part
      L_2_Error_Part_time = t
    ELSE
      L_2_Error_Part(1:AnalyticPartDim) = -1.0
    END IF
  END DO
END IF

END SUBROUTINE CalcErrorParticle


!===================================================================================================================================
!> Calculate the analytical position and velocity depending on the pre-defined function
!===================================================================================================================================
SUBROUTINE AnalyticParticleMovement(time,iter)
! MODULES
USE MOD_Preproc
USE MOD_Globals               ,ONLY: UNIT_StdOut,MPIRoot
USE MOD_Analyze_Vars          ,ONLY: OutputErrorNorms
USE MOD_Particle_Analyze_Vars ,ONLY: TrackParticlePosition
USE MOD_PICInterpolation_Vars ,ONLY: L_2_Error_Part,AnalyticPartDim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: time                        !< simulation time
INTEGER(KIND=8),INTENT(IN)    :: iter                        !< iteration
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: PartStateAnalytic(1:AnalyticPartDim)   !< analytic position and velocity
CHARACTER(LEN=40)             :: formatStr
!===================================================================================================================================

CALL CalcErrorParticle(time,iter,PartStateAnalytic)
IF(MPIRoot.AND.OutputErrorNorms) THEN
  WRITE(UNIT_StdOut,'(A13,ES16.7)')' Sim time  : ',time
  WRITE(formatStr,'(A5,I1,A7)')'(A13,',AnalyticPartDim,'ES16.7)'
  WRITE(UNIT_StdOut,formatStr)' L2_Part   : ',L_2_Error_Part
  OutputErrorNorms=.FALSE.
END IF
IF(TrackParticlePosition) CALL WriteParticleTrackingDataAnalytic(time,iter,PartStateAnalytic) ! new function

END SUBROUTINE AnalyticParticleMovement


!----------------------------------------------------------------------------------------------------------------------------------!
!> Write analytic particle info to ParticlePositionAnalytic.csv file
!> time, pos, velocity
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteParticleTrackingDataAnalytic(time,iter,PartStateAnalytic)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: MPIRoot,FILEEXISTS,unit_stdout,DOTPRODUCT
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_PICInterpolation_Vars ,ONLY: L_2_Error_Part,AnalyticPartDim
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: time
INTEGER(KIND=8),INTENT(IN)       :: iter
REAL(KIND=8),INTENT(IN)          :: PartStateAnalytic(1:AnalyticPartDim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=28),PARAMETER              :: outfile='ParticlePositionAnalytic.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=15
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    '001-time',     &
    'PartPosX_Analytic', &
    'PartPosY_Analytic', &
    'PartPosZ_Analytic', &
    'PartVelX_Analytic', &
    'PartVelY_Analytic', &
    'PartVelZ_Analytic', &
    'gamma'            , &
    'L2_PartPosX'      , &
    'L2_PartPosY'      , &
    'L2_PartPosZ'      , &
    'L2_PartVelX'      , &
    'L2_PartVelY'      , &
    'L2_PartVelZ'      , &
    'L2_gamma'           &
    /)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
LOGICAL                                  :: FileExist,CreateFile
!===================================================================================================================================
! only the root shall write this file
IF(.NOT.MPIRoot)RETURN

! check if file is to be created
CreateFile=.TRUE.
IF(iter.GT.0)CreateFile=.FALSE.                             ! don't create new file if this is not the first iteration
IF((DoRestart).AND.(FILEEXISTS(outfile)))CreateFile=.FALSE. ! don't create new file if this is a restart and the file already exists
!                                                           ! assume continued simulation and old load balance data is still needed

! check if new file with header is to be created
INQUIRE(FILE = outfile, EXIST=FileExist)
IF(.NOT.FileExist)CreateFile=.TRUE.                         ! if no file exists, create one

! create file with header
IF(CreateFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible relimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF

! Print info to file
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nOutputVar,CSVFORMAT,')'
  WRITE(tmpStr2,formatStr)&
      " ",time, &                           ! time
      delimiter,PartStateAnalytic(1), &     ! PartPosX analytic solution
      delimiter,PartStateAnalytic(2), &     ! PartPosY analytic solution
      delimiter,PartStateAnalytic(3), &     ! PartPosZ analytic solution
      delimiter,PartStateAnalytic(4), &     ! PartVelX analytic solution
      delimiter,PartStateAnalytic(5), &     ! PartVelY analytic solution
      delimiter,PartStateAnalytic(6), &     ! PartVelZ analytic solution
      delimiter,PartStateAnalytic(7), &     ! Lorentz factor
      delimiter,L_2_Error_Part(1), &     ! L2 error for PartPosX solution
      delimiter,L_2_Error_Part(2), &     ! L2 error for PartPosY solution
      delimiter,L_2_Error_Part(3), &     ! L2 error for PartPosZ solution
      delimiter,L_2_Error_Part(4), &     ! L2 error for PartVelX solution
      delimiter,L_2_Error_Part(5), &     ! L2 error for PartVelY solution
      delimiter,L_2_Error_Part(6), &     ! L2 error for PartVelZ solution
      delimiter,L_2_Error_Part(7)        ! L2 error for PartVelZ solution
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write particle tracking (analytic) info!"
END IF

END SUBROUTINE WriteParticleTrackingDataAnalytic


#endif /*defined(PARTICLES) && defined(CODE_ANALYZE)*/
END MODULE MOD_Particle_Analyze_Code
