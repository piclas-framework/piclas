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

MODULE MOD_Equation
!===================================================================================================================================
! Add comments please!
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
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE
INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE
INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE
INTERFACE DivCleaningDamping
  MODULE PROCEDURE DivCleaningDamping
END INTERFACE

INTERFACE VolInt_Pois
  MODULE PROCEDURE VolInt_weakForm
END INTERFACE

INTERFACE FillFlux_Pois
  MODULE PROCEDURE FillFlux
END INTERFACE

INTERFACE ProlongToFace_Pois
  MODULE PROCEDURE ProlongToFace_sideBased
END INTERFACE

INTERFACE SurfInt_Pois
  MODULE PROCEDURE SurfInt2
END INTERFACE

#if USE_MPI
PUBLIC::StartExchangeMPIData_Pois
#endif
PUBLIC::VolInt_Pois,FillFlux_Pois, ProlongToFace_Pois, SurfInt_Pois
PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation,DivCleaningDamping,EvalGradient,CalcSource_Pois,DivCleaningDamping_Pois
!===================================================================================================================================
PUBLIC::DefineParametersEquation
CONTAINS

!==================================================================================================================================
!> Define parameters for equation
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")

CALL prms%CreateRealOption(     'c_corr'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'multiplied with c0 results in the velocity of '//&
                                                     'introduced artificial correcting waves (HDC)' , '1.')
CALL prms%CreateRealOption(     'c0'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Velocity of light (in vacuum)' , '1.')
CALL prms%CreateRealOption(     'eps'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Electric constant (vacuum permittivity)' , '1.')
CALL prms%CreateRealOption(     'mu'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Magnetic constant (vacuum permeability = 4πE−7H/m)' &
                                                   , '1.')
CALL prms%CreateRealOption(     'fDamping'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Apply the damping factor also to PML source terms\n'//&
                                                     'but only to PML variables for Phi_E and Phi_B to prevent charge-related\n'//&
                                                     'instabilities (accumulation of divergence compensation over \n'//&
                                                     'timeU2 = U2 * fDamping' , '0.999')
CALL prms%CreateRealOption(     'fDamping_pois'    , 'TODO-DEFINE-PARAMETER' , '0.99')
CALL prms%CreateLogicalOption(  'ParabolicDamping' , 'TODO-DEFINE-PARAMETER' , '.FALSE.')
CALL prms%CreateLogicalOption(  'CentralFlux'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Flag for central or upwind flux' , '.FALSE.')
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection')

CALL prms%CreateRealArrayOption('xDipole'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Base point of electromagnetic dipole', '0. , 0. , 0.')

CALL prms%CreateIntOption(      'AlphaShape'       , 'TODO-DEFINE-PARAMETER', '2')
CALL prms%CreateRealOption(     'r_cutoff'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Modified for curved and shape-function influence'//&
                                                     ' (c*dt*SafetyFactor+r_cutoff)' , '1.0')

END SUBROUTINE DefineParametersEquation

SUBROUTINE InitEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Mesh_Vars
USE MOD_ReadInTools
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation_Vars ,ONLY: xGP
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Equation_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
#endif
!===================================================================================================================================
IF(InterpolationInitIsDone.AND.EquationInitIsDone)THEN
   LBWRITE(*,*) "InitMaxwell Poisson not ready to be called or already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT MAXWELL POISSON...'

! Read correction velocity
c_corr             = GETREAL('c_corr','1.')
c                  = GETREAL('c0','1.')
eps0               = GETREAL('eps','1.')
mu0                = GETREAL('mu','1.')
smu0               = 1./mu0
fDamping           = GETREAL('fDamping','0.99')
fDamping_pois      = GETREAL('fDamping_pois','0.99')
DoParabolicDamping = GETLOGICAL('ParabolicDamping','.FALSE.')
xDipole(1:3)       = GETREALARRAY('xDipole',3,'0.,0.,0.') ! dipole base point for CASE(4)

c_corr_c  = c_corr*c
c_corr_c2 = c_corr*c2
eta_c     = (c_corr-1.)*c

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
!WRITE(DefBCState,'(I3,A,I3,A,I3,A,I3,A,I3,A,I3)') &
!  IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc
!IF(BCType_in(1) .EQ. -999)THEN
!  BCType = GETINTARRAY('BoundaryType',6)
!ELSE
!  BCType=BCType_in
!  SWRITE(UNIT_stdOut,*)'|                   BoundaryType | -> Already read in CreateMPICart!'

!END IF
!BCState   = GETINTARRAY('BoundaryState',6,TRIM(DefBCState))
!BoundaryCondition(:,1) = BCType
!BoundaryCondition(:,2) = BCState
! Read exponent for shape function
alpha_shape = GETINT('AlphaShape','2')
rCutoff     = GETREAL('r_cutoff','1.')
! Compute factor for shape function
ShapeFuncPrefix = 1./(2. * beta(1.5, REAL(alpha_shape) + 1.) * REAL(alpha_shape) + 2. * beta(1.5, REAL(alpha_shape) + 1.)) &
                * (REAL(alpha_shape) + 1.)/(PI*(rCutoff**3))

!Init PHI
SDEALLOCATE(Phi)
ALLOCATE(Phi(4,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Phi=0.
! the time derivative computed with the DG scheme

SDEALLOCATE(Phit)
ALLOCATE(Phit(4,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
nTotalPhi=4*(PP_N+1)*(PP_N+1)*(PP_N+1)*PP_nElems

!IF(.NOT.DoRestart)THEN
!  ! U is filled with the ini solution
!  CALL FillIni()
!END IF
! Ut is set to zero because it is successively updated with DG contributions
Phit=0.

! We store the interior data at the each element face
SDEALLOCATE(Phi_master)
ALLOCATE(Phi_master(4,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))
SDEALLOCATE(Phi_slave)
ALLOCATE(Phi_slave(4,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper))
Phi_master=0.
Phi_slave=0.

! unique flux per side
SDEALLOCATE(FluxPhi)
ALLOCATE(FluxPhi(4,0:PP_N,0:PP_N,1:nSides))
FluxPhi=0.

!ElectricField as grad Phi
SDEALLOCATE(E)
ALLOCATE(E(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
E=0.
SDEALLOCATE(D)
ALLOCATE(D(0:PP_N,0:PP_N))
D=0.
CALL PolynomialDerivativeMatrix(N,xGP,D)

EquationInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT MAXWELL POISSON DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



SUBROUTINE ExactFunc(ExactFunction,t,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: c,c2,eps0
USE MOD_Equation_Vars ,ONLY: xDipole
USE MOD_Globals_Vars  ,ONLY: PI
# if (PP_TimeDiscMethod==1)
USE MOD_TimeDisc_vars ,ONLY: dt
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! state in conservative variables
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3),r,r2,zlen
REAL                            :: a, b, d, l, m, n, B0            ! aux. Variables for Resonator-Example
REAL                            :: gamma,Psi,GradPsiX,GradPsiY     !     -"-
REAL                            :: xrel(3), theta, Etheta          ! aux. Variables for Dipole
REAL,PARAMETER                  :: Q=1, dD=1, omegaD=6.28318E8     ! aux. Constants for Dipole
REAL                            :: c1,s1,b1,b2                     ! aux. Variables for Gyrotron
REAL                            :: eps,phi,z                       ! aux. Variables for Gyrotron
REAL                            :: Er,Br,Ephi,Bphi,Bz              ! aux. Variables for Gyrotron
REAL, PARAMETER                 :: B0G=1.0,g=3236.706462           ! aux. Constants for Gyrotron
REAL, PARAMETER                 :: k0=3562.936537,h=1489.378411    ! aux. Constants for Gyrotron
REAL, PARAMETER                 :: omegaG=3.562936537e+3           ! aux. Constants for Gyrotron
INTEGER, PARAMETER              :: mG=34,nG=19                     ! aux. Constants for Gyrotron
!===================================================================================================================================
Cent=x
SELECT CASE (ExactFunction)
#ifdef PARTICLES
CASE(0) ! Particles
  Resu=0.
  !resu(1:3)= x(1:3)!*x(1)
#endif
CASE(1) ! Constant
  Resu=1.
  Resu_t=0.
  Resu_tt=0.
CASE(2) ! Coaxial Waveguide
  Frequency=1.
  Amplitude=1.
  zlen=2.5
  r=0.5
  r2=(x(1)*x(1)+x(2)*x(2))/r
  omega=Frequency*2.*Pi/zlen
  resu   =0.
  resu(1)=( x(1))*sin(omega*(x(3)-c*t))/r2
  resu(2)=( x(2))*sin(omega*(x(3)-c*t))/r2
  resu(4)=(-x(2))*sin(omega*(x(3)-c*t))/(r2*c)
  resu(5)=( x(1))*sin(omega*(x(3)-c*t))/(r2*c)

  Resu_t=0.
  resu_t(1)=-omega*c*( x(1))*cos(omega*(x(3)-c*t))/r2
  resu_t(2)=-omega*c*( x(2))*cos(omega*(x(3)-c*t))/r2
  resu_t(4)=-omega*c*(-x(2))*cos(omega*(x(3)-c*t))/(r2*c)
  resu_t(5)=-omega*c*( x(1))*cos(omega*(x(3)-c*t))/(r2*c)
  Resu_tt=0.
  resu_tt(1)=-(omega*c)**2*( x(1))*sin(omega*(x(3)-c*t))/r2
  resu_tt(2)=-(omega*c)**2*( x(2))*sin(omega*(x(3)-c*t))/r2
  resu_tt(4)=-(omega*c)**2*(-x(2))*sin(omega*(x(3)-c*t))/(r2*c)
  resu_tt(5)=-(omega*c)**2*( x(1))*sin(omega*(x(3)-c*t))/(r2*c)
CASE(3) ! Resonator
  !special initial values
  !geometric perameters
  a=1.5; b=1.0; d=3.0
  !time parameters
  l=5.; m=4.; n=3.; B0=1.
  IF(a.eq.0)THEN
    CALL abort(&
      __STAMP__&
      ,' Parameter a of resonator is zero!')
  END IF
  IF(b.eq.0)THEN
    CALL abort(&
      __STAMP__&
      ,' Parameter b of resonator is zero!')
  END IF
  IF(d.eq.0)THEN
    CALL abort(&
      __STAMP__&
      ,' Parameter d of resonator is zero!')
  END IF
  omega = Pi*c*sqrt((m/a)**2+(n/b)**2+(l/d)**2)
  gamma = sqrt((omega/c)**2-(l*pi/d)**2)
  IF(gamma.eq.0)THEN
    CALL abort(&
    __STAMP__&
    ,' gamma is computed to zero!')
  END IF
  Psi      =   B0          * cos((m*pi/a)*x(1)) * cos((n*pi/b)*x(2))
  GradPsiX = -(B0*(m*pi/a) * sin((m*pi/a)*x(1)) * cos((n*pi/b)*x(2)))
  GradPsiY = -(B0*(n*pi/b) * cos((m*pi/a)*x(1)) * sin((n*pi/b)*x(2)))

  resu(1)= (-omega/gamma**2) * sin((l*pi/d)*x(3)) *(-GradPsiY)* sin(omega*t)
  resu(2)= (-omega/gamma**2) * sin((l*pi/d)*x(3)) *  GradPsiX * sin(omega*t)
  resu(3)= 0.0
  resu(4)=(1/gamma**2)*(l*pi/d) * cos((l*pi/d)*x(3)) * GradPsiX * cos(omega*t)
  resu(5)=(1/gamma**2)*(l*pi/d) * cos((l*pi/d)*x(3)) * GradPsiY * cos(omega*t)
  resu(6)= Psi                  * sin((l*pi/d)*x(3))            * cos(omega*t)
  resu(7)=0.
  resu(8)=0.

CASE(4) ! Dipole
  resu(1:8) = 0.
  RETURN
  eps=1e-10
  xrel    = x - xDipole
  r = SQRT(DOT_PRODUCT(xrel,xrel))
  IF (r.LT.eps) RETURN
  IF (xrel(3).GT.eps) THEN
    theta = ATAN(SQRT(xrel(1)**2+xrel(2)**2)/xrel(3))
  ELSE IF (xrel(3).LT.(-eps)) THEN
    theta = ATAN(SQRT(xrel(1)**2+xrel(2)**2)/xrel(3)) + pi
  ELSE
    theta = 0.5*pi
  END IF
  IF (xrel(1).GT.eps)      THEN
    phi = ATAN(xrel(2)/xrel(1))
  ELSE IF (xrel(1).LT.eps) THEN
    phi = ATAN(xrel(2)/xrel(1)) + pi
  ELSE IF (xrel(2).GT.eps) THEN
    phi = 0.5*pi
  ELSE IF (xrel(2).LT.eps) THEN
    phi = 1.5*pi
  ELSE
    phi = 0.0                                                                                     ! Vorsicht: phi ist hier undef!
  END IF

  Er = 2.*cos(theta)*Q*dD/(4.*pi*eps0) * ( 1./r**3*sin(omegaD*t-omegaD*r/c) + (omegaD/(c*r**2)*cos(omegaD*t-omegaD*r/c) ) )
  Etheta = sin(theta)*Q*dD/(4.*pi*eps0) * ( (1./r**3-omegaD**2/(c**2*r))*sin(omegaD*t-omegaD*r/c) &
          + (omegaD/(c*r**2)*cos(omegaD*t-omegaD* r/c) ) )
  Bphi = 1/(c2*eps0)*omegaD*sin(theta)*Q*dD/(4.*pi) &
       * ( - omegaD/(c*r)*sin(omegaD*t-omegaD*r/c) + 1./r**2*cos(omegaD*t-omegaD*r/c) )
  IF (ABS(phi).GT.eps) THEN
    resu(1)= sin(theta)*cos(phi)*Er + cos(theta)*cos(phi)*Etheta
    resu(2)= sin(theta)*sin(phi)*Er + cos(theta)*sin(phi)*Etheta
    resu(3)= cos(theta)         *Er - sin(theta)         *Etheta
    resu(4)=-sin(phi)*Bphi
    resu(5)= cos(phi)*Bphi
    resu(6)= 0.0
  ELSE
    resu(3)= cos(theta)         *Er - sin(theta)         *Etheta
  END IF
CASE(5) ! Initialization and BC Gyrotron Mode Converter
  eps=1e-10
  IF (x(3).GT.eps) RETURN
  r=SQRT(x(1)**2+x(2)**2)
  IF (x(1).GT.eps)      THEN
    phi = ATAN(x(2)/x(1))
  ELSE IF (x(1).LT.(-eps)) THEN
    phi = ATAN(x(2)/x(1)) + pi
  ELSE IF (x(2).GT.eps) THEN
    phi = 0.5*pi
  ELSE IF (x(2).LT.(-eps)) THEN
    phi = 1.5*pi
  ELSE
    phi = 0.0                                                                                     ! Vorsicht: phi ist hier undef!
  END IF
  z = x(3)
  Er  =-B0G*mG*omegaG/(r*g**2)*BESSEL_JN(mG,REAL(g*r))                             * &
                                                                 ( cos(h*z+mG*phi)*cos(omegaG*t)+sin(h*z+mG*phi)*sin(omegaG*t))
  Ephi= B0G*omegaG/h      *0.5*(BESSEL_JN(mG-1,REAL(g*r))-BESSEL_JN(mG+1,REAL(g*r)))* &
                                                                 (-cos(h*z+mG*phi)*sin(omegaG*t)+sin(h*z+mG*phi)*cos(omegaG*t))
  Br  =-B0G*h/g           *0.5*(BESSEL_JN(mG-1,REAL(g*r))-BESSEL_JN(mG+1,REAL(g*r)))* &
                                                                 (-cos(h*z+mG*phi)*sin(omegaG*t)+sin(h*z+mG*phi)*cos(omegaG*t))
  Bphi=-B0G*mG*h/(r*g**2)     *BESSEL_JN(mG,REAL(g*r))                             * &
                                                                 ( cos(h*z+mG*phi)*cos(omegaG*t)+sin(h*z+mG*phi)*sin(omegaG*t))
  resu(1)= cos(phi)*Er - sin(phi)*Ephi
  resu(2)= sin(phi)*Er + cos(phi)*Ephi
  resu(3)= 0.0
  resu(4)= cos(phi)*Br - sin(phi)*Bphi
  resu(5)= sin(phi)*Br + cos(phi)*Bphi
  resu(6)= B0G*BESSEL_JN(mG,REAL(g*r))*cos(h*z+mG*phi-omegaG*t)
  resu(7)= 0.0
  resu(8)= 0.0

CASE(6)
  resu   = 0.
  resu(1)= x(1)

CASE(50,51)            ! Initialization and BC Gyrotron - including derivatives
  eps=1e-10
  IF ((ExactFunction.EQ.51).AND.(x(3).GT.eps)) RETURN
  r=SQRT(x(1)**2+x(2)**2)
  IF (x(1).GT.eps)      THEN
    phi = ATAN(x(2)/x(1))
  ELSE IF (x(1).LT.(-eps)) THEN
    phi = ATAN(x(2)/x(1)) + pi
  ELSE IF (x(2).GT.eps) THEN
    phi = 0.5*pi
  ELSE IF (x(2).LT.(-eps)) THEN
    phi = 1.5*pi
  ELSE
    phi = 0.0                                                                                     ! Vorsicht: phi ist hier undef!
  END IF
  z = x(3)
  a = h*z+mG*phi
  b0 = BESSEL_JN(mG,REAL(g*r))
  b1 = BESSEL_JN(mG-1,REAL(g*r))
  b2 = BESSEL_JN(mG+1,REAL(g*r))
  SELECT CASE(MOD(tDeriv,4))
    CASE(0)
      c1  =  omegaG**tDeriv * cos(a-omegaG*t)
      s1  =  omegaG**tDeriv * sin(a-omegaG*t)
    CASE(1)
      c1  =  omegaG**tDeriv * sin(a-omegaG*t)
      s1  = -omegaG**tDeriv * cos(a-omegaG*t)
    CASE(2)
      c1  = -omegaG**tDeriv * cos(a-omegaG*t)
      s1  = -omegaG**tDeriv * sin(a-omegaG*t)
    CASE(3)
      c1  = -omegaG**tDeriv * sin(a-omegaG*t)
      s1  =  omegaG**tDeriv * cos(a-omegaG*t)
    CASE DEFAULT
      CALL abort(&
          __STAMP__&
          , &
        'What is that weird tDeriv you gave me?',999,999.)
  END SELECT

  Er  =-B0G*mG*omegaG/(r*g**2)*b0     *c1
  Ephi= B0G*omegaG/h      *0.5*(b1-b2)*s1
  Br  =-B0G*h/g           *0.5*(b1-b2)*s1
  Bphi=-B0G*mG*h/(r*g**2)     *b0     *c1
  Bz  = B0G                   *b0     *c1
  resu(1)= cos(phi)*Er - sin(phi)*Ephi
  resu(2)= sin(phi)*Er + cos(phi)*Ephi
  resu(3)= 0.0
  resu(4)= cos(phi)*Br - sin(phi)*Bphi
  resu(5)= sin(phi)*Br + cos(phi)*Bphi
  resu(6)= Bz
  resu(7)= 0.0
  resu(8)= 0.0

CASE DEFAULT
  SWRITE(*,*)'Exact function not specified'
END SELECT ! ExactFunction

# if (PP_TimeDiscMethod==1)
! For O3 RK, the boundary condition has to be adjusted
! Works only for O3 RK!!
SELECT CASE(tDeriv)
CASE(0)
  ! resu = g(t)
CASE(1)
  ! resu = g(t) + dt/3*g'(t)
  Resu=Resu + dt/3.*Resu_t
CASE(2)
  ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
  Resu=Resu + 0.75*dt*Resu_t+5./16.*dt*dt*Resu_tt
CASE DEFAULT
  ! Stop, works only for 3 Stage O3 LS RK
  CALL abort(&
      __STAMP__&
      ,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(t,coeff,Ut)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: eps0
USE MOD_Equation_Vars ,ONLY: IniExactFunc,xDipole
#ifdef PARTICLES
USE MOD_Equation_Vars ,ONLY: c_corr
USE MOD_PICDepo_Vars  ,ONLY: PartSource,DoDeposition
#endif /*PARTICLES*/
USE MOD_Mesh_Vars     ,ONLY: Elem_xGP
#ifdef LSERK
USE MOD_Equation_Vars ,ONLY: DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars ,ONLY: dt, sdtCFLOne
USE MOD_DG_Vars       ,ONLY: U
#endif /*LSERK*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
REAL,INTENT(IN)                 :: coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: eps0inv
REAL                            :: r                                                 ! for Dipole
REAL,PARAMETER                  :: Q=1, d=1, omega=6.28318E8 !2.096     ! for Dipole
!===================================================================================================================================
eps0inv = 1./eps0
#ifdef PARTICLES
IF(DoDeposition)THEN
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      !  Get PartSource from Particles
      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - coeff*eps0inv * PartSource(1:3,i,j,k,iElem)
      Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + coeff*eps0inv * PartSource(  4,i,j,k,iElem) * c_corr
    END DO; END DO; END DO
  END DO
END IF
#endif /*PARTICLES*/

SELECT CASE (IniExactFunc)
CASE(0) ! Particles
CASE(1) ! Constant          - no sources
CASE(2) ! Coaxial Waveguide - no sources
CASE(3) ! Resonator         - no sources
CASE(4) ! Dipole
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      r = SQRT(DOT_PRODUCT(Elem_xGP(:,i,j,k,iElem)-xDipole,Elem_xGP(:,i,j,k,iElem)-xDipole))
      IF (shapefunc(r) .GT. 0 ) THEN
        Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) - (shapefunc(r)) * coeff*Q*d*omega * COS(omega*t) * eps0inv
    ! dipole should be neutral
       ! Ut(8,i,j,k,iElem) = Ut(8,i,j,k,iElem) + (shapefunc(r)) * c_corr*Q * eps0inv
      END IF
    END DO; END DO; END DO
  END DO
CASE(5) ! TE_34,19 Mode     - no sources
CASE(50,51) ! TE_34,19 Mode - no sources
CASE DEFAULT
  CALL abort(&
      __STAMP__&
      ,'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
#ifdef LSERK
IF(DoParabolicDamping)THEN
  !Ut(7:8,:,:,:,:) = Ut(7:8,:,:,:,:) - (1.0-fDamping)*sdtCFL1/RK_b(iStage)*U(7:8,:,:,:,:)
  Ut(7:8,:,:,:,:) = Ut(7:8,:,:,:,:) - (1.0-fDamping)*sdtCFL1*U(7:8,:,:,:,:)
END IF
#endif /*LSERK*/

END SUBROUTINE CalcSource

SUBROUTINE CalcSource_Pois(t)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: Phit,Phi
USE MOD_DG_Vars       ,ONLY: U
USE MOD_Globals_Vars  ,ONLY: eps0,c_corr
USE MOD_Equation_Vars ,ONLY: IniExactFunc
#ifdef PARTICLES
USE MOD_PICDepo_Vars  ,ONLY: PartSource,DoDeposition
#endif /*PARTICLES*/
#ifdef LSERK
USE MOD_Equation_Vars ,ONLY: DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars ,ONLY: dt, sdtCFLOne
#endif /*LSERK*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: eps0inv
!===================================================================================================================================
eps0inv = 1./eps0
#ifdef PARTICLES
IF(DoDeposition)THEN
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      !  Get source from Particles

      Phit(  2:4,i,j,k,iElem) = Phit(  2:4,i,j,k,iElem) - U(  1:3,i,j,k,iElem)*c_corr
    END DO; END DO; END DO
  END DO
END IF
#endif /*PARTICLES*/
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
CASE(1) ! Constant          - no sources
CASE DEFAULT
  CALL abort(&
      __STAMP__&
      ,'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
#ifdef LSERK
IF(DoParabolicDamping)THEN
  !Phit(2:4,:,:,:,:) = Phit(2:4,:,:,:,:) - (1.0-fDamping)*sdtCFL1/RK_b(iStage)*Phi(2:4,:,:,:,:)
  Phit(2:4,:,:,:,:) = Phit(2:4,:,:,:,:) - (1.0-fDamping)*sdtCFL1*Phi(2:4,:,:,:,:)
END IF
#endif /*LSERK*/
END SUBROUTINE CalcSource_Pois

SUBROUTINE DivCleaningDamping()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : U
USE MOD_Equation_Vars, ONLY : fDamping,DoParabolicDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================

  IF(DoParabolicDamping) RETURN
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      !  Get source from Particles
      U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
    END DO; END DO; END DO
  END DO
END SUBROUTINE DivCleaningDamping


SUBROUTINE DivCleaningDamping_Pois()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars,       ONLY : Phi
USE MOD_Equation_Vars, ONLY : fDamping_pois,DoParabolicDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================

  IF(DoParabolicDamping) RETURN
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      !  Get source from Particles
      Phi(2:4,i,j,k,iElem) = Phi(2:4,i,j,k,iElem) * fDamping_pois
    END DO; END DO; END DO
  END DO
END SUBROUTINE DivCleaningDamping_Pois

FUNCTION shapefunc(r)
!===================================================================================================================================
! Implementation of (possibly several different) shapefunctions
!===================================================================================================================================
! MODULES
  USE MOD_Equation_Vars, ONLY : shapeFuncPrefix, alpha_shape, rCutoff
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL                 :: r         ! radius / distance to center
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL                 :: shapefunc ! sort of a weight for the source
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
   IF (r.GE.rCutoff) THEN
     shapefunc = 0.0
   ELSE
     shapefunc = ShapeFuncPrefix *(1-(r/rCutoff)**2)**alpha_shape
   END IF
END FUNCTION shapefunc

FUNCTION beta(z,w)
   IMPLICIT NONE
   REAL beta, w, z
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)
END FUNCTION beta

SUBROUTINE FinalizeEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation



SUBROUTINE EvalGradient()
!===================================================================================================================================
! Computes the gradient of the conservative variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde, sJ
USE MOD_Equation_Vars,ONLY:D,E,Phi
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)            :: gradPhi_xi,gradPhi_eta,gradPhi_zeta
INTEGER                                :: i,j,k,l,iElem
!INTEGER,SAVE                           :: N_old=0
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
          gradPhi_xi(i,j,k)  = gradPhi_xi(i,j,k)   + D(i,l) * Phi(1,l,j,k,iElem)
          gradPhi_eta(i,j,k) = gradPhi_eta(i,j,k)  + D(j,l) * Phi(1,i,l,k,iElem)
          gradPhi_zeta(i,j,k)= gradPhi_zeta(i,j,k) + D(k,l) * Phi(1,i,j,l,iElem)
       END DO ! i
     END DO ! j
    END DO ! k
  END DO ! l
  ! Transform the gradients from the reference system to the xyz-System. Only exact for cartesian mesh!
  DO k=0,N
    DO j=0,N
      DO i=0,N
        E(1,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(1,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(1,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(1,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
        E(2,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(2,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(2,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(2,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
        E(3,i,j,k,iElem) = -1*sJ(i,j,k,iElem) * (                                   &
                          Metrics_fTilde(3,i,j,k,iElem) * gradPhi_xi(i,j,k)   + &
                          Metrics_gTilde(3,i,j,k,iElem) * gradPhi_eta(i,j,k)  + &
                          Metrics_hTilde(3,i,j,k,iElem) * gradPhi_zeta(i,j,k)   )
      END DO ! i
    END DO ! j
  END DO ! k
END DO
END SUBROUTINE EvalGradient

SUBROUTINE VolInt_weakForm(Ut)
!===================================================================================================================================
! Computes the volume integral of the weak DG form a la Kopriva
! Attention 1: 1/J(i,j,k) is not yet accounted for
! Attention 2: input Ut=0. and is updated with the volume flux derivatives
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:D_hat
USE MOD_Mesh_Vars,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_PreProc
USE MOD_Flux_Pois,ONLY:EvalFlux3D_Pois                                         ! computes volume fluxes in local coordinates
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: Ut(4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4,0:PP_N,0:PP_N,0:PP_N)      :: f,g,h                ! volume fluxes at all Gauss points
REAL,DIMENSION(4)                           :: fTilde,gTilde,hTilde ! auxiliary variables needed to store the fluxes at one GP
INTEGER                                           :: i,j,k,iElem
INTEGER                                           :: l                    ! row index for matrix vector product
!===================================================================================================================================
DO iElem=1,PP_nElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3D_Pois(iElem,f,g,h)
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        fTilde=f(:,i,j,k)
        gTilde=g(:,i,j,k)
        hTilde=h(:,i,j,k)
        ! Compute the transformed fluxes with the metric terms
        ! Attention 1: we store the transformed fluxes in f,g,h again
        f(:,i,j,k) = fTilde(:)*Metrics_fTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_fTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_fTilde(3,i,j,k,iElem)
        g(:,i,j,k) = fTilde(:)*Metrics_gTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_gTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_gTilde(3,i,j,k,iElem)
        h(:,i,j,k) = fTilde(:)*Metrics_hTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_hTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_hTilde(3,i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO ! k
  DO l=0,PP_N
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! Update the time derivative with the spatial derivatives of the transformed fluxes
          Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*f(:,l,j,k) + &
                                                  D_hat(j,l)*g(:,i,l,k) + &
                                                  D_hat(k,l)*h(:,i,j,l)
        END DO !i
      END DO ! j
    END DO ! k
  END DO ! l
END DO ! iElem
END SUBROUTINE VolInt_weakForm

SUBROUTINE FillFlux(Flux,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: Phi_master,Phi_slave
USE MOD_Mesh_Vars     ,ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars     ,ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Riemann_Pois  ,ONLY: Riemann_Pois
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:4,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = nBCSides+nInnerSides+1
  lastSideID  = firstSideID-1+nMPISides_MINE
ELSE
  ! fill only InnerSides
  firstSideID = nBCSides+1
  lastSideID  = firstSideID-1+nInnerSides
END IF
!firstSideID=nBCSides+1
!lastSideID  =nBCSides+nInnerSides+nMPISides_MINE
DO SideID=firstSideID,lastSideID
  CALL Riemann_Pois(Flux(:,:,:,SideID),Phi_master(:,:,:,SideID),Phi_slave(:,:,:,SideID),NormVec(:,:,:,SideID))
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO ! SideID

END SUBROUTINE FillFlux

SUBROUTINE ProlongToFace_SideBased(Uvol,Uface_Minus,Uface_Plus,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
REAL,INTENT(IN)                 :: Uvol(4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Uface_Minus(4,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)
REAL,INTENT(INOUT)              :: Uface_Plus(4,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,l,p,q,ElemID(2),SideID,flip(2),LocSideID(2),firstSideID,lastSideID
REAL                            :: Uface(4,0:PP_N,0:PP_N)
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE+1
  lastSideID  = firstSideID-1+nMPISides_YOUR
  flip(1)      = -1
ELSE
  ! BCSides, InnerSides and MINE MPISides are filled
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
  flip(1)      = 0
END IF
DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)     = SideToElem(S2E_ELEM_ID,SideID)
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)
  DO i=1,2 !first maste then slave side
#if (PP_NodeType==1) /* for Gauss-points*/
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,q,p)=Uface(:,q,p)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,q,p)=Uface(:,q,p)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,PP_N-p,q)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    END SELECT
#else /* for Gauss-Lobatto-points*/
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,0,p,q,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      Uface(:,:,:)=Uvol(:,:,0,:,ElemID(i))
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,p,q,0,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      Uface(:,:,:)=Uvol(:,PP_N,:,:,ElemID(i))
    CASE(ETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,PP_N-p,q)=Uvol(:,p,PP_N,q,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,q,PP_N,ElemID(i))
        END DO ! p
      END DO ! q
    END SELECT
#endif
    SELECT CASE(Flip(i))
      CASE(0) ! master side
        Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,q,p)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-p,q)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-q,PP_N-p)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,p,PP_N-q)
          END DO ! p
        END DO ! q
    END SELECT
  END DO !i=1,2, masterside & slave side
END DO !SideID
END SUBROUTINE ProlongToFace_SideBased

#if USE_MPI
SUBROUTINE StartExchangeMPIData_Pois(FaceData,LowerBound,UpperBound,SendRequest,RecRequest,SendID)
!===================================================================================================================================
! Subroutine does the send and receive operations for the face data that has to be exchanged between processors.
! FaceData: the complete face data (for inner, BC and MPI sides).
! LowerBound / UpperBound: lower side index and upper side index for last dimension of FaceData
! SendRequest, RecRequest: communication handles
! SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / recieve MINE
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
INTEGER, INTENT(IN)          :: LowerBound,UpperBound
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)         :: SendRequest(nNbProcs),RecRequest(nNbProcs)
REAL, INTENT(INOUT)          :: FaceData(1:4,0:PP_N,0:PP_N,LowerBound:UpperBound)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =4*(PP_N+1)*(PP_N+1)*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,SendRequest(iNbProc),iError)
  END IF
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =4*(PP_N+1)*(PP_N+1)*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_PICLAS,RecRequest(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartExchangeMPIData_Pois
#endif

SUBROUTINE SurfInt2(Flux,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
REAL,INTENT(IN)    :: Flux(1:4,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(4,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,ElemID(2),p,q,l,Flip(2),SideID,locSideID(2)
INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType>1)
REAL            ::L_HatMinus0,L_HatPlusN
#endif
!===================================================================================================================================
IF(doMPISides)THEN
  ! surfInt only for YOUR MPISides
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE +1
  lastSideID  = firstSideID-1+nMPISides_YOUR
ELSE
  ! fill only InnerSides
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
END IF

#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif
flip(1)        = 0 !flip=0 for master side
DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! neighbor side
  ElemID(2)    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)

  DO i=1,2
  ! update DG time derivative with corresponding SurfInt contribution
#if (PP_NodeType==1)
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))+Flux(:,q,p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,p,q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT

    CASE(ETA_MINUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))+Flux(:,p,q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,q,p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,l,q
      END SELECT

    CASE(ZETA_MINUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))+Flux(:,q,p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(1) ! slave side, SideID=q,jSide=p
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,p,q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatMinus(l)
        END DO; END DO; END DO ! p,q,l
      END SELECT

    CASE(XI_PLUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))+Flux(:,p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,q,p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
          Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT

    CASE(ETA_PLUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))+Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,l,q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
          Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,l,q
      END SELECT

    CASE(ZETA_PLUS)
      SELECT CASE(flip(i))
      CASE(0) ! master side
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))+Flux(:,p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(1) ! slave side, SideID=q,jSide=p
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,q,p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,q,l
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatPlus(l)
        END DO; END DO; END DO ! p,q,l
      END SELECT
    END SELECT !locSideID
#else
    !update local grid cell
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))+Flux(:,q,p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))-Flux(:,p,q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      END SELECT

    ! switch to right hand system for ETA_PLUS direction
    CASE(ETA_MINUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))+Flux(:,p,q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))-Flux(:,q,p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      END SELECT

    ! switch to right hand system for ZETA_MINUS direction
    CASE(ZETA_MINUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))+Flux(:,q,p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))-Flux(:,p,q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatMinus0
        END DO; END DO ! p,q
      END SELECT

    CASE(XI_PLUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))+Flux(:,p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))-Flux(:,q,p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      END SELECT

    ! switch to right hand system for ETA_PLUS direction
    CASE(ETA_PLUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))+Flux(:,PP_N-p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))-Flux(:,q,PP_N-p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))-Flux(:,p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))-Flux(:,PP_N-q,p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))-Flux(:,PP_N-p,PP_N-q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      END SELECT

    ! switch to right hand system for ZETA_MINUS direction
    CASE(ZETA_PLUS)
      SELECT CASE(flip(i))
      CASE(0)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))+Flux(:,p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(1)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))-Flux(:,q,p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(2)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))-Flux(:,PP_N-p,q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(3)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))-Flux(:,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      CASE(4)
        DO q=0,PP_N; DO p=0,PP_N
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))-Flux(:,p,PP_N-q,SideID)*L_hatPlusN
        END DO; END DO ! p,q
      END SELECT
    END SELECT !locSideID
#endif
  END DO ! i=1,2 master side, slave side
END DO ! SideID=1,nSides
END SUBROUTINE SurfInt2

END MODULE MOD_Equation

