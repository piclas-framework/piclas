#include "boltzplatz.h"

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
PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation,DivCleaningDamping
!===================================================================================================================================

CONTAINS

SUBROUTINE InitEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,ONLY:PI
USE MOD_ReadInTools
#ifdef PARTICLES
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
#endif
USE MOD_Equation_Vars 
USE MOD_TimeDisc_Vars, ONLY: TEnd
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: c_test
#ifdef MPI
#endif
!===================================================================================================================================
! Read the maximum number of time steps MaxIter and the end time TEnd from ini file
TEnd=GetReal('TEnd') ! must be read in here due to DSMC_init
IF(EquationInitIsDone)THEN
#ifdef PARTICLES
  IF(InterpolationInitIsDone)THEN
    SWRITE(*,*) "InitMaxwell not ready to be called or already called."
    RETURN
  END IF
#else
  SWRITE(*,*) "InitMaxwell not ready to be called or already called."
  RETURN
#endif /*PARTICLES*/
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MAXWELL ...' 

! Read correction velocity
c_corr             = GETREAL('c_corr','1.')
c                  = GETREAL('c0','1.')
eps0               = GETREAL('eps','1.')
mu0                = GETREAL('mu','1.')
smu0               = 1./mu0
fDamping           = GETREAL('fDamping','0.99')
!scr            = 1./ GETREAL('c_r','0.18')  !constant for damping
DipoleOmega        = GETREAL('omega','6.28318E08') ! f=100 MHz default
tPulse             = GETREAL('tPulse','30e-9')     ! half length of pulse
c_test = 1./SQRT(eps0*mu0)
IF ( ABS(c-c_test)/c.GT.10E-8) THEN
  SWRITE(*,*) "ERROR: c does not equal 1/sqrt(eps*mu)!"
  SWRITE(*,*) "c:", c
  SWRITE(*,*) "mu:", mu0
  SWRITE(*,*) "eps:", eps0
  SWRITE(*,*) "1/sqrt(eps*mu):", c_test
  STOP
END IF



c2     = c*c 
c_inv  = 1./c
c2_inv = 1./c2

c_corr2   = c_corr*c_corr
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
            
EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MAXWELL DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



SUBROUTINE ExactFunc(ExactFunction,t,tDeriv,x,resu) 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE nr,only:bessj
USE nrtype,only:SP
USE MOD_Globals
USE MOD_Globals_Vars,ONLY:PI
USE MOD_Equation_Vars,ONLY:c,c2,eps0,c_inv
# if (PP_TimeDiscMethod==1)
USE MOD_TimeDisc_vars,ONLY:dt
# endif
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
REAL,PARAMETER                  :: xDipole(1:3)=(/0,0,0/)          ! aux. Constants for Dipole
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
  omega=Frequency*2.*Pi/zlen ! shift beruecksichtigen
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
    ERRWRITE(*,*)'ERROR: a eq 0!'
    STOP
  END IF
  IF(b.eq.0)THEN
    ERRWRITE(*,*)'ERROR: b eq 0!'
    STOP
  END IF
  IF(d.eq.0)THEN
    ERRWRITE(*,*)'ERROR: d eq 0!'
    STOP
  END IF
  omega = Pi*c*sqrt((m/a)**2+(n/b)**2+(l/d)**2)
  gamma = sqrt((omega/c)**2-(l*pi/d)**2)
  IF(gamma.eq.0)THEN
    ERRWRITE(*,*)'ERROR: gamma eq 0!'
    STOP
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
  Er  =-B0G*mG*omegaG/(r*g**2)*bessj(mG,REAL(g*r,SP))                             * &
                                                                 ( cos(h*z+mG*phi)*cos(omegaG*t)+sin(h*z+mG*phi)*sin(omegaG*t))
  Ephi= B0G*omegaG/h      *0.5*(bessj(mG-1,REAL(g*r,SP))-bessj(mG+1,REAL(g*r,SP)))* &
                                                                 (-cos(h*z+mG*phi)*sin(omegaG*t)+sin(h*z+mG*phi)*cos(omegaG*t))
  Br  =-B0G*h/g           *0.5*(bessj(mG-1,REAL(g*r,SP))-bessj(mG+1,REAL(g*r,SP)))* &
                                                                 (-cos(h*z+mG*phi)*sin(omegaG*t)+sin(h*z+mG*phi)*cos(omegaG*t))
  Bphi=-B0G*mG*h/(r*g**2)     *bessj(mG,REAL(g*r,SP))                             * &
                                                                 ( cos(h*z+mG*phi)*cos(omegaG*t)+sin(h*z+mG*phi)*sin(omegaG*t))
  resu(1)= cos(phi)*Er - sin(phi)*Ephi
  resu(2)= sin(phi)*Er + cos(phi)*Ephi
  resu(3)= 0.0
  resu(4)= cos(phi)*Br - sin(phi)*Bphi
  resu(5)= sin(phi)*Br + cos(phi)*Bphi
  resu(6)= B0G*bessj(mG,REAL(g*r,SP))*cos(h*z+mG*phi-omegaG*t)
  resu(7)= 0.0
  resu(8)= 0.0

CASE(7) ! Manufactured Solution
  resu(:)=0
  resu(1)=SIN(2*pi*(x(1)-t))
  resu_t(:)=0
  resu_t(1)=-2*pi*COS(2*pi*(x(1)-t))
  resu_tt(:)=0
  resu_tt(1)=-4*pi*pi*resu(1)

CASE(10) !issautier 3D test case with source (Stock et al., divcorr paper), domain [0;1]^3!!!
  resu(:)=0.
  resu(1)=x(1)*SIN(Pi*x(2))*SIN(Pi*x(3)) !*SIN(t)
  resu(2)=x(2)*SIN(Pi*x(3))*SIN(Pi*x(1)) !*SIN(t)
  resu(3)=x(3)*SIN(Pi*x(1))*SIN(Pi*x(2)) !*SIN(t)
  resu(4)=pi*SIN(Pi*x(1))*(x(3)*COS(Pi*x(2))-x(2)*COS(Pi*x(3))) !*(COS(t)-1)
  resu(5)=pi*SIN(Pi*x(2))*(x(1)*COS(Pi*x(3))-x(3)*COS(Pi*x(1))) !*(COS(t)-1)
  resu(6)=pi*SIN(Pi*x(3))*(x(2)*COS(Pi*x(1))-x(1)*COS(Pi*x(2))) !*(COS(t)-1)

  resu_t(:)=0.
  resu_t(1)= COS(t)*resu(1)
  resu_t(2)= COS(t)*resu(2)
  resu_t(3)= COS(t)*resu(3)
  resu_t(4)=-SIN(t)*resu(4)
  resu_t(5)=-SIN(t)*resu(5)
  resu_t(6)=-SIN(t)*resu(6)
  resu_tt=0.
  resu_tt(1)=-SIN(t)*resu(1)
  resu_tt(2)=-SIN(t)*resu(2)
  resu_tt(3)=-SIN(t)*resu(3)
  resu_tt(4)=-COS(t)*resu(4)
  resu_tt(5)=-COS(t)*resu(5)
  resu_tt(6)=-COS(t)*resu(6)

  resu(1)=     SIN(t)*resu(1)
  resu(2)=     SIN(t)*resu(2)
  resu(3)=     SIN(t)*resu(3)
  resu(4)=(COS(t)-1.)*resu(4)
  resu(5)=(COS(t)-1.)*resu(5)
  resu(6)=(COS(t)-1.)*resu(6)

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
  b0 = bessj(mG,REAL(g*r,SP))
  b1 = bessj(mG-1,REAL(g*r,SP))
  b2 = bessj(mG+1,REAL(g*r,SP))
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
      c1  = 0.0
      s1  = 0.0
      CALL abort(__STAMP__,'What is that weired tDeriv you gave me?',999,999.)
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

CASE(41) ! pulsed Dipole
  resu = 0.0
  RETURN

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
  CALL abort(__STAMP__,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(t)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY : abort
USE MOD_Globals_Vars,  ONLY : PI
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : Ut,U
USE MOD_Equation_Vars, ONLY : eps0,c_corr,IniExactFunc, DipoleOmega, tPulse!,scr
#ifdef PARTICLES
USE MOD_PICDepo_Vars,  ONLY : Source
#endif /*PARTICLES*/
USE MOD_Mesh_Vars,     ONLY : Elem_xGP                  ! for shape function: xyz position of the Gauss points
!USE MOD_PIC_Analyze,   ONLY : CalcDepositedCharge
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
REAL                            :: eps0inv, x(1:3)
REAL                            :: r                                                 ! for Dipole
REAL,PARAMETER                  :: xDipole(1:3)=(/0,0,0/), Q=1, d=1    ! for Dipole
!===================================================================================================================================
eps0inv = 1./eps0
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
#ifdef PARTICLES
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      !  Get source from Particles
      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv * source(1:3,i,j,k,iElem)
      Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv * source(  4,i,j,k,iElem) * c_corr 
    END DO; END DO; END DO
  END DO
#endif /*PARTICLES*/
  !CALL CalcDepositedCharge()
CASE(1) ! Constant          - no sources
CASE(2) ! Coaxial Waveguide - no sources
CASE(3) ! Resonator         - no sources
CASE(4) ! Dipole
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      r = SQRT(DOT_PRODUCT(Elem_xGP(:,i,j,k,iElem)-xDipole,Elem_xGP(:,i,j,k,iElem)-xDipole))
      IF (shapefunc(r) .GT. 0 ) THEN
        Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) - (shapefunc(r)) * Q*d*DipoleOmega * COS(DipoleOmega*t) * eps0inv
    ! dipole should be neutral
       ! Ut(8,i,j,k,iElem) = Ut(8,i,j,k,iElem) + (shapefunc(r)) * c_corr*Q * eps0inv
      END IF
    END DO; END DO; END DO
  END DO
CASE(5) ! TE_34,19 Mode     - no sources
CASE(7) ! Manufactured Solution
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - 2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * eps0inv
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + 2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * c_corr * eps0inv
    END DO; END DO; END DO
  END DO
CASE(10) !issautier 3D test case with source (Stock et al., divcorr paper), domain [0;1]^3!!!
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N  
      x(:)=Elem_xGP(:,i,j,k,iElem)
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) + (COS(t)- (COS(t)-1.)*2*pi*pi)*x(1)*SIN(Pi*x(2))*SIN(Pi*x(3))
      Ut(2,i,j,k,iElem) =Ut(2,i,j,k,iElem) + (COS(t)- (COS(t)-1.)*2*pi*pi)*x(2)*SIN(Pi*x(3))*SIN(Pi*x(1))
      Ut(3,i,j,k,iElem) =Ut(3,i,j,k,iElem) + (COS(t)- (COS(t)-1.)*2*pi*pi)*x(3)*SIN(Pi*x(1))*SIN(Pi*x(2))
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - (COS(t)-1.)*pi*COS(Pi*x(1))*(SIN(Pi*x(2))+SIN(Pi*x(3)))
      Ut(2,i,j,k,iElem) =Ut(2,i,j,k,iElem) - (COS(t)-1.)*pi*COS(Pi*x(2))*(SIN(Pi*x(3))+SIN(Pi*x(1)))
      Ut(3,i,j,k,iElem) =Ut(3,i,j,k,iElem) - (COS(t)-1.)*pi*COS(Pi*x(3))*(SIN(Pi*x(1))+SIN(Pi*x(2)))
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + c_corr*SIN(t)*( SIN(pi*x(2))*SIN(pi*x(3)) &
                                                            +SIN(pi*x(3))*SIN(pi*x(1)) &
                                                            +SIN(pi*x(1))*SIN(pi*x(2)) )
    END DO; END DO; END DO
  END DO


CASE(41) ! Dipole via temporal Gausspuls
!t0=TEnd/5, w=t0/4 ! for pulsed Dipole (t0=offset and w=width of pulse)
!TEnd=30.E-9 -> short pulse for 100ns runtime
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - 2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * eps0inv
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + 2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * c_corr * eps0inv
    END DO; END DO; END DO
  END DO

CASE(50,51) ! TE_34,19 Mode - no sources
CASE DEFAULT
  CALL abort(__STAMP__,'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
!source fo divcorr damping!
!Ut(7:8,:,:,:,:)=Ut(7:8,:,:,:,:)-(c_corr*scr)*U(7:8,:,:,:,:)
END SUBROUTINE CalcSource

SUBROUTINE DivCleaningDamping()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : U
USE MOD_Equation_Vars, ONLY : fDamping
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
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
    !  Get source from Particles
    U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
  END DO; END DO; END DO
END DO
END SUBROUTINE DivCleaningDamping

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
   USE nr
   IMPLICIT NONE
   REAL beta, w, z                                                                                                  
   beta = exp(gammln(z)+gammln(w)-gammln(z+w))                                                                    
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

END MODULE MOD_Equation

