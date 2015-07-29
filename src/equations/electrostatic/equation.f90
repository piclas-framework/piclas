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
IF(InterpolationInitIsDone.AND.EquationInitIsDone)THEN
   SWRITE(*,*) "InitElectrostatic not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ELECTROSTATIC ...'

! Read correction velocity
c_corr             = GETREAL('c_corr','1.')
c                  = GETREAL('c0','1.')
eps0               = GETREAL('eps','1.')
mu0                = GETREAL('mu','1.')
fDamping           = GETREAL('fDamping','0.99')
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
SWRITE(UNIT_stdOut,'(A)')' INIT ELECTROSTATIC DONE!'
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
USE MOD_Equation_Vars,ONLY:c,c2,eps0
USE MOD_TimeDisc_vars,ONLY:dt
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
REAL,PARAMETER                  :: Q=1, dD=1, omegaD=2.096         ! aux. Constants for Dipole
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
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : Ut
USE MOD_Equation_Vars, ONLY : eps0,c_corr,IniExactFunc
#ifdef PARTICLES
USE MOD_PICDepo_Vars,  ONLY : Source
#endif /*PARTICLES*/
USE MOD_Mesh_Vars,     ONLY : Elem_xGP                  ! for shape function: xyz position of the Gauss points
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
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
#ifdef PARTICLES
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      !  Get source from Particles
      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv * source(1:3,i,j,k,iElem)
      Ut(  4,i,j,k,iElem) = Ut(  4,i,j,k,iElem) + eps0inv * source(  4,i,j,k,iElem) * c_corr 
      !IF((t.GT.0).AND.(ABS(source(4,i,j,k,iElem)*c_corr).EQ.0))THEN
      !print*, t
     ! print*, eps0inv * source(4,i,j,k,iElem)*c_corr
      !print*, eps0inv * source(1:3,i,j,k,iElem)
      !read*
      !END IF
    END DO; END DO; END DO
  END DO
#endif /*PARTICLES*/
CASE(1) ! Constant          - no sources
CASE DEFAULT
  CALL abort(__STAMP__,'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
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
      U(4,i,j,k,iElem) = U(4,i,j,k,iElem) * fDamping
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

