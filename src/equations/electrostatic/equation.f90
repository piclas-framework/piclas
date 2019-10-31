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
                                                     'Multiplied with c0 results in the velocity of '//&
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
CALL prms%CreateLogicalOption(  'ParabolicDamping' , 'TODO-DEFINE-PARAMETER' , '.FALSE.')
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection')

CALL prms%CreateStringOption(   'BCStateFile'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Boundary Condition State File', 'no file found')
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
#if USE_MPI
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
smu0               = 1./mu0
fDamping           = GETREAL('fDamping','0.99')
DoParabolicDamping = GETLOGICAL('ParabolicDamping','.FALSE.')
c_test = 1./SQRT(eps0*mu0)
IF ( ABS(c-c_test)/c.GT.10E-8) THEN
  SWRITE(*,*) "ERROR: c does not equal 1/sqrt(eps*mu)!"
  SWRITE(*,*) "c:", c
  SWRITE(*,*) "mu:", mu0
  SWRITE(*,*) "eps:", eps0
  SWRITE(*,*) "1/sqrt(eps*mu):", c_test
  CALL abort(&
      __STAMP__&
      ,' Speed of light coefficients does not match!')
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
BCStateFile=GETSTR('BCStateFile','')
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
USE MOD_Globals,       ONLY : abort
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY : eps0,c_corr,IniExactFunc
#ifdef PARTICLES
USE MOD_PICDepo_Vars,  ONLY : PartSource,DoDeposition
#endif /*PARTICLES*/
USE MOD_Mesh_Vars,     ONLY : Elem_xGP                  ! for shape function: xyz position of the Gauss points
#ifdef LSERK
USE MOD_Equation_Vars, ONLY : DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars, ONLY : dt, sdtCFLOne!, RK_B, iStage
USE MOD_DG_Vars,       ONLY : U
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
!===================================================================================================================================
eps0inv = 1./eps0

#ifdef PARTICLES
IF(DoDeposition)THEN
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      !  Get source from Particles
      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - coeff*eps0inv * PartSource(1:3,i,j,k,iElem)
      Ut(  4,i,j,k,iElem) = Ut(  4,i,j,k,iElem) + coeff*eps0inv * PartSource(  4,i,j,k,iElem) * c_corr
    END DO; END DO; END DO
  END DO
END IF
#endif /*PARTICLES*/

SELECT CASE (IniExactFunc)
CASE(0) ! default, empty
CASE(1) ! Constant          - no sources
CASE DEFAULT
  CALL abort(&
      __STAMP__&
      ,'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction

#ifdef LSERK
IF(DoParabolicDamping)THEN
  !Ut(4,:,:,:,:) = Ut(4,:,:,:,:) - (1.0-fDamping)*sdtCFL1*/RK_b(iStage)*U(4,:,:,:,:)
  Ut(4,:,:,:,:) = Ut(4,:,:,:,:) - (1.0-fDamping)*sdtCFL1*U(4,:,:,:,:)
END IF
#endif /*LSERK*/

END SUBROUTINE CalcSource

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

END MODULE MOD_Equation

