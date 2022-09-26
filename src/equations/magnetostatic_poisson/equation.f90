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
INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation, CalcSourceHDG
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

!CALL prms%CreateRealOption(     'c_corr'          , 'TODO-DEFINE-PARAMETER multiplied with c0 results in the velocity of '//&
!                                                     'introduced artificial correcting waves (HDC)' , '1.')
CALL prms%CreateRealOption(     'c0'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Velocity of light (in vacuum)' , '1.')
CALL prms%CreateRealOption(     'eps'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Electric constant (vacuum permittivity)' , '1.')
CALL prms%CreateRealOption(     'mu'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Magnetic constant (vacuum permeability = 4πE−7H/m)' &
                                                   , '1.')
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection')
CALL prms%CreateRealArrayOption('IniWavenumber'    , 'TODO-DEFINE-PARAMETER' , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('IniCenter'        , 'TODO-DEFINE-PARAMETER' , '0. , 0. , 0.')
CALL prms%CreateRealOption(     'IniAmplitude'     , 'TODO-DEFINE-PARAMETER' , '0.1')
CALL prms%CreateRealOption(     'IniHalfwidth'     , 'TODO-DEFINE-PARAMETER' , '0.1')

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
USE MOD_Preproc
USE MOD_ReadInTools        ,ONLY: GETREALARRAY,GETREAL,GETINT
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Equation_Vars
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Mesh_Vars          ,ONLY: nSides
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
   SWRITE(*,*) "InitPoisson not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT POISSON...'

! Read the velocity vector from ini file
IniWavenumber     = GETREALARRAY('IniWavenumber',3,'1.,1.,1.')
! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
IniCenter    = GETREALARRAY('IniCenter',3,'0.,0.,0.')
IniAmplitude = GETREAL('IniAmplitude','0.1')
IniHalfwidth = GETREAL('IniHalfwidth','0.1')

ALLOCATE(chitens(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitensInv(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitens_face(3,3,0:PP_N,0:PP_N,nSides))

chitens=0.
chitens(1,1,:,:,:,:)=1.
chitens(2,2,:,:,:,:)=1.
chitens(3,3,:,:,:,:)=1.
chitensInv=chitens
chitens_face=0.
chitens_face(1,1,:,:,:)=1.
chitens_face(2,2,:,:,:)=1.
chitens_face(3,3,:,:,:)=1.

alpha_shape = GETINT('AlphaShape','2')
rCutoff     = GETREAL('r_cutoff','1.')
! Compute factor for shape function
ShapeFuncPrefix = 1./(2. * beta(1.5, REAL(alpha_shape) + 1.) * REAL(alpha_shape) + 2. * beta(1.5, REAL(alpha_shape) + 1.)) &
                * (REAL(alpha_shape) + 1.)/(PI*(rCutoff**3))

ALLOCATE(B(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(E(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
E=0.
B=0.

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT POISSON DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



SUBROUTINE ExactFunc(ExactFunction,t,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort,MPIRoot,PI
USE MOD_Equation_Vars,ONLY: IniWavenumber
USE MOD_Equation_Vars,ONLY: IniCenter,IniHalfwidth,IniAmplitude
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
REAL,INTENT(OUT)                :: Resu(1:PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3)
REAL                            :: x0(3),r1,r2
!===================================================================================================================================
SELECT CASE (ExactFunction)
CASE(0)
    Resu(:)=0.
CASE(1) !linear
    Resu(:)=0.
CASE(21) !linear
    Resu(:)=10.+SUM(x)
CASE(101) !constant
    Resu(:)=7.7

CASE(2) !sinus
  Frequency=0.5
  Amplitude=0.3
  Omega=2.*PI*Frequency
  Resu(:)=1.+Amplitude*SIN(Omega*SUM(Cent))
CASE(102) !linear: z=-1: 0, z=1, 1000
  resu(:)=(1+x(3))*1000.
CASE(103) !dipole
  r1=SQRT(SUM((x(:)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  r2=SQRT(SUM((x(:)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  resu(:)=IniAmplitude*(1/r2-1/r1)
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Exactfunction not specified!')
END SELECT ! ExactFunction


END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(t)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation_Vars,ONLY:IniCenter,IniHalfwidth,IniAmplitude
USE MOD_DG_Vars,ONLY:Ut, U
USE MOD_Mesh_Vars,ONLY:Elem_xGP
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
REAL                             :: r1,r2
REAL,DIMENSION(3)                :: dx1,dx2,dr1dx,dr2dx,dr1dx2,dr2dx2
!===================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
#ifdef PARTICLES
!  DO iElem=1,PP_nElems
!    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!      !  Get PartSource from Particles
!      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv * PartSource(1:3,i,j,k,iElem)
!      Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv * PartSource(  4,i,j,k,iElem) * c_corr
!    END DO; END DO; END DO
!  END DO
#endif /*PARTICLES*/
CASE(103)
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
         dx1=(Elem_xGP(:,i,j,k,iElem)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))
         dx2=(Elem_xGP(:,i,j,k,iElem)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))
         r1=SQRT(SUM(dx1**2))
         r2=SQRT(SUM(dx2**2))
         dr1dx(:)= r1*dx1
         dr2dx(:)= r2*dx2
         dr1dx2(:)= r1+dr1dx(:)*dx1
         dr2dx2(:)= r2+dr2dx(:)*dx2
         Ut(:,i,j,k,iElem)=Ut(:,i,j,k,iElem)- IniAmplitude*( SUM((r1*dr1dx2(:)-2*dr1dx(:)**2)/(r1*r1*r1)) &
                                 -SUM((r2*dr2dx2(:)-2*dr2dx(:)**2)/(r2*r2*r2)) )
      END DO !i
    END DO !j
  END DO !k
END DO ! iElem=1,nElems

CASE DEFAULT
!  CALL abort(__STAMP__,&
             !'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
END SUBROUTINE CalcSource


SUBROUTINE CalcSourceHDG(t,i,j,k,iElem,resu, Phi)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals            ,ONLY: Abort
USE MOD_PreProc
USE MOD_PICDepo_Vars       ,ONLY: PartSource,DoDeposition
USE MOD_Globals_Vars       ,ONLY: Pi, eps0, mu0
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Equation_Vars      ,ONLY: IniCenter,IniHalfwidth,IniAmplitude
USE MOD_DG_Vars            ,ONLY: Ut,U
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER, INTENT(IN)             :: i, j, k,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
REAL                            :: x(3)
REAL,INTENT(IN),OPTIONAL     :: Phi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: r1,r2, source_e
REAL,DIMENSION(3)                :: dx1,dx2,dr1dx,dr2dx,dr1dx2,dr2dx2
INTEGER                         :: RegionID
!===================================================================================================================================
#ifdef PARTICLES
IF(DoDeposition)THEN
  Resu(1:3) = - PartSource(1:3,i,j,k,iElem)*mu0
  Resu(4) = - PartSource(4,i,j,k,iElem)/eps0
END IF
#endif /*PARTICLES*/
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
CASE(103)
 x(1:3) = Elem_xGP(1:3,i,j,k,iElem)
 dx1=(x(:)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))
 dx2=(x(:)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))
 r1=SQRT(SUM(dx1**2))
 r2=SQRT(SUM(dx2**2))
 dr1dx(:)= r1*dx1
 dr2dx(:)= r2*dx2
 dr1dx2(:)= r1+dr1dx(:)*dx1
 dr2dx2(:)= r2+dr2dx(:)*dx2
 resu(1)=- IniAmplitude*( SUM((r1*dr1dx2(:)-2*dr1dx(:)**2)/(r1*r1*r1)) &
                         -SUM((r2*dr2dx2(:)-2*dr2dx(:)**2)/(r2*r2*r2)) )
CASE DEFAULT
  resu=0.
!  CALL abort(__STAMP__,&
             !'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
END SUBROUTINE CalcSourceHDG


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

