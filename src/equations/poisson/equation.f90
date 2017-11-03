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
INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation, CalcSourceHDG,DivCleaningDamping
!===================================================================================================================================

CONTAINS

SUBROUTINE InitEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,             ONLY:GETREALARRAY,GETREAL,GETINT
USE MOD_Interpolation_Vars,      ONLY:InterpolationInitIsDone
USE MOD_Equation_Vars
USE MOD_HDG_vars
USE MOD_Mesh_Vars,               ONLY:nSides
USE MOD_TimeDisc_Vars,           ONLY:TEnd
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: chitensValue,chitensRadius  ! depricated variables, remove in future (by the end of 2017)
INTEGER                      :: chitensWhichField           ! depricated variables, remove in future (by the end of 2017)
!===================================================================================================================================
TEnd=GetReal('TEnd') 
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
   SWRITE(*,*) "InitPoisson not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT POISSON...'

! Read the velocity vector from ini file
Pi=ACOS(-1.)
IniWavenumber     = GETREALARRAY('IniWavenumber',3,'1.,1.,1.')
c                  = GETREAL('c0','1.')
eps0               = GETREAL('eps','1.')
mu0                = GETREAL('mu','1.')
! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
IniCenter    = GETREALARRAY('IniCenter',3,'0.,0.,0.')
IniAmplitude = GETREAL('IniAmplitude','0.1')
IniHalfwidth = GETREAL('IniHalfwidth','0.1')
ACfrequency = GETREAL('ACfrequency','0.0')
ACamplitude = GETREAL('ACamplitude','0.0')
c_inv  = 1./c
c2     = c*c
smu0=1./mu0
c2_inv = 1./c2

chitensWhichField = GETINT( 'chitensWhichField','-1')
chitensValue      = GETREAL('chitensValue','-1.0')
chitensRadius     = GETREAL('chitensRadius','-1.0')
IF(chitensWhichField.GT.0.0.OR.&
   chitensValue     .GT.0.0.OR.&
   chitensRadius    .GT.0.0)THEN
  CALL abort(&
  __STAMP__&
  ,'chitensWhichField, chitensValue and chitensRadius are no longer supported. Deactivate them!')
END IF
ALLOCATE(chitens(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitensInv(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitens_face(3,3,0:PP_N,0:PP_N,nSides))

! initialize
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

ALLOCATE(E(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
E=0.

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT POISSON DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation


SUBROUTINE ExactFunc(ExactFunction,x,resu,t) 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_Equation_Vars,ONLY:Pi
USE MOD_Equation_Vars,ONLY: IniCenter,IniHalfwidth,IniAmplitude
USE MOD_Equation_Vars,ONLY: ACfrequency,ACamplitude
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)              
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
REAL,INTENT(IN),OPTIONAl        :: t ! time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3)
REAL                            :: r1,r2
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
  Omega=2.*Pi*Frequency
  Resu(:)=1.+Amplitude*SIN(Omega*SUM(Cent))
CASE(31) !sinus
  Omega=2.*Pi*ACfrequency
  Resu(:)=ACamplitude*SIN(Omega*t)
CASE(32) !sinus
  resu=0.
return
  Omega=2.*Pi*ACfrequency
  Resu(:)=ACamplitude*SIN(Omega*t-Pi)
CASE(102) !linear: z=-1: 0, z=1, 1000
  resu(:)=(1+x(3))*1000.
CASE(103) ! dipole
  r1=SQRT(SUM((x(:)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  r2=SQRT(SUM((x(:)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  resu(:)=IniAmplitude*(1/r2-1/r1)
CASE(200) ! Dielectric Sphere of Radius R in constant electric field E_0 from book: 
  ! John David Jackson, Classical Electrodynamics, 3rd edition, New York: Wiley, 1999.
  ! E_0       : constant electric field in z-direction far away from sphere
  ! R         : constant radius of the sphere
  ! eps_outer : dielectric constant of surrouding medium
  ! eps_inner : dielectric constant of sphere
  !
  !   Phi_inner = - (3 / i(2 + eps_inner / eps_outer)) * E_0 * r * cos(Theta)
  !             = - (3 / i(2 + eps_inner / eps_outer)) * E_1 * z
  !  
  !   Phi_outer = ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   - 1 ) * E_0 * z
  !
  ! E = - grad(Phi)
  !
  !   E_r,inner = 0
  !   E_z,inner = (3 / (2 + eps_inner / eps_outer)) * E_0
  !  
  !   E_r,outer = 3 * ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^4 ) ) * E_0 * z
  !   E_z,inner =   ( - (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   + 1 ) * E_0
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Exactfunction not specified!')
END SELECT ! ExactFunction


END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(Ut)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation_Vars,ONLY:IniCenter,IniHalfwidth,IniAmplitude
USE MOD_Mesh_Vars,ONLY:Elem_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
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
!      !  Get source from Particles
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


SUBROUTINE CalcSourceHDG(i,j,k,iElem,resu, Phi)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_PICDepo_Vars,ONLY:PartSource,DoDeposition
USE MOD_Equation_Vars,ONLY: eps0
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation_Vars,ONLY:IniCenter,IniHalfwidth,IniAmplitude
USE MOD_Mesh_Vars,ONLY:Elem_xGP
USE MOD_Particle_Mesh_Vars, ONLY : GEO,NbrOfRegions
USE MOD_Particle_Vars, ONLY : RegionElectronRef
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
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
  source_e=0.
  IF (PRESENT(Phi)) THEN
    RegionID=0
    IF (NbrOfRegions .GT. 0) RegionID=GEO%ElemToRegion(iElem)
      IF (RegionID .NE. 0) THEN
        source_e = Phi-RegionElectronRef(2,RegionID)
        IF (source_e .LT. 0.) THEN
          source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
                   * EXP( (source_e) / RegionElectronRef(3,RegionID) )
        ELSE
          source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
                   * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
        END IF
        !source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
        !* EXP( (Phi-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
      END IF
  END IF
  resu(1)= - (PartSource(4,i,j,k,iElem)-source_e)/eps0
END IF
#endif /*PARTICLES*/

SELECT CASE (IniExactFunc)
CASE(0) ! Particles
  ! empty
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
  ! USE nr
   IMPLICIT NONE
   REAL beta, w, z                                                                                                  
   !beta = exp(gammln(z)+gammln(w)-gammln(z+w))  ! old - kind=6
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)           ! n   - kind=8
END FUNCTION beta 

SUBROUTINE FinalizeEquation()
!===================================================================================================================================
! Deallocate the vars !!!!
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars
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
SDEALLOCATE(chitens)
SDEALLOCATE(chitensInv)
SDEALLOCATE(chitens_face)
SDEALLOCATE(E)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

