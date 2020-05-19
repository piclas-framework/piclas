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
INTERFACE CalcSourceHDG
  MODULE PROCEDURE CalcSourceHDG
END INTERFACE
INTERFACE DivCleaningDamping
  MODULE PROCEDURE DivCleaningDamping
END INTERFACE
INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation, CalcSourceHDG,DivCleaningDamping
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
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for linear scalar advection')
CALL prms%CreateRealArrayOption('IniWavenumber'    , 'TODO-DEFINE-PARAMETER' , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('IniCenter'        , 'TODO-DEFINE-PARAMETER' , '0. , 0. , 0.')
CALL prms%CreateRealOption(     'IniAmplitude'     , 'TODO-DEFINE-PARAMETER' , '0.1')
CALL prms%CreateRealOption(     'IniHalfwidth'     , 'TODO-DEFINE-PARAMETER' , '0.1')
CALL prms%CreateRealOption(     'ACfrequency'      , 'TODO-DEFINE-PARAMETER' , '0.0')
CALL prms%CreateRealOption(     'ACamplitude'      , 'TODO-DEFINE-PARAMETER' , '0.0')

CALL prms%CreateIntOption(      'chitensWhichField', 'TODO-DEFINE-PARAMETER', '-1')
CALL prms%CreateRealOption(     'chitensValue'     , 'TODO-DEFINE-PARAMETER', '-1.0')
CALL prms%CreateRealOption(     'chitensRadius'    , 'TODO-DEFINE-PARAMETER', '-1.0')

CALL prms%CreateIntOption(      'AlphaShape'       , 'TODO-DEFINE-PARAMETER', '2')
CALL prms%CreateRealOption(     'r_cutoff'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Modified for curved and shape-function influence'//&
     ' (c*dt*SafetyFactor+r_cutoff)' , '1.0')

END SUBROUTINE DefineParametersEquation

SUBROUTINE InitEquation()
!===================================================================================================================================
! Init Poisson euqation system
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Preproc
USE MOD_ReadInTools        ,ONLY: GETREALARRAY,GETREAL,GETINT
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Equation_Vars
USE MOD_HDG_vars
USE MOD_Mesh_Vars          ,ONLY: nSides
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: chitensValue,chitensRadius  ! Deprecated variables, remove in future (by the end of 2017)
INTEGER                      :: chitensWhichField           ! Deprecated variables, remove in future (by the end of 2017)
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
ACfrequency = GETREAL('ACfrequency','0.0')
ACamplitude = GETREAL('ACamplitude','0.0')

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


SUBROUTINE ExactFunc(ExactFunction,x,resu,t,ElemID)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals         ,ONLY: Abort,mpiroot
USE MOD_Globals_Vars    ,ONLY: PI
USE MOD_Equation_Vars   ,ONLY: IniCenter,IniHalfwidth,IniAmplitude
USE MOD_Equation_Vars   ,ONLY: ACfrequency,ACamplitude
USE MOD_Dielectric_Vars ,ONLY: DielectricRatio,Dielectric_E_0,DielectricRadiusValue,DielectricEpsR
USE MOD_Mesh_Vars       ,ONLY: ElemBaryNGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
INTEGER,INTENT(IN),OPTIONAL     :: ElemID           ! ElemID
REAL,INTENT(IN),OPTIONAl        :: t ! time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3)
REAL                            :: r1,r2
REAL                            :: r_2D,r_3D,r_bary
REAL                            :: cos_theta
REAL                            :: eps1,eps2
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
CASE(30) !sinus: shifted by PI into the future (ACamplitude -> -1*ACamplitude)
  Omega=2.*PI*ACfrequency
  Resu(:)=-ACamplitude*SIN(Omega*t)
CASE(31) !sinus
  Omega=2.*PI*ACfrequency
  Resu(:)=ACamplitude*SIN(Omega*t)
CASE(32) !sinus
  resu=0.
return
  Omega=2.*PI*ACfrequency
  Resu(:)=ACamplitude*SIN(Omega*t-PI)
CASE(102) !linear: z=-1: 0, z=1, 1000
  resu(:)=(1+x(3))*1000.
CASE(103) ! dipole
  r1=SQRT(SUM((x(:)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  r2=SQRT(SUM((x(:)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  resu(:)=IniAmplitude*(1/r2-1/r1)
CASE(104) ! solution to Laplace's equation: Phi_xx + Phi_yy + Phi_zz = 0
  resu(1) = ( COS(x(1))+SIN(x(1)) )*( COS(x(2))+SIN(x(2)) )*( COSH(SQRT(2.0)*x(3))+SINH(SQRT(2.0)*x(3)) )
CASE(200) ! Dielectric Sphere of Radius R in constant electric field E_0 from book:
  ! John David Jackson, Classical Electrodynamics, 3rd edition, New York: Wiley, 1999.
  ! E_0       : constant electric field in z-direction far away from sphere
  ! R         : constant radius of the sphere
  ! eps_outer : dielectric constant of surrounding medium
  ! eps_inner : dielectric constant of sphere
  ! DielectricRatio = eps_inner / eps_outer (set in dielectric init)

  ! set radius and angle for DOF position x(1:3)
  r_2D   = SQRT(x(1)**2+x(2)**2)
  r_3D   = SQRT(x(1)**2+x(2)**2+x(3)**2)
  IF(r_3D.EQ.0.0)THEN
    cos_theta = 0.0
  ELSE
    cos_theta = x(3) / r_3D
  END IF
  IF(PRESENT(ElemID))THEN ! if ElemID is present, use for bary center determination versus sphere radius
    r_bary = SQRT(ElemBaryNGeo(1,ElemID)**2+ElemBaryNGeo(2,ElemID)**2+ElemBaryNGeo(3,ElemID)**2)
  ELSE
    r_bary = r_3D
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(r_bary.LT.DielectricRadiusValue)THEN ! inside sphere: DOF and element bary center
    ! Phi_inner = - (3 / (2 + eps_inner / eps_outer)) * E_1 * z
    !resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*Dielectric_E_0*x(3)
    resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*r_3D*cos_theta*Dielectric_E_0
  ELSEIF(r_bary.GE.DielectricRadiusValue)THEN ! outside sphere
    ! Phi_outer = ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   - 1 ) * E_0 * z
    resu(1:PP_nVar) =  ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       Dielectric_E_0*(DielectricRadiusValue**3/r_3D**2)*cos_theta-Dielectric_E_0 * r_3D*cos_theta

                       !( (DielectricRadiusValue**3) / (r_3D**3) ) - 1 )*(Dielectric_E_0 * x(3))
                       !( (DielectricRadiusValue**3) / ((r_2D**2+x(3)**2)**(3./2.)) ) - 1 )*(Dielectric_E_0 * x(3))
  ELSE
    IF(PRESENT(ElemID))THEN
      SWRITE(*,*) "ElemID                ",ElemID
      SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID),r_bary
    END IF
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(&
    __STAMP__&
    ,'Dielectric sphere. Invalid radius for exact function!')
  END IF

  ! varphi = ATAN2(x(2),x(1)) ! only needed for the electric field
  !   E_r,inner = 0
  !   E_z,inner = (3 / (2 + eps_inner / eps_outer)) * E_0
  !
  !   E_r,outer = 3 * ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^4 ) ) * E_0 * z
  !   E_z,inner =   ( - (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   + 1 ) * E_0
CASE(300) ! Dielectric Slab in z-direction of half width R in constant electric field E_0: adjusted from CASE(200)
  ! R = DielectricRadiusValue
  ! DielectricRatio = eps/eps0

  ! for BC, not ElemID will be given
  IF(.NOT.PRESENT(ElemID))THEN
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue   *((DielectricRatio-1.)/(DielectricRatio))/(abs(x(3))) + 1)
    RETURN
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(ABS(ElemBaryNGeo(3,ElemID)).LT.DielectricRadiusValue)THEN ! inside box: DOF and element bary center
    ! Phi_inner = ?

    ! marcel
    !resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(1 - (DielectricRatio-1)/(DielectricRatio+2))
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-((DielectricRatio-1.)/(DielectricRatio)) + 1)

    ! linear
    !resu(1:PP_nVar) = -(1./(DielectricRatio))*Dielectric_E_0*x(3)

    ! from sphere
    ! Phi_inner = - (3 / (2 + eps_inner / eps_outer)) * E_1 * z
    !resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*x(3)*Dielectric_E_0



  ELSEIF(ABS(ElemBaryNGeo(3,ElemID)).GT.DielectricRadiusValue)THEN ! outside sphere
    ! Phi_outer = ?
    !resu(1:PP_nVar) =( ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       !( (DielectricRadiusValue**3) / (x(3)**(3.)) ) - 1 )*(Dielectric_E_0 * x(3))
                       !( (DielectricRadiusValue**3) / ((x(3)**2)**(3./2.)) ) - 1 )*(Dielectric_E_0 * x(3))

    ! marcel
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue**2*((DielectricRatio-1.)/(DielectricRatio))/(x(3)**2) + 1)

    ! linear
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue   *((DielectricRatio-1.)/(DielectricRatio))/(abs(x(3))) + 1)

    !resu(1:PP_nVar) = -Dielectric_E_0*(x(3) - sign(1.0,x(3))*((DielectricRatio-1)/(DielectricRatio+2))*((2**3)/(x(3)**2)) )
    !resu(1:PP_nVar) =( ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       !( (2.0**3) / ((x(3)**3) ) - 1 )*(Dielectric_E_0 * x(3))
  ELSE
    SWRITE(*,*) "ElemID                ",ElemID
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(&
    __STAMP__&
    ,'Dielectric sphere. Invalid radius for exact function!')
  END IF
CASE(301) ! like CASE=300, but only in positive z-direction the dielectric region is assumed
  ! R = DielectricRadiusValue
  ! DielectricRatio = eps/eps0

  ! for BC, not ElemID will be given
  IF(.NOT.PRESENT(ElemID))THEN
    IF(x(3).GT.0.0)THEN ! inside dielectric
      resu(1:PP_nVar) = -(Dielectric_E_0/DielectricRatio)*(x(3)-DielectricRadiusValue)
    ELSE
      resu(1:PP_nVar) = -(Dielectric_E_0)*(x(3)-DielectricRadiusValue/DielectricRatio)
    END IF
    RETURN
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(     (ABS(ElemBaryNGeo(3,ElemID)).LT.2.0*DielectricRadiusValue).AND.(ElemBaryNGeo(3,ElemID).GT.0.0))THEN ! inside box: DOF and element bary center
    resu(1:PP_nVar) = -(Dielectric_E_0/DielectricRatio)*(x(3)-DielectricRadiusValue)
  ELSEIF( (ABS(ElemBaryNGeo(3,ElemID)).GT.DielectricRadiusValue).OR.(ElemBaryNGeo(3,ElemID).LT.0.0) )THEN ! outside sphere
    resu(1:PP_nVar) = -(Dielectric_E_0)*(x(3)-DielectricRadiusValue/DielectricRatio)
  ELSE
    SWRITE(*,*) "ElemID                ",ElemID
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(&
    __STAMP__&
    ,'Dielectric sphere. Invalid radius for exact function!')
  END IF

CASE(400) ! Point Source in Dielectric Region with epsR_1  = 1 for x < 0 (vacuum)
  !                                                epsR_2 != 1 for x > 0 (dielectric region)
  ! DielectricRadiusValue is used as distance between dielectric interface and position of chargeed point particle
  ! set radius and angle for DOF position x(1:3)
  ! Limitations:
  ! only valid for eps_2 = 1
  ! and q = 1
  r_2D   = SQRT(x(1)**2+x(2)**2)
  r1 = SQRT(r_2D**2 + (DielectricRadiusValue-x(3))**2)
  r2 = SQRT(r_2D**2 + (DielectricRadiusValue+x(3))**2)

  eps2=1.0
  eps1=DielectricEpsR

  IF(x(3).GT.0.0)THEN
    IF(ALL((/ x(1).EQ.0.0,  x(2).EQ.0.0, x(3).EQ.DielectricRadiusValue /)))THEN
      print*, "HERE?!?!?!"
    END IF
    IF((r1.LE.0.0).OR.(r2.LE.0.0))THEN
      SWRITE(*,*) "r1=",r1
      SWRITE(*,*) "r2=",r2
      CALL abort(&
          __STAMP__&
          ,'ExactFunc=400: Point source in dielectric region. Cannot evaluate the exact function at the singularity!')
    END IF
    resu(1:PP_nVar) = (1./eps1)*(&
                                   1./r1 + ((eps1-eps2)/(eps1+eps2))*&
                                   1./r2 )/(4*PI)
  ELSE
    IF(r1.LE.0.0)THEN
      SWRITE(*,*) "r1=",r1
      CALL abort(&
          __STAMP__&
          ,'Point source in dielectric region: Cannot evaluate the exact function at the singularity!')
    END IF
    resu(1:PP_nVar) = (2./(eps2+eps1)) * 1./r1 /(4*PI)
  END IF

CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Exactfunction not specified!')
END SELECT ! ExactFunction


END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(Ut)
!===================================================================================================================================
!
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


SUBROUTINE CalcSourceHDG(i,j,k,iElem,resu, Phi, warning_linear)
!===================================================================================================================================
! Determine the right-hand-side of Poisson's equation (either by an analytic function or deposition of charge from particles)
! TODO: currently particles are enforced, which means that they over-write the exact function solution because
! the combination of both has not been specified
! How should this function work???
! for dielectric regions DO NOT apply the scaling factor Eps_R here (which is const. in HDG due to current implementation) because
! it is in the tensor "chitens"
!===================================================================================================================================
! MODULES
USE MOD_Globals            ,ONLY: Abort
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
#ifdef PARTICLES
USE MOD_PICDepo_Vars       ,ONLY: PartSource,DoDeposition,DepositionType
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,NbrOfRegions
USE MOD_Particle_Vars      ,ONLY: RegionElectronRef
USE MOD_Globals_Vars       ,ONLY: eps0
#if IMPA
USE MOD_LinearSolver_Vars  ,ONLY: ExplicitPartSource
#endif
#endif /*PARTICLES*/
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Equation_Vars      ,ONLY: IniCenter,IniHalfwidth,IniAmplitude
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: i, j, k,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
LOGICAL,INTENT(INOUT),OPTIONAL  :: warning_linear
REAL,INTENT(IN),OPTIONAL        :: Phi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: x(3)
REAL                            :: r1,r2, source_e
REAL,DIMENSION(3)               :: dx1,dx2,dr1dx,dr2dx,dr1dx2,dr2dx2
INTEGER                         :: RegionID
!===================================================================================================================================
! Calculate IniExactFunc before particles are superimposed, because the IniExactFunc might be needed by the CalcError function
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
  resu=0. ! empty
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

#ifdef PARTICLES
IF(DoDeposition .OR. (TRIM(DepositionType).EQ.'constant'))THEN
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
        warning_linear = .TRUE.
      END IF
      !source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
      !* EXP( (Phi-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
    END IF
  END IF
#if IMPA
  resu(1)= - (PartSource(4,i,j,k,iElem)+ExplicitPartSource(4,i,j,k,iElem)-source_e)/eps0
#else
  resu(1)= - (PartSource(4,i,j,k,iElem)-source_e)/eps0
#endif
END IF
#endif /*PARTICLES*/

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

