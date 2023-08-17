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

MODULE MOD_Equation_FV
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
INTERFACE InitEquation_FV
  MODULE PROCEDURE InitEquation
END INTERFACE
INTERFACE ExactFunc_FV
  MODULE PROCEDURE ExactFunc
END INTERFACE
INTERFACE CalcSource_FV
  MODULE PROCEDURE CalcSource
END INTERFACE
INTERFACE FinalizeEquation_FV
  MODULE PROCEDURE FinalizeEquation
END INTERFACE
INTERFACE DefineParametersEquation_FV
  MODULE PROCEDURE DefineParametersEquation
END INTERFACE

PUBLIC::InitEquation_FV,ExactFunc_FV,FinalizeEquation_FV,CalcSource_FV,DefineParametersEquation_FV
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for equation
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection', '-1')
CALL prms%CreateRealOption(     'DVM-omegaVHS',        "Blabla.")
CALL prms%CreateRealOption(     'DVM-T_Ref',        "Blabla.")
CALL prms%CreateRealOption(     'DVM-d_Ref',        "Blabla.")
CALL prms%CreateRealOption(     'DVM-Internal_DOF',        "Blabla.")
CALL prms%CreateRealOption(     'DVM-Mass',        "Blabla.")
CALL prms%CreateIntOption(     'DVM-BGKModel',        "1: Maxwell, 2: Shakhov", "2")
CALL prms%CreateIntOption(     'DVM-Method',        "1: Exponential differencing, 2: DUGKS", "1")
CALL prms%CreateIntOption(     'DVM-VeloDiscretization',        "1: don't use, 2: Gauss-Hermite, 3: Newton-Cotes", "2")
CALL prms%CreateRealArrayOption(     'DVM-GaussHermiteTemp',        "Blabla.","(/273.,273.,273./)")
CALL prms%CreateRealArrayOption(     'DVM-VeloMin',        "Only for Newton-Cotes velocity quadrature", "(/-1.,-1.,-1./)")
CALL prms%CreateRealArrayOption(     'DVM-VeloMax',        "Only for Newton-Cotes velocity quadrature", "(/1.,1.,1./)")
CALL prms%CreateIntOption(     'DVM-Dimension',        "Blabla.", "3")
CALL prms%CreateIntOption( 'DVM-nVelo' , "Number of velocity discretization points", '15')
CALL prms%CreateIntArrayOption(     'DVM-NewtonCotesDegree',        "Degree of the subquadrature for composite quadrature", "(/1,1,1/)")
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE., no=8 )
CALL prms%CreateRealArrayOption('DVM-Accel', "Acceleration vector for force term", "(/0., 0., 0./)")
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters (advection velocity, diffusion coeff, exact function)  from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY : BoltzmannConst
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_ReadInTools,        ONLY:GETREALARRAY,GETINTARRAY,GETREAL,GETINT, CountOption
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Basis              ,ONLY: LegendreGaussNodesAndWeights, GaussHermiteNodesAndWeights, NewtonCotesNodesAndWeights
USE MOD_Equation_Vars_FV
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i, iGH, iDim
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLinearScalarAdvection not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DVM...'

Pi=ACOS(-1.)
IniExactFunc_FV = GETINT('IniExactFunc')

DVMSpeciesData%omegaVHS = GETREAL('DVM-omegaVHS')
DVMSpeciesData%T_Ref = GETREAL('DVM-T_Ref')
DVMSpeciesData%d_Ref = GETREAL('DVM-d_Ref')
DVMSpeciesData%Internal_DOF = GETREAL('DVM-Internal_DOF')
DVMSpeciesData%Mass = GETREAL('DVM-Mass')
DVMBGKModel = GETINT('DVM-BGKModel')
DVMMethod = GETINT('DVM-Method')
DVMVeloDisc = GETINT('DVM-VeloDiscretization')
DVMSpeciesData%mu_Ref = 30.*SQRT(DVMSpeciesData%Mass*BoltzmannConst*DVMSpeciesData%T_Ref/Pi) &
                      /(4.*(4.-2.*DVMSpeciesData%omegaVHS)*(6.-2.*DVMSpeciesData%omegaVHS) &
                      * DVMSpeciesData%d_Ref**2.)
DVMSpeciesData%R_S = BoltzmannConst / DVMSpeciesData%Mass
DVMDim = GETINT('DVM-Dimension')
IF ((DVMDim+DVMSpeciesData%Internal_DOF).LT.3.) CALL abort(__STAMP__,&
           'Degrees of freedom need to be at least 3 - number of dimensions')
! DVMSpeciesData%Prandtl =2.*(DVMSpeciesData%Internal_DOF + 5.)/(2.*DVMSpeciesData%Internal_DOF + 15.)
DVMSpeciesData%Prandtl = 2./3.
DVMnVelos(1:3) = 1
IF (DVMSpeciesData%Internal_DOF.GT.0.0) THEN
  DVMnVelos(1:DVMDim) = GETINT('DVM-nVelo')
  PP_nVar_FV = (DVMnVelos(1)**DVMDim)*2
ELSE
  DVMnVelos(:) = GETINT('DVM-nVelo')
  PP_nVar_FV = (DVMnVelos(1)**3)
END IF

ALLOCATE(DVMVelos(MAXVAL(DVMnVelos),3), DVMWeights(MAXVAL(DVMnVelos),3))
DVMVelos(:,:)=0.
DVMWeights(:,:)=1.
IF (DVMVeloDisc .EQ. 1) THEN ! Gauss-Legendre Nodes and Weights
  DVMVeloMin(:) = GETREALARRAY('DVM-VeloMin',3)
  DVMVeloMax(:) = GETREALARRAY('DVM-VeloMax',3)
  DO iDim=1,DVMDim
    CALL LegendreGaussNodesAndWeights(DVMnVelos(iDim)-1,DVMVelos(:,iDim),DVMWeights(:,iDim))
    DVMVelos(:,iDim) = DVMVelos(:,iDim)*ABS(DVMVeloMax(iDim)-DVMVeloMin(iDim))/2. &
                     + ABS(DVMVeloMax(iDim)-DVMVeloMin(iDim))/2. + DVMVeloMin(iDim)
    DVMWeights(:,iDim) = DVMWeights(:,iDim)*ABS(DVMVeloMax(iDim)-DVMVeloMin(iDim))/2.
  END DO
ELSE IF (DVMVeloDisc .EQ. 2) THEN ! Newton-Cotes Nodes and Weights
  DVMVeloMin(:) = GETREALARRAY('DVM-VeloMin',3)
  DVMVeloMax(:) = GETREALARRAY('DVM-VeloMax',3)
  DVMNewtDeg(:) = GETINTARRAY('DVM-NewtonCotesDegree',3)
  DO iDim=1,DVMDim
    CALL NewtonCotesNodesAndWeights(DVMnVelos(iDim)-1,DVMNewtDeg(iDim),DVMVelos(:,iDim),DVMWeights(:,iDim))
    DVMVelos(:,iDim) = DVMVelos(:,iDim)*ABS(DVMVeloMax(iDim)-DVMVeloMin(iDim)) + DVMVeloMin(iDim)
    DVMWeights(:,iDim) = DVMWeights(:,iDim)*ABS(DVMVeloMax(iDim)-DVMVeloMin(iDim))
  END DO
ELSE IF (DVMVeloDisc .EQ. 3) THEN ! Gauss-Hermite Nodes and Weights
  DVMGHTemp(:) = GETREALARRAY('DVM-GaussHermiteTemp',3)
  DO iDim=1,DVMDim
    CALL GaussHermiteNodesAndWeights(DVMnVelos(iDim)-1,DVMVelos(:,iDim),DVMWeights(:,iDim))
    DO iGH=1,DVMnVelos(iDim)
      DVMWeights(iGH,iDim) = DVMWeights(iGH,iDim)/(EXP(-DVMVelos(iGH,iDim)**2.)) &
                           * SQRT(2.*DVMSpeciesData%R_S*DVMGHTemp(iDim))
    END DO
    DVMVelos(:,iDim) = DVMVelos(:,iDim)*SQRT(2.*DVMSpeciesData%R_S*DVMGHTemp(iDim))
  END DO
ELSE
  CALL abort(__STAMP__,&
             'Undefined velocity space discretization')
END IF ! DVMVeloDisc

! Read Boundary information / RefStates / perform sanity check
IniRefState = GETINT('IniRefState')
nRefState=CountOption('RefState')
IF(IniRefState.GT.nRefState)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))
END IF

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefState(8,nRefState))
  DO i=1,nRefState
    RefState(1:8,i)  = GETREALARRAY('RefState',8)
  END DO
END IF

DVMForce = GETREALARRAY('DVM-Accel',3)

ALLOCATE(DVMMomentSave(3,nElems))
DVMMomentSave = 0.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

SWRITE(UNIT_stdOut,*)'DVM uses ', DVMnVelos(1), 'x', DVMnVelos(2), 'x', DVMnVelos(3),' velocities!'
SELECT CASE (DVMMethod)
CASE(1)
  SWRITE(UNIT_stdOut,*)'Method: Exponential Differencing DVM'
CASE(2)
  SWRITE(UNIT_stdOut,*)'Method: DUGKS'
END SELECT


EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DVM DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquation


SUBROUTINE ExactFunc(ExactFunction,tIn,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars,  ONLY: PI
USE MOD_DistFunc,      ONLY: MaxwellDistribution, MacroValuesFromDistribution, GradDistribution
USE MOD_Equation_Vars_FV, ONLY: DVMSpeciesData, RefState
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of
                                                          !< timestep if full boundary order is used (only with RK3)
INTEGER, INTENT(IN)                :: tDeriv
REAL,INTENT(IN)                 :: x(3)                   !< coordinates to evaluate exact function
INTEGER,INTENT(IN)              :: ExactFunction          !< specifies the exact function to be used
REAL,INTENT(OUT)                :: Resu(PP_nVar_FV)          !< output state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: MacroVal(8)
REAL                            :: WallVelo1, WallTemp1, WallVelo2, WallTemp2
!==================================================================================================================================

Resu   =0.

SELECT CASE (ExactFunction)
CASE(0)
  Resu=0.

CASE(1,4) !Grad 13 uniform init (only heat flux for now); 4 = heat flux relaxation test case
  MacroVal(:) = RefState(:,1)
  CALL GradDistribution(MacroVal,Resu(:))

CASE(2) ! couette flow, L=1m between walls, RefState: 1 -> initial state, 2 -> y=-0.5 boundary, 3 -> y=0.5 boundary
  MacroVal(:) = RefState(:,1)
  CALL GradDistribution(MacroVal,Resu(:))
  ! steady state
  IF (tIn.GT.0.) THEN
    WallVelo1 = RefState(2,2)
    WallTemp1 = RefState(5,2)
    WallVelo2 = RefState(2,3)
    WallTemp2 = RefState(5,3)
    MacroVal(2) = WallVelo2 + (WallVelo1-WallVelo2)*(0.5-x(2))
    MacroVal(5) = WallTemp2 + (0.5-x(2))*(WallTemp1-WallTemp2+((WallVelo1-WallVelo2)**2)*(0.5+x(2))*4/15/DVMSpeciesData%R_S/2)
                                                                                                    ! cf Eucken's relation
    CALL MaxwellDistribution(MacroVal,Resu(:))
  END IF

CASE(3) !sod shock
  Resu=0.

  IF (x(1).LT.0.) THEN
    CALL MaxwellDistribution(RefState(:,1),Resu(:))
  ELSE
    CALL MaxwellDistribution(RefState(:,2),Resu(:))
  END IF
  ! IF (x(1).LT.(-tIn*cL)) THEN
  !   CALL MaxwellDistribution(SodMacro_L,Resu(:))
  ! ELSE IF (x(1).LT.(tIn*(SodMacro_M(2)-cM))) THEN
  !   SodMacro_LL(2) = 2./(gamma+1.) * (cL + x(1)/tIn)
  !   SodMacro_LL(1) = SodMacro_L(1) * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2./(gamma-1.))
  !   SodMacro_LL(5) = pL * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2*gamma/(gamma-1))/DVMSpeciesData%R_S/SodMacro_LL(1)
  !   CALL MaxwellDistribution(SodMacro_LL,Resu(:))
  ! ELSE IF (x(1).LT.(tIn*SodMacro_M(2))) THEN
  !   CALL MaxwellDistribution(SodMacro_M,Resu(:))
  ! ELSE IF (x(1).LT.(tIn*vs)) THEN
  !   CALL MaxwellDistribution(SodMacro_RR,Resu(:))
  ! ELSE
  !   CALL MaxwellDistribution(SodMacro_R,Resu(:))
  ! END IF

CASE(5) !Taylor-Green vortex
  MacroVal(:) = RefState(:,1)
  MacroVal(2) = MacroVal(2)*SIN(2*PI*x(1))*COS(2*PI*x(2))*COS(2*PI*x(3))
  MacroVal(3) = -MacroVal(3)*COS(2*PI*x(1))*SIN(2*PI*x(2))*COS(2*PI*x(3))
  MacroVal(4) = 0.
  CALL MaxwellDistribution(MacroVal,Resu(:))

CASE DEFAULT
  CALL abort(__STAMP__,&
             'Specified exact function not implemented!')
END SELECT ! ExactFunction

END SUBROUTINE ExactFunc

SUBROUTINE CalcSource(t,coeff,Ut)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals           ,ONLY: abort
USE MOD_Globals_Vars      ,ONLY: PI,eps0
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar_FV,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================


END SUBROUTINE CalcSource

!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars_FV,ONLY:EquationInitIsDone, DVMVelos, DVMWeights, RefState
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone = .FALSE.
SDEALLOCATE(DVMVelos)
SDEALLOCATE(DVMWeights)
SDEALLOCATE(RefState)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation_FV
