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
CALL prms%CreateIntOption(      'IniExactFunc'     , 'Define exact function necessary for '//&
                                                     'discrete velocity method', '-1')
CALL prms%CreateRealOption(     'DVM-omegaVHS',      'Variable Hard Sphere parameter')
CALL prms%CreateRealOption(     'DVM-T_Ref',         'VHS reference temperature')
CALL prms%CreateRealOption(     'DVM-d_Ref',         'VHS reference diameter')
CALL prms%CreateRealOption(     'DVM-Mass',          'Molecular mass')
CALL prms%CreateIntOption(      'DVM-Internal_DOF',  'Number of (non translational) internal degrees of freedom', '0')
CALL prms%CreateIntOption(      'DVM-Dimension',     'Number of space dimensions for velocity discretization', '3')
CALL prms%CreateIntOption(      'DVM-BGKCollModel',  'Select the BGK method:\n'//&
                                                     '1: Ellipsoidal statistical (ESBGK)\n'//&
                                                     '2: Shakov (SBGK)\n'//&
                                                     '3: Standard BGK (Maxwell)'//&
                                                     '4: Conservative Maxwell)')
CALL prms%CreateIntOption(      'DVM-Method',        'Select the DVM model:\n'//&
                                                     '1: Exponential differencing (EDDVM)\n'//&
                                                     '2: DUGKS')
CALL prms%CreateIntOption(      'DVM-VeloDiscretization',      '1: do not use, 2: Gauss-Hermite, 3: Newton-Cotes', '2')
CALL prms%CreateRealArrayOption('DVM-GaussHermiteTemp',        'Reference temperature for GH quadrature (per direction)',&
                                                               '(/273.,273.,273./)')
CALL prms%CreateRealArrayOption('DVM-VeloMin',                 'Only for Newton-Cotes velocity quadrature', '(/-1.,-1.,-1./)')
CALL prms%CreateRealArrayOption('DVM-VeloMax',                 'Only for Newton-Cotes velocity quadrature', '(/1.,1.,1./)')
CALL prms%CreateIntOption(      'DVM-nVelo' ,                  'Number of velocity discretization points', '15')
CALL prms%CreateIntArrayOption( 'DVM-NewtonCotesDegree',       'Degree of the subquadrature for composite quadrature', '(/1,1,1/)')
CALL prms%CreateIntOption(      'IniRefState',  'Refstate required for initialization.')
CALL prms%CreateRealArrayOption('RefState',     'State(s) in primitive variables (density, velo, temp, press, heatflux).',&
                                                 multiple=.TRUE., no=14 )
CALL prms%CreateRealArrayOption('DVM-Accel',    'Acceleration vector for force term', '(/0., 0., 0./)')
CALL prms%CreateLogicalOption(  'DVM-WriteMacroSurfaceValues',  'Surface output', '.FALSE.')
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY : BoltzmannConst
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_ReadInTools,        ONLY:GETREALARRAY,GETINTARRAY,GETREAL,GETINT,GETLOGICAL, CountOption
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Basis              ,ONLY: LegendreGaussNodesAndWeights, GaussHermiteNodesAndWeights, NewtonCotesNodesAndWeights
USE MOD_Equation_Vars_FV
USE MOD_DVM_Boundary_Analyze,ONLY: InitDVMBoundaryAnalyze
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
DVMSpeciesData%Internal_DOF = GETINT('DVM-Internal_DOF')
DVMSpeciesData%Mass = GETREAL('DVM-Mass')
DVMBGKModel = GETINT('DVM-BGKCollModel')
DVMMethod = GETINT('DVM-Method')
DVMVeloDisc = GETINT('DVM-VeloDiscretization')
DVMSpeciesData%mu_Ref = 30.*SQRT(DVMSpeciesData%Mass*BoltzmannConst*DVMSpeciesData%T_Ref/Pi) &
                      /(4.*(4.-2.*DVMSpeciesData%omegaVHS)*(6.-2.*DVMSpeciesData%omegaVHS) &
                      * DVMSpeciesData%d_Ref**2.)
DVMSpeciesData%R_S = BoltzmannConst / DVMSpeciesData%Mass
DVMDim = GETINT('DVM-Dimension')
IF ((DVMDim.GT.3).OR.(DVMDim.LT.1)) CALL abort(__STAMP__,'DVM error: dimension must be between 1 and 3')
DVMSpeciesData%Prandtl =2.*(DVMSpeciesData%Internal_DOF + 5.)/(2.*DVMSpeciesData%Internal_DOF + 15.)

DVMnVelos(1:3) = 1
DVMnVelos(1:DVMDim) = GETINT('DVM-nVelo')
PP_nVar_FV = (DVMnVelos(1)**DVMDim)
IF (DVMDim.LT.3) PP_nVar_FV = PP_nVar_FV*2 !double variables for reduced distributions

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
  ALLOCATE(RefState(14,nRefState))
  DO i=1,nRefState
    RefState(1:14,i)  = GETREALARRAY('RefState',14)
  END DO
END IF

DVMForce = GETREALARRAY('DVM-Accel',3)

ALLOCATE(DVMMomentSave(15,nElems))
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

WriteDVMSurfaceValues = GETLOGICAL('DVM-WriteMacroSurfaceValues')
IF (WriteDVMSurfaceValues) CALL InitDVMBoundaryAnalyze()

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
USE MOD_Equation_Vars_FV, ONLY: DVMSpeciesData, RefState, DVMBGKModel
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
REAL                            :: MacroVal(14)
REAL                            :: WallVelo1, WallTemp1, WallVelo2, WallTemp2
REAL                            :: SecondDist(PP_nVar_FV)
REAL                            :: SodMacro_L(5), SodMacro_R(5), SodMacro_LL(5), SodMacro_M(5), SodMacro_RR(5)
REAL                            :: gamma, Ggamma, beta, pL, pR, pM, cL, cR, cM, vs
REAL                            :: mu, tau
!==================================================================================================================================

Resu   =0.

SELECT CASE (ExactFunction)
CASE(0)
  Resu=0.

CASE(1) !Grad 13 uniform init
  CALL GradDistribution(RefState(:,1),Resu(:))

CASE(2) ! couette flow, L=1m between walls, RefState: 1 -> initial state, 2 -> y=-0.5 boundary, 3 -> y=0.5 boundary
  CALL GradDistribution(RefState(:,1),Resu(:))
  ! steady state
  IF (tIn.GT.0.) THEN
    MacroVal(:) = RefState(:,1)
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

  IF (tIn.EQ.0.) THEN ! initial state
    IF (x(1).LT.0.) THEN
      CALL MaxwellDistribution(RefState(:,1),Resu(:))
    ELSE
      CALL MaxwellDistribution(RefState(:,2),Resu(:))
    END IF
  ELSE ! analytical solution
    MacroVal = 0.
    SodMacro_L=RefState(1:5,1)
    SodMacro_R=RefState(1:5,2)
    SodMacro_LL = 0.         !     | L | LL | M | RR | R |
    SodMacro_M = 0.
    SodMacro_RR = 0.
    gamma = 5./3. !monatomic gas
    Ggamma = (gamma-1)/(gamma+1)
    beta = (gamma-1)/gamma/2.
    pL = DVMSpeciesData%R_S*SodMacro_L(1)*SodMacro_L(5)
    pR = DVMSpeciesData%R_S*SodMacro_R(1)*SodMacro_R(5)
    pM=(pL+pR)/2.
    CALL SecantSod(pM, pL, pR, SodMacro_L(1), SodMacro_R(1), gamma, Ggamma, beta, 1e-15, 100)
    cL = sqrt(gamma*DVMSpeciesData%R_S*SodMacro_L(5))
    cR = sqrt(gamma*DVMSpeciesData%R_S*SodMacro_R(5))
    cM = cL * (pM/pL)**beta
    SodMacro_M(1) = SodMacro_L(1)*(pM/pL)**(1./gamma)
    SodMacro_M(2) = (cL-cM)*2/(gamma-1)
    SodMacro_M(5) = pM/DVMSpeciesData%R_S/SodMacro_M(1)
    SodMacro_RR(1) = SodMacro_R(1)*(pM+Ggamma*pR)/(pR+Ggamma*pM)
    SodMacro_RR(2) = SodMacro_M(2)
    SodMacro_RR(5) = pM/DVMSpeciesData%R_S/SodMacro_RR(1)
    vs = cR*sqrt((beta/Ggamma)*(pM/pR + Ggamma))
    IF (x(1).LT.(-tIn*cL)) THEN
      MacroVal(1:5) = SodMacro_L
    ELSE IF (x(1).LT.(tIn*(SodMacro_M(2)-cM))) THEN
      SodMacro_LL(2) = 2./(gamma+1.) * (cL + x(1)/tIn)
      SodMacro_LL(1) = SodMacro_L(1) * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2./(gamma-1.))
      SodMacro_LL(5) = pL * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2*gamma/(gamma-1))/DVMSpeciesData%R_S/SodMacro_LL(1)
      MacroVal(1:5) = SodMacro_LL
    ELSE IF (x(1).LT.(tIn*SodMacro_M(2))) THEN
      MacroVal(1:5) = SodMacro_M
    ELSE IF (x(1).LT.(tIn*vs)) THEN
      MacroVal(1:5) = SodMacro_RR
    ELSE
      MacroVal(1:5) = SodMacro_R
    END IF
    CALL MaxwellDistribution(MacroVal,Resu(:))
  END IF

CASE(4) !heat flux relaxation test case
  CALL GradDistribution(RefState(:,1),Resu(:))
  IF (tIn.GT.0.) THEN ! relaxation
    MacroVal(:) = RefState(:,1)
    mu = DVMSpeciesData%mu_Ref*(MacroVal(5)/DVMSpeciesData%T_Ref)**(DVMSpeciesData%omegaVHS+0.5)
    tau = mu/(DVMSpeciesData%R_S*MacroVal(1)*MacroVal(5))
    MacroVal(12:14) = MacroVal(12:14)*EXP(-tIn*DVMSpeciesData%Prandtl/tau) !Heat flux relaxes with rate Pr/tau
    CALL GradDistribution(MacroVal(:),Resu(:))
  END IF

CASE(5) !Taylor-Green vortex
  MacroVal(:) = RefState(:,1)
  MacroVal(2) = RefState(2,1)*SIN(x(1))*COS(x(2))*COS(x(3))
  MacroVal(3) = -RefState(2,1)*COS(x(1))*SIN(x(2))*COS(x(3))
  MacroVal(4) = 0.
  MacroVal(5) = MacroVal(5)+(RefState(2,1)**2)/(16.*DVMSpeciesData%R_S)*(COS(2.*x(1))+COS(2.*x(2)))*(COS(2.*x(3))+2.)
  CALL MaxwellDistribution(MacroVal,Resu(:))

CASE(6) !Sum of 2 distributions
  CALL GradDistribution(RefState(:,1),Resu(:))
  CALL GradDistribution(RefState(:,2),SecondDist(:))
  Resu(:) = Resu(:) + SecondDist(:)

CASE DEFAULT
  CALL abort(__STAMP__,&
             'Specified exact function not implemented!')
END SELECT ! ExactFunction

END SUBROUTINE ExactFunc


SUBROUTINE SecantSod(pM, pL, pR, rhoL, rhoR, gamma, Ggamma, beta, eps, maxit)
!==================================================================================================================================
!> Secant method to get Sod shock middle state
!==================================================================================================================================
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: pL, pR, rhoL, rhoR, gamma, Ggamma, beta, eps
INTEGER,INTENT(IN)  :: maxit
REAL,INTENT(INOUT)  :: pM
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: u3, u4, pMold, pMnew, diffold
!==================================================================================================================================
pMold = 2*pM
u3 = (pL**beta-pMold**beta)*SQRT(((1.-Ggamma**2)*pL**(1./gamma))/(Ggamma**2 * rhoL))
u4 = (pMold-pR)*SQRT((1.-Ggamma)/(rhoR*(pMold+Ggamma*pR)))
diffold = u3-u4

DO i=1, maxit
  u3 = (pL**beta-pM**beta)*SQRT(((1.-Ggamma**2)*pL**(1./gamma))/(Ggamma**2 * rhoL))
  u4 = (pM-pR)*SQRT((1.-Ggamma)/(rhoR*(pM+Ggamma*pR)))
  pMnew = pM - (u3-u4)*(pM-pMold)/(u3-u4-diffold)
  pMold = pM
  pM = pMnew
  diffold = u3-u4
  IF (ABS(pM-pMold).LT.eps) RETURN
END DO

SWRITE(UNIT_stdOut,*) 'NewtonSod: max number of iterations reached'

END SUBROUTINE SecantSod


SUBROUTINE CalcSource(t,coeff,Ut)
!==================================================================================================================================
! Dummy
!==================================================================================================================================
! MODULES
USE MOD_Globals           ,ONLY: abort
USE MOD_Globals_Vars      ,ONLY: PI,eps0
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,coeff
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar_FV,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!==================================================================================================================================


END SUBROUTINE CalcSource


!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars_FV,ONLY:EquationInitIsDone, DVMVelos, DVMWeights, RefState, WriteDVMSurfaceValues
USE MOD_DVM_Boundary_Analyze,ONLY: FinalizeDVMBoundaryAnalyze
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone = .FALSE.
IF (WriteDVMSurfaceValues) CALL FinalizeDVMBoundaryAnalyze()
SDEALLOCATE(DVMVelos)
SDEALLOCATE(DVMWeights)
SDEALLOCATE(RefState)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation_FV
