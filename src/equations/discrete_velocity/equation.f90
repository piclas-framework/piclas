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
CALL prms%CreateIntOption(      'IniExactFunc-FV'     , 'Define exact function necessary for '//&
                                                     'discrete velocity method', '-1')
CALL prms%CreateIntOption(      'DVM-nSpecies',      'Number of species for DVM', '1')
CALL prms%CreateRealOption(     'DVM-Species[$]-omegaVHS',      'Variable Hard Sphere parameter', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DVM-Species[$]-T_Ref',         'VHS reference temperature', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DVM-Species[$]-d_Ref',         'VHS reference diameter', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DVM-Species[$]-Mass',          'Molecular mass', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DVM-Species[$]-Charge',          'Electrical charge', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DVM-Species[$]-Internal_DOF',  'Number of (non translational) internal degrees of freedom', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DVM-Species[$]-VeloDiscretization',      '1: do not use, 2: Gauss-Hermite, 3: Newton-Cotes', '2', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('DVM-Species[$]-GaussHermiteTemp',        'Reference temperature for GH quadrature (per direction)',&
                                                               '(/273.,273.,273./)', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('DVM-Species[$]-VeloMin',                 'Only for Newton-Cotes velocity quadrature', '(/-1.,-1.,-1./)', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('DVM-Species[$]-VeloMax',                 'Only for Newton-Cotes velocity quadrature', '(/1.,1.,1./)', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DVM-Species[$]-nVelo' ,                  'Number of velocity discretization points', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DVM-Species[$]-NewtonCotesDegree',       'Degree of the subquadrature for composite quadrature', '(/1,1,1/)', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DVM-Dimension',     'Number of space dimensions for velocity discretization', '3')
CALL prms%CreateIntOption(      'DVM-BGKCollModel',  'Select the BGK method:\n'//&
                                                     '1: Ellipsoidal statistical (ESBGK)\n'//&
                                                     '2: Shakov (SBGK)\n'//&
                                                     '3: Standard BGK (Maxwell)'//&
                                                     '4: Conservative ESBGK'//&
                                                     '5: Conservative Maxwell'//&
                                                     '6: SkewNormal BGK (SNBGK)'//&
                                                     '7: Grad 13 BGK')
CALL prms%CreateIntOption(      'DVM-Method',        'Select the DVM model:\n'//&
                                                     '1: Exponential differencing (EDDVM)\n'//&
                                                     '2: DUGKS')
CALL prms%CreateIntOption(      'IniRefState-FV',  'Refstate required for initialization.')
CALL prms%CreateRealArrayOption('RefState-FV',     'State(s) in primitive variables (density, velo, temp, press, heatflux).',&
                                                 multiple=.TRUE., no=14 )
CALL prms%CreateRealArrayOption('DVM-Accel',    'Acceleration vector for force term', '(/0., 0., 0./)')
CALL prms%CreateLogicalOption(  'DVM-Collisions',  'Activate collision (RHS BGK equation)', '.TRUE.')
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
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i, iGH, iDim, iSpec
CHARACTER(32)         :: hilf
CHARACTER(LEN=255)    :: SpecID
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone_FV)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLinearScalarAdvection not ready to be called or already called.")
END IF
LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT DVM...'

Pi=ACOS(-1.)
IniExactFunc_FV = GETINT('IniExactFunc-FV')

DVMDim = GETINT('DVM-Dimension')
IF ((DVMDim.GT.3).OR.(DVMDim.LT.1)) CALL abort(__STAMP__,'DVM error: dimension must be between 1 and 3')
DVMColl = GETLOGICAL('DVM-Collisions')
IF (DVMColl) THEN
  LBWRITE(UNIT_stdOut,'(A)')'DVM: Collisions activated!'
  DVMBGKModel = GETINT('DVM-BGKCollModel')
  DVMMethod = GETINT('DVM-Method')
ELSE
  LBWRITE(UNIT_stdOut,'(A)')'DVM: Collisions deactivated!'
  DVMBGKModel = 0
  DVMMethod = 0
END IF

DVMnSpecies = GETINT('DVM-nSpecies')
ALLOCATE(DVMSpecData(DVMnSpecies))
ALLOCATE(StrVarNames_FV(14*(DVMnSpecies+1)+1))
PP_nVar_FV = 0

ALLOCATE(DVMVeloDisc(DVMnSpecies))

DO iSpec = 1, DVMnSpecies
  LBWRITE (UNIT_stdOut,'(68(". "))')
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  ASSOCIATE(Sp => DVMSpecData(iSpec))
  Sp%omegaVHS     = GETREAL('DVM-Species'//TRIM(hilf)//'-omegaVHS')
  Sp%T_Ref        = GETREAL('DVM-Species'//TRIM(hilf)//'-T_Ref')
  Sp%d_Ref        = GETREAL('DVM-Species'//TRIM(hilf)//'-d_Ref')
  Sp%Internal_DOF = GETINT('DVM-Species'//TRIM(hilf)//'-Internal_DOF')
  Sp%Mass         = GETREAL('DVM-Species'//TRIM(hilf)//'-Mass')
  Sp%Charge       = GETREAL('DVM-Species'//TRIM(hilf)//'-Charge')
  Sp%mu_Ref       = 30.*SQRT(Sp%Mass*BoltzmannConst*Sp%T_Ref/Pi)/(4.*(4.-2.*Sp%omegaVHS)*(6.-2.*Sp%omegaVHS)*Sp%d_Ref**2.)
  Sp%R_S          = BoltzmannConst / Sp%Mass
  Sp%Prandtl      = 2.*(Sp%Internal_DOF + 5.)/(2.*Sp%Internal_DOF + 15.)

  DVMVeloDisc(iSpec)  = GETINT('DVM-Species'//TRIM(hilf)//'-VeloDiscretization')
  Sp%nVelos(1:3)      = 1
  Sp%nVelos(1:DVMDim) = GETINT('DVM-Species'//TRIM(hilf)//'-nVelo')
  Sp%nVar             = Sp%nVelos(1)**DVMDim !non reduced dimensions must have same nVelos for now
  IF (DVMDim.LT.3) Sp%nVar = Sp%nVar*2 !double variables for reduced distributions
  PP_nVar_FV          = PP_nVar_FV + Sp%nVar
  LBWRITE(UNIT_stdOut,*)'DVM species '//TRIM(hilf)//':', Sp%nVelos(1), 'x', Sp%nVelos(2), 'x', Sp%nVelos(3),' velocities!'

  ALLOCATE(Sp%Velos(MAXVAL(Sp%nVelos),3), Sp%Weights(MAXVAL(Sp%nVelos),3))
  Sp%Velos(:,:)=0.
  Sp%Weights(:,:)=1.
  SELECT CASE(DVMVeloDisc(iSpec))
  CASE(1) ! Gauss-Legendre Nodes and Weights
    Sp%VeloMin(:) = GETREALARRAY('DVM-Species'//TRIM(hilf)//'-VeloMin',3)
    Sp%VeloMax(:) = GETREALARRAY('DVM-Species'//TRIM(hilf)//'-VeloMax',3)
    DO iDim=1,DVMDim
      CALL LegendreGaussNodesAndWeights(Sp%nVelos(iDim)-1,Sp%Velos(:,iDim),Sp%Weights(:,iDim))
      Sp%Velos(:,iDim) = Sp%Velos(:,iDim)*ABS(Sp%VeloMax(iDim)-Sp%VeloMin(iDim))/2. &
                      + ABS(Sp%VeloMax(iDim)-Sp%VeloMin(iDim))/2. + Sp%VeloMin(iDim)
      Sp%Weights(:,iDim) = Sp%Weights(:,iDim)*ABS(Sp%VeloMax(iDim)-Sp%VeloMin(iDim))/2.
    END DO
  CASE(2) ! Newton-Cotes Nodes and Weights
    Sp%VeloMin(:) = GETREALARRAY('DVM-Species'//TRIM(hilf)//'-VeloMin',3)
    Sp%VeloMax(:) = GETREALARRAY('DVM-Species'//TRIM(hilf)//'-VeloMax',3)
    Sp%NewtDeg(:) = GETINTARRAY('DVM-Species'//TRIM(hilf)//'-NewtonCotesDegree',3)
    DO iDim=1,DVMDim
      CALL NewtonCotesNodesAndWeights(Sp%nVelos(iDim)-1,Sp%NewtDeg(iDim),Sp%Velos(:,iDim),Sp%Weights(:,iDim))
      Sp%Velos(:,iDim) = Sp%Velos(:,iDim)*ABS(Sp%VeloMax(iDim)-Sp%VeloMin(iDim)) + Sp%VeloMin(iDim)
      Sp%Weights(:,iDim) = Sp%Weights(:,iDim)*ABS(Sp%VeloMax(iDim)-Sp%VeloMin(iDim))
    END DO
  CASE(3) ! Gauss-Hermite Nodes and Weights
    Sp%GHTemp(:) = GETREALARRAY('DVM-Species'//TRIM(hilf)//'-GaussHermiteTemp',3)
    DO iDim=1,DVMDim
      CALL GaussHermiteNodesAndWeights(Sp%nVelos(iDim)-1,Sp%Velos(:,iDim),Sp%Weights(:,iDim))
      DO iGH=1,Sp%nVelos(iDim)
        Sp%Weights(iGH,iDim) = Sp%Weights(iGH,iDim)/(EXP(-Sp%Velos(iGH,iDim)**2.)) &
                            * SQRT(2.*Sp%R_S*Sp%GHTemp(iDim))
      END DO
      Sp%Velos(:,iDim) = Sp%Velos(:,iDim)*SQRT(2.*Sp%R_S*Sp%GHTemp(iDim))
    END DO
  CASE DEFAULT
    CALL abort(__STAMP__,&
              'Undefined velocity space discretization')
  END SELECT ! DVMVeloDisc
  END ASSOCIATE
  ! Set output variable names
  WRITE(SpecID,'(I3.3)') iSpec
  StrVarNames_FV(14*(iSpec-1)+1)  = 'Spec'//TRIM(SpecID)//'_NumberDensity'
  StrVarNames_FV(14*(iSpec-1)+2)  = 'Spec'//TRIM(SpecID)//'_VelocityX'
  StrVarNames_FV(14*(iSpec-1)+3)  = 'Spec'//TRIM(SpecID)//'_VelocityY'
  StrVarNames_FV(14*(iSpec-1)+4)  = 'Spec'//TRIM(SpecID)//'_VelocityZ'
  StrVarNames_FV(14*(iSpec-1)+5)  = 'Spec'//TRIM(SpecID)//'_Temperature'
  StrVarNames_FV(14*(iSpec-1)+6)  = 'Spec'//TRIM(SpecID)//'_PressureXX'
  StrVarNames_FV(14*(iSpec-1)+7)  = 'Spec'//TRIM(SpecID)//'_PressureYY'
  StrVarNames_FV(14*(iSpec-1)+8)  = 'Spec'//TRIM(SpecID)//'_PressureZZ'
  StrVarNames_FV(14*(iSpec-1)+9)  = 'Spec'//TRIM(SpecID)//'_PressureXY'
  StrVarNames_FV(14*(iSpec-1)+10) = 'Spec'//TRIM(SpecID)//'_PressureXZ'
  StrVarNames_FV(14*(iSpec-1)+11) = 'Spec'//TRIM(SpecID)//'_PressureYZ'
  StrVarNames_FV(14*(iSpec-1)+12) = 'Spec'//TRIM(SpecID)//'_HeatfluxX'
  StrVarNames_FV(14*(iSpec-1)+13) = 'Spec'//TRIM(SpecID)//'_HeatfluxY'
  StrVarNames_FV(14*(iSpec-1)+14) = 'Spec'//TRIM(SpecID)//'_HeatfluxZ'
END DO

StrVarNames_FV(14*DVMnSpecies+1)  = 'Total_NumberDensity'
StrVarNames_FV(14*DVMnSpecies+2)  = 'Total_VelocityX'
StrVarNames_FV(14*DVMnSpecies+3)  = 'Total_VelocityY'
StrVarNames_FV(14*DVMnSpecies+4)  = 'Total_VelocityZ'
StrVarNames_FV(14*DVMnSpecies+5)  = 'Total_Temperature'
StrVarNames_FV(14*DVMnSpecies+6)  = 'Total_PressureXX'
StrVarNames_FV(14*DVMnSpecies+7)  = 'Total_PressureYY'
StrVarNames_FV(14*DVMnSpecies+8)  = 'Total_PressureZZ'
StrVarNames_FV(14*DVMnSpecies+9)  = 'Total_PressureXY'
StrVarNames_FV(14*DVMnSpecies+10) = 'Total_PressureXZ'
StrVarNames_FV(14*DVMnSpecies+11) = 'Total_PressureYZ'
StrVarNames_FV(14*DVMnSpecies+12) = 'Total_HeatfluxX'
StrVarNames_FV(14*DVMnSpecies+13) = 'Total_HeatfluxY'
StrVarNames_FV(14*DVMnSpecies+14) = 'Total_HeatfluxZ'
StrVarNames_FV(14*DVMnSpecies+15) = 'RelaxationFactor'

! Read Boundary information / RefStates / perform sanity check
IniRefState_FV = GETINT('IniRefState-FV')
nRefState_FV=CountOption('RefState-FV')/DVMnSpecies
IF(IniRefState_FV.GT.nRefState_FV)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState_FV):',IniRefState_FV,REAL(nRefState_FV))
END IF

IF(nRefState_FV .GT. 0)THEN
  ALLOCATE(RefState_FV(14,DVMnSpecies,nRefState_FV))
  DO i=1,nRefState_FV
    DO iSpec=1,DVMnSpecies
      RefState_FV(1:14,iSpec,i)  = GETREALARRAY('RefState-FV',14)
    END DO
  END DO
END IF

DVMForce = GETREALARRAY('DVM-Accel',3)

ALLOCATE(DVMMomentSave(15,DVMnSpecies+1,nElems))
DVMMomentSave = 0.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

SELECT CASE (DVMMethod)
CASE(1)
  LBWRITE(UNIT_stdOut,*)'Method: Exponential Differencing DVM'
CASE(2)
  LBWRITE(UNIT_stdOut,*)'Method: DUGKS'
END SELECT

WriteDVMSurfaceValues = GETLOGICAL('DVM-WriteMacroSurfaceValues')
IF (WriteDVMSurfaceValues) CALL InitDVMBoundaryAnalyze()

EquationInitIsDone_FV=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT DVM DONE!'
LBWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquation


SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars,  ONLY: PI, BoltzmannConst
USE MOD_DistFunc,      ONLY: MaxwellDistribution, GradDistribution
USE MOD_Equation_Vars_FV, ONLY: DVMSpecData, RefState_FV, DVMnSpecies
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of
                                                          !< timestep if full boundary order is used (only with RK3)
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
INTEGER                         :: iSpec,vFirstID,vLastID
!==================================================================================================================================
Resu=0.

vFirstID=1
vLastID=0
DO iSpec=1,DVMnSpecies
  vLastID = vLastID + DVMSpecData(iSpec)%nVar

  SELECT CASE (ExactFunction)
  CASE(0)
    Resu(vFirstID:vLastID)=0.

  CASE(1) !Grad 13 uniform init
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)

  CASE(2) ! couette flow, L=1m between walls, RefState_FV: 1 -> initial state, 2 -> y=-0.5 boundary, 3 -> y=0.5 boundary
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
    ! steady state
    IF (tIn.GT.0.) THEN
      MacroVal(:) = RefState_FV(:,iSpec,1)
      WallVelo1 = RefState_FV(2,iSpec,2)
      WallTemp1 = RefState_FV(5,iSpec,2)
      WallVelo2 = RefState_FV(2,iSpec,3)
      WallTemp2 = RefState_FV(5,iSpec,3)
      MacroVal(2) = WallVelo2 + (WallVelo1-WallVelo2)*(0.5-x(2))
      MacroVal(5) = WallTemp2 + (0.5-x(2))*(WallTemp1-WallTemp2+((WallVelo1-WallVelo2)**2)*(0.5+x(2))*4/15/DVMSpecData(iSpec)%R_S/2)
                                                                                                      ! cf Eucken's relation
      CALL MaxwellDistribution(MacroVal,Resu(vFirstID:vLastID),iSpec)
    END IF

  CASE(3) !sod shock
    IF (tIn.EQ.0.) THEN ! initial state
      IF (x(1).LT.0.) THEN
        CALL MaxwellDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
      ELSE
        CALL MaxwellDistribution(RefState_FV(:,iSpec,2),Resu(vFirstID:vLastID),iSpec)
      END IF
    ELSE ! analytical solution
      MacroVal = 0.
      SodMacro_L=RefState_FV(1:5,iSpec,1)
      SodMacro_R=RefState_FV(1:5,iSpec,2)
      SodMacro_LL = 0.         !     | L | LL | M | RR | R |
      SodMacro_M = 0.
      SodMacro_RR = 0.
      gamma = 5./3. !monatomic gas
      Ggamma = (gamma-1)/(gamma+1)
      beta = (gamma-1)/gamma/2.
      pL = BoltzmannConst*SodMacro_L(1)*SodMacro_L(5)
      pR = BoltzmannConst*SodMacro_R(1)*SodMacro_R(5)
      pM=(pL+pR)/2.
      CALL SecantSod(pM, pL, pR, SodMacro_L(1)*DVMSpecData(iSpec)%Mass, SodMacro_R(1)*DVMSpecData(iSpec)%Mass, gamma, Ggamma, beta, 1e-15, 100)
      cL = sqrt(gamma*DVMSpecData(iSpec)%R_S*SodMacro_L(5))
      cR = sqrt(gamma*DVMSpecData(iSpec)%R_S*SodMacro_R(5))
      cM = cL * (pM/pL)**beta
      SodMacro_M(1) = SodMacro_L(1)*(pM/pL)**(1./gamma)
      SodMacro_M(2) = (cL-cM)*2/(gamma-1)
      SodMacro_M(5) = pM/BoltzmannConst/SodMacro_M(1)
      SodMacro_RR(1) = SodMacro_R(1)*(pM+Ggamma*pR)/(pR+Ggamma*pM)
      SodMacro_RR(2) = SodMacro_M(2)
      SodMacro_RR(5) = pM/BoltzmannConst/SodMacro_RR(1)
      vs = cR*sqrt((beta/Ggamma)*(pM/pR + Ggamma))
      IF (x(1).LT.(-tIn*cL)) THEN
        MacroVal(1:5) = SodMacro_L
      ELSE IF (x(1).LT.(tIn*(SodMacro_M(2)-cM))) THEN
        SodMacro_LL(2) = 2./(gamma+1.) * (cL + x(1)/tIn)
        SodMacro_LL(1) = SodMacro_L(1) * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2./(gamma-1.))
        SodMacro_LL(5) = pL * (1.-(gamma-1.)*SodMacro_LL(2)/cL/2.)**(2*gamma/(gamma-1))/BoltzmannConst/SodMacro_LL(1)
        MacroVal(1:5) = SodMacro_LL
      ELSE IF (x(1).LT.(tIn*SodMacro_M(2))) THEN
        MacroVal(1:5) = SodMacro_M
      ELSE IF (x(1).LT.(tIn*vs)) THEN
        MacroVal(1:5) = SodMacro_RR
      ELSE
        MacroVal(1:5) = SodMacro_R
      END IF
      CALL MaxwellDistribution(MacroVal,Resu(vFirstID:vLastID),iSpec)
    END IF

  CASE(4) !heat flux relaxation test case
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
    IF (tIn.GT.0.) THEN ! relaxation
      MacroVal(:) = RefState_FV(:,iSpec,1)
      mu = DVMSpecData(iSpec)%mu_Ref*(MacroVal(5)/DVMSpecData(iSpec)%T_Ref)**(DVMSpecData(iSpec)%omegaVHS+0.5)
      tau = mu/(BoltzmannConst*MacroVal(1)*MacroVal(5))
      MacroVal(12:14) = MacroVal(12:14)*EXP(-tIn*DVMSpecData(iSpec)%Prandtl/tau) !Heat flux relaxes with rate Pr/tau
      CALL GradDistribution(MacroVal(:),Resu(vFirstID:vLastID),iSpec)
    END IF

  CASE(5) !Taylor-Green vortex
    MacroVal(:) = RefState_FV(:,iSpec,1)
    MacroVal(2) = RefState_FV(2,iSpec,1)*SIN(x(1))*COS(x(2))*COS(x(3))
    MacroVal(3) = -RefState_FV(2,iSpec,1)*COS(x(1))*SIN(x(2))*COS(x(3))
    MacroVal(4) = 0.
    MacroVal(5) = MacroVal(5)+(RefState_FV(2,iSpec,1)**2)/(16.*DVMSpecData(iSpec)%R_S)*(COS(2.*x(1))+COS(2.*x(2)))*(COS(2.*x(3))+2.)
    CALL MaxwellDistribution(MacroVal,Resu(vFirstID:vLastID),iSpec)

  CASE(6) !Sum of 2 distributions
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
    CALL GradDistribution(RefState_FV(:,iSpec,2),SecondDist(vFirstID:vLastID),iSpec)
    Resu(vFirstID:vLastID) = Resu(vFirstID:vLastID) + SecondDist(vFirstID:vLastID)

  CASE(7) ! Poiseuille flow with force corresponding to 0.01 Pa/m
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
    ! steady state in continuum limit
    IF (tIn.GT.0.) THEN
      MacroVal(:) = RefState_FV(:,iSpec,1)
      mu = DVMSpecData(iSpec)%mu_Ref*(MacroVal(5)/DVMSpecData(iSpec)%T_Ref)**(DVMSpecData(iSpec)%omegaVHS+0.5)
      MacroVal(2) = 0.01*(1.*x(2)-x(2)*x(2))/mu/2.
      CALL MaxwellDistribution(MacroVal,Resu(vFirstID:vLastID),iSpec)
    END IF

  CASE(11) !Grad 13 uniform init with perturbation (two stream instability)
    CALL GradDistribution(RefState_FV(:,iSpec,1),Resu(vFirstID:vLastID),iSpec)
    IF (iSpec.GE.2.AND.tIn.EQ.0..AND.x(1).GT.0.49) THEN
      print*, 'perturbation at x = ',x(1)
      Resu(vFirstID:vLastID)=1.1*Resu(vFirstID:vLastID)
    END IF
    IF (iSpec.GE.2.AND.tIn.EQ.0..AND.x(1).LT.-0.49) THEN
      print*, 'perturbation at x = ',x(1)
      Resu(vFirstID:vLastID)=0.9*Resu(vFirstID:vLastID)
    END IF

  CASE DEFAULT
    CALL abort(__STAMP__,&
              'Specified exact function not implemented!')
  END SELECT ! ExactFunction

  vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
END DO !iSpec

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
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!==================================================================================================================================


END SUBROUTINE CalcSource


!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars_FV,ONLY:EquationInitIsDone_FV, RefState_FV, WriteDVMSurfaceValues, DVMSpecData, DVMnSpecies
USE MOD_Equation_Vars_FV,ONLY: DVMMomentSave, DVMVeloDisc, StrVarNames_FV
USE MOD_DVM_Boundary_Analyze,ONLY: FinalizeDVMBoundaryAnalyze
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iSpec
!==================================================================================================================================
EquationInitIsDone_FV = .FALSE.
IF (WriteDVMSurfaceValues) CALL FinalizeDVMBoundaryAnalyze()
SDEALLOCATE(DVMMomentSave)
SDEALLOCATE(DVMVeloDisc)
SDEALLOCATE(StrVarNames_FV)
SDEALLOCATE(RefState_FV)
DO iSpec=1,DVMnSpecies
  SDEALLOCATE(DVMSpecData(iSpec)%Velos)
  SDEALLOCATE(DVMSpecData(iSpec)%Weights)
END DO
SDEALLOCATE(DVMSpecData)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation_FV
