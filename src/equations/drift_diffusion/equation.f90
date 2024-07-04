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
CALL prms%SetSection("Equation-FV")
CALL prms%CreateIntOption(      'IniExactFunc-FV'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection', '-1')
CALL prms%CreateIntOption(      'IniRefState-FV',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState-FV',     "State(s) in primitive variables (density).",&
                                                multiple=.TRUE., no=1 )
END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters (advection velocity, diffusion coeff, exact function)  from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nSides
USE MOD_ReadInTools,        ONLY:GETREALARRAY,GETINTARRAY,GETREAL,GETINT, CountOption
USE MOD_Equation_Vars_FV
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================
IF(EquationInitIsDone_FV)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLinearScalarAdvection not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Drift-Diffusion...'

IniExactFunc_FV = GETINT('IniExactFunc-FV')


! Read Boundary information / RefStates / perform sanity check
IniRefState_FV = GETINT('IniRefState-FV')
nRefState_FV=CountOption('RefState-FV')
IF(IniRefState_FV.GT.nRefState_FV)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState_FV):',IniRefState_FV,REAL(nRefState_FV))
END IF

IF(nRefState_FV .GT. 0)THEN
  ALLOCATE(RefState_FV(1,nRefState_FV))
  DO i=1,nRefState_FV
    RefState_FV(1:1,i)  = GETREALARRAY('RefState-FV',1)
  END DO
END IF

ALLOCATE(EFluid_GradSide(nSides))
EFluid_GradSide=0.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

EquationInitIsDone_FV=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Drift-Diffusion DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquation


SUBROUTINE ExactFunc(ExactFunction,tIn,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars_FV, ONLY: RefState_FV
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

!==================================================================================================================================

Resu   =0.

SELECT CASE (ExactFunction)
CASE(0)
  Resu=0.

CASE(1)
  Resu=RefState_FV(:,1)

CASE(2) !shock
  IF (x(1).LT.0.) THEN
    Resu = RefState_FV(:,1)
  ELSE
    Resu = RefState_FV(:,2)
  END IF

CASE(3) !1D streamer (Markosyan et al. - 2013)
  Resu(1) = RefState_FV(1,1)*EXP(-((x(1)-0.8e-3)/2.9e-5)**2)

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
USE MOD_Globals           ,ONLY: abort, vecnorm
USE MOD_PreProc
USE MOD_FV_Vars           ,ONLY: U_FV
USE MOD_Equation_Vars     ,ONLY: E
USE MOD_Transport_Data    ,ONLY: CalcDriftDiffusionCoeff
USE MOD_DSMC_Vars         ,ONLY: BGGas
USE MOD_Particle_Vars     ,ONLY: nSpecies, Species
USE MOD_Interpolation_Vars,ONLY: wGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID, iSpec, i, j, k
REAL                            :: ionRate, mu, D, E_avg(1:3)
!===================================================================================================================================
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    EXIT
  END IF
END DO

DO ElemID = 1, PP_nElems
  E_avg = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    E_avg(:) = E_avg(:) + wGP(i)*wGP(j)*wGP(k)*E(1:3,i,j,k,ElemID)/((PP_N+1.)**3)
  END DO; END DO; END DO

  CALL CalcDriftDiffusionCoeff(VECNORM(E_avg),BGGas%NumberDensity(iSpec),mu,D,ionRate)
  Ut(1,:,:,:,ElemID) = Ut(1,:,:,:,ElemID) + ionRate*mu*VECNORM(E_avg)*U_FV(1,:,:,:,ElemID)
END DO

END SUBROUTINE CalcSource

!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars_FV,ONLY:EquationInitIsDone_FV, EFluid_GradSide, RefState_FV
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone_FV = .FALSE.
SDEALLOCATE(EFluid_GradSide)
SDEALLOCATE(RefState_FV)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation_FV
