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

PUBLIC::InitEquation,ExactFunc,FinalizeEquation,CalcSource
!===================================================================================================================================

PUBLIC::DefineParametersEquation
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
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density).",&
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
USE MOD_Equation_Vars
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================
IF(EquationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLinearScalarAdvection not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Drift-Diffusion...'

IniExactFunc = GETINT('IniExactFunc')


! Read Boundary information / RefStates / perform sanity check
IniRefState = GETINT('IniRefState')
nRefState=CountOption('RefState')
IF(IniRefState.GT.nRefState)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))
END IF

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefState(1,nRefState))
  DO i=1,nRefState
    RefState(1:1,i)  = GETREALARRAY('RefState',1)
  END DO
END IF

ALLOCATE(EFluid_GradSide(nSides))

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.


EquationInitIsDone=.TRUE.
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
USE MOD_Equation_Vars, ONLY: RefState
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
  Resu=RefState(:,1)

CASE(2) !shock
  IF (x(1).LT.0.) THEN
    Resu = 2.!e21
  ELSE
    Resu = 1.!e21
  END IF

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
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar_FV,0:0,0:0,0:0,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================


END SUBROUTINE CalcSource

!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone, EFluid_GradSide, RefState
IMPLICIT NONE
!==================================================================================================================================
EquationInitIsDone = .FALSE.
SDEALLOCATE(EFluid_GradSide)
SDEALLOCATE(RefState)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation
