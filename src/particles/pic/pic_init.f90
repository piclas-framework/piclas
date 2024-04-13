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

MODULE MOD_PICInit
!===================================================================================================================================
! Includes PIC Init
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
PUBLIC::InitPIC
!===================================================================================================================================
PUBLIC::DefineParametersPIC
CONTAINS

!==================================================================================================================================
!> Define parameters for PIC
!==================================================================================================================================
SUBROUTINE DefineParametersPIC()
! MODULES
USE MOD_PICInterpolation          ,ONLY: DefineParametersPICInterpolation
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
USE MOD_PICDepo_Method            ,ONLY: DefineParametersDepositionMethod
USE MOD_PICDepo                   ,ONLY: DefineParametersPICDeposition
USE MOD_InitializeBackgroundField ,ONLY: DefineParametersBGField
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
IMPLICIT NONE
!==================================================================================================================================
CALL DefineParametersPICInterpolation() ! Get PIC Interpolation parameters
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
CALL DefineParametersDepositionMethod() ! Get PIC-DoDeposition and PIC-Deposition-Type
CALL DefineParametersPICDeposition()    ! Get more PIC Deposition parameters
CALL DefineParametersBGField()          ! Get PIC Background field parameters
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
END SUBROUTINE DefineParametersPIC


SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation          ,ONLY: InitializeParticleInterpolation
USE MOD_PIC_Vars                  ,ONLY: PICInitIsDone
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
USE MOD_PICInterpolation_Vars     ,ONLY: useBGField
USE MOD_PICDepo                   ,ONLY: InitializeDeposition
USE MOD_InitializeBackgroundField ,ONLY: InitializeBackgroundField
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(PICInitIsDone)THEN
   LBWRITE(*,*) "InitPIC already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT PIC ...'

CALL InitializeParticleInterpolation()
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
CALL InitializeDeposition()
IF(useBGField) CALL InitializeBackgroundField()
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

PICInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
