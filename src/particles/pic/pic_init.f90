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
INTERFACE InitPIC
  MODULE PROCEDURE InitPIC
END INTERFACE
PUBLIC::InitPIC
!===================================================================================================================================
PUBLIC::DefineParametersPIC
CONTAINS

!==================================================================================================================================
!> Define parameters for PIC
!==================================================================================================================================
SUBROUTINE DefineParametersPIC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools      ,ONLY: prms
USE MOD_PICDepo_Method   ,ONLY: DefineParametersDepositionMethod
USE MOD_PICInterpolation ,ONLY: DefineParametersPICInterpolation
USE MOD_PICDepo          ,ONLY: DefineParametersPICDeposition
IMPLICIT NONE
!==================================================================================================================================
CALL DefineParametersPICInterpolation() ! Get PIC Interpolation parameters
CALL DefineParametersDepositionMethod() ! Get PIC-DoDeposition and PIC-Deposition-Type
CALL DefineParametersPICDeposition()    ! Get remaining PIC Deposition parameters
END SUBROUTINE DefineParametersPIC


SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation       ,ONLY: InitializeParticleInterpolation
USE MOD_PICDepo                ,ONLY: InitializeDeposition
USE MOD_PIC_Vars ,              ONLY: PICInitIsDone
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
   SWRITE(*,*) "InitPIC already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PIC ...'

CALL InitializeParticleInterpolation()
CALL InitializeDeposition()

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
