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
#if USE_QDS_DG
#include "piclas.h"

MODULE MOD_QDS
!===================================================================================================================================
!> Contains the routines to
!> - initialze the QDS method for DG + equation
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitQDS
  MODULE PROCEDURE InitQDS
END INTERFACE
INTERFACE FinalizeQDS
  MODULE PROCEDURE FinalizeQDS
END INTERFACE


PUBLIC::InitQDS
PUBLIC::FinalizeQDS
PUBLIC::DefineParametersQDS
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for QDS
!==================================================================================================================================
SUBROUTINE DefineParametersQDS()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("QDS")

CALL prms%CreateLOGICALOption(  'DoQDS'              , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntOption(      'QDSIniExactFunc'    , 'TODO-DEFINE-PARAMETER', '0')

!QDS_Species = GETINT('Particles-QDSSpecies','0')

END SUBROUTINE DefineParametersQDS


SUBROUTINE InitQDS
!===================================================================================================================================
!> Allocate all QDS variables, determine
!===================================================================================================================================
! MODULES
USE MOD_Globals,         ONLY:UNIT_stdOut,mpiroot
USE MOD_QDS_DG,          ONLY:QDS_InitDG
USE MOD_QDS_Equation,    ONLY:QDS_InitEquation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL    :: tempNorm
!INTEGER :: iWeight,i,j,k
!REAL    :: Velo(3), Temp, Dens, Mass
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT QDS ...'

CALL QDS_InitEquation()
CALL QDS_InitDG()


SWRITE(UNIT_stdOut,'(A)')' INIT QDS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitQDS


SUBROUTINE FinalizeQDS
!===================================================================================================================================
!> Allocate all QDS variables, determine
!===================================================================================================================================
! MODULES
USE MOD_QDS_DG,          ONLY:QDS_FinalizeDG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL QDS_FinalizeDG()
!CALL QDS_FinalizeEquation()

END SUBROUTINE FinalizeQDS


END MODULE MOD_QDS
#endif /*USE_QDS_DG*/
