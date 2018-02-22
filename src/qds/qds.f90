#if USE_QDS_DG
#include "boltzplatz.h"

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
