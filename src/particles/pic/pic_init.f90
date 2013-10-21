#include "boltzplatz.h"

MODULE MOD_PICInit
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
INTERFACE InitPIC
  MODULE PROCEDURE InitPIC
END INTERFACE
PUBLIC::InitPIC
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPIC()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file 
!===================================================================================================================================
! MODULES
USE MOD_Globals, ONLY: MPIRoot,UNIT_STDOUT
USE MOD_ReadInTools
USE MOD_PICInterpolation_Vars, ONLY:externalField
USE MOD_PIC_Vars!, ONLY: 
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf 
  CHARACTER(200)        :: tmpString
#ifdef MPI
#endif
!===================================================================================================================================
IF(PICInitIsDone)THEN
   SWRITE(*,*) "InitPIC already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PIC ...'

! So far, nothing to do here...
IF (externalField(6).NE.0) PIC%GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
