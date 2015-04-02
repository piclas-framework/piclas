#include "boltzplatz.h"

MODULE MOD_Particle_SFC
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitParticleSFC
  MODULE PROCEDURE InitParticleSFC
END INTERFACE

INTERFACE FinalizeParticleSFC
  MODULE PROCEDURE FinalizeParticleSFC
END INTERFACE


PUBLIC::InitParticleSFC,FinalizeParticleSFC
!===================================================================================================================================
!
CONTAINS

SUBROUTINE InitParticleSFC()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_SFC_Vars
!USE MOD_Particle_Surfaces_Vars, ONLY:neighborElemID,neighborLocSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ALLOCSTAT
INTEGER           :: iElem, ilocSide,SideID,flip,iSide
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SFC ...'
IF(ParticleSFCInitIsDone) &
    CALL abort(__STAMP__,&
    ' process local space-filling curve for particle localization already allocated!.')

ParticleSFCInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SFC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleSFC


SUBROUTINE FinalizeParticleSFC()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_SFC_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ParticleSFCInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleSFC

END MODULE MOD_Particle_SFC
