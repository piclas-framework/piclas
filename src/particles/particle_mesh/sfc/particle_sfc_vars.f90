#include "boltzplatz.h"

MODULE MOD_Particle_SFC_Vars
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
#ifdef donotcompilethis
LOGICAL                       ::ParticleSFCInitIsDone=.FALSE.                       ! flag if init is done
INTEGER                       ::whichBoundBox                                       ! select bounding box for SFC
TYPE tBox
  REAL(KIND=8) :: mini(3)
  INTEGER      :: nbits
  REAL(KIND=8) :: spacing(3)
END TYPE tBox
#endif /*donotcompilethis*/

END MODULE MOD_Particle_SFC_Vars
