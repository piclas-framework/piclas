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
LOGICAL                       ::ParticleSFCInitIsDone=.FALSE.

END MODULE MOD_Particle_SFC_Vars
