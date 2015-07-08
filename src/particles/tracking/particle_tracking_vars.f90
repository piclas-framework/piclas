#include "boltzplatz.h"

MODULE MOD_Particle_Tracking_Vars
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
LOGICAL                                 :: DoRefMapping                  ! tracking by mapping particle into reference element
REAL                                    :: tTracking                     ! Tracking time
REAL                                    :: tLocalization                 ! localization time
INTEGER                                 :: nTracks                       ! number of tracked particles
LOGICAL                                 :: MeassureTrackTime             ! switch, if tracking time is meassured
!===================================================================================================================================


END MODULE MOD_Particle_Tracking_Vars
