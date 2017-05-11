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
INTEGER                                 :: nCurrentParts                 ! current number of particles
LOGICAL                                 :: MeasureTrackTime              ! switch, if tracking time is measure
LOGICAL                                 :: FastPeriodic                  ! moves the particle along whole periodic vector, 
                                                                         ! neglecting possible reflexions
LOGICAL                                 :: CountNbOfLostParts            ! logical, to count the lost particles
LOGICAL                                 :: CartesianPeriodic             ! old periodic for refmapping and ALL bcs periocic
INTEGER                                 :: nLostParts                    ! Counter for lost particle
REAL,ALLOCATABLE                        :: Distance(:)                   ! list of distance between particle and element-origin
                                                                         ! to all elements in the same background element
INTEGER,ALLOCATABLE                     :: ListDistance(:)               ! the corresponding element id
#ifdef CODE_ANALYZE
INTEGER                                 :: PartOut
INTEGER                                 :: MPIRankOut
#endif /*CODE_ANALYZE*/
!===================================================================================================================================


END MODULE MOD_Particle_Tracking_Vars
