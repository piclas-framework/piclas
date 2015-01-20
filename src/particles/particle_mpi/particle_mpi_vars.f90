#include "boltzplatz.h"

MODULE MOD_Particle_MPI_Vars
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
INTEGER, ALLOCATABLE                     :: casematrix(:,:)                   ! matrix to compute periodic cases
INTEGER                                  :: NbrOfCases                        ! Number of periodic cases

!===================================================================================================================================

END MODULE Particle_MPI_Vars
