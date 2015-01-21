#include "boltzplatz.h"

MODULE MOD_Particle_MPI_Halo
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

#ifdef MPI
INTERFACE IdentifyHaloMPINeighborhood
  MODULE PROCEDURE IdentifyHaloMPINeighborhood
END INTERFACE


PUBLIC :: IdentifyHaloMPINeighborhood

!===================================================================================================================================

CONTAINS



SUBROUTINE IdentifyHaloMPINeighborhood(iProc,SideIndex)
!===================================================================================================================================
! Searches for sides in the neighborhood of MPI-neighbors
! mark all Sides for communication which are in range of halo_eps of own MPI_Sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)   :: iProc  ! MPI proc with which the local proc has to exchange boundary information
INTEGER, INTENT(INOUT):: SideIndex(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iBGMCellX,iBGMCellY,iBGMCellZ 
!=================================================================================================================================

! 1) Exchange Sides:
!    Each proc receives all sides that lie in FIBGM cells within an  eps-distance of FIBGM cells containing 
!    any of my own MPI sides.
! 2) Go through each FIBGM cell and if there are MPI-Neighbors, search the neighbor  surrounding FIBGM
!     cells for my own MPI sides.

! Step1: find Sides of myProc that are within halo_eps distance to iProc

  SideIndex(:)=0

  DO iBGMCellX=GEO%FIBGMimin,GEO%FIBGMimax
    DO iBGMCellY=GEO%FIBGMkmin,GEO%FIBGMkmax
      DO iBGMCellZ=GEO%FIBGMlmin,GEO%FIBGMlmax




END SUBROUTINE IdentifyHaloMPINeighborhood

#endif /*MPI*/
END MODULE Particle_MPI_Halo
