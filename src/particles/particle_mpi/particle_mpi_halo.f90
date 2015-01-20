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
INTERFACE InitFIBGM
  MODULE PROCEDURE InitFIBGM
END INTERFACE

INTERFACE IdentifyHaloMPINeighborhood
  MODULE PROCEDURE IdentifyHaloMPINeighborhood
END INTERFACE


PUBLIC :: IdentifyHaloMPINeighborhood, InitFIBGM

!===================================================================================================================================

CONTAINS

SUBROUTINE InitFIBGM()
!===================================================================================================================================
! Build Fast-Init-Background-Mesh.
! The BGM is a cartesian mesh for easier locating of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Surfaces_Vars,             ONLY:BezierControlPoints3D
USE MOD_Mesh_Vars,                          ONLY:nSides,NGeo!,ElemToSide,SideToElem,NGeo
#ifdef MPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: xmin, xmax, ymin, ymax, zmin, zmax                          !
INTEGER                          :: iSide
#ifdef MPI
INTEGER                          :: ii,jj,kk
#endif /*MPI*/
!=================================================================================================================================

! zeros
#ifdef MPI
ii=0
jj=0
kk=0
#endif /*MPI*/

!#ifdef MPI
!   !--- If this MPI process does not contain particles, step out
!   IF (PMPIVAR%GROUP.EQ.MPI_GROUP_EMPTY) RETURN
!#endif
!--- calc min and max coordinates for mesh
xmin = HUGE(1.0)
xmax =-HUGE(1.0)
ymin = HUGE(1.0)
ymax =-HUGE(1.0)
zmin = HUGE(1.0)
zmax =-HUGE(1.0)

! serch for min,max of BezierControlPoints, e.g. the convec hull of the domain
DO iSide=1,nSides
  xmin=MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,iSide)))
  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,iSide)))
  ymin=MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,iSide)))
  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,iSide)))
  zmin=MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,iSide)))
  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,iSide)))
END DO ! iSide

CALL InitPeriodic()
CALL InitializeInterpolation()
CALL InitializeDeposition()
CALL InitPIC()



END SUBROUTINE InitFIBGM


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
