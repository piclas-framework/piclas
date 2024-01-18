!==================================================================================================================================
! Copyright (c) 2023 - 2023 Marcel Pfeiffer, Stephen Copplestone
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Photon_TrackingVars
!===================================================================================================================================
! Contains the tadiation transport variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL RAY TRACING VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tPhotonProps
  REAL               :: PhotonPos(3)        !>
  REAL               :: PhotonLastPos(3)    !>
  REAL               :: PhotonDirection(3)  !>
  REAL               :: PhotonDirectionBeforeReflection(3)  !>
  REAL               :: PhotonEnergy        !>
  INTEGER            :: ElemID              !>
  INTEGER            :: WaveLength          !>
END TYPE

TYPE (tPhotonProps)  :: PhotonProps         !>

REAL,ALLOCPOINT      :: PhotonSampWallHDF5(:,:,:,:)
#if USE_MPI
INTEGER              :: PhotonSampWallHDF5_Shared_Win
REAL,POINTER         :: PhotonSampWallHDF5_Shared(:,:,:,:)
#endif

REAL,ALLOCPOINT      :: PhotonSampWall(:,:,:,:)
REAL,ALLOCATABLE     :: PhotonSampWall_loc(:,:,:)

INTEGER              :: PhotonModeBPO       !> 0: Output nothing to PartStateBoundary.h5
                                            !> 1: Output the initial position of the rays and their direction vector
                                            !> 2: Output initial position and all calculated intersection points calculated in radtrans_tracking.f90
LOGICAL              :: UsePhotonTriaTracking !> True/False: Use TriaTracking methods for photon tracking or Bilinear methods (default is True)

#if USE_MPI
INTEGER              :: PhotonSampWall_Shared_Win
REAL,POINTER         :: PhotonSampWall_Shared(:,:,:,:)
LOGICAL              :: PhotonSampWall_Shared_Win_allocated
REAL,ALLOCATABLE     :: PhotonSampWallProc(:,:,:,:)
#endif /*USE_MPI*/

REAL,ALLOCPOINT,DIMENSION(:,:,:,:) :: PhotonSurfSideSamplingMidPoints     !> Mid point of supersampled surface side
REAL,ALLOCPOINT,DIMENSION(:,:,:)   :: PhotonSurfSideArea                  !> Area of supersampled surface side

#if USE_MPI
REAL,POINTER,DIMENSION(:,:,:,:)    :: PhotonSurfSideSamplingMidPoints_Shared     !> Physical coordinate of the center of supersampled surface side
INTEGER                            :: PhotonSurfSideSamplingMidPoints_Shared_Win
REAL,POINTER,DIMENSION(:,:,:)      :: PhotonSurfSideArea_Shared                  !> Area of supersampled surface side
INTEGER                            :: PhotonSurfSideArea_Shared_Win
#endif /*USE_MPI*/


CHARACTER(LEN=255) :: RadiationSurfState,RadiationVolState !> Output file names for surface and volume sampling
!===================================================================================================================================
END MODULE MOD_Photon_TrackingVars
