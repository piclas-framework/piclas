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
  !REAL               :: nSurfSampleFac      !> Scaling factor: nSurfSampleFac= 1.0/(nSurfSample**2)
  REAL               :: PhotonStartPos(3)   !> super sampled ray path
END TYPE

TYPE (tPhotonProps)  :: PhotonProps         !>

REAL, ALLOCATABLE    :: PhotonSampWall(:,:,:,:)

INTEGER              :: PhotonModeBPO       !> 0: Output nothing to PartStateBoundary.h5
                                            !> 1: Output the initial position of the rays and their direction vector
                                            !> 2: Output initial position and all calculated intersection points calculated in radtrans_tracking.f90

#if USE_MPI
INTEGER              :: PhotonSampWall_Shared_Win
REAL,POINTER         :: PhotonSampWall_Shared(:,:,:,:)
#endif /*USE_MPI*/
!===================================================================================================================================
END MODULE MOD_Photon_TrackingVars
