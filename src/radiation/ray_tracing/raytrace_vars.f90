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

MODULE MOD_RayTracing_Vars
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
TYPE tRayTrace
  REAL    :: PulseDuration      !>
  REAL    :: tShift             !>
  REAL    :: tActive            !>
  REAL    :: Period             !>
  INTEGER :: NbrOfPulses        !>
  REAL    :: WaistRadius        !>
  REAL    :: WaveLength         !>
  REAL    :: RepetitionRate     !>
  REAL    :: Power              !>
  REAL    :: Area               !>
  REAL    :: Energy             !>
  REAL    :: IntensityAmplitude !>
END TYPE

TYPE (tRayTrace)     :: Ray                            !>

TYPE tRadTrans
  INTEGER            :: NumPhotonsPerCell              !>
  REAL               :: GlobalRadiationPower           !>
  REAL               :: ScaledGlobalRadiationPower     !>
  INTEGER            :: GlobalPhotonNum                !>
END TYPE

TYPE (tRadTrans)     :: RadTrans                       !>

TYPE tPhotonProps
  REAL               :: PhotonPos(3)                   !>
  REAL               :: PhotonLastPos(3)               !>
  REAL               :: PhotonDirection(3)             !>
  REAL               :: PhotonEnergy                   !>
  INTEGER            :: ElemID                         !>
  INTEGER            :: WaveLength                     !>
END TYPE

TYPE (tPhotonProps)  :: PhotonProps                    !>

LOGICAL              :: AdaptiveRays                   !>
INTEGER              :: NumRays                        !>
INTEGER              :: RayPosModel                    !>
INTEGER              :: RayPartBound                   !> Particle boundary ID where rays are emitted from
REAL                 :: RayDirection(3)                !>

REAL, ALLOCATABLE    :: RayElemPassedEnergy(:)         !>
REAL, ALLOCATABLE    :: RaySampWall(:,:)               !>
#if USE_MPI
INTEGER              :: RaySampWall_Shared_Win         !>
REAL,POINTER         :: RaySampWall_Shared(:,:)        !>
INTEGER              :: RayElemPassedEnergy_Shared_Win !>
REAL,POINTER         :: RayElemPassedEnergy_Shared(:)  !>
#endif
!===================================================================================================================================
END MODULE MOD_RayTracing_Vars
