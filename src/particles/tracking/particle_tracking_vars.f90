!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
INTEGER               :: TrackingMethod         ! Method that is used for tracking of particles
LOGICAL               :: DoRefMapping           ! tracking by mapping particle into reference element
LOGICAL               :: TriaTracking           ! tracking by using surface representation by triangles
LOGICAL               :: WriteTriaDebugMesh     !
REAL                  :: tTracking              ! Tracking time
REAL                  :: tLocalization          ! localization time
INTEGER               :: nTracks                ! number of tracked particles
INTEGER               :: nCurrentParts          ! current number of particles
LOGICAL               :: MeasureTrackTime       ! switch, if tracking time is measure
LOGICAL               :: FastPeriodic           ! moves the particle along whole periodic vector,
                                                ! neglecting possible reflexions
LOGICAL               :: DisplayLostParticles   ! Display position, velocity, species and hots element of particles lost during 
!                                               ! particle tracking (TrackingMethod = triatracking, tracing). Default=F
LOGICAL               :: CountNbrOfLostParts    ! logical, to count the lost particles
REAL, ALLOCATABLE     :: PartStateLost(:,:)     ! (1:14,1:NParts) 1st index: LastPartPos-X,LastPartPos-Y,LastPartPos-Z,vx,vy,vz,
!                                               !                            SpecID,MPF,time,ElemID,iPart,x,y,z
!                                               !                 2nd index: 1 to number of lost particles
INTEGER               :: PartStateLostVecLength ! Number of lost particles for each process
LOGICAL               :: CartesianPeriodic      ! old periodic for refmapping and ALL bcs periocic
INTEGER               :: NbrOfLostParticles     ! Counter for lost particle per process
INTEGER               :: NbrOfLostParticlesTotal! Counter for lost particles across all processes
REAL,ALLOCATABLE      :: Distance(:)            ! list of distance between particle and element-origin
                                                ! to all elements in the same background element
INTEGER,ALLOCATABLE   :: ListDistance(:)        ! the corresponding element id

TYPE tTrackingInfo
  INTEGER             :: CurrElem
  INTEGER             :: LocSide
  INTEGER             :: GlobSide
  INTEGER             :: flip
  INTEGER             :: TriNum
  REAL                :: xi
  REAL                :: eta
  REAL                :: alpha
  REAL                :: PartTrajectory(1:3)
  REAL                :: LengthPartTrajectory
END TYPE

TYPE(tTrackingInfo)   :: TrackInfo

#ifdef CODE_ANALYZE
INTEGER               :: PartOut
INTEGER               :: MPIRankOut
#endif /*CODE_ANALYZE*/
!===================================================================================================================================


END MODULE MOD_Particle_Tracking_Vars
