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


#if (PP_TimeDiscMethod==42)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStep_DSMC_Debug
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStep_DSMC_Debug()
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals                ,ONLY: abort
USE MOD_TimeDisc_Vars          ,ONLY: dt
#ifdef PARTICLES
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: PartState, LastPartPos, PDM,PEM
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_part_pos_and_velo      ,ONLY: SetParticleVelocity
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_Particle_Tracking_vars ,ONLY: tTracking,MeasureTrackTime
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPart
REAL                  :: timeStart, timeEnd, RandVal, dtFrac
!===================================================================================================================================
IF(DSMC%UseOctree) CALL abort(__STAMP__,'Particles-DSMC-UseOctree = T not allowed for RESERVOIR simulation')

IF (DSMC%ReservoirSimu) THEN ! fix grid should be defined for reservoir simu

  CALL UpdateNextFreePosition()

  CALL DSMC_main()

  IF(DSMC%CompareLandauTeller) THEN
    DO iPart=1,PDM%ParticleVecLength
      PDM%nextFreePosition(iPart)=iPart
    END DO
    CALL SetParticleVelocity(1,1,PDM%ParticleVecLength)
  END IF
ELSE
  IF (DoSurfaceFlux) THEN
    CALL ParticleSurfaceflux()
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%dtFracPush(iPart)) THEN
          LastPartPos(1:3,iPart)=PartState(1:3,iPart)
          PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
          PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dt
        ELSE !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RandVal
          PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtFrac
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
    END DO
  ELSE
    LastPartPos(1:3,1:PDM%ParticleVecLength)=PartState(1:3,1:PDM%ParticleVecLength)
    ! bugfix if more than 2.x mio (2000001) particle per process
    ! tested with 64bit Ubuntu 12.04 backports
    DO iPart = 1, PDM%ParticleVecLength
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
    END DO
    !PEM%LastGlobalElemID(1:PDM%ParticleVecLength)=PEM%GlobalElemID(1:PDM%ParticleVecLength)
    PartState(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength) + PartState(4:6,1:PDM%ParticleVecLength) * dt
  END IF
#if USE_MPI
  ! open receive buffer for number of particles
  CALL IRecvNbOfParticles()
#endif /*USE_MPI*/
  IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
  ! actual tracking
  CALL PerformTracking()
  IF(MeasureTrackTime) THEN
    CALL CPU_TIME(TimeEnd)
    tTracking=tTracking+TimeEnd-TimeStart
  END IF
#if USE_MPI
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
#endif /*USE_MPI*/
  CALL ParticleInserting()
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
END IF

END SUBROUTINE TimeStep_DSMC_Debug


END MODULE MOD_TimeStep
#endif
