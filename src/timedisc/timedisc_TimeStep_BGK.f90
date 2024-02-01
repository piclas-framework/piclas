!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"


#if (PP_TimeDiscMethod==400)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStep_BGK
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStep_BGK()
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt, IterDisplayStep, iter, TEnd, Time
USE MOD_Globals                ,ONLY: abort, CROSS
USE MOD_Particle_Vars          ,ONLY: PartState, LastPartPos, PDM, PEM, DoSurfaceFlux, WriteMacroVolumeValues, Symmetry
USE MOD_Particle_Vars          ,ONLY: UseRotRefFrame, RotRefFrameOmega, PartVeloRotRef
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_DSMC_Vars              ,ONLY: DSMC, CollisMode
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_Particle_Tracking_vars ,ONLY: tTracking,MeasureTrackTime
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_part_RHS               ,ONLY: CalcPartRHSRotRefFrame
USE MOD_Part_Tools             ,ONLY: InRotRefFrameCheck
USE MOD_Part_Tools             ,ONLY: CalcPartSymmetryPos
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars      ,ONLY: DoParticleLatencyHiding
USE MOD_Globals                ,ONLY: CollectiveStop
#endif /*USE_MPI*/
USE MOD_BGK                    ,ONLY: BGK_main, BGK_DSMC_main
USE MOD_BGK_Vars               ,ONLY: CoupledBGKDSMC,DoBGKCellAdaptation
USE MOD_SurfaceModel_Porous    ,ONLY: PorousBoundaryRemovalProb_Pressure
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_DSMC_ParticlePairing   ,ONLY: GeoCoordToMap2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: timeEnd, timeStart
INTEGER               :: iPart
REAL                  :: RandVal, dtVar
! Rotational frame of reference
REAL                  :: Pt_local(1:3), Pt_local_old(1:3), VeloRotRef_half(1:3), PartState_half(1:3)
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/

IF (DoSurfaceFlux) THEN
  CALL ParticleSurfaceflux()
END IF

DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Variable time step: getting the right time step for the particle (can be constant across an element)
    IF (UseVarTimeStep) THEN
      dtVar = dt * PartTimeStep(iPart)
    ELSE
      dtVar = dt
    END IF
    IF (PDM%dtFracPush(iPart)) THEN
      ! Surface flux (dtFracPush): previously inserted particles are pushed a random distance between 0 and v*dt
      !                            LastPartPos and LastElem already set!
      CALL RANDOM_NUMBER(RandVal)
      dtVar = dtVar * RandVal
      PDM%dtFracPush(iPart) = .FALSE.
    ELSE
      LastPartPos(1:3,iPart)=PartState(1:3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
    END IF
    IF(UseRotRefFrame) THEN
      IF(PDM%InRotRefFrame(iPart)) THEN
        ! Midpoint method
        ! calculate the acceleration (force / mass) at the current time step
        Pt_local_old(1:3) = CalcPartRHSRotRefFrame(PartState(1:3,iPart), PartVeloRotRef(1:3,iPart))
        ! estimate the mid-point velocity in the rotational frame
        VeloRotRef_half(1:3) = PartVeloRotRef(1:3,iPart) + 0.5*Pt_local_old(1:3)*dtVar
        ! estimate the mid-point position
        PartState_half(1:3) = PartState(1:3,iPart) + 0.5*PartVeloRotRef(1:3,iPart)*dtVar
        ! calculate the acceleration (force / mass) at the mid-point
        Pt_local(1:3) = CalcPartRHSRotRefFrame(PartState_half(1:3), VeloRotRef_half(1:3))
        ! update the position using the mid-point velocity in the rotational frame
        PartState(1:3,iPart) = PartState(1:3,iPart) + VeloRotRef_half(1:3)*dtVar
        ! update the velocity in the rotational frame using the mid-point acceleration
        PartVeloRotRef(1:3,iPart) = PartVeloRotRef(1:3,iPart) + Pt_local(1:3)*dtVar
      ELSE
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtVar
      END IF
    ELSE
      PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtVar
    END IF
    ! Axisymmetric treatment of particles: rotation of the position and velocity vector
    CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))
  END IF
END DO

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#if USE_MPI
! open receive buffer for number of particles
CALL IRecvNbOfParticles()
#endif /*USE_MPI*/
IF(MeasureTrackTime) CALL CPU_TIME(TimeStart)
! actual tracking
CALL PerformTracking()
IF (nPorousBC.GT.0) THEN
  CALL PorousBoundaryRemovalProb_Pressure()
END IF
IF(MeasureTrackTime) THEN
  CALL CPU_TIME(TimeEnd)
  tTracking=tTracking+TimeEnd-TimeStart
END IF
#if USE_MPI
! send number of particles
CALL SendNbOfParticles()
#endif
CALL ParticleInserting()
IF (CollisMode.NE.0) THEN
  CALL UpdateNextFreePosition(.TRUE.)
ELSE IF ( (MOD(iter,IterDisplayStep).EQ.0) .OR. &
          (Time.ge.(1-DSMC%TimeFracSamp)*TEnd) .OR. &
          WriteMacroVolumeValues ) THEN
  CALL UpdateNextFreePosition(.TRUE.) !postpone UNFP for CollisMode=0 to next IterDisplayStep or when needed for DSMC-Sampling
END IF

#if USE_MPI
! finish communication of number of particles and send particles
CALL MPIParticleSend(.TRUE.)
#endif /*USE_MPI*/

IF(DoBGKCellAdaptation.OR.(CoupledBGKDSMC.AND.DSMC%UseOctree)) THEN
  IF(Symmetry%Order.EQ.2) THEN
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Store reference position in LastPartPos array to reduce memory demand
        CALL GeoCoordToMap2D(PartState(1:2,iPart), LastPartPos(1:2,iPart), PEM%LocalElemID(iPart))
      END IF
    END DO
  ELSE
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Store reference position in LastPartPos array to reduce memory demand
        CALL GetPositionInRefElem(PartState(1:3,iPart),LastPartPos(1:3,iPart),PEM%GlobalElemID(iPart))
      END IF
    END DO
  END IF ! Symmetry%Order.EQ.2
END IF ! DoBGKCellAdaptation.OR.(CoupledBGKDSMC.AND.DSMC%UseOctree)

IF(UseRotRefFrame) THEN
  DO iPart = 1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(iPart)) THEN
      IF(InRotRefFrameCheck(iPart)) THEN
        ! Particle moved into the rotational frame of reference, initialize velocity
        IF(.NOT.PDM%InRotRefFrame(iPart)) THEN
          PartVeloRotRef(1:3,iPart) = PartState(4:6,iPart) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,iPart))
        END IF
      ELSE
        ! Particle left (or never was in) the rotational frame of reference
        PartVeloRotRef(1:3,iPart) = 0.
      END IF
      PDM%InRotRefFrame(iPart) = InRotRefFrameCheck(iPart)
    END IF
  END DO
END IF

#if USE_MPI
IF(DoParticleLatencyHiding)THEN
  ! IF (CoupledBGKDSMC) THEN
  !   CALL BGK_DSMC_main(1)
  ! ELSE
    CALL BGK_main(1)
  ! END IF
END IF ! DoParticleLatencyHiding

! finish communication
CALL MPIParticleRecv(.TRUE.)
#endif /*USE_MPI*/

! After MPI communication of particles, call UNFP including the MPI particles
CALL UpdateNextFreePosition()

!#ifdef EXTRAE
!CALL extrae_eventandcounters(int(9000001), int8(51))
!#endif /*EXTRAE*/

!#ifdef EXTRAE
!CALL extrae_eventandcounters(int(9000001), int8(0))
!#endif /*EXTRAE*/
#if USE_MPI
IF(DoParticleLatencyHiding)THEN
  ! IF (CoupledBGKDSMC) THEN
  !   CALL BGK_DSMC_main(2)
  ! ELSE
    CALL BGK_main(2)
  ! END IF
ELSE
#endif /*USE_MPI*/
  IF (CoupledBGKDSMC) THEN
    CALL BGK_DSMC_main()
  ELSE
    CALL BGK_main()
  END IF
#if USE_MPI
END IF ! DoParticleLatencyHiding
#endif /*USE_MPI*/

!#ifdef EXTRAE
!CALL extrae_eventandcounters(int(9000001), int8(52))
!#endif /*EXTRAE*/



!#ifdef EXTRAE
!CALL extrae_eventandcounters(int(9000001), int8(0))
!#endif /*EXTRAE*/

END SUBROUTINE TimeStep_BGK


END MODULE MOD_TimeStep
#endif /*(PP_TimeDiscMethod==400)*/
