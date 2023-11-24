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


#if (PP_TimeDiscMethod==4)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStep_DSMC
!===================================================================================================================================
CONTAINS

SUBROUTINE TimeStep_DSMC()
!===================================================================================================================================
!> Direct Simulation Monte Carlo
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_TimeDisc_Vars            ,ONLY: dt, IterDisplayStep, iter, TEnd, Time
#ifdef PARTICLES
USE MOD_Globals                  ,ONLY: abort, CROSS
USE MOD_Particle_Vars            ,ONLY: PartState, LastPartPos, PDM, PEM, DoSurfaceFlux, WriteMacroVolumeValues
USE MOD_Particle_Vars            ,ONLY: UseRotRefFrame, RotRefFrameOmega, PartVeloRotRef
USE MOD_Particle_Vars            ,ONLY: WriteMacroSurfaceValues, Symmetry, Species, PartSpecies
USE MOD_Particle_Vars            ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_Particle_Vars            ,ONLY: UseSplitAndMerge
USE MOD_DSMC_Vars                ,ONLY: DSMC, CollisMode, AmbipolElecVelo
USE MOD_DSMC                     ,ONLY: DSMC_main
USE MOD_part_tools               ,ONLY: UpdateNextFreePosition
USE MOD_part_emission            ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux        ,ONLY: ParticleSurfaceflux
USE MOD_Particle_Tracking_vars   ,ONLY: tTracking,MeasureTrackTime
USE MOD_Particle_Tracking        ,ONLY: PerformTracking
USE MOD_SurfaceModel_Porous      ,ONLY: PorousBoundaryRemovalProb_Pressure
USE MOD_SurfaceModel_Vars        ,ONLY: nPorousBC
USE MOD_vMPF                     ,ONLY: SplitAndMerge
USE MOD_part_RHS                 ,ONLY: CalcPartRHSRotRefFrame
USE MOD_part_pos_and_velo        ,ONLY: SetParticleVelocity
USE MOD_Part_Tools               ,ONLY: InRotRefFrameCheck
USE MOD_Part_Tools               ,ONLY: CalcPartSymmetryPos
#if USE_MPI
USE MOD_Particle_MPI             ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/
USE MOD_DSMC_ParticlePairing     ,ONLY: GeoCoordToMap2D
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iPart
REAL                       :: timeEnd, timeStart, dtVar, RandVal
! Rotational frame of reference
REAL                       :: Pt_local(1:3), Pt_local_old(1:3), VeloRotRef_half(1:3), PartState_half(1:3)
#if USE_LOADBALANCE
REAL                  :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
! for reservoir simulation: no surface flux, particle push, tracking, ...
IF (DSMC%ReservoirSimu) THEN ! fix grid should be defined for reservoir simu
  CALL UpdateNextFreePosition()
  CALL DSMC_main()
  IF(DSMC%CompareLandauTeller) THEN
    DO iPart=1,PDM%ParticleVecLength
      PDM%nextFreePosition(iPart)=iPart
    END DO
    CALL SetParticleVelocity(1,1,PDM%ParticleVecLength)
  END IF
  RETURN
END IF

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF (DoSurfaceFlux) THEN
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_SURF,tLBStart)
#endif /*USE_LOADBALANCE*/

  CALL ParticleSurfaceflux()
END IF

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Variable time step: getting the right time step for the particle (can be constant across an element)
    IF (UseVarTimeStep) THEN
      dtVar = dt * PartTimeStep(iPart)
    ELSE
      dtVar = dt
    END IF
    IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(PartSpecies(iPart))%TimeStepFactor
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
    IF(DSMC%DoAmbipolarDiff.AND.(Species(PartSpecies(iPart))%ChargeIC.GT.0.0)) THEN
      CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart),AmbipolElecVelo(iPart)%ElecVelo)
    ELSE
      CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))
    END IF
  END IF
END DO

#if USE_LOADBALANCE
CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
! open receive buffer for number of particles
CALL IRecvNbOfParticles()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
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
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
! send number of particles
CALL SendNbOfParticles()
! finish communication of number of particles and send particles
CALL MPIParticleSend()
! finish communication
CALL MPIParticleRecv()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL ParticleInserting()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_EMISSION,tLBStart)
#endif /*USE_LOADBALANCE*/

IF (CollisMode.NE.0) THEN
  CALL UpdateNextFreePosition()
ELSE IF ( (MOD(iter,IterDisplayStep).EQ.0) .OR. &
          (Time.ge.(1-DSMC%TimeFracSamp)*TEnd) .OR. &
          WriteMacroVolumeValues.OR.WriteMacroSurfaceValues ) THEN
  CALL UpdateNextFreePosition() !postpone UNFP for CollisMode=0 to next IterDisplayStep or when needed for DSMC-Sampling
END IF

IF(DSMC%UseOctree)THEN
  ! Case Symmetry%Order=1 is performed in DSMC main
  IF(Symmetry%Order.EQ.2)THEN
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Store reference position in LastPartPos array to reduce memory demand
        CALL GeoCoordToMap2D(PartState(1:2,iPart), LastPartPos(1:2,iPart), PEM%LocalElemID(iPart))
      END IF
    END DO
  ELSEIF(Symmetry%Order.EQ.3) THEN
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Store reference position in LastPartPos array to reduce memory demand
        CALL GetPositionInRefElem(PartState(1:3,iPart),LastPartPos(1:3,iPart),PEM%GlobalElemID(iPart))
      END IF
    END DO
  END IF ! Symmetry%Order.EQ.2
END IF ! DSMC%UseOctree

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

CALL DSMC_main()

! Split & Merge: Variable particle weighting
IF(UseSplitAndMerge) CALL SplitAndMerge()

END SUBROUTINE TimeStep_DSMC


END MODULE MOD_TimeStep
#endif
