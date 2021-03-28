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
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: PartState, LastPartPos, PDM, PEM, DoSurfaceFlux, WriteMacroVolumeValues
USE MOD_Particle_Vars          ,ONLY: VarTimeStep, Symmetry
USE MOD_DSMC_Vars              ,ONLY: DSMC_RHS, DSMC, CollisMode
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_Particle_Tracking_vars ,ONLY: tTracking,MeasureTrackTime
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
USE MOD_BGK                    ,ONLY: BGK_main, BGK_DSMC_main
USE MOD_BGK_Vars               ,ONLY: CoupledBGKDSMC
USE MOD_SurfaceModel_Porous    ,ONLY: PorousBoundaryRemovalProb_Pressure
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: timeEnd, timeStart
INTEGER               :: iPart
REAL                  :: RandVal, dtVar, NewYPart, NewYVelo
!===================================================================================================================================
IF (DoSurfaceFlux) THEN
  CALL ParticleSurfaceflux()
END IF

DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
  ! Variable time step: getting the right time step for the particle (can be constant across an element)
  IF (VarTimeStep%UseVariableTimeStep) THEN
    dtVar = dt * VarTimeStep%ParticleTimeStep(iPart)
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
  PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtVar
  ! Axisymmetric treatment of particles: rotation of the position and velocity vector
  IF(Symmetry%Axisymmetric) THEN
    IF (PartState(2,iPart).LT.0.0) THEN
      NewYPart = -SQRT(PartState(2,iPart)**2 + (PartState(3,iPart))**2)
    ELSE
      NewYPart = SQRT(PartState(2,iPart)**2 + (PartState(3,iPart))**2)
    END IF
    ! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
    !           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
    ! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
    ! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
    NewYVelo = (PartState(5,iPart)*(PartState(2,iPart))+PartState(6,iPart)*PartState(3,iPart))/NewYPart
    PartState(6,iPart) = (-PartState(5,iPart)*PartState(3,iPart)+PartState(6,iPart)*(PartState(2,iPart)))/NewYPart
    PartState(2,iPart) = NewYPart
    PartState(3,iPart) = 0.0
    PartState(5,iPart) = NewYVelo
    END IF
  END IF
END DO

! Resetting the particle positions in the second/third dimension for the 1D/2D/axisymmetric case
IF(Symmetry%Order.LT.3) THEN
  LastPartPos(Symmetry%Order+1:3,1:PDM%ParticleVecLength) = 0.0
  PartState(Symmetry%Order+1:3,1:PDM%ParticleVecLength) = 0.0
END IF

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
! finish communication of number of particles and send particles
CALL MPIParticleSend()
! finish communication
CALL MPIParticleRecv()
#endif /*USE_MPI*/
CALL ParticleInserting()
IF (CollisMode.NE.0) THEN
  CALL UpdateNextFreePosition()
ELSE IF ( (MOD(iter,IterDisplayStep).EQ.0) .OR. &
          (Time.ge.(1-DSMC%TimeFracSamp)*TEnd) .OR. &
          WriteMacroVolumeValues ) THEN
  CALL UpdateNextFreePosition() !postpone UNFP for CollisMode=0 to next IterDisplayStep or when needed for DSMC-Sampling
ELSE IF (PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).GT.PDM%maxParticleNumber .OR. &
         PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).EQ.0) THEN
  CALL abort(&
__STAMP__,&
'maximum nbr of particles reached!')  !gaps in PartState are not filled until next UNFP and array might overflow more easily!
END IF
  IF (CoupledBGKDSMC) THEN
    CALL BGK_DSMC_main()
  ELSE
    CALL BGK_main()
  END IF

PartState(4:6,1:PDM%ParticleVecLength) = PartState(4:6,1:PDM%ParticleVecLength) + DSMC_RHS(1:3,1:PDM%ParticleVecLength)

END SUBROUTINE TimeStep_BGK


END MODULE MOD_TimeStep
#endif
