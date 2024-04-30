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


#if USE_HDG
#if (PP_TimeDiscMethod==500) || (PP_TimeDiscMethod==509)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepPoisson
!===================================================================================================================================
CONTAINS

SUBROUTINE TimeStepPoisson()
!===================================================================================================================================
! Euler (500) or Leapfrog (509) -push with HDG
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort, LocalTime
USE MOD_DG_Vars                ,ONLY: U
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iter,time
!#if (PP_TimeDiscMethod==509)
!USE MOD_TimeDisc_Vars          ,ONLY: dt_old
!#endif /*(PP_TimeDiscMethod==509)*/
USE MOD_HDG                    ,ONLY: HDG
#ifdef PARTICLES
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars          ,ONLY: PartState, Pt, LastPartPos,PEM, PDM, DelayTime, Species, PartSpecies
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux, DoForceFreeSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_Particle_Analyze_Tools ,ONLY: CalcCoupledPowerPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCoupledPower,PCoupl
#if (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
#endif /*(PP_TimeDiscMethod==509)*/
USE MOD_part_RHS               ,ONLY: CalcPartRHS
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_Part_Tools             ,ONLY: CalcPartSymmetryPos
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition,isPushParticle
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_vMPF                   ,ONLY: SplitAndMerge
USE MOD_Particle_Vars          ,ONLY: UseSplitAndMerge
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iPart
REAL                       :: RandVal, dtFrac, dtVar
#if USE_LOADBALANCE
REAL                       :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#ifdef PARTICLES
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL Deposition()
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#endif /*PARTICLES*/

CALL HDG(time,U,iter)

#ifdef PARTICLES
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
! set last data already here, since surface flux moved before interpolation
LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength)
PEM%LastGlobalElemID(1:PDM%ParticleVecLength) = PEM%GlobalElemID(1:PDM%ParticleVecLength)
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
IF (time.GE.DelayTime) THEN
  IF (DoSurfaceFlux) THEN
    CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
  END IF
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  !-- get E(x(n))
  CALL InterpolateFieldToParticle()   ! forces on particles
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
  !-- get a(x(n))
  IF(DoInterpolation) CALL CalcPartRHS()
  IF (CalcCoupledPower) PCoupl = 0. ! if output of coupled power is active: reset PCoupl
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      ! Set the particle time step
      IF (UseVarTimeStep) THEN
        dtVar = dt * PartTimeStep(iPart)
      ELSE
        dtVar = dt
      END IF
      ! Set the species-specific time step
      IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(PartSpecies(iPart))%TimeStepFactor
      IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !DoSurfaceFlux for compiler-optimization if .FALSE.
        CALL RANDOM_NUMBER(RandVal)
        dtFrac = dtVar * RandVal
        PDM%IsNewPart(iPart)=.FALSE. !no IsNewPart-treatment for surffluxparts
      ELSE
        dtFrac = dtVar
#if (PP_TimeDiscMethod==509)
        IF (PDM%IsNewPart(iPart)) THEN
          ! Don't push the velocity component of neutral particles!
          IF(isPushParticle(iPart))THEN
            !-- v(n) => v(n-0.5) by a(n):
            PartState(4:6,iPart) = PartState(4:6,iPart) - Pt(1:3,iPart) * dtVar*0.5
          END IF
          PDM%IsNewPart(iPart)=.FALSE. !IsNewPart-treatment is now done
        ELSE
          !IF ((ABS(dt-dt_old).GT.1.0E-6*dt_old).AND.&
          !     isPushParticle(iPart)) THEN ! Don't push the velocity component of neutral particles!
          !  PartState(4:6,iPart)  = PartState(4:6,iPart) + Pt(1:3,iPart) * (dt_old-dt)*0.5
          !END IF
        END IF
#endif /*(PP_TimeDiscMethod==509)*/
      END IF
      ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
#if (PP_TimeDiscMethod==509)
      IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart) .AND. .NOT.DoForceFreeSurfaceFlux) THEN
        !-- x(BC) => x(n+1) by v(BC+X):
        PartState(1:3,iPart) = PartState(1:3,iPart) + ( PartState(4:6,iPart) + Pt(1:3,iPart) * dtFrac*0.5 ) * dtFrac
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          !-- v(BC) => v(n+0.5) by a(BC):
          PartState(4:6,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart) * (dtFrac - dtVar*0.5)
        END IF
        PDM%dtFracPush(iPart) = .FALSE.
      ELSE IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !DoForceFreeSurfaceFlux
        !-- x(n) => x(n+1) by v(n+0.5)=v(BC)
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtFrac
        PDM%dtFracPush(iPart) = .FALSE.
      ELSE
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          !-- v(n-0.5) => v(n+0.5) by a(n):
          PartState(4:6,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart) * dtVar
        END IF
        !-- x(n) => x(n+1) by v(n+0.5):
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtVar
      END IF
#else /*(PP_TimeDiscMethod==509)*/
        !-- x(n) => x(n+1) by v(n):
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtFrac
      IF (DoForceFreeSurfaceFlux .AND. DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN
        PDM%dtFracPush(iPart) = .FALSE.
      ELSE
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          !-- v(n) => v(n+1) by a(n):
          PartState(4:6,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart) * dtFrac
        END IF
        IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
#endif /*(PP_TimeDiscMethod==509)*/
      CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))
      ! If coupled power output is active and particle carries charge, calculate energy difference and add to output variable
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'after')
    END IF
  END DO
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  CALL PerformTracking()
  CALL ParticleInserting()
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
#endif
#if (PP_TimeDiscMethod==509)
  IF (velocityOutputAtTime) THEN
#ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
    CALL Deposition() ! because of emission and UpdateParticlePosition
#ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
    CALL HDG(time,U,iter)
#ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    CALL InterpolateFieldToParticle()   ! forces on particles
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
    IF(DoInterpolation) CALL CalcPartRHS()
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! Set the particle time step
        IF (UseVarTimeStep) THEN
          dtVar = dt * PartTimeStep(iPart)
        ELSE
          dtVar = dt
        END IF
        ! Set the species-specific time step
        IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(PartSpecies(iPart))%TimeStepFactor
        !-- v(n+0.5) => v(n+1) by a(n+1):
        velocityAtTime(1:3,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart) * dtVar*0.5
      END IF
    END DO
#ifdef EXTRAE
    CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF !velocityOutputAtTime
#endif /*(PP_TimeDiscMethod==509)*/
END IF

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/

IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL UpdateNextFreePosition()

IF (time.GE.DelayTime) THEN
  ! Direct Simulation Monte Carlo
  IF (useDSMC) THEN
    CALL DSMC_main()
  END IF
  ! Split & Merge: Variable particle weighting
  IF(UseSplitAndMerge) CALL SplitAndMerge()
END IF
#endif /*PARTICLES*/

END SUBROUTINE TimeStepPoisson


END MODULE MOD_TimeStep
#endif /*(PP_TimeDiscMethod==500) || (PP_TimeDiscMethod==509)*/
#endif /*USE_HDG*/
