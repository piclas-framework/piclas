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
#if (PP_TimeDiscMethod==508)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepPoissonByBorisLeapfrog
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStepPoissonByBorisLeapfrog()
!===================================================================================================================================
! Boris-Leapfrog (508) -push with HDG
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort, LocalTime, CROSS, DOTPRODUCT, UNITVECTOR, VECNORM, PARTISELECTRON
USE MOD_DG_Vars                ,ONLY: U
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iter,time
USE MOD_Globals_Vars           ,ONLY: c2_inv
USE MOD_HDG                    ,ONLY: HDG
#ifdef PARTICLES
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars          ,ONLY: PartState, Pt, LastPartPos,PEM, PDM, DelayTime
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_Particle_Analyze_Tools ,ONLY: CalcCoupledPowerPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCoupledPower,PCoupl
USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
USE MOD_part_RHS               ,ONLY: CalcPartRHS, CalcPartRHSSingleParticle
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
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iPart
REAL                       :: RandVal, dtFrac, gamma, gamma_minus, dtVar
#if USE_LOADBALANCE
REAL                       :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
REAL                       :: c_1
REAL, DIMENSION(3)         :: v_minus, v_plus, v_prime, t_vec, v_minus_old, v_n1
#endif /*PARTICLES*/
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
  !-- get E(x(n)) and B(x(n))
  CALL InterpolateFieldToParticle()   ! forces on particles
  !CALL InterpolateFieldToParticle() ! only needed when MPI communication changes the number of parts
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
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
        PDM%dtFracPush(iPart) = .FALSE.
      ELSE
        dtFrac = dtVar
        IF (PDM%IsNewPart(iPart)) THEN
          ! Don't push the velocity component of neutral particles!
          IF(isPushParticle(iPart).AND.DoInterpolation)THEN
              !-- Shift particle velocity back in time by half a time step dtVar/2.
              !-- get Pt(1:3,iPart) = a(x(n))
              CALL CalcPartRHSSingleParticle(iPart)

              !-- v(n) => v(n-0.5) by a(n):
              PartState(4:6,iPart) = PartState(4:6,iPart) - Pt(1:3,iPart) * dtVar*0.5
          END IF
          PDM%IsNewPart(iPart)=.FALSE. !IsNewPart-treatment is now done
        END IF
      END IF
      ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
      IF(isPushParticle(iPart).AND.DoInterpolation)THEN ! Don't push the velocity component of neutral particles!
        !-- v(n-0.5) => v(n+0.5) by a(n):
        !PartState(4:6,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart) * dt

        !-- const. factor
        c_1 =  (Species(PartSpecies(iPart))%ChargeIC * dtFrac) / (Species(PartSpecies(iPart))%MassIC * 2.)

        !-- v_minus = v(n-1/2) + q/m*E(n)*dt/2
        gamma = 1./SQRT(1-(DOTPRODUCT(PartState(4:6,iPart))*c2_inv))

        v_minus_old = PartState(4:6,iPart)*gamma

        v_minus = v_minus_old + c_1 * FieldAtParticle(1:3,iPart)

        !-- t_vec
        ! gamma_minus = (SQRT(1+DOTPRODUCT(v_minus_old)*c2_inv) + SQRT(1+DOTPRODUCT(v_minus)*c2_inv) ) / 2.
        gamma_minus = SQRT(1+DOTPRODUCT(v_minus)*c2_inv)
        t_vec = TAN(c_1/gamma_minus*VECNORM(FieldAtParticle(4:6,iPart))) * UNITVECTOR(FieldAtParticle(4:6,iPart))

        !-- v_prime = v_minus + v_minus x (q*B/m)*dt/2
        v_prime = v_minus + CROSS(v_minus, t_vec )

        !-- v_plus = v_minus + v_prime x 2*t_vec/(1+t_vec^2) where t_vec = c_1 * B
        v_plus = v_minus + (2.0/(1.+DOTPRODUCT(t_vec))) * CROSS(v_prime, t_vec)

        !-- v(n+1/2) = v_plus + c_1 * E
        v_n1 = v_plus + c_1 * FieldAtParticle(1:3,iPart)

        !-- v(n+1) = u(n+1)/gamma(n+1)
        PartState(4:6,iPart) = v_n1 / SQRT(1+DOTPRODUCT(v_n1)*c2_inv)
      END IF

      !-- x(n) => x(n+1) by v(n+0.5):
      PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * dtFrac

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
    !CALL InterpolateFieldToParticle() ! only needed when MPI communication changes the number of parts
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

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#endif /*PARTICLES*/

END SUBROUTINE TimeStepPoissonByBorisLeapfrog


END MODULE MOD_TimeStep
#endif /*(PP_TimeDiscMethod==508)*/
#endif /*USE_HDG*/
