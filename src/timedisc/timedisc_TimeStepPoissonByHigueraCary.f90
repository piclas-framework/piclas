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
#if (PP_TimeDiscMethod==507)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepPoissonByHigueraCary
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStepPoissonByHigueraCary()
!===================================================================================================================================
! Higuera-Cary (507) -push with HDG: relativistic push 2nd order
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: c_inv,c2_inv
USE MOD_Globals                ,ONLY: Abort, LocalTime, CROSS, DOTPRODUCT, UNITVECTOR, VECNORM
USE MOD_DG_Vars                ,ONLY: U
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iter,time
USE MOD_HDG                    ,ONLY: HDG
#ifdef PARTICLES
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars          ,ONLY: PartState, LastPartPos,PEM, PDM, DelayTime
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies
USE MOD_Particle_Analyze_Tools ,ONLY: CalcCoupledPowerPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCoupledPower,PCoupl
!#if (PP_TimeDiscMethod==509)
!USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
!#endif /*(PP_TimeDiscMethod==509)*/
USE MOD_part_RHS               ,ONLY: CalcPartRHS, CalcPartRHSSingleParticle
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_Part_Tools             ,ONLY: CalcPartSymmetryPos
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange
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
#if USE_LOADBALANCE
REAL                       :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
REAL                       :: RandVal, dtFrac
REAL                       :: gamma1,c1,gammaMinus,gammaPlus,tau2,sigma,s,uStar
REAL, DIMENSION(3)         :: tauVec, tVec, uMinus, uPlus
#endif /*PARTICLES*/
!===================================================================================================================================
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
  IF (DoSurfaceFlux) CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! 1.  Update position by half step before calculating the fields
  IF (CalcCoupledPower) PCoupl = 0. ! if output of coupled power is active: reset PCoupl
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
      ! 1st part of the position update (also neutral particles)
      !-- x(n) => x(n+1/2) by v(n) = u(n)/gamma(n):
      PartState(1:3,iPart) = PartState(1:3,iPart) + 0.5 * PartState(4:6,iPart) * dt
    END IF ! PDM%ParticleInside(iPart)
  END DO ! iPart=1,PDM%ParticleVecLength
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
#ifdef EXTRAE
  CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/

! Tracking after 1st half push
#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  CALL PerformTracking()
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
#endif
END IF ! time.GE.DelayTime

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL Deposition()
#ifdef EXTRAE
  CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#endif /*PARTICLES*/

! HDG at x(n+1/2) gives E(x(n+1/2))
CALL HDG(time,U,iter)

#ifdef PARTICLES
IF (time.GE.DelayTime) THEN
  LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength)
  PEM%LastGlobalElemID(1:PDM%ParticleVecLength) = PEM%GlobalElemID(1:PDM%ParticleVecLength)
#ifdef EXTRAE
  CALL extrae_eventandcounters(int(9000001), int8(5))
#endif /*EXTRAE*/
  !-- get E(x(n+1/2)) and B(x(n+1/2))
  CALL InterpolateFieldToParticle()   ! forces on particles
  !CALL InterpolateFieldToParticle() ! only needed when MPI communication changes the number of parts
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/

  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !DoSurfaceFlux for compiler-optimization if .FALSE.
        CALL RANDOM_NUMBER(RandVal)
        dtFrac = dt * RandVal
        PDM%dtFracPush(iPart) = .FALSE.
      ELSE
        dtFrac = dt
      END IF

      ! Don't push the velocity component of neutral particles!
      IF(isPushParticle(iPart).AND.DoInterpolation)THEN
        !-- get Lorentz factor gamma1(n)
        gamma1=DOTPRODUCT(PartState(4:6,iPart))*c2_inv
        IF(gamma1.GE.1.0) CALL abort(__STAMP__,'Velocity is geater than c',RealInfoOpt=SQRT(DOTPRODUCT(PartState(4:6,iPart))))
        gamma1=1./SQRT(1.0-gamma1)

        !-- u(n) = v(n)*gamma(n) (store vector u in v)
        PartState(4:6,iPart) = PartState(4:6,iPart)*gamma1

        !-- const. factor
        c1 = (Species(PartSpecies(iPart))%ChargeIC * dtFrac) / (Species(PartSpecies(iPart))%MassIC * 2.)

        !-- u- = u(n) + q/m*E(n+1/2)*dt/2
        uMinus = PartState(4:6,iPart) + c1 * FieldAtParticle(1:3,iPart)

        !-- Auxiliary variables
        ! gamma- = SQRT(1+(|u-/c|)^2)
        gammaMinus = SQRT(1.0 + DOTPRODUCT(uMinus)*c2_inv)
        ! tau = q/m*B(n+1/2)*dt/2
        tauVec = c1 * FieldAtParticle(4:6,iPart)
        ! u* = (u- * tau)/c
        uStar = DOT_PRODUCT(uMinus, tauVec)*c_inv
        ! sigma = (gamma-)^2 - tau^2
        tau2 = DOTPRODUCT(tauVec)
        sigma = gammaMinus**2 - tau2
        ! gamma+ = SQRT(...)
        gammaPlus = SQRT(sigma**2 + 4.0*(tau2 + uStar**2))
        gammaPlus = SQRT(0.5*(sigma+gammaPlus))
        ! t = tau/gamma+
        tVec = tauVec / gammaPlus
        ! s = 1/(1 + t^2)
        s = 1.0/(1.0 + DOTPRODUCT(tVec))

        !-- u+ = s[u- + (u- * t)t + u- x t)]
        uPlus = s*(uMinus + DOT_PRODUCT(uMinus, tVec)*tVec + CROSS(uMinus, tVec))

        !-- u(n+1) = u+ + q/m*E(n+1/2)*dt/2 + u- x t
        PartState(4:6,iPart) = uPlus + c1 * FieldAtParticle(1:3,iPart) + CROSS(uPlus, tVec)

        !-- v(n+1) = u(n+1)/gamma(n+1)
        PartState(4:6,iPart) = PartState(4:6,iPart)/SQRT(1+DOTPRODUCT(PartState(4:6,iPart))*c2_inv)
      END IF

      ! 2nd part of the position update (also neutral particles)
      !-- x(n) => x(n+1) by v(n+0.5):
      PartState(1:3,iPart) = PartState(1:3,iPart) + 0.5 * PartState(4:6,iPart) * dtFrac

      CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))

      ! If coupled power output is active and particle carries charge, calculate energy difference and add to output variable
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'after')
    END IF ! PDM%ParticleInside(iPart)
  END DO ! iPart=1,PDM%ParticleVecLength
#ifdef EXTRAE
  CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

! Tracking after 2nd half push
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
END IF ! time.GE.DelayTime

#if USE_MPI
PartMPIExchange%nMPIParticles=0 ! and set number of received particles to zero for deposition
#endif

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

END SUBROUTINE TimeStepPoissonByHigueraCary


END MODULE MOD_TimeStep
#endif /*(PP_TimeDiscMethod==507)*/
#endif /*USE_HDG*/
