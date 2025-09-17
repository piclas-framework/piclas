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


#if (PP_TimeDiscMethod==701)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStep_ExplicitFV
!===================================================================================================================================

CONTAINS


SUBROUTINE TimeStep_ExplicitFV()
!===================================================================================================================================
! Explicit timestep with finite volumes
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars               ,ONLY: U_FV,Ut_FV
USE MOD_TimeDisc_Vars         ,ONLY: dt,time,iter
USE MOD_FV                    ,ONLY: FV_main
USE MOD_Equation_Vars_FV      ,ONLY:IniExactFunc_FV
#ifdef PARTICLES
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars          ,ONLY: PartState, Pt, LastPartPos, DelayTime,  PEM, PDM, Species, PartSpecies
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_PICModels              ,ONLY: FieldIonization
USE MOD_part_RHS               ,ONLY: CalcPartRHS
USE MOD_Ionization             ,ONLY: InsertNewIons
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition,isPushParticle
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_vMPF                   ,ONLY: SplitAndMerge
USE MOD_Particle_Vars          ,ONLY: UseSplitAndMerge
USE MOD_Particle_Vars          ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_Particle_Analyze_Tools ,ONLY: CalcCoupledPowerPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCoupledPower,PCoupl
USE MOD_Part_Tools             ,ONLY: CalcPartSymmetryPos
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux, DoForceFreeSurfaceFlux, DoFieldIonization
#endif /*PARTICLES*/
#if USE_HDG
USE MOD_HDG                    ,ONLY: HDG
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iPart, iElem
REAL                       :: RandVal, dtFrac, dtVar
#if USE_LOADBALANCE
REAL                       :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#ifdef PARTICLES
IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL Deposition()
IF ((IniExactFunc_FV.EQ.3).AND.(iter.EQ.0)) CALL InsertNewIons(init=IniExactFunc_FV)
#endif /*PARTICLES*/

! Electric field calculation
#if USE_HDG
CALL HDG(time,iter)
#endif

#ifdef PARTICLES
! set last data already here, since surfaceflux moved before interpolation
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
LastPartPos(1,1:PDM%ParticleVecLength)=PartState(1,1:PDM%ParticleVecLength)
LastPartPos(2,1:PDM%ParticleVecLength)=PartState(2,1:PDM%ParticleVecLength)
LastPartPos(3,1:PDM%ParticleVecLength)=PartState(3,1:PDM%ParticleVecLength)
PEM%LastGlobalElemID(1:PDM%ParticleVecLength)=PEM%GlobalElemID(1:PDM%ParticleVecLength)
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
  CALL InterpolateFieldToParticle()   ! forces on particles
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
  IF(DoFieldIonization) CALL FieldIonization()
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
      END IF
      ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
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
      CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))
      ! If coupled power output is active and particle carries charge, calculate energy difference and add to output variable
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'after')
    END IF
  END DO
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  CALL PerformTracking()
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  CALL ParticleInserting()
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
END IF

#endif /*PARTICLES*/

! Calculation of the electron density time derivative
CALL FV_main(time,time,doSource=.TRUE.)

#ifdef PARTICLES
IF ((time.GE.DelayTime).OR.(iter.EQ.0)) CALL UpdateNextFreePosition()

IF (time.GE.DelayTime) THEN
  ! Direct Simulation Monte Carlo
  IF (useDSMC) THEN
    CALL DSMC_main()
  END IF

  ! Particle emission for ionization
  CALL InsertNewIons()

  ! Split & Merge: Variable particle weighting
  IF(UseSplitAndMerge) CALL SplitAndMerge()
END IF
#endif /*PARTICLES*/

! Electron density update
U_FV = U_FV + Ut_FV*dt

! Prevent negative densities
DO iElem=1,PP_nElems
  U_FV(1,iElem) = MAX(U_FV(1,iElem),0.)
END DO

END SUBROUTINE TimeStep_ExplicitFV


END MODULE MOD_TimeStep
#endif /*PP_TimeDiscMethod==701*/