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
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepPoissonByLSERK
!===================================================================================================================================
CONTAINS


SUBROUTINE TimeStepPoissonByLSERK()
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort, LocalTime
USE MOD_PreProc
USE MOD_TimeDisc_Vars          ,ONLY: dt,iStage,time,iter
USE MOD_TimeDisc_Vars          ,ONLY: RK_c,nRKStages
USE MOD_DG_Vars                ,ONLY: U
#ifdef PARTICLES
USE MOD_TimeDisc_Vars          ,ONLY: RK_a,RK_b,dt_Min,dtWeight
USE MOD_TimeDisc_Vars          ,ONLY: RKdtFracTotal,RKdtFrac
USE MOD_PICDepo                ,ONLY: Deposition
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars          ,ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime,  PEM, PDM
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux, DoForceFreeSurfaceFlux, DoFieldIonization
USE MOD_PICModels              ,ONLY: FieldIonization
USE MOD_part_RHS               ,ONLY: CalcPartRHS
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolation
USE MOD_part_emission          ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux      ,ONLY: ParticleSurfaceflux
USE MOD_DSMC                   ,ONLY: DSMC_main
USE MOD_DSMC_Vars              ,ONLY: useDSMC
USE MOD_Particle_Analyze_Tools ,ONLY: CalcCoupledPowerPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCoupledPower,PCoupl
USE MOD_Part_Tools             ,ONLY: CalcPartSymmetryPos
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition,isPushParticle
USE MOD_Particle_Tracking      ,ONLY: PerformTracking
USE MOD_vMPF                   ,ONLY: SplitAndMerge
USE MOD_Particle_Vars          ,ONLY: UseSplitAndMerge
#endif /*PARTICLES*/
USE MOD_HDG                    ,ONLY: HDG
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL           :: tStage
#ifdef PARTICLES
INTEGER        :: iPart, iStage_loc
REAL           :: RandVal,b_dt(1:nRKStages)
REAL           :: Pa_rebuilt_coeff(1:nRKStages),Pa_rebuilt(1:3,1:nRKStages),Pv_rebuilt(1:3,1:nRKStages),v_rebuilt(1:3,0:nRKStages-1)
#endif /*PARTICLES*/
#if defined(PARTICLES) && USE_LOADBALANCE
REAL           :: tLBStart
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
!===================================================================================================================================

#ifdef PARTICLES
DO iStage_loc=1,nRKStages
  ! RK coefficients
  b_dt(iStage_loc)=RK_b(iStage_loc)*dt
  ! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RK_a(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO
#endif /*PARTICLES*/

! first RK step
iStage=1
tStage=time
#ifdef PARTICLES
CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
RKdtFrac = RK_c(2)
dtWeight = dt/dt_Min(DT_MIN) * RKdtFrac
RKdtFracTotal=RKdtFrac

IF ((time.GE.DelayTime).OR.(iter.EQ.0)) THEN
  ! communicate shape function particles
  CALL Deposition() ! because of emission and UpdateParticlePosition
END IF
#endif /*PARTICLES*/

CALL HDG(tStage,U,iter)

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
      ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
      IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
      !-- Pt is not known only for new Surfaceflux-Parts -> change IsNewPart back to F for other Parts
      IF (.NOT.DoSurfaceFlux) THEN
        PDM%IsNewPart(iPart)=.FALSE.
      ELSE
        IF (.NOT.PDM%dtFracPush(iPart)) PDM%IsNewPart(iPart)=.FALSE.
      END IF
      !-- Particle Push
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(1:3,iPart) = PartState(4:6,iPart)
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart) * b_dt(1)
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          Pt_temp(4:6,iPart) = Pt(1:3,iPart)
          PartState(4:6,iPart) = PartState(4:6,iPart) + Pt(1:3,iPart)*b_dt(1)
        END IF
      ELSE !IsNewPart: no Pt_temp history available!
        IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage
          CALL RANDOM_NUMBER(RandVal)
          IF (DoForceFreeSurfaceFlux) Pt(1:3,iPart)=0.
        ELSE
          CALL abort(__STAMP__,'Error in LSERK-HDG-Timedisc: This case should be impossible...')
        END IF
        Pa_rebuilt(:,:)=0.
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          DO iStage_loc=1,iStage
            Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(1:3,iPart)
          END DO
        END IF
        v_rebuilt(:,:)=0.
        DO iStage_loc=iStage-1,0,-1
          IF (iStage_loc.EQ.iStage-1) THEN
            v_rebuilt(1:3,iStage_loc) = PartState(4:6,iPart) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          ELSE
            v_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc+1) - b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          END IF
        END DO
        Pv_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          IF (iStage_loc.EQ.1) THEN
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,0)
          ELSE
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc-1) - RK_a(iStage_loc)*Pv_rebuilt(1:3,iStage_loc-1)
          END IF
        END DO
        Pt_temp(1:3,iPart) = Pv_rebuilt(1:3,iStage)
        PartState(1:3,iPart) = PartState(1:3,iPart) + Pt_temp(1:3,iPart)*b_dt(iStage)*RandVal
        ! Don't push the velocity component of neutral particles!
        IF(isPushParticle(iPart))THEN
          Pt_temp(4:6,iPart) = Pa_rebuilt(1:3,iStage)
          PartState(4:6,iPart) = PartState(4:6,iPart) + Pt_temp(4:6,iPart)*b_dt(iStage)*RandVal
        END IF
        PDM%dtFracPush(iPart) = .FALSE.
        IF (.NOT.DoForceFreeSurfaceFlux) PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
      END IF !IsNewPart

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
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_EMISSION,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
END IF

#endif /*PARTICLES*/

! perform RK steps
DO iStage=2,nRKStages
  tStage=time+dt*RK_c(iStage)
#ifdef PARTICLES
  CALL CountPartsPerElem(ResetNumberOfParticles=.FALSE.) !for scaling of tParts of LB
  IF (iStage.NE.nRKStages) THEN
    RKdtFrac = RK_c(iStage+1)-RK_c(iStage)
    RKdtFracTotal=RKdtFracTotal+RKdtFrac
  ELSE
    RKdtFrac = 1.-RK_c(nRKStages)
    RKdtFracTotal=1.
  END IF
  dtWeight = dt/dt_Min(DT_MIN) * RKdtFrac

  ! deposition
  IF (time.GE.DelayTime) THEN
    CALL Deposition() ! because of emission and UpdateParticlePosition
  END IF
#endif /*PARTICLES*/

  CALL HDG(tStage,U,iter)

#ifdef PARTICLES
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! set last data already here, since surfaceflux moved before interpolation
  LastPartPos(1:3,1:PDM%ParticleVecLength)=PartState(1:3,1:PDM%ParticleVecLength)
  PEM%LastGlobalElemID(1:PDM%ParticleVecLength)=PEM%GlobalElemID(1:PDM%ParticleVecLength)
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  IF (time.GE.DelayTime) THEN
    IF (DoSurfaceFlux)THEN
      CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
    END IF
    ! forces on particle
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    CALL InterpolateFieldToParticle()   ! forces on particles
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
    IF(DoInterpolation) CALL CalcPartRHS()
    ! particle step
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        ! If coupled power output is active and particle carries charge, determine its kinetic energy and store in EDiff
        IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'before')
        IF (.NOT.PDM%IsNewPart(iPart)) THEN
          Pt_temp(1:3,iPart) = PartState(4:6,iPart) - RK_a(iStage) * Pt_temp(1:3,iPart)
          PartState(1:3,iPart) = PartState(1:3,iPart) + Pt_temp(1:3,iPart)*b_dt(iStage)
          ! Don't push the velocity component of neutral particles!
          IF(isPushParticle(iPart))THEN
            Pt_temp(4:6,iPart) = Pt(1:3,iPart) - RK_a(iStage) * Pt_temp(4:6,iPart)
            PartState(4:6,iPart) = PartState(4:6,iPart) + Pt_temp(4:6,iPart)*b_dt(iStage)
          END IF
        ELSE !IsNewPart: no Pt_temp history available!
          IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage
            CALL RANDOM_NUMBER(RandVal)
            PDM%dtFracPush(iPart) = .FALSE.
          ELSE !new but without SF in current RKStage (i.e., from ParticleInserting or diffusive wall reflection)
               ! -> rebuild Pt_tmp-coefficients assuming F=const. (value at last Pos) in previous stages
            RandVal=1. !"normal" particles (i.e. not from SurfFlux) are pushed with whole timestep!
          END IF
          IF (DoForceFreeSurfaceFlux) Pt(1:3,iPart)=0.
          Pa_rebuilt(:,:)=0.
          ! Don't push the velocity component of neutral particles!
          IF(isPushParticle(iPart))THEN
            DO iStage_loc=1,iStage
              Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(1:3,iPart)
            END DO
          END IF
          v_rebuilt(:,:)=0.
          DO iStage_loc=iStage-1,0,-1
            IF (iStage_loc.EQ.iStage-1) THEN
              v_rebuilt(1:3,iStage_loc) = PartState(4:6,iPart) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
            ELSE
              v_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc+1) - b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
            END IF
          END DO
          Pv_rebuilt(:,:)=0.
          DO iStage_loc=1,iStage
            IF (iStage_loc.EQ.1) THEN
              Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,0)
            ELSE
              Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc-1) - RK_a(iStage_loc)*Pv_rebuilt(1:3,iStage_loc-1)
            END IF
          END DO
          Pt_temp(1:3,iPart) = Pv_rebuilt(1:3,iStage)
          PartState(1:3,iPart) = PartState(1:3,iPart) + Pt_temp(1:3,iPart)*b_dt(iStage)*RandVal
          ! Don't push the velocity component of neutral particles!
          IF(isPushParticle(iPart))THEN
            Pt_temp(4:6,iPart) = Pa_rebuilt(1:3,iStage)
            PartState(4:6,iPart) = PartState(4:6,iPart) + Pt_temp(4:6,iPart)*b_dt(iStage)*RandVal
          END IF
          IF (.NOT.DoForceFreeSurfaceFlux .OR. iStage.EQ.nRKStages) PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
        END IF !IsNewPart

        CALL CalcPartSymmetryPos(PartState(1:3,iPart),PartState(4:6,iPart))
      
        ! If coupled power output is active and particle carries charge, calculate energy difference and add to output variable
        IF (CalcCoupledPower) CALL CalcCoupledPowerPart(iPart,'after')
      END IF
    END DO
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

    ! particle tracking
#if USE_MPI
    CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
    CALL PerformTracking()
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    CALL ParticleInserting()
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_EMISSION,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
    CALL SendNbOfParticles() ! send number of particles
    CALL MPIParticleSend()   ! finish communication of number of particles and send particles
    CALL MPIParticleRecv()   ! finish communication
#endif
  END IF
#endif /*PARTICLES*/
END DO

#ifdef PARTICLES
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

END SUBROUTINE TimeStepPoissonByLSERK

END MODULE MOD_TimeStep
#endif /*(PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)*/
#endif /*USE_HDG*/
