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
USE MOD_Globals                  ,ONLY: abort
USE MOD_Particle_Vars            ,ONLY: PartState, LastPartPos, PDM, PEM, DoSurfaceFlux, WriteMacroVolumeValues
USE MOD_Particle_Vars            ,ONLY: WriteMacroSurfaceValues, Symmetry, VarTimeStep, Species, PartSpecies
USE MOD_Particle_Vars            ,ONLY: UseSplitAndMerge
USE MOD_DSMC_Vars                ,ONLY: DSMC_RHS, DSMC, CollisMode, AmbipolElecVelo
USE MOD_DSMC                     ,ONLY: DSMC_main
USE MOD_part_tools               ,ONLY: UpdateNextFreePosition
USE MOD_part_emission            ,ONLY: ParticleInserting
USE MOD_Particle_SurfFlux        ,ONLY: ParticleSurfaceflux
USE MOD_Particle_SurfChemFlux
USE MOD_Particle_Tracking_vars   ,ONLY: tTracking,MeasureTrackTime
USE MOD_Particle_Tracking        ,ONLY: PerformTracking
USE MOD_SurfaceModel_Porous      ,ONLY: PorousBoundaryRemovalProb_Pressure
USE MOD_SurfaceModel_Vars        ,ONLY: nPorousBC, DoChemSurface
USE MOD_Particle_Boundary_Vars   ,ONLY: PartBound
USE MOD_vMPF                     ,ONLY: SplitAndMerge
#if USE_MPI
USE MOD_Particle_MPI             ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Boundary_Sampling, ONLY: ExchangeChemSurfData
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                       :: timeEnd, timeStart, dtVar, RandVal, NewYPart, NewYVelo
INTEGER                    :: iPart
#if USE_LOADBALANCE
REAL                  :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

  IF (DoSurfaceFlux) THEN
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURF,tLBStart)
#endif /*USE_LOADBALANCE*/

    CALL ParticleSurfaceflux()
  END IF

  IF (DoChemSurface) THEN
    CALL ParticleSurfChemFlux()
    CALL ExchangeChemSurfData()
    ! Write the adsorption and desorption values into separate files for the two reactants
    OPEN(10, file='adsorbtion_O2.txt', position="APPEND")
    OPEN(30, file='adsorbtion_CO.txt', position="APPEND")
    OPEN(20, file='desorption_O2.txt', position="APPEND")
    OPEN(40, file='desorption_CO.txt', position="APPEND")
    OPEN(50, file='reaction_LH.txt', position="APPEND")
    OPEN(60, file='reaction_ER.txt', position="APPEND")

      WRITE(10,*) PartBound%AdCount(1,4)
      WRITE(20,*) PartBound%DesCount(1,4)
      WRITE(30,*) PartBound%AdCount(1,2)
      WRITE(40,*) PartBound%DesCount(1,2)
      WRITE(50,*) PartBound%LHCount(1,1)
      WRITE(60,*) PartBound%ERCount(1,1)

    CLOSE(10)
    CLOSE(20)
    CLOSE(30)
    CLOSE(40)
    CLOSE(50)
    CLOSE(60)
  END IF

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
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
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
          NewYVelo = (AmbipolElecVelo(iPart)%ElecVelo(2)*(PartState(2,iPart))+AmbipolElecVelo(iPart)%ElecVelo(3)*PartState(3,iPart))/NewYPart
          AmbipolElecVelo(iPart)%ElecVelo(3)= (-AmbipolElecVelo(iPart)%ElecVelo(2)*PartState(3,iPart) &
            + AmbipolElecVelo(iPart)%ElecVelo(3)*(PartState(2,iPart)))/NewYPart         
          AmbipolElecVelo(iPart)%ElecVelo(2) = NewYVelo
        END IF
      END IF
      NewYVelo = (PartState(5,iPart)*(PartState(2,iPart))+PartState(6,iPart)*PartState(3,iPart))/NewYPart
      PartState(6,iPart) = (-PartState(5,iPart)*PartState(3,iPart)+PartState(6,iPart)*(PartState(2,iPart)))/NewYPart
      PartState(2,iPart) = NewYPart
      PartState(3,iPart) = 0.0
      PartState(5,iPart) = NewYVelo
      END IF
    END IF
  END DO
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

  ! Resetting the particle positions in the second/third dimension for the 1D/2D/axisymmetric case
  IF(Symmetry%Order.LT.3) THEN
    LastPartPos(Symmetry%Order+1:3,1:PDM%ParticleVecLength) = 0.0
    PartState(Symmetry%Order+1:3,1:PDM%ParticleVecLength) = 0.0
  END IF

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
  ELSE IF (PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).GT.PDM%maxParticleNumber .OR. &
           PDM%nextFreePosition(PDM%CurrentNextFreePosition+1).EQ.0) THEN
    CALL abort(&
    __STAMP__&
    ,'maximum nbr of particles reached!')  !gaps in PartState are not filled until next UNFP and array might overflow more easily!
  END IF

  CALL DSMC_main()

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  PartState(4:6,1:PDM%ParticleVecLength) = PartState(4:6,1:PDM%ParticleVecLength) + DSMC_RHS(1:3,1:PDM%ParticleVecLength)

  IF(UseSplitAndMerge) CALL SplitAndMerge()

#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE TimeStep_DSMC


END MODULE MOD_TimeStep
#endif
