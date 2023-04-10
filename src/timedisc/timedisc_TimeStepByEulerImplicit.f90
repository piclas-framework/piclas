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


#if (PP_TimeDiscMethod==100)
MODULE MOD_TimeStep
!===================================================================================================================================
! Module for the Temporal discretization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: TimeStepByEulerImplicit
!===================================================================================================================================

CONTAINS


SUBROUTINE TimeStepByEulerImplicit()
!===================================================================================================================================
! Euler Implicit method:
! U^n+1 = U^n + dt*R(U^n+1)
! (I -dt*R)*U^n+1 = U^n
! Solve Linear System
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars               ,ONLY: U,Ut
USE MOD_DG                    ,ONLY: DGTimeDerivative_weakForm
USE MOD_TimeDisc_Vars         ,ONLY: dt,time
USE MOD_LinearSolver          ,ONLY: LinearSolver
USE MOD_LinearOperator        ,ONLY: EvalResidual
#ifdef PARTICLES
USE MOD_PICDepo               ,ONLY: Deposition                                      ! , DepositionMPF
USE MOD_PICInterpolation      ,ONLY: InterpolateFieldToParticle
USE MOD_PIC_Vars              ,ONLY: PIC
USE MOD_Particle_Vars         ,ONLY: PartState, Pt, LastPartPos, DelayTime, PEM, PDM ! , usevMPF
USE MOD_part_RHS              ,ONLY: CalcPartRHS
USE MOD_PICInterpolation_Vars ,ONLY: DoInterpolation
USE MOD_part_emission         ,ONLY: ParticleInserting
USE MOD_DSMC                  ,ONLY: DSMC_main
USE MOD_DSMC_Vars             ,ONLY: useDSMC, DSMC
USE MOD_part_tools            ,ONLY: UpdateNextFreePosition
#endif
USE MOD_Particle_Tracking     ,ONLY: PerformTracking
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: tstage,coeff,Norm_R0
!===================================================================================================================================

! one Euler implicit step
! time for source is time + dt
tstage = time + dt
coeff  = dt*1.

#ifdef maxwell
IF(PrecondType.GT.0)THEN
  IF (iter==0) CALL BuildPrecond(time,time,0,1.,dt)
END IF
#endif /*maxwell*/

IF (time.GE.DelayTime) CALL ParticleInserting()

IF ((time.GE.DelayTime).OR.(time.EQ.0)) CALL Deposition()

IF (time.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle()
  IF(DoInterpolation) CALL CalcPartRHS()
END IF
! particles
LastPartPos(1:3,1:PDM%ParticleVecLength)=PartState(1:3,1:PDM%ParticleVecLength)
PEM%LastGlobalElemID(1:PDM%ParticleVecLength)=PEM%GlobalElemID(1:PDM%ParticleVecLength)
IF (time.GE.DelayTime) THEN ! Euler-Explicit only for Particles
  PartState(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength) + dt * PartState(4:6,1:PDM%ParticleVecLength)
  PartState(4:6,1:PDM%ParticleVecLength) = PartState(4:6,1:PDM%ParticleVecLength) + dt * Pt(1:3,1:PDM%ParticleVecLength)
END IF


#if USE_MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
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


! EM field
! U predict
!U = U
! b
LinSolverRHS = U
ImplicitSource=0.
CALL EvalResidual(time,Coeff,Norm_R0)
CALL LinearSolver(tstage,coeff,Norm_R0=Norm_R0)
CALL DivCleaningDamping()
CALL UpdateNextFreePosition()
IF (useDSMC) THEN
  CALL DSMC_main()
END IF

END SUBROUTINE TimeStepByEulerImplicit


END MODULE MOD_TimeStep
#endif /*PP_TimeDiscMethod==100*/
