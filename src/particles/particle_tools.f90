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

MODULE MOD_part_tools
!===================================================================================================================================
! Contains tools for particles
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, GetParticleWeight
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars        ,ONLY: PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem, PartPressureCell
USE MOD_Particle_Vars        ,ONLY: KeepWallParticles, PartState, VarTimeStep, usevMPF
USE MOD_DSMC_Vars            ,ONLY: useDSMC, CollInf
USE MOD_Particle_VarTimeStep ,ONLY: CalcVarTimeStep
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools    ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: counter1,i,n
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF(PDM%maxParticleNumber.EQ.0) RETURN
counter1 = 1
IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell.OR.usevMPF) THEN
  PEM%pNumber(:) = 0
  IF (KeepWallParticles) PEM%wNumber(:) = 0
END IF
n = PDM%ParticleVecLength !PDM%maxParticleNumber
PDM%ParticleVecLength = 0
PDM%insideParticleNumber = 0
IF (doParticleMerge) vMPF_SpecNumElem = 0
IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell.OR.usevMPF) THEN
  DO i=1,n
    IF (.NOT.PDM%ParticleInside(i)) THEN
      IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
      PDM%nextFreePosition(counter1) = i
      counter1 = counter1 + 1
    ELSE
      IF (PEM%pNumber(PEM%Element(i)).EQ.0) THEN
        PEM%pStart(PEM%Element(i)) = i                    ! Start of Linked List for Particles in Elem
      ELSE
        PEM%pNext(PEM%pEnd(PEM%Element(i))) = i ! Next Particle of same Elem (Linked List)
      END IF
      PEM%pEnd(PEM%Element(i)) = i
      PEM%pNumber(PEM%Element(i)) = &                      ! Number of Particles in Element
          PEM%pNumber(PEM%Element(i)) + 1
      IF (VarTimeStep%UseVariableTimeStep) THEN
        VarTimeStep%ParticleTimeStep(i) = CalcVarTimeStep(PartState(i,1),PartState(i,2),PEM%Element(i))
      END IF
      IF (KeepWallParticles) THEN
        IF (PDM%ParticleAtWall(i)) THEN
          PEM%wNumber(PEM%Element(i)) = PEM%wNumber(PEM%Element(i)) + 1
        END IF
      END IF
      PDM%ParticleVecLength = i
      IF(doParticleMerge) vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) = vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) + 1
    END IF
  END DO
ELSE ! no DSMC
  DO i=1,n
    IF (.NOT.PDM%ParticleInside(i)) THEN
      PDM%nextFreePosition(counter1) = i
      counter1 = counter1 + 1
    ELSE
      PDM%ParticleVecLength = i
    END IF
  END DO
ENDIF
PDM%insideParticleNumber = PDM%ParticleVecLength - counter1+1
PDM%CurrentNextFreePosition = 0
DO i = n+1,PDM%maxParticleNumber
  IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
  PDM%nextFreePosition(counter1) = i
  counter1 = counter1 + 1
END DO
PDM%nextFreePosition(counter1:PDM%MaxParticleNumber)=0 ! exists if MaxParticleNumber is reached!!!
IF (counter1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber)=0

#if USE_LOADBALANCE
CALL LBPauseTime(LB_UNFP,tLBStart)
#endif /*USE_LOADBALANCE*/

  RETURN
END SUBROUTINE UpdateNextFreePosition

FUNCTION DiceUnitVector()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals_Vars ,ONLY: Pi
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: DiceUnitVector(3)
REAL                     :: iRan, bVec, aVec
!===================================================================================================================================
CALL RANDOM_NUMBER(iRan)
bVec              = 1. - 2.*iRan
aVec              = SQRT(1. - bVec**2.)
DiceUnitVector(3) = bVec
CALL RANDOM_NUMBER(iRan)
bVec              = Pi *2. * iRan
DiceUnitVector(1) = aVec * COS(bVec)
DiceUnitVector(2) = aVec * SIN(bVec)

END FUNCTION DiceUnitVector


PURE REAL FUNCTION GetParticleWeight(iPart)
!===================================================================================================================================
!> Determines the appropriate particle weighting for the axisymmetric case with radial weighting and the variable time step. For
!> radial weighting, the radial factor is multiplied by the regular weighting factor. If only a variable time step is used, at the
!> moment, the regular weighting factor is not included.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Particle_Vars           ,ONLY: usevMPF, VarTimeStep, PartMPF
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  IF (VarTimeStep%UseVariableTimeStep) THEN
    GetParticleWeight = PartMPF(iPart) * VarTimeStep%ParticleTimeStep(iPart)
  ELSE
    GetParticleWeight = PartMPF(iPart)
  END IF
ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
  GetParticleWeight = VarTimeStep%ParticleTimeStep(iPart)
ELSE
  GetParticleWeight = 1.
END IF

END FUNCTION GetParticleWeight

END MODULE MOD_part_tools
