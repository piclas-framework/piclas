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
USE MOD_DSMC_Vars, ONLY : useDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

INTERFACE INSIDEMACROPART
  MODULE PROCEDURE INSIDEMACROPART
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, INSIDEMACROPART
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Particle_Vars, ONLY : PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem, PartPressureCell
  USE MOD_Particle_Vars, ONLY : KeepWallParticles
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                          :: counter1,i,n
!===================================================================================================================================

  IF(PDM%maxParticleNumber.EQ.0) RETURN
  counter1 = 1
  IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell) THEN
    PEM%pNumber(:) = 0
    IF (KeepWallParticles) PEM%wNumber(:) = 0
  END IF
  n = PDM%ParticleVecLength !PDM%maxParticleNumber
  PDM%ParticleVecLength = 0
  PDM%insideParticleNumber = 0
  IF (doParticleMerge) vMPF_SpecNumElem = 0
  IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell) THEN
   DO i=1,n
     IF (.NOT.PDM%ParticleInside(i)) THEN
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
   PDM%nextFreePosition(counter1) = i
   counter1 = counter1 + 1
  END DO 
  PDM%nextFreePosition(counter1:PDM%MaxParticleNumber)=0 ! exists if MaxParticleNumber is reached!!!
  IF (counter1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber)=0

  RETURN
END SUBROUTINE UpdateNextFreePosition

FUNCTION DiceUnitVector()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  USE MOD_Globals_Vars,           ONLY : Pi
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


LOGICAL FUNCTION INSIDEMACROPART(Particle_pos)
!===================================================================================================================================
!> Function for checking if particle position would be inside of any macro-particle in the local domain
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars ,ONLY: MacroPart, nMacroParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: refPosSphere(1:3), distance
INTEGER :: iMP
!===================================================================================================================================
INSIDEMACROPART = .FALSE.
DO iMP=1,nMacroParticle
  !IF (ElemHasMacroPart(ElemID,:))) THEN
  refPosSphere(1:3) = MacroPart(iMP)%center(1:3)
  distance=SQRT(DOT_PRODUCT((Particle_pos-refPosSphere),(Particle_pos-refPosSphere)))
  IF (distance.LE.MacroPart(iMP)%radius) THEN
    INSIDEMACROPART=.TRUE.
    RETURN
  END IF
END DO
END FUNCTION INSIDEMACROPART


END MODULE MOD_part_tools
