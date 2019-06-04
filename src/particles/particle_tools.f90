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

INTERFACE diceCollVector
  MODULE PROCEDURE diceCollVector
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector,diceCollVector
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


SUBROUTINE diceCollVector(diceVector,alpha)
!===================================================================================================================================
! Calculates deflection angle and resulting deflection vector after bird1994/Hawk. Variables of functions must not be optional.
! Therefore, this is a subroutine.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  USE MOD_Globals_Vars,           ONLY : Pi
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN), OPTIONAL  :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL,INTENT(OUT)           :: diceVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                       :: iRan, eps, cos_chi, sin_chi,alpha_aux
!===================================================================================================================================
! im Code die alphas hinterlegen #datenbank.
! macht es unterschied ob a b oder b a?
IF(.NOT.PRESENT(alpha)) THEN
   !VHS isotropic scattering angle with random value between [-1,1]
   alpha_aux=1.
   WRITE(*,*) "VHS" 
   cos_chi         = 2.*iRan**(1./alpha_aux)-1. !chi - scattering angle   
ELSE
   
   IF (alpha.GE.1) THEN
      cos_chi         = 2.*iRan**(1./alpha)-1. !chi - scattering angle   
      IF (alpha.GT.1) WRITE(*,*) "VSS"
      IF (alpha.EQ.1) WRITE(*,*) "VHS" 
   ELSE !Error
      WRITE (*,*) "Alpha must not be less than 1"
   END IF
END IF
CALL RANDOM_NUMBER(iRan)
sin_chi         = SQRT(1. - cos_chi**2.) 
diceVector(3)   = cos_chi
CALL RANDOM_NUMBER(iRan)
eps             = 2.*PI*iRan ! azimuthal impact angle epsilon
diceVector(1)   = sin_chi*cos(eps)
diceVector(2)   = sin_chi*sin(eps)
! pre-collision relative velocity is scaled, to have the magnitude of post-collision relative speed
END SUBROUTINE diceCollVector

FUNCTION DiceUnitVector()
!===================================================================================================================================
! Calculates random unit vector.
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
  REAL                     :: iRan,z,prefactor,phi
!===================================================================================================================================
  CALL RANDOM_NUMBER(iRan)
  z                 = 1. - 2.*iRan ! z random value between [-1,1]
  prefactor         = SQRT(1. - z**2.) 
  DiceUnitVector(3) = z
  CALL RANDOM_NUMBER(iRan)
  phi               = Pi *2. * iRan ! phi random value between [0,2*pi]
  DiceUnitVector(1) = prefactor * COS(phi)
  DiceUnitVector(2) = prefactor * SIN(phi)

END FUNCTION DiceUnitVector 

END MODULE MOD_part_tools
