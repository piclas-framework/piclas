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

INTERFACE DiceDeflectedVector
  MODULE PROCEDURE DiceDeflectedVector
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector,DiceDeflectedVector
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


FUNCTION DiceDeflectedVector(alpha)
!===================================================================================================================================
! Calculates deflection angle and resulting deflection vector after bird1994.
! Subsequent coordinate transformation from independent coordinate system to Center of Mass system of colliding particles.
! Matrix-wise coordinate transformation A*b=(2.22) since it is faster executed, than running transformation line-wise
!(bird1994,p.36)
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
  REAL                       :: DiceDeflectedVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
 REAL                   :: u_r_new,v_r_new,w_r_new
 REAL                   ::  cr ! absolute value of pre-collision relative velocity CRela=c(ipart)-c(jpart)
 ! pre-collision relative velocity is scaled, to have the magnitude of post-collision relative speed
 REAL                   :: iRan, eps, cos_chi, sin_chi
 INTEGER,DIMENSION(3,3) :: trafoMatrix
 INTEGER,DIMENSION(3)   :: deflectedVector ! post-collision relative velocity vector in an independent vector system

!===================================================================================================================================
   ! im Code die alphas hinterlegen #datenbank.
   ! macht es unterschied ob a b oder b a?
   ! oder readin
   CALL RANDOM_NUMBER(iRan)
   IF((.NOT.PRESENT(alpha)).OR.(alpha.EQ.1)) THEN
      !VHS 
      cos_chi         = 2.*iRan-1.! isotropic scattering angle chi between [-1,1]
      WRITE(*,*) "VHS - default" 
     ! einbauen in particle init und vars 
   ELSEIF (alpha.GT.1) THEN
      WRITE(*,*) "VSS - alpha greater than 1"
      cos_chi         = 2.*iRan**(1./alpha)-1. ! deflected (anisotrop) scattering angle chi 
   ELSE !Error
         WRITE (*,*) "alpha must not be less than 1."
         !either abort or vhs default.
   END IF
   sin_chi                  = SQRT(1. - cos_chi**2.) 
   DiceDeflectedVector(3)   = cos_chi
   CALL RANDOM_NUMBER(iRan)
   eps                      = 2.*PI*iRan ! azimuthal impact angle epsilon between [0,2*pi]
   DiceDeflectedVector(1)   = sin_chi*cos(eps)
   DiceDeflectedVector(2)   = sin_chi*sin(eps)



!Initialization rotation matrix
trafoMatrix(1,1)=u_r_new/cr
trafoMatrix(1,2)=0
trafoMatrix(1,3)=sqrt(v_r_new**2+w_r_new**2)/cr
trafoMatrix(2,1)=v_r_new/cr
trafoMatrix(2,2)=w_r_new/sqrt(v_r_new**2+w_r_new**2)
trafoMatrix(2,3)=-u_r_new*v_r_new/(u_r_new*sqrt(v_r_new**2+w_r_new**2))
trafoMatrix(3,1)=w_r_new/cr
trafoMatrix(3,2)=-v_r_new/sqrt(v_r_new**2+w_r_new**2)
trafoMatrix(3,3)=-u_r_new*v_r_new/(u_r_new*sqrt(v_r_new**2+w_r_new**2))


! Diced deflected vector in the independent coordinate system
deflectedVector(1)=cr*cos_chi
deflectedVector(2)=cr*sin_chi*cos(eps)
deflectedVector(3)=cr*sin_chi*sin(eps)

!DiceDeflectedVetor (=CRelaN,bird) in the original coordinate system
DiceDeflectedVector(:)=MATMUL(trafoMatrix,deflectedVector)
!write(*,*) trafoMatrix !debug

END FUNCTION DiceDeflectedVector

FUNCTION DiceUnitVector()
!===================================================================================================================================
! Calculates random unit vector, which expresses isotropic scattering as featured in VHS.
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
  REAL                     :: iRan,cos_chi,sin_chi,eps
!===================================================================================================================================
  CALL RANDOM_NUMBER(iRan)
  cos_chi           = 2.*iRan-1. ! z random value between [-1,1]
  sin_chi           = SQRT(1. - cos_chi**2.) 
  DiceUnitVector(3) = cos_chi
  CALL RANDOM_NUMBER(iRan)
  eps               = 2.*Pi* iRan ! phi random value between [0,2*pi]
  DiceUnitVector(1) = sin_chi * COS(eps)
  DiceUnitVector(2) = sin_chi * SIN(eps)

END FUNCTION DiceUnitVector 

END MODULE MOD_part_tools
