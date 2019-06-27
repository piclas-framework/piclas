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

FUNCTION DiceDeflectedVector(absCRela,ur,vr,wr,alpha)
!===================================================================================================================================
! Calculates deflection angle and resulting deflection vector after Bird 1994.
! VHS: isotropic scattering vector
! VSS: anisotropic scattering vector
! Subsequent coordinate transformation from independent collision coordinate system to original coordinate system.
! Matrix-wise coordinate transformation A*b=(2.22) since it is faster executed, than running transformation line-wise
!(Bird1994,p.36)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  USE MOD_Globals_Vars,           ONLY : Pi
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN)            :: absCRela ! absolute value of pre-coll relative velocity. Bird (2.3),(2.8) 
  REAL,INTENT(IN)            :: ur,vr,wr ! pre-coll relative velocity CRela=(/ur,vr,wr/)
  REAL,INTENT(IN), OPTIONAL  :: alpha    ! VSS parameter

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL                       :: DiceDeflectedVector(3) ! post-coll relative velocity vector DiceDeflectedVector. Bird (CRelaN)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
 REAL                        :: iRan, eps, cos_chi, sin_chi
 REAL,DIMENSION(3,3)         :: trafoMatrix
!===================================================================================================================================
   ! im Code die alphas hinterlegen #datenbank.
   ! macht es unterschied ob a b oder b a?
   ! oder readin
absCRela=SQRT(absCRela)    
   ! determination of DiceDeflectedVector in the independent coordinate system
   CALL RANDOM_NUMBER(iRan)
   IF((.NOT.PRESENT(alpha)).OR.(alpha.EQ.1)) THEN
      !VHS 
      cos_chi         = 2.*iRan-1.! isotropic scattering angle chi between [-1,1]
      WRITE(*,*) "VHS - default" 
     ! einbauen in particle init und vars 
   ELSEIF (alpha.GT.1) THEN
      WRITE(*,*) "VSS - anisotropic scattering (alpha greater than 1)"
      cos_chi         = 2.*iRan**(1./alpha)-1. ! deflected (anisotrop) scattering angle chi 
   ELSE !Error
         WRITE (*,*) "alpha must not be less than 1."
         !either abort or vhs default.
   END IF
   sin_chi                  = SQRT(1. - cos_chi**2.)
   ! DiceDeflectedVector(x,y,z) order according to Bird 1994, p.36  
   DiceDeflectedVector(1)   = absCRela*cos_chi
   CALL RANDOM_NUMBER(iRan)
   eps                      = 2.*PI*iRan ! azimuthal impact angle epsilon between [0,2*pi]
   DiceDeflectedVector(2)   = absCRela*sin_chi*cos(eps)
   DiceDeflectedVector(3)   = absCRela*sin_chi*sin(eps)

   !Transformation to original coordinate system via Bird (2.22)
   IF ((vr.EQ.0) .AND. (wr.EQ.0)) THEN
      ! In case the impact plane system points into the same direction as the original coordinate system the
      ! DiceDeflectedVector needs no change
      ! Initializing scaling matrix
      trafoMatrix(1,1)=1
      trafoMatrix(1,2)=0
      trafoMatrix(1,3)=0
      trafoMatrix(2,1)=0
      trafoMatrix(2,2)=1
      trafoMatrix(2,3)=0
      trafoMatrix(3,1)=0
      trafoMatrix(3,2)=0
      trafoMatrix(3,3)=1
      ! ANDY - to be solved 
   ELSE
      !Initializing rotation matrix
      trafoMatrix(1,1)=ur/absCRela
      trafoMatrix(1,2)=0
      trafoMatrix(1,3)=sqrt(vr**2+wr**2)/absCRela
      trafoMatrix(2,1)=vr/absCRela
      trafoMatrix(2,2)=wr/sqrt(vr**2+wr**2)
      trafoMatrix(2,3)=-ur*vr/(ur*sqrt(vr**2+wr**2))
      trafoMatrix(3,1)=wr/absCRela
      trafoMatrix(3,2)=-vr/sqrt(vr**2+wr**2)
      trafoMatrix(3,3)=-ur*vr/(ur*sqrt(vr**2+wr**2))
   END IF
   ! Transformation and scaling
   DiceDeflectedVector(:)=MATMUL(trafoMatrix,DiceDeflectedVector)
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
