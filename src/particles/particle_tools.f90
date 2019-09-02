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

INTERFACE DiceDeflectedVector
  MODULE PROCEDURE DiceDeflectedVector
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, GetParticleWeight, DiceDeflectedVector
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Particle_Vars,          ONLY: PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem, PartPressureCell
  USE MOD_Particle_Vars,          ONLY: KeepWallParticles, PartState, VarTimeStep
  USE MOD_DSMC_Vars,              ONLY: useDSMC, CollInf
  USE MOD_Particle_VarTimeStep,   ONLY: CalcVarTimeStep
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
       IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
       PDM%nextFreePosition(counter1) = i
       counter1 = counter1 +  1
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

  RETURN
END SUBROUTINE UpdateNextFreePosition

FUNCTION DiceDeflectedVector(CRela2,ur,vr,wr,alpha)
!===================================================================================================================================
! Calculates deflection angle and resulting deflection vector after Bird 1994.
! VHS: isotropic scattering vector
! VSS: anisotropic scattering vector
! Subsequent coordinate transformation from independent collision coordinate system to original coordinate system.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  USE MOD_Globals_Vars,           ONLY : Pi
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN)            :: CRela2                 ! instead of reading in the whole container  
  REAL,INTENT(IN)            :: ur,vr,wr               ! pre-coll relative velocity CRela=(/ur,vr,wr/)
  REAL,INTENT(IN), OPTIONAL  :: alpha                  ! VSS parameter

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL                       :: DiceDeflectedVector(3) ! post-coll relative velocity vector DiceDeflectedVector. Bird (CRelaN)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES   
 REAL                        :: CRela                  ! absolute value of pre-coll relative velocity. Bird (2.3),(2.8) 
 REAL                        :: iRan, eps, cos_chi, sin_chi
 REAL,DIMENSION(3,3)         :: trafoMatrix
!===================================================================================================================================
  ! determination of DiceDeflectedVector in independent coordinate system and scaling
  CRela=SQRT(CRela2)
  CALL RANDOM_NUMBER(iRan)
  cos_chi                  = 2.*iRan**(1./alpha)-1.     ! deflected (anisotrop) scattering angle chi 
                                                        ! if alpha=1 VHS isotropic scattering angle chi between [-1,1]
  sin_chi                  = SQRT(1. - cos_chi**2.)
  DiceDeflectedVector(1)   = CRela*cos_chi              ! DiceDeflectedVector(x,y,z) order according to Bird 1994, p.36  
  CALL RANDOM_NUMBER(iRan)
  eps                      = 2.*PI*iRan                 ! azimuthal impact angle epsilon between [0,2*pi]
  DiceDeflectedVector(2)   = CRela*sin_chi*cos(eps)
  DiceDeflectedVector(3)   = CRela*sin_chi*sin(eps)
  IF (alpha.GT.1) THEN ! VSS
    IF ((vr.EQ.0) .AND. (wr.EQ.0)) THEN
      ! In case the impact plane system points into the same direction as the
      ! original coordinate system the DiceDeflectedVector needs no change.
    ELSE   
      ! Transformation to original coordinate system via Bird1994 p.36
      ! Matrix-wise coordinate transformation A*b=(2.22) since it is faster
      ! Initializing rotation matrix
      trafoMatrix(1,1)=ur/CRela
      trafoMatrix(1,2)=0
      trafoMatrix(1,3)=sqrt(vr**2+wr**2)/CRela
      trafoMatrix(2,1)=vr/CRela
      trafoMatrix(2,2)=wr/sqrt(vr**2+wr**2)
      trafoMatrix(2,3)=-ur*vr/(CRela*sqrt(vr**2+wr**2))
      trafoMatrix(3,1)=wr/CRela
      trafoMatrix(3,2)=-vr/sqrt(vr**2+wr**2)
      trafoMatrix(3,3)=-ur*wr/(CRela*sqrt(vr**2+wr**2))
      ! Transformation 
      DiceDeflectedVector(:)=MATMUL(trafoMatrix,DiceDeflectedVector)
    END IF
  END IF
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
