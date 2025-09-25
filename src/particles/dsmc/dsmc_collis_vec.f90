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

MODULE MOD_DSMC_CollisVec
!===================================================================================================================================
! Module for chemical reactions including calculation of probabilities and collisions
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DiceDeflectedVelocityVector4Coll
  MODULE PROCEDURE DiceDeflectedVelocityVector4Coll
END INTERFACE

INTERFACE DiceVelocityVector4Coll
  MODULE PROCEDURE DiceVelocityVector4Coll
END INTERFACE

ABSTRACT INTERFACE
  !FUNCTION PostCollisionVeloVec(iPair)
  FUNCTION PostCollisionVeloVec(iPair,ForceUnitVector) RESULT(VeloVec)
    INTEGER,INTENT(IN)          :: iPair               ! index of collision pair
    LOGICAL,INTENT(IN),OPTIONAL :: ForceUnitVector
    REAL :: VeloVec(3)
  END FUNCTION
END INTERFACE

PROCEDURE(PostCollisionVeloVec),POINTER :: PostCollVec    !< pointer defining the function called for random vector after collision
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DiceDeflectedVelocityVector4Coll
PUBLIC :: DiceVelocityVector4Coll
PUBLIC :: PostCollVec
PUBLIC :: VelocityCOMBackscatter
!===================================================================================================================================

CONTAINS


FUNCTION DiceDeflectedVelocityVector4Coll(iPair,ForceUnitVector) RESULT(VeloVec)
!===================================================================================================================================
!> Calculation of post collision velocity vector
!>
!> Calculates deflection angle and resulting deflection relative velocity vector including the coordinate transformation
!> from the reduced mass system back to the COM frame - see Bird 1994 p.36
!> VHS: isotropic    scattering vector for alphaVSS = 1
!> VSS: anisotropic  scattering vector     alphaVSS e [1,2] see collision parameters in dsmc_init for sources
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals       ,ONLY: VECNORM
USE MOD_Globals_Vars  ,ONLY: Pi
USE MOD_DSMC_Vars     ,ONLY: Coll_pData, CollInf
USE MOD_Part_Tools    ,ONLY: DiceUnitVector
USE MOD_Particle_Vars ,ONLY: PartSpecies, PartState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iPair
LOGICAL,INTENT(IN),OPTIONAL :: ForceUnitVector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                       :: VeloVec(3) ! post-collision relative velocity vector cRela*
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                       :: cRela               !> absolute value of pre-coll relative velocity abs(cRela), Bird1994 (2.3),(2.8)
REAL                       :: VeloRel(3),VeloMag,VeloYZMag
REAL                       :: rRan, rotAngle, cos_scatAngle, sin_scatAngle
REAL,DIMENSION(3,3)        :: trafoMatrix
INTEGER                    :: iPart1, iPart2, iSpec1, iSpec2
LOGICAL                    :: ForceUnitVector_loc
!===================================================================================================================================
IF (PRESENT(ForceUnitVector)) THEN
  ForceUnitVector_loc = ForceUnitVector
ELSE
  ForceUnitVector_loc = .FALSE.
END IF

iPart1 = Coll_pData(iPair)%iPart_p1
iPart2 = Coll_pData(iPair)%iPart_p2
iSpec1 = PartSpecies(iPart1)
iSpec2 = PartSpecies(iPart2)

cRela = SQRT(Coll_pData(iPair)%cRela2)  ! absolute value of post-collision relative velocity
IF (CollInf%alphaVSS(iSpec1,iSpec2).NE.1.0 .AND. .NOT.ForceUnitVector_loc) THEN ! VSS

  CALL RANDOM_NUMBER(rRan) ! rRan = (b / d) ^ 2  : dice impact parameter b to distance d relation in y-direction
                           ! 0                   : frontal collision
                           ! 1                   : brush without change of direction

  cos_scatAngle = 2. * rRan ** ( 1. / CollInf%alphaVSS(iSpec1,iSpec2) ) - 1. ! deflection x-component in collision plane  (chi e [-1,1], away from center)
  sin_scatAngle = SQRT ( 1. - cos_scatAngle ** 2. )   ! deflection y-component in collision plane  (                      -of-mass)

  ! transfer 2D collision vector to 3D space through relation of collision to reference plane
  CALL RANDOM_NUMBER(rRan) ! dice rotation angle between collision and reference plane :  epsilon e [0,2*pi]
  rotAngle = 2. * Pi * rRan

  VeloVec(1) = cRela * cos_scatAngle                 ! x-component in collision plane
  VeloVec(2) = cRela * sin_scatAngle * COS(rotAngle) ! y-component between collision and reference plane
  VeloVec(3) = cRela * sin_scatAngle * SIN(rotAngle) ! z-component between collision and reference plane

  VeloRel(1:3) = PartState(4:6,iPart1) - PartState(4:6,iPart2)

  ! if no radial component: collision plane and laboratory identical-> no transformation
  IF ((VeloRel(2).NE.0.) .AND. (VeloRel(3).NE.0.)) THEN
    VeloMag = VECNORM(VeloRel(1:3))
    VeloYZMag = SQRT(VeloRel(2)**2 + VeloRel(3)**2)
    ! axis transformation from reduced- mass frame back to center-of-mass frame via Bird1994 p.36 (2.22)=A*b MATMUL for performance reasons
    ! initializing rotation matrix
    trafoMatrix(1,1) = VeloRel(1) / VeloMag
    trafoMatrix(1,2) = 0.
    trafoMatrix(1,3) = VeloYZMag / VeloMag
    trafoMatrix(2,1) = VeloRel(2) / VeloMag
    trafoMatrix(2,2) = VeloRel(3) / VeloYZMag
    trafoMatrix(2,3) = -VeloRel(1) * VeloRel(2) / (VeloMag * VeloYZMag)
    trafoMatrix(3,1) = VeloRel(3) / VeloMag
    trafoMatrix(3,2) = -VeloRel(2) / VeloYZMag
    trafoMatrix(3,3) = -VeloRel(1) * VeloRel(3) / (VeloMag * VeloYZMag)
    ! relative post collision velocity transformation from reduced mass to COM frame
    VeloVec(:) = MATMUL(trafoMatrix , VeloVec)
  END IF ! transformation
ELSE
  VeloVec(:) = CRela * DiceUnitVector()
END IF ! VSS

END FUNCTION DiceDeflectedVelocityVector4Coll


FUNCTION DiceVelocityVector4Coll(iPair,ForceUnitVector) RESULT(VeloVec)
!===================================================================================================================================
!> Calculates random unit vector with a shift towards one axis
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_DSMC_Vars  ,ONLY: Coll_pData
USE MOD_Part_Tools ,ONLY: DiceUnitVector
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iPair
LOGICAL,INTENT(IN),OPTIONAL :: ForceUnitVector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: VeloVec(3)
!===================================================================================================================================
VeloVec(:) = SQRT(Coll_pData(iPair)%cRela2) * DiceUnitVector()
RETURN
IF(PRESENT(ForceUnitVector)) THEN
END IF
END FUNCTION DiceVelocityVector4Coll


!===================================================================================================================================
!> Backscattering in the centre of mass plane
!===================================================================================================================================
PPURE FUNCTION VelocityCOMBackscatter(iPair) RESULT(VeloVec)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals         ,ONLY: VECNORM
USE MOD_Particle_Vars   ,ONLY: PartState
USE MOD_DSMC_Vars       ,ONLY: Coll_pData
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: VeloVec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iPart1, iPart2
REAL                :: VeloRel(3),VeloMag,VeloYZMag,trafoMatrix(3,3)
!===================================================================================================================================

iPart1 = Coll_pData(iPair)%iPart_p1
iPart2 = Coll_pData(iPair)%iPart_p2

VeloVec(1) = -SQRT(Coll_pData(iPair)%cRela2)                  ! x-component in collision plane
VeloVec(2) = 0.0
VeloVec(3) = 0.0

VeloRel(1:3) = PartState(4:6,iPart1) - PartState(4:6,iPart2)

! if no radial component: collision plane and laboratory identical-> no transformation
IF ((VeloRel(2).NE.0.) .AND. (VeloRel(3).NE.0.)) THEN
  VeloMag = VECNORM(VeloRel(1:3))
  VeloYZMag = SQRT(VeloRel(2)**2 + VeloRel(3)**2)
  ! axis transformation from reduced- mass frame back to center-of-mass frame via Bird1994 p.36 (2.22)=A*b MATMUL for performance reasons
  ! initializing rotation matrix
  trafoMatrix(1,1) = VeloRel(1) / VeloMag
  trafoMatrix(1,2) = 0.
  trafoMatrix(1,3) = VeloYZMag / VeloMag
  trafoMatrix(2,1) = VeloRel(2) / VeloMag
  trafoMatrix(2,2) = VeloRel(3) / VeloYZMag
  trafoMatrix(2,3) = -VeloRel(1) * VeloRel(2) / (VeloMag * VeloYZMag)
  trafoMatrix(3,1) = VeloRel(3) / VeloMag
  trafoMatrix(3,2) = -VeloRel(2) / VeloYZMag
  trafoMatrix(3,3) = -VeloRel(1) * VeloRel(3) / (VeloMag * VeloYZMag)
  ! relative post collision velocity transformation from reduced mass to COM frame
  VeloVec(:) = MATMUL(trafoMatrix , VeloVec)
END IF ! transformation

END FUNCTION VelocityCOMBackscatter

END MODULE MOD_DSMC_CollisVec
