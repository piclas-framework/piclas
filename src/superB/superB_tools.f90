!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_SuperB_Tools
!===================================================================================================================================
!> Contains the routines and algorithms required for the calculation of magnetic fields with SuperB
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: FindLinIndependentVectors, GramSchmidtAlgo
!===================================================================================================================================
INTERFACE FindLinIndependentVectors
  MODULE PROCEDURE FindLinIndependentVectors
END INTERFACE FindLinIndependentVectors

INTERFACE GramSchmidtAlgo
  MODULE PROCEDURE GramSchmidtAlgo
END INTERFACE GramSchmidtAlgo
!===================================================================================================================================

CONTAINS

SUBROUTINE FindLinIndependentVectors(NormalVector, Vector1, Vector2)
!===================================================================================================================================
! Finds two linear vectors of a normal vector around a base point
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: NormalVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: Vector1(3), Vector2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Find the second vector which is in the normal plane
IF (NormalVector(1).NE.0) THEN
  Vector1(1) = (0 - NormalVector(2) - NormalVector(3)) / NormalVector(1)
  Vector1(2) = 1
  Vector1(3) = 1
ELSE IF (NormalVector(2).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = (0 - NormalVector(1) - NormalVector(3)) / NormalVector(2)
  Vector1(3) = 1
ELSE IF (NormalVector(3).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = 1
  Vector1(3) = (0 - NormalVector(1) - NormalVector(2)) / NormalVector(3)
ELSE
  CALL abort(__STAMP__&
      ,'The coil direction vector can not be (0,0,)')
END IF

! Find the third vecord vector with the cross product
Vector2(1) = NormalVector(2)*Vector1(3) - NormalVector(3)*Vector1(2)
Vector2(2) = NormalVector(3)*Vector1(1) - NormalVector(1)*Vector1(3)
Vector2(3) = NormalVector(1)*Vector1(2) - NormalVector(2)*Vector1(1)

END SUBROUTINE FindLinIndependentVectors


SUBROUTINE GramSchmidtAlgo(Vector1, Vector2, Vector3)
!===================================================================================================================================
! Contains the Gram Schmidt algorithm for an orthonormal basis
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT) :: Vector1(3), Vector2(3), Vector3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! v1 = w1/||w1||
Vector1(:) = Vector1(:) / SQRT(Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2)

! v2 = w2 - <v1,w2>*v1
Vector2(:) = Vector2(:) - DOT_PRODUCT(Vector1, Vector2) * Vector1(:)
! v2 = v2/||v2||
Vector2(:) = Vector2(:) / SQRT(Vector2(1)**2 + Vector2(2)**2 + Vector2(3)**2)

! v3 = w3 - <v1,w3>*v1 - <v2,w3>*v2
Vector3(:) = Vector3(:) - DOT_PRODUCT(Vector1, Vector3) * Vector1(:) -&
                          DOT_PRODUCT(Vector2,Vector3) *  Vector3(:)
! v3 = v3/||v3||
Vector3(:) = Vector3(:) / SQRT(Vector3(1)**2 + Vector3(2)**2 + Vector3(3)**2)

END SUBROUTINE GramSchmidtAlgo

END MODULE MOD_SuperB_Tools