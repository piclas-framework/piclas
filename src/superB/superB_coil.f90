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

MODULE MOD_SuperB_Coil
!===================================================================================================================================
! Contains the calculation of the magnetic field of coils
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: SetUpCoil, SetUpCircleCoil, SetUpRectangleCoil, SetUpLinearConductor, FindLinIndependendVectors, GramSchmidtAlgo, &
          BiotSavart, Jefimenko, WriteCoilVTK, FinalizeCoil
!===================================================================================================================================

INTERFACE SetUpCoil
  MODULE PROCEDURE SetUpCoil
END INTERFACE SetUpCoil

INTERFACE SetUpCircleCoil
  MODULE PROCEDURE SetUpCircleCoil
END INTERFACE SetUpCircleCoil

INTERFACE SetUpRectangleCoil
  MODULE PROCEDURE SetUpRectangleCoil
END INTERFACE SetUpRectangleCoil

INTERFACE SetUpLinearConductor
  MODULE PROCEDURE SetUpLinearConductor
END INTERFACE SetUpLinearConductor

INTERFACE FindLinIndependendVectors
  MODULE PROCEDURE FindLinIndependendVectors
END INTERFACE FindLinIndependendVectors

INTERFACE GramSchmidtAlgo
  MODULE PROCEDURE GramSchmidtAlgo
END INTERFACE GramSchmidtAlgo

INTERFACE BiotSavart
  MODULE PROCEDURE BiotSavart
END INTERFACE BiotSavart

INTERFACE Jefimenko
  MODULE PROCEDURE Jefimenko
END INTERFACE Jefimenko

INTERFACE WriteCoilVTK
  MODULE PROCEDURE WriteCoilVTK
END INTERFACE WriteCoilVTK

INTERFACE FinalizeCoil
  MODULE PROCEDURE FinalizeCoil
END INTERFACE FinalizeCoil
!===================================================================================================================================

CONTAINS

SUBROUTINE SetUpCoil(iCoil)
!===================================================================================================================================
! Sets up a coil defined by segments
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Coil_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: AxisVec2(3), LastSegPoint(3)
REAL                :: TrafoMatrix(3,3)
INTEGER             :: iLoop, iSegment, iPoint
INTEGER             :: PointNumber
REAL                :: LengthAxisVec1, LengthAxisVec2, Phi1, Phi2, Radius
!===================================================================================================================================

! Calculate the second axis vector with the cross product
AxisVec2(1) = CoilInfo(iCoil)%LengthVector(2) * CoilInfo(iCoil)%AxisVec1(3) -&
              CoilInfo(iCoil)%LengthVector(3) * CoilInfo(iCoil)%AxisVec1(2)
AxisVec2(2) = CoilInfo(iCoil)%LengthVector(3) * CoilInfo(iCoil)%AxisVec1(1) -&
              CoilInfo(iCoil)%LengthVector(1) * CoilInfo(iCoil)%AxisVec1(3)
AxisVec2(3) = CoilInfo(iCoil)%LengthVector(1) * CoilInfo(iCoil)%AxisVec1(2) -&
              CoilInfo(iCoil)%LengthVector(2) * CoilInfo(iCoil)%AxisVec1(1)
! Calculate the length of the axis vectors in order to get the normalized base vectors
LengthAxisVec1 = SQRT(CoilInfo(iCoil)%AxisVec1(1)**2 + CoilInfo(iCoil)%AxisVec1(2)**2 +&
                      CoilInfo(iCoil)%AxisVec1(3)**2)
LengthAxisVec2 = SQRT(AxisVec2(1)**2 + AxisVec2(2)**2 + AxisVec2(3)**2)
CoilInfo(iCoil)%LengthVector = CoilInfo(iCoil)%LengthVector / CoilInfo(iCoil)%Length
CoilInfo(iCoil)%AxisVec1     = CoilInfo(iCoil)%AxisVec1 / LengthAxisVec1
AxisVec2                     = AxisVec2 / LengthAxisVec2
! Transformation matrix
TrafoMatrix(:,1) = CoilInfo(iCoil)%LengthVector
TrafoMatrix(:,2) = CoilInfo(iCoil)%AxisVec1
TrafoMatrix(:,3) = AxisVec2

ALLOCATE(CoilNodes(1:3,CoilInfo(iCoil)%NumNodes))
! Start at (0,0,0) in the untransformed coordinate system
CoilNodes(:,1)  = 0
LastSegPoint(:) = 0
PointNumber = 1
DO iLoop=1,CoilInfo(iCoil)%LoopNum
  DO iSegment=1,CoilInfo(iCoil)%NumOfSegments
    DO iPoint=1,CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints-1
      PointNumber = PointNumber+1
      CoilNodes(1,PointNumber) = CoilInfo(iCoil)%Length * (PointNumber - 1) / (CoilInfo(iCoil)%NumNodes - 1)
      IF (CoilInfo(iCoil)%SegmentInfo(iSegment)%SegmentType.EQ.1) THEN
        CoilNodes(2,PointNumber) = CoilInfo(iCoil)%SegmentInfo(iSegment)%LineVector(1) * iPoint /&
                                  (CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1) + LastSegPoint(2)
        CoilNodes(3,PointNumber) = CoilInfo(iCoil)%SegmentInfo(iSegment)%LineVector(2) * iPoint /&
                                  (CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1) + LastSegPoint(3)
      ELSE IF (CoilInfo(iCoil)%SegmentInfo(iSegment)%SegmentType.EQ.2) THEN
        Phi1   = CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi1
        Phi2   = CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi2
        Radius = CoilInfo(iCoil)%SegmentInfo(iSegment)%Radius
        CoilNodes(2,PointNumber) = LastSegPoint(2) - COS(Phi1)*Radius +&
                                   COS((Phi2 - Phi1)/(CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1)*&
                                        iPoint + Phi1)*Radius
        CoilNodes(3,PointNumber) = LastSegPoint(3) - SIN(Phi1)*Radius +&
                                   SIN((Phi2 - Phi1)/(CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1)*&
                                        iPoint + Phi1)*Radius
      ENDIF
    ENDDO !iPoint
    LastSegPoint(:) = CoilNodes(:,PointNumber)
  ENDDO !iSegment
ENDDO !iLoop

! Transform the nodes to the coil coordinate system
DO iPoint=1,CoilInfo(iCoil)%NumNodes
  CoilNodes(:,iPoint) = MATMUL(TrafoMatrix, CoilNodes(:,iPoint))
  CoilNodes(:,iPoint) = CoilNodes(:,iPoint) + CoilInfo(iCoil)%BasePoint(:)
ENDDO

END SUBROUTINE SetUpCoil

SUBROUTINE SetUpCircleCoil(iCoil)
!===================================================================================================================================
! Sets up a coil consisting of a circle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Coil_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: AxisVec1(3), AxisVec2(3)
REAL                :: TrafoMatrix(3,3)
INTEGER             :: iPoint
!===================================================================================================================================

CALL FindLinIndependendVectors(CircleCoilInfo(iCoil)%LengthVector, AxisVec1, AxisVec2)
CALL GramSchmidtAlgo(CircleCoilInfo(iCoil)%LengthVector, AxisVec1, AxisVec2)
TrafoMatrix(:,1) = CircleCoilInfo(iCoil)%LengthVector
TrafoMatrix(:,2) = AxisVec1
TrafoMatrix(:,3) = AxisVec2

ALLOCATE(CoilNodes(1:3,CircleCoilInfo(iCoil)%NumNodes))

DO iPoint=0,CircleCoilInfo(iCoil)%NumNodes - 1
  CoilNodes(1, iPoint+1) = CircleCoilInfo(iCoil)%Length * iPoint / (CircleCoilInfo(iCoil)%NumNodes - 1)
  CoilNodes(2, iPoint+1) = CircleCoilInfo(iCoil)%Radius * COS(2*PI*iPoint/CircleCoilInfo(iCoil)%PointsPerLoop)
  CoilNodes(3, iPoint+1) = CircleCoilInfo(iCoil)%Radius * SIN(2*PI*iPoint/CircleCoilInfo(iCoil)%PointsPerLoop)
  CoilNodes(:, iPoint+1) = MATMUL(TrafoMatrix, CoilNodes(:, iPoint+1))
  CoilNodes(:, iPoint+1) = CoilNodes(:, iPoint+1) + CircleCoilInfo(iCoil)%BasePoint(:)
ENDDO

END SUBROUTINE SetUpCircleCoil

SUBROUTINE SetUpRectangleCoil(iCoil)
!===================================================================================================================================
! Sets up a coil consisting of a rectangle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Coil_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: AxisVec2(3)
REAL                :: LengthAxisVec1, LengthAxisVec2, RectVecLength1, RectVecLength2
INTEGER             :: PointsPerSide, PointNumber, iLoop, iPoint
REAL                :: TrafoMatrix(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------

! Calculate the second axis vector with the cross product
AxisVec2(1) = RectangleCoilInfo(iCoil)%LengthVector(2) * RectangleCoilInfo(iCoil)%AxisVec1(3) -&
              RectangleCoilInfo(iCoil)%LengthVector(3) * RectangleCoilInfo(iCoil)%AxisVec1(2)
AxisVec2(2) = RectangleCoilInfo(iCoil)%LengthVector(3) * RectangleCoilInfo(iCoil)%AxisVec1(1) -&
              RectangleCoilInfo(iCoil)%LengthVector(1) * RectangleCoilInfo(iCoil)%AxisVec1(3)
AxisVec2(3) = RectangleCoilInfo(iCoil)%LengthVector(1) * RectangleCoilInfo(iCoil)%AxisVec1(2) -&
              RectangleCoilInfo(iCoil)%LengthVector(2) * RectangleCoilInfo(iCoil)%AxisVec1(1)              
! Calculate the length of the axis vectors in order to get the normalized base vectors
LengthAxisVec1 = SQRT(RectangleCoilInfo(iCoil)%AxisVec1(1)**2 + RectangleCoilInfo(iCoil)%AxisVec1(2)**2 +&
                      RectangleCoilInfo(iCoil)%AxisVec1(3)**2)
LengthAxisVec2 = SQRT(AxisVec2(1)**2 + AxisVec2(2)**2 + AxisVec2(3)**2)
RectangleCoilInfo(iCoil)%LengthVector = RectangleCoilInfo(iCoil)%LengthVector / RectangleCoilInfo(iCoil)%Length
RectangleCoilInfo(iCoil)%AxisVec1     = RectangleCoilInfo(iCoil)%AxisVec1 / LengthAxisVec1
AxisVec2                              = AxisVec2 / LengthAxisVec2
! Transformation matrix
TrafoMatrix(:,1) = RectangleCoilInfo(iCoil)%LengthVector
TrafoMatrix(:,2) = RectangleCoilInfo(iCoil)%AxisVec1
TrafoMatrix(:,3) = AxisVec2

ALLOCATE(CoilNodes(1:3,RectangleCoilInfo(iCoil)%NumNodes))

PointsPerSide = (RectangleCoilInfo(iCoil)%PointsPerLoop - 1) / 4

! Start at (0,0,0) in the untransformed coordinate system
CoilNodes(1,1) = 0
CoilNodes(2,1) = - RectangleCoilInfo(iCoil)%RectVec1(1)/2 - RectangleCoilInfo(iCoil)%RectVec2(1)/2
CoilNodes(3,1) = - RectangleCoilInfo(iCoil)%RectVec1(2)/2 - RectangleCoilInfo(iCoil)%RectVec2(2)/2
PointNumber = 1
DO iLoop=1,RectangleCoilInfo(iCoil)%LoopNum
  DO iPoint=2,RectangleCoilInfo(iCoil)%PointsPerLoop
    PointNumber = PointNumber + 1
    CoilNodes(1,PointNumber) = RectangleCoilInfo(iCoil)%Length * (PointNumber - 1) / (RectangleCoilInfo(iCoil)%NumNodes - 1)
    IF ((iPoint.GE.2).AND.(iPoint.LE.(1 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = - RectangleCoilInfo(iCoil)%RectVec1(1) / 2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec1(1) * (iPoint - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = - RectangleCoilInfo(iCoil)%RectVec1(2) / 2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec1(2) * (iPoint - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(1 * PointsPerSide + 2)).AND.(iPoint.LE.(2 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = RectangleCoilInfo(iCoil)%RectVec1(1) / 2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) /2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) * (iPoint - PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = RectangleCoilInfo(iCoil)%RectVec1(2) / 2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) /2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) * (iPoint - PointsPerSide - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(2 * PointsPerSide + 2)).AND.(iPoint.LE.(3 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = RectangleCoilInfo(iCoil)%RectVec1(1) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) /2 -&
                                 RectangleCoilInfo(iCoil)%RectVec1(1) * (iPoint - 2 * PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = RectangleCoilInfo(iCoil)%RectVec1(2) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) /2 -&
                                 RectangleCoilInfo(iCoil)%RectVec1(2) * (iPoint - 2 * PointsPerSide - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(3 * PointsPerSide + 2)).AND.(iPoint.LE.RectangleCoilInfo(iCoil)%PointsPerLoop)) THEN
      CoilNodes(2,PointNumber) = - RectangleCoilInfo(iCoil)%RectVec1(1) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) /2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(1) * (iPoint - 3 * PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = - RectangleCoilInfo(iCoil)%RectVec1(2) / 2 +&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) /2 -&
                                 RectangleCoilInfo(iCoil)%RectVec2(2) * (iPoint - 3 * PointsPerSide - 1) / PointsPerSide
    ENDIF
  ENDDO
ENDDO !iLoop

! Transform the nodes to the coil coordinate system
DO iPoint=1,RectangleCoilInfo(iCoil)%NumNodes
  CoilNodes(:,iPoint) = MATMUL(TrafoMatrix, CoilNodes(:,iPoint))
  CoilNodes(:,iPoint) = CoilNodes(:,iPoint) + RectangleCoilInfo(iCoil)%BasePoint(:)
ENDDO

END SUBROUTINE SetUpRectangleCoil

SUBROUTINE SetUpLinearConductor(iCoil)
!===================================================================================================================================
! Sets up a linear conductor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Coil_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iPoint
!-----------------------------------------------------------------------------------------------------------------------------------

ALLOCATE(CoilNodes(1:3,LinearConductorInfo(iCoil)%NumNodes))
DO iPoint=0,LinearConductorInfo(iCoil)%NumNodes - 1
  CoilNodes(:,iPoint+1) = LinearConductorInfo(iCoil)%BasePoint(:) + LinearConductorInfo(iCoil)%LengthVector(:) *&
                          iPoint / (LinearConductorInfo(iCoil)%NumNodes - 1)
ENDDO

END SUBROUTINE SetUpLinearConductor

SUBROUTINE FindLinIndependendVectors(NormalVector, Vector1, Vector2)
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
ENDIF

! Find the third vecord vector with the cross product
Vector2(1) = NormalVector(2)*Vector1(3) - NormalVector(3)*Vector1(2)
Vector2(2) = NormalVector(3)*Vector1(1) - NormalVector(1)*Vector1(3)
Vector2(3) = NormalVector(1)*Vector1(2) - NormalVector(2)*Vector1(1)

END SUBROUTINE FindLinIndependendVectors

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

SUBROUTINE RotateCoordinateSystem(Vector1, Vector2, Vector3, AxisVector)
!===================================================================================================================================
! Rotates the coordinate system around vector 1 till vector 2 is parallel to the base vector
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)    :: AxisVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT) :: Vector1(3),Vector2(3),Vector3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Phi
REAL                :: LengthVector2, LengthAxisVector
REAL                :: UnitVector1(3)
REAL                :: RotationMatrix(3,3)
!===================================================================================================================================

write(*,*) Vector1
write(*,*) Vector2
write(*,*) Vector3
read *
LengthVector2    = SQRT(Vector2(1)**2 + Vector2(2)**2 + Vector2(3)**2)
LengthAxisVector = SQRT(AxisVector(1)**2 + AxisVector(2)**2 + AxisVector(3)**2)
Phi = ACOS(DOT_PRODUCT(Vector2,AxisVector)/LengthVector2/LengthAxisVector)
write(*,*) Phi
read *
UnitVector1(:) = Vector1(:)/SQRT(Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2)

RotationMatrix(1,1) = UnitVector1(1)*UnitVector1(1) * (1 - COS(Phi)) + COS(Phi)
RotationMatrix(1,2) = UnitVector1(1)*UnitVector1(2) * (1 - COS(Phi)) - UnitVector1(3) * SIN(Phi)
RotationMatrix(1,3) = UnitVector1(1)*UnitVector1(3) * (1 - COS(Phi)) + UnitVector1(2) * SIN(Phi)
RotationMatrix(2,1) = UnitVector1(2)*UnitVector1(1) * (1 - COS(Phi)) + UnitVector1(3) * SIN(Phi)
RotationMatrix(2,2) = UnitVector1(2)*UnitVector1(2) * (1 - COS(Phi)) + COS(Phi)
RotationMatrix(2,3) = UnitVector1(2)*UnitVector1(3) * (1 - COS(Phi)) - UnitVector1(1) * SIN(Phi)
RotationMatrix(3,1) = UnitVector1(3)*UnitVector1(1) * (1 - COS(Phi)) - UnitVector1(2) * SIN(Phi)
RotationMatrix(3,2) = UnitVector1(3)*UnitVector1(2) * (1 - COS(Phi)) + UnitVector1(1) * SIN(Phi)
RotationMatrix(3,3) = UnitVector1(3)*UnitVector1(3) * (1 - COS(Phi)) + COS(Phi)

Vector2 = MATMUL(RotationMatrix,Vector2)
Vector3 = MATMUL(RotationMatrix,Vector3)

END SUBROUTINE RotateCoordinateSystem

SUBROUTINE BiotSavart(iCoil, coilType)
!===================================================================================================================================
! Calculates the magnetic field of a coil with Biot-Savarts law
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Preproc
USE MOD_Mesh_Vars,             ONLY: nElems, Elem_xGP
USE MOD_Coil_Vars
USE MOD_PICInterpolation_Vars, ONLY: BGField
USE MOD_Equation_Vars,  ONLY: mu0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
INTEGER, INTENT(IN) :: coilType
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, i, j, k
INTEGER             :: nNodes, iPoint
REAL                :: current, segmentLength, segmentToDOFLengthP3
REAL                :: currentDirection(3), segmentPosition(3), segmentToDOF(3), BFieldTemp(3)
!===================================================================================================================================

IF (coilType.EQ.1) THEN
  nNodes = CoilInfo(iCoil)%NumNodes
  current = CoilInfo(iCoil)%Current
ELSE IF (coilType.EQ.2) THEN
  nNodes = CircleCoilInfo(iCoil)%NumNodes
  current = CircleCoilInfo(iCoil)%Current
ELSE IF (coilType.EQ.3) THEN
  nNodes = RectangleCoilInfo(iCoil)%NumNodes
  current = RectangleCoilInfo(iCoil)%Current
ELSE IF (coilType.EQ.4) THEN
  nNodes = LinearConductorInfo(iCoil)%NumNodes
  current = LinearConductorInfo(iCoil)%Current
ENDIF

! Iterate over all DOFs
DO iElem=1,nElems
  DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
    ! Iterate over the coil segements
    DO iPoint=1,nNodes-1
      currentDirection(:)  = CoilNodes(:,iPoint + 1) - CoilNodes(:,iPoint)
      segmentPosition(:)   = CoilNodes(:,iPoint) + currentDirection(:)/2
      segmentLength        = SQRT(currentDirection(1)**2 + currentDirection(2)**2 + currentDirection(3)**2)
      currentDirection(:)  = currentDirection(:) / segmentLength

      segmentToDOF(:)      = Elem_xGP(:,i,j,k,iElem) - segmentPosition(:)
      segmentToDOFLengthP3 = SQRT(segmentToDOF(1)**2 + segmentToDOF(2)**2 + segmentToDOF(3)**2)**3

      ! Calculate the cross product between the currentDirection and the segmentToDOF vectors
      BFieldTemp(1)        = currentDirection(2)*segmentToDOF(3) - currentDirection(3)*segmentToDOF(2)
      BFieldTemp(2)        = currentDirection(3)*segmentToDOF(1) - currentDirection(1)*segmentToDOF(3)
      BFieldTemp(3)        = currentDirection(1)*segmentToDOF(2) - currentDirection(2)*segmentToDOF(1)

      BFieldTemp(:)        = BFieldTemp(:) * mu0 * current * segmentLength /&
                             (4 * PI * segmentToDOFLengthP3)

      BGField(1:3,i,j,k,iElem) = BGField(1:3,i,j,k,iElem) + BFieldTemp(:)
    ENDDO
  ENDDO; ENDDO; ENDDO
ENDDO

END SUBROUTINE BiotSavart

SUBROUTINE Jefimenko(iCoil, t, OutFlag)
!===================================================================================================================================
! Calculates the magnetic field of a coil with a time-dependent current with Jefimenkos Law
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars, ONLY: PI
USE MOD_Preproc
USE MOD_Mesh_Vars, ONLY: nElems, Elem_xGP
USE MOD_Equation_Vars, ONLY: mu0, c
USE MOD_Coil_Vars, ONLY: CoilInfo, CircleCoilInfo, RectangleCoilInfo, LinearConductorInfo, &
                        CurrentInfo, NumOfCoils, NumOfCircleCoils, NumOfRectangleCoils, CoilNodes
USE MOD_PICInterpolation_Vars, ONLY: BGField
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil, OutFlag
REAL, INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: current, timedercurrent, segmentLength, segmentToDOFLength
INTEGER             :: nNodes, iElem, i, j, k, iPoint
REAL                :: currentDirection(3), segmentPosition(3), segmentToDOF(3), BFieldTemp(3)
!===================================================================================================================================

IF (OutFlag.EQ.1) THEN
  nNodes = CoilInfo(iCoil)%NumNodes
  current = CurrentInfo(iCoil)%CurrentAmpl*SIN(2*PI*CurrentInfo(iCoil)%CurrentFreq*t + CurrentInfo(iCoil)%CurrentPhase)
  timedercurrent = 2*Pi*CurrentInfo(iCoil)%CurrentFreq*CurrentInfo(iCoil)%CurrentAmpl *&
                   COS(2*PI*CurrentInfo(iCoil)%CurrentFreq*t + CurrentInfo(iCoil)%CurrentPhase)
ELSE IF (OutFlag.EQ.2) THEN
  nNodes = CircleCoilInfo(iCoil)%NumNodes
  current = CurrentInfo(iCoil+NumOfCoils)%CurrentAmpl*SIN(2*PI*CurrentInfo(iCoil+NumOfCoils)%CurrentFreq*t +&
            CurrentInfo(iCoil+NumOfCoils)%CurrentPhase)
  timedercurrent = 2*Pi*CurrentInfo(iCoil+NumOfCoils)%CurrentFreq*CurrentInfo(iCoil+NumOfCoils)%CurrentAmpl *&
                   COS(2*PI*CurrentInfo(iCoil+NumOfCoils)%CurrentFreq*t + CurrentInfo(iCoil+NumOfCoils)%CurrentPhase)
ELSE IF (OutFlag.EQ.3) THEN
  nNodes = RectangleCoilInfo(iCoil)%NumNodes
  current = CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentAmpl *&
            SIN(2*PI*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentFreq*t +&
            CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentPhase)
  timedercurrent = 2*Pi*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentFreq *&
                   CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentAmpl *&
                   COS(2*PI*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentFreq*t +&
                   CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils)%CurrentPhase)
ELSE IF (OutFlag.EQ.4) THEN
  nNodes = LinearConductorInfo(iCoil)%NumNodes
  current = CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentAmpl *&
            SIN(2*PI*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentFreq*t +&
            CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentPhase)
  timedercurrent = 2*Pi*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentFreq *&
                   CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentAmpl *&
                   COS(2*PI*CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentFreq*t +&
                   CurrentInfo(iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils)%CurrentPhase)
ENDIF

! Iterate over all DOFs
DO iElem=1,nElems
  DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
    ! Iterate over the coil segements
    DO iPoint=1,nNodes-1
      currentDirection(:)  = CoilNodes(:,iPoint + 1) - CoilNodes(:,iPoint)
      segmentPosition(:)   = CoilNodes(:,iPoint) + currentDirection(:)/2
      segmentLength        = SQRT(currentDirection(1)**2 + currentDirection(2)**2 + currentDirection(3)**2)
      currentDirection(:)  = currentDirection(:) / segmentLength

      segmentToDOF(:)      = Elem_xGP(:,i,j,k,iElem) - segmentPosition(:)
      segmentToDOFLength   = SQRT(segmentToDOF(1)**2 + segmentToDOF(2)**2 + segmentToDOF(3)**2)

      ! Calculate the cross product between the currentDirection and the segmentToDOF vectors
      BFieldTemp(1)        = currentDirection(2)*segmentToDOF(3) - currentDirection(3)*segmentToDOF(2)
      BFieldTemp(2)        = currentDirection(3)*segmentToDOF(1) - currentDirection(1)*segmentToDOF(3)
      BFieldTemp(3)        = currentDirection(1)*segmentToDOF(2) - currentDirection(2)*segmentToDOF(1)

      BFieldTemp(:)        = BFieldTemp(:) * mu0 * (current / segmentToDOFLength**3 +&
                             timedercurrent / (segmentToDOFLength**2 * c)) * segmentLength / (4 * PI)

      BGField(:,i,j,k,iElem) = BGField(:,i,j,k,iElem) + BFieldTemp(:)
    ENDDO
  ENDDO; ENDDO; ENDDO
ENDDO

END SUBROUTINE Jefimenko

SUBROUTINE WriteCoilVTK(iCoil, OutFlag)
!===================================================================================================================================
! Write Coil to VTK
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Coil_Vars, ONLY: CoilNodes, CoilInfo, CircleCoilInfo, RectangleCoilInfo, LinearConductorInfo,&
                         NumOfCoils, NumOfCircleCoils, NumOfRectangleCoils
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil, OutFlag ! 1:coil, 2:circleCoil, 3: rectangleCoil, 4: linearConductor
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=26)  :: myFileName
INTEGER            :: iNode, iElem, nNodes
!===================================================================================================================================

IF (OutFlag.EQ.1) THEN
  WRITE(myFileName,'(A9,I3.3,A4)')'CoilMesh_',iCoil,'.vtk'
  nNodes = CoilInfo(iCoil)%NumNodes
ELSE IF (OutFlag.EQ.2) THEN
  WRITE(myFileName,'(A9,I3.3,A4)')'CoilMesh_',iCoil+NumOfCoils,'.vtk'
  nNodes = CircleCoilInfo(iCoil)%NumNodes
ELSE IF (OutFlag.EQ.3) THEN
  WRITE(myFileName,'(A9,I3.3,A4)')'CoilMesh_',iCoil+NumOfCoils+NumOfCircleCoils,'.vtk'
  nNodes = RectangleCoilInfo(iCoil)%NumNodes
ELSE IF (OutFlag.EQ.4) THEN
  WRITE(myFileName,'(A9,I3.3,A4)')'CoilMesh_',iCoil+NumOfCoils+NumOfCircleCoils+NumOfRectangleCoils,'.vtk'
  nNodes = LinearConductorInfo(iCoil)%NumNodes
ENDIF

OPEN(1112,FILE=myFileName,STATUS='replace')
WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
WRITE(1112,'(A)')'Debug Mesh'
WRITE(1112,'(A)')'ASCII'
WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
WRITE(1112,'(A)')''
WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
DO iNode=1,nNodes
  WRITE(1112,*) CoilNodes(1:3, iNode)
ENDDO
WRITE(1112,*)''
WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nNodes - 1,3*(nNodes - 1)
DO iElem=1,nNodes - 1
  WRITE(1112,'(I0)',ADVANCE="NO")2
  WRITE(1112,'(1X,I0)',ADVANCE="NO") iElem - 1
  WRITE(1112,'(1X,I0)',ADVANCE="NO") iElem
  WRITE(1112,*)''
ENDDO
WRITE(1112,*)''
WRITE(1112,'(A,I0)')'CELL_TYPES ',nNodes - 1
DO iElem=1,nNodes - 1
  WRITE(1112,'(1X,I0)',ADVANCE="NO")3
ENDDO
WRITE(1112,*)''

CLOSE(1112)

END SUBROUTINE WriteCoilVTK

SUBROUTINE FinalizeCoil()
!===================================================================================================================================
! Finalize coil
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Coil_Vars, ONLY: CoilNodes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( CoilNodes)

END SUBROUTINE FinalizeCoil

END MODULE MOD_SuperB_Coil




