!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_SuperB_Coil
!===================================================================================================================================
!> Contains the calculation of the magnetic field of coils
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: SetUpCoils, SetUpCoil, SetUpCircleCoil, SetUpRectangleCoil, SetUpLinearConductor, BiotSavart, Jefimenko, WriteCoilVTK
!===================================================================================================================================

INTERFACE SetUpCoils
  MODULE PROCEDURE SetUpCoils
END INTERFACE SetUpCoils

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

INTERFACE BiotSavart
  MODULE PROCEDURE BiotSavart
END INTERFACE BiotSavart

INTERFACE Jefimenko
  MODULE PROCEDURE Jefimenko
END INTERFACE Jefimenko

!===================================================================================================================================

CONTAINS

SUBROUTINE SetUpCoils(iCoil)
!===================================================================================================================================
!> Sets up a coil defined by segments (line and circle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_SuperB_Vars
USE MOD_Interpolation_Vars    ,ONLY: BGFieldVTKOutput
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,I2)',ADVANCE='NO') ' Set up coil: ', iCoil

! Select different coil types
SELECT CASE(TRIM(CoilInfo(iCoil)%Type))
CASE('custom')
  CALL SetUpCoil(iCoil)
CASE('circle')
  CALL SetUpCircleCoil(iCoil)
CASE('rectangle')
  CALL SetUpRectangleCoil(iCoil)
CASE('linear')
  CALL SetUpLinearConductor(iCoil)
CASE DEFAULT
  CALL abort(__STAMP__,'Unknown (time-dependent) coil type ['//TRIM(CoilInfo(iCoil)%Type)//']')
END SELECT

! Output to VTK
IF(BGFieldVTKOutput) CALL WriteCoilVTK(iCoil)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') 'DONE'

END SUBROUTINE SetUpCoils


SUBROUTINE SetUpCoil(iCoil)
!===================================================================================================================================
!> Sets up a coil defined by segments (line and circle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_SuperB_Vars
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
      SELECT CASE(TRIM(CoilInfo(iCoil)%SegmentInfo(iSegment)%SegmentType))
      CASE('line')
        CoilNodes(2,PointNumber) = CoilInfo(iCoil)%SegmentInfo(iSegment)%LineVector(1) * iPoint /&
            (CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1) + LastSegPoint(2)
        CoilNodes(3,PointNumber) = CoilInfo(iCoil)%SegmentInfo(iSegment)%LineVector(2) * iPoint /&
            (CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1) + LastSegPoint(3)
      CASE('circle')
        Phi1   = CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi1
        Phi2   = CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi2
        Radius = CoilInfo(iCoil)%SegmentInfo(iSegment)%Radius
        CoilNodes(2,PointNumber) = LastSegPoint(2) - COS(Phi1)*Radius +&
            COS((Phi2 - Phi1)/(CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1)*&
            iPoint + Phi1)*Radius
        CoilNodes(3,PointNumber) = LastSegPoint(3) - SIN(Phi1)*Radius +&
            SIN((Phi2 - Phi1)/(CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1)*&
            iPoint + Phi1)*Radius
      END SELECT
    END DO !iPoint
    LastSegPoint(:) = CoilNodes(:,PointNumber)
  END DO !iSegment
END DO !iLoop

! Transform the nodes to the coil coordinate system
DO iPoint=1,CoilInfo(iCoil)%NumNodes
  CoilNodes(:,iPoint) = MATMUL(TrafoMatrix, CoilNodes(:,iPoint))
  CoilNodes(:,iPoint) = CoilNodes(:,iPoint) + CoilInfo(iCoil)%BasePoint(:)
END DO

END SUBROUTINE SetUpCoil


SUBROUTINE SetUpCircleCoil(iCoil)
!===================================================================================================================================
!> Sets up a coil consisting of a circle cross-section
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_SuperB_Vars
USE MOD_SuperB_Tools  ,ONLY: FindLinIndependentVectors, GramSchmidtAlgo
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

CALL FindLinIndependentVectors(CoilInfo(iCoil)%LengthVector, AxisVec1, AxisVec2)
CALL GramSchmidtAlgo(CoilInfo(iCoil)%LengthVector, AxisVec1, AxisVec2)
TrafoMatrix(:,1) = CoilInfo(iCoil)%LengthVector
TrafoMatrix(:,2) = AxisVec1
TrafoMatrix(:,3) = AxisVec2

ALLOCATE(CoilNodes(1:3,CoilInfo(iCoil)%NumNodes))

DO iPoint=0,CoilInfo(iCoil)%NumNodes - 1
  CoilNodes(1, iPoint+1) = CoilInfo(iCoil)%Length * iPoint / (CoilInfo(iCoil)%NumNodes - 1)
  CoilNodes(2, iPoint+1) = CoilInfo(iCoil)%Radius * COS(2*PI*iPoint/CoilInfo(iCoil)%PointsPerLoop)
  CoilNodes(3, iPoint+1) = CoilInfo(iCoil)%Radius * SIN(2*PI*iPoint/CoilInfo(iCoil)%PointsPerLoop)
  CoilNodes(:, iPoint+1) = MATMUL(TrafoMatrix, CoilNodes(:, iPoint+1))
  CoilNodes(:, iPoint+1) = CoilNodes(:, iPoint+1) + CoilInfo(iCoil)%BasePoint(:)
END DO

END SUBROUTINE SetUpCircleCoil

SUBROUTINE SetUpRectangleCoil(iCoil)
!===================================================================================================================================
! Sets up a coil consisting of a rectangle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_SuperB_Vars
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
REAL                :: LengthAxisVec1, LengthAxisVec2
INTEGER             :: PointsPerSide, PointNumber, iLoop, iPoint
REAL                :: TrafoMatrix(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------

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
AxisVec2                              = AxisVec2 / LengthAxisVec2
! Transformation matrix
TrafoMatrix(:,1) = CoilInfo(iCoil)%LengthVector
TrafoMatrix(:,2) = CoilInfo(iCoil)%AxisVec1
TrafoMatrix(:,3) = AxisVec2

ALLOCATE(CoilNodes(1:3,CoilInfo(iCoil)%NumNodes))

PointsPerSide = (CoilInfo(iCoil)%PointsPerLoop - 1) / 4

! Start at (0,0,0) in the untransformed coordinate system
CoilNodes(1,1) = 0
CoilNodes(2,1) = - CoilInfo(iCoil)%RectVec1(1)/2 - CoilInfo(iCoil)%RectVec2(1)/2
CoilNodes(3,1) = - CoilInfo(iCoil)%RectVec1(2)/2 - CoilInfo(iCoil)%RectVec2(2)/2
PointNumber = 1
DO iLoop=1,CoilInfo(iCoil)%LoopNum
  DO iPoint=2,CoilInfo(iCoil)%PointsPerLoop
    PointNumber = PointNumber + 1
    CoilNodes(1,PointNumber) = CoilInfo(iCoil)%Length * (PointNumber - 1) / (CoilInfo(iCoil)%NumNodes - 1)
    IF ((iPoint.GE.2).AND.(iPoint.LE.(1 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = - CoilInfo(iCoil)%RectVec1(1) / 2. -&
          CoilInfo(iCoil)%RectVec2(1) / 2. +&
          CoilInfo(iCoil)%RectVec1(1) * (iPoint - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = - CoilInfo(iCoil)%RectVec1(2) / 2. -&
          CoilInfo(iCoil)%RectVec2(2) / 2. +&
          CoilInfo(iCoil)%RectVec1(2) * (iPoint - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(1 * PointsPerSide + 2)).AND.(iPoint.LE.(2 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = CoilInfo(iCoil)%RectVec1(1) / 2. -&
          CoilInfo(iCoil)%RectVec2(1) /2 +&
          CoilInfo(iCoil)%RectVec2(1) * (iPoint - PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = CoilInfo(iCoil)%RectVec1(2) / 2. -&
          CoilInfo(iCoil)%RectVec2(2) /2 +&
          CoilInfo(iCoil)%RectVec2(2) * (iPoint - PointsPerSide - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(2 * PointsPerSide + 2)).AND.(iPoint.LE.(3 * PointsPerSide + 1))) THEN
      CoilNodes(2,PointNumber) = CoilInfo(iCoil)%RectVec1(1) / 2. +&
          CoilInfo(iCoil)%RectVec2(1) /2 -&
          CoilInfo(iCoil)%RectVec1(1) * (iPoint - 2 * PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = CoilInfo(iCoil)%RectVec1(2) / 2. +&
          CoilInfo(iCoil)%RectVec2(2) /2 -&
          CoilInfo(iCoil)%RectVec1(2) * (iPoint - 2 * PointsPerSide - 1) / PointsPerSide
    ELSE IF ((iPoint.GE.(3 * PointsPerSide + 2)).AND.(iPoint.LE.CoilInfo(iCoil)%PointsPerLoop)) THEN
      CoilNodes(2,PointNumber) = - CoilInfo(iCoil)%RectVec1(1) / 2. +&
          CoilInfo(iCoil)%RectVec2(1) /2 -&
          CoilInfo(iCoil)%RectVec2(1) * (iPoint - 3 * PointsPerSide - 1) / PointsPerSide
      CoilNodes(3,PointNumber) = - CoilInfo(iCoil)%RectVec1(2) / 2. +&
          CoilInfo(iCoil)%RectVec2(2) /2 -&
          CoilInfo(iCoil)%RectVec2(2) * (iPoint - 3 * PointsPerSide - 1) / PointsPerSide
    END IF
  END DO
END DO !iLoop

! Transform the nodes to the coil coordinate system
DO iPoint=1,CoilInfo(iCoil)%NumNodes
  CoilNodes(:,iPoint) = MATMUL(TrafoMatrix, CoilNodes(:,iPoint))
  CoilNodes(:,iPoint) = CoilNodes(:,iPoint) + CoilInfo(iCoil)%BasePoint(:)
END DO

END SUBROUTINE SetUpRectangleCoil


SUBROUTINE SetUpLinearConductor(iCoil)
!===================================================================================================================================
!> Sets up a linear conductor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_SuperB_Vars
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

ALLOCATE(CoilNodes(1:3,CoilInfo(iCoil)%NumNodes))
DO iPoint=0,CoilInfo(iCoil)%NumNodes - 1
  CoilNodes(:,iPoint+1) = CoilInfo(iCoil)%BasePoint(:) + CoilInfo(iCoil)%LengthVector(:) * iPoint / (CoilInfo(iCoil)%NumNodes - 1)
END DO

END SUBROUTINE SetUpLinearConductor


SUBROUTINE BiotSavart(iCoil)
!===================================================================================================================================
!> Calculates the magnetic field of a coil with Biot-Savart's law
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Preproc
USE MOD_SuperB_Vars
USE MOD_Mesh_Vars          ,ONLY: nElems, Elem_xGP
USE MOD_Interpolation_Vars ,ONLY: BGField
USE MOD_Globals_Vars       ,ONLY: mu0
USE MOD_SuperB_Tools       ,ONLY: CalcErrorSuperB
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, i, j, k, iPoint
REAL                :: segmentLength, segmentToDOFLengthP3, currentDirection(3), segmentPosition(3), segmentToDOF(3), BFieldTemp(3)
CHARACTER(LEN=40)   :: formatStr
INTEGER             :: ExactFunctionNumber    ! Number of exact function to be used for the calculation of the analytical solution
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') '...Calculation of the B-Field'

! Iterate over all DOFs
DO iElem=1,nElems
  DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
    ! Iterate over the coil segements
    DO iPoint=1,CoilInfo(iCoil)%NumNodes-1
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

      BFieldTemp(:)        = BFieldTemp(:) * mu0 * CoilInfo(iCoil)%Current * segmentLength /&
          (4 * PI * segmentToDOFLengthP3)

      BGField(1:3,i,j,k,iElem) = BGField(1:3,i,j,k,iElem) + BFieldTemp(:)
    END DO
  END DO; END DO; END DO
END DO

SDEALLOCATE(CoilNodes)

IF(DoCalcErrorNormsSuperB)THEN

  ! Check pre-defined cases
  SELECT CASE(TRIM(CoilInfo(iCoil)%Type))
  !CASE('custom')
  !CASE('circle')
  !CASE('rectangle')
  CASE('linear')
    ExactFunctionNumber = 10
  CASE DEFAULT
    CALL abort(&
        __STAMP__&
        ,'Cannot calculate L2/LInf error for coil type ['//TRIM(CoilInfo(iCoil)%Type)//']')
  END SELECT

  ! Get L2 errors
  CALL CalcErrorSuperB(L_2_ErrorSuperB,L_Inf_ErrorSuperB,ExactFunctionNumber,iCoil)

  ! Graphical output
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A25,',4,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' L_2_ErrorSuperB       : ',L_2_ErrorSuperB
    WRITE(UNIT_StdOut,formatStr)' L_Inf_ErrorSuperB     : ',L_Inf_ErrorSuperB
  END IF
END IF ! DoCalcErrorNormsSuperB

END SUBROUTINE BiotSavart


SUBROUTINE Jefimenko(iCoil, timestep, iTimePoint)
!===================================================================================================================================
!> Calculates the magnetic field of a coil with a time-dependent current with Jefimenko's Law
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nElems, Elem_xGP
USE MOD_Globals_Vars       ,ONLY: mu0, c
USE MOD_SuperB_Vars        ,ONLY: CoilInfo, CurrentInfo, CoilNodes, BGFieldTDep
!USE MOD_Interpolation_Vars  ,ONLY: BGField
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
REAL, INTENT(IN)    :: timestep
INTEGER, INTENT(IN) :: iTimePoint
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: current, timedercurrent, segmentLength, segmentToDOFLength
INTEGER             :: iElem, i, j, k, iPoint
REAL                :: currentDirection(3), segmentPosition(3), segmentToDOF(3), BFieldTemp(3)
REAL                :: omega
!===================================================================================================================================

ASSOCIATE( f   => CurrentInfo(iCoil)%CurrentFreq   ,&
           phi => CurrentInfo(iCoil)%CurrentPhase  ,&
           A   => CoilInfo(iCoil)%Current          ,&
           t   => timestep*REAL(iTimePoint-1)      )
  IF(f.GT.0.)THEN
    omega          = 2.*PI*f
    current        = A*SIN(omega*t+phi)
    timedercurrent = omega*A*COS(omega*t+phi)
  ELSE
    current = A
    timedercurrent = 0.
  END IF ! f.LE.0
END ASSOCIATE

! Iterate over all DOFs
DO iElem=1,nElems
  DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
    ! Iterate over the coil segments
    DO iPoint=1,CoilInfo(iCoil)%NumNodes-1
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

      BGFieldTDep(:,i,j,k,iElem,iTimePoint) = BGFieldTDep(:,i,j,k,iElem,iTimePoint) + BFieldTemp(:)
    END DO
  END DO; END DO; END DO
END DO


END SUBROUTINE Jefimenko


SUBROUTINE WriteCoilVTK(iCoil)
!===================================================================================================================================
!> Write Coil to VTK
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SuperB_Vars ,ONLY: CoilNodes, CoilInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCoil
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=26)  :: myFileName
INTEGER            :: iNode, iElem, nNodes
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,I0,A)',ADVANCE='NO') ' ... Write coil #',iCoil,' to VTK file ...'
WRITE(myFileName,'(A9,I3.3,A4)')'CoilMesh_',iCoil,'.vtk'
nNodes = CoilInfo(iCoil)%NumNodes

OPEN(1112,FILE=myFileName,STATUS='replace')
WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
WRITE(1112,'(A)')'Debug Mesh'
WRITE(1112,'(A)')'ASCII'
WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
WRITE(1112,'(A)')''
WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
DO iNode=1,nNodes
  WRITE(1112,*) CoilNodes(1:3, iNode)
END DO
WRITE(1112,*)''
WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nNodes - 1,3*(nNodes - 1)
DO iElem=1,nNodes - 1
  WRITE(1112,'(I0)',ADVANCE="NO")2
  WRITE(1112,'(1X,I0)',ADVANCE="NO") iElem - 1
  WRITE(1112,'(1X,I0)',ADVANCE="NO") iElem
  WRITE(1112,*)''
END DO
WRITE(1112,*)''
WRITE(1112,'(A,I0)')'CELL_TYPES ',nNodes - 1
DO iElem=1,nNodes - 1
  WRITE(1112,'(1X,I0)',ADVANCE="NO")3
END DO
WRITE(1112,*)''

CLOSE(1112)

END SUBROUTINE WriteCoilVTK


END MODULE MOD_SuperB_Coil
