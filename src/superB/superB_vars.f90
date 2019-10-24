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

MODULE MOD_SuperB_Vars
!===================================================================================================================================
! Contains the global variables for different coils and permanent magnets
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

! === Coils
INTEGER :: NumOfCoils
INTEGER :: NumOfCircleCoils
INTEGER :: NumOfRectangleCoils
INTEGER :: NumOfLinearConductors

TYPE tSegmentInfo
  INTEGER :: SegmentType ! 1: Line, 2: Circle Segment
  INTEGER :: NumOfPoints
  REAL    :: LineVector(2)
  REAL    :: Radius
  REAL    :: Phi1
  REAL    :: Phi2
END TYPE tSegmentInfo

TYPE tCoilInfo
  INTEGER                        :: NumNodes
  REAL                           :: Current
  INTEGER                        :: LoopNum
  INTEGER                        :: PointsPerLoop
  INTEGER                        :: NumOfSegments
  REAL                           :: LengthVector(3)
  REAL                           :: Length
  REAL                           :: AxisVec1(3)
  REAL                           :: BasePoint(3)
  TYPE(tSegmentInfo),ALLOCATABLE :: SegmentInfo(:)
  REAL                           :: LoopLength
END TYPE tCoilInfo

TYPE tCircleCoil
  INTEGER                        :: NumNodes
  REAL                           :: Current
  INTEGER                        :: LoopNum
  INTEGER                        :: PointsPerLoop
  REAL                           :: Radius
  REAL                           :: BasePoint(3)
  REAL                           :: LengthVector(3)
  REAL                           :: Length
END TYPE tCircleCoil

TYPE tRectangleCoil
  INTEGER                        :: NumNodes
  REAL                           :: Current
  INTEGER                        :: LoopNum
  INTEGER                        :: PointsPerLoop
  REAL                           :: RectVec1(2)
  REAL                           :: RectVec2(2)
  REAL                           :: BasePoint(3)
  REAL                           :: LengthVector(3)
  REAL                           :: Length
  REAL                           :: AxisVec1(3)
END TYPE tRectangleCoil

TYPE tLinearConductor
  INTEGER                        :: NumNodes
  REAL                           :: Current
  REAL                           :: BasePoint(3)
  REAL                           :: LengthVector(3)
END TYPE tLinearConductor

TYPE(tCoilInfo),ALLOCATABLE      :: CoilInfo(:)
TYPE(tCircleCoil),ALLOCATABLE    :: CircleCoilInfo(:)
TYPE(tRectangleCoil),ALLOCATABLE :: RectangleCoilInfo(:)
TYPE(tLinearConductor),ALLOCATABLE :: LinearConductorInfo(:)

REAL, ALLOCATABLE                :: CoilNodes(:,:)

TYPE tCurrentInfo
  REAL                           :: CurrentAmpl
  REAL                           :: CurrentFreq
  REAL                           :: CurrentPhase
END TYPE tCurrentInfo

TYPE(tCurrentInfo),ALLOCATABLE   :: CurrentInfo(:)

LOGICAL, ALLOCATABLE             :: TimeDepCoil(:)
INTEGER                           :: nTimePoints
REAL, ALLOCATABLE                 :: BGFieldTDep(:,:,:,:,:,:)      !< Time dep. BGField data (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems,1:nTime)

! === Permanent Magnets

! INPUT
INTEGER :: NumOfCuboidMagnets
INTEGER :: NumOfSphericMagnets
INTEGER :: NumOfCylindricMagnets
INTEGER :: NumOfConicMagnets

TYPE tCuboidMagnetInfo
  REAL    :: BasePoint(3)
  REAL    :: BaseVector1(3)
  REAL    :: BaseVector2(3)
  REAL    :: BaseVector3(3)
  INTEGER :: NumNodes
  REAL    :: Magnetisation(3)
END TYPE tCuboidMagnetInfo

TYPE(tCuboidMagnetInfo),ALLOCATABLE :: CuboidMagnetInfo(:)

TYPE tSphericMagnetInfo
  REAL    :: BasePoint(3)
  REAL    :: Radius
  INTEGER :: NumNodes
  REAL    :: Magnetisation(3)
END TYPE tSphericMagnetInfo

TYPE(tSphericMagnetInfo),ALLOCATABLE :: SphericMagnetInfo(:)

TYPE tCylindricMagnetInfo
  REAL    :: BasePoint(3)
  REAL    :: HeightVector(3)
  REAL    :: Radius
  INTEGER :: NumNodes
  REAL    :: Magnetisation(3)
END TYPE tCylindricMagnetInfo

TYPE(tCylindricMagnetInfo),ALLOCATABLE :: CylindricMagnetInfo(:)

TYPE tConicMagnetInfo
  REAL    :: BasePoint(3)
  REAL    :: HeightVector(3)
  REAL    :: Radius1
  REAL    :: Radius2
  INTEGER :: NumNodes
  REAL    :: Magnetisation(3)
END TYPE tConicMagnetInfo

TYPE(tConicMagnetInfo),ALLOCATABLE :: ConicMagnetInfo(:)

REAL, ALLOCATABLE                   :: PsiMag(:,:,:,:)
INTEGER, ALLOCATABLE                :: MagnetFlag(:,:,:,:)

END MODULE MOD_SuperB_Vars