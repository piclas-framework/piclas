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

MODULE MOD_SuperB_Vars
!===================================================================================================================================
! Contains the global variables for different coils and permanent magnets
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

LOGICAL           :: DoCalcErrorNormsSuperB                     !< perform L2, LInf error calculation
REAL              :: L_2_ErrorSuperB(4)                         !< L2 error for magnetic field Bx, By, Bz
REAL              :: L_Inf_ErrorSuperB(4)                       !< LInf error for magnetic field Bx, By, Bz

! === Coils
INTEGER :: NumOfCoils                                           !< Total number of coils

TYPE tSegmentInfo
  CHARACTER(LEN=255)              :: SegmentType                !< 1: Line, 2: Circle Segment
  INTEGER                         :: NumOfPoints                !< Number of discretization points per segment
  REAL                            :: LineVector(2)              !< 2D vector defining the line segment in the cross-section plane
  REAL                            :: Radius                     !< Radial of circle segment [m]
  REAL                            :: Phi1                       !< Initial angle of circle segment [deg]
  REAL                            :: Phi2                       !< Final angle of circle segment [deg]
END TYPE tSegmentInfo

TYPE tCoilInfo
  CHARACTER(LEN=255)              :: Type                       !< "custom", "circle", "rectangle" cross-section, "linear" conductor
  INTEGER                         :: NumNodes                   !< Discretization of the linear conductor
  REAL                            :: Current                    !< Current in [A]
  INTEGER                         :: LoopNum                    !< Number of coils loops
  INTEGER                         :: PointsPerLoop              !< Number of discretization points per loop
  REAL                            :: BasePoint(3)               !< Origin vector for the coil
  REAL                            :: LengthVector(3)            !< Length vector normal to the coil cross-section
  REAL                            :: AxisVec1(3)                !< Axial vector orthogonal to the length vector (custom, cuboid)
  INTEGER                         :: NumOfSegments              !< Custom coil: Number of segments
  TYPE(tSegmentInfo),ALLOCATABLE  :: SegmentInfo(:)             !< Custom coil: Container for segment information [1:NumOfSegments]
  REAL                            :: Radius                     !< Circular coil-specific
  REAL                            :: RectVec1(2)                !< Cuboid coil-specific vector in the cross-section plane
  REAL                            :: RectVec2(2)                !< Cuboid coil-specific vector in the cross-section plane
  REAL                            :: Length                     !< Length of coil, calculated from the length vector
END TYPE tCoilInfo

TYPE(tCoilInfo),ALLOCATABLE       :: CoilInfo(:)                !< Container for the coil information [1:NumOfCoils]

REAL, ALLOCATABLE                 :: CoilNodes(:,:)             !< Geometric information of the coils [1:3,CoilInfo(iCoil)%NumNodes]

! === Time-dependent coils
TYPE tCurrentInfo
  REAL                            :: CurrentFreq                !< Current frequency [1/s]
  REAL                            :: CurrentPhase               !< Current phase [rad]
END TYPE tCurrentInfo

TYPE(tCurrentInfo),ALLOCATABLE    :: CurrentInfo(:)             !< Container for the time-dependent coil information [1:NumOfCoils]

REAL                              :: BGFieldFrequency           !< Frequency of time-dependent coil (written as attribute to h5)
REAL                              :: BGFieldCurrent             !< Current maximum of time-dependent coil (written as attribute to h5)

LOGICAL, ALLOCATABLE              :: TimeDepCoil(:)             !< Flag if the coil has a time-dependent current [1:NumOfCoils]
LOGICAL                           :: UseTimeDepCoil             !< Flag if any coil has a time-dependent current
INTEGER                           :: nTimePoints                !< Number of time discretization points for the sinusoidal curve
REAL, ALLOCATABLE                 :: BGFieldTDep(:,:,:,:,:,:)   !< Time-dependent background field
                                                                !< [1:3,0:NBG,0:NBG,0:NBG,1:PP_nElems,1:nTimePoints]

! === Permanent Magnets

INTEGER                 :: NumOfPermanentMagnets                !< Total number of permanent magnets
INTEGER, ALLOCATABLE    :: MagnetFlag(:,:,:,:)                  !< Number of the magnet that occupies the point, otherwise zero
                                                                !< [0:PP_N,0:PP_N,0:PP_N,1:nElems]

TYPE tPermanentMagnetInfo
  CHARACTER(LEN=255)    :: Type                                 !< Cuboid, sphere, cylinder, conic
  REAL                  :: BasePoint(3)                         !< Origin vector for the permanent magnet
  REAL                  :: BaseVector1(3)                       !< Vector 1 spanning the cuboid permanent magnet
  REAL                  :: BaseVector2(3)                       !< Vector 2 spanning the cuboid permanent magnet
  REAL                  :: BaseVector3(3)                       !< Vector 3 spanning the cuboid permanent magnet
  INTEGER               :: NumNodes                             !< Number of nodes for the discretization:
                                                                !< Cuboid: N points in each direction (total number: 6N^2)
                                                                !< Sphere: N divisions in the zenith direction with 2*N points in
                                                                !<          the azimuthal direction
                                                                !< Cylinder: N divisions along height vector, 2*N points in the
                                                                !<            azimuthal direction, N points in radial direction on
                                                                !<            the top and bottom face
                                                                !< Conic: analogous to the cylinder
  REAL                  :: Magnetisation(3)                     !< Magnetisation vector in [A/m]
  REAL                  :: Radius                               !< Radius in [m] for sphere, cylinder, conic magnets
  REAL                  :: HeightVector(3)                      !< Height vector [m] for cylinder, conic magnets
  REAL                  :: Radius2                              !< Second radius for conic magnets or inner radius for hollow cylinders
END TYPE tPermanentMagnetInfo

TYPE(tPermanentMagnetInfo),ALLOCATABLE :: PermanentMagnetInfo(:)!< Container for the permanent magnet info [1:NumOfPermanentMagnets]

END MODULE MOD_SuperB_Vars
