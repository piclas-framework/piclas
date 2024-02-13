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

MODULE MOD_Particle_InterSection
!===================================================================================================================================
! Provides routines to calculate the intersection of the particle trajectory with a side depending on the side type
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ComputePlanarCurvedIntersection
  MODULE PROCEDURE ComputePlanarCurvedIntersection
END INTERFACE

INTERFACE ComputePlanarRectInterSection
  MODULE PROCEDURE ComputePlanarRectInterSection
END INTERFACE

INTERFACE ComputeBilinearIntersection
  MODULE PROCEDURE ComputeBilinearIntersection
END INTERFACE

INTERFACE ComputeCurvedIntersection
  MODULE PROCEDURE ComputeCurvedIntersection
END INTERFACE

#ifdef CODE_ANALYZE
INTERFACE OutputTrajectory
  MODULE PROCEDURE OutputTrajectory
END INTERFACE
#endif /*CODE_ANALYZE*/

PUBLIC :: IntersectionWithWall
PUBLIC :: ComputePlanarRectInterSection
PUBLIC :: ComputePlanarCurvedIntersection
PUBLIC :: ComputeBilinearIntersection
PUBLIC :: ComputeCurvedIntersection
PUBLIC :: ParticleThroughSideCheck3DFast
PUBLIC :: ParticleThroughSideLastPosCheck
#ifdef CODE_ANALYZE
PUBLIC :: OutputTrajectory
#endif /*CODE_ANALYZE*/
PUBLIC :: ComputeXi,ComputeSurfaceDistance2
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


SUBROUTINE IntersectionWithWall(iPart,iLocSide,Element,TriNum)
!===================================================================================================================================
! Compute the Intersection with bilinear surface by approximating the surface with two triangles
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Vars           ,ONLY: lastPartPos
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iPart
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: Node1,Node2
REAL                             :: PoldX, PoldY, PoldZ, nx, ny, nz, nVal
REAL                             :: bx,by,bz,ax,ay,az,dist
REAL                             :: xNod, yNod, zNod
REAL                             :: Vector1(1:3),Vector2(1:3)
!===================================================================================================================================

CNElemID = GetCNElemID(Element)

PoldX = LastPartPos(1,iPart)
PoldY = LastPartPos(2,iPart)
PoldZ = LastPartPos(3,iPart)

xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)

!---- Calculate normal vector:
Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - xNod
Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - yNod
Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - zNod

Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - xNod
Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - yNod
Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - zNod


nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

nVal = SQRT(nx*nx + ny*ny + nz*nz)

nx = nx/nVal
ny = ny/nVal
nz = nz/nVal

!---- Calculate Intersection
bx = PoldX - xNod
by = PoldY - yNod
bz = PoldZ - zNod

ax = bx - nx * (bx * nx + by * ny + bz * nz)
ay = by - ny * (bx * nx + by * ny + bz * nz)
az = bz - nz * (bx * nx + by * ny + bz * nz)

dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
      (az * bx - ax * bz) * (az * bx - ax * bz) +   &
      (ax * by - ay * bx) * (ax * by - ay * bx))/   &
      (ax * ax + ay * ay + az * az))

! If vector from old point to new point goes through the node, a will be zero
! dist is then simply length of vector b instead of |axb|/|a|
IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

TrackInfo%alpha = TrackInfo%PartTrajectory(1) * nx + TrackInfo%PartTrajectory(2) * ny + TrackInfo%PartTrajectory(3) * nz
IF(ABS(TrackInfo%alpha).GT.0.) TrackInfo%alpha = dist / TrackInfo%alpha

RETURN

END SUBROUTINE IntersectionWithWall


PPURE SUBROUTINE ParticleThroughSideCheck3DFast(PartID,iLocSide,Element,ThroughSide,TriNum, IsMortar)
!===================================================================================================================================
!> Routine to check whether a particle crossed the given triangle of a side. The determinant between the normalized trajectory
!> vector and the vectors from two of the three nodes to the old particle position is calculated. If the determinants for the three
!> possible combinations are greater than zero, then the particle went through this triangle of the side.
!> Note that if this is a mortar side, the side of the small neighbouring mortar element has to be checked. Thus, the orientation
!> is reversed.
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars              ,ONLY: EpsMach
USE MOD_Particle_Vars             ,ONLY: lastPartPos
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: ThroughSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3)
!===================================================================================================================================

CNElemID = GetCNElemID(Element)

ThroughSide = .FALSE.

Px = lastPartPos(1,PartID)
Py = lastPartPos(2,PartID)
Pz = lastPartPos(3,PartID)

! Normalized particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = TrackInfo%PartTrajectory(1)
Vy = TrackInfo%PartTrajectory(2)
Vz = TrackInfo%PartTrajectory(3)

! Get the coordinates of the first node and the vector from the particle position to the node
xNode(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNode(1) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNode(1) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)

Ax(1) = xNode(1) - Px
Ay(1) = yNode(1) - Py
Az(1) = zNode(1) - Pz
! Get the vectors to the other two nodes, depending on the triangle number
IF(PRESENT(IsMortar)) THEN
  ! Note: reverse orientation in the mortar case, as the side is treated from the perspective of the smaller neighbouring element
  !       (TriNum=1: NodeID=3,2; TriNum=2: NodeID=4,3)
  xNode(2) = NodeCoords_Shared(1,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  yNode(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  zNode(2) = NodeCoords_Shared(3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)

  Ax(2) = xNode(2) - Px
  Ay(2) = yNode(2) - Py
  Az(2) = zNode(2) - Pz

  xNode(3) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  yNode(3) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  zNode(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)

  Ax(3) = xNode(3) - Px
  Ay(3) = yNode(3) - Py
  Az(3) = zNode(3) - Pz
ELSE
  DO n = 2,3
    NodeID = n+TriNum-1       ! m = true node number of the sides (TriNum=1: NodeID=2,3; TriNum=2: NodeID=3,4)
    xNode(n) = NodeCoords_Shared(1,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    yNode(n) = NodeCoords_Shared(2,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    zNode(n) = NodeCoords_Shared(3,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)

    Ax(n) = xNode(n) - Px
    Ay(n) = yNode(n) - Py
    Az(n) = zNode(n) - Pz
  END DO
END IF
!--- check whether v and the vectors from the particle to the two edge nodes build
!--- a right-hand-system. If yes for all edges: vector goes potentially through side
det(1) = (Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + (Ax(1) * Vy - Ay(1) * Vx) * Az(3)
det(2) = (Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + (Ax(2) * Vy - Ay(2) * Vx) * Az(1)
det(3) = (Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + (Ax(3) * Vy - Ay(3) * Vx) * Az(2)

! Comparison of the determinants with eps due to machine precision
IF ((det(1).ge.-epsMach).AND.(det(2).ge.-epsMach).AND.(det(3).ge.-epsMach)) THEN
  ThroughSide = .TRUE.
END IF

RETURN

END SUBROUTINE ParticleThroughSideCheck3DFast


SUBROUTINE ParticleThroughSideLastPosCheck(i,iLocSide,Element,InElementCheck,TriNum,det,isMortarSide,detPartPos)
!===================================================================================================================================
!> Routine used in the case of a particle crossing multipe sides. Calculates the determinant of the three vectors from the last
!> (and current in the case of mortar sides) particle position to the nodes of the triangle in order to determine if the particle
!> is inside the element. Output of determinant is used to determine which of the sides was crossed first.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Tools         ,         ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: i, Element, iLocSide, TriNum
LOGICAL, INTENT(IN),OPTIONAL     :: isMortarSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: InElementCheck
REAL   ,INTENT(OUT)              :: det
REAL   ,INTENT(OUT), OPTIONAL    :: detPartPos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: NodeNum, ind, iNode
REAL                             :: Ax(3),Ay(3),Az(3)
REAL                             :: NodeCoord(1:3,1:3)
!===================================================================================================================================

CNElemID = GetCNElemID(Element)

InElementCheck = .TRUE.

!--- coords of first node:
DO ind = 1,3
  NodeCoord(ind,1) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
END DO

!--- coords of other two nodes (depending on triangle):
IF(PRESENT(isMortarSide)) THEN
  ! Note: reversed orientation as the triangle is treated from the perspective of the smaller neighbouring mortar element
  NodeCoord(1:3,2) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  NodeCoord(1:3,3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
ELSE
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(NodeNum,iLocSide,CNElemID)+1)
    END DO
  END DO
END IF
!--- vector from lastPos(!) to triangle nodes
DO ind = 1,3
  Ax(ind) = NodeCoord(1,ind) - lastPartPos(1,i)
  Ay(ind) = NodeCoord(2,ind) - lastPartPos(2,i)
  Az(ind) = NodeCoord(3,ind) - lastPartPos(3,i)
END DO

!--- determine whether particle is on inner side (rel. to element) of triangle
!--- set corresponding "flag" (see below)
det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
       (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
       (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

IF ((det.lt.0).OR.(det.NE.det)) THEN
  InElementCheck = .FALSE.
END IF

IF(PRESENT(isMortarSide).AND.PRESENT(detPartPos)) THEN
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - PartState(1,i)
    Ay(ind) = NodeCoord(2,ind) - PartState(2,i)
    Az(ind) = NodeCoord(3,ind) - PartState(3,i)
  END DO

  detPartPos = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
                (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
                (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))
END IF

RETURN

END SUBROUTINE ParticleThroughSideLastPosCheck


SUBROUTINE ComputePlanarRectIntersection(isHit                       &
                                        ,PartTrajectory              &
                                        ,lengthPartTrajectory        &
                                        ,alpha                       &
                                        ,xi                          &
                                        ,eta                         &
                                        ,PartID                       &
                                        ,flip                        &
                                        ,SideID                      &
                                        ,opt_CriticalParallelInSide   )
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,            ONLY:epsMach
USE MOD_Mesh_Tools,              ONLY:GetCNSideID
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,epsilontol,SideDistance
USE MOD_Particle_Surfaces_Vars,  ONLY:BaseVectors0,BaseVectors1,BaseVectors2
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
USE MOD_Mesh_Vars,               ONLY:NGeo
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: PartID,SideID
INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
LOGICAL,INTENT(OUT)               :: isHit
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_CriticalParallelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2
REAL                              :: NormVec(1:3),locDistance,Inter1(1:3), alphaNorm
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,locSideDistance
REAL                              :: sdet
REAL                              :: epsLoc
LOGICAL                           :: CriticalParallelInSide
INTEGER                           :: CNSideID
!INTEGER                           :: flip
!===================================================================================================================================

CNSideID = GetCNSideID(SideID)

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(110("-"))')
      WRITE(UNIT_stdout,'(A,I0)')       '     | Output of planar face constants for Side: ',SideID
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | SideNormVec  : ',SideNormVec(1:3,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint1: ',BezierControlPoints3D(:,0,0,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint2: ',BezierControlPoints3D(:,NGeo,0,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint3: ',BezierControlPoints3D(:,0,NGeo,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint4: ',BezierControlPoints3D(:,NGeo,NGeo,SideID)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! set alpha to minus 1, assume no intersection
alpha    = -1.0
xi       = -2.
eta      = -2.
isHit    = .FALSE.

! new with flip
IF(flip.EQ.0)THEN
  NormVec     =  SideNormVec(1:3,CNSideID)
  locDistance =  SideDistance(CNSideID)
ELSE
  NormVec     = -SideNormVec(1:3,CNSideID)
  locDistance = -SideDistance(CNSideID)
END IF

coeffA=DOT_PRODUCT(NormVec,PartTrajectory)

!! corresponding to particle starting in plane
!! interaction should be computed in last step
CriticalParallelInSide=.FALSE.
IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.

! extension for periodic sides
locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(1:3,PartID),NormVec)

IF (CriticalParallelInSide) THEN ! particle parallel to side
  IF (ALMOSTZERO(locSideDistance)) THEN ! particle on/in side
    IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide=.TRUE.
    ! move particle eps into interior
    alpha=-1.
    RETURN
  END IF
  IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide=.FALSE.
  alpha=-1.
  RETURN
ELSE
  IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide=.FALSE.
  alpha=locSideDistance/coeffA
END IF

IF (locSideDistance.LT.-100*epsMach) THEN
  ! particle is located outside of element, THEREFORE, an intersection were not detected
  alpha = -1. ! here, alpha was set to zero? why?
  isHit = .FALSE.
  RETURN
  ! do I have to compute the xi and eta value? first try: do not re-check new element!
END IF

alphaNorm=alpha/lengthPartTrajectory

! found intersection further than normalized alpha or within negative machine accuracy. Move particle back inside
IF((alphaNorm.GT.1.0) .OR.(alphaNorm.LT.-epsilontol))THEN
  alpha = -1.0
  ishit = .FALSE.
  RETURN
END IF

! Calculate intersection point and initial base vectors
Inter1=LastPartPos(1:3,PartID)+alpha*PartTrajectory
P0 =-0.25*BaseVectors0(:,SideID)+Inter1
P1 = 0.25*BaseVectors1(:,SideID)
P2 = 0.25*BaseVectors2(:,SideID)

A1=P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3)
B1=P2(1)*P1(1)+P2(2)*P1(2)+P2(3)*P1(3)
C1=P1(1)*P0(1)+P1(2)*P0(2)+P1(3)*P0(3)

A2=B1
B2=P2(1)*P2(1)+P2(2)*P2(2)+P2(3)*P2(3)
C2=P2(1)*P0(1)+P2(2)*P0(2)+P2(3)*P0(3)

sdet=A1*B2-A2*B1
IF(ABS(sdet).EQ.0)THEN
  CALL abort(&
  __STAMP__&
  ,' ABS(sdet).EQ.0!')
END IF
sdet=1.0/sdet
epsLoc=1.0+100.*epsMach


xi=(B2*C1-B1*C2)*sdet
! xi outside of reference element, no intersection
IF(ABS(xi).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF

eta=(-A2*C1+A1*C2)*sdet
! eta outside of reference element, no intersection
IF(ABS(eta).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF

! every check passed, the particle crossed the side. alpha, xi and eta are still on correct values
isHit=.TRUE.

END SUBROUTINE ComputePlanarRectIntersection


SUBROUTINE ComputePlanarCurvedIntersection(isHit                        &
                                           ,PartTrajectory              &
                                           ,lengthPartTrajectory        &
                                           ,alpha                       &
                                           ,xi                          &
                                           ,eta                         &
                                           ,PartID                      &
                                           ,flip                        &
                                           ,SideID                      &
                                           ,opt_CriticalParallelInSide)
!===================================================================================================================================
! Compute the intersection with a planar non rectangular face
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Globals                ,ONLY: Cross,abort,UNIT_stdOut,CROSSNORM,UNITVECTOR
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Mesh_Tools             ,ONLY: GetCNSideID
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Particle_Surfaces_Vars ,ONLY: SideNormVec,SideSlabNormals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars ,ONLY: locXi,locEta,locAlpha,SideDistance
USE MOD_Utils                  ,ONLY: InsertionSort
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: rBoundingBoxChecks
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Globals                ,ONLY: myrank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)           :: PartTrajectory
REAL,INTENT(IN)                          :: lengthPartTrajectory
INTEGER,INTENT(IN)                       :: PartID,SideID
INTEGER,INTENT(IN)                       :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                         :: alpha,xi,eta
LOGICAL,INTENT(OUT)                      :: isHit
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParallelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
REAL                                     :: NormVec(1:3),locDistance
INTEGER                                  :: CNSideID,nInterSections,p,q
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
LOGICAL                                  :: CriticalParallelInSide
REAL                                     :: XiNewton(2)
REAL                                     :: coeffA,locSideDistance
! fallback algorithm
LOGICAL                                  :: failed
INTEGER(KIND=2)                          :: ClipMode
REAL                                     :: LineNormVec(1:2,1:2)
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: PartFaceAngle
!===================================================================================================================================

! set alpha to minus 1, assume no intersection
alpha    = -1.0
xi       = 2.0
eta      = 2.0
isHit    = .FALSE.
CNSideID = GetCNSideID(SideID)

#ifdef CODE_ANALYZE
rBoundingBoxChecks = rBoundingBoxChecks+1.
#endif /*CODE_ANALYZE*/

CriticalParallelInSide=.FALSE.

! new with flip
IF(flip.EQ.0)THEN
  NormVec     =  SideNormVec(1:3,CNSideID)
  locDistance =  SideDistance(CNSideID)
ELSE
  NormVec     = -SideNormVec(1:3,CNSideID)
  locDistance = -SideDistance(CNSideID)
END IF

! Calculate distance from particle to planar side face
!> 1) check if particle is moving in other direction or exactly parallel, no intersection
!> 2) difference between SideDistance (distance from origin to sice) and the dot product is the distance of the particle to the side
!> 3) check if distance from particle to side is longer than the particle vector, no intersection
coeffA = DOT_PRODUCT(NormVec,PartTrajectory)
IF (coeffA.LE.0.) RETURN

! difference between SideDistance (distance from origin to side) and the dot product is the distance of the particle to the side
locSideDistance = locDistance-DOT_PRODUCT(LastPartPos(1:3,PartID),NormVec)
locSideDistance = locSideDistance/coeffA
IF (locSideDistance.GT.lengthPartTrajectory) RETURN

IF (TrackingMethod.NE.REFMAPPING) THEN
  IF (ALMOSTZERO(coeffA)) CriticalParallelInSide = .TRUE.
END IF

! Check if the particle intersects the bounding box of the side. If not, we can eliminate the side without doing more checking
IF(.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN

! 2.) Bezier intersection: transformation of bezier patch 3D->2D
!PartTrajectory = PartTrajectoryOrig + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF(ABS(PartTrajectory(3)).LT.0.)THEN
  n1 = (/ -PartTrajectory(2) - PartTrajectory(3),  PartTrajectory(1),  PartTrajectory(1) /)
ELSE
  n1 = (/  PartTrajectory(3),  PartTrajectory(3), -PartTrajectory(1) - PartTrajectory(2) /)
END IF

n1 = UNITVECTOR(n1)
n2 = CROSSNORM(PartTrajectory,n1)

DO q = 0,NGeo
  DO p = 0,NGeo
    BezierControlPoints2D(1,p,q) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n1)
    BezierControlPoints2D(2,p,q) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n2)
  END DO
END DO

XiNewton=0.
CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,PartID,SideID,failed)
! write Xinewton to locXi and locEta
locXi (1) = XiNewton(1)
locEta(1) = XiNewton(2)

IF (failed) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,2,CNSideID))))
  IPWRITE(UNIT_stdout,*) ' Intersection-angle-of-BezierNetwon: ',PartFaceAngle*180./PI

  iClipIter   = 0
  nXiClip     = 0
  nEtaClip    = 0
  nInterSections = 0
  ClipMode    = 1
  LineNormVec = 0.
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                ,iClipIter,nXiClip,nEtaClip,nInterSections,PartID,SideID)
  IF (nInterSections.GT.1) nInterSections = 1
END IF

nInterSections = 0
IF (locAlpha(1).GT.-1) nInterSections = 1

IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide = .FALSE.
IF (CriticalParallelInSide) THEN
  IF (ALMOSTZERO(locAlpha(1))) THEN
    IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide = .TRUE.
  END IF
END IF

SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha = locAlpha(1)
  xi    = locXi (1)
  eta   = loceta(1)
  isHit = .TRUE.
  RETURN
CASE DEFAULT
  CALL abort(__STAMP__,' The code should never go here')
END SELECT


END SUBROUTINE ComputePlanarCurvedIntersection


SUBROUTINE ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                      ,PartID,SideID,ElemCheck_Opt,alpha2)
!===================================================================================================================================
! Compute the Intersection with planar surface, improved version by
! Haselbacher, A.; Najjar, F. M. & Ferry, J. P., An efficient and robust particle-localization algorithm for unstructured grids
! Journal of Computational Physics, Elsevier BV, 2007, 225, 2198-2213
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Utils                  ,ONLY: QuadraticSolver
USE MOD_Mesh_Tools             ,ONLY: GetCNSideID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,SideNormVec,epsilonTol!,BaseVectorsScale
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBilinear
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos!,PEM
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
USE MOD_Mesh_Vars              ,ONLY: NGeo
#endif /*CODE_ANALYZE*/
#if USE_MPI
!USE MOD_Mesh_Vars              ,ONLY: BC
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: PartID,SideID
LOGICAL,INTENT(IN),OPTIONAL       :: ElemCheck_Opt
REAL,INTENT(IN),OPTIONAL          :: alpha2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff,NormalCoeff
REAL                              :: A,B,C,alphaNorm
REAL                              :: xi(2),eta(2),t(2),scaleFac
INTEGER                           :: CNSideID,InterType,nRoot
LOGICAL                           :: ElemCheck
!===================================================================================================================================

! set alpha to minus one // no intersection
alpha    = -1.0
xitild   = -2.0
etatild  = -2.0
isHit    = .FALSE.
CNSideID = GetCNSideID(SideID)

! compute initial vectors
BiLinearCoeff(:,1) = 0.25*BaseVectors3(:,SideID)
BiLinearCoeff(:,2) = 0.25*BaseVectors1(:,SideID)
BiLinearCoeff(:,3) = 0.25*BaseVectors2(:,SideID)
BiLinearCoeff(:,4) = 0.25*BaseVectors0(:,SideID)

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(110("-"))')
      WRITE(UNIT_stdout,'(A)') '     | Output of bilinear intersection equation constants: '
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | SideNormVec  : ',SideNormVec(1:3,CNSideID)
      WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(1,1:4)
      WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(2,1:4)
      WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(3,1:4)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint1: ',BezierControlPoints3D(:,0,0,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint2: ',BezierControlPoints3D(:,NGeo,0,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint3: ',BezierControlPoints3D(:,0,NGeo,SideID)
      WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint4: ',BezierControlPoints3D(:,NGeo,NGeo,SideID)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! Check if the site can be encountered. Both vectors are already normalized
scaleFac = DOT_PRODUCT(PartTrajectory,SideNormVec(1:3,CNSideID))
IF (ABS(scaleFac).LT.epsilontol) RETURN

! Haselbacher et al. define d = d - r_p
BiLinearCoeff(:,4) = BiLinearCoeff(:,4) - LastPartPos(:,PartID)

! Calculate component normal to ray
NormalCoeff(:,1) = BiLinearCoeff(:,1) - SUM(BiLinearCoeff(:,1)*PartTrajectory(:))*PartTrajectory
NormalCoeff(:,2) = BiLinearCoeff(:,2) - SUM(BiLinearCoeff(:,2)*PartTrajectory(:))*PartTrajectory
NormalCoeff(:,3) = BiLinearCoeff(:,3) - SUM(BiLinearCoeff(:,3)*PartTrajectory(:))*PartTrajectory
NormalCoeff(:,4) = BiLinearCoeff(:,4) - SUM(BiLinearCoeff(:,4)*PartTrajectory(:))*PartTrajectory

! A1 is X_xz = X_z - X_x
A1(:) = NormalCoeff(3,:) - NormalCoeff(1,:)
! A2 is X_yz = X_z - X_y
A2(:) = NormalCoeff(3,:) - NormalCoeff(2,:)

! Bring into quadratic form
A = a1(1)*a2(3) - a2(1)*a1(3)
B = a1(1)*a2(4) - a2(1)*a1(4) + a1(2)*a2(3) - a2(2)*a1(3)
C = a1(2)*a2(4) - a2(2)*a1(4)

! Scale with <PartTraj.,NormVec>^2 and cell-scale (~area) for getting coefficients at least approx. in the order of 1
!scaleFac = scaleFac**2 * BaseVectorsScale(SideID) !<...>^2 * cell-scale
!scaleFac = 1./scaleFac
!A = A * scaleFac
!B = B * scaleFac
!C = C * scaleFac

CALL QuadraticSolver(A,B,C,nRoot,Eta(1),Eta(2))

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | Output after QuadraticSolver: '
      WRITE(UNIT_stdout,'(A,I0,A,2(1X,G0))') '     | number of root: ',nRoot,' | Eta: ',Eta(1:2)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! nRoot equals the number of possible intersections with the bilinear surface. However, only values between [-1,1] are valid
SELECT CASE(nRoot)
  ! No intersection
  CASE(0) ! nRoot = 0
    RETURN

  ! One possible intersection
  CASE(1) ! nRoot = 1
#ifdef CODE_ANALYZE
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
      IF(PartID.EQ.PARTOUT)THEN
        WRITE(UNIT_stdout,'(A)') '     | nRoot = 1 '
      END IF
    END IF
#endif /*CODE_ANALYZE*/
    ! Check if eta is valid
    IF (ABS(eta(1)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(1) = ComputeXi(eta(1),A1=A1,A2=A2)

      IF (Xi(1).EQ.HUGE(1.)) THEN
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartID:             ', PartID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPartPos:   ', LastPartPos(1:3,PartID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PartPos:       ', PartState  (1:3,PartID)
        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(1)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(1) = ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory)

        IF (PRESENT(alpha2)) THEN
          IF (alpha2.GT.-1.0 .AND. ALMOSTEQUAL(t(1),alpha2)) THEN
            t(1) = -1.0
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(PartID.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A)') 'changed t1'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
          END IF
        END IF

        ! Normalize alpha to unitLength
        alphaNorm = t(1)/lengthPartTrajectory

        IF ((alphaNorm.LE.1.0) .AND.(alphaNorm.GE.0.)) THEN
          alpha   = t(1)
          xitild  = xi(1)
          etatild = eta(1)
          isHit   = .TRUE.
#ifdef CODE_ANALYZE
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
            IF(PartID.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(A,G0,A,G0)') '     | alphanorm: ',alphaNorm,' | epsilonTolerance: ',epsilontol
            END IF
          END IF
#endif /*CODE_ANALYZE*/
          ! This is the only possible intersection, so we are done
          RETURN
        ELSE ! t is not in range
          RETURN
        END IF
      ELSE ! xi not in range
        RETURN
      END IF ! ABS(xi(1)).LE.1.0
    ELSE ! eta not in range
      RETURN
    END IF ! ABS(eta(1)).LE.1.0

  CASE(2) ! nRoot = 2
#ifdef CODE_ANALYZE
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
      IF(PartID.EQ.PARTOUT)THEN
        WRITE(UNIT_stdout,'(A)') '     | nRoot = 2 '
      END IF
    END IF
#endif /*CODE_ANALYZE*/
    InterType = 0
    t(:)      =-1.

    ! Check if eta(1)) is valid
    IF (ABS(eta(1)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(1) = ComputeXi(eta(1),A1=A1,A2=A2)

      IF (Xi(1).EQ.HUGE(1.)) THEN
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartID:             ', PartID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPartPos:   ', LastPartPos(1:3,PartID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PartPos:       ', PartState  (1:3,PartID)
        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(1)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(1) = ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory)

        IF (PRESENT(alpha2)) THEN
          IF (alpha2.GT.-1.0 .AND. ALMOSTEQUAL(t(1),alpha2)) THEN
#ifdef CODE_ANALYZE
              IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
                IF(PartID.EQ.PARTOUT)THEN
                  WRITE(UNIT_stdout,'(A)') 'changed t1'
                END IF
              END IF
#endif /*CODE_ANALYZE*/
            t(1) = -1.0
          END IF
        END IF

        ! Normalize alpha to unitLength
        alphaNorm = t(1)/lengthPartTrajectory

#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(A,G0,A,G0,A,G0)') '     | xi: ',xi(1),' | t: ',t(1),' | alphaNorm: ',alphaNorm
          END IF
        END IF
#endif /*CODE_ANALYZE*/

        IF ((alphaNorm.LE.1.0) .AND.(alphaNorm.GE.0.)) THEN
          InterType = InterType+1
          isHit     = .TRUE.
#ifdef CODE_ANALYZE
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
            IF(PartID.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | alphanorm1: ',alphaNorm,' | epsilonTolerance: ',epsilontol
            END IF
          END IF
#endif /*CODE_ANALYZE*/
        END IF
      END IF ! ABS(xi(1)).LE.1.0
    END IF ! ABS(eta(1)).LE.1.0

    ! Check if eta(2) is valid
    IF (ABS(eta(2)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(2) = ComputeXi(eta(2),A1=A1,A2=A2)

      IF (Xi(2).EQ.HUGE(1.)) THEN
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' PartID:             ', PartID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPartPos:   ', LastPartPos(1:3,PartID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PartPos:       ', PartState  (1:3,PartID)
        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(2)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(2) = ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory)

        IF (PRESENT(alpha2)) THEN
          IF (alpha2.GT.-1.0 .AND. ALMOSTEQUAL(t(2),alpha2)) THEN
            t(2) = -1.0
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(PartID.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A)') 'changed t2'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
          END IF
        END IF

        ! Normalize alpha to unitLength
        alphaNorm = t(2)/lengthPartTrajectory

#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(A,G0,A,G0,A,G0)') '     | xi: ',xi(2),' | t: ',t(2),' | alphaNorm: ',alphaNorm
          END IF
        END IF
#endif /*CODE_ANALYZE*/

        IF ((alphaNorm.LE.1.0) .AND.(alphaNorm.GE.0.)) THEN
          InterType = InterType+2
          isHit     = .TRUE.
#ifdef CODE_ANALYZE
          IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
            IF(PartID.EQ.PARTOUT)THEN
              WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | alphanorm2: ',alphaNorm,' | epsilonTolerance: ',epsilontol
            END IF
          END IF
#endif /*CODE_ANALYZE*/
        END IF
      END IF ! ABS(xi(2)).LE.1.0
    END IF ! ABS(eta(2)).LE.1.0

    SELECT CASE(InterType)
      ! No intersection found, return
      CASE(0)
        RETURN

      ! First intersection is only hit
      CASE(1)
        alpha  =t  (1)
        xitild =xi (1)
        etatild=eta(1)

      ! Second intersection is only hit
      CASE(2)
        alpha  =t  (2)
        xitild =xi (2)
        etatild=eta(2)

      ! Two intersections found, decide on the correct one
      CASE(3)
        ! If side is a BC side, take only the intersection encountered first
        IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
          SELECT CASE(TrackingMethod)
            ! Take the one encountered first
            CASE(REFMAPPING)
              IF(t(1).LT.t(2))THEN
                alpha  =t  (1)
                xitild =xi (1)
                etatild=eta(1)
              ELSE
                alpha  =t  (2)
                xitild =xi (2)
                etatild=eta(2)
              END IF

            CASE(TRACING)
              ! Check if the element is supposed to be checked
              ElemCheck = .FALSE.
              IF(PRESENT(ElemCheck_Opt))THEN
                ElemCheck = ElemCheck_Opt
              END IF

              IF(ElemCheck)THEN
                alpha  =-1
                xitild =-2
                etatild=-2
              ELSE
                ! Apparently we don't care about the direction of the PartTrajectory
                IF(ABS(t(1)).LT.ABS(t(2)))THEN
                  alpha  =t  (1)
                  xitild =xi (1)
                  etatild=eta(1)
                ELSE
                  alpha  =t  (2)
                  xitild =xi (2)
                  etatild=eta(2)
                END IF
              END IF
          END SELECT ! TrackingMethod
        ! Inner side with double intersection, particle leaves and enters element
        ELSE
          alpha  =-1
          xitild = 0.
          etatild= 0.
          isHit  = .FALSE.
        END IF
    END SELECT ! InterType
END SELECT ! nRoot

END SUBROUTINE ComputeBiLinearIntersection


SUBROUTINE ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID &
                ,flip,SideID,opt_CriticalParallelInSide,ElemCheck_Opt)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals,                 ONLY:Cross,abort,CROSSNORM,UNITVECTOR
USE MOD_Globals_Vars,            ONLY:PI
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Mesh_Tools,              ONLY:GetCNSideID
USE MOD_Particle_Mesh_Vars,      ONLY:SideInfo_Shared
USE MOD_Particle_Surfaces,       ONLY:CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,BezierNewtonAngle
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,  ONLY:locXi,locEta,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipLocalTol
USE MOD_Particle_Tracking_Vars,  ONLY:TrackingMethod
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Utils,                   ONLY:InsertionSort
#ifdef CODE_ANALYZE
USE MOD_Globals,                 ONLY:MyRank,UNIT_stdOut
USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints
USE MOD_Particle_Surfaces_Vars,  ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)           :: PartTrajectory
REAL,INTENT(IN)                          :: lengthPartTrajectory
INTEGER,INTENT(IN)                       :: PartID,SideID
INTEGER,INTENT(IN)                       :: flip
LOGICAL,INTENT(IN),OPTIONAL              :: ElemCheck_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                         :: alpha,xi,eta
LOGICAL,INTENT(OUT)                      :: isHit
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParallelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3),NormVec(1:3)
INTEGER                                  :: CNSideID,nInterSections,iInter,p,q
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
#ifdef CODE_ANALYZE
REAL                                     :: BezierControlPoints2D_tmp(2,0:NGeo,0:NGeo)
#endif /*CODE_ANALYZE*/
INTEGER,ALLOCATABLE,DIMENSION(:)         :: locID,realInterID
INTEGER(KIND=2)                          :: ClipMode
REAL                                     :: LineNormVec(1:2,1:2)
INTEGER                                  :: realnInter,isInter
REAL                                     :: XiNewton(2)
REAL                                     :: PartFaceAngle,dXi,dEta
LOGICAL                                  :: CriticalParallelInSide,failed
!REAL                                     :: Interval1D,dInterVal1D
!===================================================================================================================================
! set alpha to minus 1, asume no intersection
alpha    = -1.0
Xi       = 2.0
Eta      = 2.0
isHit    = .FALSE.
CNSideID = GetCNSideID(SideID)

#ifdef CODE_ANALYZE
rBoundingBoxChecks = rBoundingBoxChecks + 1.
#endif /*CODE_ANALYZE*/

! new with flip
IF(flip.EQ.0)THEN
  NormVec     =  SideNormVec(1:3,CNSideID)
ELSE
  NormVec     = -SideNormVec(1:3,CNSideID)
END IF

CriticalParallelInSide = .FALSE.
IF (BoundingBoxIsEmpty(CNSideID)) THEN
  IF (TrackingMethod.EQ.REFMAPPING) THEN
    IF (DOT_PRODUCT(NormVec,PartTrajectory).LT.0.) RETURN
  ELSE
    IF (ALMOSTZERO(DOT_PRODUCT(NormVec,PartTrajectory))) CriticalParallelInSide=.TRUE.
  END IF
  IF (.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the
                                                                                                  ! bounding box
ELSE
  ! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
  IF (.NOT.InsideBoundingBox(LastPartPos(1:3,PartID),SideID)) THEN ! the old particle position is not inside the bounding box
    IF (.NOT.InsideBoundingBox(PartState(1:3,PartID),SideID)) THEN ! the new particle position is not inside the bounding box
      IF (.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the
                                                                                                  ! bounding box
    END IF
  END IF
END IF

! 2.) Bezier intersection: transformation of bezier patch 3D->2D
! PartTrajectoryOrig = PartTrajectory
! PartTrajectory     = PartTrajectory + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF (ABS(PartTrajectory(3)).LT.0.) THEN
  n1=(/-PartTrajectory(2) - PartTrajectory(3), PartTrajectory(1),  PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3),  PartTrajectory(3),-PartTrajectory(1) - PartTrajectory(2) /)
END IF

! check angle to boundingbox (height normal vector)
PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,2,CNSideID))))
IF (ALMOSTZERO(PartFaceAngle*180/ACOS(-1.))) THEN
  n1 = n1 ! +epsilontol
END IF

! n1=n1/SQRT(DOT_PRODUCT(n1,n1))
! n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
!        , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
!        , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
! n2=n2/SQRT(DOT_PRODUCT(n2,n2))
n1 = UNITVECTOR(n1)
n2 = CROSSNORM(PartTrajectory,n1)
!PartTrajectory = PartTrajectoryOrig !set back for preventing angles > 90 deg (0.5pi+eps)

! projection like Nishita
! plane 1 with n1 becomes y-axis and plane 2 with n2 becomes the x-axis
DO q = 0,NGeo
  DO p = 0,NGeo
    ! n2 is perpendicular to x-axis => gives distance to new x-axis
    BezierControlPoints2D(1,p,q) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID) - LastPartPos(1:3,PartID),n2)
    ! n1 is perpendicular to y-axis => gives distance to new y-axis
    BezierControlPoints2D(2,p,q) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID) - LastPartPos(1:3,PartID),n1)
  END DO
END DO

! calculate angle between particle path and slab normal plane of face
IF ((PartFaceAngle.LT.BezierNewtonAngle)) THEN ! 1° = 0.01745rad: critical side at the moment need: 0.57° angle
#ifdef CODE_ANALYZE
rPerformBezierClip = rPerformBezierClip+1.
#endif /*CODE_ANALYZE*/
  !  this part in a new function or subroutine
  locAlpha  = -1.0
  iClipIter = 0
  nXiClip   = 0
  nEtaClip  = 0
  nInterSections = 0
  ! check extend in Xi  and eta direction
  dXi  = MAXVAL(BezierControlPoints2D(1,:,:)) - MINVAL(BezierControlPoints2D(1,:,:))
  dEta = MAXVAL(BezierControlPoints2D(2,:,:)) - MINVAL(BezierControlPoints2D(2,:,:))
  IF(dXi.GT.dEta)THEN
    ! first Clip is in XI diretion
    ClipMode = 1
  ELSE
    ! first Clip is in ETA diretion
    ClipMode = 2
  END IF
  BezierClipLocalTol = MIN(dXi,dEta)*BezierClipTolerance
  LineNormVec = 0.
  ! CALL recursive Bezier clipping algorithm
#ifdef CODE_ANALYZE
  IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
    IF (PartID.EQ.PARTOUT) THEN
      IPWRITE(UNIT_stdout,*) ' --------------------------------------------- '
      IPWRITE(UNIT_stdout,*) ' clipping '
      IPWRITE(UNIT_stdout,*) ' BezierClipTolerance   ', BezierClipTolerance
      IPWRITE(UNIT_stdout,*) ' BezierClipLocalTol    ', BezierClipLocalTol
      IPWRITE(UNIT_stdout,*) ' ClipMode    ', ClipMode
      IPWRITE(UNIT_stdout,*) ' n1    ', n1
      IPWRITE(UNIT_stdout,*) ' n2    ', n2
      IPWRITE(UNIT_stdout,*) ' dXi,dEta    ', dXi,dEta
      IPWRITE(UNIT_stdout,*) ' BezierControlpoints3D '
      CALL OutputBezierControlPoints(BezierControlPoints3D_in=BezierControlPoints3D(:,:,:,SideID))
    END IF
  END IF
#endif /*CODE_ANALYZE*/
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory&
                ,iClipIter,nXiClip,nEtaClip,nInterSections,PartID,SideID)

ELSE !BezierNewtonAngle
#ifdef CODE_ANALYZE
  rPerformBezierNewton      = rPerformBezierNewton+1.
  BezierControlPoints2D_tmp = BezierControlPoints2D
  locAlpha    = -1.0
  iClipIter   = 0
  nXiClip     = 0
  nEtaClip    = 0
  nInterSections = 0
  ClipMode    = 1
  LineNormVec = 0.
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_tmp,LineNormVec,PartTrajectory,lengthPartTrajectory&
                ,iClipIter,nXiClip,nEtaClip,nInterSections,PartID,SideID)
  IF (nInterSections.GT.1) &
    CALL abort(__STAMP__,' More then one intersection! Cannot use Newton!' ,nInterSections)
#endif /*CODE_ANALYZE*/

  XiNewton = 0.
  CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,PartID,SideID,failed)
  nInterSections = 0
  IF (locAlpha(1).GT.-1) nInterSections=1
  IF (failed) CALL abort(__STAMP__,' Bezier-Newton does not yield root! ')

#ifdef CODE_ANALYZE
  IF(nInterSections.EQ.1)THEN
    dXi  = ABS(locXi(1)-XiNewton(1)) !/(400.*BezierClipTolerance)
    dEta = ABS(locEta(1)-XiNewton(2))!/(400.*BezierClipTolerance)
    dXi  = dXi*dXi+dEta*dEta
    dXi  = SQRT(dXi)/(400.*BezierClipTolerance)
    IF (dXi.GT.1.0) THEN
      IPWRITE(UNIT_stdout,*) ': Difference between Intersections > Tolerance'
      IPWRITE(UNIT_stdout,*) ': xi-clip,   xi-newton', locXi(1), XiNewton(1)
      IPWRITE(UNIT_stdout,*) ': eta-clip, eta-newton', loceta(1), XiNewton(2)
      !CALL abort(__STAMP__ &
      ! ' Wrong intersection in Xi! Clip/Newton=',nInterSections,dXi)
    END IF
    !IF(dXi.GT.1.0)THEN
    !  IPWRITE(UNIT_stdout,*) ' eta-clip, eta-newton', loceta(1), XiNewton(2)
    !  CALL abort(__STAMP__ &
    !   ' Wrong intersection in Eta! Clip/Newton=',nInterSections, dXi)
    !END IF
  END IF
#endif /*CODE_ANALYZE*/

  locXi (1) = XiNewton(1)
  locEta(1) = XiNewton(2)
END IF

IF (PRESENT(ElemCheck_Opt)) THEN
  IF (ElemCheck_Opt) THEN
    DO iInter = 1,nInterSections
      CALL CalcNormAndTangBezier(nVec=n1,xi=locXi(iInter),eta=locEta(iInter),SideID=SideID)
      IF (ABS(DOT_PRODUCT(PartTrajectory,n1)).LT.0.125) THEN
        locAlpha(iInter) = -1
        nInterSections   = nInterSections-1
      END IF
    END DO ! iInter=2,nInterSections
  END IF
END IF

#ifdef CODE_ANALYZE
IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
  IF (PartID.EQ.PARTOUT) THEN
    IPWRITE(UNIT_stdout,*)'----------------------------------------------'
    IPWRITE(UNIT_stdout,*)' PARTOUT        = ',PARTOUT
    IPWRITE(UNIT_stdout,*)' nInterSections = ',nInterSections
  END IF
END IF
#endif /*CODE_ANALYZE*/

SELECT CASE(nInterSections)
  CASE(0)
    RETURN

  CASE(1)
    alpha = locAlpha(1)
    xi    = locXi (1)
    eta   = loceta(1)
    isHit = .TRUE.
    IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide=.FALSE.
    IF (CriticalParallelInSide) THEN
      IF (ALMOSTZERO(alpha)) THEN
        IF (PRESENT(opt_CriticalParallelInSide)) opt_CriticalParallelInSide=.TRUE.
      END IF
    END IF
    RETURN

  CASE DEFAULT
    ! more than one intersection
    ALLOCATE(locID(nInterSections))
    DO iInter = 1,nInterSections
      locID(iInter) = iInter
    END DO ! iInter
    ! sort intersection distance
    CALL InsertionSort(locAlpha(1:nIntersections),locID,nIntersections)
#ifdef CODE_ANALYZE
    IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
      IF (PartID.EQ.PARTOUT) THEN
        IPWRITE(UNIT_stdout,*) ' locAlpha-sorted ',locAlpha(1:nIntersections)
      END IF
    END IF
#endif /*CODE_ANALYZE*/
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      DO iInter = 1,nInterSections
        IF (locAlpha(iInter).GT.-1.0) THEN
          alpha = locAlpha(iInter)
          xi    = locXi (locID(iInter))
          eta   = loceta(locID(iInter))
          DEALLOCATE(locID)
          isHit = .TRUE.
          RETURN
        END IF
      END DO ! iInter
    ELSE
      ! no ref mapping
      ! get real number of intersections
      realnInter  = 1
      ALLOCATE(realInterID(1:nInterSections))
      realInterID = 0
      realInterID(1) = 1
      ! PO & CS:
      ! we used the approach to check the previous (i-1)  with the current (i) alpha, if they
      ! are almost identically, it is ignored (multiple intersections are reduced to one)
      ! second possibility:
      ! check only to the accepted alphas
      DO iInter = 2,nInterSections
        IF (.NOT.ALMOSTEQUALRELATIVE(locAlpha(iInter-1),locAlpha(iInter),0.002)) THEN
          realNInter = realNInter + 1
          realInterID(realNInter) = iInter
          isInter   = iInter
        END IF
      END DO ! iInter=2,nInterSections
#ifdef CODE_ANALYZE
       IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN
         IF (PartID.EQ.PARTOUT) THEN
           IPWRITE(UNIT_stdout,*) ' realnInter ',realnInter
         END IF
       END IF
#endif /*CODE_ANALYZE*/
      IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
        IF (PRESENT(ElemCheck_Opt)) THEN
          IF (ElemCheck_Opt) THEN
            IF (MOD(realNInter,2).EQ.0) THEN
              alpha = -1
              nInterSections = 0
              RETURN
            END IF
          END IF
        END IF

        ! boundary side, take first intersection
        alpha = locAlpha(1)
        xi    = locXi (locID(1))
        eta   = loceta(locID(1))
        DEALLOCATE(locID)
        DEALLOCATE(realInterID)
        isHit = .TRUE.
        IF(PRESENT(opt_CriticalParallelInSide)) THEN
          opt_CriticalParallelInSide = .FALSE.
          IF(CriticalParallelInSide) THEN
            IF(ALMOSTZERO(alpha)) THEN
              opt_CriticalParallelInSide=.TRUE.
            END IF
          END IF
        END IF
        RETURN

      ELSE
        IF (MOD(realNInter,2).EQ.0) THEN
          ! particle leaves and enters cell multiple times, however, remain
          ! still inside of the element
          DEALLOCATE(locID)
          DEALLOCATE(realInterID)
          alpha = -1.0
          isHit = .FALSE.
          IF (PRESENT(opt_CriticalParallelInSide)) THEN
            opt_CriticalParallelInSide = .FALSE.
            IF (CriticalParallelInSide) THEN
              IF (ALMOSTZERO(alpha)) THEN
                opt_CriticalParallelInSide = .TRUE.
              END IF
            END IF
          END IF
          RETURN

        ELSE
          ! particle leaves and enters, take the LAST intersection
          alpha = locAlpha(realInterID(realNInter))
          xi    = locXi (locID(realInterID(realNInter)))
          eta   = loceta(locID(realInterID(realNInter)))
          isHit = .TRUE.
          DEALLOCATE(locID)
          DEALLOCATE(realInterID)
          IF(PRESENT(opt_CriticalParallelInSide)) THEN
            opt_CriticalParallelInSide=.FALSE.
            IF(CriticalParallelInSide)THEN
              IF(ALMOSTZERO(alpha))THEN
                opt_CriticalParallelInSide=.TRUE.
              END IF
            END IF
          END IF
          RETURN
        END IF
      END IF ! BC or no BC side
    END IF
    SDEALLOCATE(locID)
END SELECT

CALL abort(__STAMP__,' The code should never go here')

! remove compiler warning
IF (PRESENT(ElemCheck_Opt)) ClipMode = 0

END SUBROUTINE ComputeCurvedIntersection


RECURSIVE SUBROUTINE BezierClipRecursive(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                               ,iClipIter,nXiClip,nEtaClip&
                               ,nInterSections,PartID,SideID)
!================================================================================================================================
! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
! original article:
!   author = {Nishita, Tomoyuki and Sederberg, Thomas W. and Kakimoto, Masanori},
!   title = {Ray Tracing Trimmed Rational Surface Patches},
!   year = {1990},
! book:
!   author = {Farin, Gerald},
!   title = {Curves and Surfaces for CAGD: A Practical Guide},
!   year = {2002},
!================================================================================================================================
USE MOD_Globals,                 ONLY:Abort,UNIT_stdOut
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipMaxIter
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipLineVectorMethod
#ifdef CODE_ANALYZE
USE MOD_Globals,                 ONLY:MyRank
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: lengthPartTrajectory
REAL,INTENT(INOUT)                   :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: SideID,PartID
REAL,INTENT(IN),DIMENSION(1:3)       :: PartTrajectory
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                  :: iClipIter,nXiClip,nEtaClip,nInterSections
INTEGER(KIND=2),INTENT(INOUT)          :: ClipMode
REAL,DIMENSION(2,2),INTENT(INOUT)      :: LineNormVec
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: PatchDOF2D
!================================================================================================================================

PatchDOF2D=1.0/REAL((NGeo+1)*(NGeo+1))

! 3.) Bezier intersection: solution Newton's method or Bezier clipping
! outcome: no intersection, single intersection, multiple intersection with patch
DO WHILE(iClipIter.LE.BezierClipMaxIter)
  IF(iClipIter.EQ.0)THEN
    IF(BezierClipLineVectorMethod.EQ.0) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
  END IF
  iClipIter=iClipIter+1
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A,I0,1X,I0)') ' iClipIter,ClipMode ', iClipIter, ClipMode
    END IF
  END IF
#endif /*CODE_ANALYZE*/
  SELECT CASE(ClipMode)
  CASE(-1)
    ! no intersection possible
    RETURN
  CASE(1)
    ! LineNormVec is only computed, if a Xi and Eta Clip is performed.
    ! we compute LineNormVecs only until one direction is converged, than we keep the vector to report the correct
    ! results, see. Efremov 2005
    IF(BezierClipLineVectorMethod.EQ.1) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
    IF(BezierClipLineVectorMethod.EQ.2) THEN
      !IF(MOD(iClipIter,2).EQ.1) CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec(:,:),NGeo,0)
      IF(MOD(iClipIter,2).EQ.1) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
    END IF
    CALL CheckXiClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                       ,iClipIter,nXiClip,nEtaClip&
                       ,nInterSections,PartID,SideID)
  CASE(2)
    ! LineNormVec is only computed, if a Xi and Eta Clip is performed.
    ! we compute LineNormVecs only until one direction is converged, than we keep the vector to report the correct
    ! results, see. Efremov 2005
    IF(BezierClipLineVectorMethod.EQ.1) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
    IF(BezierClipLineVectorMethod.EQ.2) THEN
      IF(MOD(iClipIter,2).EQ.1) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
    END IF
    CALL CheckEtaClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                       ,iClipIter,nXiClip,nEtaClip&
                       ,nInterSections,PartID,SideID)
  CASE(3)
    CALL CheckXiClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                       ,iClipIter,nXiClip,nEtaClip&
                       ,nInterSections,PartID,SideID)
  CASE(4)
    CALL CheckEtaClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                       ,iClipIter,nXiClip,nEtaClip&
                       ,nInterSections,PartID,SideID)
  CASE(5)
    ! validate found intersection
    CALL ComputeBezierIntersectionPoint(nXiClip,nEtaClip,PartID,SideID,nInterSections,PartTrajectory,lengthPartTrajectory)
    RETURN ! leave, because convergence
  CASE DEFAULT

    CALL abort(&
__STAMP__ &
      ,' ClipMode is defined in [1-5]! ')
  END SELECT

END DO ! iClipIter=iClipIter,BezierClipMaxIter

IF(iClipIter.GE.BezierClipMaxIter)THEN
  WRITE(UNIT_stdout,'(A,I0)') 'Iter   ',iClipIter
  WRITE(UNIT_stdout,'(A,I0)') 'PartID ',PartID
  WRITE(UNIT_stdout,'(A)') 'Bezier Clipping not converged!'
  STOP
END IF

END SUBROUTINE BezierClipRecursive


SUBROUTINE BezierNewton(alpha,Xi,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,PartID,SideID,failed)
!===================================================================================================================================
! Newton to find root in projected plane for curved tracking
! whole Newton operates in [-1,1]
! output: [-1,1]
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierNewtonTolerance2,BezierNewtonMaxIter,BezierControlPoints3D,epsilontol &
                                     ,BezierNewtonGuess
USE MOD_Particle_Vars,           ONLY:LastPartPos!,PartState
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars,  ONLY:D_Bezier
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierNewtonHit
USE MOD_Globals,                 ONLY:MyRank,UNIT_stdOut
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)    :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
REAL,INTENT(IN)    :: PartTrajectory(3),lengthPartTrajectory
INTEGER,INTENT(IN) :: PartID
INTEGER,INTENT(IN) :: SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Xi(2)
REAL,INTENT(OUT)   :: alpha
LOGICAL,INTENT(OUT):: failed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: dXi(2),sdet,dXi2,alphaNorm
REAL               :: dBezierControlPoints2D(2,2,0:NGeo,0:NGeo)
REAL               :: P(2),gradXi(2,2),InterP(3)
REAL               :: IntersectionVector(3)
INTEGER            :: nIter
INTEGER            :: dd,nn,l,i,j
INTEGER            :: CycIJ(2),Cyc(2)
INTEGER            :: iArmijo
REAL               :: lambda,Xi_old(2)
REAL               :: Norm_P, Norm_P_old
LOGICAL            :: hasInter
REAL               :: MinMax(1:2,1:2)
!===================================================================================================================================

failed=.FALSE.
! check if intersection is possible
hasInter=.TRUE.
DO l=1,2
  MinMax(1,l)=MINVAL(BezierControlPoints2D(l,:,:))
  MinMax(2,l)=MAXVAL(BezierControlPoints2D(l,:,:))
  IF(MinMax(1,l)*MinMax(2,l).GT.0.) hasInter=.FALSE.
END DO ! l=1,2

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,*) ' Can intersect? ', hasInter
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF(.NOT.hasInter) THEN
  alpha=-1.
  RETURN
END IF

! inital guess
SELECT CASE(BezierNewtonGuess)
CASE(1)
  ! assume:_minvalue at -1, maxvalue at 1
  DO l=1,2
    Xi(l) = -1.0 - 2*MinMax(1,l) / ( MinMax(2,l)-MinMax(1,l) )
  END DO ! l=1,2
CASE(2)
  ! take the absolute value of control points to get init-value
  Norm_P=SQRT(DOT_PRODUCT(BezierControlPoints2D(:,0,0),BezierControlPoints2D(:,0,0)))
  Norm_P_old=Norm_P
  Xi(:) = (/0.,0./)
  DO j=0,NGeo
    DO i=0,NGeo
      dXi(1) = ABS(BezierControlPoints2D(1,i,j))
      IF(dXi(1).GT.Norm_P_old) CYCLE
      dXi(2) = ABS(BezierControlPoints2D(2,i,j))
      IF(dXi(2).GT.Norm_P_old) CYCLE
      Norm_P=SQRT(dXi(1)*dXi(1)+dXi(2)*dXi(2))
      IF(Norm_P.LT.Norm_P_Old)THEN
        Norm_P_old=Norm_P
        Xi(:) = (/REAL(i),REAL(j)/)
      END IF
    END DO ! i=0,NGeo
  END DO ! j=0,NGeo
  ! compute actual position
  DO l=1,2
    Xi(l) = 2./REAL(NGeo) * Xi(l) -1.
  END DO ! l=1,2
CASE(4)
  ! trival guess
  Xi(:)=(/0.,0./)
CASE DEFAULT
  CALL abort(&
__STAMP__ &
    ,' BezierNewtonGuess is not implemented! ', BezierNewtonGuess)
END SELECT

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,*) ' Initial-guess in BezierNewton ', Xi
  END IF
END IF
#endif /*CODE_ANALYZE*/

nIter=1
alpha=-1.0
dXi2=1.0
! compute gradient at each control point and use D-matrix
dBezierControlPoints2D=0.
DO nn=1,2
  DO dd=1,2 !iSize
    DO j=0,NGeo
      CycIJ(2)=j
      DO i=0,NGeo
        CycIJ(1)=i
        ! Matrix-vector multiplication
        Cyc=CycIJ
        DO l=0,NGeo
          Cyc(dd)=l  ! d/dxi_dd
          dBezierControlPoints2D(dd,nn,i,j) = dBezierControlPoints2D(dd,nn,i,j) + &
                                              D_Bezier(CycIJ(dd),l)*BezierControlPoints2D(nn,Cyc(1),Cyc(2))
        END DO ! l=0,NGeo
      END DO ! j=0,NGeo
    END DO ! i=0,NGeo
  END DO ! dd=1,2
END DO ! nn=1,2


! compute f(xi) and df(xi)/dxi
CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,2, BezierControlPoints2D(1:2,    0:NGeo,0:NGeo) &
                                                  ,dBezierControlPoints2D(1:2,1:2,0:NGeo,0:NGeo),Point=P,Gradient=gradXi)
! and norm
Norm_P=P(1)*P(1)+P(2)*P(2)

DO WHILE((dXi2.GT.BezierNewtonTolerance2).AND.(nIter.LE.BezierNewtonMaxIter))

  ! caution with index
  sdet=gradXi(1,1)*gradXi(2,2)-gradXi(1,2)*gradXi(2,1)
  !CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,2, BezierControlPoints2D(1:2,    0:NGeo,0:NGeo) &
  !                                                  , Point=P,Gradient=gradXi)
  ! caution with index

  IF(ABS(sdet).GT.epsilon(0.)) THEN
    sdet=1./sdet
  ELSE !shit
    alpha=-1.0
    Xi=1.5
    EXIT
    !CALL abort(__STAMP__ &
    !   'Bezier-Netwton singular. iter,sdetJac',nIter,sDet)
  END IF

  ! build 2x2 inverse and multiply by vector
  dXi(1) = gradXi(2,2)*P(1)-gradXi(2,1)*P(2)
  dXi(2) =-gradXi(1,2)*P(1)+gradXi(1,1)*P(2)

  dXi    = sdet*dXi
  dXi2   = dXi(1)*dXi(1)+dXi(2)*dXi(2)
  Xi_old = Xi

  ! Armijo method to enforce global convergence
  lambda=1.0
  iArmijo=1
  Norm_P_old=Norm_P
  Norm_P=Norm_P*2.
  DO WHILE(Norm_P.GT.Norm_P_old*(1.0-0.0001*lambda) .AND.iArmijo.LE.6)
    ! update to new position
    Xi = Xi_old - lambda*dXi
    ! if a particle hit the edge, then the solution process can produce an overshoot
    ! currently, the overshoot is only accounted in the first iteration
    CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,2, BezierControlPoints2D(1:2,    0:NGeo,0:NGeo) &
                                                      , Point=P,Gradient=gradXi)
    ! compute Norm_P
    Norm_P=P(1)*P(1)+P(2)*P(2)
    lambda=0.3*lambda
    iArmijo=iArmijo+1
  END DO
  IF((nIter.GE.4).AND.(ANY(ABS(Xi).GT.1.5))) THEN
    ! no intersection of ray and bezier patch
    Xi=1.5
    EXIT
  END IF
  nIter=nIter+1
END DO

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,*) ' Newton converget to ', Xi
    IPWRITE(UNIT_stdout,*) ' Tolarance Vaule ', BezierNewtonHit
    IPWRITE(UNIT_stdout,*) ' Check if it is a zero: ',P
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF(nIter.GT.BezierNewtonMaxIter) THEN
  IPWRITE(UNIT_stdout,*) ' WARNING: Bezier-Newton not converged!'
  failed=.TRUE.
  RETURN
END IF

! check if found Xi,Eta are in parameter range
IF(ABS(xi(1)).GT.1.002) RETURN
IF(ABS(xi(2)).GT.1.002) RETURN

! compute 3D intersection
CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=InterP)
IntersectionVector=InterP-LastPartPos(1:3,PartID)

alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)

alphaNorm=alpha/lengthPartTrajectory

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,*) ' Intersection Point ', InterP
    IPWRITE(UNIT_stdout,*) ' alphanorm  ',alphaNorm
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF((alphaNorm.LE.1.0).AND.(alphaNorm.GT.-epsilontol)) RETURN
alpha=-1.0

END SUBROUTINE BezierNewton


FUNCTION InsideBoundingBox(ParticlePosition,SideID)
!================================================================================================================================
! check is the particles is inside the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Mesh_Tools,              ONLY: GetCNSideID
USE MOD_Particle_Surfaces_Vars,  ONLY: SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: ParticlePosition
INTEGER,INTENT(IN)                   :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: InsideBoundingBox
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: CNSideID
REAL                                 :: x,y,z,P(3)
!================================================================================================================================

CNSideID = GetCNSideID(SideID)
P = ParticlePosition - BezierControlPoints3D(1:3,0,0,SideID)

! y is perpendicular to xi & eta directions --> check first, smallest interval
y = DOT_PRODUCT(P,SideSlabNormals(:,2,CNSideID))
IF ((y.LT.SideSlabIntervals(3,CNSideID)-100.*epsMach).OR.(y.GT.SideSlabIntervals(4,CNSideID)+100.*epsMach)) THEN
  InsideBoundingBox = .FALSE.
  RETURN
END IF
! than xi
x = DOT_PRODUCT(P,SideSlabNormals(:,1,CNSideID))
IF ((x.LT.SideSlabIntervals(1,CNSideID)-100.*epsMach).OR.(x.GT.SideSlabIntervals(2,CNSideID)+100.*epsMach)) THEN
  InsideBoundingBox = .FALSE.
  RETURN
END IF
! than eta
z = DOT_PRODUCT(P,SideSlabNormals(:,3,CNSideID))
IF ((z.LT.SideSlabIntervals(5,CNSideID)-100.*epsMach).OR.(z.GT.SideSlabIntervals(6,CNSideID)+100.*epsMach)) THEN
  InsideBoundingBox = .FALSE.
  RETURN
END IF

InsideBoundingBox = .TRUE.

END FUNCTION InsideBoundingBox


FUNCTION BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals,                 ONLY: ABORT
USE MOD_Globals_Vars
USE MOD_Mesh_Tools,              ONLY: GetCNSideID
USE MOD_Particle_Surfaces_Vars,  ONLY: SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
USE MOD_Particle_Vars,           ONLY: LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: lengthPartTrajectory
INTEGER,INTENT(IN)                   :: PartID,SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: BoundingBoxIntersection
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: dnk,alpha(2,3)
REAL                                 :: maxvalue,minvalue
INTEGER                              :: CNSideID
INTEGER                              :: i
!================================================================================================================================

CNSideID = GetCNSideID(SideID)

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
DO i=1,3!x,y,z direction
  dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,CNSideID))

  !IF(ABS(dnk).LT.epsilontol)THEN
  IF(ABS(dnk).LT.100.*epsMach)THEN
    dnk=100.*epsMach ! ÜBERPRÜFEN OB SIGN sinn macht
  END IF

  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
  END IF
END DO!i
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAXVAL(alpha(1,:)) ! taken the maxvalue of the minima
minvalue=MINVAL(alpha(2,:)) ! taken the minvalue of the maxima

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
    IF((maxvalue.LT.lengthPartTrajectory+100*epsMach).AND.(maxvalue+100*epsMach.GT.0.))THEN
      !the first intersection is less than lengthPartTrajectory and greater 0
      BoundingBoxIntersection=.TRUE.
    ELSE
      BoundingBoxIntersection=.FALSE.
    END IF
ELSE
  BoundingBoxIntersection=.FALSE.
END IF
END FUNCTION BoundingBoxIntersection


FUNCTION FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Mesh_Tools,              ONLY: GetCNSideID
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
USE MOD_Particle_Vars,           ONLY:LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: lengthPartTrajectory
INTEGER,INTENT(IN)                   :: PartID,SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: FlatBoundingBoxIntersection
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: dnk,alpha(2,3)! alpha(2,2): dummy because we are lazy
REAL                                 :: maxvalue,minvalue
INTEGER                              :: CNSideID
INTEGER                              :: i
!================================================================================================================================

CNSideID = GetCNSideID(SideID)

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
i   = 1
dnk = DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,CNSideID))

IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0. ! ÜBERPRÜFEN OB SIGN sinn macht
  alpha(1,1) = -HUGE(1.0)
  alpha(2,1) =  HUGE(1.0)
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
  END IF
END IF
i=3
dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,CNSideID))

IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0.
  alpha(1,3) = -HUGE(1.0)
  alpha(2,3) =  HUGE(1.0)
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(2*i-1,CNSideID))/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,CNSideID))&
                                                                                    +SideSlabIntervals(  2*i,CNSideID))/dnk!t_max
  END IF
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAX(alpha(1,1),alpha(1,3))!only x,z directions due to flat surface
minvalue=MIN(alpha(2,1),alpha(2,3))!only x,z directions due to flat surface

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
  IF((maxvalue.LT.0).AND.(minvalue.GT.0))THEN
    FlatBoundingBoxIntersection=.TRUE.
  ELSE
    IF((maxvalue.LT.lengthPartTrajectory+100.*epsMach).AND.(maxvalue+100.*epsMach.GT.0.))THEN
    !the first intersection is less than lengthPartTrajectory and greater 0
      FlatBoundingBoxIntersection=.TRUE.
    ELSE
      FlatBoundingBoxIntersection=.FALSE.
    END IF
  END IF
ELSE
  FlatBoundingBoxIntersection=.FALSE.
END IF

END FUNCTION FlatBoundingBoxIntersection


PPURE FUNCTION ComputeSurfaceDistance2(BiLinearCoeff,xi,eta,PartTrajectory)
!================================================================================================================================
! compute the required vector length to intersection
! ramsey paper algorithm 3.4
!================================================================================================================================
USE MOD_Preproc
USE MOD_Globals
!#ifdef CODE_ANALYZE
!USE MOD_Particle_Tracking_Vars,      ONLY:PartOut,MPIRankOut
!#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: BiLinearCoeff(1:3,4)
REAL,INTENT(IN)                      :: xi,eta
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeSurfaceDistance2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
!================================================================================================================================

t = SUM((xi*eta*BiLinearCoeff(:,1) + xi*BiLinearCoeff(:,2) + eta*BiLinearCoeff(:,3) + BiLinearCoeff(:,4))*PartTrajectory(:))

ComputeSurfaceDistance2 = t

END FUNCTION ComputeSurfaceDistance2


PPURE FUNCTION ComputeXi(eta,A1,A2)
!================================================================================================================================
! compute the xi value with algorithm 3.3 of Ramsey paper
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: eta
REAL,DIMENSION(4),INTENT(IN)         :: A1,A2
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeXi
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: a,b
!================================================================================================================================

a = eta* A2(1)+A2(2)
b = eta*(A2(1)-A1(1))+A2(2)-A1(2)

IF (ABS(B).GE.ABS(A)) THEN
  ! Both denominators are zero, no possible methods to calculate xi left
  IF (ALMOSTZERO(ABS(B))) THEN
    ComputeXi = HUGE(1.)
    RETURN
  END IF
  ComputeXi = (-eta*(A2(3)-A1(3))-(A2(4)-A1(4)))/b
ELSE
  ComputeXi = (-eta* A2(3)-A2(4))/a
END IF

END FUNCTION ComputeXi


SUBROUTINE ComputeBezierIntersectionPoint(nXiClip,nEtaClip,PartID,SideID,nInterSections,PartTrajectory,lengthPartTrajectory)
!===================================================================================================================================
! Compute the BezierInterSectionPoint in 3D, and verify if this intersection is satisfies
! a) alpha in [0,lenghtPartTrajectrory]  ! intersection point is between LastPartPos and PartPos
! b) alpha is not a multiple intersection
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Surfaces_Vars ,ONLY: XiArray,EtaArray,locAlpha,locXi,locEta
USE MOD_Particle_Surfaces_Vars ,ONLY: epsilontol,Beziercliphit
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierClipTolerance,BezierClipLocalTol,BezierClipMaxIntersec
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Globals                ,ONLY: UNIT_stdout,abort
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBezier
#endif /*CODE_ANALYZE*/
#if defined(CODE_ANALYZE) || USE_MPI
USE MOD_Globals                ,ONLY: myrank
#endif /*defined(CODE_ANALYZE)*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: nXiClip
INTEGER,INTENT(IN)                 :: nEtaClip
INTEGER,INTENT(IN)                 :: PartID
INTEGER,INTENT(IN)                 :: SideID
REAL,INTENT(IN)                    :: lengthPartTrajectory
REAL,INTENT(IN)                    :: PartTrajectory(3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)              :: nInterSections
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: tmpXi,tmpEta,Xi,Eta,alpha,alphaNorm,MinusXi,MinusEta
REAL,DIMENSION(3,0:NGeo,0:NGeo)    :: ReducedBezierControlPoints
REAL,DIMENSION(3)                  :: IntersectionVector
REAL                               :: deltaXi,deltaEta
INTEGER                            :: p,q,l,iDeCasteljau,iClip
LOGICAL                            :: isNewIntersection
INTEGER                            :: iInter
!===================================================================================================================================

! back transformation of sub-level clipping values to original bezier surface: ximean, etamean
!   xi-direction
IF(nXiClip.EQ.0)THEN
  Xi=0.
ELSE
  Xi=0.5*SUM(XiArray(:,nXiClip))
  DO iClip=nXiClip-1,1,-1
    Xi=XiArray(1,iClip)+0.5*(Xi+1)*(XiArray(2,iClip)-XiArray(1,iClip))
  END DO
END IF ! nXIClip
!   eta-direction
IF(nEtaClip.EQ.0)THEN
  Eta=0.
ELSE
  Eta=0.5*SUM(EtaArray(:,nEtaClip))
  DO iClip=nEtaClip-1,1,-1
    Eta=EtaArray(1,iClip)+0.5*(Eta+1)*(EtaArray(2,iClip)-EtaArray(1,iClip))
  END DO
END IF ! nEtaclip

! Calculate intersection value in 3D (De Casteljau)
tmpXi=XI
tmpEta=Eta

! shift from [-1,1] to [0,1]
Xi=0.5*(Xi+1)
Eta=0.5*(Eta+1)

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,*) ' xi,eta ',xi,eta
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF((ABS(eta).GT.BezierClipHit).OR.(ABS(xi).GT.BezierClipHit))THEN
  RETURN
END IF

MinusXi =1.0-Xi
MinusEta=1.0-Eta
! BEGIN DECASTELJAU ------------------------------------
! Wikipedia: "Although the algorithm is slower for most architectures
!             when compared with the direct approach, it is more numerically stable."
! DEBUG: keep decastejau or implement horner for direct evaluation
ReducedBezierControlPoints=BezierControlPoints3D(:,:,:,SideID)
l=NGeo-1
DO iDeCasteljau=1,NGeo
  DO q=0,l
    DO p=0,l
      ReducedBezierControlPoints(:,p,q)=MinusXi*ReducedBezierControlPoints(:,p,q  )  *MinusEta & ! A
                                       +MinusXi*ReducedBezierControlPoints(:,p,q+1)  *Eta      & ! B
                                       +     Xi*ReducedBezierControlPoints(:,p+1,q)  *MinusEta & ! C
                                       +     Xi*ReducedBezierControlPoints(:,p+1,q+1)*Eta        ! D

    END DO
  END DO
  l=l-1
END DO

! resulting point is ReducedBezierControlPoints(:,1,1)
IntersectionVector=ReducedBezierControlPoints(:,0,0)-LastPartPos(1:3,PartID)
! END DECASTELJAU ------------------------------------
! verify alpha and alphanorm
alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)
alphaNorm=alpha/lengthPartTrajectory

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    CALL CalcNormAndTangBezier(nVec=IntersectionVector,xi=tmpXi,eta=tmpeta,SideID=SideID)
    IPWRITE(UNIT_stdout,*) ' nVec   ',IntersectionVector
    IPWRITE(UNIT_stdout,*) ' PartTrajectory   ',PartTrajectory
    IPWRITE(UNIT_stdout,*) '<nVec,PartTrajectory>  ',DOT_PRODUCT(PartTrajectory,IntersectionVector)
    IPWRITE(UNIT_stdout,*) ' alpha,alphanorm ',alpha,alphaNorm
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF((alphaNorm.LE.1.0).AND.(alphaNorm.GT.-epsilontol))THEN
  ! found additional intersection point
  IF(nInterSections.GE.BezierClipMaxIntersec)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') ' nInterSections > BezierClipMaxIntersec'
    IPWRITE(UNIT_stdout,'(I0,A,I0)') ' PartID ', PartID
    IPWRITE(UNIT_stdout,'(I0,A,I0)') ' SideID ', SideID
    IPWRITE(UNIT_stdout,'(I0,A,E24.12)') ' BezierClipTolerance  ', BezierClipTolerance
    IPWRITE(UNIT_stdout,'(I0,A,E24.12)') ' BezierClipLocalTol ', BezierClipLocalTol
    IPWRITE(UNIT_stdout,'(I0,A,E18.12,1X,E18.12)') ' critical error! ',alpha,alphaNorm
    IPWRITE(UNIT_stdout,'(I0,A,E18.12,1X,E18.12)') ' critical error! ',alpha,alphaNorm
    IPWRITE(UNIT_stdout,'(I0,A)') ' locAlpha, locXi,locEta ' !/ lengthPartTrajectory '
    DO iInter=1,nInterSections
      WRITE(UNIT_stdout,'(I0,3(1X,E18.12))') iInter,locAlpha(iInter),locXi(iInter),locEta(iInter)
    END DO
    STOP
    RETURN
  END IF
  IF(nInterSections.EQ.0)THEN ! first intersection is always taken
    nInterSections=nIntersections+1
    locAlpha(nInterSections)=alpha
    locXi (nInterSections)=tmpXi
    locEta(nInterSections)=tmpEta
  ELSE
   ! check if new intersection does not coincidence with with any other intersection
   isNewIntersection=.TRUE.
   DO iInter=1,nInterSections
     ! remove multiple found intersections
     deltaXi =ABS(locXI(iInter)-tmpXi)
     IF(deltaXi.GT.0.01) CYCLE
     deltaEta=ABS(locEta(iInter)-tmpEta)
     IF(deltaEta.GT.0.01) CYCLE
     isNewIntersection=.FALSE.
     EXIT
     ! compare the found alpha value absolute and relative to find duplicate intersections
     !IF((locAlpha(iInter)-alpha).LT.SQRT(BezierClipTolerance)) isNewIntersection=.FALSE.
     !IF(ALMOSTEQUALRELATIVE(locAlpha(iInter),alpha,SQRT(BezierClipTolerance))) isNewIntersection=.FALSE.
#ifdef CODE_ANALYZE
     IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
       IF(PartID.EQ.PARTOUT)THEN
         IPWRITE(UNIT_stdout,*) ' locAlpha,alpha,tol ', locAlpha(iInter),alpha,SQRT(BezierClipTolerance)
         IPWRITE(UNIT_stdout,*) ' locXi,ETa,,...    ,', locXi(iInter),locEta(iInter),tmpXi,tmpEta
       END IF
     END IF
#endif /*CODE_ANALYZE*/
   END DO ! iInter=1,nInterSections
   IF(isNewIntersection)THEN
     nInterSections=nIntersections+1
     locAlpha(nInterSections)=alpha
     locXi (nInterSections)=tmpXi
     locEta(nInterSections)=tmpEta
   END IF
 END IF
END IF

END SUBROUTINE ComputeBezierIntersectionPoint


SUBROUTINE CheckXiClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                               ,iClipIter,nXiClip,nEtaClip&
                               ,nInterSections,iPart,SideID)
!================================================================================================================================
! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
! original article:
!   author = {Nishita, Tomoyuki and Sederberg, Thomas W. and Kakimoto, Masanori},
!   title = {Ray Tracing Trimmed Rational Surface Patches},
!   year = {1990},
! book:
!   author = {Farin, Gerald},
!   title = {Curves and Surfaces for CAGD: A Practical Guide},
!   year = {2002},
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:XiArray!,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipLocalTol,FacNchooseK!,BezierClipTolerance
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierSplitLimit
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
#ifdef CODE_ANALYZE
USE MOD_Globals,                 ONLY:UNIT_stdOut
USE MOD_Globals,                 ONLY:MyRank
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: lengthPartTrajectory
REAL,INTENT(INOUT)                   :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: SideID,iPart
REAL,INTENT(IN),DIMENSION(1:3)       :: PartTrajectory
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                  :: iClipIter
INTEGER,INTENT(INOUT)                  :: nXiClip,nEtaClip,nInterSections
INTEGER(KIND=2),INTENT(INOUT)          :: ClipMode
REAL,DIMENSION(2,2),INTENT(INOUT)      :: LineNormVec
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:NGeo,0:NGeo)        :: BezierControlPoints1D
REAL                                 :: minmax(1:2,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp(2,0:NGeo,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp2(2,0:NGeo,0:NGeo)
INTEGER                              :: p,q,l
REAL                                 :: XiMin,XiMax,XiSplit,XiTmp
REAL                                 :: PlusXi,MinusXi
INTEGER                              :: tmpnClip,tmpnXi,tmpnEta
REAL                                 :: xiup(0:NGeo),xidown(0:NGeo)
REAL                                 :: XiBuf(0:NGeo,0:NGeo)
REAL                                 :: dmin,dmax
INTEGER(KIND=2)                      :: tmpClipMode
REAL,DIMENSION(2,2)                  :: tmpLineNormVec
!================================================================================================================================

! Bezier Clip (and Split) in xi
DO q=0,NGeo
  DO p=0,NGeo
    BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec(:,1))
  END DO
END DO
DO l=0,NGeo
  minmax(1,l)=MINVAL(BezierControlPoints1D(l,:))
  minmax(2,l)=MAXVAL(BezierControlPoints1D(l,:))
END DO ! l
dmin=MINVAL(minmax(1,:))
dmax=MAXVAL(minmax(2,:))
#ifdef CODE_ANALYZE
 IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
   IF(iPart.EQ.PARTOUT)THEN
     IPWRITE(UNIT_stdOut,*) ' minval-xi',minmax(1,:)
     IPWRITE(UNIT_stdOut,*) ' maxval-xi',minmax(2,:)
     IPWRITE(UNIT_stdout,*) ' dmax-dmin-xi ',(dmax-dmin)
     IPWRITE(UNIT_stdout,*) ' dmax,dmin-xi ',dmax,dmin
   END IF
 END IF
#endif /*CODE_ANALYZE*/

!print*,'dmax,min,...',dmax,dmin,ABS(dmax-dmin),BezierClipTolerance
! 1D abort criterion from
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},
!      YEAR = {2005},
IF(dMax*dMin.GT.0.)THEN ! no sign change with dMax,dMin, hence, no intersection
  ClipMode=-1
  RETURN
END IF
IF(ABS(dmax-dmin).LT.BezierClipLocalTol)THEN ! current patch is converged in xi, then skip xi
  IF(ClipMode.EQ.3) THEN  ! eta has already converged
    ClipMode=5            ! no more clipping, we have converged
    RETURN                ! or stop
  ELSE
    ClipMode=4            ! xi is converged, but not eta
    RETURN
  END IF
ELSE ! xi not converged, next clip should be in eta
  ClipMode=2 ! after clipping in xi, clip in eta
END IF

! we perform the clip  (with XiMin and XiMax) in Xi direction
!             or split (if (XiMax-XiMin).GT.BezierSplitLimit) in Xi direction

! calc Smin and Smax and check boundaries
CALL CalcSminSmax2(minmax,XiMin,XiMax,nXiClip)
#ifdef CODE_ANALYZE
 IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
   IF(iPart.EQ.PARTOUT)THEN
     IPWRITE(UNIT_stdout,*) ' XiMin,XiMax ',XiMin,XiMax,nXiClip
   END IF
 END IF
#endif /*CODE_ANALYZE*/

! check if diverged
IF((XiMin.EQ.1.5).OR.(XiMax.EQ.-1.5))THEN
  ClipMode=-1
  RETURN
END IF

IF(XiMin.GT.XiMax)THEN
  ! output, should never ever happen
!  print*,'swwwaaaaaaaap xi',XiMin,XiMax
  XiTmp = XiMax
  Ximax = XiMin
  XiMin = XiTmp
END IF


! count number of clips in xi and eta direction
nXiClip=nXiClip+1

! 1.) CLIPPING xi
IF((XiMax-XiMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
  XiSplit=0.5*(XiMax+XiMin)
  ! first split: xi upper
  ! set mapping array
  XiArray(:,nXiClip)=(/XiSplit,XiMax/)
  BezierControlPoints2D_temp=0.
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  IF(XiMax.NE.1.0)THEN
    PlusXi=1.0+XiMax
    MinusXi=1.0-XiMax
    ! compute the required stuff || pseudo Horner or precomputation
    xiup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    xidown(0)=1.0
    DO l=1,NGeo
      xiup   (l)     =xiup   (l-1)*PlusXi
      xidown (l)     =xidown (l-1)*MinusXi
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                           +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
        END DO
      END DO
    END DO

  ELSE
    BezierControlPoints2D_temp=BezierControlPoints2D
  END IF

  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  BezierControlPoints2D_temp2=0.
  PlusXi=(XiSplit+1.0)/(XiMax+1.0)
  ! MinusXi= 2.0*(1-PlusXi)-1
  ! MinusXi= 1.0-2.0*(1.0-s)+1.0
  MinusXi=2.0*PlusXi
  ! PlusXi=1+ 2.0*(1-s)-1
  PlusXi=2.0-2.0*PlusXi
  ! compute the required stuff || pseudo Horner or precomputation
  xiup(0)=1.0
  ! caution, here the indicies are switched from n-j to j  for **down
  xidown(0)=1.0
  DO l=1,NGeo
    xiup   (l)     =xiup   (l-1)*PlusXi
    xidown (l)     =xidown (l-1)*MinusXi
  END DO ! l=0,NGeo
  DO p=0,NGeo
    DO l=0,p
      XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
    END DO ! l=0,p
  END DO ! p=0,NGeo

  DO q=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
                                               +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
      END DO
    END DO
  END DO

  ! Bezier Split
  tmpnClip =iClipIter
  ! backup current split level to compute correct intersection, required for back-trafo of intervals
  tmpnXi   =nXiClip
  tmpnEta  =nEtaClip
  tmpLineNormVec=LineNormVec
  tmpClipMode   =ClipMode
  ! MAYBE set ClipMode for NEXT clip
  ! HERE, ClipMode currently set above
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
          IPWRITE(UNIT_stdout,*) ' split xi-upper '
          IPWRITE(UNIT_stdout,*) ' XiMin,XiMax ',XiArray(:,nXiClip)
  !        CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
  ! HERE, ClipMode currently set above
  ! Perform split xi-upper
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)


  ! second split: xi lower
  ! restore values to allow for correct back-trafo of intervals (required for intersectionpoint)
  ! backup current split level to compute correct intersection, required for back-trafo of intervals
  iClipIter   =tmpnClip
  nXiClip     =tmpnXi
  nEtaClip    =tmpnEta
  LineNormVec =tmpLineNormVec
  ClipMode    =tmpClipMode

  ! set mapping array
  XiArray(:,nXiClip)=(/XiMin,XiSplit/)
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  BezierControlPoints2D_temp=0.
  PlusXi=1.0+XiSplit
  MinusXi=1.0-XiSplit
  ! compute the required stuff || pseudo Horner or precomputation
  xiup(0)=1.0
  ! caution, here the indicies are switched from n-j to j  for **down
  xidown(0)=1.0
  DO l=1,NGeo
    xiup   (l)     =xiup   (l-1)*PlusXi
    xidown (l)     =xidown (l-1)*MinusXi
  END DO ! l=0,NGeo
  DO p=0,NGeo
    DO l=0,p
      XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
    END DO ! l=0,p
  END DO ! p=0,NGeo

  DO q=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                         +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
      END DO
    END DO
  END DO
  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  IF(XiMin.NE.-1.0)THEN
    BezierControlPoints2D_temp2=0.
    PlusXi=(XiMin+1.0)/(XiSplit+1.0)
    ! MinusXi= 2.0*(1-PlusXi)-1
    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
    MinusXi=2.0*PlusXi
    ! PlusXi=1+ 2.0*(1-s)-1
    PlusXi=2.0-2.0*PlusXi
    ! compute the required stuff || pseudo Horner or precomputation
    xiup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    xidown(0)=1.0
    DO l=1,NGeo
      xiup   (l)     =xiup   (l-1)*PlusXi
      xidown (l)     =xidown (l-1)*MinusXi
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
                                                 +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
        END DO
      END DO
    END DO
  ELSE
    BezierControlPoints2D_temp2=BezierControlPoints2D_temp
  END IF
  ! MAYBE set ClipMode for NEXT clip
  ! HERE, ClipMode currently set above
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
          IPWRITE(UNIT_stdout,*) ' split xi-lower '
          IPWRITE(UNIT_stdout,*) ' XiMin,XiMax ',XiArray(:,nXiClip)
  !        CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
  ! HERE, ClipMode currently set above
  ! Perform split xi-lower
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
                          ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
  ! and we are done
  ClipMode=-1
  ! after recursive steps, we are done!
ELSE  ! no split necessary, only a clip

  ! set mapping array
  XiArray(:,nXiClip)=(/XiMin,XiMax/)
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  IF(XiMax.NE.1.0)THEN
    BezierControlPoints2D_temp=0.
    PlusXi=1.0+XiMax
    MinusXi=1.0-XiMax
    ! compute the required stuff || pseudo Horner or precomputation
    xiup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    xidown(0)=1.0
    DO l=1,NGeo
      xiup   (l)     =xiup   (l-1)*PlusXi
      xidown (l)     =xidown (l-1)*MinusXi
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                           +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
        END DO
      END DO
    END DO
    BezierControlPoints2D=BezierControlPoints2D_temp
  END IF
  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  IF(XiMin.NE.-1.0)THEN
    BezierControlPoints2D_temp=0.
    PlusXi=(XiMin+1.0)/(XiMax+1.0)
    ! MinusXi= 2.0*(1-PlusXi)-1
    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
    MinusXi=2.0*PlusXi
    ! PlusXi=1+ 2.0*(1-s)-1
    PlusXi=2.0-2.0*PlusXi
    ! compute the required stuff || pseudo Horner or precomputation
    xiup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    xidown(0)=1.0
    DO l=1,NGeo
      xiup   (l)     =xiup   (l-1)*PlusXi
      xidown (l)     =xidown (l-1)*MinusXi
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)                  &
                                                +BezierControlPoints2D(:,NGeo-l,q)*XiBuf(p,l)
        END DO
      END DO
    END DO
    BezierControlPoints2D=BezierControlPoints2D_temp
  END IF
  !CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
  !               ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)

END IF ! decision between Clip or Split

END SUBROUTINE CheckXiClip

SUBROUTINE CheckEtaClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
                               ,iClipIter,nXiClip,nEtaClip&
                               ,nInterSections,iPart,SideID)
!================================================================================================================================
! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
! original article:
!   author = {Nishita, Tomoyuki and Sederberg, Thomas W. and Kakimoto, Masanori},
!   title = {Ray Tracing Trimmed Rational Surface Patches},
!   year = {1990},
! book:
!   author = {Farin, Gerald},
!   title = {Curves and Surfaces for CAGD: A Practical Guide},
!   year = {2002},
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:EtaArray
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipLocalTol,FacNchooseK
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierSplitLimit
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
#ifdef CODE_ANALYZE
USE MOD_Globals,                 ONLY:UNIT_stdOut
USE MOD_Globals,                 ONLY:MyRank
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: lengthPartTrajectory
REAL,INTENT(INOUT)                   :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: SideID,iPart
REAL,INTENT(IN),DIMENSION(1:3)       :: PartTrajectory
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                  :: iClipIter
INTEGER,INTENT(INOUT)                  :: nXiClip,nEtaClip,nInterSections
INTEGER(KIND=2),INTENT(INOUT)          :: ClipMode
REAL,DIMENSION(2,2),INTENT(INOUT)      :: LineNormVec
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:NGeo,0:NGeo)        :: BezierControlPoints1D
REAL                                 :: minmax(1:2,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp(2,0:NGeo,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp2(2,0:NGeo,0:NGeo)
INTEGER                              :: p,q,l
REAL                                 :: EtaMin,EtaMax,EtaSplit,EtaTmp
REAL                                 :: PlusEta,MinusEta
INTEGER                              :: tmpnClip,tmpnXi,tmpnEta
REAL                                 :: etaup(0:NGeo),etadown(0:NGeo)
REAL                                 :: EtaBuf(0:NGeo,0:NGeo)
REAL                                 :: dmin,dmax
INTEGER(KIND=2)                      :: tmpClipMode
REAL,DIMENSION(2,2)                  :: tmpLineNormVec
!================================================================================================================================

! Bezier Clip (and Split) in eta
DO q=0,NGeo
  DO p=0,NGeo
    BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec(:,2))
  END DO
END DO
DO l=0,NGeo
  minmax(1,l)=MINVAL(BezierControlPoints1D(:,l))
  minmax(2,l)=MAXVAL(BezierControlPoints1D(:,l))
END DO ! l
dmin=MINVAL(minmax(1,:))
dmax=MAXVAL(minmax(2,:))
#ifdef CODE_ANALYZE
 IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
   IF(iPart.EQ.PARTOUT)THEN
     IPWRITE(UNIT_stdOut,*) ' minval-eta',minmax(1,:)
     IPWRITE(UNIT_stdOut,*) ' maxval-eta',minmax(2,:)
     IPWRITE(UNIT_stdout,*) ' dmax-dmin-eta ',(dmax-dmin)
     IPWRITE(UNIT_stdout,*) ' dmax,dmin-eta ',dmax,dmin
   END IF
 END IF
#endif /*CODE_ANALYZE*/

!print*,'dmax,min,...',dmax,dmin,ABS(dmax-dmin),BezierClipTolerance
! 1D abort criterion from
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},
!      YEAR = {2005},
IF(dMax*dMin.GT.0.)THEN ! no sign change with dMax,dMin, hence, no intersection
  ClipMode=-1
  RETURN
END IF
IF(ABS(dmax-dmin).LT.BezierClipLocalTol)THEN ! current patch is converged in eta, then skip eta
  IF(ClipMode.EQ.4) THEN  ! xi has already converged
    ClipMode=5            ! no more clipping, we have converged
    RETURN                ! or stop
  ELSE
    ClipMode=3            ! eta is converged, but not xi
    RETURN
  END IF
ELSE ! eta not converged, next clip should be in xi
  ClipMode=1 ! after clipping in eta, clip in xi
END IF

! we perform the clip  (with EtaMin and EtaMax) in Eta direction
!             or split (if (EtaMax-EtaMin).GT.BezierSplitLimit) in Eta direction

! calc Smin and Smax and check boundaries
CALL CalcSminSmax2(minmax,Etamin,Etamax,nEtaClip)
#ifdef CODE_ANALYZE
 IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
   IF(iPart.EQ.PARTOUT)THEN
     IPWRITE(UNIT_stdout,*) ' EtaMin,EtaMax ',EtaMin,EtaMax,nEtaClip
   END IF
 END IF
#endif /*CODE_ANALYZE*/

! check if diverged
IF((EtaMin.EQ.1.5).OR.(EtaMax.EQ.-1.5)) THEN
  ClipMode=-1
  RETURN
END IF

IF(EtaMin.GT.EtaMax)THEN
!  print*,'swwwaaaaaaaap etta',etamin,etamax
  EtaTmp = EtaMax
  Etamax = EtaMin
  EtaMin = EtaTmp
END IF

nEtaClip=nEtaClip+1
! 2.) CLIPPING eta
IF((EtaMax-EtaMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
  EtaSplit=0.5*(EtaMax+EtaMin)
  ! first split: eta upper
  ! set mapping array
  EtaArray(:,nEtaClip)=(/EtaSplit,EtaMax/)
  BezierControlPoints2D_temp=0.
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  IF(EtaMax.NE.1.0)THEN
    PlusEta=1.0+EtaMax
    MinusEta=1.0-EtaMax
    ! compute the required stuff || pseudo Horner or precomputation
    Etaup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    Etadown(0)=1.0
    DO l=1,NGeo
      Etaup   (l)     =Etaup   (l-1)*PlusEta
      Etadown (l)     =Etadown (l-1)*MinusEta
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                           +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
        END DO
      END DO
    END DO

  ELSE
    BezierControlPoints2D_temp=BezierControlPoints2D
  END IF

  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  BezierControlPoints2D_temp2=0.
  PlusEta=(EtaSplit+1.0)/(EtaMax+1.0)
  ! MinusXi= 2.0*(1-PlusXi)-1
  ! MinusXi= 1.0-2.0*(1.0-s)+1.0
  MinusEta=2.0*PlusEta
  ! PlusXi=1+ 2.0*(1-s)-1
  PlusEta=2.0-2.0*PlusEta
  ! compute the required stuff || pseudo Horner or precomputation
  Etaup(0)=1.0
  ! caution, here the indicies are switched from n-j to j  for **down
  Etadown(0)=1.0
  DO l=1,NGeo
    Etaup   (l)     =Etaup   (l-1)*PlusEta
    Etadown (l)     =Etadown (l-1)*MinusEta
  END DO ! l=0,NGeo
  DO p=0,NGeo
    DO l=0,p
      EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
    END DO ! l=0,p
  END DO ! p=0,NGeo

  DO q=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                               +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
      END DO
    END DO
  END DO

  ! Bezier Split
  tmpnClip =iClipIter
  ! backup current split level to compute correct intersection, required for back-trafo of intervals
  tmpnXi   =nXiClip
  tmpnEta  =nEtaClip
  tmpLineNormVec=LineNormVec
  tmpClipMode   =ClipMode
  ! MAYBE set ClipMode for NEXT clip
  ! HERE, ClipMode currently set above
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
          IPWRITE(UNIT_stdout,*) ' split eta-upper '
          IPWRITE(UNIT_stdout,*) ' EtaMin,EtaMax ',EtaArray(:,nEtaClip)
          !CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
  ! HERE, ClipMode currently set above
  ! Perform split eta-upper
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)

  ! second split: eta lower
  ! restore values to allow for correct back-trafo of intervals (required for intersectionpoint)
  iClipIter   =tmpnClip
  nXiClip     =tmpnXi
  nEtaClip    =tmpnEta
  LineNormVec =tmpLineNormVec
  ClipMode    =tmpClipMode

  ! set mapping array
  EtaArray(:,nEtaClip)=(/EtaMin,EtaSplit/)
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  BezierControlPoints2D_temp=0.
  PlusEta=1.0+EtaSplit
  MinusEta=1.0-EtaSplit
  ! compute the required stuff || pseudo Horner or precomputation
  Etaup(0)=1.0
  ! caution, here the indicies are switched from n-j to j  for **down
  Etadown(0)=1.0
  DO l=1,NGeo
    Etaup   (l)     =Etaup   (l-1)*PlusEta
    Etadown (l)     =Etadown (l-1)*MinusEta
  END DO ! l=0,NGeo
  DO p=0,NGeo
    DO l=0,p
      EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
    END DO ! l=0,p
  END DO ! p=0,NGeo

  DO q=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                         +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
      END DO
    END DO
  END DO
  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  IF(EtaMin.NE.-1.0)THEN
    BezierControlPoints2D_temp2=0.
    PlusEta=(EtaMin+1.0)/(EtaSplit+1.0)
    ! MinusXi= 2.0*(1-PlusXi)-1
    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
    MinusEta=2.0*PlusEta
    ! PlusXi=1+ 2.0*(1-s)-1
    PlusEta=2.0-2.0*PlusEta
    ! compute the required stuff || pseudo Horner or precomputation
    Etaup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    Etadown(0)=1.0
    DO l=1,NGeo
      Etaup   (l)     =Etaup   (l-1)*PlusEta
      Etadown (l)     =Etadown (l-1)*MinusEta
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                                 +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
        END DO
      END DO
    END DO
  ELSE
    BezierControlPoints2D_temp2=BezierControlPoints2D_temp
  END IF
  tmpnClip      =iClipIter
  tmpnXi        =nXiClip
  tmpnEta       =nEtaClip
  tmpLineNormVec=LineNormVec
  tmpClipMode   =ClipMode
  ! MAYBE set ClipMode for NEXT clip
  ! HERE, ClipMode currently set above
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(iPart.EQ.PARTOUT)THEN
          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
          IPWRITE(UNIT_stdout,*) ' split eta-lower '
          IPWRITE(UNIT_stdout,*) ' EtaMin,EtaMax ',EtaArray(:,nEtaClip)
          !CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
  ! HERE, ClipMode currently set above
  ! Perform split eta-lower
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
  ! we are done
  ClipMode=-1
  !! after recursive steps, we are done!
ELSE  ! no split necessary, only a clip

  ! set mapping array
  EtaArray(:,nEtaClip)=(/EtaMin,EtaMax/)
  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
  IF(EtaMax.NE.1.0)THEN
    BezierControlPoints2D_temp=0.
    PlusEta=1.0+EtaMax
    MinusEta=1.0-EtaMax
    ! compute the required stuff || pseudo Horner or precomputation
    Etaup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    Etadown(0)=1.0
    DO l=1,NGeo
      Etaup   (l)     =Etaup   (l-1)*PlusEta
      Etadown (l)     =Etadown (l-1)*MinusEta
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                           +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
        END DO
      END DO
    END DO
    BezierControlPoints2D=BezierControlPoints2D_temp
  END IF
  ! BOTTOM (mirrored Bernstein Basis evaluation)
  ! s = (smin+1)/(smax+1) for [-1, +1]
  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
  IF(EtaMin.NE.-1.0)THEN
    BezierControlPoints2D_temp=0.
    PlusEta=(EtaMin+1.0)/(EtaMax+1.0)
    ! MinusXi= 2.0*(1-PlusXi)-1
    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
    MinusEta=2.0*PlusEta
    ! PlusXi=1+ 2.0*(1-s)-1
    PlusEta=2.0-2.0*PlusEta
    ! compute the required stuff || pseudo Horner or precomputation
    Etaup(0)=1.0
    ! caution, here the indicies are switched from n-j to j  for **down
    Etadown(0)=1.0
    DO l=1,NGeo
      Etaup   (l)     =Etaup   (l-1)*PlusEta
      Etadown (l)     =Etadown (l-1)*MinusEta
    END DO ! l=0,NGeo
    DO p=0,NGeo
      DO l=0,p
        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
      END DO ! l=0,p
    END DO ! p=0,NGeo

    DO q=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)             &
                                                 +BezierControlPoints2D     (:,q,NGeo-l)*EtaBuf(p,l)
        END DO
      END DO
    END DO
    BezierControlPoints2D=BezierControlPoints2D_temp
  END IF
  ! after recursive steps, we are done!
END IF ! decision between Clip or Split

END SUBROUTINE CheckEtaClip

#ifdef CODE_ANALYZE
SUBROUTINE OutputTrajectory(PartID,PartPos,PartTrajectory,lengthPartTrajectory)
!===================================================================================================================================
! subroutine to print particle trajectory, lengthPartTrajectory and Last and Current position in matlab format
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars, ONLY: LastPartPos
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
REAL, INTENT(IN)                 :: PartTrajectory(1:3), lengthPartTrajectory
REAL, INTENT(IN)                 :: PartPos(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

WRITE(UNIT_stdout,'(A,3(E24.12,A))') ' LastPartPos = [ ',LastpartPos(1,PartID), ','  &
                                                        ,LastpartPos(2,PartID), ','  &
                                                        ,LastpartPos(3,PartID), '];'

WRITE(UNIT_stdout,'(A,3(E24.12,A))') ' PartPosition   = [ ',PartPos(1), ','  &
                                                         ,PartPos(2), ','  &
                                                         ,PartPos(3), '];'

WRITE(UNIT_stdout,'(A,3(E24.12,A))') ' PartTrajectory = [ ',PartTrajectory(1) ,','  &
                                                           ,PartTrajectory(2) ,','  &
                                                           ,PartTrajectory(3) ,'];'

WRITE(UNIT_stdout,*) ' lengthPartTrajectory = ', lengthPartTrajectory

END SUBROUTINE OutputTrajectory
#endif /*CODE_ANALYZE*/


SUBROUTINE calcLineNormVec3(BezierControlPoints2D,LineNormVec)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
! Ls is the search direction. Ls is constructed via the combination of the left and right bounding vector
! e.g. clip in xi direction: 1) Ls vector is constructed in eta direction, because the left and right border are to be
!                               clipped away ("reducing the patch by moving the left and right boundary")
!                            2) in order to compute the distance between Ls and each control point, the vector has to be
!                               normalized. If the length is zero, then the unit normal vector in eta direction is used.
!                            3) and right rotated. The rotation constructs the vector which is perpendicular to Ls,
!                               the scalar product <LineNormVec,ControlPoint> gives the distance
!================================================================================================================================
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Mesh_Vars,               ONLY:NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: LineNormVec(1:2,1:2)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length,doPro
REAL,DIMENSION(2)                    :: LXi, Leta,Mbar,Mbar2
!================================================================================================================================

! compute Lxi vector
! 1) the initial Lxi vector is the combination of the two bounding vectors  which point in eta direction
!     to get the 1D distances of each point via scalar product, we have to right-rotate the vector
LXi=(BezierControlPoints2D(:,   0,NGeo)-BezierControlPoints2D(:,   0,   0))+&
    (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,NGeo,   0))

Mbar =(BezierControlPoints2D(:,   0,NGeo)-BezierControlPoints2D(:,   0,   0))
MBar2=(BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,NGeo,   0))

! 2) normalization
Length=SQRT(DOT_PRODUCT(LXi,LXi))
IF(Length.EQ.0.)THEN
  print*,'AAARG'
  STOP
  LXi=(/0.,1./)
ELSE
  LXi=LXi/Length
END IF

!print*,'Lxi',Lxi

! compute Leta vector
! 1) the initial Leta vector is the combination of the two bounding vectors  which point in xi direction
!    to get the 1D distances of each point via scalar product, we have to right-rotate the vector
Leta=(BezierControlPoints2D(:,NGeo,   0)-BezierControlPoints2D(:,0,   0))+&
     (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,0,NGeo))


Mbar =(BezierControlPoints2D(:,Ngeo,   0)-BezierControlPoints2D(:,   0,   0))
MBar2=(BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,   0,NGeo))
!print*,'DoProEta',DOT_PRODUCT(MBar,MBar2)

! 2) normalization
Length=SQRT(DOT_PRODUCT(Leta,Leta))
IF(Length.EQ.0)THEN
  LXi=(/1.,0./)
  print*,'BBBBRG'
  STOP
ELSE
  Leta=Leta/Length
END IF

doPro=DOT_PRODUCT(Lxi,Leta)

! 3) rotate  both vectors by -90 degree
LineNormVec(1,1) =-LXi (2)
LineNormVec(2,1) = LXi (1)
LineNormVec(1,2) =-Leta(2) ! stephen meint minus1
LineNormVec(2,2) = Leta(1)

END SUBROUTINE calcLineNormVec3


SUBROUTINE CalcSminSmax2(minmax,Smin,Smax,iter)
!================================================================================================================================
! find upper and lower intersection with convex hull (or no intersection)
! find the largest and smallest roots of the convex hull, pre-sorted values minmax(:,:) are required
! % convex-hull check by
! 1) check upper line for an intersection with x-axis
! 2) check lower line for an intersection with x-axis
! 3) check if both ends have a sign change
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo,Xi_NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance!,BezierClipHit
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: minmax(1:2,0:NGeo)
INTEGER,INTENT(IN)                   :: iter
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Smin,Smax
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: tmp,m
INTEGER                              :: p,q
!================================================================================================================================

Smin=1.5
Smax=-1.5

! check upper/lower line for an intersection with the x-axis
DO p=0,NGeo
  DO q=0,NGeo
    IF(p.EQ.q) CYCLE
    ! 1) check for upper line
    IF(minmax(2,p)*minmax(2,q).LE.0.)THEN
      m=(minmax(2,q)-minmax(2,p))/(Xi_NGeo(q)-Xi_NGeo(p));
      IF(m.LT.0.)THEN ! move right boundary
        tmp=Xi_NGeo(q)-minmax(2,q)/m;
        smax=MAX(smax,tmp);
      END IF
      IF(m.GT.0.)THEN ! move left boundary
        tmp=Xi_NGeo(q)-minmax(2,q)/m;
        smin=MIN(smin,tmp);
      END IF
    END IF
    ! 2) check for lower line
    IF(minmax(1,p)*minmax(1,q).LE.0.)THEN
      m=(minmax(1,q)-minmax(1,p))/(Xi_NGeo(q)-Xi_NGeo(p));
      IF(m.GT.0.)THEN ! move right boundary
        tmp=Xi_NGeo(q)-minmax(1,q)/m;
        smax=MAX(smax,tmp);
      END IF
      IF(m.LT.0.)THEN ! move right boundary
        tmp=Xi_NGeo(q)-minmax(1,q)/m;
        smin=MIN(smin,tmp);
      END IF
    END IF
  END DO ! p=0,PP_N
END DO ! q=0,PP_N

! 3) check left and right boundary
IF(minmax(1,   0)*minmax(2,   0).LT.0.) smin=-1.
IF(minmax(1,NGeo)*minmax(2,NGeo).LT.0.) smax= 1.

! adjust Smin and Smax to increase the current range
!! adapted from: 1997, Campagna, Ray tracing of spline surfaces
! modification. initial method works with smax=smax+(smax-smax,0)*eps*f
IF(Smax.GT.-1.5)THEN
  Smax=MIN(Smax+20.*BezierClipTolerance,1.0)
  !Smax=MIN(Smax+100.*BezierClipTolerance,BezierClipHit)
END IF
IF(Smin.LT.1.5)THEN
  Smin=MAX(Smin-20.*BezierClipTolerance,-1.0)
  !Smin=MAX(Smin-100.*BezierClipTolerance,-BezierClipHit)
END IF

! in first iteration direction
! due to tolerance issues in first clip, it is not allowed to diverge
! example: particle intersects close to the edge,corner, the NEXT patch
! has to be increased slightly
IF(iter.EQ.0)THEN
  IF(Smin.EQ.1.5) SMin=-1. !BezierClipHit ! BezierClipHit=1+BezierClipTolerance
  IF(Smax.EQ.-1.5)SMax=1.  !BezierClipHit
END IF

END SUBROUTINE calcSminSmax2

END MODULE MOD_Particle_Intersection
