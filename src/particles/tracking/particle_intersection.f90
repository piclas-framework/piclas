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

MODULE MOD_Particle_InterSection
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ComputePlanarCurvedIntersection
  MODULE PROCEDURE ComputePlanarCurvedIntersection
END INTERFACE

INTERFACE ComputePlanarRectIntersection
  MODULE PROCEDURE ComputePlanarRectIntersection
END INTERFACE

INTERFACE ComputeBilinearIntersection
  MODULE PROCEDURE ComputeBilinearIntersection
END INTERFACE

INTERFACE ComputeCurvedIntersection
  MODULE PROCEDURE ComputeCurvedIntersection
END INTERFACE

INTERFACE ComputeAuxBCIntersection
  MODULE PROCEDURE ComputeAuxBCIntersection
END INTERFACE

#ifdef CODE_ANALYZE
INTERFACE OutputTrajectory
  MODULE PROCEDURE OutputTrajectory
END INTERFACE
#endif /*CODE_ANALYZE*/
PUBLIC::ComputePlanarRectIntersection
PUBLIC::ComputePlanarCurvedIntersection
PUBLIC::ComputeBilinearIntersection
PUBLIC::ComputeCurvedIntersection
PUBLIC::ComputeAuxBCIntersection
#ifdef CODE_ANALYZE
PUBLIC::OutputTrajectory
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


SUBROUTINE IntersectionWithWall(PartTrajectory,alpha,iPart,iLocSide,Element,TriNum)!, IntersectionPos)
!===================================================================================================================================
! Compute the Intersection with bilinear surface by approximating the surface with two triangles
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,          ONLY : lastPartPos,PartState
!USE MOD_Particle_Mesh_Vars,     ONLY : GEO
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iPart
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
REAL, INTENT(IN)                 :: PartTrajectory(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)               :: alpha !,IntersectionPos(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: Node1, Node2
REAL                             :: PoldX, PoldY, PoldZ, PnewX, PnewY, PnewZ, nx, ny, nz, nVal
REAL                             :: bx,by,bz, ax,ay,az, dist!, PoldStarX, PoldStarY, PoldStarZ
REAL                             :: xNod, yNod, zNod
REAL                             :: Vector1(1:3), Vector2(1:3)!, VectorShift(1:3)
!===================================================================================================================================

PoldX = LastPartPos(1,iPart)
PoldY = LastPartPos(2,iPart)
PoldZ = LastPartPos(3,iPart)
PnewX = PartState(1,iPart)
PnewY = PartState(2,iPart)
PnewZ = PartState(3,iPart)

!xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
!yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
!zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,Element)+1)
yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,Element)+1)
zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,Element)+1)

!---- Calculate normal vector:

Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

!Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
!Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
!Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod
Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,iLocSide,Element)+1) - xNod
Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,iLocSide,Element)+1) - yNod
Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,iLocSide,Element)+1) - zNod

!Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
!Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
!Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod
Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,iLocSide,Element)+1) - xNod
Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,iLocSide,Element)+1) - yNod
Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,iLocSide,Element)+1) - zNod


nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

nVal = SQRT(nx*nx + ny*ny + nz*nz)

nx = nx/nVal
ny = ny/nVal
nz = nz/nVal

!---- Calculate Intersection

!--- the following has been used for impulse computations, not implemented yet?
!   IF (nx.NE.0) PIC%InverseImpulseX(iPart) = .NOT.(PIC%InverseImpulseX(iPart))
!   IF (ny.NE.0) PIC%InverseImpulseY(iPart) = .NOT.(PIC%InverseImpulseY(iPart))
!   IF (nz.NE.0) PIC%InverseImpulseZ(iPart) = .NOT.(PIC%InverseImpulseZ(iPart))

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

!   PoldStarX = PoldX + 2 * dist * nx
!   PoldStarY = PoldY + 2 * dist * ny
!   PoldStarZ = PoldZ + 2 * dist * nz

!VectorShift(1) = PnewX - PoldX
!VectorShift(2) = PnewY - PoldY
!VectorShift(3) = PnewZ - PoldZ

alpha = PartTrajectory(1) * nx + PartTrajectory(2) * ny + PartTrajectory(3) * nz
IF(ABS(alpha).GT.0.) alpha = dist / alpha


!IntersectionPos(1) = PoldX + alpha * PartTrajectory(1)
!IntersectionPos(2) = PoldY + alpha * PartTrajectory(2)
!IntersectionPos(3) = PoldZ + alpha * PartTrajectory(3)

RETURN

END SUBROUTINE IntersectionWithWall


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
INTEGER,INTENT(IN)                :: PartID,SideID!,ElemID,locSideID
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
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,locSideDistance
REAL                              :: sdet
REAL                              :: epsLoc
LOGICAL                           :: CriticalParallelInSide
!INTEGER                           :: flip
!===================================================================================================================================

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(110("-"))')
      WRITE(UNIT_stdout,'(A)') '     | Output of planar face constants: '
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | SideNormVec  : ',SideNormVec(1:3,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint1: ',BezierControlPoints3D(:,0,0,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint2: ',BezierControlPoints3D(:,NGeo,0,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint3: ',BezierControlPoints3D(:,0,NGeo,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint4: ',BezierControlPoints3D(:,NGeo,NGeo,SideID)
    END IF
  END IF
#endif /*CODE_ANALYZE*/
! set alpha to minus 1, asume no intersection
alpha = -1.0
xi    = -2.
eta   = -2.
isHit = .FALSE.

! new with flip
IF(flip.EQ.0)THEN
  NormVec     = SideNormVec(1:3,SideID)
  locDistance = SideDistance(SideID)
ELSE
  NormVec     = -SideNormVec(1:3,SideID)
  locDistance = -SideDistance(SideID)
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

!IF((alphaNorm.GT.OnePlusEps) .OR.(alphaNorm.LT.-epsilontol))THEN
IF((alphaNorm.GT.1.0) .OR.(alphaNorm.LT.-epsilontol))THEN
  ishit = .FALSE.
  alpha = -1.0
  RETURN
END IF

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

IF(ABS(xi).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF

eta=(-A2*C1+A1*C2)*sdet
IF(ABS(eta).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF

isHit=.TRUE.

END SUBROUTINE ComputePlanarRectIntersection


SUBROUTINE ComputePlanarCurvedIntersection(isHit                       &
                                           ,PartTrajectory              &
                                           ,lengthPartTrajectory        &
                                           ,alpha                       &
                                           ,xi                          &
                                           ,eta                         &
                                           ,PartID                       &
                                           ,flip                        &
                                           ,SideID                      &
                                           ,opt_CriticalParllelInSide)
!===================================================================================================================================
! Compute the intersection with a planar non rectangular face
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Globals                ,ONLY: Cross,abort,UNIT_stdOut,CROSSNORM,UNITVECTOR
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Particle_Surfaces_Vars ,ONLY: SideNormVec,SideSlabNormals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars ,ONLY: locXi,locEta,locAlpha,SideDistance
USE MOD_Utils                  ,ONLY: InsertionSort
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: rBoundingBoxChecks
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Globals                ,ONLY: myrank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
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
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,p,q
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
LOGICAL                                  :: CriticalParallelInSide
REAL                                     :: XiNewton(2)
REAL                                     :: coeffA,locSideDistance
!REAL                                     :: Interval1D,dInterVal1D
! fallback algorithm
LOGICAL                                  :: failed
INTEGER(KIND=2)                          :: ClipMode
REAL                                     :: LineNormVec(1:2,1:2)
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: PartFaceAngle
!===================================================================================================================================
!PartTrajectory = PartTrajectory
! set alpha to minus 1, asume no intersection
alpha=-1.0
Xi   = 2.0
Eta  = 2.0
isHit=.FALSE.

#ifdef CODE_ANALYZE
rBoundingBoxChecks=rBoundingBoxChecks+1.
#endif /*CODE_ANALYZE*/

CriticalParallelInSide=.FALSE.

IF(DoRefMapping)THEN
  coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
  IF(coeffA.LE.0.)RETURN
  locSideDistance=SideDistance(SideID)-DOT_PRODUCT(LastPartPos(1:3,PartID),SideNormVec(1:3,SideID))
  locSideDistance=locSideDistance/coeffA
  IF(locSideDistance.GT.lengthPartTrajectory) RETURN
ELSE
  coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
  IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.
  IF(flip.EQ.0)THEN
    IF(coeffA.LE.0.)RETURN
    locSideDistance=SideDistance(SideID)-DOT_PRODUCT(LastPartPos(1:3,PartID),SideNormVec(1:3,SideID))
    locSideDistance=locSideDistance/coeffA
    IF(locSideDistance.GT.lengthPartTrajectory) RETURN
  ELSE
    IF(coeffA.GE.0.)RETURN
    locSideDistance=-SideDistance(SideID)+DOT_PRODUCT(LastPartPos(1:3,PartID),SideNormVec(1:3,SideID))
    locSideDistance=locSideDistance/coeffA
    IF(locSideDistance.GT.lengthPartTrajectory) RETURN
  END IF
END IF
IF(.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the

!IF(DoRefMapping)THEN
!  IF(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory).LT.0.)RETURN
!ELSE
!  ! dependend on master/slave flip
!  coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!  IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.
!  IF(flip.EQ.0)THEN
!    IF(coeffA.LT.0.)RETURN
!  ELSE
!    IF(coeffA.GT.0.)RETURN
!  END IF
!END IF
!! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
!IF(.NOT.InsideBoundingBox(LastPartPos(1:3,PartID),SideID))THEN ! the old particle position is not inside the bounding box
!  IF(.NOT.InsideBoundingBox(PartState(1:3,PartID),SideID))THEN ! the new particle position is not inside the bounding box
!    IF(.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the
!                                                                                              ! bounding box
!  END IF
!END IF

! 2.) Bezier intersection: transformation of bezier patch 3D->2D
!PartTrajectory = PartTrajectoryOrig + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF(ABS(PartTrajectory(3)).LT.0.)THEN
  n1=(/ -PartTrajectory(2)-PartTrajectory(3)  , PartTrajectory(1) ,PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3) , PartTrajectory(3) , -PartTrajectory(1)-PartTrajectory(2) /)
END IF

n1=UNITVECTOR(n1)
!n1=n1/SQRT(DOT_PRODUCT(n1,n1))
n2=CROSSNORM(PartTrajectory,n1)
!n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
!       , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
!       , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
!n2=n2/SQRT(DOT_PRODUCT(n2,n2))

DO q=0,NGeo
  DO p=0,NGeo
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n1)
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n2)
  END DO
END DO

XiNewton=0.
CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,PartID,SideID,failed)
! write Xinewton to locXi and locEta
locXi (1)=XiNewton(1)
locEta(1)=XiNewton(2)
IF(failed)THEN
  PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,2,SideID))))
  IPWRITE(UNIT_stdout,*) ' Intersection-angle-of-BezierNetwon: ',PartFaceAngle*180./PI
  iClipIter=0
  nXiClip=0
  nEtaClip=0
  nInterSections=0
  ClipMode=1
  LineNormVec=0.
  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory&
                ,iClipIter,nXiClip,nEtaClip,nInterSections,PartID,SideID)
  IF(nInterSections.GT.1)THEN
    nInterSections=1
  END IF
END IF

nInterSections=0
IF(locAlpha(1).GT.-1) nInterSections=1

IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
IF(CriticalParallelInSide)THEN
  IF(ALMOSTZERO(locAlpha(1))) THEN
    IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
  END IF
END IF

SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha=locAlpha(1)
  xi =locXi (1)
  eta=loceta(1)
  isHit=.TRUE.
  RETURN
CASE DEFAULT
  CALL abort(&
__STAMP__&
,' The code should never go here')
END SELECT


END SUBROUTINE ComputePlanarCurvedIntersection


SUBROUTINE ComputePlanarNonRectIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                                   ,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
! robust version
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_Utils                  ,ONLY: QuadraticSolver
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Particle_Surfaces_Vars ,ONLY: epsilontol,Beziercliphit
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID
!LOGICAL,INTENT(IN),OPTIONAL       :: ElemCheck_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: B,C,alphaNorm!,A
REAL                              :: xi(2),eta(2),t(2), scaleFac
INTEGER                           :: nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0
isHit=.FALSE.

! compute initial vectors
BiLinearCoeff(:,1) = 0.25*BaseVectors3(:,SideID)
BiLinearCoeff(:,2) = 0.25*BaseVectors1(:,SideID)
BiLinearCoeff(:,3) = 0.25*BaseVectors2(:,SideID)
BiLinearCoeff(:,4) = 0.25*BaseVectors0(:,SideID)

! compute product with particle trajectory
a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,iPart))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(3,iPart))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
a2(4)=(BilinearCoeff(2,4)-LastPartPos(2,iPart))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(3,iPart))*PartTrajectory(2)

!A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)

!scale with <PartTraj.,NormVec>^2 and cell-scale (~area) for getting coefficients at least approx. in the order of 1
scaleFac = DOT_PRODUCT(PartTrajectory,SideNormVec(1:3,SideID)) !both vectors are already normalized
IF(scaleFac.NE.0.)THEN
  scaleFac = scaleFac**2 * BaseVectorsScale(SideID) !<...>^2 * cell-scale
  !A = A d/ scaleFac
  B = B / scaleFac
  C = C / scaleFac
END IF

IF(ABS(B).GT.0.)THEN
  nRoot=1
  Eta(1)=-C/B
  Eta(2)=0.
ELSE
  nRoot=0
  Eta(1)=0.
  Eta(2)=0.
END IF
!CALL QuadraticSolver(A,B,C,nRoot,Eta(1),Eta(2))

IF(nRoot.EQ.0)THEN
  RETURN
END IF
IF (nRoot.EQ.1) THEN
  IF(ABS(eta(1)).LT.BezierClipHit)THEN
    ! check for Xi only, if eta is possible
    xi(1)=ComputeXi(a1,a2,eta(1))
    IF(ABS(xi(1)).LT.BezierClipHit)THEN
      ! compute alpha only with valid xi and eta
      t(1)=ComputeSurfaceDistance2(SideNormVec(1:3,SideID),BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      alphaNorm=t(1)/lengthPartTrajectory
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
      IF((alphaNorm.LE.1.0) .AND.(alphaNorm.GT.-epsilontol))THEN
        alpha=t(1)!/LengthPartTrajectory
        xitild=xi(1)
        etatild=eta(1)
        isHit=.TRUE.
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. OnePlusEps
  ELSE ! eta not in reange
    RETURN
  END IF ! eta .lt. OnePlusEps
END IF

END SUBROUTINE ComputePlanarNonRectIntersection


SUBROUTINE ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                      ,PartID,SideID,ElemCheck_Opt,alpha2)
!===================================================================================================================================
! Compute the Intersection with planar surface
! robust version
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Utils                  ,ONLY: QuadraticSolver
USE MOD_Particle_Vars          ,ONLY: LastPartPos
!USE MOD_Mesh_Vars              ,ONLY: nBCSides,nSides
!USE MOD_Particle_Surfaces_Vars  ,ONLY: Beziercliphit
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale,SideNormVec
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBilinear
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Surfaces_Vars ,ONLY: epsilontol
#endif /*CODE_ANALYZE*/
#if USE_MPI
!USE MOD_Mesh_Vars              ,ONLY: BC
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: PartID,SideID!,flip
LOGICAL,INTENT(IN),OPTIONAL       :: ElemCheck_Opt
REAL,INTENT(IN),OPTIONAL          :: alpha2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C,alphaNorm
REAL                              :: xi(2),eta(2),t(2), scaleFac!, n_loc(1:3)
INTEGER                           :: InterType,nRoot
LOGICAL                           :: ElemCheck
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0
isHit=.FALSE.

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
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | SideNormVec  : ',SideNormVec(1:3,SideID)
      WRITE(UNIT_stdout,'(A,4(X,G0))') '     | BilinearCoeff: ',BilinearCoeff(1,1:4)
      WRITE(UNIT_stdout,'(A,4(X,G0))') '     | BilinearCoeff: ',BilinearCoeff(2,1:4)
      WRITE(UNIT_stdout,'(A,4(X,G0))') '     | BilinearCoeff: ',BilinearCoeff(3,1:4)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint1: ',BezierControlPoints3D(:,0,0,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint2: ',BezierControlPoints3D(:,NGeo,0,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint3: ',BezierControlPoints3D(:,0,NGeo,SideID)
      WRITE(UNIT_stdout,'(A,3(X,G0))') '     | Beziercontrolpoint4: ',BezierControlPoints3D(:,NGeo,NGeo,SideID)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! compute product with particle trajectory
!a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
!a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
!a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
!a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(3) &
!     -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(1)

!a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
!a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
!a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
!a2(4)=(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(3) &
!     -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(2)
!-- use the following instead of the previous lines, since otherwise coeffs can be cancelled out for trajectories with 0-components
IF((ABS(SideNormVec(1,SideID)).GE.ABS(SideNormVec(2,SideID))) .AND.(ABS(SideNormVec(1,SideID)).LE.ABS(SideNormVec(3,SideID))) &
  .AND. .NOT.ALMOSTZERO(PartTrajectory(1)))THEN
  a1(1)= BilinearCoeff(2,1)*PartTrajectory(1) - BilinearCoeff(1,1)*PartTrajectory(2)
  a1(2)= BilinearCoeff(2,2)*PartTrajectory(1) - BilinearCoeff(1,2)*PartTrajectory(2)
  a1(3)= BilinearCoeff(2,3)*PartTrajectory(1) - BilinearCoeff(1,3)*PartTrajectory(2)
  a1(4)=(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(1) &
       -(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(2)

  a2(1)= BilinearCoeff(3,1)*PartTrajectory(1) - BilinearCoeff(1,1)*PartTrajectory(3)
  a2(2)= BilinearCoeff(3,2)*PartTrajectory(1) - BilinearCoeff(1,2)*PartTrajectory(3)
  a2(3)= BilinearCoeff(3,3)*PartTrajectory(1) - BilinearCoeff(1,3)*PartTrajectory(3)
  a2(4)=(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(1) &
       -(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(3)
ELSE IF(ABS(SideNormVec(2,SideID)).LE.ABS(SideNormVec(3,SideID)) &
  .AND. .NOT.ALMOSTZERO(PartTrajectory(1)))THEN
  a1(1)= BilinearCoeff(1,1)*PartTrajectory(2) - BilinearCoeff(2,1)*PartTrajectory(1)
  a1(2)= BilinearCoeff(1,2)*PartTrajectory(2) - BilinearCoeff(2,2)*PartTrajectory(1)
  a1(3)= BilinearCoeff(1,3)*PartTrajectory(2) - BilinearCoeff(2,3)*PartTrajectory(1)
  a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(2) &
       -(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(1)

  a2(1)= BilinearCoeff(3,1)*PartTrajectory(2) - BilinearCoeff(2,1)*PartTrajectory(3)
  a2(2)= BilinearCoeff(3,2)*PartTrajectory(2) - BilinearCoeff(2,2)*PartTrajectory(3)
  a2(3)= BilinearCoeff(3,3)*PartTrajectory(2) - BilinearCoeff(2,3)*PartTrajectory(3)
  a2(4)=(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(2) &
       -(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(3)
ELSE IF(.NOT.ALMOSTZERO(PartTrajectory(3)))THEN
  a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
  a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
  a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
  a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(3) &
       -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(1)

  a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
  a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
  a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
  a2(4)=(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(3) &
       -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(2)
!if PartTrajectory should be zero in largest component of SideNormVec, decide based on original check:
ELSE IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(3))))THEN
  a1(1)= BilinearCoeff(2,1)*PartTrajectory(1) - BilinearCoeff(1,1)*PartTrajectory(2)
  a1(2)= BilinearCoeff(2,2)*PartTrajectory(1) - BilinearCoeff(1,2)*PartTrajectory(2)
  a1(3)= BilinearCoeff(2,3)*PartTrajectory(1) - BilinearCoeff(1,3)*PartTrajectory(2)
  a1(4)=(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(1) &
       -(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(2)

  a2(1)= BilinearCoeff(3,1)*PartTrajectory(1) - BilinearCoeff(1,1)*PartTrajectory(3)
  a2(2)= BilinearCoeff(3,2)*PartTrajectory(1) - BilinearCoeff(1,2)*PartTrajectory(3)
  a2(3)= BilinearCoeff(3,3)*PartTrajectory(1) - BilinearCoeff(1,3)*PartTrajectory(3)
  a2(4)=(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(1) &
       -(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(3)
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  a1(1)= BilinearCoeff(1,1)*PartTrajectory(2) - BilinearCoeff(2,1)*PartTrajectory(1)
  a1(2)= BilinearCoeff(1,2)*PartTrajectory(2) - BilinearCoeff(2,2)*PartTrajectory(1)
  a1(3)= BilinearCoeff(1,3)*PartTrajectory(2) - BilinearCoeff(2,3)*PartTrajectory(1)
  a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(2) &
       -(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(1)

  a2(1)= BilinearCoeff(3,1)*PartTrajectory(2) - BilinearCoeff(2,1)*PartTrajectory(3)
  a2(2)= BilinearCoeff(3,2)*PartTrajectory(2) - BilinearCoeff(2,2)*PartTrajectory(3)
  a2(3)= BilinearCoeff(3,3)*PartTrajectory(2) - BilinearCoeff(2,3)*PartTrajectory(3)
  a2(4)=(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(2) &
       -(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(3)
ELSE
  a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
  a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
  a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
  a1(4)=(BilinearCoeff(1,4)-LastPartPos(1,PartID))*PartTrajectory(3) &
       -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(1)

  a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
  a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
  a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
  a2(4)=(BilinearCoeff(2,4)-LastPartPos(2,PartID))*PartTrajectory(3) &
       -(BilinearCoeff(3,4)-LastPartPos(3,PartID))*PartTrajectory(2)
END IF

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A,4(X,G0))') '     | a1: ',a1
      WRITE(UNIT_stdout,'(A,4(X,G0))') '     | a2: ',a2
      WRITE(UNIT_stdout,'(A)') '     | Quadratic equation constants: '
      WRITE(UNIT_stdout,'(3(A,G0))') '     | A: ',A,' | B: ',B,' | C: ',C
    END IF
  END IF
#endif /*CODE_ANALYZE*/

!scale with <PartTraj.,NormVec>^2 and cell-scale (~area) for getting coefficients at least approx. in the order of 1
scaleFac = DOT_PRODUCT(PartTrajectory,SideNormVec(1:3,SideID)) !both vectors are already normalized
IF(scaleFac.NE.0.)THEN
  scaleFac = scaleFac**2 * BaseVectorsScale(SideID) !<...>^2 * cell-scale
  scaleFac = 1./scaleFac
  A = A * scaleFac
  B = B * scaleFac
  C = C * scaleFac
END IF

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | Quadratic equation constants (after scaling): '
      WRITE(UNIT_stdout,'(3(A,G0))') '     | A: ',A,' | B: ',B,' | C: ',C
    END IF
  END IF
#endif /*CODE_ANALYZE*/

CALL QuadraticSolver(A,B,C,nRoot,Eta(1),Eta(2))

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | Output after QuadraticSolver: '
      WRITE(UNIT_stdout,'(A,I0,A,2(X,G0))') '     | number of root: ',nRoot,' | Eta: ',Eta(1:2)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | nRoot = 1 '
    END IF
  END IF
#endif /*CODE_ANALYZE*/

  IF(ABS(eta(1)).LE.1.0) THEN!.LT.BezierClipHit)THEN
    ! check for Xi only, if eta is possible
    xi(1)=ComputeXi(a1,a2,eta(1))
    IF(ABS(xi(1)).LE.1.0) THEN!.LT.BezierClipHit)THEN
      ! compute alpha only with valid xi and eta
      t(1)=ComputeSurfaceDistance2(SideNormVec(1:3,SideID),BiLinearCoeff,xi(1),eta(1),PartTrajectory,PartID)
      IF (PRESENT(alpha2)) THEN
        IF (alpha2.GT.-1.0) THEN
          IF (ALMOSTEQUAL(t(1),alpha2)) THEN
            t(1)=-1.0
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(PartID.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A)') 'changed t1'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
          END IF
        END IF
      END IF
      alphaNorm=t(1)/lengthPartTrajectory
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,G0,A,G0,A,G0)') '     | xi: ',xi(1),' | t: ',t(1),' | alphaNorm: ',alphaNorm
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
      IF((alphaNorm.LE.1.0) .AND.(alphaNorm.GE.0.))THEN!.GT.-epsilontol))THEN
        alpha=t(1)!/LengthPartTrajectory
        xitild=xi(1)
        etatild=eta(1)
        isHit=.TRUE.
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,G0,A,G0)') '     | alphanorm: ',alphaNorm,' | epsilonTolerance: ',epsilontol
        END IF
      END IF
#endif /*CODE_ANALYZE*/
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. OnePlusEps
  ELSE ! eta not in reange
    RETURN
  END IF ! eta .lt. OnePlusEps
ELSE
  InterType=0
  t(:)=-1.

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | nRoot = 2 '
    END IF
  END IF
#endif /*CODE_ANALYZE*/

  ! check if intersection is possible
  ! t(1) has to be nullified if intersection is NOT possible
  ! else, the selection scheme is WRONG
  IF(ABS(eta(1)).LE.1.0) THEN!.LT.BezierClipHit)THEN
    ! check for Xi only, if eta is possible
    xi(1)=ComputeXi(a1,a2,eta(1))
    IF(ABS(xi(1)).LE.1.0) THEN!.LT.BezierCliphit)THEN
      ! compute alpha only with valid xi and eta
      t(1)=ComputeSurfaceDistance2(SideNormVec(1:3,SideID),BiLinearCoeff,xi(1),eta(1),PartTrajectory,PartID)
      IF (PRESENT(alpha2)) THEN
        IF (alpha2.GT.-1.0) THEN
          IF (ALMOSTEQUAL(t(1),alpha2)) THEN
            t(1)=-1.0
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(PartID.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A)') 'changed t1'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
          END IF
        END IF
      END IF
      alphaNorm=t(1)/lengthPartTrajectory
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,G0,A,G0,A,G0)') '     | xi: ',xi(1),' | t: ',t(1),' | alphaNorm: ',alphaNorm
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GE.0.))THEN
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
      IF((alphaNorm.LE.1.0) .AND.(alphaNorm.GE.0.))THEN!.GT.-epsilontol))THEN
        InterType=InterType+1
        isHit=.TRUE.
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | alphanorm1: ',alphaNorm,' | epsilonTolerance: ',epsilontol
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      END IF
    END IF
  END IF ! eta(1)


  ! check if intersection is possible
  ! t(2) has to be nullified if intersection is NOT possible
  ! else, the selection scheme is WRONG
  IF(ABS(eta(2)).LE.1.0) THEN!.LT.BezierClipHit)THEN
    ! check for Xi only, if eta is possible
    xi(2)=ComputeXi(a1,a2,eta(2))
    IF(ABS(xi(2)).LE.1.0) THEN!.LT.BezierClipHit)THEN
      ! compute alpha only with valid xi and eta
      t(2)=ComputeSurfaceDistance2(SideNormVec(1:3,SideID),BiLinearCoeff,xi(2),eta(2),PartTrajectory,PartID)
      IF (PRESENT(alpha2)) THEN
        IF (alpha2.GT.-1.0) THEN
          IF (ALMOSTEQUAL(t(2),alpha2)) THEN
            t(2)=-1.0
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(PartID.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A)') 'changed t2'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
          END IF
        END IF
      END IF
      alphaNorm=t(2)/lengthPartTrajectory
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(PartID.EQ.PARTOUT)THEN
          WRITE(UNIT_stdout,'(A,G0,A,G0,A,G0)') '     | xi: ',xi(2),' | t: ',t(2),' | alphaNorm: ',alphaNorm
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      IF((alphaNorm.LT.1.0) .AND.(alphaNorm.GE.0.))THEN!.GT.-epsilontol))THEN
        ! Two solutions can be correspond to one unique intersection (?!)
        IF(InterType.EQ.1)THEN
          IF(.NOT.ALMOSTEQUALRELATIVE(t(2),t(1),1e-8))THEN
            isHit=.TRUE.
            InterType=InterType+2
          END IF
        ELSE
          isHit=.TRUE.
          InterType=InterType+2
        END IF
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | alphanorm2: ',alphaNorm,' | epsilonTolerance: ',epsilontol
          END IF
        END IF
#endif /*CODE_ANALYZE*/
      END IF
    END IF
  END IF

  IF(InterType.EQ.0) THEN
    RETURN
  END IF
  isHit=.TRUE.
  IF(DoRefMapping) THEN
    SELECT CASE(InterType)
    CASE(1)
      alpha  =t  (1)
      xitild =xi (1)
      etatild=eta(1)
    CASE(2)
       alpha  =t  (2)
       xitild =xi (2)
       etatild=eta(2)
    CASE DEFAULT
     ! two intersections
      IF(t(1).LT.t(2))THEN
        alpha  =t  (1)
        xitild =xi (1)
        etatild=eta(1)
      ELSE
        alpha  =t  (2)
        xitild =xi (2)
        etatild=eta(2)
      END IF
    END SELECT
    RETURN
  END IF
  ! no refmapping

  ! TODO: this is obsolete with new halo region
!  IF(SideID.LE.nSides)THEN
!    IF(SideID.LE.nBCSides)THEN
      ! take closest
      SELECT CASE(InterType)
      CASE(1)
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      CASE(2)
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      CASE(3)
        ElemCheck = .FALSE.
        IF(PRESENT(ElemCheck_Opt))THEN
          ElemCheck = ElemCheck_Opt
        END IF
        IF(ElemCheck)THEN
          alpha = -1
          xitild = -2
          etatild = -2
        ELSE
          IF(ABS(t(1)).LT.ABS(t(2)))THEN
            alpha=t(1)
            xitild=xi(1)
            etatild=eta(1)
          ELSE
            alpha=t(2)
            xitild=xi(2)
            etatild=eta(2)
          END IF
        END IF
      END SELECT
!    ELSE
!      SELECT CASE(InterType)
!      CASE(1)
!        alpha=t(1)
!        xitild=xi(1)
!        etatild=eta(1)
!      CASE(2)
!        alpha=t(2)
!        xitild=xi(2)
!        etatild=eta(2)
!      CASE(3) ! double intersection leaves and entries element
!!        IF(ABS(t(1)).LT.ABS(t(2)))THEN
!!          CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi(2),eta=eta(2),SideID=SideID)
!!          IF(flip.NE.0) n_loc=-n_loc
!!          IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0)THEN
!!            alpha=t(2)
!!            xitild=xi(2)
!!            etatild=eta(2)
!!          ELSE
!            alpha=-1.0
!            xitild=0.
!            etatild=0.
!            isHit=.FALSE.
!!          END IF
!!        ELSE
!!          CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi(1),eta=eta(1),SideID=SideID)
!!          IF(flip.NE.0) n_loc=-n_loc
!!          IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0)THEN
!!            alpha=t(1)
!!            xitild=xi(1)
!!            etatild=eta(1)
!!          ELSE
!!            alpha=-1.0
!!            xitild=0.
!!            etatild=0.
!!            isHit=.FALSE.
!!          END IF
!!        END IF
!      END SELECT
!    END IF
!#if USE_MPI
!  ELSE
!    ! halo side
!    IF(BC(SideID).GT.0)THEN ! BC Sides
!      ! take closest
!      SELECT CASE(InterType)
!      CASE(1)
!        alpha=t(1)
!        xitild=xi(1)
!        etatild=eta(1)
!      CASE(2)
!        alpha=t(2)
!        xitild=xi(2)
!        etatild=eta(2)
!      CASE(3)
!        ElemCheck = .FALSE.
!        IF(PRESENT(ElemCheck_Opt))THEN
!          ElemCheck = ElemCheck_Opt
!        END IF
!        IF(ElemCheck)THEN
!          alpha = -1
!          xitild = -2
!          etatild = -2
!        ELSE
!          IF(ABS(t(1)).LT.ABS(t(2)))THEN
!            alpha=t(1)
!            xitild=xi(1)
!            etatild=eta(1)
!          ELSE
!            alpha=t(2)
!            xitild=xi(2)
!            etatild=eta(2)
!          END IF
!        END IF
!      END SELECT
!    ELSE
!      SELECT CASE(InterType)
!      CASE(1)
!        alpha=t(1)
!        xitild=xi(1)
!        etatild=eta(1)
!      CASE(2)
!        alpha=t(2)
!        xitild=xi(2)
!        etatild=eta(2)
!      CASE(3) ! double intersection leaves and entries element
!        alpha=-1.0
!        xitild=0.
!        etatild=0.
!        isHit=.FALSE.
!      END SELECT
!    END IF
!#endif /*USE_MPI*/
!  END IF
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersection


SUBROUTINE ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID &
                ,SideID,opt_CriticalParllelInSide,ElemCheck_Opt)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,            ONLY:PI
USE MOD_Globals,                 ONLY:Cross,abort,CROSSNORM,UNITVECTOR
USE MOD_Mesh_Vars,               ONLY:NGeo,BC
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,BezierNewtonAngle
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,  ONLY:locXi,locEta,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals!,epsilonTol
USE MOD_Utils,                   ONLY:InsertionSort !BubbleSortID
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipLocalTol
USE MOD_Particle_Surfaces,       ONLY:CalcNormAndTangBezier
#ifdef CODE_ANALYZE
USE MOD_Globals,                 ONLY:MyRank,UNIT_stdOut
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces_Vars,  ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints

#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)           :: PartTrajectory
REAL,INTENT(IN)                          :: lengthPartTrajectory
INTEGER,INTENT(IN)                       :: PartID,SideID
LOGICAL,INTENT(IN),OPTIONAL              :: ElemCheck_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                         :: alpha,xi,eta
LOGICAL,INTENT(OUT)                      :: isHit
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,iInter,p,q
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
alpha=-1.0
Xi   = 2.0
Eta  = 2.0
isHit=.FALSE.

#ifdef CODE_ANALYZE
rBoundingBoxChecks=rBoundingBoxChecks+1.
#endif /*CODE_ANALYZE*/

CriticalParallelInSide=.FALSE.
IF(BoundingBoxIsEmpty(SideID))THEN
  IF(DoRefMapping)THEN
    IF(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory).LT.0.)RETURN
  ELSE
    IF(ALMOSTZERO(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory))) CriticalParallelInSide=.TRUE.
  END IF
  IF(.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the
                                                                                                ! bounding box
ELSE
  ! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
  IF(.NOT.InsideBoundingBox(LastPartPos(1:3,PartID),SideID))THEN ! the old particle position is not inside the bounding box
    IF(.NOT.InsideBoundingBox(PartState(1:3,PartID),SideID))THEN ! the new particle position is not inside the bounding box
      IF(.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)) RETURN ! the particle does not intersect the
                                                                                                ! bounding box
    END IF
  END IF
END IF
! 2.) Bezier intersection: transformation of bezier patch 3D->2D
!PartTrajectoryOrig=PartTrajectory
!PartTrajectory = PartTrajectory + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF(ABS(PartTrajectory(3)).LT.0.)THEN
  n1=(/ -PartTrajectory(2)-PartTrajectory(3)  , PartTrajectory(1) ,PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3) , PartTrajectory(3) , -PartTrajectory(1)-PartTrajectory(2) /)
END IF

! check angle to boundingbox (height normal vector)
PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,2,SideID))))
IF(ALMOSTZERO(PartFaceAngle*180/ACOS(-1.)))THEN
  n1=n1! +epsilontol
END IF

!n1=n1/SQRT(DOT_PRODUCT(n1,n1))
!n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
!       , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
!       , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
!n2=n2/SQRT(DOT_PRODUCT(n2,n2))
n1=UNITVECTOR(n1)
n2=CROSSNORM(PartTrajectory,n1)
!PartTrajectory = PartTrajectoryOrig !set back for preventing angles > 90 deg (0.5pi+eps)

! projection like Nishita
! plane 1 with n1 becomes y-axis and plane 2 with n2 becomes the x-axis
DO q=0,NGeo
  DO p=0,NGeo
    ! n2 is perpendicular to x-axis => gives distance to new x-axis
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n2)
    ! n1 is perpendicular to y-axis => gives distance to new y-axis
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(1:3,PartID),n1)
  END DO
END DO

! calculate angle between particle path and slab normal plane of face
! angle2=abs(90-RadianToDegree*acos(scalar_product/(VECTOR_LENGTH(t)*VECTOR_LENGTH(n2))))
!IF(BoundingBoxIsEmpty(SideID))THEN
!  PartFaceAngle=BezierNewtonAngle+1.0
!ELSE
!END IF
!IF(.NOT.BezierNewtonAngle)THEN
IF((PartFaceAngle.LT.BezierNewtonAngle))THEN ! 1° = 0.01745rad: critical side at the moment need: 0.57° angle
#ifdef CODE_ANALYZE
rPerformBezierClip=rPerformBezierClip+1.
#endif /*CODE_ANALYZE*/
  !  this part in a new function or subroutine
  locAlpha=-1.0
  iClipIter=0
  nXiClip=0
  nEtaClip=0
  nInterSections=0
  ! check extend in Xi  and eta direction
  dXi =MAXVAL(BezierControlPoints2D(1,:,:))-MINVAL(BezierControlPoints2D(1,:,:))
  dEta=MAXVAL(BezierControlPoints2D(2,:,:))-MINVAL(BezierControlPoints2D(2,:,:))
  IF(dXi.GT.dEta)THEN
    ! first Clip is in XI diretion
    ClipMode=1
  ELSE
    ! first Clip is in ETA diretion
    ClipMode=2
  END IF
  BezierClipLocalTol=MIN(dXi,dEta)*BezierClipTolerance
  LineNormVec=0.
  ! CALL recursive Bezier clipping algorithm
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
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
ELSE!BezierNewtonAngle
#ifdef CODE_ANALYZE
rPerformBezierNewton=rPerformBezierNewton+1.
BezierControlPoints2D_tmp=BezierControlPoints2D
locAlpha=-1.0
iClipIter=0
nXiClip=0
nEtaClip=0
nInterSections=0
ClipMode=1
LineNormVec=0.
CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_tmp,LineNormVec,PartTrajectory,lengthPartTrajectory&
              ,iClipIter,nXiClip,nEtaClip,nInterSections,PartID,SideID)
IF(nInterSections.GT.1)THEN
  CALL abort(&
  __STAMP__&
  ,' More then one intersection! Cannot use Newton!' ,nInterSections)
END IF
#endif /*CODE_ANALYZE*/
  !dInterVal1D =MINVAL(BezierControlPoints2D(1,:,:))
  !InterVal1D  =MAXVAL(BezierControlPoints2D(1,:,:))-dInterVal1D
  !XiNewton(1) =MIN(-1.0+ABS(dInterVal1D)/InterVal1D,1.0)
  !InterVal1D  =MINVAL(BezierControlPoints2D(2,:,:))
  !InterVal1D  =MAXVAL(BezierControlPoints2D(2,:,:))-dInterVal1D
  !XiNewton(2) =MIN(-1.0+ABS(dInterVal1D)/InterVal1D,1.0)
  XiNewton=0.
  !CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D_tmp,PartTrajectory,lengthPartTrajectory,PartID,SideID)
  CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,PartID,SideID,failed)
  nInterSections=0
  IF(locAlpha(1).GT.-1) nInterSections=1
  IF(failed) CALL abort(&
    __STAMP__&
    ,' Bezier-Newton does not yield root! ')
#ifdef CODE_ANALYZE
  IF(nInterSections.EQ.1)THEN
    dXi =ABS(locXi(1)-XiNewton(1)) !/(400.*BezierClipTolerance)
    dEta=ABS(locEta(1)-XiNewton(2))!/(400.*BezierClipTolerance)
    dXi =dXi*dXi+dEta*dEta
    dXi =SQRT(dXi)/(400.*BezierClipTolerance)
    IF(dXi.GT.1.0) THEN
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
  locXi (1)=XiNewton(1)
  locEta(1)=XiNewton(2)
END IF

IF(PRESENT(ElemCheck_Opt))THEN
  IF(ElemCheck_Opt)THEN
    DO iInter=1,nInterSections
      CALL CalcNormAndTangBezier(nVec=n1,xi=locXi(iInter),eta=locEta(iInter),SideID=SideID)
      IF(ABS(DOT_PRODUCT(PartTrajectory,n1)).LT.0.125) THEN
        locAlpha(iInter)=-1
        nInterSections=nInterSections-1
      END IF
    END DO ! iInter=2,nInterSections
  END IF
END IF

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
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
  alpha=locAlpha(1)
  xi =locXi (1)
  eta=loceta(1)
  isHit=.TRUE.
  IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
  IF(CriticalParallelInSide)THEN
    IF(ALMOSTZERO(alpha))THEN
      IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
    END IF
  END IF
  RETURN
CASE DEFAULT
  ! more than one intersection
  !ALLOCATE(locID(nInterSections))
  ALLOCATE(locID(nInterSections))
  DO iInter=1,nInterSections
    locID(iInter)=iInter
  END DO ! iInter
  ! sort intersection distance
!  CALL BubbleSortID(locAlpha,locID,nIntersections)
  CALL InsertionSort(locAlpha(1:nIntersections),locID,nIntersections)
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      IPWRITE(UNIT_stdout,*) ' locAlpha-sorted ',locAlpha(1:nIntersections)
    END IF
  END IF
#endif /*CODE_ANALYZE*/
  IF(DoRefMapping)THEN
    DO iInter=1,nInterSections
      IF(locAlpha(iInter).GT.-1.0)THEN
        alpha=locAlpha(iInter)
        xi =locXi (locID(iInter))
        eta=loceta(locID(iInter))
        DEALLOCATE(locID)
        isHit=.TRUE.
        RETURN
      END IF
    END DO ! iInter
  ELSE
    ! no ref mapping
    ! get real number of intersections
    realnInter=1
    ALLOCATE(realInterID(1:nInterSections))
    realInterID=0
    realInterID(1)=1
    ! PO & CS:
    ! we used the approach to check the previous (i-1)  with the current (i) alpha, if they
    ! are almost identically, it is ignored (multiple intersections are reduced to one)
    ! second possibility:
    ! check only to the accepted alphas
    DO iInter=2,nInterSections
      IF(.NOT.ALMOSTEQUALRELATIVE(locAlpha(iInter-1),locAlpha(iInter),0.002))THEN
        realNInter=realNInter+1
        realInterID(realNInter)=iInter
        isInter=iInter
      END IF
    END DO ! iInter=2,nInterSections
#ifdef CODE_ANALYZE
     IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
       IF(PartID.EQ.PARTOUT)THEN
         IPWRITE(UNIT_stdout,*) ' realnInter ',realnInter
       END IF
     END IF
#endif /*CODE_ANALYZE*/
    IF(BC(SideID).GT.0)THEN
      IF(PRESENT(ElemCheck_Opt))THEN
        IF(ElemCheck_Opt)THEN
          IF(MOD(realNInter,2).EQ.0) THEN
            alpha=-1
            nInterSections=0
            RETURN
          END IF
        END IF
      END IF
      ! boundary side, take first intersection
      alpha=locAlpha(1)
      xi =locXi (locID(1))
      eta=loceta(locID(1))
      DEALLOCATE(locID)
      DEALLOCATE(realInterID)
      isHit=.TRUE.
      IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
      IF(CriticalParallelInSide)THEN
        IF(ALMOSTZERO(alpha))THEN
          IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
        END IF
      END IF
      RETURN
    ELSE
      IF(MOD(realNInter,2).EQ.0) THEN
        ! particle leaves and enters cell multiple times, however, remain
        ! still inside of the element
        DEALLOCATE(locID)
        DEALLOCATE(realInterID)
        alpha=-1.0
        isHit=.FALSE.
        IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
        IF(CriticalParallelInSide)THEN
          IF(ALMOSTZERO(alpha))THEN
            IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
          END IF
        END IF
        RETURN
      ELSE
        ! particle leaves and enters, take the LAST intersection
        alpha=locAlpha(realInterID(realNInter))
        xi =locXi (locID(realInterID(realNInter)))
        eta=loceta(locID(realInterID(realNInter)))
        isHit=.TRUE.
        DEALLOCATE(locID)
        DEALLOCATE(realInterID)
        IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
        IF(CriticalParallelInSide)THEN
          IF(ALMOSTZERO(alpha))THEN
            IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
          END IF
        END IF
        RETURN
      END IF
    END IF ! BC or no BC side
  END IF
  SDEALLOCATE(locID)
END SELECT

CALL abort(&
__STAMP__&
,' The code should never go here')

! remove compiler warning
IF(PRESENT(ElemCheck_Opt)) ClipMode=0

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
USE MOD_Globals,                 ONLY:Abort
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipMaxIter!,BezierClipTolerance
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipLineVectorMethod
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Globals,                 ONLY:UNIT_stdOut
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
! INTEGER,INTENT(INOUT),DIMENSION(:)   :: locID
INTEGER,INTENT(INOUT)                  :: iClipIter,nXiClip,nEtaClip,nInterSections
INTEGER(KIND=2),INTENT(INOUT)          :: ClipMode
REAL,DIMENSION(2,2),INTENT(INOUT)      :: LineNormVec
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: PatchDOF2D
!REAL                                 :: ZeroDistance,BezierClipTolerance2
!================================================================================================================================

PatchDOF2D=1.0/REAL((NGeo+1)*(NGeo+1))

! 3.) Bezier intersection: solution Newton's method or Bezier clipping
! outcome: no intersection, single intersection, multiple intersection with patch
DO WHILE(iClipIter.LE.BezierClipMaxIter)
  IF(iClipIter.EQ.0)THEN
    !IF(BezierClipLineVectorMethod.EQ.0) CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec(:,:),NGeo,0)
    IF(BezierClipLineVectorMethod.EQ.0) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
  END IF
  iClipIter=iClipIter+1
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A,I0,X,I0)') ' iClipIter,ClipMode ', iClipIter, ClipMode
      !read*
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
    !IF(BezierClipLineVectorMethod.EQ.1) CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec(:,:),NGeo,0)
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
    !IF(BezierClipLineVectorMethod.EQ.1) CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec(:,:),NGeo,0)
    IF(BezierClipLineVectorMethod.EQ.1) CALL CalcLineNormVec3(BezierControlPoints2D(:,:,:),LineNormVec(:,:))
    IF(BezierClipLineVectorMethod.EQ.2) THEN
      !IF(MOD(iClipIter,2).EQ.1) CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec(:,:),NGeo,0)
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
!USE MOD_Particle_Mesh_Vars,      ONLY:PartSideToElem
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
!  IPWRITE(UNIT_stdout,*) ' SideId      : ', SideID
!  IPWRITE(UNIT_stdout,*) ' PartID      : ', PartID
!  IPWRITE(UNIT_stdout,*) ' ElemID      : ', PartSideToElem(S2E_ELEM_ID,SideID)
!  IPWRITE(UNIT_stdout,*) ' Norm_P      : ', Norm_P
!  IPWRITE(UNIT_stdout,*) ' minmax-1    : ', MinMax(:,1)
!  IPWRITE(UNIT_stdout,*) ' minmax-2    : ', MinMax(:,2)
!  IPWRITE(UNIT_stdout,*) ' xi, eta     : ', xi
!  IPWRITE(UNIT_stdout,*) ' dxi, dxi2   : ', dXi, dXi2
!  IPWRITE(UNIT_stdout,*) ' PartState   : ', PartState(1:3,PartID)
!  IPWRITE(UNIT_stdout,*) ' lastPos     : ', LastPartPos(1:3,PartID)
!  IPWRITE(UNIT_stdout,*) ' Trajectory  : ', PartTrajectory
!  IPWRITE(UNIT_stdout,*) ' Calling-Bezier-Clipping  '
  failed=.TRUE.
  RETURN
!  CALL abort(&
!    __STAMP__&
!    ,' Bezier-Newton does not yield root! ')
END IF

! check if found Xi,Eta are in parameter range
!IF(ABS(xi(1)).GT.BezierNewtonHit) RETURN
!IF(ABS(xi(2)).GT.BezierNewtonHit) RETURN
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

!IF((alphaNorm.LE.BezierNewtonHit).AND.(alphaNorm.GT.-epsilontol)) RETURN
IF((alphaNorm.LE.1.0).AND.(alphaNorm.GT.-epsilontol)) RETURN
alpha=-1.0

END SUBROUTINE BezierNewton


SUBROUTINE calcLineNormVec2(BezierControlPoints2D,LineNormVec,a,b)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
!================================================================================================================================
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Mesh_Vars,               ONLY:NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: a ! NGeo
INTEGER,INTENT(IN)                   :: b ! 0
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: LineNormVec(1:2,1:2)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length,doPro,dcorr!,alpha(2),dalpha
REAL,DIMENSION(2)                    :: LXi, Leta, MBar
REAL                                 :: LineNormVecOld(1:2,1:2)
!================================================================================================================================

! backup old linenormvec
!LineNormVecOld=LineNormVec
LXi=(BezierControlPoints2D(:,a,b)-BezierControlPoints2D(:,0,0))+&
    (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,b,a))
Length=SQRT(DOT_PRODUCT(LXi,LXi))

IF(Length.EQ.0)THEN
  LineNormVec(:,1)=(/1.,0./)
  LineNormVec(:,2)=(/0.,1./)
  RETURN
END IF
LXi=LXi/Length

Leta=(BezierControlPoints2D(:,b,a)-BezierControlPoints2D(:,0,0))+&
     (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,a,b))
Length=SQRT(DOT_PRODUCT(Leta,Leta))
IF(Length.EQ.0)THEN
  LineNormVec(:,1)=(/1.,0./)
  LineNormVec(:,2)=(/0.,1./)
  RETURN
END IF
Leta=Leta/Length

LineNormVecOld(:,1)=Lxi
LineNormVecOld(:,2)=Leta

! if Lxi and Leta are orientated in opposite directions
doPro=DOT_PRODUCT(Lxi,Leta) ! can be negative

! efremov200501 step1
!IF(doPro.LT.0) Lxi=Lxi

IF(ABS(doPro).GT.0.6)THEN
!  print*,'doPro',doPro,ACOS(doPro)*180/pi

!  ! only required here
!  alpha(1)=atan2(Lxi (2),Lxi (1))
!  alpha(2)=atan2(Leta(2),Leta(1))
!
!  IF(alpha(1).LT.0) alpha(1)=alpha(1)+2.*PI
!  IF(alpha(2).LT.0) alpha(2)=alpha(2)+2.*PI
!
!  dalpha=alpha(1)-alpha(2)
!
!  !print*,'dalpha',dalpha
!  IF(dalpha.LT.0)THEN
!    print*,'dalpha<0'
!    ! alpha2.GT.alpha1
!    IF(dalpha.LT.-PI)THEN
!      print*,'<-pi'
!      IF(dalpha.LT.-1.5*PI)THEN
!        print*,'<-1.5pi'
!        ! dcorr=2*PI+dalpha
!        dcorr=2*PI+dalpha
!
!      ELSE
!        print*,'ese'
!        dcorr=dalpha+PI
!      END IF
!    ELSE
!      IF(dalpha.LT.-0.5*PI)THEN
!        ! dcorr=2*PI+dalpha
!        dcorr=PI+dalpha
!      ELSE
!        dcorr=dalpha !+0.5*PI
!      END IF
!    END IF
!  ELSE
!    ! angle angle3.LT.angle1
!    print*,'dalpha>0'
!    ! how can the angle be neg.
!    IF(dalpha.GT.PI)THEN
!      IF(dalpha.GT.1.5*PI)THEN
!        dcorr=-2*PI+dalpha
!      ELSE
!        dcorr=-PI+dalpha
!      END IF
!    ELSE
!      IF(dalpha.GT.0.5*PI)THEN
!        dcorr=-PI+dalpha
!      ELSE
!        dcorr=dalpha
!      END IF
!    END IF
!  END IF
!  dcorr=dcorr*0.5
!  alpha(1)=alpha(1)+dcorr
!  LineNormVec(1,1)=COS(alpha(1))
!  LineNormVec(2,1)=SIN(alpha(1))
!  alpha(2)=alpha(2)-dcorr
!  LineNormVec(1,2)=COS(alpha(2))
!  LineNormVec(2,2)=SIN(alpha(2))
  MBar=0.5*(Lxi+Leta)
  dCorr=0.5*SQRT(3.)
  IF(doPro.LT.0.)THEN
    ! Lxi
    LineNormVec(1,1)= 0.5*Mbar(2) + dCorr*Mbar(1)
    LineNormVec(2,1)=-0.5*Mbar(1) + dCorr*Mbar(2)
    ! Leta
    LineNormVec(1,2)=-0.5*Mbar(2) + dCorr*Mbar(1)
    LineNormVec(2,2)= 0.5*Mbar(1) + dCorr*Mbar(2)
  ELSE
    ! Lxi
    LineNormVec(1,1)=-0.5*Mbar(2) + dCorr*Mbar(1)
    LineNormVec(2,1)= 0.5*Mbar(1) + dCorr*Mbar(2)
    ! Leta
    LineNormVec(1,2)= 0.5*Mbar(2) + dCorr*Mbar(1)
    LineNormVec(2,2)=-0.5*Mbar(1) + dCorr*Mbar(2)
  END IF
!  print*,'Lxi-new' ,LineNormVec(:,1)
!  print*,'Leta-new',LineNormVec(:,2)
  doPro=DOT_PRODUCT(LineNormVec(:,1),LineNormVec(:,2))
!  print*,'doprofixed',doPro, ACOS(doPro)*180/pi
!  read*
ELSE
  ! do not change the line vectors
  LineNormVec(:,1)=Lxi
  LineNormVec(:,2)=Leta
END IF

IF(DOT_PRODUCT(LineNormVec(:,1),LineNormVecOld(:,1)).LT.0.)THEN
  LineNormVec(:,1)=-LineNormVec(:,1)
!   IPWRITE(UNIT_stdout,'(I0,A)')            ' LineNormVec-switched: Xi '
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVecOld-Xi        ', LineNormVecOld(:,1)
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVec-Xi           ', LineNormVec   (:,1)
!   IPWRITE(UNIT_stdout,'(I0,A,1(E24.12))')  ' DotProduct               ', DOT_PRODUCT(LineNormVec(:,1),LineNormVecOld(:,1))
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVecOld-Eta       ', LineNormVecOld(:,2)
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVec-Eta          ', LineNormVec   (:,2)
END IF
!
IF(DOT_PRODUCT(LineNormVec(:,2),LineNormVecOld(:,2)).LT.0.)THEN
  LineNormVec(:,2)=-LineNormVec(:,2)
!   IPWRITE(UNIT_stdout,'(I0,A)')            ' LineNormVec-switched: Eta '
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVecOld-Eta       ', LineNormVecOld(:,2)
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVec-Eta          ', LineNormVec   (:,2)
!   IPWRITE(UNIT_stdout,'(I0,A,1(E24.12))')  ' DotProduct               ', DOT_PRODUCT(LineNormVec(:,2),LineNormVecOld(:,2))
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVecOld-Xi        ', LineNormVecOld(:,1)
!   IPWRITE(UNIT_stdout,'(I0,A,2(E24.12))')  ' LineNormVec-Xi           ', LineNormVec   (:,1)
!   ! stop
END IF
! DEBUG: fix from (could become zero)
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},
!      YEAR = {2005},
END SUBROUTINE calcLineNormVec2



SUBROUTINE calcLineNormVec(BezierControlPoints2D,LineNormVec,a,b)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
!================================================================================================================================
USE MOD_Globals
USE MOD_Mesh_Vars,               ONLY:NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: a,b
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: LineNormVec(1:2)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length
!================================================================================================================================
LineNormVec=(BezierControlPoints2D(:,a,b)-BezierControlPoints2D(:,0,0))+&
            (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,b,a))
Length=SQRT(DOT_PRODUCT(LineNormVec,LineNormVec))
! DEBUG: fix from (could become zero)
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},
!      YEAR = {2005},
IF(Length.EQ.0)THEN
  ! DEBUG: is the complete IF statement dispensable?
  CALL abort(&
  __STAMP__&
  ,'Bezier Clipping -> LineNormVec is Null vector!')
  RETURN
END IF
LineNormVec=LineNormVec/Length
END SUBROUTINE calcLineNormVec


SUBROUTINE CalcSminSmax(minmax,Smin,Smax,iter)
!================================================================================================================================
! find upper and lower intersection with convex hull (or no intersection)
! find the largest and smallest roots of the convex hull, pre-sorted values minmax(:,:) are required
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo,Xi_NGeo,DeltaXi_NGeo
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
INTEGER                              :: l
!================================================================================================================================

Smin=1.5
Smax=-1.5
DO l=0,NGeo-1
  ! 1.) check traverse line UPPER/LOWER
  ! upper line, max values
  IF(minmax(2,l)*minmax(2,l+1).LE.0.)THEN
    m    = (minmax(2,l+1)-minmax(2,l))/DeltaXi_NGeo
    tmp  = Xi_NGeo(l)-minmax(2,l)/m
    Smin = MIN(tmp,Smin)
  END IF
  ! lower line, min values
  IF(minmax(1,l)*minmax(1,l+1).LE.0.)THEN
    m    = (minmax(1,l+1)-minmax(1,l))/DeltaXi_NGeo
    tmp  = Xi_NGeo(l)-minmax(1,l)/m
    Smax = MAX(tmp,Smax)
  END IF
END DO

! 2.) check BEGINNING/END upper convex hull
DO l=1,NGeo
  IF(minmax(2,0)*minmax(2,l) .LE.0.)THEN
    ! interval is the whole parameter space
    m    = (minmax(2,l)-minmax(2,0))/(DeltaXi_NGeo*l)
    tmp  = -1.0-minmax(2,0)/m
    Smin = MIN(tmp,Smin)
  END IF
  IF(minmax(1,0)*minmax(1,l) .LE.0.)THEN
    ! interval is the whole parameter space
    m    = (minmax(1,l)-minmax(1,0))/(DeltaXi_NGeo*l)
    tmp  = -1.0-minmax(1,0)/m
    Smax = MAX(tmp,Smax)
  END IF
END DO ! l
DO l=0,NGeo-1
  IF(minmax(2,l)*minmax(2,NGeo) .LE.0.)THEN
    ! interval is the whole parameter space
    m    = (minmax(2,NGeo)-minmax(2,l))/(DeltaXi_NGeo*(NGeo-l))
    tmp  = Xi_NGeo(l)-minmax(2,l)/m
    Smin = MIN(tmp,Smin)
  END IF
  IF(minmax(1,l)*minmax(1,NGeo) .LE.0.)THEN
    ! interval is the whole parameter space
    m    = (minmax(1,NGeo)-minmax(1,l))/(DeltaXi_NGeo*(NGeo-l))
    tmp  = Xi_NGeo(l)-minmax(1,l)/m
    Smax = MAX(tmp,Smax)
  END IF
END DO ! l

! 3.) check vertical line LEFT/RIGHT of convex hull
IF(minmax(1,0)*minmax(2,0)    .LE.0.)THEN
  tmp = -1.0
  Smin=MIN(tmp,Smin)
END IF
IF(minmax(1,NGeo)*minmax(2,NGeo)    .LE.0.)THEN
  tmp =  1.0
  Smax=MAX(tmp,Smax)
END IF

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
  print*,'initia shrink', smin,smax
  IF(Smin.EQ.1.5) SMin=-1. !BezierClipHit ! BezierClipHit=1+BezierClipTolerance
  IF(Smax.EQ.-1.5)SMax=1.  !BezierClipHit
END IF

END SUBROUTINE calcSminSmax


FUNCTION InsideBoundingBox(ParticlePosition,SideID)
!================================================================================================================================
! check is the particles is inside the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
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
REAL                                 :: x,y,z,P(3)
!================================================================================================================================
P=ParticlePosition-BezierControlPoints3D(1:3,0,0,SideID)
! y is perpendicular to xi & eta directions --> check first, smallest intervall
y=DOT_PRODUCT(P,SideSlabNormals(:,2,SideID))
!IF((y.LT.SideSlabIntervals(3,SideID)-epsilontol).OR.(y.GT.SideSlabIntervals(4,SideID)+epsilontol))THEN
IF((y.LT.SideSlabIntervals(3,SideID)-100.*epsMach).OR.(y.GT.SideSlabIntervals(4,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
! than xi
x=DOT_PRODUCT(P,SideSlabNormals(:,1,SideID))
!IF((x.LT.SideSlabIntervals(1,SideID)-epsilontol).OR.(x.GT.SideSlabIntervals(2,SideID)+epsilontol))THEN
IF((x.LT.SideSlabIntervals(1,SideID)-100.*epsMach).OR.(x.GT.SideSlabIntervals(2,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
! than eta
z=DOT_PRODUCT(P,SideSlabNormals(:,3,SideID))
!IF((z.LT.SideSlabIntervals(5,SideID)-epsilontol).OR.(z.GT.SideSlabIntervals(6,SideID)+epsilontol))THEN
IF((z.LT.SideSlabIntervals(5,SideID)-100.*epsMach).OR.(z.GT.SideSlabIntervals(6,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
InsideBoundingBox=.TRUE.
END FUNCTION InsideBoundingBox


FUNCTION BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Globals,                  ONLY:abort
USE MOD_Particle_Vars,            ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
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
INTEGER                              :: i
!================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
DO i=1,3!x,y,z direction
  !dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(i,:,SideID))
  dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))
  !IF(ABS(dnk).LT.epsilontol)THEN
  IF(ABS(dnk).LT.100.*epsMach)THEN
    dnk=100.*epsMach ! ÜBERPRÜFEN OB SIGN sinn macht
  END IF
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
  END IF
END DO!i
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAXVAL(alpha(1,:)) ! taken the maxvalue of the minima
minvalue=MINVAL(alpha(2,:)) ! taken the minvalue of the maxima

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
!  IF((maxvalue.LT.0).AND.(minvalue.GT.0))THEN
!    CALL abort(&
!  __STAMP__&
!  ,' BoundingBox check failed!')
!  ELSE
    IF((maxvalue.LT.lengthPartTrajectory+100*epsMach).AND.(maxvalue+100*epsMach.GT.0.))THEN
      !the first intersection is less than lengthPartTrajectory and greater 0
      BoundingBoxIntersection=.TRUE.
    ELSE
      BoundingBoxIntersection=.FALSE.
    END IF
  !END IF
ELSE
  BoundingBoxIntersection=.FALSE.
END IF
END FUNCTION BoundingBoxIntersection


FUNCTION FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,PartID,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Particle_Vars,            ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
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
INTEGER                              :: i
!================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
i=1
dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))

!IF(ABS(dnk).LT.epsilontol)THEN
IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0. ! ÜBERPRÜFEN OB SIGN sinn macht
  alpha(1,1) = -HUGE(1.0)
  alpha(2,1) =  HUGE(1.0)
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
  END IF
END IF
i=3
dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))
!IF(ABS(dnk).LT.epsilontol)THEN
IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0.
  alpha(1,3) = -HUGE(1.0)
  alpha(2,3) =  HUGE(1.0)
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(1:3,PartID),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
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


FUNCTION ComputeSurfaceDistance2(SideNormVec,BiLinearCoeff,xi,eta,PartTrajectory,PartID)
!================================================================================================================================
! compute the required vector length to intersection
! ramsey paper algorithm 3.4
!================================================================================================================================
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol
USE MOD_Particle_Vars,            ONLY:LastPartPos
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars,      ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,DIMENSION(3),INTENT(IN)         :: SideNormVec !non-oriented, averaged normal vector based on all four edges
REAL,DIMENSION(3),INTENT(IN)         :: BiLinearCoeff(1:3,4)
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: PartID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeSurfaceDistance2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
!================================================================================================================================

#ifdef CODE_ANALYZE
  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -LastPartPos(1,PartID)
  IF (PartTrajectory(1).EQ.0.) THEN
    t=SIGN(HUGE(t),t)
  ELSE
    t = t/ PartTrajectory(1)-epsilontol
  END IF
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(2(A,E15.8))') '     -- t1: ',t
   END IF
  END IF
  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -LastPartPos(2,PartID)
  IF (PartTrajectory(2).EQ.0.) THEN
    t=SIGN(HUGE(t),t)
  ELSE
    t = t/ PartTrajectory(2)-epsilontol
  END IF
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(2(A,E15.8))') '     -- t2: ',t
    END IF
  END IF
  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -LastPartPos(3,PartID)
  IF (PartTrajectory(3).EQ.0.) THEN
    t=SIGN(HUGE(t),t)
  ELSE
    t = t/ PartTrajectory(3)-epsilontol
  END IF
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(2(A,E15.8))') '     -- t3: ',t
    END IF
  END IF
#endif /*CODE_ANALYZE*/
!in ramsey paper the direction was chosen based on the largest component of PartTrajectory for preventing a division by zero.
!however, by this significant floating point inaccuracies can occur if this direction is approx. orthogonal to side normal vec.
!solution: chose direction based on SideNormVec and additionally check that no division by zero occurs.
!IF((ABS(SideNormVec(1)).GE.ABS(SideNormVec(2))).AND.(ABS(SideNormVec(1)).GE.ABS(SideNormVec(3))) &
!  .AND. .NOT.ALMOSTZERO(PartTrajectory(1)))THEN
IF((ABS(SideNormVec(1)).GE.ABS(SideNormVec(2))) .AND.(ABS(SideNormVec(1)).GE.ABS(SideNormVec(3))) &
  .AND. .NOT.ALMOSTZERO(PartTrajectory(1)))THEN
  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -LastPartPos(1,PartID)
  t = t/ PartTrajectory(1)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t1: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
ELSE IF(ABS(SideNormVec(2)).GE.ABS(SideNormVec(3)) &
  .AND. .NOT.ALMOSTZERO(PartTrajectory(2)))THEN
  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -LastPartPos(2,PartID)
  t = t/ PartTrajectory(2)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t2: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
ELSE IF(.NOT.ALMOSTZERO(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -LastPartPos(3,PartID)
  t = t/ PartTrajectory(3)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t3: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
!if PartTrajectory should be zero in largest component of SideNormVec, decide based on original check:
ELSE IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(3))))THEN
  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -LastPartPos(1,PartID)
  t = t/ PartTrajectory(1)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t1: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -LastPartPos(2,PartID)
  t = t/ PartTrajectory(2)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t2: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
ELSE
  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -LastPartPos(3,PartID)
  t = t/ PartTrajectory(3)-epsilontol
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(PartID.EQ.PARTOUT)THEN
            WRITE(UNIT_stdout,'(2(A,E15.8))') '     >> t3: ',t
          END IF
        END IF
#endif /*CODE_ANALYZE*/
END IF

!IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -LastPartPos(1,PartID)
!  t = t/ PartTrajectory(1)!-epsilontol
!ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -LastPartPos(2,PartID)
!  t = t/ PartTrajectory(2)!-epsilontol
!ELSE
!  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -LastPartPos(3,PartID)
!  t = t/ PartTrajectory(3)!-epsilontol
!END IF

ComputeSurfaceDistance2=t

END FUNCTION ComputeSurfaceDistance2


#ifdef CODE_ANALYZE
FUNCTION ComputeXi(A1,A2,eta)
#else
PURE FUNCTION ComputeXi(A1,A2,eta)
#endif
!================================================================================================================================
! compute the xi value with algorithm 3.3 of Ramsey paper
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
#ifdef CODE_ANALYZE
USE MOD_Globals, ONLY: abort,MyRank
#endif /*CODE_ANALYZE*/
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(4),INTENT(IN)         :: A1,A2
REAL,INTENT(IN)                      :: eta
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeXi
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: a,b
!================================================================================================================================

a=eta*A2(1)+A2(2)
b=eta*(A2(1)-A1(1))+A2(2)-A1(2)

IF(ABS(B).GE.ABS(A))THEN
  IF(ABS(B).LE.0.)THEN
#ifdef CODE_ANALYZE
  !ComputeXi=10
  !RETURN
    IPWRITE(*,*) 'eta',eta
    IPWRITE(*,*) 'A',a2(1)*a1(3)-a1(1)*a2(3)
    IPWRITE(*,*) 'B',a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
    IPWRITE(*,*) 'C', a1(4)*a2(2)-a1(2)*a2(4)
    IPWRITE(*,*) 'a1', a1(:)
    IPWRITE(*,*) 'a2', a2(:)
    IPWRITE(*,*) A2(1),A2(2)
    IPWRITE(*,*) A1(1),A1(2)
    CALL abort(&
    __STAMP__&
    ,' Division by zero. Invalid b')
#else
  ComputeXi=10
  RETURN
#endif /*CODE_ANALYZE*/
  END IF
  ComputeXi=(-eta*(A2(3)-A1(3))-(A2(4)-A1(4)))/b
ELSE
#ifdef CODE_ANALYZE
  IF(ABS(A).LE.0.)THEN
    CALL abort(&
    __STAMP__&
    ,' Division by zero. Invalid a')
  END IF
#endif /*CODE_ANALYZE*/
  ComputeXi=(-eta*A2(3)-A2(4))/a
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
USE MOD_Globals                ,ONLY: myrank
#endif /*CODE_ANALYZE*/
USE MOD_Globals                ,ONLY: MyRank
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
      !ReducedBezierControlPoints_temp(:,p,q)=DOT_PRODUCT((/1-.Smean(1), Smean(1)/),MATMUL(
                                       ![ReducedBezierControlPoints(p,q  ,:),ReducedBezierControlPoints(p  ,q+1,:);
                                       ! ReducedBezierControlPoints(p,q+1,:),ReducedBezierControlPoints(p+1,q+1,:)]
                                                        !,(/1-.Smean(2), Smean(2)/))
      !ReducedBezierControlPoints(:,p,q)=MinusXi*ReducedBezierControlPoints(:,p,q  )          &
      !                                 +    Xi *ReducedBezierControlPoints(:,p,q+1)*MinusEta &
      !                                 +MinusXi*ReducedBezierControlPoints(:,p  ,q+1)        &
      !                                 +    Xi *ReducedBezierControlPoints(:,p+1,q+1)*Eta
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
    IPWRITE(UNIT_stdout,'(I0,A,E18.12,x,E18.12)') ' critical error! ',alpha,alphaNorm
    IPWRITE(UNIT_stdout,'(I0,A,E18.12,x,E18.12)') ' critical error! ',alpha,alphaNorm
    IPWRITE(UNIT_stdout,'(I0,A)') ' locAlpha, locXi,locEta ' !/ lengthPartTrajectory '
    DO iInter=1,nInterSections
      WRITE(UNIT_stdout,'(I0,3(X,E18.12))') iInter,locAlpha(iInter),locXi(iInter),locEta(iInter)
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
! REAL,INTENT(INOUT),DIMENSION(:)      :: locAlpha
! INTEGER,INTENT(INOUT),DIMENSION(:)   :: locID
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
!REAL                                 :: ZeroDistance,BezierClipTolerance2
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
!              BezierControlPoints2D_temp(:,p,q)=&
!              BezierControlPoints2D_temp(:,p,q)+&
!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!              BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
!                                               *arrayNchooseK(p,l) &
!                                               *PlusXi**l        &
!                                               *MinusXi**(p-l)
!              DEBUG: optimize this !
!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!                                               *(PlusXi**l)*(MinusXi**(p-l))
!               not Horner!
!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!                                               *xiup(l)*XiDown(p-l)

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
        !BezierControlPoints2D_temp2(:,NGeo-p,q)=&
        !BezierControlPoints2D_temp2(:,NGeo-p,q)+&
        !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
        !BezierControlPoints2D_temp  (:,NGeo-l,q)*(1./(2.**p))                        &
        !                                      *arrayNchooseK(p,l)                  &
        !                                      *PlusXi**l &
        !                                      *MinusXi**(p-l)
        !DEBUG: optimize this !
        !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
        !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
        !                                 *(PlusXi**l)*(MinusXi**(p-l))
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
        !BezierControlPoints2D_temp(:,p,q)=&
        !BezierControlPoints2D_temp(:,p,q)+&
        !!BezierControlPoints2D(:,l,q)*B(p,l,Smax)
        !BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
        !                                 *arrayNchooseK(p,l) &
        !                                 *(1+XiSplit)**l        &
        !                                 *(1-XiSplit)**(p-l)
        !DEBUG: optimize this !
        !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
        !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
        !                                 *(PlusXi**l)*(MinusXi**(p-l))

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
    !      BezierControlPoints2D_temp2(:,NGeo-p,q)=&
    !      BezierControlPoints2D_temp2(:,NGeo-p,q)+&
    !      !BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
    !      BezierControlPoints2D_temp (:,NGeo-l,q)*(1./(2.**p))                        &
    !                                            *arrayNchooseK(p,l)                  &
    !                                            *(1.+2*((XiMin+1.)/(XiSplit+1.)))**(l-1) &
    !                                            *(1.-2*((XiMin+1.)/(XiSplit+1.)))**(p-l)
          !DEBUG: optimize this !
          !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
          !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
          !                                 *(PlusXi**l)*(MinusXi**(p-l))
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
  !iClipIter   =tmpnClip
  !nXiClip     =tmpnXi
  !nEtaClip    =tmpnEta
  !LineNormVec =tmpLineNormVec
  !ClipMode    =tmpClipMode
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
!              BezierControlPoints2D_temp(:,p,q)=&
!              BezierControlPoints2D_temp(:,p,q)+&
!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!              BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
!                                               *arrayNchooseK(p,l) &
!                                               *(1+XiMax)**l        &
!                                               *(1-XiMax)**(p-l)
!             !DEBUG: optimize this !
!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!                                               *(PlusXi**l)*(MinusXi**(p-l))
!
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
          !BezierControlPoints2D_temp(:,NGeo-p,q)=&
          !BezierControlPoints2D_temp(:,NGeo-p,q)+&
          !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
          !BezierControlPoints2D     (:,NGeo-l,q)*(1./(2.**p))                        &
          !                                      *arrayNchooseK(p,l)                  &
          !                                      *(1.+2*((XiMin+1.)/(XiMax+1.)))**(l-1) &
          !                                      *(1.-2*((XiMin+1.)/(XiMax+1.)))**(p-l)
          !DEBUG: optimize this !
          !BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)             &
          !                                +BezierControlPoints2D  (:,NGeo-l,q)*FacNchooseK(p,l) &
          !                                *(PlusXi**l)*(MinusXi**(p-l))
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
USE MOD_Particle_Surfaces_Vars,  ONLY:EtaArray!,locAlpha
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
! REAL,INTENT(INOUT),DIMENSION(:)      :: locAlpha
! INTEGER,INTENT(INOUT),DIMENSION(:)   :: locID
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
!REAL                                 :: ZeroDistance,BezierClipTolerance2
REAL                                 :: PlusEta,MinusEta!,MinusXi,PlusXi
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
!              BezierControlPoints2D_temp(:,q,p)=&
!              BezierControlPoints2D_temp(:,q,p)+&
!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!              BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!                                               *arrayNchooseK(p,l) &
!                                               *(1.+Etamax)**l       &
!                                               *(1.-Etamax)**(p-l)
!              DEBUG: optimize this !
!              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                               *(PlusEta**l)*(MinusEta**(p-l))
!
!
!
!

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
        !BezierControlPoints2D_temp2(:,q,p)=&
        !BezierControlPoints2D_temp2(:,q,p)+&
        !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
        !BezierControlPoints2D_temp(:,q,NGeo-l)*(1./(2.**p))                     &
        !                                      *arrayNchooseK(p,l)               &
        !                                      *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
        !                                      *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
        !  !DEBUG: optimize this !
        !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
        !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
        !                                 *(PlusEta**l)*(MinusEta**(p-l))
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
        !BezierControlPoints2D_temp(:,q,p)=&
        !BezierControlPoints2D_temp(:,q,p)+&
        !!BezierControlPoints2D(:,l,q)*B(p,l,Smax)
        !BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
        !                                 *arrayNchooseK(p,l) &
        !                                 *(1.+EtaSplit)**l       &
        !                                 *(1.-EtaSplit)**(p-l)
        !DEBUG: optimize this !
!            BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                             +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                             *(PlusEta**l)*(MinusEta**(p-l))
!
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
        !  BezierControlPoints2D_temp2(:,q,p)=&
        !  BezierControlPoints2D_temp2(:,q,p)+&
        !  !BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
        !  BezierControlPoints2D_temp(:,q,NGeo-l)*(1./(2.**p))                     &
        !                                        *arrayNchooseK(p,l)               &
        !                                        *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
        !                                        *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
        !DEBUG: optimize this !
          !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
          !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
          !                                 *(PlusEta**l)*(MinusEta**(p-l))
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
  !iClipIter   =tmpnClip
  !! restore values to allow for correct back-trafo of intervals (required for intersectionpoint)
  !iClipIter   =tmpnClip
  !nXiClip     =tmpnXi
  !nEtaClip    =tmpnEta
  !LineNormVec =tmpLineNormVec
  !ClipMode    =tmpClipMode
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
!              BezierControlPoints2D_temp(:,q,p)=&
!              BezierControlPoints2D_temp(:,q,p)+&
!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!              BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!                                               *arrayNchooseK(p,l) &
!                                               *(1.+Etamax)**l       &
!                                               *(1.-Etamax)**(p-l)
!             !DEBUG: optimize this !
!              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                               *(PlusEta**l)*(MinusEta**(p-l))
!
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
          !BezierControlPoints2D_temp(:,q,p)=&
          !BezierControlPoints2D_temp(:,q,p)+&
          !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
          !BezierControlPoints2D     (:,q,NGeo-l)*(1./(2.**p))                     &
          !                                      *arrayNchooseK(p,l)               &
          !                                      *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
          !                                      *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
          !DEBUG: optimize this !
          !BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)        &
          !                                 +BezierControlPoints2D(:,q,NGeo-l)*FacNchooseK(p,l) &
          !                                 *(PlusEta**l)*(MinusEta**(p-l))
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
REAL                                 :: Length,doPro!,dCorr
REAL,DIMENSION(2)                    :: LXi, Leta,Mbar,Mbar2
!================================================================================================================================

! compute Lxi vector
! 1) the initial Lxi vector is the combination of the two bounding vectors  which point in eta direction
!     to get the 1D distances of each point via scalar product, we have to right-rotate the vector
LXi=(BezierControlPoints2D(:,   0,NGeo)-BezierControlPoints2D(:,   0,   0))+&
    (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,NGeo,   0))

Mbar =(BezierControlPoints2D(:,   0,NGeo)-BezierControlPoints2D(:,   0,   0))
MBar2=(BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,NGeo,   0))
!print*,'DoProXi',DOT_PRODUCT(MBar,MBar2)


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
!print*,'leta',Leta
doPro=DOT_PRODUCT(Lxi,Leta)
!IF(doPro.LT.0) Lxi=-Lxi
!print*,'doPro',doPro

! HERE some fixes for the line vectors
!IF(ABS(doPro).GT.0.5)THEN
!  MBar=0.5*(Lxi+Leta)
!  dCorr=0.5*SQRT(3.)
!  IF(doPro.GT.0.)THEN
!    ! Lxi
!    Lxi(1)= 0.5*Mbar(2) + dCorr*Mbar(1)
!    Lxi(2)=-0.5*Mbar(1) + dCorr*Mbar(2)
!    ! Leta
!    Leta(1)=-0.5*Mbar(2) + dCorr*Mbar(1)
!    Leta(2)= 0.5*Mbar(1) + dCorr*Mbar(2)
!  ELSE
!    ! Lxi
!    Lxi(1)=-0.5*Mbar(2) + dCorr*Mbar(1)
!    Lxi(2)= 0.5*Mbar(1) + dCorr*Mbar(2)
!    ! Leta
!    Leta(1)= 0.5*Mbar(2) + dCorr*Mbar(1)
!    Leta(2)=-0.5*Mbar(1) + dCorr*Mbar(2)
!  END IF
!END IF

! 3) rotate  both vectors by -90 degree
LineNormVec(1,1) =-LXi (2)
LineNormVec(2,1) = LXi (1)
LineNormVec(1,2) =-Leta(2) ! stephen meint minus1
LineNormVec(2,2) = Leta(1)
!print*,'linenormvec1',LineNormVec(:,1)
!print*,'linenormvec2',LineNormVec(:,2)


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


SUBROUTINE ComputeAuxBCIntersection     (isHit                       &
                                        ,PartTrajectory              &
                                        ,lengthPartTrajectory        &
                                        ,AuxBCIdx                    &
                                        ,alpha                       &
                                        ,iPart                       &
                                        ,opt_CriticalParllelInSide   )
!===================================================================================================================================
! Compute the Intersection with auxBC. (based partly on PlanarRect)
! Implemtented types:
! - plane
! - cylinder
! - cone
! - parabol(oid)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Utils                  ,ONLY: QuadraticSolver
!USE MOD_Globals_Vars           ,ONLY:  epsMach
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Particle_Surfaces_Vars ,ONLY: epsilontol
USE MOD_Particle_Boundary_Vars ,ONLY: AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: AuxBCIdx,iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha
LOGICAL,INTENT(OUT)               :: isHit
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: r_vec(3),n_vec(3),locSideDistance,coeffA,alphaNorm,radius,lmin,lmax,halfangle
REAL                              :: axis(3),tang1(3),tang2(3),geomatrix(3,3),matU(3,1),matLambda(3,1),A(1,1),B(1,1),C(1,1),cos2inv
REAL                              :: geomatrix4(4,4),rotmatrix(3,3),matU4(4,1),matLambda4(4,1),zfac
REAL                              :: trajTang(2),originTang(2),roots(2),intersec(3),alphadir(2),origindist(2) !,roots2(2)
INTEGER                           :: nRoot !,nRoot2
LOGICAL                           :: CriticalParallelInSide,inwards
!===================================================================================================================================
isHit=.FALSE.
SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
CASE ('plane')
  r_vec=AuxBC_plane(AuxBCMap(AuxBCIdx))%r_vec
  n_vec=AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
  radius=AuxBC_plane(AuxBCMap(AuxBCIdx))%radius
  coeffA=DOT_PRODUCT(n_vec,PartTrajectory)
  CriticalParallelInSide=.FALSE.
  IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.
  locSideDistance = DOT_PRODUCT(n_vec,r_vec) - DOT_PRODUCT(LastPartPos(1:3,iPart),n_vec)
  IF(CriticalParallelInSide)THEN ! particle parallel to side
    IF(ALMOSTZERO(locSideDistance))THEN ! particle on/in side
      IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
      ! move particle eps into interior (?!)
      alpha=-1.
      RETURN
    END IF
    IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
    alpha=-1.
    RETURN
  ELSE
    IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
    alpha=locSideDistance/coeffA
  END IF
  alphaNorm=alpha/lengthPartTrajectory
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory - r_vec !intersec from basepoint, not origin!
  ! check besides alpha and radius already here the dir. of trajectory since no inner auxBCs possible (can happen due to tolerances)
  IF((alphaNorm.GT.1.0) .OR.(alphaNorm.LT.-epsilontol) .OR. SQRT(DOT_PRODUCT(intersec,intersec)).GT.radius &
    .OR. DOT_PRODUCT(n_vec,PartTrajectory).LT.0.)THEN
    ishit=.FALSE.
    alpha=-1.0
    RETURN
  END IF
!  epsLoc=1.0+100.*epsMach
!  xi=...
!  IF(ABS(xi).GT.epsLoc)THEN
!    alpha=-1.0
!    RETURN
!  END IF
!  IF(ABS(eta).GT.epsLoc)THEN
!    alpha=-1.0
!    RETURN
!  END IF
  isHit=.TRUE.
CASE ('cylinder','cone','parabol')
  IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE. !not used for cylinder and cone
  IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cylinder') THEN
    r_vec=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
    axis=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
    lmin=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%lmin
    lmax=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%lmax
    IF (axis(3).NE.0.) THEN
      tang1(1) = 1.0
      tang1(2) = 1.0
      tang1(3) = -(axis(1)+axis(2))/axis(3)
    ELSE
      IF (axis(2).NE.0.) THEN
        tang1(1) = 1.0
        tang1(3) = 1.0
        tang1(2) = -(axis(1)+axis(3))/axis(2)
      ELSE
        IF (axis(1).NE.0.) THEN
          tang1(2) = 1.0
          tang1(3) = 1.0
          tang1(1) = -(axis(2)+axis(3))/axis(1)
        ELSE
          CALL abort(&
__STAMP__&
,'Error in ComputeAuxBCIntersection, axis is zero for AuxBC',AuxBCIdx)
        END IF
      END IF
    END IF
    tang1=UNITVECTOR(tang1)
    tang2=CROSSNORM(axis,tang1)
    radius=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%radius
    inwards=AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards
    !- project trajectory and origin into circle-area of cylinder
    trajTang(1)=DOT_PRODUCT(tang1,PartTrajectory)
    trajTang(2)=DOT_PRODUCT(tang2,PartTrajectory)
    originTang(1)=DOT_PRODUCT(tang1,LastPartPos(1:3,iPart)-r_vec)
    originTang(2)=DOT_PRODUCT(tang2,LastPartPos(1:3,iPart)-r_vec)
    !- solve quadratic equation from trajectory inserted in circle-equation
    CALL QuadraticSolver(trajTang(1)*trajTang(1)+trajTang(2)*trajTang(2) &
      ,2.*originTang(1)*trajTang(1)+2.*originTang(2)*trajTang(2) &
      ,originTang(1)*originTang(1)+originTang(2)*originTang(2)-radius*radius &
      ,nRoot,roots(1),roots(2))
  ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cone') THEN
    r_vec=AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
    axis=AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
    lmin=AuxBC_cone(AuxBCMap(AuxBCIdx))%lmin
    lmax=AuxBC_cone(AuxBCMap(AuxBCIdx))%lmax
    halfangle=AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle
    cos2inv=1./COS(halfangle)**2
    inwards=AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards
    !- coefficients and matrices according to "Intersection of a Line and a Cone" by David Eberly 2000/2014, Geometric Tools, CC
    geomatrix=AuxBC_cone(AuxBCMap(AuxBCIdx))%geomatrix
    !geomatrix2=AuxBC_cone(AuxBCMap(AuxBCIdx))%geomatrix2
    !rotmatrix=AuxBC_cone(AuxBCMap(AuxBCIdx))%rotmatrix
    matU(:,1)=PartTrajectory
    matLambda(:,1)=LastPartPos(1:3,iPart)-r_vec
    A=MATMUL(MATMUL(TRANSPOSE(matU),geomatrix),matU)
    B=2.*MATMUL(MATMUL(TRANSPOSE(matU),geomatrix),matLambda)
    C=MATMUL(MATMUL(TRANSPOSE(matLambda),geomatrix),matLambda)
    !- solve quadratic equation from trajectory inserted in cone-equation
    CALL QuadraticSolver(A(1,1),B(1,1),C(1,1),nRoot,roots(1),roots(2))
  ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'parabol') THEN
    r_vec=AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
    axis=AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
    lmin=AuxBC_parabol(AuxBCMap(AuxBCIdx))%lmin
    lmax=AuxBC_parabol(AuxBCMap(AuxBCIdx))%lmax
    zfac=AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac
    inwards=AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards
    geomatrix4=AuxBC_parabol(AuxBCMap(AuxBCIdx))%geomatrix4
    rotmatrix=AuxBC_parabol(AuxBCMap(AuxBCIdx))%rotmatrix
    matU(:,1)=PartTrajectory
    matLambda(:,1)=LastPartPos(1:3,iPart)-r_vec
    matU=MATMUL(rotmatrix,matU)
    matLambda=MATMUL(rotmatrix,matLambda)
    matU4(1:3,1)=matU(1:3,1)
    matU4(4,1)=0.
    matLambda4(1:3,1)=matLambda(1:3,1)
    matLambda4(4,1)=1.
    A=MATMUL(MATMUL(TRANSPOSE(matU4),geomatrix4),matU4)
    B=2.*MATMUL(MATMUL(TRANSPOSE(matU4),geomatrix4),matLambda4)
    C=MATMUL(MATMUL(TRANSPOSE(matLambda4),geomatrix4),matLambda4)
    !- solve quadratic equation from trajectory inserted in parabol-equation
    CALL QuadraticSolver(A(1,1),B(1,1),C(1,1),nRoot,roots(1),roots(2))
  ELSE
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END IF !cylinder, cone, or paraboloid
  SELECT CASE (nRoot)
  CASE (1)
    alpha=roots(1)
    !- check for normal vec / trajectory direction
    ! (already here since no inner auxBCs possible (can happen due to tolerances)
    intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
    origindist(1) = DOT_PRODUCT(intersec-r_vec,axis)
    IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cylinder') THEN
      n_vec = intersec - ( r_vec + axis*origindist(1) )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cone') THEN
      n_vec = intersec - ( r_vec + axis*origindist(1)*cos2inv )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'parabol') THEN
      n_vec = intersec - ( r_vec + axis*(origindist(1)+0.5*zfac) )
    ELSE
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist')
    END IF
    IF (.NOT.inwards) n_vec=-n_vec
    IF(DOT_PRODUCT(n_vec,PartTrajectory).LT.0.)THEN
      ishit=.FALSE.
      alpha=-1.0
      RETURN
    END IF
    !- check for lmin and lmax
    IF (origindist(1).LT.lmin .OR. origindist(1).GT.lmax) THEN
      ishit=.FALSE.
      alpha=-1.0
      RETURN
    END IF
  CASE (2)
    !- 2 roots: check for smallest alpha>-eps
    IF (roots(1).LT.roots(2)) THEN
      IF (roots(1).GE.-epsilontol*lengthPartTrajectory) THEN
        alpha=roots(1)
      ELSE
        alpha=roots(2)
        roots(2)=roots(1)
        roots(1)=alpha
      END IF
    ELSE
      IF (roots(2).GE.-epsilontol*lengthPartTrajectory) THEN
        alpha=roots(2)
        roots(2)=roots(1)
        roots(1)=alpha
      ELSE
        alpha=roots(1)
      END IF
    END IF
    !- check for lmin and lmax of cylinder and normal vec / trajectory direction
    ! (already here since no inner auxBCs possible (can happen due to tolerances)
    intersec = LastPartPos(1:3,iPart) + roots(1)*PartTrajectory
    origindist(1) = DOT_PRODUCT(intersec-r_vec,axis)
    IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cylinder') THEN
      n_vec = intersec - ( r_vec + axis*origindist(1) )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cone') THEN
      n_vec = intersec - ( r_vec + axis*origindist(1)*cos2inv )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'parabol') THEN
      n_vec = intersec - ( r_vec + axis*(origindist(1)+0.5*zfac) )
    ELSE
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist')
    END IF
    alphadir(1)=DOT_PRODUCT(n_vec,PartTrajectory)
    intersec = LastPartPos(1:3,iPart) + roots(2)*PartTrajectory
    origindist(2) = DOT_PRODUCT(intersec-r_vec,axis)
    IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cylinder') THEN
      n_vec = intersec - ( r_vec + axis*origindist(2) )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'cone') THEN
      n_vec = intersec - ( r_vec + axis*origindist(2)*cos2inv )
    ELSE IF (TRIM(AuxBCType(AuxBCIdx)).EQ.'parabol') THEN
      n_vec = intersec - ( r_vec + axis*(origindist(2)+0.5*zfac) )
    ELSE
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist')
    END IF
    alphadir(2)=DOT_PRODUCT(n_vec,PartTrajectory)
    IF (.NOT.inwards) alphadir=-alphadir
    IF (alphadir(1).GE.0. .AND. origindist(1).GE.lmin .AND. origindist(1).LE.lmax) THEN
      ! alpha stays alpha
    ELSE IF (alphadir(2).GE.0. .AND. origindist(2).GE.lmin .AND. origindist(2).LE.lmax) THEN
      alpha=roots(2)
    ELSE
      ishit=.FALSE.
      alpha=-1.0
      RETURN
    END IF
  CASE DEFAULT
    ishit=.FALSE.
    alpha=-1.0
    RETURN
  END SELECT
  alphaNorm=alpha/lengthPartTrajectory
  IF((alphaNorm.GT.1.0) .OR.(alphaNorm.LT.-epsilontol))THEN
    ishit=.FALSE.
    alpha=-1.0
    RETURN
  END IF
  isHit=.TRUE.
CASE DEFAULT
  SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(AuxBCIdx))
  CALL abort(&
    __STAMP__&
    ,'AuxBC does not exist')
END SELECT

END SUBROUTINE ComputeAuxBCIntersection


END MODULE MOD_Particle_Intersection
