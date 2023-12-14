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

MODULE MOD_Photon_TrackingTools
!===================================================================================================================================
! Routines for photon tracking in radiave transfer solver
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE PhotonThroughSideCheck3DFast
  MODULE PROCEDURE PhotonThroughSideCheck3DFast
END INTERFACE

PUBLIC :: PhotonThroughSideCheck3DFast, PhotonIntersectionWithSide, CalcAbsoprtion, PerfectPhotonReflection, DiffusePhotonReflection
PUBLIC :: CalcWallAbsoprtion, PointInObsCone, PhotonIntersectSensor, PhotonThroughSideCheck3DDir, PhotonIntersectionWithSide2DDir
PUBLIC :: PhotonIntersectionWithSide2D, RotatePhotonIn2DPlane, PerfectPhotonReflection2D,DiffusePhotonReflection2D, PhotonOnLineOfSight
PUBLIC :: PeriodicPhotonBC
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


SUBROUTINE PhotonThroughSideCheck3DFast(iLocSide,Element,ThroughSide,TriNum, IsMortar)
!===================================================================================================================================
!> Routine to check whether a particle crossed the given triangle of a side. The determinant between the normalix_photon_startd trajectory
!> vector and the vectors from two of the three nodes to the old particle position is calculated. If the determinants for the three
!> possible combinations are greater than x_photon_startro, then the particle went through this triangle of the side.
!> Note that if this is a mortar side, the side of the small neighbouring mortar element has to be checked. Thus, the orientation
!> is reversed.
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars        ,ONLY: EpsMach
USE MOD_Globals             ,ONLY: abort
USE MOD_Particle_Mesh_Vars  ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
LOGICAL,INTENT(OUT)              :: ThroughSide
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz!, Vall
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3),tolerance
!===================================================================================================================================
CNElemID = GetCNElemID(Element)

! Sanity check
IF(CNElemID.LE.0) CALL abort(__STAMP__,'PhotonThroughSideCheck3DFast() found CNElemID<=0')

ThroughSide = .FALSE.

Px = PhotonProps%PhotonLastPos(1)
Py = PhotonProps%PhotonLastPos(2)
Pz = PhotonProps%PhotonLastPos(3)

! Normalix_photon_startd particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = PhotonProps%PhotonDirection(1)
Vy = PhotonProps%PhotonDirection(2)
Vz = PhotonProps%PhotonDirection(3)
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
!--- a right-hand-szstem. If yes for all edges: vector goes potentially through side
det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
          (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
          (Ax(1) * Vy - Ay(1) * Vx) * Az(3))

det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
          (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
          (Ax(2) * Vy - Ay(2) * Vx) * Az(1))

det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
          (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
          (Ax(3) * Vy - Ay(3) * Vx) * Az(2))

! Changed tolerance from -epsMach to -0.1*epsMach because this check finds intersections, but when PhotonIntersectionWithSide() is
! called, the photon vector and the side normal vector are perfectly parallel, e.g.,  v = (/0,0,-1/) and n=(/0,-1,0./), hence, the
! scalar product is exactly zero and the routine breaks.
tolerance = -0.1*epsMach
! Comparison of the determinants with eps, where a x_photon_startro is stored (due to machine precision)
IF(det(1).LT.tolerance) RETURN
IF(det(2).LT.tolerance) RETURN
IF(det(3).LT.tolerance) RETURN
ThroughSide = .TRUE.

RETURN
END SUBROUTINE PhotonThroughSideCheck3DFast


SUBROUTINE PhotonThroughSideCheck3DDir(iLocSide,CNElemID,ThroughSide,TriNum,StartPoint,Dir)
!===================================================================================================================================
!> Routine to check whether a particle crossed the given triangle of a side. The determinant between the normalix_photon_startd trajectory
!> vector and the vectors from two of the three nodes to the old particle position is calculated. If the determinants for the three
!> possible combinations are greater than x_photon_startro, then the particle went through this triangle of the side.
!> Note that if this is a mortar side, the side of the small neighbouring mortar element has to be checked. Thus, the orientation
!> is reversed.
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars              ,ONLY: EpsMach
USE MOD_Particle_Mesh_Vars, ONLY : NodeCoords_Shared, ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: CNElemID
INTEGER,INTENT(IN)               :: TriNum
LOGICAL,INTENT(OUT)              :: ThroughSide
REAL, INTENT(IN)                 :: StartPoint(3), Dir(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz!, Vall
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3)
!===================================================================================================================================
ThroughSide = .FALSE.

Px = StartPoint(1)
Py = StartPoint(2)
Pz = StartPoint(3)

! Normalix_photon_startd particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = Dir(1)
Vy = Dir(2)
Vz = Dir(3)
! Get the coordinates of the first node and the vector from the particle position to the node
xNode(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNode(1) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNode(1) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
Ax(1) = xNode(1) - Px
Ay(1) = yNode(1) - Py
Az(1) = zNode(1) - Pz
! Get the vectors to the other two nodes, depending on the triangle number

DO n = 2,3
  NodeID = n+TriNum-1       ! m = true node number of the sides (TriNum=1: NodeID=2,3; TriNum=2: NodeID=3,4)
  xNode(n) = NodeCoords_Shared(1,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
  yNode(n) = NodeCoords_Shared(2,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
  zNode(n) = NodeCoords_Shared(3,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)

  Ax(n) = xNode(n) - Px
  Ay(n) = yNode(n) - Py
  Az(n) = zNode(n) - Pz
END DO

!--- check whether v and the vectors from the particle to the two edge nodes build
!--- a right-hand-szstem. If yes for all edges: vector goes potentially through side
det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
          (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
          (Ax(1) * Vy - Ay(1) * Vx) * Az(3))

det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
          (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
          (Ax(2) * Vy - Ay(2) * Vx) * Az(1))

det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
          (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
          (Ax(3) * Vy - Ay(3) * Vx) * Az(2))

! Comparison of the determinants with eps, where a x_photon_startro is stored (due to machine precision)
IF ((det(1).ge.-epsMach).AND.(det(2).ge.-epsMach).AND.(det(3).ge.-epsMach)) THEN
  ThroughSide = .TRUE.
END IF

RETURN

END SUBROUTINE PhotonThroughSideCheck3DDir


!===================================================================================================================================
!> Routine to check whether a photon crossed the given side.
!===================================================================================================================================
SUBROUTINE PhotonIntersectionWithSide2D(iLocSide,Element,ThroughSide,IntersectionPos,isLastSide,Distance)
! MODULES
USE MOD_Globals             ,ONLY: abort
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemSideNodeID2D_Shared, NodeCoords_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: ThroughSide
INTEGER,INTENT(IN)               :: iLocSide, Element
REAL, INTENT(OUT)                :: IntersectionPos(3)
REAL, INTENT(OUT)                :: Distance
LOGICAL, INTENT(IN)              :: isLastSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
REAL                             :: y_photon_start,x_photon_start,yNode1,xNode1,yNode2,xNode2,sy,sz,sx
REAL                             :: l1,S1,l2,S2,l,S
REAL                             :: beta, alpha, deltay, a, b, c, tmpsqrt
!===================================================================================================================================
CNElemID = GetCNElemID(Element)

! Sanity check
IF(CNElemID.LE.0) CALL abort(__STAMP__,'PhotonIntersectionWithSide2D() found CNElemID<=0')

ThroughSide = .FALSE.

xNode1 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(1,iLocSide, CNElemID))
yNode1 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(1,iLocSide, CNElemID))
xNode2 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(2,iLocSide, CNElemID))
yNode2 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(2,iLocSide, CNElemID))

x_photon_start=PhotonProps%PhotonLastPos(1)
y_photon_start=PhotonProps%PhotonLastPos(2)

sx=PhotonProps%PhotonDirection(1)
sy=PhotonProps%PhotonDirection(2)
sz=PhotonProps%PhotonDirection(3)

IF (sx .EQ. 0.0) THEN
  IF (xNode1.EQ.xNode2) THEN
    l = (y_photon_start-yNode1)/(yNode2-yNode1)
  ELSE
    l = (x_photon_start-xNode1)/(xNode2-xNode1)
  END IF
  a = sy*sy + sz*sz
  b = 2*sy*y_photon_start
  c = y_photon_start*y_photon_start - yNode1*yNode1 + 2.*l*yNode1*yNode1 - l*l*yNode1*yNode1 &
      - 2.*yNode1*yNode2*l + 2.*yNode1*yNode2*l*l - yNode2*yNode2*l*l
  tmpsqrt = b*b - 4.*a*c
  IF (tmpsqrt.LE.0.0) THEN
    RETURN
  END IF
  S1 = (-b+SQRT(tmpsqrt))/(2.*a)
  S2 = (-b-SQRT(tmpsqrt))/(2.*a)

  IF(isLastSide) THEN
    IF (ALMOSTEQUAL(S1,S2)) THEN
      RETURN ! TODO
    ELSE IF (ABS(S1).GT.ABS(S2)) THEN
      S=S1
    ELSE
      S=S2
    END IF
  ELSE
    IF (S1.LE.0.0) THEN
      S = S2
    ELSE
      IF (S2.GT.0.0) THEN
        IF(S2.GT.S1) THEN
          S = S1
        ELSE
          S = S2
        END IF
      ELSE
        S = S1
      END IF
    END IF
  END IF


ELSE
  alpha = (xNode1 - x_photon_start) / sx
  beta = (xNode2 - xNode1) / sx
  deltay = (yNode2 - yNode1)
  a = beta*beta*sy*sy - deltay*deltay + beta*beta*sz*sz
  b = 2.*beta*sy*y_photon_start + 2.*alpha*beta*sy*sy - 2.*deltay*yNode1 + 2.*alpha*beta*sz*sz
  c = y_photon_start*y_photon_start - yNode1*yNode1 + 2.*alpha*sy*y_photon_start + alpha*alpha*sy*sy + sz*sz*alpha*alpha
  tmpsqrt = b*b - 4.*a*c
  IF (tmpsqrt.LE.0.0) THEN
    RETURN
  END IF
  l1 = (-b + SQRT(tmpsqrt))/(2.*a)
  S1 = (xNode1-x_photon_start+(xNode2-xNode1)*l1)/sx
  l2 = (-b - SQRT(tmpsqrt))/(2.*a)
  S2 = (xNode1-x_photon_start+(xNode2-xNode1)*l2)/sx

  IF (isLastSide) THEN
    IF (ALMOSTEQUAL(S1,S2).AND.ALMOSTEQUAL(ABS(l1),ABS(l2))) THEN
      RETURN
    ELSE IF (ALMOSTEQUAL(S1,S2)) THEN
      IF (ABS(l1).GT.ABS(l2)) THEN
        l=l1; S=S1
      ELSE
        l=l2; S=S2
      END IF
    ELSE IF (ALMOSTZERO(S1).AND.ALMOSTZERO(S2)) THEN
      IF (ABS(l1).GT.ABS(l2)) THEN
        l=l1; S=S1
      ELSE
        l=l2; S=S2
      END IF
    ELSE IF (ABS(S1).GT.ABS(S2)) THEN       !though same spot again, caused by numerical inaccuray (discard shorter solution)
      l=l1; S=S1
    ELSE
      l=l2; S=S2
    END IF
  ELSE IF ((l1.LE.0.0).OR.(l1.GE.1.0)) THEN !if 1 is not a valid intersection -> 2
    l = l2; S = S2
  ELSE                                      !1 is valid intersection
    IF ((S1.LE.0.0)) THEN                   !1 would be moving backwards -> 2
      l = l2; S = S2
    ELSE
      IF ((l2.GT.0.0).AND.(l2.LT.1.0).AND.(S2.GT.0.0)) THEN !1 and 2 valid -> chose shorter one
        IF (S2.GT.S1) THEN
          l=l1; S=S1
        ELSE
          l=l2; S=S2
        END IF
      ELSE                                  !1 is only valid intersection -> 1
        l=l1; S=S1
      END IF
    END IF
  END IF

END IF

IF((S .GT. 0.0) .AND. (0.0 .LE. l) .AND. (l .LE. 1.0)) THEN
  ThroughSide = .TRUE.
  IntersectionPos(1) = PhotonProps%PhotonLastPos(1) + S*sx
  IntersectionPos(2) = PhotonProps%PhotonLastPos(2) + S*sy
  IntersectionPos(3) = S*sz
  Distance = S
END IF

END SUBROUTINE PhotonIntersectionWithSide2D


!===================================================================================================================================
!> Routine to check whether a photon crossed the given side.
!===================================================================================================================================
SUBROUTINE PhotonIntersectionWithSide2DDir(iLocSide,CNElemID,ThroughSide,StartPoint, Dir, IntersectionPos, Distance)
! MODULES
USE MOD_Particle_Mesh_Vars,          ONLY : ElemSideNodeID2D_Shared, NodeCoords_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: ThroughSide
INTEGER,INTENT(IN)               :: iLocSide, CNElemID
REAL, INTENT(OUT), OPTIONAL      :: IntersectionPos(3), Distance
REAL,INTENT(IN)                  :: StartPoint(3), Dir(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: y_photon_start,x_photon_start,yNode1,xNode1,yNode2,xNode2,sy,sz,sx
REAL                             :: l1,S1,l2,S2,l,S
REAL                             :: beta, alpha, deltay, a, b, c, tmpsqrt
!===================================================================================================================================
ThroughSide = .FALSE.

xNode1 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(1,iLocSide, CNElemID))
yNode1 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(1,iLocSide, CNElemID))
xNode2 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(2,iLocSide, CNElemID))
yNode2 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(2,iLocSide, CNElemID))

x_photon_start=StartPoint(1)
y_photon_start=StartPoint(2)

sx=Dir(1)
sy=Dir(2)
sz=Dir(3)

IF (sx .EQ. 0.0) THEN
  IF (xNode1.EQ.xNode2) THEN
    l = (y_photon_start-yNode1)/(yNode2-yNode1)
  ELSE
    l = (x_photon_start-xNode1)/(xNode2-xNode1)
  END IF
  a = sy*sy + sz*sz
  b = 2*sy*y_photon_start
  c = y_photon_start*y_photon_start - yNode1*yNode1 + 2.*l*yNode1*yNode1 - l*l*yNode1*yNode1 &
      - 2.*yNode1*yNode2*l + 2.*yNode1*yNode2*l*l - yNode2*yNode2*l*l
  tmpsqrt = b*b - 4.*a*c
  IF (tmpsqrt.LE.0.0) THEN
    RETURN
  END IF
  S1 = (-b+SQRT(tmpsqrt))/(2.*a)
  S2 = (-b-SQRT(tmpsqrt))/(2.*a)

  IF (S1.LE.0.0) THEN
    S = S2
  ELSE
    IF (S2.GT.0.0) THEN
      IF(S2.GT.S1) THEN
        S = S1
      ELSE
        S = S2
      END IF
    ELSE
      S = S1
    END IF
  END IF
ELSE
  alpha = (xNode1 - x_photon_start) / sx
  beta = (xNode2 - xNode1) / sx
  deltay = (yNode2 - yNode1)
  a = beta*beta*sy*sy - deltay*deltay + beta*beta*sz*sz
  b = 2.*beta*sy*y_photon_start + 2.*alpha*beta*sy*sy - 2.*deltay*yNode1 + 2.*alpha*beta*sz*sz
  c = y_photon_start*y_photon_start - yNode1*yNode1 + 2.*alpha*sy*y_photon_start + alpha*alpha*sy*sy + sz*sz*alpha*alpha
  tmpsqrt = b*b - 4.*a*c
  IF (tmpsqrt.LE.0.0) THEN
    RETURN
  END IF
  l1 = (-b + SQRT(tmpsqrt))/(2.*a)
  S1 = (xNode1-x_photon_start+(xNode2-xNode1)*l1)/sx
  l2 = (-b - SQRT(tmpsqrt))/(2.*a)
  S2 = (xNode1-x_photon_start+(xNode2-xNode1)*l2)/sx

  IF ((l1.LE.0.0).OR.(l1.GE.1.0)) THEN !if 1 is not a valid intersection -> 2
    l = l2; S = S2
  ELSE                                      !1 is valid intersection
    IF ((S1.LE.0.0)) THEN                   !1 would be moving backwards -> 2
      l = l2; S = S2
    ELSE
      IF ((l2.GT.0.0).AND.(l2.LT.1.0).AND.(S2.GT.0.0)) THEN !1 and 2 valid -> chose shorter one
        IF (S2.GT.S1) THEN
          l=l1; S=S1
        ELSE
          l=l2; S=S2
        END IF
      ELSE                                  !1 is only valid intersection -> 1
        l=l1; S=S1
      END IF
    END IF
  END IF

END IF

IF((S .GT. 0.0) .AND. (0.0 .LE. l) .AND. (l .LE. 1.0)) THEN
  ThroughSide = .TRUE.
  IF (PRESENT(IntersectionPos)) THEN
    IntersectionPos(1) = StartPoint(1) + S*sx
    IntersectionPos(2) = StartPoint(2) + S*sy
    IntersectionPos(3) = S*sz
    Distance = S
  END IF
END IF
END SUBROUTINE PhotonIntersectionWithSide2DDir


!===================================================================================================================================
!> Routine to check whether a photon crossed the given side.
!===================================================================================================================================
SUBROUTINE RotatePhotonIn2DPlane(IntersectionPos)
! MODULES
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(OUT)                :: IntersectionPos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: NewYPho, NewYVelo
!===================================================================================================================================
PhotonProps%PhotonLastPos(1:3) = IntersectionPos(1:3)
NewYPho = SQRT(PhotonProps%PhotonLastPos(2)**2 + PhotonProps%PhotonLastPos(3)**2)
! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
!           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
NewYVelo = (PhotonProps%PhotonDirection(2)*PhotonProps%PhotonLastPos(2) &
          + PhotonProps%PhotonDirection(3)*PhotonProps%PhotonLastPos(3))/NewYPho
PhotonProps%PhotonDirection(3) = (-PhotonProps%PhotonDirection(2)*PhotonProps%PhotonLastPos(3) &
          + PhotonProps%PhotonDirection(3)*PhotonProps%PhotonLastPos(2))/NewYPho
PhotonProps%PhotonLastPos(2) = NewYPho
PhotonProps%PhotonLastPos(3) = 0.0
PhotonProps%PhotonDirection(2) = NewYVelo
PhotonProps%PhotonPos(1:3) = PhotonProps%PhotonLastPos(1:3)

END SUBROUTINE RotatePhotonIn2DPlane


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE PhotonIntersectionWithSide(iLocSide,Element,TriNum, IntersectionPos, PhotonLost, IsMortar)
USE MOD_Globals             ,ONLY: abort,UNIT_StdOut
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemSideNodeID_Shared, NodeCoords_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_Globals             ,ONLY: myRank
#endif /*USE_MPI*/
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER,INTENT(IN)                          :: iLocSide
INTEGER,INTENT(IN)                          :: Element
INTEGER,INTENT(IN)                          :: TriNum
LOGICAL,INTENT(OUT)                         :: PhotonLost
REAL,INTENT(OUT)                            :: IntersectionPos(1:3)
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
! Local variable declaration
INTEGER                          :: CNElemID
INTEGER                          :: Node1, Node2
REAL                             :: PoldX, PoldY, PoldZ, nx, ny, nz, nVal
REAL                             :: bx,by,bz, ax,ay,az, dist
REAL                             :: xNod, yNod, zNod,IntersecPara
REAL                             :: Vector1(1:3), Vector2(1:3), VectorShift(1:3)
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

CNElemID = GetCNElemID(Element)

! Sanity check
IF(CNElemID.LE.0) THEN
  IPWRITE(UNIT_StdOut,*) "Element  =", Element
  IPWRITE(UNIT_StdOut,*) "CNElemID =", CNElemID
  CALL abort(__STAMP__,'PhotonIntersectionWithSide() found CNElemID<=0')
END IF

PoldX = PhotonProps%PhotonLastPos(1)
PoldY = PhotonProps%PhotonLastPos(2)
PoldZ = PhotonProps%PhotonLastPos(3)

xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)

!---- Calculate normal vector:
IF(PRESENT(IsMortar)) THEN
  Node1 = TriNum+2     ! normal = cross product of 1-2 and 1-3 for first triangle
  Node2 = TriNum+1     !          and 1-3 and 1-4 for second triangle
ELSE
  Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
  Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
END IF

Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - xNod
Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - yNod
Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - zNod

Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - xNod
Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - yNod
Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - zNod

! Normal vector of side (one of the two triangles)
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

! If vector from old point to new point goes through the node, a will be x_photon_startro
! dist is then simply length of vector b instead of |axb|/|a|
IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

VectorShift(1) = PhotonProps%PhotonDirection(1)
VectorShift(2) = PhotonProps%PhotonDirection(2)
VectorShift(3) = PhotonProps%PhotonDirection(3)

! IntersecPara is zero when photon is perfectly perpendicular to side (cannot actually happen)
PhotonLost=.FALSE.
IntersecPara = VectorShift(1) * nx + VectorShift(2) * ny + VectorShift(3) * nz
IF(ABS(IntersecPara).GT.0.0)THEN
  IntersecPara = dist / IntersecPara
ELSE
  ! Mark photon as lost so that it is written to PartStateLost.h5
  PhotonLost      = .TRUE.
  IntersectionPos = HUGE(1.0)
  RETURN
END IF ! ABS(IntersecPara).GT.0.0

IntersectionPos(1) = PoldX + IntersecPara * VectorShift(1)
IntersectionPos(2) = PoldY + IntersecPara * VectorShift(2)
IntersectionPos(3) = PoldZ + IntersecPara * VectorShift(3)

RETURN
END SUBROUTINE PhotonIntersectionWithSide


!===================================================================================================================================
!> Calculate the absorbed energy in the (sub-)volume of each element. For higher-order sampling, the ray path between element entry
!> and exit is sampled NbrOfSamples = MAX(30,(Nloc+1)**2) and a nearest neighbour search finds the nearest sub-volume element on
!> which the energy is deposited. The sum of all sub-volume energies must equal the element-constant value, which is also
!> determined. Not that only the energies are conserved and not the density, which is the ratio of energy and (sub-)volume when
!> a change basis is used to switch between point sets (polynomial representations).
!===================================================================================================================================
SUBROUTINE CalcAbsorptionRayTrace(IntersectionPos,GlobalElemID,PhotonDir)
USE MOD_Globals             ,ONLY: VECNORM
USE MOD_RayTracing_Vars     ,ONLY: RayElemPassedEnergy,Ray,U_N_Ray,N_DG_Ray
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER, INTENT(IN) :: GlobalElemID
REAL, INTENT(IN)    :: PhotonDir(3)
REAL, INTENT(IN)    :: IntersectionPos(3)
! Local variable declaration
INTEGER           :: k,l,m,Nloc,NbrOfSamples,iIntersec,idx
REAL              :: SamplePos(3)
REAL              :: direction(3),subdirection(3),length,sublength
REAL              :: RandVal(3)
!--------------------------------------------------------------------------------------------------!
! Calculate the direction and length of the path of the ray through the element
direction(1:3) = IntersectionPos(1:3)-PhotonProps%PhotonPos(1:3)
length = VECNORM(direction(1:3))

! Check primary or secondary direction
IF(DOT_PRODUCT(PhotonDir,Ray%Direction).GT.0.0)THEN
  ! 1st energy direction
  idx = 1
ELSE
  ! 2nd energy direction
  idx = 2
  RayElemPassedEnergy(3:5,GlobalElemID) = RayElemPassedEnergy(3:5,GlobalElemID) + PhotonDir(1:3)
  RayElemPassedEnergy(6,GlobalElemID)   = RayElemPassedEnergy(6,GlobalElemID)   + 1.0
END IF
RayElemPassedEnergy(idx,GlobalElemID)   = RayElemPassedEnergy(idx,GlobalElemID) + PhotonProps%PhotonEnergy * length

! High-order sampling: Use nearest Gauss point (NGP) from PIC deposition
Nloc = N_DG_Ray(GlobalElemID)

! Loop over number of sub-samples
NbrOfSamples = MAX(30,(Nloc+1)**2) ! Nloc+1 ! must be at least 3 for this sampling method (one point between the two intersections of the element)!
subdirection(1:3) = direction(1:3)/REAL(NbrOfSamples-1)
sublength = VECNORM(subdirection(1:3))
! Loop over the number of sub lengths and assign them to the nearest DOF. Choose the intersection points at random to prevent artifacts
IF(ABS(direction(3)).GT.1e6*(ABS(direction(1))+ABS(direction(2))))THEN
  ! only in z-dir
  DO iIntersec = 1, NbrOfSamples-1
    CALL RANDOM_NUMBER(RandVal(1))
    SamplePos(1:3) = PhotonProps%PhotonPos(1:3) + (/1.0 , 1.0 , RandVal(1)/) * direction(1:3)
    CALL GetNestestDOFInRefElem(Nloc,SamplePos(1:3),GlobalElemID,k,l,m)
    U_N_Ray(GlobalElemID)%U(idx,k,l,m) = U_N_Ray(GlobalElemID)%U(idx,k,l,m) + sublength*PhotonProps%PhotonEnergy
  END DO
ELSE
  DO iIntersec = 1, NbrOfSamples-1
    CALL RANDOM_NUMBER(RandVal(1:3))
    SamplePos(1:3) = PhotonProps%PhotonPos(1:3) + RandVal(1:3) * direction(1:3)
    CALL GetNestestDOFInRefElem(Nloc,SamplePos(1:3),GlobalElemID,k,l,m)
    U_N_Ray(GlobalElemID)%U(idx,k,l,m) = U_N_Ray(GlobalElemID)%U(idx,k,l,m) + sublength*PhotonProps%PhotonEnergy
  END DO
END IF ! ABS(direction(3).GT.1e6*(ABS(direction(1))+ABS(direction(2))))
PhotonProps%PhotonPos(1:3) = IntersectionPos(1:3)

! Store intersection point as new starting point
PhotonProps%PhotonPos(1:3) = IntersectionPos(1:3)

END SUBROUTINE CalcAbsorptionRayTrace


!===================================================================================================================================
!> Find the nearest DOF of a particle position (or sampled particle path coordinate) and return k,l,m
!===================================================================================================================================
SUBROUTINE GetNestestDOFInRefElem(Nloc,SamplePos,GlobalElemID,k,l,m)
! MODULES
USE MOD_globals
USE MOD_Eval_xyz        ,ONLY: GetPositionInRefElem
USE MOD_RayTracing_Vars ,ONLY: N_Inter_Ray
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: SamplePos(1:3)
INTEGER,INTENT(IN)  :: Nloc
INTEGER,INTENT(IN)  :: GlobalElemID
INTEGER,INTENT(OUT) :: k,l,m
!INTEGER             :: kOld,lOld,mOld
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: IntersectionPosRef(3)
!INTEGER :: ii,a,b
!===================================================================================================================================
! OLD METHOD FROM NEAREST GAUSS POINT

!IF(MOD(Nloc,2).EQ.0) THEN
!  a = Nloc/2
!  b = a
!ELSE
!  a = (Nloc+1)/2
!  b = a-1
!END IF
!
!! Get position in reference element
!CALL GetPositionInRefElem(SamplePos(1:3),IntersectionPosRef(1:3),GlobalElemID)
!
!k = a
!DO ii = 0,b-1
!  IF(ABS(IntersectionPosRef(1)).GE.N_Inter_Ray(Nloc)%GaussBorder(Nloc-ii))THEN
!    k = Nloc-ii
!    EXIT
!  END IF
!END DO
!k = NINT((Nloc+SIGN(2.0*k-Nloc,IntersectionPosRef(1)))/2)
!!! y-direction
!l = a
!DO ii = 0,b-1
!  IF(ABS(IntersectionPosRef(2)).GE.N_Inter_Ray(Nloc)%GaussBorder(Nloc-ii))THEN
!    l = Nloc-ii
!    EXIT
!  END IF
!END DO
!l = NINT((Nloc+SIGN(2.0*l-Nloc,IntersectionPosRef(2)))/2)
!!! z-direction
!m = a
!DO ii = 0,b-1
!  IF(ABS(IntersectionPosRef(3)).GE.N_Inter_Ray(Nloc)%GaussBorder(Nloc-ii))THEN
!    m = Nloc-ii
!    EXIT
!  END IF
!END DO
!m = NINT((Nloc+SIGN(2.0*m-Nloc,IntersectionPosRef(3)))/2)
!
!kOld = k
!lOld = l
!mOld = m

! Get reference position
CALL GetPositionInRefElem(SamplePos(1:3),IntersectionPosRef(1:3),GlobalElemID)

k = MINLOC(ABS(N_Inter_Ray(Nloc)%xGP(:) - IntersectionPosRef(1)),DIM=1) - 1
l = MINLOC(ABS(N_Inter_Ray(Nloc)%xGP(:) - IntersectionPosRef(2)),DIM=1) - 1
m = MINLOC(ABS(N_Inter_Ray(Nloc)%xGP(:) - IntersectionPosRef(3)),DIM=1) - 1
!IPWRITE(UNIT_StdOut,*) "Nloc =", Nloc
!IPWRITE(UNIT_StdOut,*) "N_Inter_Ray(Nloc)%xGP(:) =", N_Inter_Ray(Nloc)%xGP(:)
!IPWRITE(UNIT_StdOut,*) "IntersectionPosRef =", IntersectionPosRef
!
!IF(kOld.ne.k .or. lOld.ne.l .or. mOld.ne.m)THEN
!  write(*,*) ""
!  IPWRITE(UNIT_StdOut , *) "k    , l    , m    =" , k    , l    , m
!  IPWRITE(UNIT_StdOut , *) "kOld , lOld , mOld =" , kOld , lOld , mOld
!  CALL abort(__STAMP__,'not equal!')
!END IF ! kOld.eq.k .or. lOld.ne.l .or. mOld.ne.m
!
!  IPWRITE(UNIT_StdOut , *) "k    , l    , m    ="    , k    , l    , m
!  !IF(myrank.eq.0) read*; CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)

END SUBROUTINE GetNestestDOFInRefElem


!===================================================================================================================================
!> Calculates absorbed energy of photons along their paths stochastically
!===================================================================================================================================
SUBROUTINE CalcAbsoprtionMC(IntersectionPos,Element, DONE)
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergy
USE MOD_Radiation_Vars      ,ONLY: Radiation_Absorption_spec,Radiation_Absorption_SpecPercent
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergySpec
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)              :: Element
REAL, INTENT(IN)                 :: IntersectionPos(3)
LOGICAL, INTENT(OUT)             :: DONE
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
REAL                            :: iRan, DistanceVec(3), Distance, opticalPath
!--------------------------------------------------------------------------------------------------!
IF ((Radiation_Absorption_Spec(PhotonProps%WaveLength,Element).GT.0.0)&
    .AND.(SUM(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element)).GT.0.0)) THEN
  DistanceVec(1:3) = PhotonProps%PhotonPos(1:3) - IntersectionPos(1:3)
  Distance = SQRT(DistanceVec(1)*DistanceVec(1) + DistanceVec(2)*DistanceVec(2) + DistanceVec(3)*DistanceVec(3))
  CALL RANDOM_NUMBER(iRan)
  opticalPath = Distance*Radiation_Absorption_Spec(PhotonProps%WaveLength,Element)
  IF (-LOG(iRan).LT.opticalPath) THEN
    RadiationElemAbsEnergySpec(:,Element) = RadiationElemAbsEnergySpec(:,Element) &
    + PhotonProps%PhotonEnergy*(REAL(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element))&
    /SUM(REAL(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element))))
    DONE = .TRUE.
  ELSE
    PhotonProps%PhotonPos(1:3) = IntersectionPos(1:3)
  END IF
ELSE
  PhotonProps%PhotonPos(1:3) = IntersectionPos(1:3)
END IF
RadiationElemAbsEnergy(1,Element) = RadiationElemAbsEnergy(1,Element) + opticalPath
RadiationElemAbsEnergy(2,Element) = RadiationElemAbsEnergy(2,Element) + 1.0

END SUBROUTINE CalcAbsoprtionMC


!===================================================================================================================================
!> Calculates absorbed energy of photons along their paths analytically
!===================================================================================================================================
SUBROUTINE CalcAbsoprtionAnalytic(IntersectionPos,Element, DONE)
USE MOD_Globals
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergy, RadiationElemAbsEnergySpec
USE MOD_Radiation_Vars      ,ONLY: Radiation_Absorption_spec, Radiation_Absorption_SpecPercent
!--------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER, INTENT(IN)              :: Element
REAL, INTENT(IN)                 :: IntersectionPos(3)
LOGICAL, INTENT(INOUT)             :: DONE
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
REAL                            :: DistanceVec(3), Distance, LostEnergy, opticalPath
!--------------------------------------------------------------------------------------------------!
IF ((Radiation_Absorption_Spec(PhotonProps%WaveLength,Element).GT.0.0)&
    .AND.(SUM(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element)).GT.0)) THEN
  DistanceVec(1:3) = PhotonProps%PhotonPos(1:3) - IntersectionPos(1:3)
  Distance = SQRT(DistanceVec(1)*DistanceVec(1) + DistanceVec(2)*DistanceVec(2) + DistanceVec(3)*DistanceVec(3))
  opticalPath = Distance*Radiation_Absorption_Spec(PhotonProps%WaveLength,Element)
  IF (CHECKEXP(opticalPath)) THEN
    LostEnergy = PhotonProps%PhotonEnergy*(1.-EXP(-opticalPath))
  ELSE
    LostEnergy = PhotonProps%PhotonEnergy
    DONE = .TRUE.
  END IF
  IF (SUM(REAL(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element))).EQ.0.0) THEN
      print*,'arg',Element,PhotonProps%WaveLength,  Radiation_Absorption_Spec(PhotonProps%WaveLength,Element)
      print*, 'percent', Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element)
  END IF
  PhotonProps%PhotonEnergy = PhotonProps%PhotonEnergy - LostEnergy
  RadiationElemAbsEnergySpec(:,Element) = RadiationElemAbsEnergySpec(:,Element) &
    + LostEnergy*(REAL(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element))&
    /SUM(REAL(Radiation_Absorption_SpecPercent(PhotonProps%WaveLength,:,Element))))
ELSE
  opticalPath = 0.0
END IF
!  IF (PhotonProps%PhotonEnergy.LE.(RadTrans%GlobalRadiationPower/(1000.*RadTrans%GlobalPhotonNum))) THEN
!    DONE = .TRUE.
!  ELSE

!  END IF
PhotonProps%PhotonPos(1:3) = IntersectionPos(1:3)
RadiationElemAbsEnergy(1,Element) = RadiationElemAbsEnergy(1,Element) + opticalPath
RadiationElemAbsEnergy(2,Element) = RadiationElemAbsEnergy(2,Element) + 1.0

END SUBROUTINE CalcAbsoprtionAnalytic


!===================================================================================================================================
!> Calculates absorbed energy of photons along their paths
!===================================================================================================================================
SUBROUTINE CalcAbsoprtion(IntersectionPos, Element, DONE, before)
USE MOD_Globals
USE MOD_RadiationTrans_Vars ,ONLY: RadiationAbsorptionModel
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER, INTENT(IN)           :: Element
REAL, INTENT(IN)              :: IntersectionPos(3)
LOGICAL, INTENT(INOUT)        :: DONE
LOGICAL, INTENT(IN), OPTIONAL :: before !> before=T (before reflection), before=F (after reflection), before=NOT PRESENT (volume)
! Local variable declaration
!--------------------------------------------------------------------------------------------------!
IF (RadiationAbsorptionModel.EQ.0) THEN
  ! For ray tracing, check if routine is called be
  IF(PRESENT(before))THEN
    IF(before)THEN
      ! Before reflection: Nothing to do as it is not yet known if the ray is absorbed at the surface
      RETURN
    ELSE
      ! After reflection: Use old or new ray direction depending on whether the ray was absorbed
      IF(DONE)THEN
        ! Ray was absorbed at the wall
        CALL CalcAbsorptionRayTrace(IntersectionPos, Element,PhotonProps%PhotonDirectionBeforeReflection)
      ELSE
        ! Ray was reflected at the wall
        ! TODO: Not sure which ray vector should be used
        !CALL CalcAbsorptionRayTrace(Element,PhotonProps%PhotonDirection)
        CALL CalcAbsorptionRayTrace(IntersectionPos, Element,PhotonProps%PhotonDirectionBeforeReflection)
      END IF ! DONE
    END IF ! before
  ELSE
    CALL CalcAbsorptionRayTrace(IntersectionPos, Element,PhotonProps%PhotonDirection)
  END IF ! PRESENT(before)
ELSEIF (RadiationAbsorptionModel.EQ.1) THEN
  CALL CalcAbsoprtionAnalytic(IntersectionPos,Element, DONE)
ELSEIF (RadiationAbsorptionModel.EQ.2) THEN
  CALL CalcAbsoprtionMC(IntersectionPos,Element, DONE)
ELSE
  CALL Abort(__STAMP__,'AbsorptionModel must be 1 or 2!')
END IF

END SUBROUTINE CalcAbsoprtion


!===================================================================================================================================
!> Determines velocity vectors of photons after a perfect reflection at a boundary
!===================================================================================================================================
SUBROUTINE PerfectPhotonReflection(iLocSide,Element,TriNum, IntersectionPos, IntersecAlreadyCalc)
USE MOD_Particle_Mesh_Vars  ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER,INTENT(IN)              :: iLocSide
INTEGER,INTENT(IN)              :: Element
INTEGER,INTENT(IN)              :: TriNum
REAL, INTENT(INOUT)             :: IntersectionPos(1:3)
LOGICAL, INTENT(IN)             :: IntersecAlreadyCalc
! Local variable declaration
INTEGER                          :: CNElemID
INTEGER                          :: Node1, Node2
REAL                             :: PoldX, PoldY, PoldZ, nx, ny, nz, nVal
REAL                             :: xNod, yNod, zNod
REAL                             :: VelX, VelY, VelZ
REAL                             :: Vector1(1:3), Vector2(1:3), POI_fak, ProjVel
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

CNElemID = GetCNElemID(Element)
PoldX = PhotonProps%PhotonLastPos(1)
PoldY = PhotonProps%PhotonLastPos(2)
PoldZ = PhotonProps%PhotonLastPos(3)

VelX = PhotonProps%PhotonDirection(1)
VelY = PhotonProps%PhotonDirection(2)
VelZ = PhotonProps%PhotonDirection(3)

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

!---- Calculate Point of Intersection (POI)
IF (.NOT.IntersecAlreadyCalc) THEN
  POI_fak = (Vector2(2)*(Vector1(1)*(zNod-PoldZ)+Vector1(3)*(PoldX-xNod)) &
          +Vector1(2)*(Vector2(1)*(PoldZ-zNod)+Vector2(3)*(xNod-PoldX)) &
          +yNod*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)) &
          +PoldY*(Vector1(1)*Vector2(3)-Vector1(3)*Vector2(1))) &
          /(Vector1(2)*(Vector2(3)*VelX-Vector2(1)*VelZ) &
          + Vector2(2)*(Vector1(1)*VelZ-Vector1(3)*VelX) &
          + VelY*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)))

  IntersectionPos(1) = PoldX + POI_fak * VelX
  IntersectionPos(2) = PoldY + POI_fak * VelY
  IntersectionPos(3) = PoldZ + POI_fak * VelZ
END IF

!---- Calculate new velocity vector
ProjVel = nx*PhotonProps%PhotonDirection(1)+ny*PhotonProps%PhotonDirection(2) &
     +nz*PhotonProps%PhotonDirection(3)
VelX=PhotonProps%PhotonDirection(1)-2.*ProjVel*nx
VelY=PhotonProps%PhotonDirection(2)-2.*ProjVel*ny
VelZ=PhotonProps%PhotonDirection(3)-2.*ProjVel*nz

!---- Assign new values to "old" variables to continue loop

PhotonProps%PhotonLastPos(1) = IntersectionPos(1)
PhotonProps%PhotonLastPos(2) = IntersectionPos(2)
PhotonProps%PhotonLastPos(3) = IntersectionPos(3)

PhotonProps%PhotonDirection(1) = VelX
PhotonProps%PhotonDirection(2) = VelY
PhotonProps%PhotonDirection(3) = VelZ
RETURN
END SUBROUTINE PerfectPhotonReflection


SUBROUTINE PerfectPhotonReflection2D(iLocSide,Element, IntersectionPos)
!--------------------------------------------------------------------------------------------------!
!Based on PerfectReflection3D
!--------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars  ,ONLY: SideNormalEdge2D_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER,INTENT(IN)              :: iLocSide
INTEGER,INTENT(IN)              :: Element
REAL, INTENT(INOUT)             :: IntersectionPos(1:3)
! Local variable declaration
INTEGER                          :: CNElemID
REAL                             :: nx, ny, nz, nValIntersec
REAL                             :: VelX, VelY, VelZ
REAL                             :: ProjVel
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
CNElemID = GetCNElemID(Element)

nx = SideNormalEdge2D_Shared(1,iLocSide, CNElemID)
nValIntersec = SQRT(IntersectionPos(2)*IntersectionPos(2) + IntersectionPos(3)*IntersectionPos(3))
ny = IntersectionPos(2)/nValIntersec * SideNormalEdge2D_Shared(2,iLocSide, CNElemID)
nz = IntersectionPos(3)/nValIntersec * SideNormalEdge2D_Shared(2,iLocSide, CNElemID)

!---- Calculate new velocity vector
ProjVel = nx*PhotonProps%PhotonDirection(1)+ny*PhotonProps%PhotonDirection(2) &
         +nz*PhotonProps%PhotonDirection(3)
VelX=PhotonProps%PhotonDirection(1)-2.*ProjVel*nx
VelY=PhotonProps%PhotonDirection(2)-2.*ProjVel*ny
VelZ=PhotonProps%PhotonDirection(3)-2.*ProjVel*nz

!---- Assign new values to "old" variables to continue loop

PhotonProps%PhotonLastPos(1) = IntersectionPos(1)
PhotonProps%PhotonLastPos(2) = IntersectionPos(2)
PhotonProps%PhotonLastPos(3) = IntersectionPos(3)

PhotonProps%PhotonDirection(1)   = VelX
PhotonProps%PhotonDirection(2)   = VelY
PhotonProps%PhotonDirection(3)   = VelZ
END SUBROUTINE PerfectPhotonReflection2D


SUBROUTINE DiffusePhotonReflection(iLocSide,Element,TriNum, IntersectionPos, IntersecAlreadyCalc)
!--------------------------------------------------------------------------------------------------!
!Based on PerfectReflection3D
!--------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemSideNodeID_Shared, NodeCoords_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE Ziggurat
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER,INTENT(IN)              :: iLocSide
INTEGER,INTENT(IN)              :: Element
INTEGER,INTENT(IN)              :: TriNum
REAL, INTENT(INOUT)             :: IntersectionPos(1:3)
LOGICAL, INTENT(IN)             :: IntersecAlreadyCalc
! Local variable declaration
INTEGER                          :: CNElemID
INTEGER                          :: Node1, Node2
REAL                             :: PoldX, PoldY, PoldZ, nx, ny, nz, nVal
REAL                             :: xNod, yNod, zNod, VecX, VecY, VecZ
REAL                             :: VelX, VelY, VelZ, VeloCx, VeloCy, VeloCz, NormVec, RanNum
REAL                             :: Vector1(1:3), Vector2(1:3), POI_fak
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
CNElemID = GetCNElemID(Element)
PoldX = PhotonProps%PhotonLastPos(1)
PoldY = PhotonProps%PhotonLastPos(2)
PoldZ = PhotonProps%PhotonLastPos(3)

VelX = PhotonProps%PhotonDirection(1)
VelY = PhotonProps%PhotonDirection(2)
VelZ = PhotonProps%PhotonDirection(3)

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

!---- Calculate Point of Intersection (POI)
!---- Calculate Point of Intersection (POI)
IF (.NOT.IntersecAlreadyCalc) THEN
  POI_fak = (Vector2(2)*(Vector1(1)*(zNod-PoldZ)+Vector1(3)*(PoldX-xNod)) &
          +Vector1(2)*(Vector2(1)*(PoldZ-zNod)+Vector2(3)*(xNod-PoldX)) &
          +yNod*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)) &
          +PoldY*(Vector1(1)*Vector2(3)-Vector1(3)*Vector2(1))) &
          /(Vector1(2)*(Vector2(3)*VelX-Vector2(1)*VelZ) &
          + Vector2(2)*(Vector1(1)*VelZ-Vector1(3)*VelX) &
          + VelY*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)))

  IntersectionPos(1) = PoldX + POI_fak * VelX
  IntersectionPos(2) = PoldY + POI_fak * VelY
  IntersectionPos(3) = PoldZ + POI_fak * VelZ
END IF
!---- Calculate new velocity vector (Extended Maxwellian Model)

VeloCx  = rnor()                !normal distri
VeloCy  = rnor()                !normal distri
CALL RANDOM_NUMBER(RanNum)
VeloCz  = SQRT(-2.*LOG(RanNum)) ! rayleigh distri

!---- Transformation local distribution -> global coordinates
VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
VecZ = Vector1(3) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )

VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz
!---- Assign new values to "old" variables to continue loop

PhotonProps%PhotonLastPos(1) = IntersectionPos(1)
PhotonProps%PhotonLastPos(2) = IntersectionPos(2)
PhotonProps%PhotonLastPos(3) = IntersectionPos(3)

!----  saving new particle velocity
NormVec = SQRT(VelX*VelX + VelY*VelY + VelZ*VelZ)
PhotonProps%PhotonDirection(1)   = VelX / NormVec
PhotonProps%PhotonDirection(2)   = VelY / NormVec
PhotonProps%PhotonDirection(3)   = VelZ / NormVec

RETURN
END SUBROUTINE DiffusePhotonReflection


SUBROUTINE DiffusePhotonReflection2D(iLocSide,Element, IntersectionPos)
!--------------------------------------------------------------------------------------------------!
!Based on PerfectReflection3D
!--------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars  ,ONLY: SideNormalEdge2D_Shared
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
USE Ziggurat
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
INTEGER,INTENT(IN)              :: iLocSide
INTEGER,INTENT(IN)              :: Element
REAL, INTENT(IN)                :: IntersectionPos(1:3)
! Local variable declaration
INTEGER                          :: CNElemID
REAL                             :: nx, ny, nz, nValIntersec, VecX, VecY, VecZ
REAL                             :: VelX, VelY, VelZ, VeloCx, VeloCy, VeloCz, NormVec, RanNum
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
CNElemID = GetCNElemID(Element)
nx = SideNormalEdge2D_Shared(1,iLocSide, CNElemID)
nValIntersec = SQRT(IntersectionPos(2)*IntersectionPos(2) + IntersectionPos(3)*IntersectionPos(3))
ny = IntersectionPos(2)/nValIntersec * SideNormalEdge2D_Shared(2,iLocSide, CNElemID)
nz = IntersectionPos(3)/nValIntersec * SideNormalEdge2D_Shared(2,iLocSide, CNElemID)

VecX = SideNormalEdge2D_Shared(3,iLocSide, CNElemID)
VecY = IntersectionPos(2)/nValIntersec * SideNormalEdge2D_Shared(4,iLocSide, CNElemID)
VecZ = IntersectionPos(3)/nValIntersec * SideNormalEdge2D_Shared(4,iLocSide, CNElemID)
!---- Calculate new velocity vector (Extended Maxwellian Model)

VeloCx  = rnor()                !normal distri
VeloCy  = rnor()                !normal distri
CALL RANDOM_NUMBER(RanNum)
VeloCz  = SQRT(-2.*LOG(RanNum)) ! rayleigh distri

VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz
!---- Assign new values to "old" variables to continue loop

PhotonProps%PhotonLastPos(1) = IntersectionPos(1)
PhotonProps%PhotonLastPos(2) = IntersectionPos(2)
PhotonProps%PhotonLastPos(3) = IntersectionPos(3)

!----  saving new particle velocity
NormVec = SQRT(VelX*VelX + VelY*VelY + VelZ*VelZ)
PhotonProps%PhotonDirection(1)   = VelX / NormVec
PhotonProps%PhotonDirection(2)   = VelY / NormVec
PhotonProps%PhotonDirection(3)   = VelZ / NormVec

END SUBROUTINE DiffusePhotonReflection2D


SUBROUTINE PeriodicPhotonBC(iLocSide, Element, TriNum, IntersectionPos, IntersecAlreadyCalc, SideID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the periodic movement of a photon and updates PhotonLastPos, PhotonPos and IntersectionPos to the periodically moves
! location
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemSideNodeID_Shared, NodeCoords_Shared
USE MOD_Mesh_Vars              ,ONLY: BoundaryType
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Photon_TrackingVars    ,ONLY: PhotonProps
!USE MOD_Particle_Tracking_Vars ,ONLY: TrackInfo
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: TriNum
INTEGER,INTENT(IN)               :: SideID
LOGICAL,INTENT(IN)               :: IntersecAlreadyCalc
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: Element
REAL, INTENT(INOUT)              :: IntersectionPos(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: PVID
!INTEGER                              :: moved(2),locSideID
! Local variable declaration
INTEGER                          :: CNElemID
INTEGER                          :: Node1, Node2
REAL                             :: PoldX, PoldY, PoldZ
REAL                             :: xNod, yNod, zNod
REAL                             :: VelX, VelY, VelZ
REAL                             :: Vector1(1:3), Vector2(1:3), POI_fak
!===================================================================================================================================
CNElemID = GetCNElemID(Element)
PoldX = PhotonProps%PhotonLastPos(1)
PoldY = PhotonProps%PhotonLastPos(2)
PoldZ = PhotonProps%PhotonLastPos(3)

VelX = PhotonProps%PhotonDirection(1)
VelY = PhotonProps%PhotonDirection(2)
VelZ = PhotonProps%PhotonDirection(3)

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

!---- Calculate Point of Intersection (POI)
IF (.NOT.IntersecAlreadyCalc) THEN
  POI_fak = (Vector2(2)*(Vector1(1)*(zNod-PoldZ)+Vector1(3)*(PoldX-xNod)) &
          +Vector1(2)*(Vector2(1)*(PoldZ-zNod)+Vector2(3)*(xNod-PoldX)) &
          +yNod*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)) &
          +PoldY*(Vector1(1)*Vector2(3)-Vector1(3)*Vector2(1))) &
          /(Vector1(2)*(Vector2(3)*VelX-Vector2(1)*VelZ) &
          + Vector2(2)*(Vector1(1)*VelZ-Vector1(3)*VelX) &
          + VelY*(Vector1(3)*Vector2(1)-Vector1(1)*Vector2(3)))

  IntersectionPos(1) = PoldX + POI_fak * VelX
  IntersectionPos(2) = PoldY + POI_fak * VelY
  IntersectionPos(3) = PoldZ + POI_fak * VelZ
END IF

! Set last particle position on face
PhotonProps%PhotonLastPos = IntersectionPos
! Perform the periodic movement
PVID = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)
PhotonProps%PhotonLastPos = PhotonProps%PhotonLastPos + SIGN( GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
! Update particle positon after periodic BC
!PartState(1:3,PartID) = PhotonProps%PhotonLastPos + (TrackInfo%lengthPartTrajectory-TrackInfo%alpha)*TrackInfo%PartTrajectory
!TrackInfo%lengthPartTrajectory  = TrackInfo%lengthPartTrajectory - TrackInfo%alpha

! refmapping and tracing
! Move particle from old element to new element
Element = SideInfo_Shared(SIDE_NBELEMID,SideID)

! Periodic movement the intersection point as it will be considered the starting point in the next iteration
IntersectionPos = PhotonProps%PhotonLastPos

! Periodic movement the photon position
PhotonProps%PhotonPos = PhotonProps%PhotonLastPos

END SUBROUTINE PeriodicPhotonBC


!===================================================================================================================================
!> Photon interactions with walls:
!> Compare random number with photon accommodation and absorb/reflect depending on the outcome
!>
!>   ForceWallSample (OPTIONAL): When true, the sampling is performed independent of the actual absorption/reflection outcome
!===================================================================================================================================
SUBROUTINE CalcWallAbsoprtion(IntersectionPos, GlobSideID, DONE, ForceWallSample)
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Photon_TrackingVars    ,ONLY: PhotonProps,PhotonSurfSideSamplingMidPoints
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, GlobalSide2SurfSide
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
#if USE_MPI
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSampWallProc
#else
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSampWall
#endif /*USE_MPI*/
USE MOD_RayTracing_Vars        ,ONLY: Ray
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! argument list declaration
REAL, INTENT(IN)                 :: IntersectionPos(3)
INTEGER, INTENT(IN)              :: GlobSideID
LOGICAL, INTENT(OUT)             :: DONE
LOGICAL, INTENT(IN), OPTIONAL    :: ForceWallSample !>
! Local variable declaration
!--------------------------------------------------------------------------------------------------!
REAL                            :: iRan,PhotonEnACC,distance,distanceMin
INTEGER                         :: SurfSideID,p,q,pp,qq
LOGICAL                         :: ForceWallSampleLoc
!--------------------------------------------------------------------------------------------------!
#if USE_MPI
ASSOCIATE( PhotonSampWall => PhotonSampWallProc )
#endif /*USE_MPI*/
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,GlobSideID)
! Check if photon is to be added to PhotonSampWall independent of the actual absorption/reflection
IF(PRESENT(ForceWallSample))THEN
  ForceWallSampleLoc = ForceWallSample
  ! Sample impact
  IF(ForceWallSampleLoc)THEN
    IF(Ray%nSurfSample.GT.1)THEN
      distanceMin = HUGE(1.)
      DO pp = 1, Ray%nSurfSample
        DO qq = 1, Ray%nSurfSample
          distance = VECNORM(IntersectionPos(1:3) - PhotonSurfSideSamplingMidPoints(1:3,pp,qq,SurfSideID))
          IF(distance.LT.distanceMin)THEN
            p = pp
            q = qq
            distanceMin = distance
          END IF ! distance.LT.distanceMin
        END DO ! q = 1, Ray%nSurfSample
      END DO ! p = 1, Ray%nSurfSample
      PhotonSampWall(1,p,q,SurfSideID) = PhotonSampWall(1,p,q,SurfSideID) + 1.0
      PhotonSampWall(2,p,q,SurfSideID) = PhotonSampWall(2,p,q,SurfSideID) + PhotonProps%PhotonEnergy
    ELSE
      PhotonSampWall(1,1,1,SurfSideID) = PhotonSampWall(1,1,1,SurfSideID) + 1.0
      PhotonSampWall(2,1,1,SurfSideID) = PhotonSampWall(2,1,1,SurfSideID) + PhotonProps%PhotonEnergy
    END IF ! Ray%nSurfSample.GT.1
  END IF ! ForceWallSampleLoc
ELSE
  ForceWallSampleLoc = .FALSE.
END IF ! PRESENT(ForceWallSample)
DONE = .FALSE. ! initialize
PhotonEnACC = PartBound%PhotonEnACC(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobSideID)))
IF(PhotonEnACC.LE.0.0) RETURN ! Skip sides without absorption (pure reflection)
CALL RANDOM_NUMBER(iRan)
IF (PhotonEnACC.GT.iRan) THEN
  DONE = .TRUE.
  ! Do not sample twice
  IF(.NOT.ForceWallSampleLoc)THEN
    IF(Ray%nSurfSample.GT.1)THEN
      distanceMin = HUGE(1.)
      DO pp = 1, Ray%nSurfSample
        DO qq = 1, Ray%nSurfSample
          distance = VECNORM(IntersectionPos(1:3) - PhotonSurfSideSamplingMidPoints(1:3,pp,qq,SurfSideID))
          IF(distance.LT.distanceMin)THEN
            p = pp
            q = qq
          END IF ! distance.LT.distanceMin
        END DO ! q = 1, Ray%nSurfSample
      END DO ! p = 1, Ray%nSurfSample
      PhotonSampWall(1,p,q,SurfSideID) = PhotonSampWall(1,p,q,SurfSideID) + 1.0
      PhotonSampWall(2,p,q,SurfSideID) = PhotonSampWall(2,p,q,SurfSideID) + PhotonProps%PhotonEnergy
    ELSE
      PhotonSampWall(1,1,1,SurfSideID) = PhotonSampWall(1,1,1,SurfSideID) + 1.0
      PhotonSampWall(2,1,1,SurfSideID) = PhotonSampWall(2,1,1,SurfSideID) + PhotonProps%PhotonEnergy
    END IF ! Ray%nSurfSample.GT.1
  END IF ! .NOT.ForceWallSampleLoc
END IF
#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/
END SUBROUTINE CalcWallAbsoprtion


LOGICAL FUNCTION PointInObsCone(Point)
!===================================================================================================================================
! Checks if a point is in the opening cone of an external observer
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_RadiationTrans_Vars,    ONLY: RadObservationPoint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
REAL, INTENT(IN)             :: Point(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: ConeDist, ConeRadius, orthoDist
!===================================================================================================================================
PointInObsCone = .FALSE.
ConeDist = DOT_PRODUCT(Point(1:3) - RadObservationPoint%StartPoint(1:3), RadObservationPoint%ViewDirection(1:3))
ConeRadius = TAN(RadObservationPoint%AngularAperture/2.) * ConeDist
orthoDist = VECNORM(Point(1:3) - RadObservationPoint%StartPoint(1:3) - ConeDist*RadObservationPoint%ViewDirection(1:3))
IF (orthoDist.LE.ConeRadius) PointInObsCone = .TRUE.

END FUNCTION PointInObsCone


LOGICAL FUNCTION PhotonIntersectSensor(Point, Direction)
!===================================================================================================================================
! Checks if the photon's apth intersect with the opening cone of an external observer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RadiationTrans_Vars,    ONLY: RadObservationPoint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
REAL, INTENT(IN)             :: Point(3), Direction(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: projectedDist, DirectionVec(3)
!===================================================================================================================================
PhotonIntersectSensor = .FALSE.

projectedDist = DOT_PRODUCT(RadObservationPoint%ViewDirection(1:3), Direction(1:3))
IF (projectedDist.LT.0.0) THEN
  DirectionVec(1:3) =  RadObservationPoint%MidPoint(1:3) - Point(1:3)
  !distance to travel
  projectedDist = DOT_PRODUCT(DirectionVec(1:3), RadObservationPoint%ViewDirection(1:3))/projectedDist
  ! actual intersection point
  DirectionVec(1:3) = Point(1:3) + projectedDist*Direction(1:3)
  !Vector from midpoint of sensor
  DirectionVec(1:3) = DirectionVec(1:3) - RadObservationPoint%MidPoint(1:3)
  !distance to midpoint
  projectedDist = VECNORM(DirectionVec(1:3))
  IF (projectedDist.LE.RadObservationPoint%Diameter/2.) PhotonIntersectSensor = .TRUE.
END IF

END FUNCTION PhotonIntersectSensor


LOGICAL FUNCTION PhotonOnLineOfSight(Direction)
!===================================================================================================================================
! Checks if a photon is on the simulated line of sight
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RadiationTrans_Vars,    ONLY: RadObservationPoint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
REAL, INTENT(IN)             :: Direction(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: SkalarFactors(3)
INTEGER                       :: iDir, jDir
!===================================================================================================================================
PhotonOnLineOfSight = .FALSE.
DO iDir = 1, 3
  IF (Direction(iDir).EQ.0.0) THEN
    IF (RadObservationPoint%ViewDirection(iDir).NE.0.0) THEN
      RETURN
    ELSE
      SkalarFactors(iDir) = 0.0
    END IF
  ELSE
    IF (RadObservationPoint%ViewDirection(iDir).EQ.0.0) THEN
      RETURN
    ELSE
      SkalarFactors(iDir) = Direction(iDir)/ RadObservationPoint%ViewDirection(iDir)
    END IF
  END IF
END DO
PhotonOnLineOfSight = .TRUE.
DO iDir = 1, 2
  DO jDir = iDir+1 , 3
    IF (SkalarFactors(iDir).EQ.0.0) CYCLE
    IF (SkalarFactors(jDir).EQ.0.0) CYCLE
    IF (.NOT.ALMOSTEQUAL(SkalarFactors(iDir),SkalarFactors(jDir))) THEN
      PhotonOnLineOfSight = .FALSE.
      RETURN
    END IF
  END DO
END DO

END FUNCTION PhotonOnLineOfSight

END MODULE MOD_Photon_TrackingTools
