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

MODULE MOD_Particle_Mesh_Tools
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE BoundsOfElement
  MODULE PROCEDURE BoundsOfElement
END INTERFACE

INTERFACE ParticleInsideQuad3D
  MODULE PROCEDURE ParticleInsideQuad3D
END INTERFACE

INTERFACE ParticleInsideQuad3D_MortarMPI
  MODULE PROCEDURE ParticleInsideQuad3D_MortarMPI
END INTERFACE

INTERFACE GetGlobalElemID
  PROCEDURE GetGlobalElemID
END INTERFACE

INTERFACE GetCNElemID
  PROCEDURE GetCNElemID
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: GetGlobalElemID
PUBLIC :: GetCNElemID
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalElemIDInterface),POINTER :: GetGlobalElemID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNElemIDInterface),POINTER :: GetCNElemID    !< pointer defining the mapping: global element ID -> compute-node element ID

! Initialization routines
INTERFACE InitGetGlobalElemID
  MODULE PROCEDURE InitGetGlobalElemID
END INTERFACE

INTERFACE InitGetCNElemID
  MODULE PROCEDURE InitGetCNElemID
END INTERFACE

PUBLIC::InitGetGlobalElemID
PUBLIC::InitGetCNElemID
PUBLIC::BoundsOfElement, ParticleInsideQuad3D, ParticleInsideQuad3D_MortarMPI
!===================================================================================================================================
CONTAINS

SUBROUTINE BoundsOfElement(ElemID,Bounds)
!===================================================================================================================================
! computes the min/max of element in xyz (Based on BGMIndexOfElement)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Particle_Mesh_Vars     ,ONLY: PartElemToSide
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: Bounds(1:2,1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ilocSide, SideID
REAL                      :: xmin,xmax,ymin,ymax,zmin,zmax
REAL                      :: BezierControlPoints3D_tmp(1:3,0:NGeo,0:NGeo)
!===================================================================================================================================

xmin = HUGE(1.0)
xmax =-HUGE(1.0)
ymin = HUGE(1.0)
ymax =-HUGE(1.0)
zmin = HUGE(1.0)
zmax =-HUGE(1.0)

! get min,max of BezierControlPoints of Element
DO iLocSide = 1,6
  SideID = PartElemToSide(E2S_SIDE_ID, ilocSide, ElemID)
  IF(DoRefMapping)THEN
    IF(SideID.GT.0)THEN
      IF(PartElemToSide(E2S_FLIP,ilocSide,ElemID).EQ.0)THEN
        BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
      ELSE
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(XI_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_MINUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
        CASE(ZETA_PLUS)
          CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
        END SELECT
      END IF
    ELSE
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(XI_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0,ElemID),BezierControlPoints3D_tmp)
      CASE(ZETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo,ElemID),BezierControlPoints3D_tmp)
      END SELECT
    END IF
  ELSE ! pure tracing
    BezierControlPoints3d_tmp=BezierControlPoints3D(:,:,:,SideID)
  END IF
  xmin=MIN(xmin,MINVAL(BezierControlPoints3D_tmp(1,:,:)))
  xmax=MAX(xmax,MAXVAL(BezierControlPoints3D_tmp(1,:,:)))
  ymin=MIN(ymin,MINVAL(BezierControlPoints3D_tmp(2,:,:)))
  ymax=MAX(ymax,MAXVAL(BezierControlPoints3D_tmp(2,:,:)))
  zmin=MIN(zmin,MINVAL(BezierControlPoints3D_tmp(3,:,:)))
  zmax=MAX(zmax,MAXVAL(BezierControlPoints3D_tmp(3,:,:)))
END DO ! ilocSide
Bounds(:,1)=(/xmin,xmax/)
Bounds(:,2)=(/ymin,ymax/)
Bounds(:,3)=(/zmin,zmax/)

END SUBROUTINE BoundsOfElement


!PURE SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
!===================================================================================================================================
!> Checks if particle is inside of a linear element with triangulated faces, compatible with mortars
!> Regular element: The determinant of a 3x3 matrix, where the three vectors point from the particle to the nodes of a triangle, is
!>                  is used to determine whether the particle is inside the element. The geometric equivalent is the triple product
!>                  A*(B x C), spanning a signed volume. If the volume/determinant is positive, then the particle is inside.
!> Element with neighbouring mortar elements: Additional checks of the smaller sides are required if the particle is in not in the
!>                                       concave part of the element but in the convex. Analogous procedure using the determinants.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO,PartElemToSide,PartElemToElemAndSide
USE MOD_Mesh_Vars             ,ONLY: MortarType
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   ,INTENT(OUT)           :: Det(6,2)
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum, SideID, SideIDMortar, ind, NbElemID, nNbMortars
LOGICAL                       :: PosCheck, NegCheck, InElementCheckMortar, InElementCheckMortarNb
REAL                          :: A(1:3,1:4), crossP(3)
!===================================================================================================================================
print *,PartStateLoc
DO NodeNum = 1,8
  print *,NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+NodeNum),ElemID
END DO

InElementCheck = .TRUE.
InElementCheckMortar = .TRUE.
!--- Loop over the 6 sides of the element
DO iLocSide = 1,6
  DO NodeNum = 1,4
    !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,iLocSide,ElemID)+1) - PartStateLoc(1:3)
  END DO
  IF (ElemID.EQ.1 .AND. iLocSide.EQ.1) print *, NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,ElemID)+1)
  IF (ElemID.EQ.1 .AND. iLocSide.EQ.1) print *, NodeCoords_Shared(:,ElemSideNodeID_Shared(2,iLocSide,ElemID)+1)
  IF (ElemID.EQ.1 .AND. iLocSide.EQ.1) print *, NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,ElemID)+1)
  IF (ElemID.EQ.1 .AND. iLocSide.EQ.1) print *, NodeCoords_Shared(:,ElemSideNodeID_Shared(4,iLocSide,ElemID)+1)
  IF (ElemID.EQ.1 .AND. iLocSide.EQ.1) STOP

  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
  SideIDMortar = MortarMapping_Shared(SideID)
  IF(ElemID.EQ.1) print *, SideIDMortar
  !--- Treatment of sides which are adjacent to mortar elements
  IF (SideIDMortar.GT.0) THEN
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- Checking the concave part of the side
    IF (ConcaveElemSide_Shared(iLocSide,ElemID)) THEN
      ! If the element is actually concave, CalcDetOfTrias determines its determinants
      Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar = .FALSE.
    ELSE
      ! If its a convex element, CalcDetOfTrias determines the concave determinants
      Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar= .FALSE.
    END IF
    !--- Checking the convex part of the side
    IF (.NOT.InElementCheckMortar) THEN
      InElementCheckMortar = .TRUE.
      IF (ConcaveElemSide_Shared(iLocSide,ElemID)) THEN
        Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar = .FALSE.
      ELSE
        Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar= .FALSE.
      END IF
      !--- Particle is in a convex elem but not in concave, checking additionally the mortar neighbors. If particle is not inside
      !    the mortar elements, it has to be in the original element.
      IF (InElementCheckMortar) THEN
        nNbMortars = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)
        DO ind = 1, nNbMortars
          InElementCheckMortarNb = .TRUE.
          NbElemID = SideInfo_Shared(SIDE_ELEMID,MortarInfo_Shared(ind,MortarMapping_Shared(SideID)))
          ! If small mortar element not defined, skip it for now, likely not inside the halo region (additional check is performed
          ! after the MPI communication: ParticleInsideQuad3D_MortarMPI)
          IF (NbElemID.LT.1) CYCLE
          CALL ParticleInsideNbMortar(PartStateLoc,NbElemID,InElementCheckMortarNb)
          IF (InElementCheckMortarNb) THEN
            InElementCheck = .FALSE.
            EXIT
          END IF
        END DO
      ELSE
        InElementCheck = .FALSE.
      END IF
    END IF
  ELSE ! Treatment of regular elements without mortars
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- compute cross product for vector 1 and 3
    crossP(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
    crossP(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
    crossP(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    !--- negative determinant of triangle 1 (points 1,3,2):
    Det(iLocSide,1) = crossP(1) * A(1,2) + &
                      crossP(2) * A(2,2) + &
                      crossP(3) * A(3,2)
    Det(iLocSide,1) = -det(iLocSide,1)
    !--- determinant of triangle 2 (points 1,3,4):
    Det(iLocSide,2) = crossP(1) * A(1,4) + &
                      crossP(2) * A(2,4) + &
                      crossP(3) * A(3,4)
    IF (Det(iLocSide,1).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    IF (Det(iLocSide,2).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    !--- final determination whether particle is in element
    IF (ElemID.EQ.1) print *, InElementCheck,posCheck,NegCheck,ConcaveElemSide_Shared(iLocSide,ElemID)
    IF (ConcaveElemSide_Shared(iLocSide,ElemID)) THEN
      IF (.NOT.PosCheck) InElementCheck = .FALSE.
    ELSE
      IF (NegCheck) InElementCheck = .FALSE.
    END IF
  END IF ! Mortar element or regular element

!  SideID =PartElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
!  SideIDMortar=MortarType(2,SideID)
!  !--- Treatment of sides which are adjacent to mortar elements
!  IF (SideIDMortar.GT.0) THEN
!    PosCheck = .FALSE.
!    NegCheck = .FALSE.
!    !--- Checking the concave part of the side
!    IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
!      ! If the element is actually concave, CalcDetOfTrias determines its determinants
!      Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
!      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
!      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
!      !--- final determination whether particle is in element
!      IF (.NOT.PosCheck) InElementCheckMortar = .FALSE.
!    ELSE
!      ! If its a convex element, CalcDetOfTrias determines the concave determinants
!      Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
!      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
!      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
!      !--- final determination whether particle is in element
!      IF (.NOT.PosCheck) InElementCheckMortar= .FALSE.
!    END IF
!    !--- Checking the convex part of the side
!    IF (.NOT.InElementCheckMortar) THEN
!      InElementCheckMortar = .TRUE.
!      IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
!        Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
!        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
!        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
!        !--- final determination whether particle is in element
!        IF (NegCheck) InElementCheckMortar = .FALSE.
!      ELSE
!        Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
!        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
!        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
!        !--- final determination whether particle is in element
!        IF (NegCheck) InElementCheckMortar= .FALSE.
!      END IF
!      !--- Particle is in a convex elem but not in concave, checking additionally the mortar neighbors. If particle is not inside
!      !    the mortar elements, it has to be in the original element.
!      IF (InElementCheckMortar) THEN
!        IF (MortarType(1,SideID).EQ.1) THEN
!          nNbMortars = 4
!        ELSE
!          nNbMortars = 2
!        END IF
!        !--- Loop over the number of neighbouring mortar elements, leave the routine if the particle is found within one of the
!        !    mortar elements
!        DO ind = 1, nNbMortars
!          InElementCheckMortarNb = .TRUE.
!          NbElemID = PartElemToElemAndSide(ind,iLocSide,ElemID)
!          ! If small mortar element not defined, skip it for now, likely not inside the halo region (additional check is performed
!          ! after the MPI communication: ParticleInsideQuad3D_MortarMPI)
!          IF (NbElemID.LT.1) CYCLE
!          CALL ParticleInsideNbMortar(PartStateLoc,NbElemID,InElementCheckMortarNb)
!          IF (InElementCheckMortarNb) THEN
!            InElementCheck = .FALSE.
!            EXIT
!          END IF
!        END DO
!      ELSE
!        InElementCheck = .FALSE.
!      END IF
!    END IF
!  ELSE ! Treatment of regular elements without mortars
!    PosCheck = .FALSE.
!    NegCheck = .FALSE.
!    !--- compute cross product for vector 1 and 3
!    crossP(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
!    crossP(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
!    crossP(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
!    !--- negative determinant of triangle 1 (points 1,3,2):
!    Det(iLocSide,1) = crossP(1) * A(1,2) + &
!                      crossP(2) * A(2,2) + &
!                      crossP(3) * A(3,2)
!    Det(iLocSide,1) = -det(iLocSide,1)
!    !--- determinant of triangle 2 (points 1,3,4):
!    Det(iLocSide,2) = crossP(1) * A(1,4) + &
!                      crossP(2) * A(2,4) + &
!                      crossP(3) * A(3,4)
!    IF (Det(iLocSide,1).LT.0) THEN
!      NegCheck = .TRUE.
!    ELSE
!      PosCheck = .TRUE.
!    END IF
!    IF (Det(iLocSide,2).LT.0) THEN
!      NegCheck = .TRUE.
!    ELSE
!      PosCheck = .TRUE.
!    END IF
!    !--- final determination whether particle is in element
!    IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
!      IF (.NOT.PosCheck) InElementCheck = .FALSE.
!    ELSE
!      IF (NegCheck) InElementCheck = .FALSE.
!    END IF
!  END IF ! Mortar element or regular element
END DO ! iLocSide = 1,6

RETURN

END SUBROUTINE ParticleInsideQuad3D


PURE SUBROUTINE ParticleInsideNbMortar(PartStateLoc,ElemID,InElementCheck)
!===================================================================================================================================
!> Routines checks if the particle is inside the neighbouring mortar element. Used for the regular ParticleInsideQuad3D routine
!> after it was determined that the particle is not in the concave part but in the convex part of the element.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO, PartElemToSide
USE MOD_Mesh_Vars             ,ONLY: MortarType
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum, SideID
LOGICAL                       :: PosCheck, NegCheck
REAL                          :: A(1:3,1:4), cross(3)
REAL                          :: Det(2)
!===================================================================================================================================
InElementCheck = .TRUE.
DO iLocSide = 1,6                 ! for all 6 sides of the element
  DO NodeNum = 1,4
  !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,iLocSide,ElemID)+1) - PartStateLoc(1:3)
  END DO
  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
  IF (MortarMapping_Shared(SideID).GT.0) THEN ! Mortar side
    !--- initialize flags for side checks
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- Check if the particle is inside the convex element. If its outside, it has to be inside the original element
    IF (ConcaveElemSide_Shared(iLocSide,ElemID)) THEN
      Det(1:2) = CalcDetOfTrias(A,2)
      IF (Det(1).LT.0) NegCheck = .TRUE.
      IF (Det(2).LT.0) NegCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    ELSE
      Det(1:2) = CalcDetOfTrias(A,1)
      IF (Det(1).LT.0) NegCheck = .TRUE.
      IF (Det(2).LT.0) NegCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    END IF
  ELSE ! Regular side
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- compute cross product for vector 1 and 3
    cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
    cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
    cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    !--- negative determinant of triangle 1 (points 1,3,2):
    Det(1) = cross(1) * A(1,2) + &
                      cross(2) * A(2,2) + &
                      cross(3) * A(3,2)
    Det(1) = -det(1)
    !--- determinant of triangle 2 (points 1,3,4):
    Det(2) = cross(1) * A(1,4) + &
                      cross(2) * A(2,4) + &
                      cross(3) * A(3,4)
    IF (Det(1).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    IF (Det(2).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    !--- final determination whether particle is in element
    IF (ConcaveElemSide_Shared(iLocSide,ElemID)) THEN
      IF (.NOT.PosCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    ELSE
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    END IF
  END IF  ! Mortar or regular side
END DO  ! iLocSide = 1,6

END SUBROUTINE ParticleInsideNbMortar


PURE FUNCTION CalcDetOfTrias(A,bending)
!================================================================================================================================
!> Calculates the determinant A*(B x C) for both triangles of a side. bending = 1 gives the determinant considering the actual
!> orientation of the side (concave/convex), 2 gives the opposite of the saved form (e.g. a concave side gets the convex analog)
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                     :: A(3,4)
INTEGER,INTENT(IN)                  :: bending
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                :: CalcDetOfTrias(2)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: cross(3)
!================================================================================================================================
IF (bending.EQ.1) THEN
  !--- compute cross product for vector 1 and 3
  cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
  cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
  cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

  !--- negative determinant of triangle 1 (points 1,3,2):
  CalcDetOfTrias(1) = cross(1) * A(1,2) + &
                   cross(2) * A(2,2) + &
                   cross(3) * A(3,2)
  CalcDetOfTrias(1)  = -CalcDetOfTrias(1)
  !--- determinant of triangle 2 (points 1,3,4):
  CalcDetOfTrias(2)  = cross(1) * A(1,4) + &
                   cross(2) * A(2,4) + &
                   cross(3) * A(3,4)
ELSE
  !--- compute cross product for vector 2 and 4
  cross(1) = A(2,2) * A(3,4) - A(3,2) * A(2,4)
  cross(2) = A(3,2) * A(1,4) - A(1,2) * A(3,4)
  cross(3) = A(1,2) * A(2,4) - A(2,2) * A(1,4)

  !--- negative determinant of triangle 1 (points 2,4,1):
  CalcDetOfTrias(1) = cross(1) * A(1,1) + &
                   cross(2) * A(2,1) + &
                   cross(3) * A(3,1)
  !--- determinant of triangle 2 (points 2,4,3):
  CalcDetOfTrias(2) = cross(1) * A(1,3) + &
                   cross(2) * A(2,3) + &
                   cross(3) * A(3,3)
  CalcDetOfTrias(2) = -CalcDetOfTrias(2)
END IF

END FUNCTION CalcDetOfTrias


SUBROUTINE ParticleInsideQuad3D_MortarMPI(PartStateLoc,ElemID,PartInside)
!===================================================================================================================================
!> Analogous to the original routine, checks particles which were communicated if they are in the correct element.
!> Checks if particle is inside of a linear element with triangulated faces, compatible with mortars
!> Regular element: The determinant of a 3x3 matrix, where the three vectors point from the particle to the nodes of a triangle, is
!>                  is used to determine whether the particle is inside the element. The geometric equivalent is the triple product
!>                  A*(B x C), spanning a signed volume. If the volume/determinant is positive, then the particle is inside.
!> Element with neighbouring mortar elements: Additional checks of the smaller sides are required if the particle is in not in the
!>                                       concave part of the element but in the convex. Analogous procedure using the determinants.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO,PartElemToSide,PartElemToElemAndSide
USE MOD_Mesh_Vars             ,ONLY: MortarType
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL   ,INTENT(IN)            :: PartStateLoc(3)
INTEGER,INTENT(INOUT)         :: ElemID
LOGICAL,INTENT(INOUT)         :: PartInside
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum, SideID, SideIDMortar, ind, NbElemID, nNbMortars
LOGICAL                       :: PosCheck, NegCheck, InElementCheckMortar, InElementCheckMortarNb
REAL                          :: A(1:3,1:4), Det(2)
!===================================================================================================================================
InElementCheckMortar = .TRUE.
!--- Loop over the 6 sides of the element
DO iLocSide = 1,6
  DO NodeNum = 1,4
    !--- A = vector from particle to node coords
    A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,ElemID)) - PartStateLoc(1:3)
  END DO
  SideID =PartElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
  SideIDMortar=MortarType(2,SideID)
  IF (SideIDMortar.LE.0) CYCLE
  PosCheck = .FALSE.
  NegCheck = .FALSE.
  !--- Checking the concave part of the side
  IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
    ! If the element is actually concave, CalcDetOfTrias determines its determinants
    Det(1:2) = CalcDetOfTrias(A,1)
    IF (Det(1).GE.0) PosCheck = .TRUE.
    IF (Det(2).GE.0) PosCheck = .TRUE.
    !--- final determination whether particle is in element
    IF (.NOT.PosCheck) InElementCheckMortar = .FALSE.
  ELSE
    ! If its a convex element, CalcDetOfTrias determines the concave determinants
    Det(1:2) = CalcDetOfTrias(A,2)
    IF (Det(1).GE.0) PosCheck = .TRUE.
    IF (Det(2).GE.0) PosCheck = .TRUE.
    !--- final determination whether particle is in element
    IF (.NOT.PosCheck) InElementCheckMortar= .FALSE.
  END IF
  IF (InElementCheckMortar) THEN
    IF (MortarType(1,SideID).EQ.1) THEN
      nNbMortars = 4
    ELSE
      nNbMortars = 2
    END IF
    !--- Loop over the number of neighbouring mortar elements, leave the routine if the particle is found within one of the
    !    mortar elements
    DO ind = 1, nNbMortars
      InElementCheckMortarNb = .TRUE.
      NbElemID = PartElemToElemAndSide(ind,iLocSide,ElemID)
      IF (NbElemID.LT.1) THEN
        IPWRITE(*,*) 'PartState:', PartStateLoc(1:3)
        IPWRITE(*,*) 'ElemID:', ElemID
        IPWRITE(*,*) 'WARNING: Particle deleted! Think about increasing the halo region (HaloEpsVelo)!'
        PartInside = .FALSE.
      END IF
      CALL ParticleInsideNbMortar(PartStateLoc,NbElemID,InElementCheckMortarNb)
      IF (InElementCheckMortarNb) THEN
        ElemID = NbElemID
        RETURN
      END IF
    END DO
  END IF
END DO ! iLocSide = 1,6

END SUBROUTINE ParticleInsideQuad3D_MortarMPI


!==================================================================================================================================!
!> Initialize GetGlobalElemID function (mapping of compute-node element ID to global element ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalElemID => GetGlobalElemID_iElem
ELSE
  GetGlobalElemID => GetGlobalElemID_fromTotalElem
END IF
#else
GetGlobalElemID => GetGlobalElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElemID_fromTotalElem(1)
#endif
dummy=GetGlobalElemID_iElem(1)
END SUBROUTINE InitGetGlobalElemID


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalElemID_iElem
!===================================================================================================================================
GetGlobalElemID_iElem = iElem
END FUNCTION GetGlobalElemID_iElem


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalElemID_fromTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:CNTotalElem2GlobalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElemID_fromTotalElem = CNTotalElem2GlobalElem(iElem)
END FUNCTION GetGlobalElemID_fromTotalElem
#endif /*USE_MPI*/


!==================================================================================================================================!
!> Initialize GetCNElemID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNElemID => GetCNElemID_iElem
ELSE
  GetCNElemID => GetCNElemID_fromTotalElem
END IF
#else
GetCNElemID => GetCNElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetCNElemID_fromTotalElem(1)
#endif
dummy=GetCNElemID_iElem(1)
END SUBROUTINE InitGetCNElemID


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION GetCNElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetCNElemID_iElem
!===================================================================================================================================
GetCNElemID_iElem = iElem
END FUNCTION GetCNElemID_iElem


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE FUNCTION GetCNElemID_fromTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:GlobalElem2CNTotalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetCNElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetCNElemID_fromTotalElem = GlobalElem2CNTotalElem(iElem)
END FUNCTION GetCNElemID_fromTotalElem
#endif /*USE_MPI*/


END MODULE MOD_Particle_Mesh_Tools
