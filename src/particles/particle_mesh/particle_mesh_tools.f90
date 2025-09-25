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

MODULE MOD_Particle_Mesh_Tools
!===================================================================================================================================
! Contains subroutines to build (curvilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Initialization routines
INTERFACE GetGlobalNonUniqueSideID
  MODULE PROCEDURE GetGlobalNonUniqueSideID
END INTERFACE

INTERFACE GetSideBoundingBoxTria
  MODULE PROCEDURE GetSideBoundingBoxTria
END INTERFACE

! Define an interface for the function pointer
ABSTRACT INTERFACE
  SUBROUTINE ParticleInsideQuadInterface(PartStateLoc,GlobalElemID,InElementCheck,Det_Out)
    INTEGER,INTENT(IN)            :: GlobalElemID
    REAL   ,INTENT(IN)            :: PartStateLoc(3)
    LOGICAL,INTENT(OUT)           :: InElementCheck
    REAL   ,INTENT(OUT),OPTIONAL  :: Det_Out(6,2)
  END SUBROUTINE
END INTERFACE

!> Pointer defining the localization routine based on the symmetry order
PROCEDURE(ParticleInsideQuadInterface),POINTER :: ParticleInsideQuad => NULL()

PUBLIC :: ParticleInsideQuad3D, InitPEM_LocalElemID, InitPEM_CNElemID, GetGlobalNonUniqueSideID, GetSideBoundingBoxTria
PUBLIC :: GetMeshMinMax, IdentifyElemAndSideType, WeirdElementCheck, CalcParticleMeshMetrics, CalcXCL_NGeo
PUBLIC :: CalcBezierControlPoints, InitParticleGeometry, ComputePeriodicVec
PUBLIC :: InitParticleInsideQuad, ParticleInsideQuad
PUBLIC :: InitVolumes_2D, InitVolumes_1D, DSMC_2D_CalcSymmetryArea, DSMC_1D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
!===================================================================================================================================
CONTAINS

!==================================================================================================================================!
!> Initialize ParticleInsideQuad depending on symmetry dimension using a function pointer
!==================================================================================================================================!
SUBROUTINE InitParticleInsideQuad()
! MODULES
USE MOD_Globals
USE MOD_Symmetry_Vars            ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================

IF (Symmetry%Order.EQ.3) THEN
  ParticleInsideQuad => ParticleInsideQuad3D
ELSE IF (Symmetry%Order.EQ.2) THEN
  ParticleInsideQuad => ParticleInsideQuad2D
ELSE IF (Symmetry%Order.EQ.1) THEN
  ParticleInsideQuad => ParticleInsideQuad1D
ELSE
  CALL abort(__STAMP__,'ERROR in InitParticleInsideQuad: Function pointer could not be properly defined!')
END IF

END SUBROUTINE InitParticleInsideQuad

!PPURE SUBROUTINE ParticleInsideQuad3D(PartStateLoc,GlobalElemID,InElementCheck,Det)
SUBROUTINE ParticleInsideQuad3D(PartStateLoc,GlobalElemID,InElementCheck,Det_Out)
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
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars    ,ONLY :ConcaveElemSide_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: GlobalElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
REAL   ,INTENT(OUT),OPTIONAL  :: Det_Out(6,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum, SideID, SideIDMortar, ind, NbElemID, nNbMortars, nlocSides, localSideID
INTEGER                       :: CNElemID
LOGICAL                       :: PosCheck, NegCheck, InElementCheckMortar, InElementCheckMortarNb
REAL                          :: A(1:3,1:4), crossP(3), Det(6,2)
!===================================================================================================================================
InElementCheck = .TRUE.
InElementCheckMortar = .TRUE.
!--- Loop over the 6 sides of the element
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
CNElemID = GetCNElemID(GlobalElemID)
DO iLocSide = 1,nlocSides
  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
  localSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
  IF (localSideID.LE.0) CYCLE
  DO NodeNum = 1,4
    !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,CNElemID)+1) - PartStateLoc(1:3)
  END DO

  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
  !--- Treatment of sides which are adjacent to mortar elements
  IF (NbElemID.LT.0) THEN
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- Checking the concave part of the side
    IF (ConcaveElemSide_Shared(localSideID,CNElemID)) THEN
      ! If the element is actually concave, CalcDetOfTrias determines its determinants
      Det(localSideID,1:2) = CalcDetOfTrias(A,1)
      IF (Det(localSideID,1).GE.0) PosCheck = .TRUE.
      IF (Det(localSideID,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar = .FALSE.
    ELSE
      ! If its a convex element, CalcDetOfTrias determines the concave determinants
      Det(localSideID,1:2) = CalcDetOfTrias(A,2)
      IF (Det(localSideID,1).GE.0) PosCheck = .TRUE.
      IF (Det(localSideID,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar= .FALSE.
    END IF
    !--- Checking the convex part of the side
    IF (.NOT.InElementCheckMortar) THEN
      InElementCheckMortar = .TRUE.
      IF (ConcaveElemSide_Shared(localSideID,CNElemID)) THEN
        Det(localSideID,1:2) = CalcDetOfTrias(A,2)
        IF (Det(localSideID,1).LT.0) NegCheck = .TRUE.
        IF (Det(localSideID,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar = .FALSE.
      ELSE
        Det(localSideID,1:2) = CalcDetOfTrias(A,1)
        IF (Det(localSideID,1).LT.0) NegCheck = .TRUE.
        IF (Det(localSideID,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar= .FALSE.
      END IF
      !--- Particle is in a convex elem but not in concave, checking additionally the mortar neighbors. If particle is not inside
      !    the mortar elements, it has to be in the original element.
      IF (InElementCheckMortar) THEN
        nNbMortars = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)
        DO ind = 1, nNbMortars
          InElementCheckMortarNb = .TRUE.
          SideIDMortar = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide + ind
          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideIDMortar)
          ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
          ! no way to recover it during runtime
          IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',GlobalElemID)

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
    Det(localSideID,1) = crossP(1) * A(1,2) + &
                      crossP(2) * A(2,2) + &
                      crossP(3) * A(3,2)
    Det(localSideID,1) = -det(localSideID,1)
    !--- determinant of triangle 2 (points 1,3,4):
    Det(localSideID,2) = crossP(1) * A(1,4) + &
                      crossP(2) * A(2,4) + &
                      crossP(3) * A(3,4)
    IF (Det(localSideID,1).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    IF (Det(localSideID,2).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF
    !--- final determination whether particle is in element
    IF (ConcaveElemSide_Shared(localSideID,CNElemID)) THEN
      IF (.NOT.PosCheck) InElementCheck = .FALSE.
    ELSE
      IF (NegCheck) InElementCheck = .FALSE.
    END IF
  END IF ! Mortar element or regular element
END DO ! iLocSide = 1,6

IF(PRESENT(Det_Out)) Det_Out = Det

RETURN

END SUBROUTINE ParticleInsideQuad3D


SUBROUTINE ParticleInsideQuad2D(PartStateLoc,GlobalElemID,InElementCheck,Det_Out)
!===================================================================================================================================
!> Checks if particle is inside of a linear 2D element with 4  faces, compatible with mortars. The "Ray Casting Algorithm" is used.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared, SideIsSymSide
USE MOD_Particle_Mesh_Vars    ,ONLY :ElemSideNodeID2D_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: GlobalElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
REAL   ,INTENT(OUT),OPTIONAL  :: Det_Out(6,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, TempSideID, nlocSides, localSideID
INTEGER                       :: CNElemID
REAL                          :: x_int, xNode1, xNode2, yNode1, yNode2
!===================================================================================================================================
CNElemID = GetCNElemID(GlobalElemID)
InElementCheck = .FALSE.
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)

DO iLocSide=1,nlocSides
  TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
  localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
  ! Side is not one of the 6 local sides
  IF (localSideID.LE.0) CYCLE
  IF (SideIsSymSide(TempSideID)) CYCLE
  xNode1 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(1,localSideID, CNElemID))
  yNode1 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(1,localSideID, CNElemID))
  xNode2 = NodeCoords_Shared(1,ElemSideNodeID2D_Shared(2,localSideID, CNElemID))
  yNode2 = NodeCoords_Shared(2,ElemSideNodeID2D_Shared(2,localSideID, CNElemID))
  IF ( (yNode1 >= PartStateLoc(2) .AND. yNode2 < PartStateLoc(2)) .OR. &
       (yNode1 < PartStateLoc(2) .AND. yNode2 >= PartStateLoc(2)) ) THEN
    ! Compute x-coordinate of the intersection point
    x_int = (PartStateLoc(2)- yNode1) * (xNode2 - xNode1) / (yNode2 - yNode1) + xNode1
    ! Check if the ray crosses the edge to the right
    IF (x_int > PartStateLoc(1)) THEN
      InElementCheck = .NOT.InElementCheck
    END IF
  END IF
END DO

! Surpress compiler warning
RETURN
IF(PRESENT(Det_Out)) Det_Out = 0.0

END SUBROUTINE ParticleInsideQuad2D

SUBROUTINE ParticleInsideQuad1D(PartStateLoc,GlobalElemID,InElementCheck,Det_Out)
!===================================================================================================================================
!> Checks if particle is inside of a 1D element.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared, SideIsSymSide
USE MOD_Particle_Mesh_Vars    ,ONLY :ElemSideNodeID1D_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: GlobalElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
REAL   ,INTENT(OUT),OPTIONAL  :: Det_Out(6,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, TempSideID, nlocSides, localSideID, DiffSign(2), iSide
INTEGER                       :: CNElemID
REAL                          :: xNode
!===================================================================================================================================
CNElemID = GetCNElemID(GlobalElemID)
InElementCheck = .FALSE.
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
iSide = 0
DO iLocSide=1,nlocSides
  TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
  localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
  ! Side is not one of the 6 local sides
  IF (localSideID.LE.0) CYCLE
  IF (SideIsSymSide(TempSideID)) CYCLE
  xNode = NodeCoords_Shared(1,ElemSideNodeID1D_Shared(localSideID, CNElemID))
  iSide = iSide + 1
  DiffSign(iSide) = NINT(SIGN(1.,PartStateLoc(1) - xNode))
  IF (iSide.EQ.2) EXIT
END DO
IF (DiffSign(1).NE.DiffSign(2)) InElementCheck = .TRUE.
! Surpress compiler warning
RETURN
IF(PRESENT(Det_Out)) Det_Out = 0.0
END SUBROUTINE ParticleInsideQuad1D


PPURE SUBROUTINE ParticleInsideNbMortar(PartStateLoc,GlobalElemID,InElementCheck)
!===================================================================================================================================
!> Routines checks if the particle is inside the neighbouring mortar element. Used for the regular ParticleInsideQuad3D routine
!> after it was determined that the particle is not in the concave part but in the convex part of the element.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars    ,ONLY :ConcaveElemSide_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: GlobalElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide, NodeNum, SideID, nlocSides, localSideID, NbElemID, CNElemID
LOGICAL                       :: PosCheck, NegCheck
REAL                          :: A(1:3,1:4), cross(3)
REAL                          :: Det(2)
!===================================================================================================================================
InElementCheck = .TRUE.
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
CNElemID = GetCNElemID(GlobalElemID)
DO iLocSide = 1,nlocSides
  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
  localSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
  IF (localSideID.LE.0) CYCLE
  DO NodeNum = 1,4
  !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,CNElemID)+1) - PartStateLoc(1:3)
  END DO
  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
  !--- Treatment of sides which are adjacent to mortar elements
  IF (NbElemID.LT.0) THEN
    !--- initialize flags for side checks
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- Check if the particle is inside the convex element. If its outside, it has to be inside the original element
    IF (ConcaveElemSide_Shared(localSideID,CNElemID)) THEN
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
    IF (ConcaveElemSide_Shared(localSideID,CNElemID)) THEN
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


PPURE FUNCTION CalcDetOfTrias(A,bending)
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


!==================================================================================================================================!
!> Initialize PEM%LocalElemID(iPart) function (mapping of global element ID, which is first obtained from PEM%GlobalElemID(iPart) to
!> the core local element ID)
!==================================================================================================================================!
SUBROUTINE InitPEM_LocalElemID()
! MODULES
USE MOD_Particle_Vars ,ONLY: PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
PEM%LocalElemID => GetGlobalID_offset
#else
PEM%LocalElemID => GetGlobalID
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalID_offset(1)
#else
dummy=GetGlobalID(1)
#endif
END SUBROUTINE InitPEM_LocalElemID


!==================================================================================================================================!
!> Get the global element ID from PEM%GlobalElemID(iPart)
!==================================================================================================================================!
PPURE FUNCTION GetGlobalID(iPart)
! MODULES
USE MOD_Particle_Vars ,ONLY: PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iPart ! Particle ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalID
!===================================================================================================================================
GetGlobalID = PEM%GlobalElemID(iPart)
END FUNCTION GetGlobalID


#if USE_MPI
!==================================================================================================================================!
!> Get the local element ID from PEM%GlobalElemID(iPart) - offsetElem
!==================================================================================================================================!
PPURE FUNCTION GetGlobalID_offset(iPart)
! MODULES
USE MOD_Mesh_Vars     ,ONLY: offSetElem
USE MOD_Particle_Vars ,ONLY: PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iPart ! Particle ID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalID_offset
!===================================================================================================================================
GetGlobalID_offset = PEM%GlobalElemID(iPart) - offsetElem
END FUNCTION GetGlobalID_offset
#endif /*USE_MPI*/


!==================================================================================================================================!
!> Initialize PEM%CNElemID(iPart) function (mapping of global element ID, which is first obtained from PEM%GlobalElemID(iPart) to
!> compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitPEM_CNElemID()
! MODULES
USE MOD_Particle_Vars   ,ONLY: PEM
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  PEM%CNElemID => GetGlobalID
ELSE
  PEM%CNElemID => GetGlobalElem2CNTotalElem_iPart
END IF
#else
PEM%CNElemID => GetGlobalID
#endif /*USE_MPI*/

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElem2CNTotalElem_iPart(1)
#endif /*USE_MPI*/
dummy=GetGlobalID(1)
END SUBROUTINE InitPEM_CNElemID


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID from the global element ID, which is first obtained from PEM%GlobalElemID(iPart) in case of MPI=ON for
!> single or multiple compute nodes (CN)
!==================================================================================================================================!
PPURE FUNCTION GetGlobalElem2CNTotalElem_iPart(iPart)
! MODULES
USE MOD_Particle_Mesh_Vars ,ONLY: GlobalElem2CNTotalElem
USE MOD_Particle_Vars      ,ONLY: PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iPart ! Particle ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalElem2CNTotalElem_iPart
!===================================================================================================================================
GetGlobalElem2CNTotalElem_iPart = GlobalElem2CNTotalElem(PEM%GlobalElemID(iPart))
END FUNCTION GetGlobalElem2CNTotalElem_iPart
#endif /*USE_MPI*/


FUNCTION GetGlobalNonUniqueSideID(GlobalElemID,localSideID)
!===================================================================================================================================
!> Determines the non-unique global side ID of the local side in global element ID
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared, SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: GlobalElemID                        !< global element ID
INTEGER,INTENT(IN) :: localSideID                         !< local side id of an element (1:6)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER :: GetGlobalNonUniqueSideID
INTEGER :: iSide,firstSide,lastSide
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + 1
lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND, GlobalElemID)

! Small mortar sides are added after
DO iSide = firstSide,lastSide
  IF (SideInfo_Shared(SIDE_LOCALID,iSide).EQ.localSideID) THEN
    GetGlobalNonUniqueSideID = iSide
    RETURN
  END IF
END DO

! We should never arrive here
GetGlobalNonUniqueSideID=-1
CALL ABORT(__STAMP__,'GlobalSideID not found for Elem',GlobalElemID)
END FUNCTION GetGlobalNonUniqueSideID

!==================================================================================================================================!
!> Initialize GetGlobalElemID function (mapping of compute-node element ID to global element ID)
!==================================================================================================================================!
SUBROUTINE GetSideBoundingBoxTria(SideID, BoundingBox)
! MODULES
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: SideID
REAL, INTENT(OUT)             :: BoundingBox(1:3,1:8)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iLocSide, CNElemID, iNode
REAL                      :: NodePoints(1:3,1:4)
REAL                      :: xMin, xMax, yMin, yMax, zMin, zMax
!==================================================================================================================================
CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
DO iNode = 1, 4
  NodePoints(1:3,iNode) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)
END DO
xMin = MINVAL(NodePoints(1,:))
yMin = MINVAL(NodePoints(2,:))
zMin = MINVAL(NodePoints(3,:))
xMax = MAXVAL(NodePoints(1,:))
yMax = MAXVAL(NodePoints(2,:))
zMax = MAXVAL(NodePoints(3,:))
BoundingBox(1:3,1) = (/xMin,yMin,zMin/)
BoundingBox(1:3,2) = (/xMax,yMin,zMin/)
BoundingBox(1:3,3) = (/xMax,yMax,zMin/)
BoundingBox(1:3,4) = (/xMin,yMax,zMin/)
BoundingBox(1:3,5) = (/xMin,yMin,zMax/)
BoundingBox(1:3,6) = (/xMax,yMin,zMax/)
BoundingBox(1:3,7) = (/xMax,yMax,zMax/)
BoundingBox(1:3,8) = (/xMin,yMax,zMax/)

END SUBROUTINE GetSideBoundingBoxTria


SUBROUTINE GetMeshMinMax()
!===================================================================================================================================
! computes the minimum and maximum value of the mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: offsetElem,nElems
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED,myComputeNodeRank
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: offsetLocalSide,nLocalSides
INTEGER                        :: offsetLocalNode,nLocalNodes
!===================================================================================================================================

SELECT CASE(TrackingMethod)
  ! Build mesh min/max on BezierControlPoints for possibly curved elements
  CASE(REFMAPPING,TRACING)
    ! calculate all offsets
    offsetLocalSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)
    nLocalSides     = ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%xmax     = MAXVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymin     = MINVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymax     = MAXVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmin     = MINVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmax     = MAXVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))

#if USE_MPI
    ! compute-node local
    CALL MPI_ALLREDUCE(GEO%xmin,GEO%CNxmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%xmax,GEO%CNxmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%ymin,GEO%CNymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%ymax,GEO%CNymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%zmin,GEO%CNzmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%zmax,GEO%CNzmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)

    ! global
    IF (myComputeNodeRank.EQ.0) THEN
      CALL MPI_ALLREDUCE(GEO%CNxmin,GEO%xminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNxmax,GEO%xmaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNymin,GEO%yminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNymax,GEO%ymaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNzmin,GEO%zminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNzmax,GEO%zmaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
    END IF
    CALL MPI_BCAST(GEO%xminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%xmaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%yminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%ymaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%zminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%zmaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/
  ! TriaTracking does not have curved elements, nodeCoords are sufficient
  CASE(TRIATRACKING)
    ! calculate all offsets
    offsetLocalNode = ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)
    nLocalNodes     = ElemInfo_Shared(ELEM_LASTNODEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%xmax     = MAXVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymin     = MINVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymax     = MAXVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmin     = MINVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmax     = MAXVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))

#if USE_MPI
    ! compute-node local
    CALL MPI_ALLREDUCE(GEO%xmin,GEO%CNxmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%xmax,GEO%CNxmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%ymin,GEO%CNymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%ymax,GEO%CNymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%zmin,GEO%CNzmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_SHARED,iError)
    CALL MPI_ALLREDUCE(GEO%zmax,GEO%CNzmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_SHARED,iError)

    ! global
    IF (myComputeNodeRank.EQ.0) THEN
      CALL MPI_ALLREDUCE(GEO%CNxmin,GEO%xminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNxmax,GEO%xmaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNymin,GEO%yminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNymax,GEO%ymaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNzmin,GEO%zminglob,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_LEADERS_SHARED,iError)
      CALL MPI_ALLREDUCE(GEO%CNzmax,GEO%zmaxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LEADERS_SHARED,iError)
    END IF
    CALL MPI_BCAST(GEO%xminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%xmaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%yminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%ymaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%zminglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
    CALL MPI_BCAST(GEO%zmaxglob,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/
END SELECT

#if !USE_MPI
! compute-node local (dummy)
GEO%CNxmin   = GEO%xmin
GEO%CNxmax   = GEO%xmax
GEO%CNymin   = GEO%ymin
GEO%CNymax   = GEO%ymax
GEO%CNzmin   = GEO%zmin
GEO%CNzmax   = GEO%zmax

! global (dummy)
GEO%xminglob = GEO%xmin
GEO%xmaxglob = GEO%xmax
GEO%yminglob = GEO%ymin
GEO%ymaxglob = GEO%ymax
GEO%zminglob = GEO%zmin
GEO%zmaxglob = GEO%zmax
#endif /*USE_MPI*/

LBWRITE(UNIT_stdOut,'(A,E18.8,A,E18.8,A,E18.8)') ' | Total MESH   Dim (x,y,z): '                                     &
                                                , MAXVAL(NodeCoords_Shared(1,:))-MINVAL(NodeCoords_Shared(1,:)),', '&
                                                , MAXVAL(NodeCoords_Shared(2,:))-MINVAL(NodeCoords_Shared(2,:)),', '&
                                                , MAXVAL(NodeCoords_Shared(3,:))-MINVAL(NodeCoords_Shared(3,:))
IF (TrackingMethod.EQ.REFMAPPING .OR. TrackingMethod.EQ. TRACING) THEN
  LBWRITE(UNIT_stdOut,'(A,E18.8,A,E18.8,A,E18.8)') ' | Total BEZIER Dim (x,y,z): '                                   &
                                                  , GEO%xmaxglob-GEO%xminglob,', '                                  &
                                                  , GEO%ymaxglob-GEO%yminglob,', '                                  &
                                                  , GEO%zmaxglob-GEO%zminglob
END IF

END SUBROUTINE GetMeshMinMax


SUBROUTINE IdentifyElemAndSideType()
!===================================================================================================================================
!> get the element and side type of each element
!> 1) Get Elem Type (curved_elem)
!> 2) Get Side Type (planar_rect, planar_nonrect, bilineard, curved, planar_curved)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_ChangeBasis            ,ONLY: changeBasis3D
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNSideID
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides
USE MOD_Mesh_Vars              ,ONLY: Vdm_CLNGeo1_CLNGeo,NGeo,Vdm_CLNGeo1_CLNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared,ElemBaryNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared,ElemCurved
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Surfaces_Vars ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,SideType,SideNormVec,SideDistance
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems,nComputeNodeTotalSides
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideDistance_Shared,SideDistance_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideType_Shared,SideType_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideNormVec_Shared,SideNormVec_Shared_Win
#else
USE MOD_Particle_Mesh_Vars     ,ONLY: nComputeNodeElems
#endif /* USE_MPI */
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,SideID,CNSideID,flip
INTEGER                                  :: iCNElem,CNElemID,firstElem,lastElem,iGlobalElem,GlobalElemID
REAL,DIMENSION(1:3)                      :: v1,v2,v3
LOGICAL,ALLOCATABLE                      :: SideIsDone(:)
REAL                                     :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                                     :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: BezierControlPoints_loc(1:3,0:NGeo,0:NGeo)
INTEGER                                  :: NGeo3,NGeo2
REAL                                     :: XCL_NGeoSideNew(1:3,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoSideOld(1:3,0:NGeo,0:NGeo)
LOGICAL                                  :: isCurvedSide,isRectangular
! output and sanity check
INTEGER                                  :: nPlanarRectangular,   nPlanarNonRectangular,   nPlanarCurved,   nBilinear,   nCurved
INTEGER                                  :: nPlanarRectangularTot,nPlanarNonRectangularTot,nPlanarCurvedTot,nBilinearTot,nCurvedTot
INTEGER                                  :: nLinearElems,   nCurvedElems
INTEGER                                  :: nLinearElemsTot,nCurvedElemsTot
!#if USE_MPI
!INTEGER                                  :: nDummy
!#endif /* USE_MPI */
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_StdOut,'(A)') ' Identifying side types (planar, bilinear) and whether elements are curved ...'

! elements
#if USE_MPI
CALL Allocate_Shared((/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
ElemCurved => ElemCurved_Shared
#else
ALLOCATE(ElemCurved(1:nComputeNodeElems))
#endif /*USE_MPI*/

! sides
#if USE_MPI
CALL Allocate_Shared((/nComputeNodeTotalSides/),  SideType_Shared_Win,    SideType_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideType_Shared_Win,IERROR)
SideType => SideType_Shared
CALL Allocate_Shared((/nComputeNodeTotalSides/),  SideDistance_Shared_Win,SideDistance_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideDistance_Shared_Win,IERROR)
SideDistance => SideDistance_Shared
CALL Allocate_Shared((/3,nComputeNodeTotalSides/),SideNormVec_Shared_Win, SideNormVec_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideNormVec_Shared_Win,IERROR)
SideNormVec => SideNormVec_Shared
#else
ALLOCATE(SideType(       nNonUniqueGlobalSides))
ALLOCATE(SideDistance(   nNonUniqueGlobalSides))
ALLOCATE(SideNormVec(1:3,nNonUniqueGlobalSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemCurved   = .FALSE.
  SideType     = -1
  SideDistance = -0.
  SideNormVec  = 0.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(ElemCurved_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideType_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideDistance_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideNormVec_Shared_Win ,MPI_COMM_SHARED)
#endif /* USE_MPI*/

ALLOCATE(SideIsDone(nNonUniqueGlobalSides))
SideIsDone = .FALSE.

NGeo2 = (NGeo+1)*(NGeo+1)
NGeo3 = NGeo2   *(NGeo+1)

! decide if element is (bi-)linear or curved
! decide if sides are planar-rect, planar-nonrect, planar-curved, bilinear or curved
#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

DO iCNElem=firstElem,lastElem
  GlobalElemID = GetGlobalElemID(iCNElem)

  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,GlobalElemID)
  ! 1) check if elem is curved
  !   a) get the coordinates of the eight nodes of the hexahedral
  XCL_NGeo1(1:3,0,0,0) = XCL_NGeoLoc(1:3, 0  , 0  , 0  )
  XCL_NGeo1(1:3,1,0,0) = XCL_NGeoLoc(1:3,NGeo, 0  , 0  )
  XCL_NGeo1(1:3,0,1,0) = XCL_NGeoLoc(1:3, 0  ,NGeo, 0  )
  XCL_NGeo1(1:3,1,1,0) = XCL_NGeoLoc(1:3,NGeo,NGeo, 0  )
  XCL_NGeo1(1:3,0,0,1) = XCL_NGeoLoc(1:3, 0  , 0  ,NGeo)
  XCL_NGeo1(1:3,1,0,1) = XCL_NGeoLoc(1:3,NGeo, 0  ,NGeo)
  XCL_NGeo1(1:3,0,1,1) = XCL_NGeoLoc(1:3, 0  ,NGeo,NGeo)

  XCL_NGeo1(1:3,1,1,1) = XCL_NGeoLoc(1:3,NGeo,NGeo,NGeo)
  !  b) interpolate from the nodes to NGeo
  !     Compare the bi-linear mapping with the used mapping
  !     For NGeo=1, this should always be true, because the mappings are identical
  CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
  ! check the coordinates of all Chebychev-Lobatto geometry points between the bi-linear and used mapping
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo),ElemCurved(iCNElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID   = GetGlobalNonUniqueSideID(GlobalElemID,iLocSide)
    CNSideID = GetCNSideID(SideID)
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    IF(.NOT.ElemCurved(iCNElem))THEN
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! linear element
      IF(BoundingBoxIsEmpty(CNSideID))THEN
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
        v1=0.25*(BezierControlPoints_loc(:,0,0      )     &
                +BezierControlPoints_loc(:,NGeo,0   )  &
                +BezierControlPoints_loc(:,0,NGeo   )  &
                +BezierControlPoints_loc(:,NGeo,NGeo))
        ! check if normal vector points outwards
        v2=v1-ElemBaryNGeo(:,iCNElem)
        IF(flip.EQ.0)THEN
          IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
        ELSE
          IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
        END IF
        SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
        ! check if it is rectangular
        isRectangular=.TRUE.
        v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
        v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
        v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        IF(isRectangular)THEN
          v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        END IF
        IF(isRectangular)THEN
          SideType(CNSideID)=PLANAR_RECT
        ELSE
          SideType(CNSideID)=PLANAR_NONRECT
        END IF
      ELSE
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,CNSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
        SideType(CNSideID)=BILINEAR
      END IF
    ELSE
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! possible curved face
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0,0:NGeo,0:NGeo)
      CASE(XI_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,NGeo,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,NGeo,0:NGeo,0:NGeo)
      CASE(ETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0,0:NGeo)
      CASE(ETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,NGeo,0:NGeo)
      CASE(ZETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0)
      CASE(ZETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,NGeo)
        XCL_NGeoSideNew=XCL_NGeoNEw(1:3,0:NGeo,0:NGeo,NGeo)
      END SELECT
      CALL PointsEqual(NGeo2,XCL_NGeoSideNew,XCL_NGeoSideOld,isCurvedSide)
      IF(isCurvedSide)THEn
        IF(BoundingBoxIsEmpty(CNSideID))THEN
          SideType(CNSideID)=PLANAR_CURVED
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0   )  &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iCNElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          END IF
          SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
        ELSE
          SideType(CNSideID)=CURVED
        END IF
      ELSE
        IF (BoundingBoxIsEmpty(CNSideID)) THEN
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0)     &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iCNElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          END IF
          SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
          ! check if it is rectangular
          isRectangular=.TRUE.
          v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
          v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
          v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          IF(isRectangular)THEN
            v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          END IF
          IF(isRectangular)THEN
            SideType(CNSideID)=PLANAR_RECT
          ELSE
            SideType(CNSideID)=PLANAR_NONRECT
          END IF
        ELSE
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
          SideType(CNSideID)=BILINEAR
        END IF
      END IF
    END IF
    SideIsDone(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iCNElem = firstElem, lastElem (nComputeNodeTotalElems)

DEALLOCATE(SideIsDone)

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemCurved_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI */

! zero counter for side types
nPlanarRectangular         = 0
nPlanarNonRectangular      = 0
nPlanarCurved              = 0
nBilinear                  = 0
nCurved                    = 0
! zero counter for elem types
nCurvedElems               = 0
nLinearElems               = 0

#if USE_MPI
firstElem = offsetElem + 1
lastElem  = offsetElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif


DO iGlobalElem = firstElem,lastElem
  CNElemID = GetCNElemID(iGlobalElem)
  IF (ElemCurved(CNElemID)) THEN
    nCurvedElems = nCurvedElems+1
  ELSE
    nLinearElems = nLinearElems+1
  END IF

  DO ilocSide = 1,6
    ! ignore small mortar sides attached to big mortar sides
    SideID   = GetGlobalNonUniqueSideID(iGlobalElem,ilocSide)
    CNSideID = GetCNSideID(SideID)
    SELECT CASE(SideType(CNSideID))
      CASE (PLANAR_RECT)
        nPlanarRectangular    = nPlanarRectangular   +1
      CASE (PLANAR_NONRECT)
        nPlanarNonRectangular = nPlanarNonRectangular+1
      CASE (BILINEAR)
        nBilinear             = nBilinear            +1
      CASE (PLANAR_CURVED)
        nPlanarCurved         = nPlanarCurved        +1
      CASE (CURVED)
        nCurved               = nCurved              +1
    END SELECT
  END DO
END DO

#if USE_MPI
CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nLinearElems         ,nLinearElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
CALL MPI_REDUCE(nCurvedElems         ,nCurvedElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
#else
nPlanarRectangularTot    = nPlanarRectangular
nPlanarNonRectangularTot = nPlanarNonRectangular
nBilinearTot             = nBilinear
nPlanarCurvedTot         = nPlanarCurved
nCurvedTot               = nCurved
nLinearElemsTot          = nLinearElems
nCurvedElemsTot          = nCurvedElems
#endif /* USE_MPI */

LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-rectangular     faces: ', nPlanarRectangulartot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of bi-linear              faces: ', nBilineartot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-curved          faces: ', nPlanarCurvedtot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 faces: ', nCurvedtot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of (bi-)linear            elems: ', nLinearElemsTot
LBWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 elems: ', nCurvedElemsTot

END SUBROUTINE IdentifyElemAndSideType


SUBROUTINE PointsEqual(N,Points1,Points2,IsNotEqual)
!===================================================================================================================================
! compute the distance between two data sets
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: N
REAL,INTENT(IN)           :: Points1(1:3,1:N)
REAL,INTENT(IN)           :: Points2(1:3,1:N)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL                   :: IsNotEqual
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i
!===================================================================================================================================

IsNotEqual=.FALSE.

DO i=1,N
  IF( ABS(Points1(1,i)-Points2(1,i)).GT.1e-14 .OR. &
      ABS(Points1(2,i)-Points2(2,i)).GT.1e-14 .OR. &
      ABS(Points1(3,i)-Points2(3,i)).GT.1e-14 ) THEN
    IsNotEqual=.TRUE.
    RETURN
  END IF
END DO ! i=0,N

END SUBROUTINE PointsEqual


SUBROUTINE WeirdElementCheck()
!===================================================================================================================================
! Calculate whether element edges intersect other sides
! If this is the case it means that part of the element is turned inside-out
! which results in a warning so the user can decide whether it is a problem that
! necessitates a new mesh.
! Fixing the problem would involve defining the bilinear edge between nodes 2 and 4
! (instead of 1 and 3). This information would need to be stored and used throughout
! the particle treatment. Additionally, since the edge would need to be changed
! for both neighboring elements, it is possible that both elements might have the problem
! hence no solution exists.
! tl;dr: Hard/maybe impossible to fix, hence only a warning is given so the user can decide
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: NodeCoords_Shared,ConcaveElemSide_Shared,ElemSideNodeID_Shared,WeirdElems
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, iLocSide, kLocSide, iNode
INTEGER,ALLOCATABLE :: WeirdElemNbrs(:)
REAL              :: vec(1:3), Node(1:3,1:4),det(1:3)
LOGICAL           :: WEIRD, TRICHECK, TRIABSCHECK
INTEGER           :: firstElem,lastElem
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' CHECKING FOR WEIRD ELEMENTS...'

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

ALLOCATE(WeirdElemNbrs(1:lastElem-firstElem+1))

WeirdElems = 0

! go through all CN elements
DO iElem = firstElem,lastElem
  WEIRD = .FALSE.
  DO iLocSide = 1,5  ! go through local sides
    IF (.not.WEIRD) THEN  ! if one is found there is no need to continue
      IF (ConcaveElemSide_Shared(iLocSide,iElem)) THEN  ! only concave elements need to be checked
        ! build vector from node 1 to node 3
        vec(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
               - NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1)
        ! check all other sides
        DO kLocSide = iLocSide + 1, 6
          IF (ConcaveElemSide_Shared(kLocSide,iElem)) THEN  ! only concave elements need to be checked
            ! build 4 vectors from point 1 of edge to 4 nodes of kLocSide
            DO iNode = 1,4
              Node(:,iNode) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                            - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
            END DO
            ! Compute whether any of the triangle intersects with the vector vec:
            ! If all three volumes built by the vector vec and the vectors Node
            ! are either positive or negative then there is an intersection

            ! Triangle 1 (Nodes 1,2,3)
            ! Only check this if neither point of vec is part of the triangle.
            ! If points of vec correspont to point 1 or 3 or triangle then both
            ! triangles can be skipped (triabscheck = true), else point 4 needs to be checked
            ! separately for triangle 2 (see below)
            TRICHECK = .FALSE.
            TRIABSCHECK = .FALSE.
            DO iNode = 1,3
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
            END DO
            IF (.not.TRICHECK) THEN
              det(1) = ((Node(2,1) * Node(3,2) - Node(3,1) * Node(2,2)) * vec(1)  + &
                        (Node(3,1) * Node(1,2) - Node(1,1) * Node(3,2)) * vec(2)  + &
                        (Node(1,1) * Node(2,2) - Node(2,1) * Node(1,2)) * vec(3))
              det(2) = ((Node(2,2) * Node(3,3) - Node(3,2) * Node(2,3)) * vec(1)  + &
                        (Node(3,2) * Node(1,3) - Node(1,2) * Node(3,3)) * vec(2)  + &
                        (Node(1,2) * Node(2,3) - Node(2,2) * Node(1,3)) * vec(3))
              det(3) = ((Node(2,3) * Node(3,1) - Node(3,3) * Node(2,1)) * vec(1)  + &
                        (Node(3,3) * Node(1,1) - Node(1,3) * Node(3,1)) * vec(2)  + &
                        (Node(1,3) * Node(2,1) - Node(2,3) * Node(1,1)) * vec(3))
              IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
              IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
            END IF

            ! Triangle 2 (Nodes 1,3,4)
            TRICHECK = .FALSE.
            IF (.not.TRIABSCHECK) THEN
              ! Node 4 needs to be checked separately (see above)
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              IF (.not.TRICHECK) THEN
                det(1) = ((Node(2,1) * Node(3,3) - Node(3,1) * Node(2,3)) * vec(1)  + &
                          (Node(3,1) * Node(1,3) - Node(1,1) * Node(3,3)) * vec(2)  + &
                          (Node(1,1) * Node(2,3) - Node(2,1) * Node(1,3)) * vec(3))
                det(2) = ((Node(2,3) * Node(3,4) - Node(3,3) * Node(2,4)) * vec(1)  + &
                          (Node(3,3) * Node(1,4) - Node(1,3) * Node(3,4)) * vec(2)  + &
                          (Node(1,3) * Node(2,4) - Node(2,3) * Node(1,4)) * vec(3))
                det(3) = ((Node(2,4) * Node(3,1) - Node(3,4) * Node(2,1)) * vec(1)  + &
                          (Node(3,4) * Node(1,1) - Node(1,4) * Node(3,1)) * vec(2)  + &
                          (Node(1,4) * Node(2,1) - Node(2,4) * Node(1,1)) * vec(3))
                IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
                IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END IF
  END DO
  IF (WEIRD) THEN
    WeirdElems = WeirdElems + 1
    WeirdElemNbrs(WeirdElems) = GetGlobalElemID(iElem)
  END IF
END DO

IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  END DO
  IPWRITE(Unit_StdOut,*) 'This check is optional. You can disable it by setting meshCheckWeirdElements = F'
  CALL abort(__STAMP__,'Weird elements found: it means that part of the element is turned inside-out')
END IF

LBWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'

DEALLOCATE(WeirdElemNbrs)

LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WeirdElementCheck


SUBROUTINE CalcXCL_NGeo()
!===================================================================================================================================
!> calculates Chebyshev-Lobatto basis XCL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis                  ,ONLY: BuildBezierVdm,BuildBezierDMat
USE MOD_Basis                  ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
USE MOD_Mesh_Vars              ,ONLY: NGeo,XCL_NGeo,wBaryCL_NGeo,XiCL_NGeo,Xi_NGeo
USE MOD_Mesh_Vars              ,ONLY: wBaryCL_NGeo1,Vdm_CLNGeo1_CLNGeo,XiCL_NGeo1
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: Vdm_Bezier,sVdm_Bezier,D_Bezier
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems,offsetElem,nElems
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
#if USE_MPI
INTEGER                :: iElem
#endif /*USE_MPI*/
REAL                   :: Vdm_NGeo_CLNGeo(0:NGeo,0:NGeo)
INTEGER                :: i
REAL,DIMENSION(0:NGeo) :: wBary_NGeo
!===================================================================================================================================

! Equidistant-Lobatto NGeo
ALLOCATE(Xi_NGeo(0:NGeo))
DO i=0,NGeo
  Xi_NGeo(i) = 2./REAL(NGeo) * REAL(i) - 1.
END DO
CALL BarycentricWeights(NGeo,Xi_NGeo,wBary_NGeo)

! small wBaryCL_NGEO
ALLOCATE(Vdm_CLNGeo1_CLNGeo(0:NGeo,0:1))

! new for curved particle sides
ALLOCATE( Vdm_Bezier(0:NGeo,0:NGeo))
ALLOCATE(sVdm_Bezier(0:NGeo,0:NGeo))
ALLOCATE(   D_Bezier(0:NGeo,0:NGeo))

! XiCL_NGeo1, wBaryCL_NGeo1
ALLOCATE(wBaryCL_NGeo1(0:1))
ALLOCATE(XiCL_NGeo1(0:1))
CALL ChebyGaussLobNodesAndWeights(1,XiCL_NGeo1)
CALL BarycentricWeights(1,XiCL_NGeo1,wBaryCL_NGeo1)
! Vandermonde: Vdm_CLNGeo1_CLNGeo
CALL InitializeVandermonde(1,NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo,Vdm_CLNGeo1_CLNGeo)
! initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm(NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier) ! Required for FacNchooseK (requires BezierElevation)
CALL BuildBezierDMat(NGeo,Xi_NGeo,D_Bezier) ! Required for DMat=D_Bezier (requires FacNchooseK)
CALL InitializeVandermonde(NGeo,NGeo,wBaryCL_NGeo,Xi_NGeo,XiCL_NGeo,Vdm_NGeo_CLNGeo)

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared((/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => XCL_NGeo_Array

DO iElem = 1,nElems
  XCL_NGeo_Shared(:,:,:,:,offsetElem+iElem) = XCL_NGeo(:,:,:,:,iElem)
END DO ! iElem = 1, nElems

! Communicate XCL and dXCL between compute node roots instead of calculating globally
CALL BARRIER_AND_SYNC(XCL_NGeo_Shared_Win ,MPI_COMM_SHARED)

IF (nComputeNodeProcessors.NE.nProcessors_Global .AND. myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
      , 0                             &
      , MPI_DATATYPE_NULL             &
      , XCL_NGeo_Shared               &
      , 3*(NGeo+1)**3*recvcountElem   &
      , 3*(NGeo+1)**3*displsElem      &
      , MPI_DOUBLE_PRECISION          &
      , MPI_COMM_LEADERS_SHARED       &
      , IERROR)
END IF

CALL BARRIER_AND_SYNC(XCL_NGeo_Shared_Win ,MPI_COMM_SHARED)
#else
XCL_NGeo_Shared  => XCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcXCL_NGeo


SUBROUTINE CalcParticleMeshMetrics()
!===================================================================================================================================
!> calculates dXCL_Ngeo for global mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars              ,ONLY: N_VolMesh
USE MOD_Mesh_Vars              ,ONLY: dXCL_NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_DG_Vars                ,ONLY: nDofsMapping, N_DG_Mapping
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nElems
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: NGeo,nGlobalElems
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_DG_Vars                ,ONLY: nDofsMapping, N_DG_Mapping, recvcountDofs, displsDofs
!#else
!USE MOD_Mesh_Vars              ,ONLY: N_VolMesh_Shared
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!#if USE_MPI
INTEGER                        :: iElem, Nloc, i,j,k,r, offSetDof
!#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared((/3*  nDofsMapping/), Elem_xGP_Shared_Win,Elem_xGP_Array)
CALL Allocate_Shared((/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
Elem_xGP_Shared (1:3    ,1:nDofsMapping) => Elem_xGP_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => dXCL_NGeo_Array

DO iElem = 1,nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  offSetDof = N_DG_Mapping(1,iElem+offSetElem)
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
        Elem_xGP_Shared (:,r+offSetDof) = N_VolMesh(iElem)%Elem_xGP (:,i,j,k)
      END DO
    END DO
  END DO
  dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
END DO ! iElem = 1, nElems

! Communicate XCL and dXCL between compute node roots instead of calculating globally
CALL BARRIER_AND_SYNC(Elem_xGP_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(dXCL_NGeo_Shared_Win,MPI_COMM_SHARED)

IF (nComputeNodeProcessors.NE.nProcessors_Global .AND. myComputeNodeRank.EQ.0) THEN

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , Elem_xGP_Shared               &
                     , 3*recvcountDofs   &
                     , 3*displsDofs      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , dXCL_NGeo_Shared              &
                     , 3*3*(NGeo+1)**3*recvcountElem &
                     , 3*3*(NGeo+1)**3*displsElem    &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
END IF

CALL BARRIER_AND_SYNC(Elem_xGP_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(dXCL_NGeo_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(Elem_xGP_Array(3*nDofsMapping))
Elem_xGP_Shared(1:3,1:nDofsMapping) => Elem_xGP_Array
DO iElem = 1,nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  offSetDof = N_DG_Mapping(1,iElem+offSetElem)
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
        Elem_xGP_Shared(:,r+offSetDof) = N_VolMesh(iElem)%Elem_xGP(:,i,j,k)
      END DO
    END DO
  END DO
END DO ! iElem = 1, nElems
dXCL_NGeo_Shared => dXCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcParticleMeshMetrics


SUBROUTINE CalcBezierControlPoints()
!===================================================================================================================================
!> calculate the Bezier control point (+elevated) for shared compute-node mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoElevated
USE MOD_Particle_Mesh_Vars     ,ONLY: nNonUniqueGlobalSides,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Surfaces      ,ONLY: GetBezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3DElevated,BezierElevation
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
USE MOD_Particle_Mesh_Vars     ,ONLY: BezierControlPoints3D_Shared,BezierControlPoints3D_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: BezierControlPoints3DElevated_Shared,BezierControlPoints3DElevated_Shared_Win
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iSide,ilocSide
INTEGER                        :: SideID
INTEGER                        :: firstElem,lastElem,firstSide,lastSide
INTEGER                        :: p,q,pq(2)
REAL                           :: tmp(3,0:NGeo,0:NGeo)
REAL                           :: tmp2(3,0:NGeo,0:NGeo)
#if !USE_MPI
INTEGER                        :: ALLOCSTAT
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! BezierControlPoints are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

LBWRITE(UNIT_stdOut,'(A)') ' CALCULATING BezierControlPoints ...'

! Build BezierControlPoints3D (compute-node local+halo)
#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared((/3*(NGeo+1)*(NGeo+1)*nNonUniqueGlobalSides/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3D_Shared_Win,IERROR)
BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) => BezierControlPoints3D_Shared
IF (myComputeNodeRank.EQ.0) THEN
  BezierControlPoints3D         = 0.
END IF
IF (BezierElevation.GT.0) THEN
  ! Sanity check
  IF(NGeoElevated.LT.0) CALL abort(__STAMP__,'NGeoElevated<0 is not allowed. Check correct initialisation of NGeoElevated.')
  CALL Allocate_Shared((/3*(NGeoElevated+1)*(NGeoElevated+1)*nNonUniqueGlobalSides/),BezierControlPoints3DElevated_Shared_Win,BezierControlPoints3DElevated_Shared)
  CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3DElevated_Shared_Win,IERROR)
  BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides) => BezierControlPoints3DElevated_Shared
  IF (myComputeNodeRank.EQ.0) THEN
    BezierControlPoints3DElevated = 0.
  END IF
END IF
#else
ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) &
        ,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3D!')
BezierControlPoints3D         = 0.

IF (BezierElevation.GT.0) THEN
  ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3DElevated!')
  BezierControlPoints3DElevated = 0.
END IF
#endif /*USE_MPI*/

#if USE_MPI
CALL BARRIER_AND_SYNC(BezierControlPoints3D_Shared_Win,MPI_COMM_SHARED)
IF (BezierElevation.GT.0) THEN
  CALL BARRIER_AND_SYNC(BezierControlPoints3DElevated_Shared_Win,MPI_COMM_SHARED)
END IF

firstElem = INT(REAL( myComputeNodeRank   )*REAL(nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nGlobalElems)/REAL(nComputeNodeProcessors))
firstSide = INT(REAL (myComputeNodeRank   )*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! iElem is CN elem
DO iElem = firstElem,lastElem
  DO ilocSide=1,6
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , 0    , :    , :   ,iElem )
    CASE(XI_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , NGeo , :    , :   ,iElem )
    CASE(ETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , 0    , :   ,iElem )
    CASE(ETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , NGeo , :   ,iElem )
    CASE(ZETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , 0   ,iElem )
    CASE(ZETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , NGeo,iElem )
    END SELECT
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,tmp,tmp2)

    ! get global SideID of local side
    SideID = GetGlobalNonUniqueSideID(iElem,iLocSide)

    DO q=0,NGeo; DO p=0,NGeo
      ! turn into right hand system of side
      pq = CGNS_SideToVol2(NGeo,p,q,iLocSide)
      BezierControlPoints3D(1:3,p,q,SideID) = tmp2(1:3,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! ilocSide=1,6
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL BARRIER_AND_SYNC(BezierControlPoints3D_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

! Calculate elevated BezierControlPoints
IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,LastSide
    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
    ! Indices in shared arrays are shifted by 1
    CALL GetBezierControlPoints3DElevated( NGeo,NGeoElevated                                                       &
                                         , BezierControlPoints3D        (1:3,0:NGeo        ,0:NGeo        ,iSide)  &
                                         , BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide))
  END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(BezierControlPoints3DElevated_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
END IF

END SUBROUTINE CalcBezierControlPoints


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
!> Subroutine for particle geometry initialization (GEO container):
!> - ConcaveElemSide_Shared(   1:6,1:nComputeNodeTotalElems)
!> - ElemSideNodeID_Shared(1:4,1:6,1:nComputeNodeTotalElems)
!> - ElemMidPoint_Shared(      1:3,1:nComputeNodeTotalElems)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID, GetCornerNodeMapCGNS
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif
USE MOD_ReadInTools        ,ONLY: GETLOGICAL
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,FirstElem,LastElem,GlobalElemID
INTEGER            :: GlobalSideID,nlocSides,localSideID,iLocSide
INTEGER            :: iNode, NodeMap(1:4,1:6), nStart, NodeNum
REAL               :: A(3,3),detcon
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION...'

! Get the node map to convert from the CGNS format (as given by HOPR)
CALL GetCornerNodeMapCGNS(NGeo,NodeMapCGNS=NodeMap(1:4,1:6))

#if USE_MPI
CALL Allocate_Shared((/6,nComputeNodeTotalElems/),ConcaveElemSide_Shared_Win,ConcaveElemSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,ConcaveElemSide_Shared_Win,IERROR)
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

CALL Allocate_Shared((/4,6,nComputeNodeTotalElems/),ElemSideNodeID_Shared_Win,ElemSideNodeID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemSideNodeID_Shared_Win,IERROR)

CALL Allocate_Shared((/3,nComputeNodeTotalElems/),ElemMidPoint_Shared_Win,ElemMidPoint_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemMidPoint_Shared_Win,IERROR)
#else
ALLOCATE(ConcaveElemSide_Shared(   1:6,1:nElems))
ALLOCATE(ElemSideNodeID_Shared(1:4,1:6,1:nElems))
ALLOCATE(ElemMidPoint_Shared(      1:3,1:nElems))
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  ElemSideNodeID_Shared  = 0
  ConcaveElemSide_Shared = .FALSE.
#if USE_MPI
END IF
CALL BARRIER_AND_SYNC(ElemSideNodeID_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ConcaveElemSide_Shared_Win,MPI_COMM_SHARED)
#endif

! iElem is CNElemID
DO iElem = firstElem,lastElem
  GlobalElemID = GetGlobalElemID(iElem)

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    ! Get global SideID
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide

    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    IF (localSideID.LE.0) CYCLE
    ! Find start of CGNS mapping from flip
    IF (SideInfo_Shared(SIDE_ID,GlobalSideID).GT.0) THEN
      nStart = 0
    ELSE
      nStart = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1)
    END IF
    ! Shared memory array starts at 1, but NodeID at 0
    ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart  ,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)

  END DO ! iLocSide = 1,nlocSides
END DO ! iElem = firstElem,lastElem

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemSideNodeID_Shared_Win,MPI_COMM_SHARED)
#endif

!--- Save whether Side is concave or convex
DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    IF (localSideID.LE.0) CYCLE
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
       A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,iElem)+1) &
                    - NodeCoords_Shared(:,ElemSideNodeID_Shared(4      ,localSideID,iElem)+1)
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
    IF (detcon.LT.0) THEN
      ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    ELSE IF (detcon.EQ.0.0) THEN
      IF (GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    END IF
  END DO
END DO

DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  ElemMidPoint_Shared(:,iElem) = 0.
  DO iNode = 1,8
    ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+iNode)
  END DO
  ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) / 8.
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(ConcaveElemSide_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemMidPoint_Shared_Win   ,MPI_COMM_SHARED)
#endif

!--- check for elements with intersecting sides (e.g. very flat elements)
meshCheckWeirdElements = GETLOGICAL('meshCheckWeirdElements','.TRUE.')
IF(meshCheckWeirdElements) CALL WeirdElementCheck()

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GEOMETRY INFORMATION DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleGeometry


SUBROUTINE ComputePeriodicVec()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: NGeo,offsetElem,BoundaryType
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Vars          ,ONLY: PartMeshHasPeriodicBCs
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Mesh_Tools             ,ONLY: GetCornerNodeMapCGNS
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: iNode=1
INTEGER                        :: iVec,iBC,iPartBC
INTEGER                        :: firstElem,lastElem,NbSideID,BCALPHA,flip
INTEGER                        :: SideID,iGlobalElem,NbElemID,localSideID,localSideNbID,nStart,NodeMap(4,6)
REAL,DIMENSION(3)              :: MasterCoords,SlaveCoords,Vec
LOGICAL,ALLOCATABLE            :: PeriodicFound(:)
#if USE_MPI
REAL                           :: sendbuf
REAL,ALLOCATABLE               :: recvbuf(:)
#endif
INTEGER                        :: nPeriodicVectorsParameterIni
!-----------------------------------------------------------------------------------------------------------------------------------

! Find number of periodic vectors
nPeriodicVectorsParameterIni = GEO%nPeriodicVectors
GEO%nPeriodicVectors = MERGE(MAXVAL(BoundaryType(:,BC_ALPHA)),0,PartMeshHasPeriodicBCs)
IF(nPeriodicVectorsParameterIni.GT.GEO%nPeriodicVectors)THEN
  SWRITE (*,*) "Number of periodic vectors in parameter file: ", nPeriodicVectorsParameterIni
  SWRITE (*,*) "Number of periodic vectors in mesh      file: ", GEO%nPeriodicVectors
  CALL CollectiveStop(__STAMP__,'Wrong number of periodic vectors!')
END IF
IF (GEO%nPeriodicVectors.EQ.0) RETURN

firstElem = offsetElem+1
lastElem  = offsetElem+nElems

ALLOCATE(PeriodicFound(1:GEO%nPeriodicVectors))
PeriodicFound(:) = .FALSE.

ALLOCATE(GEO%PeriodicVectors(1:3,GEO%nPeriodicVectors))
GEO%PeriodicVectors = 0.

CALL GetCornerNodeMapCGNS(NGeo,NodeMapCGNS=NodeMap)

DO iGlobalElem = firstElem,lastElem
  ! Every periodic vector already found
  IF (ALL(PeriodicFound(:))) EXIT

SideLoop: DO SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobalElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobalElem)
    ! Get BC
    iBC = SideInfo_Shared(SIDE_BCID,SideID)
    IF(iBC.EQ.0) CYCLE

    ! Get particle BC
    iPartBC = PartBound%MapToPartBC(iBC)
    IF (iPartBC.EQ.0) CYCLE

    ! Boundary is a periodic boundary
    IF (PartBound%TargetBoundCond(iPartBC).NE.3) CYCLE

    ! Check if side is master side
    BCALPHA = BoundaryType(iBC,BC_ALPHA)

    IF (BCALPHA.GT.0) THEN
      ! Periodic vector already found
      IF (PeriodicFound(BCALPHA)) CYCLE

      ! Periodic slave side has same ID, but negative sign
      flip          = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)
      NbElemID      = SideInfo_Shared(SIDE_NBELEMID,SideID)
      localSideID   = SideInfo_Shared(SIDE_LOCALID,SideID)
      localSideNbID = (SideInfo_Shared(SIDE_FLIP,SideID)-flip)/10
      NbSideID      = GetGlobalNonUniqueSideID(NbElemID,localSideNbID)
      nStart        = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,NbSideID),10)-1)

      ! Only take the first node into account, no benefit in accuracy if running over others as well
      MasterCoords  = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,iGlobalElem)+NodeMap(iNode                  ,localSideID))
      SlaveCoords   = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,NbElemID)+NodeMap(MOD(nStart+5-iNode,4)+1,localSideNbID))
      Vec           = SlaveCoords-MasterCoords

      ! Might consider aborting here, malformed periodic sides
      IF (VECNORM(Vec).EQ.0) CYCLE

      ! Check if the periodic vector is ALMOST aligned with a Cartesian direction
      DO iVec = 1,3
        ! IF (ABS(Vec(iVec)).GT.0 .AND. ABS(Vec(iVec))*VECNORM(Vec).LT.1E-12) CYCLE SideLoop
        IF (ABS(Vec(iVec)).GT.0 .AND. ABS(Vec(iVec)).LT.1E-12*VECNORM(Vec)) Vec(iVec) = 0.
      END DO

      GEO%PeriodicVectors(:,BCALPHA) = Vec
      PeriodicFound(BCALPHA) = .TRUE.
    END IF
  END DO SideLoop
END DO

#if USE_MPI
ALLOCATE(recvbuf(0:nProcessors-1))
sendbuf = 0.
recvbuf = 0.

DO iVec = 1,GEO%nPeriodicVectors
  sendbuf = MERGE(VECNORM(GEO%PeriodicVectors(:,iVec)),HUGE(1.),PeriodicFound(iVec))

! Do it by hand, MPI_ALLREDUCE seems problematic with MPI_2DOUBLE_PRECISION and MPI_MINLOC
! https://stackoverflow.com/questions/56307320/mpi-allreduce-not-synchronizing-properly
!CALL MPI_ALLREDUCE(MPI_IN_PLACE,sendbuf,GEO%nPeriodicVectors,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_SHARED,iERROR)

  CALL MPI_ALLGATHER(sendbuf,1,MPI_DOUBLE_PRECISION,recvbuf(:),1,MPI_DOUBLE_PRECISION,MPI_COMM_PICLAS,iERROR)
  IF (ALL(recvbuf(:).EQ.HUGE(1.))) CALL CollectiveStop(__STAMP__,'No periodic vector for BC_ALPHA found!',IntInfo=iVec)

  ! MINLOC does not follow array bounds, so root rank = 1
  CALL MPI_BCAST(GEO%PeriodicVectors(:,iVec),3,MPI_DOUBLE_PRECISION,MINLOC(recvbuf(:),1)-1,MPI_COMM_PICLAS,iError)
END DO
#endif /*USE_MPI*/

SDEALLOCATE(PeriodicFound)

#if USE_MPI
IF (myRank.EQ.0) THEN
#endif /*USE_MPI*/
  LBWRITE(UNIT_stdOut,'(A,I0,A)') ' Found ',GEO%nPeriodicVectors,' periodic vectors for particle tracking'
  DO iVec = 1,GEO%nPeriodicVectors
    LBWRITE(UNIT_stdOut,'(A,I1,A,F12.8,2(", ",F12.8))') ' | Periodic vector ',iVec,': ', GEO%PeriodicVectors(:,iVec)
  END DO

  ! Sanity check
  DO iVec = 1,GEO%nPeriodicVectors
    IF(VECNORM(GEO%PeriodicVectors(:,iVec)).LE.0.)THEN
      CALL abort(__STAMP__,'Norm of GEO%PeriodicVectors(:,iVec) <= 0 for iVec =',IntInfoOpt=iVec)
    END IF
  END DO
#if USE_MPI
END IF
#endif /*USE_MPI*/

END SUBROUTINE ComputePeriodicVec


SUBROUTINE InitVolumes_2D()
!===================================================================================================================================
!> Routine determines the symmetry sides and calculates the 2D (area faces in symmetry plane) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side (SymmetrySide array) will be used later on to determine in which two
!> directions the quadtree shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement. Additionally,
!> symmetry sides will be skipped during tracking (SideIsSymSide_Shared array).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: Pi
USE MOD_PreProc
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem,nBCSides,SideToElem
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO,LocalVolume,MeshVolume
USE MOD_Particle_Mesh_Vars      ,ONLY: SymmetrySide
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared, SideIsSymSide
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle
#if USE_MPI
USE MOD_Mesh_Vars               ,ONLY: ELEM_HALOFLAG
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: nNonUniqueGlobalSides, offsetComputeNodeElem, ElemInfo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared_Win,ElemCharLength_Shared_Win,SideIsSymSide_Shared,SideIsSymSide_Shared_Win
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors,MPI_COMM_LEADERS_SHARED
#else
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeSides
#endif /*USE_MPI*/
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, iLocSide, iNode, BCSideID, locElemID, CNElemID, iSide
REAL                            :: radius, triarea(2)
#if USE_MPI
REAL                            :: CNVolume
INTEGER                         :: offsetElemCNProc
#endif /*USE_MPI*/
LOGICAL                         :: SymmetryBCExists
INTEGER                         :: firstElem, lastElem, firstSide, lastSide
!===================================================================================================================================

#if USE_MPI
CALL Allocate_Shared((/nNonUniqueGlobalSides/),SideIsSymSide_Shared_Win,SideIsSymSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideIsSymSide_Shared_Win,IERROR)
SideIsSymSide => SideIsSymSide_Shared
! only CN root nullifies
IF(myComputeNodeRank.EQ.0) SideIsSymSide = .FALSE.
! This sync/barrier is required as it cannot be guaranteed that the zeros have been written to memory by the time the MPI_REDUCE
! is executed (see MPI specification). Until the Sync is complete, the status is undefined, i.e., old or new value or utter nonsense.
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(SideIsSymSide(nComputeNodeSides))
SideIsSymSide = .FALSE.
#endif  /*USE_MPI*/

! Flag of symmetry sides to be skipped during tracking
#if USE_MPI
  firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
  lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
  firstSide = 1
  lastSide  = nComputeNodeSides
#endif

DO iSide = firstSide, lastSide
  SideIsSymSide(iSide) = .FALSE.
  ! ignore non-BC sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
#if USE_MPI
  ! ignore sides outside of halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.0) CYCLE
#endif /*USE_MPI*/
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,iSide))).EQ.PartBound%SymmetryDim) &
    SideIsSymSide(iSide) = .TRUE.
END DO

SymmetryBCExists = .FALSE.
ALLOCATE(SymmetrySide(1:nElems,1:2))                ! 1: GlobalSide, 2: LocalSide
SymmetrySide = -1

! Sanity check: mesh has to be centered at z = 0
IF(.NOT.ALMOSTEQUALRELATIVE(GEO%zmaxglob,ABS(GEO%zminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in z:', GEO%zmaxglob
  SWRITE(*,*) 'Minimum dimension in z:', GEO%zminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%zmaxglob)-ABS(GEO%zminglob))/ABS(GEO%zminglob), ' > 1e-5'
  CALL abort(__STAMP__,'ERROR: Please orient your mesh with one cell in z-direction around 0, |z_min| = z_max !')
END IF

! Calculation of the correct volume and characteristic length
DO BCSideID=1,nBCSides
  locElemID = SideToElem(S2E_ELEM_ID,BCSideID)
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  SideID=GetGlobalNonUniqueSideID(offsetElem+locElemID,iLocSide)
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%SymmetryDim) THEN
    CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
    iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
    ! Exclude the symmetry axis (y=0)
    IF(Symmetry%Axisymmetric) THEN
      IF(MAXVAL(NodeCoords_Shared(2,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1)).LE.0.0) CYCLE
    END IF
    ! The z-plane with the positive z component is chosen
    IF(MINVAL(NodeCoords_Shared(3,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1)).GT.(GEO%zmaxglob+GEO%zminglob)/2.) THEN
      IF(SymmetrySide(locElemID,1).GT.0) THEN
        CALL abort(__STAMP__&
          ,'ERROR: PICLas could not determine a unique symmetry surface for 2D/axisymmetric calculation!'//&
          ' Please orient your mesh with x as the symmetry axis and positive y as the second/radial direction!')
      END IF
      SymmetrySide(locElemID,1) = BCSideID
      SymmetrySide(locElemID,2) = iLocSide
      ! The volume calculated at this point (final volume for the 2D case) corresponds to the cell face area (z-dimension=1) in
      ! the xy-plane.
      ElemVolume_Shared(CNElemID) = 0.0
      CALL CalcNormAndTangTriangle(area=triarea(1),TriNum=1, SideID=SideID)
      CALL CalcNormAndTangTriangle(area=triarea(2),TriNum=2, SideID=SideID)
      ElemVolume_Shared(CNElemID) = triarea(1) + triarea(2)
      ! Characteristic length is compared to the mean free path as the condition to refine the mesh. For the 2D/axisymmetric case
      ! the third dimension is not considered as particle interaction occurs in the xy-plane, effectively reducing the refinement
      ! requirement.
      ElemCharLength_Shared(CNElemID) = SQRT(ElemVolume_Shared(CNElemID))
      ! Axisymmetric case: The volume is multiplied by the circumference to get the volume of the ring. The cell face in the
      ! xy-plane is rotated around the x-axis. The radius is the middle point of the cell face.
      IF (Symmetry%Axisymmetric) THEN
        radius = 0.
        DO iNode = 1, 4
          radius = radius + NodeCoords_Shared(2,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)
        END DO
        radius = radius / 4.
        ElemVolume_Shared(CNElemID) = ElemVolume_Shared(CNElemID) * 2. * Pi * radius
      END IF
      SymmetryBCExists = .TRUE.
    END IF      ! Greater z-coord
  END IF
END DO

IF(.NOT.SymmetryBCExists) CALL abort(__STAMP__,'At least one symmetric BC (in the xy-plane) has to be defined for 2D simulations')

! LocalVolume & MeshVolume: Recalculate the volume of the mesh of a single process and the total mesh volume
#if USE_MPI
FirstElem = offsetElem - offsetComputeNodeElem + 1
LastElem  = offsetElem - offsetComputeNodeElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

LocalVolume = SUM(ElemVolume_Shared(FirstElem:LastElem))

#if USE_MPI
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
! Compute-node mesh volume
offsetElemCNProc = offsetElem - offsetComputeNodeElem
CNVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))
CALL MPI_ALLREDUCE(MPI_IN_PLACE,CNVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SHARED,iError)
IF (myComputeNodeRank.EQ.0) THEN
  ! All-reduce between node leaders
  CALL MPI_ALLREDUCE(CNVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
! Broadcast from node leaders to other processors on the same node
CALL MPI_BCAST(MeshVolume,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iERROR)
#else
MeshVolume = LocalVolume
#endif /*USE_MPI*/

END SUBROUTINE InitVolumes_2D


SUBROUTINE InitVolumes_1D()
!===================================================================================================================================
!> Routine determines a symmetry side and calculates the 1D (area faces at x axis) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side will be used later on to determine in which two directions the quadtree
!> shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: nElems, offsetElem
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO,LocalVolume,MeshVolume, SideIsSymSide
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_Mesh_Vars               ,ONLY: ELEM_HALOFLAG
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeElem,SideIsSymSide_Shared,SideIsSymSide_Shared_Win,nNonUniqueGlobalSides
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared_Win,ElemCharLength_Shared_Win, ElemInfo_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors,MPI_COMM_LEADERS_SHARED
#else
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iLocSide, iElem, SideID, iOrder, CNElemID, iSide
REAL                            :: X(2), MaxCoord, MinCoord
LOGICAL                         :: SideInPlane, X1Occupied
INTEGER                         :: firstElem,lastElem, firstSide, lastSide
#if USE_MPI
REAL                            :: CNVolume
INTEGER                         :: offsetElemCNProc
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
firstElem = offsetElem - offsetComputeNodeElem + 1
lastElem  = offsetElem - offsetComputeNodeElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

ALLOCATE(GEO%XMinMax(2,nElems))

IF(.NOT.ALMOSTEQUALRELATIVE(GEO%ymaxglob,ABS(GEO%yminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in y:', GEO%ymaxglob
  SWRITE(*,*) 'Minimum dimension in y:', GEO%yminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%ymaxglob)-ABS(GEO%yminglob))/ABS(GEO%yminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in y-direction around 0, |y_min| = y_max !')
END IF

IF(.NOT.ALMOSTEQUALRELATIVE(GEO%zmaxglob,ABS(GEO%zminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in z:', GEO%zmaxglob
  SWRITE(*,*) 'Minimum dimension in z:', GEO%zminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%zmaxglob)-ABS(GEO%zminglob))/ABS(GEO%zminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in z-direction around 0, |z_min| = z_max !')
END IF

DO iElem = 1,nElems
  ! Check if all sides of the element are parallel to xy-, xz-, or yz-plane and Sides parallel to xy-,and xz-plane are symmetric
  ! And determine xmin and xmax of the element
  X = 0
  X1Occupied = .FALSE.
  DO iLocSide = 1,6
    SideInPlane=.FALSE.
    SideID = GetGlobalNonUniqueSideID(offsetElem+iElem,iLocSide)
    CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
    DO iOrder=1,3
      MaxCoord = MAXVAL(NodeCoords_Shared(iOrder,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1))
      MinCoord = MINVAL(NodeCoords_Shared(iOrder,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1))
      IF(ALMOSTALMOSTEQUAL(MaxCoord,MinCoord).OR.(ALMOSTZERO(MaxCoord).AND.ALMOSTZERO(MinCoord))) THEN
        IF(SideInPlane) CALL abort(__STAMP__&
          ,'ERROR: Please orient your mesh with all element sides parallel to xy-,xz-,or yz-plane')
        SideInPlane=.TRUE.
        IF(iOrder.GE.2)THEN
          IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).NE.PartBound%SymmetryDim) &
            CALL abort(__STAMP__,&
              'ERROR: Sides parallel to xy-,and xz-plane has to be the symmetric boundary condition')
        END IF
        IF(iOrder.EQ.1) THEN
          IF(X1Occupied) THEN
            X(2) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
          ELSE
            X(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
            X1Occupied = .TRUE.
          END IF
        END IF !iOrder.EQ.1
      END IF
    END DO !iOrder=1,3
    IF(.NOT.SideInPlane) THEN
      IPWRITE(*,*) 'ElemID:',iElem,'SideID:',iLocSide
      DO iOrder=1,4
        IPWRITE(*,*) 'Node',iOrder,'x:',NodeCoords_Shared(1,ElemSideNodeID_Shared(iOrder,iLocSide,CNElemID)+1), &
                                  'y:',NodeCoords_Shared(2,ElemSideNodeID_Shared(iOrder,iLocSide,CNElemID)+1), &
                                  'z:',NodeCoords_Shared(3,ElemSideNodeID_Shared(iOrder,iLocSide,CNElemID)+1)
      END DO
      CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with all element sides parallel to xy-,xz-,or yz-plane')
    END IF
  END DO ! iLocSide = 1,6
  GEO%XMinMax(1,iElem) = MINVAL(X)
  GEO%XMinMax(2,iElem) = MAXVAL(X)
  ElemVolume_Shared(CNElemID) = ABS(X(1)-X(2))
  ElemCharLength_Shared(CNElemID) = ElemVolume_Shared(CNElemID)
END DO
! LocalVolume & MeshVolume: Recalculate the volume of the mesh of a single process and the total mesh volume
LocalVolume = SUM(ElemVolume_Shared(firstElem:lastElem))
#if USE_MPI
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
! Compute-node mesh volume
offsetElemCNProc = offsetElem - offsetComputeNodeElem
CNVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))
CALL MPI_ALLREDUCE(MPI_IN_PLACE,CNVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SHARED,iError)
IF (myComputeNodeRank.EQ.0) THEN
  ! All-reduce between node leaders
  CALL MPI_ALLREDUCE(CNVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
! Broadcast from node leaders to other processors on the same node
CALL MPI_BCAST(MeshVolume,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iERROR)
#else
MeshVolume = SUM(ElemVolume_Shared)
#endif /*USE_MPI*/

#if USE_MPI
CALL Allocate_Shared((/nNonUniqueGlobalSides/),SideIsSymSide_Shared_Win,SideIsSymSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideIsSymSide_Shared_Win,IERROR)
SideIsSymSide => SideIsSymSide_Shared
! only CN root nullifies
IF(myComputeNodeRank.EQ.0) SideIsSymSide = .FALSE.
! This sync/barrier is required as it cannot be guaranteed that the zeros have been written to memory by the time the MPI_REDUCE
! is executed (see MPI specification). Until the Sync is complete, the status is undefined, i.e., old or new value or utter nonsense.
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(SideIsSymSide(nComputeNodeSides))
SideIsSymSide = .FALSE.
#endif  /*USE_MPI*/

! Flag of symmetry sides to be skipped during tracking
#if USE_MPI
  firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
  lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
  firstSide = 1
  lastSide  = nComputeNodeSides
#endif

DO iSide = firstSide, lastSide
  SideIsSymSide(iSide) = .FALSE.
  ! ignore non-BC sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
#if USE_MPI
  ! ignore sides outside of halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.0) CYCLE
#endif /*USE_MPI*/
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,iSide))).EQ.PartBound%SymmetryDim) &
    SideIsSymSide(iSide) = .TRUE.
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win ,MPI_COMM_SHARED)
#endif

END SUBROUTINE InitVolumes_1D


REAL FUNCTION DSMC_2D_CalcSymmetryArea(iLocSide,iElem, ymin, ymax)
!===================================================================================================================================
!> Calculates the actual area of an element for 2D simulations (plane/axisymmetric) regardless of the mesh dimension in z
!> Utilized in the particle emission (surface flux) and boundary sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem,iLocSide           !> iElem is the compute-node element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, OPTIONAL, INTENT(OUT)   :: ymax,ymin
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode
REAL                          :: P(1:2,1:4), Pmin(2), Pmax(2), Length, MidPoint
!===================================================================================================================================

Pmin = HUGE(Pmin)
Pmax = -HUGE(Pmax)

DO iNode = 1,4
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1)
END DO

Pmax(1) = MAXVAL(P(1,:))
Pmax(2) = MAXVAL(P(2,:))
Pmin(1) = MINVAL(P(1,:))
Pmin(2) = MINVAL(P(2,:))

IF (PRESENT(ymax).AND.PRESENT(ymin)) THEN
  ymin = Pmin(2)
  ymax = Pmax(2)
END IF

Length = SQRT((Pmax(1)-Pmin(1))**2 + (Pmax(2)-Pmin(2))**2)

MidPoint = (Pmax(2)+Pmin(2)) / 2.
IF(Symmetry%Axisymmetric) THEN
  DSMC_2D_CalcSymmetryArea = Length * MidPoint * Pi * 2.
  ! Area of the cells on the rotational symmetry axis is set to one
  IF(.NOT.(DSMC_2D_CalcSymmetryArea.GT.0.0)) DSMC_2D_CalcSymmetryArea = 1.
ELSE
  DSMC_2D_CalcSymmetryArea = Length
END IF
RETURN

END FUNCTION DSMC_2D_CalcSymmetryArea


REAL FUNCTION DSMC_1D_CalcSymmetryArea(iLocSide,iElem)
!===================================================================================================================================
!> Calculates the actual area of an element for 1D simulations regardless of the mesh dimension in z and y
!> Utilized in the particle emission (surface flux) and boundary sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem,iLocSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Pmin, Pmax, Length
!===================================================================================================================================

Pmax = MAXVAL(NodeCoords_Shared(1,ElemSideNodeID_Shared(:,iLocSide,iElem)+1))
Pmin = MINVAL(NodeCoords_Shared(1,ElemSideNodeID_Shared(:,iLocSide,iElem)+1))

Length = ABS(Pmax-Pmin)

! IF(Symmetry%Axisymmetric) THEN
!   MidPoint = (Pmax(2)+Pmin(2)) / 2.
!   DSMC_1D_CalcSymmetryArea = Length * MidPoint * Pi * 2.
!   ! Area of the cells on the rotational symmetry axis is set to one
!   IF(.NOT.(DSMC_1D_CalcSymmetryArea.GT.0.0)) DSMC_1D_CalcSymmetryArea = 1.
! ELSE
  DSMC_1D_CalcSymmetryArea = Length
  IF (DSMC_1D_CalcSymmetryArea.EQ.0.) DSMC_1D_CalcSymmetryArea = 1.
! END IF
RETURN

END FUNCTION DSMC_1D_CalcSymmetryArea


FUNCTION DSMC_2D_CalcSymmetryAreaSubSides(iLocSide,iElem)
!===================================================================================================================================
!> Calculates the area of the subsides for the insertion with the surface flux
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: Pi
USE MOD_DSMC_Vars                 ,ONLY: ParticleWeighting
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iLocSide,iElem           !> iElem is the compute-node element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                              :: DSMC_2D_CalcSymmetryAreaSubSides(ParticleWeighting%nSubSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iNode
REAL                              :: P(1:2,1:4), Pmin(2), Pmax(2), MidPoint, PminTemp, PmaxTemp, Length
!===================================================================================================================================

Pmin = HUGE(Pmin)
Pmax = -HUGE(Pmax)

DO iNode = 1,4
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1)
END DO

Pmax(1) = MAXVAL(P(1,:))
Pmax(2) = MAXVAL(P(2,:))
Pmin(1) = MINVAL(P(1,:))
Pmin(2) = MINVAL(P(2,:))
Length = SQRT((Pmax(1)-Pmin(1))**2 + (Pmax(2)-Pmin(2))**2)

DO iNode = 1, ParticleWeighting%nSubSides
  PminTemp = Pmin(2) + (Pmax(2) - Pmin(2))/ParticleWeighting%nSubSides*(iNode-1.)
  PmaxTemp = Pmin(2) + (Pmax(2) - Pmin(2))/ParticleWeighting%nSubSides*iNode
  MidPoint = (PmaxTemp+PminTemp) / 2.
  DSMC_2D_CalcSymmetryAreaSubSides(iNode) = Length/ParticleWeighting%nSubSides * MidPoint * Pi * 2.
END DO

RETURN

END FUNCTION DSMC_2D_CalcSymmetryAreaSubSides


END MODULE MOD_Particle_Mesh_Tools
