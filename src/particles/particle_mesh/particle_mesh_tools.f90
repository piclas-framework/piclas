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

PUBLIC :: ParticleInsideQuad3D, InitPEM_LocalElemID, InitPEM_CNElemID, GetGlobalNonUniqueSideID, GetSideBoundingBoxTria
PUBLIC :: GetMeshMinMax, IdentifyElemAndSideType, WeirdElementCheck, CalcParticleMeshMetrics, CalcXCL_NGeo
PUBLIC :: CalcBezierControlPoints, InitParticleGeometry, ComputePeriodicVec
!===================================================================================================================================
CONTAINS


!PPURE SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det_Out)
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
INTEGER,INTENT(IN)            :: ElemID
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
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
CNElemID = GetCNElemID(ElemID)
DO iLocSide = 1,nlocSides
  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
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
          SideIDMortar = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideIDMortar)
          ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
          ! no way to recover it during runtime
          IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)

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


PPURE SUBROUTINE ParticleInsideNbMortar(PartStateLoc,ElemID,InElementCheck)
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
INTEGER,INTENT(IN)            :: ElemID
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
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
CNElemID = GetCNElemID(ElemID)
DO iLocSide = 1,nlocSides
  SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
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


FUNCTION GetGlobalNonUniqueSideID(ElemID,localSideID)
!===================================================================================================================================
!> Determines the non-unique global side ID of the local side in global element ElemID
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared, SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID                              !< global element ID
INTEGER,INTENT(IN) :: localSideID                         !< local side id of an element (1:6)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER :: GetGlobalNonUniqueSideID
INTEGER :: iSide,firstSide,lastSide
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + 1
lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND, ElemID)

! Small mortar sides are added after
DO iSide = firstSide,lastSide
  IF (SideInfo_Shared(SIDE_LOCALID,iSide).EQ.localSideID) THEN
    GetGlobalNonUniqueSideID = iSide
    RETURN
  END IF
END DO

! We should never arrive here
GetGlobalNonUniqueSideID=-1
CALL ABORT(__STAMP__,'GlobalSideID not found for Elem',ElemID)
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
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeSide,nComputeNodeSides
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeNode,nComputeNodeNodes
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
    GEO%CNxmin   = MINVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNxmax   = MAXVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymin   = MINVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymax   = MAXVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmin   = MINVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmax   = MAXVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))

    ! global
    GEO%xminglob = MINVAL(BezierControlPoints3D(1,:,:,:))
    GEO%xmaxglob = MAXVAL(BezierControlPoints3D(1,:,:,:))
    GEO%yminglob = MINVAL(BezierControlPoints3D(2,:,:,:))
    GEO%ymaxglob = MAXVAL(BezierControlPoints3D(2,:,:,:))
    GEO%zminglob = MINVAL(BezierControlPoints3D(3,:,:,:))
    GEO%zmaxglob = MAXVAL(BezierControlPoints3D(3,:,:,:))
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
    GEO%CNxmin   = MINVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNxmax   = MAXVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymin   = MINVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymax   = MAXVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmin   = MINVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmax   = MAXVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))

    ! global
    GEO%xminglob = MINVAL(NodeCoords_Shared(1,:))
    GEO%xmaxglob = MAXVAL(NodeCoords_Shared(1,:))
    GEO%yminglob = MINVAL(NodeCoords_Shared(2,:))
    GEO%ymaxglob = MAXVAL(NodeCoords_Shared(2,:))
    GEO%zminglob = MINVAL(NodeCoords_Shared(3,:))
    GEO%zmaxglob = MAXVAL(NodeCoords_Shared(3,:))
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
INTEGER                                  :: iElem,firstElem,lastElem,ElemID
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
LBWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved ...'

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

DO iElem=firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)

  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID)
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
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo),ElemCurved(iElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
    CNSideID = GetCNSideID(SideID)
    flip = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    IF(.NOT.ElemCurved(iElem))THEN
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
        v2=v1-ElemBaryNGeo(:,iElem)
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
          v2=v1-ElemBaryNGeo(:,iElem)
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
          v2=v1-ElemBaryNGeo(:,iElem)
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
END DO ! iElem=1,nTotalElems

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


DO iElem = firstElem,lastElem
  ElemID = GetCNElemID(iElem)
  IF (ElemCurved(ElemID)) THEN
    nCurvedElems = nCurvedElems+1
  ELSE
    nLinearElems = nLinearElems+1
  END IF

  DO ilocSide = 1,6
    ! ignore small mortar sides attached to big mortar sides
    SideID   = GetGlobalNonUniqueSideID(iElem,ilocSide)
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
INTEGER                        :: iElem
#endif /*USE_MPI*/
REAL                           :: Vdm_NGeo_CLNGeo(0:NGeo,0:NGeo)
!===================================================================================================================================

! small wBaryCL_NGEO
ALLOCATE( wBaryCL_NGeo1(            0:1)    ,&
    XiCL_NGeo1(               0:1)    ,&
    Vdm_CLNGeo1_CLNGeo(0:NGeo,0:1)    ,&
    ! new for curved particle sides
Vdm_Bezier(        0:NGeo,0:NGeo) ,&
    sVdm_Bezier(       0:NGeo,0:NGeo) ,&
    D_Bezier(          0:NGeo,0:NGeo))

CALL ChebyGaussLobNodesAndWeights(1,XiCL_NGeo1)
CALL BarycentricWeights(1,XiCL_NGeo1,wBaryCL_NGeo1)
CALL InitializeVandermonde(1, NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)
! initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm(NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier) !CHANGETAG
CALL BuildBezierDMat(NGeo,Xi_NGeo,D_Bezier)
CALL InitializeVandermonde(NGeo,NGeo,wBaryCL_NGeo,Xi_NGeo,XiCL_NGeo,Vdm_NGeo_CLNGeo)

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared((/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
XCL_NGeo_Shared (1:3    ,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => XCL_NGeo_Array

DO iElem = 1,nElems
  XCL_NGeo_Shared (:  ,:,:,:,offsetElem+iElem) = XCL_NGeo (:  ,:,:,:,iElem)
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
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP
USE MOD_Mesh_Vars              ,ONLY: dXCL_NGeo
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: NGeo,nGlobalElems,offsetElem,nElems
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
INTEGER                        :: iElem
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared((/3*  (PP_N+1)*(PP_N+1)*(PP_N+1)*nGlobalElems/), Elem_xGP_Shared_Win,Elem_xGP_Array)
CALL Allocate_Shared((/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
Elem_xGP_Shared (1:3    ,0:PP_N,0:PP_N,0:PP_N,1:nGlobalElems) => Elem_xGP_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => dXCL_NGeo_Array

DO iElem = 1,nElems
  Elem_xGP_Shared (:  ,:,:,:,offsetElem+iElem) = Elem_xGP (:  ,:,:,:,iElem)
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
                     , 3*(PP_N+1)**3*recvcountElem   &
                     , 3*(PP_N+1)**3*displsElem      &
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
Elem_xGP_Shared  => Elem_xGP
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
INTEGER                        :: SideID,ElemID,NbElemID,localSideID,localSideNbID,nStart,NodeMap(4,6)
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

DO ElemID = firstElem,lastElem
  ! Every periodic vector already found
  IF (ALL(PeriodicFound(:))) EXIT

SideLoop: DO SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)
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
      MasterCoords  = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)  +NodeMap(iNode                  ,localSideID))
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


END MODULE MOD_Particle_Mesh_Tools
