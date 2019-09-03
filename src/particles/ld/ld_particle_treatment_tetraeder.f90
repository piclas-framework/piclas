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
MODULE MOD_LD_part_treat

!===================================================================================================================================
! module for determination of Lagrangian cell velocity
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC :: LDPartTreament
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE LDPartTreament(iElem)

USE MOD_LD_Vars
USE MOD_TimeDisc_Vars,         ONLY : dt
USE MOD_Particle_Vars,         ONLY : PEM, PartState!,GEO
USE MOD_Basis,                 ONLY : GetInverse
!USE MOD_part_MPFtools,         ONLY : MapToGeo, GeoCoordToMap

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER                       :: iNode, iLocSide1, iLocSide2, iLocSide3, iLocSide, iVertex
  INTEGER                       :: nPart, iPart, iPartIndx, trinum
  REAL                          :: PartShift(3)
  REAL                          :: Matrix(3,3)
  REAL                          :: MatrixInv(3,3)
  REAL                          :: Vector(3,1)
  REAL                          :: NewNodePos(3,8)
  REAL                          :: ChosenMeanBaseD1, ChosenMeanBaseD2, ChosenMeanBaseD3
  REAL                          :: vLAG1, vLAG2, vLAG3
  REAL                          :: NVec1(3),NVec2(3),NVec3(3), xi_Out(3),NewHexCentroid(3),OldHexCentroid(3)
  LOGICAL, ALLOCATABLE        :: PartFound(:)
  REAL                          :: Wert_min, Wert_sum
  REAL                          :: DistVertex, DistPart

  INTEGER, INTENT(IN)          :: iElem
!--------------------------------------------------------------------------------------------------!

  DO iNode = 1, 8
    ! ----which three Sides belongs to Node (cgns stards)
    IF(iNode.eq.1) THEN
      iLocSide1 = 1
      iLocSide2 = 2
      iLocSide3 = 5
    ELSE IF(iNode.eq.2) THEN
      iLocSide1 = 1
      iLocSide2 = 2
      iLocSide3 = 3
    ELSE IF(iNode.eq.3) THEN
      iLocSide1 = 1
      iLocSide2 = 3
      iLocSide3 = 4
    ELSE IF(iNode.eq.4) THEN
      iLocSide1 = 1
      iLocSide2 = 4
      iLocSide3 = 5
    ELSE IF(iNode.eq.5) THEN
      iLocSide1 = 2
      iLocSide2 = 5
      iLocSide3 = 6
    ELSE IF(iNode.eq.6) THEN
      iLocSide1 = 2
      iLocSide2 = 3
      iLocSide3 = 6
    ELSE IF(iNode.eq.7) THEN
      iLocSide1 = 3
      iLocSide2 = 4
      iLocSide3 = 6
    ELSE
      iLocSide1 = 4
      iLocSide2 = 5
      iLocSide3 = 6
    END IF
    NVec1(1:3) = MeanSurfValues(iLocSide1,iElem)%MeanNormVec(1:3)
    vLAG1 = MeanSurfValues(iLocSide1,iElem)%MeanLagVelo
    NVec2(1:3) = MeanSurfValues(iLocSide2,iElem)%MeanNormVec(1:3)
    vLAG2 = MeanSurfValues(iLocSide2,iElem)%MeanLagVelo
    NVec3(1:3) = MeanSurfValues(iLocSide3,iElem)%MeanNormVec(1:3)
    vLAG3 = MeanSurfValues(iLocSide3,iElem)%MeanLagVelo
    ChosenMeanBaseD1 = MeanSurfValues(iLocSide1,iElem)%MeanBaseD
    ChosenMeanBaseD2 = MeanSurfValues(iLocSide2,iElem)%MeanBaseD
    ChosenMeanBaseD3 = MeanSurfValues(iLocSide3,iElem)%MeanBaseD
    Matrix(1,1) = NVec1(1)
    Matrix(2,1) = NVec2(1)
    Matrix(3,1) = NVec3(1)
! ----
    Matrix(1,2) = NVec1(2)
    Matrix(2,2) = NVec2(2)
    Matrix(3,2) = NVec3(2)
! ----
    Matrix(1,3) = NVec1(3)
    Matrix(2,3) = NVec2(3)
    Matrix(3,3) = NVec3(3)
! ----
    Vector(1,1) = ChosenMeanBaseD1 &
                + vLAG1 * dt

    Vector(2,1) = ChosenMeanBaseD2 &
                + vLAG2 * dt

    Vector(3,1) = ChosenMeanBaseD3 &
                + vLAG3 * dt
    MatrixInv=getInverse(3,Matrix)
    Vector(:,1)=MATMUL(MatrixInv,Vector(:,1))
    !CALL gaussj(Matrix,Vector)
    NewNodePos(1,iNode) = Vector(1,1)
    NewNodePos(2,iNode) = Vector(2,1)
    NewNodePos(3,iNode) = Vector(3,1)
    ! PO old
    !NewNodePosIndx(1,GEO%ElemToNodeID(iNode,iElem)) = Vector(1,1)
    !NewNodePosIndx(2,GEO%ElemToNodeID(iNode,iElem)) = Vector(2,1)
    !NewNodePosIndx(3,GEO%ElemToNodeID(iNode,iElem)) = Vector(3,1)
    !OldNodePos(1,iNode) = GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,iElem))
    !OldNodePos(2,iNode) = GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,iElem))
    !OldNodePos(3,iNode) = GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,iElem))
  END DO
  xi_Out = 0.0
!  NewHexCentroid = MapToGeo(xi_Out,NewNodePos)
!  xi_Out = 0.0
!  OldHexCentroid = MapToGeo(xi_Out,OldNodePos)

  nPart = PEM%pNumber(iElem)
  ALLOCATE(PartFound(nPart))
  PartFound(1:nPart)=.FALSE.

  DO iLocSide = 1, 6
    DO trinum = 1, 2
!      IF (trinum.EQ.1) THEN
!        TetraPoint(1) = GEO%ElemSideNodeID(1,iLocSide,iElem)
!        TetraPoint(2) = GEO%ElemSideNodeID(2,iLocSide,iElem)
!        TetraPoint(3) = GEO%ElemSideNodeID(3,iLocSide,iElem)
!      ELSE
!        TetraPoint(1) = GEO%ElemSideNodeID(1,iLocSide,iElem)
!        TetraPoint(2) = GEO%ElemSideNodeID(4,iLocSide,iElem)
!        TetraPoint(3) = GEO%ElemSideNodeID(3,iLocSide,iElem)
!      END IF
      iPartIndx = PEM%pStart(iElem)
      DO ipart = 1, nPart
        IF (.NOT.PartFound(ipart)) THEN
          ! Localisation
          !Matrix(1,1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(1))
          !Matrix(2,1) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(1))
          !Matrix(3,1) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(1))

          !Matrix(1,2) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(2))
          !Matrix(2,2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(2))
          !Matrix(3,2) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(2))

          !Matrix(1,3) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(3))
          !Matrix(2,3) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(3))
          !Matrix(3,3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(3))

          Vector(1,1) = - OldHexCentroid(1) + PartState(1,iPartIndx)
          Vector(2,1) = - OldHexCentroid(2) + PartState(2,iPartIndx)
          Vector(3,1) = - OldHexCentroid(3) + PartState(3,iPartIndx)

          !CALL gaussj(Matrix,Vector)
          MatrixInv=getInverse(3,Matrix)
          Vector(:,1)=MATMUL(MatrixInv,Vector(:,1))


          Wert_min = MIN(Vector(1,1),Vector(2,1),Vector(3,1)) + 1e-14 !!!
          Wert_sum = Vector(1,1) + Vector(2,1) + Vector(3,1)  - 1e-14 !!!
          IF ((Wert_min .GE. 0.0).and.(Wert_sum .LE. 1.0)) THEN
            PartFound(ipart)=.TRUE.
            PartShift(1:3) = 0.0
            DO iVertex = 1, 4
              IF (iVertex.EQ.4) THEN
                !CALL CalcDistOpposite(iVertex,OldHexCentroid,DistVertex,OldHexCentroid,TetraPoint)
                !CALL CalcDistOpposite(iVertex,OldHexCentroid,DistPart,PartState(1:3,iPartIndx),TetraPoint)
                PartShift(1) = PartShift(1) + DistPart / DistVertex * &
                             (- OldHexCentroid(1) + NewHexCentroid(1))
                PartShift(2) = PartShift(2) + DistPart / DistVertex * &
                             (- OldHexCentroid(2) + NewHexCentroid(2))
                PartShift(3) = PartShift(3) + DistPart / DistVertex * &
                             (- OldHexCentroid(3) + NewHexCentroid(3))
              ELSE
                !CALL CalcDistOpposite(iVertex,OldHexCentroid,DistVertex,GEO%NodeCoords(1:3,TetraPoint(iVertex)),TetraPoint)
                !CALL CalcDistOpposite(iVertex,OldHexCentroid,DistPart,PartState(1:3,iPartIndx),TetraPoint)
                !PartShift(1) = PartShift(1) + DistPart / DistVertex * &
                !             (- GEO%NodeCoords(1,TetraPoint(iVertex)) &
                !              + NewNodePosIndx(1,TetraPoint(iVertex)))
                !PartShift(2) = PartShift(2) + DistPart / DistVertex * &
                !             (- GEO%NodeCoords(2,TetraPoint(iVertex)) &
                !              + NewNodePosIndx(2,TetraPoint(iVertex)))
                !PartShift(3) = PartShift(3) + DistPart / DistVertex * &
                !             (- GEO%NodeCoords(3,TetraPoint(iVertex)) &
                !              + NewNodePosIndx(3,TetraPoint(iVertex)))
              END IF
            END DO ! iVertex
            LD_RHS(iPartIndx,1) = PartShift(1) / dt - PartState(4,iPartIndx)
            LD_RHS(iPartIndx,2) = PartShift(2) / dt - PartState(5,iPartIndx)
            LD_RHS(iPartIndx,3) = PartShift(3) / dt - PartState(6,iPartIndx)
            IF (ABS(PartShift(1)).LE. 1E-14) THEN
              LD_RHS(iPartIndx,1) = 0.0
            END IF
            IF (ABS(PartShift(2)).LE. 1E-14) THEN
              LD_RHS(iPartIndx,2) = 0.0
            END IF
            IF (ABS(PartShift(3)).LE. 1E-14) THEN
              LD_RHS(iPartIndx,3) = 0.0
            END IF
          END IF ! PartLocalisation
        END IF ! PartFound
        iPartIndx = PEM%pNext(iPartIndx)
      END DO ! END iPart
    END DO ! END trinum
  END DO ! END iLocSide

  DEALLOCATE(PartFound)

END SUBROUTINE LDPartTreament
!    !--------------------------------------------------------------------------------------------------!
!    !--------------------------------------------------------------------------------------------------!
!    SUBROUTINE CalcDistOpposite(TetraVertex,OldHexCentroid,DistOpposite,PointCoords,TetraPoint)
!    !--------------------------------------------------------------------------------------------------!
!      USE MOD_LD_Vars
!      !USE MOD_Particle_Vars,         ONLY : GEO
!      IMPLICIT NONE                                                                                   !
!    !--------------------------------------------------------------------------------------------------!
!      INTEGER, INTENT(IN)                :: TetraVertex
!      REAL, INTENT(IN)                   :: OldHexCentroid(3),PointCoords(3)
!      INTEGER, INTENT(IN)                :: TetraPoint(3)
!      REAL, INTENT(OUT)                  :: DistOpposite
!
!      REAL                                :: Fak_a, Fak_b, Fak_c
!    !--------------------------------------------------------------------------------------------------!
!    !  IF (TetraVertex .EQ. 1) THEN
!    !    Vec1(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(2))
!    !    Vec1(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(2))
!    !    Vec1(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(2))
!    !
!    !    Vec2(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(3))
!    !    Vec2(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(3))
!    !    Vec2(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(3))
!    !
!    !    BaseVec(1:3) = OldHexCentroid(1:3)
!    !  ELSE IF (TetraVertex .EQ. 2) THEN
!    !    Vec1(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(1))
!    !    Vec1(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(1))
!    !    Vec1(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(1))
!    !
!    !    Vec2(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(3))
!    !    Vec2(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(3))
!    !    Vec2(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(3))
!    !
!    !    BaseVec(1:3) = OldHexCentroid(1:3)
!    !  ELSE IF (TetraVertex .EQ. 3) THEN
!    !    Vec1(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(1))
!    !    Vec1(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(1))
!    !    Vec1(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(1))
!    !
!    !    Vec2(1) = - OldHexCentroid(1) + GEO%NodeCoords(1,TetraPoint(2))
!    !    Vec2(2) = - OldHexCentroid(2) + GEO%NodeCoords(2,TetraPoint(2))
!    !    Vec2(3) = - OldHexCentroid(3) + GEO%NodeCoords(3,TetraPoint(2))
!    !
!    !    BaseVec(1:3) = OldHexCentroid(1:3)
!    !  ELSE ! TetraVertex .EQ. OldHexCentroid
!    !    Vec1(1) = - GEO%NodeCoords(1,TetraPoint(1)) + GEO%NodeCoords(1,TetraPoint(2))
!    !    Vec1(2) = - GEO%NodeCoords(2,TetraPoint(1)) + GEO%NodeCoords(2,TetraPoint(2))
!    !    Vec1(3) = - GEO%NodeCoords(3,TetraPoint(1)) + GEO%NodeCoords(3,TetraPoint(2))
!    !
!    !    Vec2(1) = - GEO%NodeCoords(1,TetraPoint(1)) + GEO%NodeCoords(1,TetraPoint(3))
!    !    Vec2(2) = - GEO%NodeCoords(2,TetraPoint(1)) + GEO%NodeCoords(2,TetraPoint(3))
!    !    Vec2(3) = - GEO%NodeCoords(3,TetraPoint(1)) + GEO%NodeCoords(3,TetraPoint(3))
!    !
!    !    BaseVec(1:3) = GEO%NodeCoords(1:3,TetraPoint(1))
!    !  END IF
!
!      ! -> Koordinatenform
!    !  Fak_a = Vec1(2)*Vec2(3) - Vec2(2)*Vec1(3)
!    !  Fak_b = Vec1(3)*Vec2(1) - Vec2(3)*Vec1(1)
!    !  Fak_c = Vec1(1)*Vec2(2) - Vec2(1)*Vec1(2)
!    !  DistOpposite = ABS((Fak_a*(PointCoords(1)-BaseVec(1))  &
!    !             + Fak_b*(PointCoords(2)-BaseVec(2))  &
!    !             + Fak_c*(PointCoords(3)-BaseVec(3))) &
!    !             / SQRT(Fak_a*Fak_a + Fak_b*Fak_b + Fak_c*Fak_c))
!
!    END SUBROUTINE CalcDistOpposite

!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_part_treat
