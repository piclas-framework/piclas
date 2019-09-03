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
! module for determination of particle push velocity
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LDPartTreament
  MODULE PROCEDURE LDPartTreament
END INTERFACE

PUBLIC :: LDPartTreament
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE LDPartTreament(iElem)
!===================================================================================================================================
! determination of particle push velocity
!===================================================================================================================================
! MODULES
USE MOD_LD_Vars
USE MOD_TimeDisc_Vars,         ONLY : dt
USE MOD_Particle_Vars,         ONLY : PEM, PartState,PartPosRef
USE MOD_Eval_xyz,              ONLY : TensorProductInterpolation
USE MOD_Mesh_Vars,             ONLY : NGeo
USE MOD_Basis,                 ONLY:GetInverse
USE MOD_Mesh_Vars,             ONLY: wBaryCL_NGeo,XiCL_NGeo
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER                       :: iNode, iLocSide1, iLocSide2, iLocSide3
  INTEGER                       :: nPart, iPart, iPartIndx
  REAL                          :: PartNewPos(3)
  REAL                          :: Matrix(3,3)
  REAL                          :: MatrixInv(3,3)
  REAL                          :: Vector(3,1)
  REAL                          :: NewNodePos(3,8),XCL_NGEO_tmp(1:3,0:NGeo,0:NGeo,0:NGeo)
  REAL                          :: ChosenMeanBaseD1, ChosenMeanBaseD2, ChosenMeanBaseD3
  REAL                          :: vLAG1, vLAG2, vLAG3
  REAL                          :: NVec1(3),NVec2(3),NVec3(3)!, xi_Out(3)
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
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
    !print*,'Node',iNode
    !print*,'iloside,vlag,nvec',ilocside1,vLag1,nVec1
    !print*,'iloside,vlag,nvec',ilocside2,vLag2,nVec2
    !print*,'iloside,vlag,nvec',ilocside3,vLag3,nVec3
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
    !print*,'matrix'
    !print*,matrix(1,:)
    !print*,matrix(2,:)
    !print*,matrix(3,:)
    MatrixInv=getInverse(3,Matrix)
    Vector(:,1)=MATMUL(MatrixInv,Vector(:,1))
    !CALL gaussj(Matrix,Vector)
    !print*,'vector',vector(:,1)
    !read*
    NewNodePos(1,iNode) = Vector(1,1)
    NewNodePos(2,iNode) = Vector(2,1)
    NewNodePos(3,iNode) = Vector(3,1)
  END DO
  ! map from vector in CGNS-Standard into polynomial order
  ! z-minus || zeta_minus
  XCL_NGeo_tmp(1:3,0   ,0   ,0)=NewNodePos(1:3,1)
  XCL_NGeo_tmp(1:3,NGeo,0   ,0)=NewNodePos(1:3,2)
  XCL_NGeo_tmp(1:3,NGeo,NGeo,0)=NewNodePos(1:3,3)
  XCL_NGeo_tmp(1:3,0   ,NGeo,0)=NewNodePos(1:3,4)
  ! z-plus || zeta_plus
  XCL_NGeo_tmp(1:3,0   ,0   ,NGeo)=NewNodePos(1:3,5)
  XCL_NGeo_tmp(1:3,NGeo,0   ,NGeo)=NewNodePos(1:3,6)
  XCL_NGeo_tmp(1:3,NGeo,NGeo,NGeo)=NewNodePos(1:3,7)
  XCL_NGeo_tmp(1:3,0   ,NGeo,NGeo)=NewNodePos(1:3,8)

  nPart = PEM%pNumber(iElem)
  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    CALL TensorProductInterpolation(PartPosRef(1:3,iPartIndx),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo_tmp,PartNewPos)

    LD_RHS(iPartIndx,1) = (PartNewPos(1) - PartState(1,iPartIndx)) / dt - PartState(4,iPartIndx)
    LD_RHS(iPartIndx,2) = (PartNewPos(2) - PartState(2,iPartIndx)) / dt - PartState(5,iPartIndx)
    LD_RHS(iPartIndx,3) = (PartNewPos(3) - PartState(3,iPartIndx)) / dt - PartState(6,iPartIndx)

    IF (LD_RHS(iPartIndx,1).NE.LD_RHS(iPartIndx,1)) THEN
      LD_RHS(iPartIndx,1) = 0.0
    END IF
    IF (LD_RHS(iPartIndx,2).NE.LD_RHS(iPartIndx,2)) THEN
      LD_RHS(iPartIndx,2) = 0.0
    END IF
    IF (LD_RHS(iPartIndx,3).NE.LD_RHS(iPartIndx,3)) THEN
      LD_RHS(iPartIndx,3) = 0.0
    END IF
#if (PP_TimeDiscMethod==1001)
    LD_DSMC_RHS(1,iPartIndx) = LD_RHS(iPartIndx,1)
    LD_DSMC_RHS(2,iPartIndx) = LD_RHS(iPartIndx,2)
    LD_DSMC_RHS(3,iPartIndx) = LD_RHS(iPartIndx,3)
#endif

    iPartIndx = PEM%pNext(iPartIndx)
  END DO

END SUBROUTINE LDPartTreament

!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_part_treat
