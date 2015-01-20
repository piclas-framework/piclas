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
USE MOD_Particle_Vars,         ONLY : PEM, PartState
USE nr,                        ONLY : gaussj 
USE MOD_part_MPFtools,         ONLY : MapToGeo, GeoCoordToMap
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
  REAL                          :: Vector(3,1)
  REAL                          :: NewNodePos(3,8)
  REAL                          :: ChosenMeanBaseD1, ChosenMeanBaseD2, ChosenMeanBaseD3
  REAL                          :: vLAG1, vLAG2, vLAG3
  REAL                          :: NVec1(3),NVec2(3),NVec3(3), xi_Out(3)
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
    CALL gaussj(Matrix,Vector)
    NewNodePos(1,iNode) = Vector(1,1)
    NewNodePos(2,iNode) = Vector(2,1)
    NewNodePos(3,iNode) = Vector(3,1)
  END DO

  nPart = PEM%pNumber(iElem)
  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    CALL GeoCoordToMap(PartState(iPartIndx,1:3),xi_Out,iElem)
    PartNewPos = MapToGeo(xi_Out,NewNodePos)

    LD_RHS(iPartIndx,1) = (PartNewPos(1) - PartState(iPartIndx,1)) / dt - PartState(iPartIndx,4)
    LD_RHS(iPartIndx,2) = (PartNewPos(2) - PartState(iPartIndx,2)) / dt - PartState(iPartIndx,5)
    LD_RHS(iPartIndx,3) = (PartNewPos(3) - PartState(iPartIndx,3)) / dt - PartState(iPartIndx,6)

    IF (LD_RHS(iPartIndx,1).NE.LD_RHS(iPartIndx,1)) THEN
      LD_RHS(iPartIndx,1) = 0.0
    END IF
    IF (LD_RHS(iPartIndx,2).NE.LD_RHS(iPartIndx,2)) THEN
      LD_RHS(iPartIndx,2) = 0.0
    END IF
    IF (LD_RHS(iPartIndx,3).NE.LD_RHS(iPartIndx,3)) THEN
      LD_RHS(iPartIndx,3) = 0.0
    END IF


!    IF ((ABS(PartNewPos(1) - PartState(iPartIndx,1)).LE. 1E-14) .OR. (dt.LE. 1E-14)) THEN
!      LD_RHS(iPartIndx,1) = 0.0
!    ELSE
!      LD_RHS(iPartIndx,1) = (PartNewPos(1) - PartState(iPartIndx,1)) / dt - PartState(iPartIndx,4)
!    END IF
!    IF ((ABS(PartNewPos(2) - PartState(iPartIndx,2)).LE. 1E-14) .OR. (dt.LE. 1E-14)) THEN
!      LD_RHS(iPartIndx,2) = 0.0
!    ELSE
!      LD_RHS(iPartIndx,2) = (PartNewPos(2) - PartState(iPartIndx,2)) / dt - PartState(iPartIndx,5)
!    END IF
!    IF ((ABS(PartNewPos(3) - PartState(iPartIndx,3)).LE. 1E-14) .OR. (dt.LE. 1E-14)) THEN
!      LD_RHS(iPartIndx,3) = 0.0
!    ELSE
!      LD_RHS(iPartIndx,3) = (PartNewPos(3) - PartState(iPartIndx,3)) / dt - PartState(iPartIndx,6)
!    END IF
#if (PP_TimeDiscMethod==1001)
    LD_DSMC_RHS(iPartIndx,1) = LD_RHS(iPartIndx,1)
    LD_DSMC_RHS(iPartIndx,2) = LD_RHS(iPartIndx,2)
    LD_DSMC_RHS(iPartIndx,3) = LD_RHS(iPartIndx,3)
#endif

    iPartIndx = PEM%pNext(iPartIndx)
  END DO
 
END SUBROUTINE LDPartTreament

!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_part_treat
