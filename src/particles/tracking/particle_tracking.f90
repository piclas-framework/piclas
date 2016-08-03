#include "boltzplatz.h"

MODULE MOD_Particle_Tracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

!INTERFACE ParticleTracking
!  MODULE PROCEDURE ParticleTracking
!END INTERFACE

INTERFACE ParticleTrackingCurved
  MODULE PROCEDURE ParticleTrackingCurved
END INTERFACE

INTERFACE ParticleRefTracking
  MODULE PROCEDURE ParticleRefTrackingFast
  !MODULE PROCEDURE ParticleRefTrackingSLOW
END INTERFACE

!PUBLIC::ParticleTracking,ParticleTrackingCurved
PUBLIC::ParticleTrackingCurved
PUBLIC::ParticleRefTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTrackingCurved(doParticle_In)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Tracking_vars,      ONLY:ntracks,nCurrentParts
USE MOD_Particle_Mesh,               ONLY:SingleParticleToExactElementNoMap
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,PartInElemCheck
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
#ifdef MPI
USE MOD_LoadBalance_Vars,            ONLY:ElemTime
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: doParticle(1:PDM%ParticleVecLength)
INTEGER                       :: iPart,ElemID,flip,OldElemID
INTEGER                       :: ilocSide,SideID, locSideList(1:6), hitlocSide,nInterSections,nLoc
LOGICAL                       :: PartisDone,dolocSide(1:6),isHit,markTol
REAL                          :: localpha(1:6),xi(1:6),eta(1:6)
!INTEGER                       :: lastlocSide
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory,xNodes(1:3,1:4)
#ifdef MPI
REAL                          :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

DO iPart=1,PDM%ParticleVecLength
  IF(DoParticle(iPart))THEN
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    nTracks=nTracks+1
    nCurrentParts=nCurrentParts+1
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )
    IF(ALMOSTZERO(lengthPartTrajectory))THEN
      PEM%Element(iPart)=ElemID
      PartisDone=.TRUE.
      CYCLE
    END IF
    PartTrajectory=PartTrajectory/lengthPartTrajectory!+epsilontol
    !PartTrajectory=PartTrajectory/(lengthPartTrajectory+epsilontol)
    lengthPartTrajectory=lengthPartTrajectory!+epsilontol
    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    DO WHILE (.NOT.PartisDone)
      locAlpha=-1.
      nInterSections=0
      markTol=.FALSE.
      DO ilocSide=1,6
        locSideList(ilocSide)=ilocSide
        IF(.NOT.dolocSide(ilocSide)) CYCLE
        !SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
        flip  = PartElemToSide(E2S_FLIP,ilocSide,ElemID)
        SELECT CASE(SideType(SideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                        ,xi (ilocSide)      &
                                                                                        ,eta(ilocSide)   ,iPart,flip,SideID)

!          CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                        ,xi (ilocSide)      &
!                                                                                        ,eta(ilocSide)   ,iPart,flip,SideID)
!                                                                              !,doTest=.TRUE.)
!                                                                                  !,eta(ilocSide)   ,iPart,ilocSide,SideID,ElemID)

        CASE(BILINEAR,PLANAR_NONRECT)
          xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,SideID)
          xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,SideID)
          xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,SideID)
          xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,SideID)
          !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,xNodes &
          !                                                    ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
          !                                                                                  ,xi (ilocSide)      &
          !                                                                                  ,eta(ilocSide)      &
          !                                                                                  ,iPart,flip,SideID)
          CALL ComputeBiLinearIntersectionRobust(isHit,xNodes &
                                                ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                            ,xi (ilocSide)      &
                                                                                            ,eta(ilocSide)      &
                                                                                            ,iPart,flip,SideID)

          !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,[BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
          !                                                    ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
          !                                                    ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
          !                                                    ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
          !                                                    ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
          !                                                                                  ,xi (ilocSide)      &
          !                                                                                  ,eta(ilocSide)      &
          !                                                                                  ,iPart,flip,SideID)
!                                                                                            !ilocSide,iPart,SideID)
!          CALL ComputeBezierIntersection(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                            ,xi (ilocSide)      &
!                                                                            ,eta(ilocSide)      ,iPart,SideID)
!

        CASE(CURVED)
          CALL ComputeBezierIntersection(ishit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)      &
                                                                                  ,eta(ilocSide)      ,iPart,SideID)
        CASE DEFAULT
          CALL abort(&
          __STAMP__ &
          ,' Missing required side-data. Please increase halo region. ',SideID)
        END SELECT
        IF(isHit) THEN
          nInterSections=nInterSections+1
          IF((ABS(xi(ilocSide)).GE.1.0).OR.(ABS(eta(ilocSide)).GE.1.0)) markTol=.TRUE.
          IF(ALMOSTZERO(locAlpha(ilocSide))) markTol=.TRUE.
          !IF((ABS(xi(ilocSide)).GE.0.99).OR.(ABS(eta(ilocSide)).GE.0.99)) markTol=.TRUE.
          !IF(locAlpha(ilocSide)/lengthPartTrajectory.GE.0.99) markTol=.TRUE.
        END IF
      END DO ! ilocSide
      !IF((ipart.eq.40).AND.(iter.GE.68)) print*,' nIntersections',nIntersections
      SELECT CASE(nInterSections)
      CASE(0) ! no intersection
        PEM%Element(iPart)=ElemID
        PartisDone=.TRUE.
      CASE(1) ! one intersection
        ! get intersection side
        DO ilocSide=1,6
          IF(locAlpha(ilocSide).GT.-1.0) THEN
            hitlocSide=ilocSide
            EXIT
          END IF
        END DO ! ilocSide
        SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
        OldElemID=ElemID
        CALL SelectInterSectionType(PartIsDone,doLocSide,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                         ,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,ElemID)
        IF(ElemID.NE.OldElemID)THEN
#ifdef MPI
          tLBEnd = LOCALTIME() ! LB Time End
          ElemTime(OldELemID)=ElemTime(OldElemID)+tLBEnd-tLBStart
          tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
        END IF
      CASE DEFAULT ! two or more hits
        ! take last possible intersection, furthest
        !CALL BubbleSortID(locAlpha,locSideList,6)
        CALL InsertionSort(locAlpha,locSideList,6)
!        IF((ipart.eq.40).AND.(iter.GE.68)) THEN
!        END IF
        nloc=0
        DO ilocSide=6,1,-1
          IF(locAlpha(ilocSide).GT.-1.0)THEN
            !nloc=nloc+1
            hitlocSide=locSideList(ilocSide)
            EXIT
            !IF(nloc.EQ.nInterSections) EXIT
          END IF
        END DO ! ilocSide
        SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
        OldElemID=ElemID
        CALL SelectInterSectionType(PartIsDone,doLocSide,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                         ,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,ElemID)
        IF(ElemID.NE.OldElemID)THEN
#ifdef MPI
          tLBEnd = LOCALTIME() ! LB Time End
          ElemTime(OldELemID)=ElemTime(OldElemID)+tLBEnd-tLBStart
          tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
        END IF
        markTol=.TRUE.
      END SELECT
      IF(markTol)THEN
        CALL PartInElemCheck(iPart,ElemID,isHit)
        !print*,'partid',ipart,'in elem',isHit
        !print*,'old elem',ElemID
        !print*,'ipart,loc',ipart,localpha
        PEM%Element(iPart)=ElemID
        IF(.NOT.isHit) CALL SingleParticleToExactElementNoMap(iPart,doHALO=.TRUE.)!debug=.TRUE.)
!        print*,'new elem',PEM%Element(ipart)
        PartIsDone=.TRUE.
        IF(.NOT.PDM%ParticleInside(iPart))THEN
          IPWRITE(UNIT_stdOut,*) 'lost particle with id', ipart
        END IF
      END IF ! markTol
    END DO ! PartisDone=.FALSE.
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    IF(PEM%Element(iPart).LE.PP_nElems) ElemTime(PEM%Element(iPart))=ElemTime(PEM%Element(iPart))+tLBEnd-tLBStart
#endif /*MPI*/
!    IF(markTol)THEN
!      CALL PartInElemCheck(iPart,ElemID,isHit)
!      PEM%Element(iPart)=ElemID
!      IF(.NOT.isHit) CALL SingleParticleToExactElementNoMap(iPart,debug=.TRUE.)
!      PartIsDone=.TRUE.
!      IF(.NOT.PDM%ParticleInside(iPart))THEN
!        IPWRITE(UNIT_stdOut,*) 'lost particle with id', ipart
!      END IF
!    END IF
  END IF ! Part inside
END DO ! iPart

END SUBROUTINE ParticleTrackingCurved


SUBROUTINE ParticleBCTrackingfast(ElemID,firstSide,LastSide,nlocSides,PartId,PartisDone,PartisMoved)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Mesh_Vars,          ONLY:BCElem,GEO
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,ComputePlanarIntersectionBezierRobust2
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
LOGICAL,INTENT(INOUT)         :: PartisDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES!
LOGICAL,INTENT(INOUT)         :: PartisMoved
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip,BCSideID
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory,xNodes(1:3,1:4)
LOGICAL                       :: DoTracing,PeriMoved
!===================================================================================================================================


PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )

IF(ALMOSTZERO(lengthPartTrajectory))THEN
  PEM%Element(PartID)=ElemID
  PartIsDone=.TRUE.
  RETURN
END IF

PartTrajectory=PartTrajectory/lengthPartTrajectory

PartisMoved=.FALSE.
DoTracing=.TRUE.
DO WHILE(DoTracing)
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! call here function for mapping of partpos and lastpartpos
    CALL PeriodicMovement(PartID,PeriMoved)
  ELSE
    PeriMoved=.FALSE.
  END IF
  locAlpha=-1.0
  nInter=0
  !nlocSides=lastSide-firstSide+1
  DO iLocSide=firstSide,LastSide
    ! track particle vector until the final particle position is achieved
    SideID=BCElem(ElemID)%BCSideID(ilocSide)
    BCSideID=PartBCSideList(SideID)
    locSideList(ilocSide)=ilocSide
    ! get correct flip
    flip  = 0 
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT)

      !CALL ComputePlanarIntersectionBezier(ishit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
      !                                                                              ,xi (ilocSide)            &
      !                                                                              ,eta(ilocSide)   ,PartID,flip,BCSideID)

      !CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
      CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                    ,xi (ilocSide)            &
                                                                                    ,eta(ilocSide)   ,PartID,flip,BCSideID)
    CASE(BILINEAR,PLANAR_NONRECT)
      xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,BCSideID)
      xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,BCSideID)
      xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,BCSideID)
      xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,BCSideID)
      CALL ComputeBiLinearIntersectionRobust(isHit,xNodes &
                                                   ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                       ,xi (ilocSide)      &
                                                                                       ,eta(ilocSide)      &
                                                                                       ,PartID,flip,BCSideID)
    CASE(CURVED)
      CALL ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                              ,xi (ilocSide)      &
                                                                              ,eta(ilocSide)      ,PartID,BCSideID)
    END SELECT
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      nInter=nInter+1
    END IF
  END DO ! ilocSide
  
  IF(nInter.EQ.0)THEN
    IF(.NOT.PeriMoved) DoTracing=.FALSE. 
  ELSE
    ! take first possible intersection
    !CALL BubbleSortID(locAlpha,locSideList,6)
    PartIsMoved=.TRUE.
    CALL InsertionSort(locAlpha,locSideList,nlocSides)
    ilocSide=LastSide-nInter+1
    hitlocSide=locSideList(ilocSide)
    SideID=BCElem(ElemID)%BCSideID(hitlocSide)
    BCSideID=PartBCSideList(SideID)
    CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                      ,xi(hitlocSide)     &
                                                                      ,eta(hitlocSide)    &
                                                                      ,PartId,SideID,ElemID)
    IF(.NOT.PDM%ParticleInside(PartID)) THEN
      PartisDone = .TRUE.
       RETURN
     END IF
  END IF ! nInter>0
END DO

END SUBROUTINE ParticleBCTrackingfast


SUBROUTINE ParticleBCTracking(ElemID,firstSide,LastSide,nlocSides,PartId,PartisDone)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Mesh_Vars,          ONLY:BCElem
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,ComputePlanarIntersectionBezierRobust2
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
LOGICAL,INTENT(INOUT)         :: PartisDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: dolocSide(firstSide:lastSide),ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip,BCSideID
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory,xNodes(1:3,1:4)
!===================================================================================================================================

IF(LastSide.EQ.0) RETURN

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )

IF(ALMOSTZERO(lengthPartTrajectory))THEN
  PEM%Element(PartID)=ElemID
  PartIsDone=.TRUE.
  RETURN
END IF

PartTrajectory=PartTrajectory/lengthPartTrajectory
!lengthPartTrajectory=lengthPartTrajectory!+epsilontol

locAlpha=-1.0
nInter=0
dolocSide=.TRUE.
!nlocSides=lastSide-firstSide+1
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID=BCElem(ElemID)%BCSideID(ilocSide)
  BCSideID=PartBCSideList(SideID)
  locSideList(ilocSide)=ilocSide
  ! get correct flip
  flip  = 0 
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT)
    !CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
    CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)            &
                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)

!    CALL ComputePlanarIntersectionBezierRobust2(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                  ,xi (ilocSide)      &
!                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)
!
!                                                                            !,eta(ilocSide)   ,PartID,ilocSide,SideID,ElemID)
  CASE(BILINEAR,PLANAR_NONRECT)
    xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,BCSideID)
    xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,BCSideID)
    xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,BCSideID)
    xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,BCSideID)
    !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,xNodes &
    CALL ComputeBiLinearIntersectionRobust(isHit,xNodes &
                                                 ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                     ,xi (ilocSide)      &
                                                                                     ,eta(ilocSide)      &
                                                                                     ,PartID,flip,BCSideID)

    !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,[BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
    !                                                    ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
    !                                                                                        ,xi (ilocSide)      &
    !                                                                                        ,eta(ilocSide)  ,PartID,flip,SideID)
!    CALL ComputeBezierIntersection(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                      ,xi (ilocSide)      &
!                                                                      ,eta(ilocSide)      ,PartID,SideID)

  CASE(CURVED)
    CALL ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                            ,xi (ilocSide)      &
                                                                            ,eta(ilocSide)      ,PartID,BCSideID)
  END SELECT
  IF(locAlpha(ilocSide).GT.-1.0)THEN
    nInter=nInter+1
  END IF
END DO ! ilocSide

IF(nInter.EQ.0) THEN
  RETURN
ELSE
  ! take first possible intersection
  !CALL BubbleSortID(locAlpha,locSideList,6)
  CALL InsertionSort(locAlpha,locSideList,nlocSides)
  ilocSide=LastSide-nInter+1
  hitlocSide=locSideList(ilocSide)
  SideID=BCElem(ElemID)%BCSideID(hitlocSide)
  BCSideID=PartBCSideList(SideID)
  CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                    ,xi(hitlocSide)     &
                                                                    ,eta(hitlocSide)    &
                                                                    ,PartId,SideID,ElemID)
  IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  RETURN
  !DO ilocSide=firstSide,LastSide
  !DO ilocSide=LastSide
  !  IF(locAlpha(ilocSide).GT.-1.0)THEN
  !    hitlocSide=locSideList(ilocSide)
  !    !SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
  !    SideID=BCElem(ElemID)%BCSideID(hitlocSide)
  !    BCSideID=PartBCSideList(SideID)
  !    CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
  !                                                                      ,xi(hitlocSide)     &
  !                                                                      ,eta(hitlocSide)    &
  !                                                                      ,PartId,SideID,ElemID)
  !    IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  !    RETURN
  !  END IF ! locAlpha>-1.0
  !END DO ! ilocSide
END IF ! nInter>0

END SUBROUTINE ParticleBCTracking


SUBROUTINE ParticleRefTrackingSlow(doParticle_In)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Particle_Vars,           ONLY:PDM,PEM,PartState,PartPosRef,LastpartPOs
USE MOD_Mesh_Vars,               ONLY:OffSetElem
USE MOD_Eval_xyz,                ONLY:eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars,  ONLY:nTracks,nCurrentParts
USE MOD_Particle_Mesh_Vars,      ONLY:Geo,IsBCElem,BCElem,epsInCell,epsOneCell
USE MOD_Utils,                   ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Mesh_Vars,      ONLY:ElemBaryNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElement
USE MOD_Eval_xyz,                ONLY:Eval_XYZ_Poly
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloElemToProc
USE MOD_LoadBalance_Vars,        ONLY:ElemTime,nTracksPerElem,tTracking
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: doParticle(1:PDM%ParticleVecLength)
INTEGER                          :: iPart, ElemID,oldElemID,iElem, newElemID
INTEGER                          :: CellX,CellY,CellZ,iBGMElem,nBGMElems
REAL,ALLOCATABLE                 :: Distance(:)
REAL                             :: oldXi(3),newXi(3), LastPos(3),epsLowOne
INTEGER,ALLOCATABLE              :: ListDistance(:)
!REAL                             :: epsOne
#ifdef MPI
INTEGER                          :: InElem
#endif
INTEGER                          :: TestElem
LOGICAL                          :: ParticleFound(1:PDM%ParticleVecLength),PartisDone
!LOGICAL                          :: HitBC(1:PDM%ParticleVecLength)
! load balance
#ifdef MPI
REAL                               :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

ParticleFound=.FALSE.
!HitBC=.FALSE.
epsLowOne=1.0-2.0*epsInCell
! first step, reuse Elem cache, therefore, check if particle are still in element, if not, search later
DO iElem=1,PP_nElems ! loop only over internal elems, if particle is already in HALO, it shall not be found here
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  DO iPart=1,PDM%ParticleVecLength
    IF(DoParticle(iPart))THEN
      ElemID = PEM%lastElement(iPart)
      IF(ElemID.NE.iElem) CYCLE
      nTracks=nTracks+1
      ! sanity check
      !IF(PartState(iPart,3).GE.0.089)THEN
      !  IPWRITE(UNIT_stdOut,*) ' Part out of area, z,ipart', PartState(iPart,3),iPart
      !END IF
      IF(IsBCElem(ElemID))THEN
#if defined(LSERK)
        CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,iPart,ParticleFound(iPart))
        IF(ParticleFound(iPart)) CYCLE
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
        IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle is inside 
        !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle is inside 
          PEM%lastElement(iPart)=ElemID
          ParticleFound(iPart)=.TRUE.
          CYCLE
        END IF
#else
        ! simple and stupid
        LastPos=PartState(iPart,1:3)
        IF(GEO%nPeriodicVectors.GT.0)THEN
          ! call here function for mapping of partpos and lastpartpos
          CALL PeriodicMovement(iPart)
        END IF
        CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,iPart,ParticleFound(iPart))
        IF(ParticleFound(iPart)) CYCLE
        DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
          LastPos=PartState(iPart,1:3)
          IF(GEO%nPeriodicVectors.GT.0)THEN
            ! call here function for mapping of partpos and lastpartpos
            CALL PeriodicMovement(iPart)
          END IF
          ! unfortunately, here all sides
          CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%lastSide &
              ,BCElem(ElemID)%lastSide,iPart,ParticleFound(iPart))
          IF(ParticleFound(iPart)) EXIT
        END DO
        IF(ParticleFound(iPart)) CYCLE
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
        IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle is inside 
        !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle is inside 
          ParticleFound(iPart)=.TRUE.
          CYCLE
        END IF

#endif /*TIMEDISCS*/
      ELSE ! no bc elem, therefore, no bc ineraction possible
        IF(GEO%nPeriodicVectors.GT.0)THEN
          ! call here function for mapping of partpos and lastpartpos
          LastPos=PartState(iPart,1:3)
          CALL PeriodicMovement(iPart)
          IF(.NOT.IsBCElem(ElemID))THEN
            DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
                .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
                .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
              LastPos=PartState(iPart,1:3)
              ! call here function for mapping of partpos and lastpartpos
              CALL PeriodicMovement(iPart)
            END DO
          END IF
        END IF
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
#endif
        IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle inside
        !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
          PEM%Element(iPart)  = ElemID
          ParticleFound(iPart)=.TRUE.
        !  IPWRITE(UNIT_stdOut,*) ' partposref to large!',iPart
        END IF
      END IF ! initial check
    ELSE
      ! caution: dummy, because particle is not inside and such a particle sould not be searched for in 
      ! the next loop
      ParticleFound(iPart)=.TRUE.
    END IF
  END DO ! iPart
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
END DO ! iElem

! now, locate not all found particle
#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
DO iPart=1,PDM%ParticleVecLength
  IF(ParticleFound(iPart)) CYCLE
  ! relocate particle
  oldElemID = PEM%lastElement(iPart) ! this is not!  a possible elem
  ! get background mesh cell of particle
  CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
  CellX = MAX(MIN(GEO%FIBGMimax,CellX),GEO%FIBGMimin)
  CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
  CellY = MAX(MIN(GEO%FIBGMjmax,CellY),GEO%FIBGMjmin)
  CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
  CellZ = MAX(MIN(GEO%FIBGMkmax,CellZ),GEO%FIBGMkmin)
        
  ! check all cells associated with this beckground mesh cell
  nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem
  !SDEALLOCATE( Distance)
  !SDEALLOCATE( ListDistance)
  ALLOCATE( Distance(1:nBGMElems) &
          , ListDistance(1:nBGMElems) )
 
  ! get closest element barycenter
  Distance=0.
  ListDistance=0
  DO iBGMElem = 1, nBGMElems
    ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    ListDistance(iBGMElem)=ElemID
    IF(ElemID.EQ.-1)CYCLE
    IF(ElemID.EQ.OldElemID)THEN
      Distance(iBGMElem)=-1.0
    ELSE
      Distance(iBGMElem)=    ((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
                             +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                             +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )

      IF(Distance(iBGMElem).GT.ElemRadius2NGeo(ElemID))THEN
        Distance(iBGMElem)=-1.0
      END IF
    END IF
  END DO ! nBGMElems

  !CALL BubbleSortID(Distance,ListDistance,nBGMElems)
  CALL InsertionSort(Distance,ListDistance,nBGMElems)

  OldXi=PartPosRef(1:3,iPart)
  newXi=HUGE(1.0)
  newElemID=-1
  ! loop through sorted list and start by closest element  
  DO iBGMElem=1,nBGMElems
    IF(ALMOSTEQUAL(Distance(iBGMELem),-1.0)) CYCLE
    ElemID=ListDistance(iBGMElem)
#ifdef MPI
    IF(ElemID.LE.PP_nElems) nTracksPerElem(ElemID)=nTracksPerElem(ElemID)+1
#endif /*MPI*/
    CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.BezierClipHit) THEN ! particle inside
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.epsOneCell) THEN ! particle inside
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
      PEM%Element(iPart) = ElemID
      ParticleFound(iPart)=.TRUE.
      EXIT
    END IF
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
      newXi=PartPosRef(1:3,iPart)
      newElemID=ElemID
    END IF
  END DO ! iBGMElem

  IF(.NOT.ParticleFound(iPart))THEN
    ! use best xi
    !IPWRITE(UNIT_stdOut,*) ' recover particle', iPart
    IF(MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi)))THEN
      PartPosRef(1:3,iPart)=OldXi
      PEM%Element(iPart)   =oldElemID
      ElemID               =oldElemID
    ELSE
      PartPosRef(1:3,iPart)=NewXi
      PEM%Element(iPart)   =NewElemID
      oldElemID            =NewElemID
      ElemID               =NewElemID
    END IF
  
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.05) THEN
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.0) THEN
      PartIsDone=.FALSE.
      TestElem=ElemID
      IF(.NOT.IsBCElem(TestElem))THEN
        ! ausgabe
        IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with internal element '
        IPWRITE(UNIT_stdOut,*) ' xi          ', PartPosRef(1:3,iPart)
        IPWRITE(UNIT_stdOut,*) ' epsOneCell  ', epsOneCell
        IPWRITE(UNIT_stdOut,*) ' oldxi       ', oldXi
        IPWRITE(UNIT_stdOut,*) ' newxi       ', newXi
        IPWRITE(UNIT_stdOut,*) ' ParticlePos ', PartState(iPart,1:3)
#ifdef MPI
        InElem=PEM%Element(iPart)
        IF(InElem.LE.PP_nElems)THEN
          IPWRITE(UNIT_stdOut,*) ' ElemID       ', InElem+offSetElem
        ELSE
          IPWRITE(UNIT_stdOut,*) ' ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                 + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
        END IF
#else
        IPWRITE(UNIT_stdOut,*) ' ElemID       ', PEM%Element(iPart)+offSetElem
#endif
        CALL abort(&
        __STAMP__ &
        ,'Particle Not inSide of Element, iPart',iPart)
      ELSE ! BCElem
        !CALL ComputeFaceIntersection(TestElem,1,BCElem(TestElem)%nInnerSides,BCElem(TestElem)%nInnerSides,iPart,PartIsDone)
        CALL ComputeFaceIntersection(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone)
        LastPos=PartState(iPart,1:3)
        CALL ParticleBCTracking(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone)
        ! check inner sides
        DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
          LastPos=PartState(iPart,1:3)
          ! unfortunately, here all sides
          CALL ParticleBCTracking(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone)
          IF(PartIsDone) EXIT
          IF(GEO%nPeriodicVectors.GT.0)THEN
            ! call here function for mapping of partpos and lastpartpos
            CALL PeriodicMovement(iPart)
          END IF
        END DO ! While
        IF(PartIsDone) THEN
          DEALLOCATE( Distance)
          DEALLOCATE( ListDistance)
          CYCLE
        END IF
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),TestElem)
        ! false, reallocate particle
        !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.0)THEN
        IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsOneCell)THEN
          IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with BC element, relocating!! '
          CALL SingleParticleToExactElement(iPart,doHalo=.TRUE.)                                                             
          IF(.NOT.PDM%ParticleInside(iPart)) THEN
            IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with BC element '
            IPWRITE(UNIT_stdOut,*) ' xi          ', partposref(1:3,ipart)
            IPWRITE(UNIT_stdOut,*) ' epsonecell  ', epsonecell
            IPWRITE(UNIT_stdOut,*) ' oldxi       ', oldxi
            IPWRITE(UNIT_stdOut,*) ' newxi       ', newxi
            IPWRITE(UNIT_stdOut,*) ' particlepos ', partstate(ipart,1:3)
            IPWRITE(UNIT_stdOut,*) ' velocity    ', partstate(ipart,4:6)
            IPWRITE(UNIT_stdOut,*) ' lastpartpos ', LastPartPos(ipart,1:3)
#ifdef MPI
            inelem=PEM%Element(ipart)
            IF(inelem.LE.PP_nElems)THEN
              IPWRITE(UNIT_stdout,*) ' elemid       ', inelem+offsetelem
            ELSE
              IPWRITE(UNIT_stdOut,*) ' elemid       ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
                                                       + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
            END IF
#else
            IPWRITE(UNIt_stdOut,*) ' elemid       ', pem%element(ipart)+offsetelem
#endif
CALL abort(&
__STAMP__ &
,' Particle not inside of element, ipart',ipart)
          END IF ! inside
        END IF ! epsCell
      END IF ! BCElem
    END IF ! inner eps to large
  END IF ! not found
  ParticleFound(iPart)=.TRUE.
  DEALLOCATE( Distance)
  DEALLOCATE( ListDistance)

END DO ! iPart
#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tTracking = tTracking +tLBEnd-tLBStart
#endif /*MPI*/

END SUBROUTINE ParticleRefTrackingSlow


SUBROUTINE ParticleRefTrackingFast(doParticle_In)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Particle_Vars,           ONLY:PDM,PEM,PartState,PartPosRef, PartSpecies
USE MOD_Mesh_Vars,               ONLY:OffSetElem
USE MOD_Eval_xyz,                ONLY:eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars,  ONLY:nTracks,nCurrentParts
USE MOD_Particle_Mesh_Vars,      ONLY:Geo,IsBCElem,BCElem,epsInCell,epsOneCell
USE MOD_Utils,                   ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Mesh_Vars,      ONLY:ElemBaryNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElement
USE MOD_Eval_xyz,                ONLY:Eval_XYZ_Poly
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloElemToProc
USE MOD_LoadBalance_Vars,        ONLY:ElemTime,nTracksPerElem,tTracking
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: doParticle(1:PDM%ParticleVecLength)
INTEGER                           :: iPart, ElemID,oldElemID,newElemID
INTEGER                           :: CellX,CellY,CellZ,iBGMElem,nBGMElems
REAL,ALLOCATABLE                  :: Distance(:)
REAL                              :: oldXi(3),newXi(3), LastPos(3)
INTEGER,ALLOCATABLE               :: ListDistance(:)
!REAL                              :: epsOne
#ifdef MPI
INTEGER                           :: InElem
#endif
INTEGER                           :: TestElem
!LOGICAL                           :: ParticleFound(1:PDM%ParticleVecLength),PartisDone
LOGICAL                           :: PartisDone,PartIsMoved
!LOGICAL                           :: HitBC(1:PDM%ParticleVecLength)
! load balance
#ifdef MPI
REAL                                :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

DO iPart=1,PDM%ParticleVecLength
  IF(DoParticle(iPart))THEN
    ElemID = PEM%lastElement(iPart)
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    nTracks=nTracks+1
    ! sanity check
    PartIsDone=.FALSE.
    IF(IsBCElem(ElemID))THEN
      CALL ParticleBCTrackingfast(ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,iPart,PartIsDone,PartIsMoved)
      IF(PartIsDone) CYCLE
      IF(PartIsMoved)THEN
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
      ELSE
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* only LSERK */
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
#endif
      END IF
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle is inside 
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle is inside 
        PEM%lastElement(iPart)=ElemID
#ifdef MPI
         tLBEnd = LOCALTIME() ! LB Time End
         ElemTime(ElemID)=ElemTime(ElemID)+tLBEnd-tLBStart
#endif /*MPI*/
        CYCLE
      END IF
    ELSE ! no bc elem, therefore, no bc ineraction possible
      IF(GEO%nPeriodicVectors.GT.0)THEN
        ! call here function for mapping of partpos and lastpartpos
        LastPos=PartState(iPart,1:3)
        CALL PeriodicMovement(iPart)
        IF(.NOT.IsBCElem(ElemID))THEN
          DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
              .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
            LastPos=PartState(iPart,1:3)
            ! call here function for mapping of partpos and lastpartpos
            CALL PeriodicMovement(iPart)
          END DO
        END IF
      END IF
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
#else
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
#endif
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
        PEM%Element(iPart)  = ElemID
#ifdef MPI
         tLBEnd = LOCALTIME() ! LB Time End
         ElemTime(ElemID)=ElemTime(ElemID)+tLBEnd-tLBStart
#endif /*MPI*/
        CYCLE
      !ELSE IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.5) THEN
      !  IPWRITE(UNIT_stdOut,*) ' partposref to large!',iPart
      END IF
    END IF ! initial check
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(ElemID)=ElemTime(ElemID)+tLBEnd-tLBStart
#endif /*MPI*/
  ! still not located
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    ! relocate particle
    oldElemID = PEM%lastElement(iPart) ! this is not!  a possible elem
    ! get background mesh cell of particle
    CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
    CellX = MAX(MIN(GEO%FIBGMimax,CellX),GEO%FIBGMimin)
    CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
    CellY = MAX(MIN(GEO%FIBGMjmax,CellY),GEO%FIBGMjmin)
    CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
    CellZ = MAX(MIN(GEO%FIBGMkmax,CellZ),GEO%FIBGMkmin)
          
    ! check all cells associated with this beckground mesh cell
    nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem
    !SDEALLOCATE( Distance)
    !SDEALLOCATE( ListDistance)
    ALLOCATE( Distance(1:nBGMElems) &
            , ListDistance(1:nBGMElems) )
 
    ! get closest element barycenter
    Distance=0.
    ListDistance=0
    DO iBGMElem = 1, nBGMElems
      ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
      ListDistance(iBGMElem)=ElemID
      IF(ElemID.EQ.-1)CYCLE
      IF(ElemID.EQ.OldElemID)THEN
        Distance(iBGMElem)=-1.0
      ELSE
        !Distance(iBGMElem)=SQRT((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))  &
        !                       +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
        !                       +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )
        Distance(iBGMElem)=    ((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
                               +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                               +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )

        IF(Distance(iBGMElem).GT.ElemRadius2NGeo(ElemID))THEN
          Distance(iBGMElem)=-1.0
        END IF
      END IF
    END DO ! nBGMElems

    !CALL BubbleSortID(Distance,ListDistance,nBGMElems)
    CALL InsertionSort(Distance,ListDistance,nBGMElems)

    OldXi=PartPosRef(1:3,iPart)
    newXi=HUGE(1.0)
    newElemID=-1
    ! loop through sorted list and start by closest element  
    DO iBGMElem=1,nBGMElems
      IF(ALMOSTEQUAL(Distance(iBGMELem),-1.0)) CYCLE
      ElemID=ListDistance(iBGMElem)
#ifdef MPI
      IF(ElemID.LE.PP_nElems) nTracksPerElem(ElemID)=nTracksPerElem(ElemID)+1
#endif /*MPI*/
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN ! particle inside
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.epsOneCell) THEN ! particle inside
        PEM%Element(iPart) = ElemID
        PartIsDone=.TRUE.
        EXIT
      END IF
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
        newXi=PartPosRef(1:3,iPart)
        newElemID=ElemID
      END IF
    END DO ! iBGMElem
    DEALLOCATE( Distance     &
              , ListDistance )
    IF(.NOT.PartIsDone)THEN
      ! use best xi
      IF(MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi)))THEN
        PartPosRef(1:3,iPart)=OldXi
        PEM%Element(iPart)=oldElemID
      ELSE
        PartPosRef(1:3,iPart)=NewXi
        PEM%Element(iPart)=NewElemID
        oldElemID=NewElemID
      END IF
    
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsOneCell) THEN
        PartIsDone=.FALSE.
        TestElem=PEM%Element(iPart)
        IF(.NOT.IsBCElem(TestElem))THEN
          ! ausgabe
          IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with internal element '
          IPWRITE(UNIT_stdOut,*) ' xi          ', PartPosRef(1:3,iPart)
          IPWRITE(UNIT_stdOut,*) ' epsOneCell  ', epsOneCell
          IPWRITE(UNIT_stdOut,*) ' oldxi       ', oldXi
          IPWRITE(UNIT_stdOut,*) ' newxi       ', newXi
          IPWRITE(UNIT_stdOut,*) ' ParticlePos ', PartState(iPart,1:3)
#ifdef MPI
          InElem=PEM%Element(iPart)
          IF(InElem.LE.PP_nElems)THEN
            IPWRITE(UNIT_stdOut,*) ' ElemID       ', InElem+offSetElem
          ELSE
            IPWRITE(UNIT_stdOut,*) ' ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                                   + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
          END IF
#else
          IPWRITE(UNIT_stdOut,*) ' ElemID       ', PEM%Element(iPart)+offSetElem
#endif
CALL abort(&
__STAMP__ &
,'Particle Not inSide of Element, iPart',iPart)
        ELSE ! BCElem
          !CALL ComputeFaceIntersection(ElemID,1,BCElem(ElemID)%nInnerSides,BCElem(ElemID)%nInnerSides,iPart,PartIsDone)
          CALL ComputeFaceIntersection(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone)
          LastPos=PartState(iPart,1:3)
          CALL ParticleBCTrackingfast(TestElem,1,BCElem(TestElem)%lastSide,BCElem(TestElem)%lastSide,iPart,PartIsDone,PartIsMoved)
          IF(PartIsDone) CYCLE
          CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),TestElem)
          ! false, reallocate particle
          IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsOneCell)THEN
            IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with BC element, relocating!! '
            CALL SingleParticleToExactElement(iPart,doHalo=.TRUE.)                                                             
            IF(.NOT.PDM%ParticleInside(iPart)) THEN
              IPWRITE(UNIT_stdOut,*) ' Tolerance Issue with BC element '
              IPWRITE(UNIT_stdOut,*) ' xi          ', partposref(1:3,ipart)
              IPWRITE(UNIT_stdOut,*) ' epsonecell  ', epsonecell
              IPWRITE(UNIT_stdOut,*) ' oldxi       ', oldxi
              IPWRITE(UNIT_stdOut,*) ' newxi       ', newxi
              IPWRITE(UNIT_stdOut,*) ' particlepos ', partstate(ipart,1:3)
#ifdef MPI
              inelem=PEM%Element(ipart)
              IF(inelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,*) ' elemid       ', inelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdOut,*) ' elemid       ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,inelem)) &
                                                         + PartHaloElemToProc(NATIVE_ELEM_ID,inelem)
              END IF
              IF(testelem.LE.PP_nElems)THEN
                IPWRITE(UNIT_stdout,*) ' testelem       ', testelem+offsetelem
              ELSE
                IPWRITE(UNIT_stdOut,*) ' testelem       ', offsetelemmpi(PartHaloElemToProc(NATIVE_PROC_ID,testelem)) &
                                                         + PartHaloElemToProc(NATIVE_ELEM_ID,testelem)
              END IF

#else
              IPWRITE(UNIt_stdOut,*) ' elemid       ', pem%element(ipart)+offsetelem
#endif
              CALL abort(&
    __STAMP__ &
    ,'particle noT inside of element, ipart',ipart)
            END IF ! inside
          ELSE
            PEM%Element(iPart)=TestElem
          END IF ! epsCell
        END IF ! BCElem
      END IF ! inner eps to large
    END IF
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    tTracking = tTracking +tLBEnd-tLBStart
#endif /*MPI*/
  END IF
END DO ! iPart

END SUBROUTINE ParticleRefTrackingFast


SUBROUTINE SelectInterSectionType(PartIsDone,doLocSide,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                 ,xi,eta,alpha,PartID,SideID,ElemID)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Mesh_Vars,                   ONLY:nBCSides
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToElem
USE MOD_Particle_Mesh_Vars,          ONLY:SidePeriodicType, SidePeriodicDisplacement
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos,PDM
#ifdef MPI
USE MOD_Mesh_Vars,                   ONLY:BC,nSides
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: PartID,SideID,hitlocSide,ilocSide
REAL,INTENT(INOUT)                :: Xi,Eta,Alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: PartIsDone
LOGICAL,INTENT(INOUT)             :: DoLocSide(1:6)
INTEGER,INTENT(INOUT)             :: ElemID
REAL,INTENT(INOUT),DIMENSION(1:3) :: PartTrajectory
REAL,INTENT(INOUT)                :: lengthPartTrajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: flip
!===================================================================================================================================
IF(SideID.LE.nBCSides)THEN
  ! check if interesction is possible and take first intersection
  CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha &
                                                                 ,xi    &
                                                                 ,eta   ,PartID,SideID,ElemID)

  IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  dolocSide=.TRUE.
  dolocSide(hitlocSide)=.FALSE.
ELSE ! no BC Side
  ! check for periodic sides
  IF(SidePeriodicType(SideID).GT.0)THEN
    flip   = PartElemToSide(E2S_FLIP,hitlocSide,ElemID)
    ! update particle position of intersection
    IF(flip.EQ.0)THEN
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+SidePeriodicDisplacement(:,SidePeriodicType(SideID))
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    ELSE
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-SidePeriodicDisplacement(:,SidePeriodicType(SideID))
      PartState(PartID,1:3)  =PartState(PartID,1:3)  -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    END IF
    ! recompute particle trajectory
    PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )
    PartTrajectory=PartTrajectory/lengthPartTrajectory
    lengthPartTrajectory=lengthPartTrajectory
    ! update particle element
    dolocSide=.TRUE.
    dolocSide(PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,ElemID))=.FALSE.
    ElemID=PartElemToElem(E2E_NB_ELEM_ID,ilocSide,ElemID)
    !lastlocSide=-1
  ELSE ! no periodic side
#ifdef MPI
    IF(SideID.GT.nSides)THEN
      IF(BC(SideID).NE.0)THEN
        ! encountered a bc side
        CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,Alpha &
                                                                       ,xi    &
                                                                       ,eta   ,PartID,SideID,ElemID)
        dolocSide=.TRUE.
        IF(SideType(SideID).EQ.PLANAR_RECT) THEN !also for PLANAR_NONRECT?
          dolocSide(hitlocSide)=.FALSE.
        END IF
      ELSE
        ! inner side
        dolocSide=.TRUE.
        dolocSide(PartElemToElem(E2E_NB_LOC_SIDE_ID,hitlocSide,ElemID))=.FALSE.
        ElemID=PartElemToElem(E2E_NB_ELEM_ID,hitlocSide,ElemID)
        IF(ElemID.LE.0) CALL abort(&
__STAMP__&
,' HaloRegion too small or critical error during halo region reconstruction!')
      END IF ! BC?
    ELSE
      dolocSide=.TRUE.
      dolocSide(PartElemToElem(E2E_NB_LOC_SIDE_ID,hitlocSide,ElemID))=.FALSE.
      ElemID=PartElemToElem(E2E_NB_ELEM_ID,hitlocSide,ElemID)
      !lastlocSide=-1
    END IF ! SideID.GT.nSides
#else
    dolocSide=.TRUE.
    dolocSide(PartElemToElem(E2E_NB_LOC_SIDE_ID,hitlocSide,ElemID))=.FALSE.
    ElemID=PartElemToElem(E2E_NB_ELEM_ID,hitlocSide,ElemID)
#endif /* MP!!I */
    END IF ! SidePeriodicType
END IF ! SideID>nCBSides

END SUBROUTINE SelectInterSectionType


SUBROUTINE PeriodicMovement(PartID,isMovedOut)
!----------------------------------------------------------------------------------------------------------------------------------!
! move particle in the periodic direction, if particle is outside of the box
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,          ONLY:GEO
USE MOD_Particle_Tracking_Vars,      ONLY:FastPeriodic
#ifdef MPI
USE MOD_Particle_MPI_Vars,           ONLY:PartShiftVector
#endif /*MPI*/
#ifdef IMEX
USE MOD_TimeDisc_Vars,               ONLY: dt,iStage
USE MOD_TimeDisc_Vars,               ONLY: ERK_a,RK_b,RK_c
USE MOD_Particle_Vars,               ONLY: PartStateN,PartStage
#endif /*IMEX*/
#ifdef IMPA
USE MOD_Particle_Vars,               ONLY: PartQ
USE MOD_TimeDisc_Vars,               ONLY: dt,iStage
USE MOD_TimeDisc_Vars,               ONLY: ERK_a,ESDIRK_a,RK_b,RK_c,RK_bs
USE MOD_Particle_Vars,               ONLY: PartStateN,PartStage
USE MOD_LinearSolver_Vars,           ONLY: PartXK
#endif /*IMPA*/
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
USE MOD_Particle_Vars,              ONLY: PartIsImplicit
#endif 
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)              :: PartID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL   :: isMovedOut
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPV
REAL                            :: MoveVector(1:3)
LOGICAL                         :: isMoved
#ifdef IMPA
INTEGER                         :: iCounter
REAL                            :: DeltaP(6),Norm
#endif /*IMPA*/
#ifdef IMEX
INTEGER                         :: iCounter
#endif /*IMEX*/
!===================================================================================================================================

#ifdef MPI
PartShiftVector(1:3,PartID)=PartState(PartID,1:3)
#endif /*MPI*/
isMoved=.FALSE.
IF(FastPeriodic)THEN
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,1)-GEO%xmaxglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,1)-GEO%xminglob)/ABS(GEO%PeriodicVectors(1,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,2)-GEO%ymaxglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,2)-GEO%yminglob)/ABS(GEO%PeriodicVectors(2,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,3)-GEO%zmaxglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      MoveVector=CEILING(ABS(PartState(PartID,3)-GEO%zminglob)/ABS(GEO%PeriodicVectors(3,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+MoveVector
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -MoveVector
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-MoveVector
        isMoved=.TRUE.
      END IF
    END IF
  END IF


  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x+, PartID',PartID)
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside x-, PartID',PartID)
    END IF
  END IF
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y+, PartID',PartID)
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside y-, PartID',PartID)
    END IF
  END IF
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z+, PartID',PartID)
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      IPWRITE(*,*) 'PartPos', PartState(PartID,:)
      CALL abort(&
      __STAMP__ &
      ,' particle outside z-, PartID',PartID)
    END IF
  END IF 
ELSE
  ! x direction
  IF(GEO%directions(1)) THEN
    IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,1).LT.GEO%xminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
      END DO
      IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  
  ! y direction
  IF(GEO%directions(2)) THEN
    IF(PartState(PartID,2).GT.GEO%ymaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,2).LT.GEO%yminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
      END DO
      IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF
  
  ! z direction
  IF(GEO%directions(3)) THEN
    IF(PartState(PartID,3).GT.GEO%zmaxglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
    IF(PartState(PartID,3).LT.GEO%zminglob) THEN
      DO iPV=1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
      END DO
      IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
        PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      ELSE
        PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
        LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
        isMoved=.TRUE.
      END IF
    END IF
  END IF
END IF

#ifdef MPI
PartShiftVector(1:3,PartID)=-PartState(PartID,1:3)+PartShiftvector(1:3,PartID)
#endif /*MPI*/

IF(isMoved)THEN
#ifdef IMEX 
  ! recompute PartStateN to kill jump in integration through periodic BC
  IF(iStage.GT.0)THEN
    PartStateN(PartID,1:6) = PartState(PartID,1:6)
    DO iCounter=1,iStage-1
      PartStateN(PartID,1:6) = PartStateN(PartID,1:6)   &
                        - ERK_a(iStage,iCounter)*dt*PartStage(PartID,1:6,iCounter)
    END DO
  END IF
#endif /*IMEX*/

#ifdef IMPA 
  ! recompute PartStateN to kill jump in integration through periodic BC
  IF(iStage.GT.0)THEN
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
    IF(PartIsImplicit(PartID))THEN
#endif
      ! partshift-vector is pointing from parallel-pos to old pos
      PartStateN(PartID,1:6) = PartState(PartID,1:6)
      ! explicit particle
      DeltaP=0.
      DO iCounter=1,iStage-1
        DeltaP=DeltaP - ESDIRK_a(iStage,iCounter)*dt*PartStage(PartID,1:6,iCounter)
      END DO
      PartStateN(PartID,1:6) = PartStateN(PartID,1:6) + DeltaP ! plus, because DeltaP is defined neg
      PartQ(1:3,PartID) = PartQ(1:3,PartID) - PartShiftVector(1:3,PartID)
      ! and move all the functions
      ! F_PartX0 is not changing, because of difference
      PartXK(1:3,PartID) = PartXK(1:3,PartID) - PartShiftVector(1:3,PartID)
      ! brainfuck, does F_PartXK(:,PartID) is changed?
      ! init: F_PartXK=F_PartXK0 
      ! sanity check for parallel...
      !! old norm
      !DeltaP(1)=PartState(PartID,4)
      !DeltaP(2)=PartState(PartID,5)
      !DeltaP(3)=PartState(PartID,6)
      !DeltaP(4)=Pt(PartID,1)
      !DeltaP(5)=Pt(PartID,2)
      !DeltaP(6)=Pt(PartID,3)
      !IF(iStage.GT.0)THEN
      !DeltaP=PartState(PartID,:) - PartQ(:,PartID)-ESDIRK_A(iStage,iStage)*dt*DeltaP
      !CALL PartVectorDotProduct(DeltaP,DeltaP,Norm)
      !ELSE
      !Norm=Norm2_F_PartXK(PartID)
      !END IF
      !IF(.NOT.ALMOSTEQUAL(Norm,Norm2_F_PartXK(PartID)))THEN
      !  print*,'norm diff',PartID
      !END IF
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
   ELSE
     PartStateN(PartID,1:6) = PartState(PartID,1:6)
     ! explicit particle
     DO iCounter=1,iStage-1
       PartStateN(PartID,1:6) = PartStateN(PartID,1:6)   &
                              - ERK_a(iStage,iCounter)*dt*PartStage(PartID,1:6,iCounter)
     END DO
   END IF
#endif
  END IF
#endif /*IMPA*/
END IF


IF(PRESENT(isMovedOut)) isMovedOut=isMoved

END SUBROUTINE PeriodicMovement


SUBROUTINE ReComputeParticleBCInteraction(xi,eta,locSideID,SideID,BCSideID,PartID) 
!----------------------------------------------------------------------------------------------------------------------------------!
! The particle BC intersection is ignored. therefore, the particle is mapped onto the BC at the lost position and a wrong
! particle BC interaction is performed. FALLBACK!!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Vars,           ONLY:PDM,PEM,PartState,PartPosRef,lastpartpos
USE MOD_Eval_xyz,                ONLY:Eval_XYZ_Poly
USE MOD_Mesh_Vars,               ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Particle_Mesh_Vars,      ONLY:Geo,epsOneCell,ElemRadius2NGeo,ElemBaryNGeo,BCElem
USE MOD_Particle_Mesh_Vars,      ONLY:GEO
USE MOD_Utils,                   ONLY:BubbleSortID,InsertionSort
USE MOD_Eval_xyz,                ONLY:eval_xyz_elemcheck
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloElemToProc
#endif
USE MOD_Mesh_Vars,               ONLY:OffSetElem
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloElemToProc
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)                    :: locSideID,PartID
INTEGER,INTENT(INOUT)                 :: SideID,BCSideID
REAL,INTENT(IN)                       :: xi,eta
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: ElemID,oldElemID,newElemID
REAL                        :: hit
REAL                        :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                     :: ParticleFound
INTEGER                     :: CellX,CellY,CellZ,iBGMElem,nBGMElems
REAL,ALLOCATABLE            :: Distance(:)
REAL                        :: oldXi(3),newXi(3), LastPos(3),OldestPos(3),blubb(3), n_loc(3)
INTEGER,ALLOCATABLE         :: ListDistance(:)
#ifdef MPI
INTEGER                     :: inElem
#endif /*MPI*/
!===================================================================================================================================

!SELECT CASE(locSideID)
!CASE(XI_MINUS)
!  !BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
!  Xi =PartPosRef(3,PartID)
!  Eta=PartPosRef(2,PartID)
!  PartPosRef(1,PartID)=-0.9999
!CASE(XI_PLUS)
!  Xi =PartPosRef(2,PartID)
!  Eta=PartPosRef(3,PartID)
!  PartPosRef(1,PartID)= 0.9999
!CASE(ETA_MINUS)
!  Xi =PartPosRef(1,PartID)
!  Eta=PartPosRef(3,PartID)
!  PartPosRef(2,PartID)=-0.9999
!CASE(ETA_PLUS)
!  !BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,NGeo-p,q)
!  ! hopefully correct
!  Xi =-PartPosRef(1,PartID)
!  Eta=PartPosRef(3,PartID)
!  PartPosRef(2,PartID)= 0.9999
!CASE(ZETA_MINUS)
!  Xi =PartPosRef(2,PartID)
!  Eta=PartPosRef(1,PartID)
!  PartPosRef(3,PartID)=-0.9999
!CASE(ZETA_PLUS)
!  Xi =PartPosRef(1,PartID)
!  Eta=PartPosRef(2,PartID)
!  PartPosRef(3,PartID)= 0.9999
!END SELECT

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory

SELECT CASE(SideType(BCSideID))
CASE(PLANAR_RECT,PLANAR_NONRECT)
  n_loc=SideNormVec(1:3,BCSideID)
CASE(BILINEAR)
  CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
CASE(CURVED)
  CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
!   CALL abort(__STAMP__'nvec for bezier not implemented!',999,999.)
END SELECT 

IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.)THEN
  ElemID=PEM%Element(PartID)
  CALL Eval_xyz_Poly(PartPosRef(:,PartID),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),PartState(PartID,1:3))
  RETURN
ELSE
  ! compute temporary last particle position
  ElemID=PEM%Element(PartID)
  oldestPos=PartState(PartID,1:3)
  CALL Eval_xyz_Poly(PartPosRef(:,PartID),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),LastPartPos(PartID,1:3))
  blubb=LastPartPos(PartID,1:3)
  PartState(PartID,1:3)=LastPartPos(PartID,1:3)+PartTrajectory*lengthPartTrajectory
END IF

LastPos=PartState(PartID,1:3)
hit=0.
CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,hit,Xi,Eta,PartID,SideID,ElemID)

! if boundary condition is an open boundary condition, particle is deleted
IF(.NOT.PDM%ParticleInside(PartID)) RETURN

IF(GEO%nPeriodicVectors.GT.0)THEN
  CALL PeriodicMovement(PartID)
END IF

DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(PartID,1)) &
    .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(PartID,2)) &
    .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(PartID,3)) )
  LastPos=PartState(PartID,1:3)
  ! unfortunately, here all sides
  CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,PartID,ParticleFound)
  IF(ParticleFound) RETURN
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! call here function for mapping of partpos and lastpartpos
    CALL PeriodicMovement(PartID)
  END IF
END DO

CALL Eval_xyz_ElemCheck(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)

! if particle is found
IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).LE.1.0) THEN ! particle inside
  PEM%Element(PartID)  = ElemID
  RETURN
ELSE
  ! particle has to to located in its final cell
  oldElemID=ElemID
  PDM%ParticleInside(PartID)=.FALSE.
  ! BGM cell
  CellX = CEILING((PartState(PartID,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
  !CellX = MIN(GEO%FIBGMimax,CellX)
  CellX = MAX(MIN(GEO%FIBGMimax,CellX),GEO%FIBGMimin)
  CellY = CEILING((PartState(PartID,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
  !CellY = MIN(GEO%FIBGMjmax,CellY)
  CellY = MAX(MIN(GEO%FIBGMjmax,CellY),GEO%FIBGMjmin)
  CellZ = CEILING((PartState(PartID,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
  CellZ = MAX(MIN(GEO%FIBGMkmax,CellZ),GEO%FIBGMkmin)
  !CellZ = MIN(GEO%FIBGMkmax,CellZ)

  ! check all cells associated with this beckground mesh cell
  nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem
  ALLOCATE( Distance(1:nBGMElems) &
          , ListDistance(1:nBGMElems) )

  ! get closest element barycenter
  Distance=0.
  ListDistance=0
  DO iBGMElem = 1, nBGMElems
    ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    ListDistance(iBGMElem)=ElemID
    IF(ElemID.EQ.-1)CYCLE
    IF(ElemID.EQ.OldElemID)THEN
      Distance(iBGMElem)=-1.0
    ELSE
      Distance(iBGMElem)=    ((PartState(PartID,1)-ElemBaryNGeo(1,ElemID))*(PartState(PartID,1)-ElemBaryNGeo(1,ElemID)) &
                             +(PartState(PartID,2)-ElemBaryNGeo(2,ElemID))*(PartState(PartID,2)-ElemBaryNGeo(2,ElemID)) &
                             +(PartState(PartID,3)-ElemBaryNGeo(3,ElemID))*(PartState(PartID,3)-ElemBaryNGeo(3,ElemID)) )
      IF(Distance(iBGMElem).GT.ElemRadius2NGeo(ElemID))THEN
        Distance(iBGMElem)=-1.0
      END IF
    END IF
  END DO ! nBGMElems
  !CALL BubbleSortID(Distance,ListDistance,nBGMElems)
  CALL InsertionSort(Distance,ListDistance,nBGMElems)

  OldXi=PartPosRef(1:3,PartID)
  newXi=HUGE(1.0)
  newElemID=-1
  ! loop through sorted list and start by closest element  
  DO iBGMElem=1,nBGMElems
    IF(ALMOSTEQUAL(Distance(iBGMELem),-1.0)) CYCLE
    ElemID=ListDistance(iBGMElem)
    CALL Eval_xyz_ElemCheck(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.BezierClipHit) THEN ! particle inside
    IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).LE.1.0) THEN ! particle inside
      PEM%Element(PartID) = ElemID
      PDM%ParticleInside(PartID)=.TRUE.
      EXIT
    END IF
    IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).LT.MAXVAL(ABS(newXi))) THEN
      newXi=PartPosRef(1:3,PartID)
      newElemID=ElemID
    END IF
  END DO ! iBGMElem
  IF(PDM%ParticleInside(PartID))THEN
    RETURN
  ELSE
    ! use best xi
    !IPWRITE(UNIT_stdOut,*) ' recover particle', iPart
    IF(MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi)))THEN
      PartPosRef(1:3,PartID)=OldXi
      PEM%Element(PartID)=oldElemID
    ELSE
      PartPosRef(1:3,PartID)=NewXi
      PEM%Element(PartID)=NewElemID
    END IF
    PDM%ParticleInside(PartID)=.TRUE.
    IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).GT.epsOneCell) THEN
      IPWRITE(UNIT_stdOut,*) ' xi          ', PartPosRef(1:3,PartID)
      IPWRITE(UNIT_stdOut,*) ' newxi       ', newXi
      IPWRITE(UNIT_stdOut,*) ' oldxi       ', oldXi
      IPWRITE(UNIT_stdOut,*) ' ParticlePos ', PartState(PartID,1:3)
      IPWRITE(UNIT_stdOut,*) ' oldPartPosi ', LastPartPos(PartID,1:3)
      IPWRITE(UNIT_stdOut,*) ' initoutpos  ', oldestpos
      IPWRITE(UNIT_stdOut,*) ' correctedpos', blubb
      IPWRITE(UNIT_stdOut,*) ' Trajectory  ', PartTrajectory
      IPWRITE(UNIT_stdOut,*) ' lengthT     ', LengthPartTrajectory
#ifdef MPI
      InElem=PEM%Element(PartID)
      IF(InElem.LE.PP_nElems)THEN
      IPWRITE(UNIT_stdOut,*) ' ElemID       ', InElem+offSetElem
      ELSE
        IPWRITE(UNIT_stdOut,*) ' ElemID       ', offSetElemMPI(PartHaloElemToProc(NATIVE_PROC_ID,InElem)) &
                                               + PartHaloElemToProc(NATIVE_ELEM_ID,InElem)
      END IF
#else
      IPWRITE(UNIT_stdOut,*) ' ElemID       ', PEM%Element(PartID)+offSetElem
#endif
      CALL abort(&
      __STAMP__ &
      ,'Particle Not inSide of Element, PartID',PartID)

    END IF
  END IF
END IF

END SUBROUTINE ReComputeParticleBCInteraction


SUBROUTINE ComputeFaceIntersection(ElemID,firstSide,LastSide,nlocSides,PartID,PartIsDone)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:SideType
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
USE MOD_Particle_Mesh_Vars,          ONLY:ElemBaryNGeo
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Mesh_Vars,          ONLY:BCElem
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,ComputePlanarIntersectionBezierRobust2
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
USE MOD_Particle_Vars,               ONLY:PartPosRef
USE MOD_Eval_xyz,                    ONLY:Eval_XYZ_Poly
USE MOD_Mesh_Vars,                   ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
LOGICAL,INTENT(INOUT)         :: PartisDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID, locSideList(firstSide:lastSide), hitlocSide
LOGICAL                       :: dolocSide(firstSide:lastSide),ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip,BCSideID
REAL                          :: tmpPos(3), tmpLastPartPos(3),tmpVec(3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory,xNodes(1:3,1:4)
!===================================================================================================================================


tmpPos=PartState(PartID,1:3)
tmpLastPartPos(1:3)=LastPartPos(PartID,1:3)
PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
tmpVec=PartTrajectory

LastPartPos(PartID,1:3)=PartState(PartID,1:3)
!PartState(PartID,1:3)=ElemBaryNGeo(:,ElemID)
LastPartPos(PartID,1:3)=ElemBaryNGeo(:,ElemID)

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory


locAlpha=-1.0
nInter=0
dolocSide=.TRUE.
!nlocSides=lastSide-firstSide+1
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID=BCElem(ElemID)%BCSideID(ilocSide)
  BCSideID=PartBCSideList(SideID)
  locSideList(ilocSide)=ilocSide
  ! get correct flip
  flip  = 0 
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT)
    !CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
    CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)            &
                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)

!    CALL ComputePlanarIntersectionBezierRobust2(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                  ,xi (ilocSide)      &
!                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)
!
!                                                                            !,eta(ilocSide)   ,PartID,ilocSide,SideID,ElemID)
  CASE(BILINEAR,PLANAR_NONRECT)
    xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,BCSideID)
    xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,BCSideID)
    xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,BCSideID)
    xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,BCSideID)
    !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,xNodes &
    CALL ComputeBiLinearIntersectionRobust(isHit,xNodes &
                                                 ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                     ,xi (ilocSide)      &
                                                                                     ,eta(ilocSide)      &
                                                                                     ,PartID,flip,BCSideID)

    !CALL ComputeBiLinearIntersectionSuperSampled2(isHit,[BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
    !                                                    ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
    !                                                    ,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
    !                                                                                        ,xi (ilocSide)      &
    !                                                                                        ,eta(ilocSide)  ,PartID,flip,SideID)
!    CALL ComputeBezierIntersection(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                      ,xi (ilocSide)      &
!                                                                      ,eta(ilocSide)      ,PartID,SideID)

  CASE(CURVED)
    CALL ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                            ,xi (ilocSide)      &
                                                                            ,eta(ilocSide)      ,PartID,BCSideID)
  END SELECT
  IF(locAlpha(ilocSide).GT.-1.0)THEN
    nInter=nInter+1
  END IF
END DO ! ilocSide

IF(nInter.EQ.0) THEN
  !IPWRITE(*,*) 'not found',PartID
  !IPWRITE(*,*) 'ElemBary',LastPartPos(PartID,1:3)
  !IPWRITE(*,*) 'Part-Pos',tmpPos
  !IPWRITE(*,*) 'LastPart-Pos',tmpLastPartPos
  PartState(PartID,1:3)=tmpPos
  LastPartPos(PartID,1:3)=tmpLastPartPos(1:3)
  IF(PartPosRef(1,PartID).GT. 1.) PartPosRef(1,PartID)= 0.99
  IF(PartPosRef(1,PartID).LT.-1.) PartPosRef(1,PartID)=-0.99
  IF(PartPosRef(2,PartID).GT. 1.) PartPosRef(2,PartID)= 0.99
  IF(PartPosRef(2,PartID).LT.-1.) PartPosRef(2,PartID)=-0.99
  IF(PartPosRef(3,PartID).GT. 1.) PartPosRef(3,PartID)= 0.99
  IF(PartPosRef(3,PartID).LT.-1.) PartPosRef(3,PartID)=-0.99
  CALL Eval_xyz_Poly(PartPosRef(:,PartID),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),PartState(PartID,1:3))
  ! crash
  RETURN
ELSE
  ! take first possible intersection
  !CALL BubbleSortID(locAlpha,locSideList,6)
  CALL InsertionSort(locAlpha,locSideList,nlocSides)
  DO ilocSide=firstSide,LastSide
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      hitlocSide=locSideList(ilocSide)
      !SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
      SideID=BCElem(ElemID)%BCSideID(hitlocSide)
      BCSideID=PartBCSideList(SideID)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+0.97*locAlpha(ilocSide)*PartTrajectory
      PartState(PartID,1:3)  =LastPartPos(PartID,1:3)!+tmpVec
      !PartState(PartID,1:3)  =PartState(PartID,1:3)+locAlpha(ilocSide)*PartTrajectory
      !PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
      !lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
      !                         +PartTrajectory(2)*PartTrajectory(2) &
      !                         +PartTrajectory(3)*PartTrajectory(3) )
      !PartTrajectory=PartTrajectory/lengthPartTrajectory
      !locAlpha(ilocSide)=0.
      !CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
      !                                                                  ,xi(hitlocSide)     &
      !                                                                  ,eta(hitlocSide)    &
      !                                                                  ,PartId,SideID,ElemID)
      !IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
    END IF ! locAlpha>-1.0
  END DO ! ilocSide
END IF ! nInter>0

END SUBROUTINE ComputeFaceIntersection


END MODULE MOD_Particle_Tracking
