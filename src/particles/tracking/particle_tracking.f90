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
  MODULE PROCEDURE ParticleRefTracking
END INTERFACE

!PUBLIC::ParticleTracking,ParticleTrackingCurved
PUBLIC::ParticleTrackingCurved
PUBLIC::ParticleRefTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTrackingCurved()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:nBCSides,NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:epsilontol,SideType,epsilonOne,epsilonbilinear
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D,BoundingBoxIsEmpty
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide,PartSideToElem,nTotalSides,nTotalElems
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToElem,GEO
USE MOD_TimeDisc_Vars,               ONLY:iter
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Particle_Vars,               ONLY:time
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Tracking_vars,      ONLY:ntracks
USE MOD_Particle_Mesh,               ONLY:SingleParticleToExactElementNoMap
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,PartInElemCheck
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
#ifdef MPI
USE MOD_Mesh_Vars,                   ONLY:BC,nSides
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID,flip
INTEGER                       :: ilocSide,SideID, locSideList(1:6), hitlocSide,nInterSections,nLoc
LOGICAL                       :: PartisDone,dolocSide(1:6),isHit,markTol
REAL                          :: localpha(1:6),xi(1:6),eta(1:6)
!INTEGER                       :: lastlocSide
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory,xNodes(1:3,1:4)
!===================================================================================================================================

DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    nTracks=nTracks+1
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )
    PartTrajectory=PartTrajectory/lengthPartTrajectory!+epsilontol
    !PartTrajectory=PartTrajectory/(lengthPartTrajectory+epsilontol)
    lengthPartTrajectory=lengthPartTrajectory!+epsilontol
    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    DO WHILE (.NOT.PartisDone)
      !IF((iPart.EQ.238).AND.(iter.GE.182)) WRITE(*,*) 'ElemID',ElemID
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
        CASE(PLANAR)
          !IF(iPart.EQ.40) WRITE(*,*) 'planar'
          CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                        ,xi (ilocSide)      &
                                                                                        ,eta(ilocSide)   ,iPart,flip,SideID)

!          CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                        ,xi (ilocSide)      &
!                                                                                        ,eta(ilocSide)   ,iPart,flip,SideID)
!                                                                              !,doTest=.TRUE.)
!                                                                                  !,eta(ilocSide)   ,iPart,ilocSide,SideID,ElemID)
        CASE(BILINEAR)
          !IF(iPart.EQ.40) WRITE(*,*) 'bilinear'
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
        CALL SelectInterSectionType(PartIsDone,doLocSide,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                         ,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,ElemID)
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
        CALL SelectInterSectionType(PartIsDone,doLocSide,hitlocSide,ilocSide,PartTrajectory,lengthPartTrajectory &
                                         ,xi(hitlocSide),eta(hitlocSide),localpha(ilocSide),iPart,SideID,ElemID)
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


SUBROUTINE ParticleBCTracking(ElemID,firstSide,LastSide,nlocSides,PartId,PartisDone)
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY:nBCSides,NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:epsilontol,SideType,epsilonOne,epsilonbilinear
USE MOD_Particle_Surfaces_Vars,      ONLY:BezierControlPoints3D,BoundingBoxIsEmpty
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide,PartSideToElem,nTotalSides,nTotalElems,PartBCSideList,nTotalBCSides
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToElem
USE MOD_TimeDisc_Vars,               ONLY:iter
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Vars,               ONLY:time
USE MOD_Particle_Mesh_Vars,          ONLY:SidePeriodicType, SidePeriodicDisplacement,BCElem
USE MOD_Utils,                       ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Tracking_Vars,      ONLY:ntracks
USE MOD_Particle_Intersection,       ONLY:ComputeBezierIntersection,ComputeBiLinearIntersectionSuperSampled2 &
                                         ,ComputePlanarIntersectionBezier,ComputePlanarIntersectionBezierRobust2
USE MOD_Particle_Intersection,       ONLY:ComputePlanarIntersectionBezierRobust,ComputeBiLinearIntersectionRobust
#ifdef MPI
USE MOD_Mesh_Vars,                   ONLY:BC,nSides
#endif /*MPI*/
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
PartTrajectory=PartTrajectory/lengthPartTrajectory
lengthPartTrajectory=lengthPartTrajectory!+epsilontol

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
  CASE(PLANAR)
    !CALL ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
    CALL ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                                  ,xi (ilocSide)      &
                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)

!    CALL ComputePlanarIntersectionBezierRobust2(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
!                                                                                  ,xi (ilocSide)      &
!                                                                                  ,eta(ilocSide)   ,PartID,flip,BCSideID)
!
!                                                                            !,eta(ilocSide)   ,PartID,ilocSide,SideID,ElemID)
  CASE(BILINEAR)
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
  DO ilocSide=firstSide,LastSide
    IF(locAlpha(ilocSide).GT.-1.0)THEN
      hitlocSide=locSideList(ilocSide)
      !SideID=PartElemToSide(E2S_SIDE_ID,hitlocSide,ElemID)
      SideID=BCElem(ElemID)%BCSideID(hitlocSide)
      BCSideID=PartBCSideList(SideID)
      CALL GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                                                        ,xi(hitlocSide)     &
                                                                        ,eta(hitlocSide)    &
                                                                        ,PartId,SideID,ElemID)
      IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
    END IF ! locAlpha>-1.0
  END DO ! ilocSide
END IF ! nInter>0

END SUBROUTINE ParticleBCTracking


SUBROUTINE ParticleRefTracking()
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Particle_Vars,           ONLY:PDM,PEM,PartState,PartPosRef,lastpartpos
USE MOD_Mesh_Vars,               ONLY:OffSetElem
USE MOD_Eval_xyz,                ONLY:eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars,  ONLY:nTracks
USE MOD_Particle_Surfaces_Vars,  ONLY:ClipHit
USE MOD_Particle_Mesh_Vars,      ONLY:Geo,IsBCElem,BCElem,epsInCell,epsOneCell
USE MOD_Particle_Mesh_Vars,      ONLY:PartElemToSide,PartSideToElem,nTotalBCSides,PartBCSideLIst
USE MOD_Utils,                   ONLY:BubbleSortID,InsertionSort
USE MOD_Particle_Mesh_Vars,      ONLY:ElemBaryNGeo,ElemRadiusNGeo
USE MOD_TimeDisc_Vars,           ONLY:iter
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloToProc
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart, ElemID,oldElemID,iElem, newElemID
INTEGER                     :: CellX,CellY,CellZ,iBGMElem,nBGMElems,nLocSides,locSideID,tmplocSideID
REAL,ALLOCATABLE            :: Distance(:)
REAL                        :: oldXi(3),newXi(3), LastPos(3),epsLowOne,xi,eta
INTEGER,ALLOCATABLE         :: ListDistance(:)
!REAL                        :: epsOne
#ifdef MPI
INTEGER                     :: InElem
#endif
INTEGER                     :: SideID,BCSideID,ilocSide
LOGICAL                     :: ParticleFound(1:PDM%ParticleVecLength)
!===================================================================================================================================


! sanity check
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF(PartState(iPart,1).NE.PartState(iPart,1)) CALL abort(&
        __STAMP__,&
        ' x-pos is nan of ipart',iPart)
    IF(PartState(iPart,2).NE.PartState(iPart,2)) CALL abort(&
        __STAMP__,&
        ' y-pos is nan of ipart',iPart)
    IF(PartState(iPart,3).NE.PartState(iPart,3)) CALL abort(&
        __STAMP__,&
        ' z-pos is nan of ipart',iPart)
  END IF
END DO ! iPart=1,,PDM%ParticleVecLength

!epsOne=1.0+epsInCell
ParticleFound=.FALSE.
epsLowOne=1.0-2.0*epsInCell
! first step, reuse Elem cache, therefore, check if particle are still in element, if not, search later
DO iElem=1,PP_nElems ! loop only over internal elems, if particle is already in HALO, it shall not be found here
  DO iPart=1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(iPart))THEN
      ElemID = PEM%lastElement(iPart)
      !IF(ElemID.GT.PP_nElems) IPWRITE(*,*) 'too large',ElemID,PP_nElems
      !IF(ElemID.LT.1)         IPWRITE(*,*) 'too small',ElemID,1
      !IF(ElemID.EQ.-888)      IPWRITE(*,*) 'not set',ElemID
      IF(ElemID.NE.iElem) CYCLE
      nTracks=nTracks+1
      ! sanity check
      !IF(PartState(iPart,3).GE.0.089)THEN
      !  IPWRITE(UNIT_stdOut,*) ' Part out of area, z,ipart', PartState(iPart,3),iPart
      !END IF
      IF(GEO%nPeriodicVectors.GT.0)THEN
        ! call here function for mapping of partpos and lastpartpos
        CALL PeriodicMovement(iPart)
      END IF
      IF(IsBCElem(ElemID))THEN
        !  multiple BC intersections
        !LastPos=PartState(iPart,1:3)
        ! check inner sides
        CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%nInnerSides,BCElem(ElemID)%nInnerSides,iPart,ParticleFound(iPart))
        IF(ParticleFound(iPart)) CYCLE
        ! check if particle is still in element
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
        ! particle can hit outside bc due to tolerance
        LastPos=PartState(iPart,1:3)
        IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.epsLowOne) THEN ! particle can be outside
          CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%nInnerSides,BCElem(ElemID)%nInnerSides,iPart,ParticleFound(iPart))
          IF(ParticleFound(iPart)) CYCLE
          IF(GEO%nPeriodicVectors.GT.0)THEN
            CALL PeriodicMovement(iPart)
          END IF
        ELSE ! only external faces
          nlocSides=BCElem(ElemID)%lastSide-BCElem(ElemID)%nInnerSides
          CALL ParticleBCTracking(ElemID,BCElem(ElemID)%nInnerSides+1,BCElem(ElemID)%lastSide,nlocSides,iPart,ParticleFound(iPart))
          IF(ParticleFound(iPart)) CYCLE
          IF(GEO%nPeriodicVectors.GT.0)THEN
            CALL PeriodicMovement(iPart)
          END IF
        END IF
        DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(iPart,1)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(iPart,2)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(iPart,3)) )
          LastPos=PartState(iPart,1:3)
          ! unfortunately, here all sides
          CALL ParticleBCTracking(ElemID,1,BCElem(ElemID)%lastSide,BCElem(ElemID)%lastSide,iPart,ParticleFound(iPart))
          IF(ParticleFound(iPart)) EXIT
          IF(GEO%nPeriodicVectors.GT.0)THEN
            ! call here function for mapping of partpos and lastpartpos
            CALL PeriodicMovement(iPart)
          END IF
        END DO
        IF(ParticleFound(iPart)) CYCLE
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
      ELSE ! no bc elem, therefore, no bc ineraction possible
        CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
      END IF ! initial check
      !IF(iPart.EQ.1) print*,'pos,elem',PartPosRef(:,ipart),ElemID
      !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.ClipHit) THEN ! particle inside
      IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.epsOneCell) THEN ! particle inside
        PEM%Element(iPart)  = ElemID
        ParticleFound(iPart)=.TRUE.
      !ELSE IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.5) THEN
      !  IPWRITE(UNIT_stdOut,*) ' partposref to large!',iPart
      END IF
      !IF((iPart.EQ.85).AND.(iter.GE.30))THEN
      !  WRITE(*,*) ' PartPosRef', PartPosRef(1:3,iPart)
      !  WRITE(*,*) ' max val',MAXVAL(ABS(PartPosRef(1:3,iPart)))
      !  WRITE(*,*) ' ElemID', ElemID
      !  WRITE(*,*) ' Found', ParticleFound(iPart)
      !END IF
    ELSE
      ! caution: dummy, because particle is not inside and such a particle sould not be searched for in 
      ! the next loop
      ParticleFound(iPart)=.TRUE.
    END IF
  END DO ! iPart
END DO ! iElem

! now, locate not all found particle
DO iPart=1,PDM%ParticleVecLength
  IF(ParticleFound(iPart)) CYCLE
  ! relocate particle
  oldElemID = PEM%lastElement(iPart) ! this is not!  a possible elem
  !IF(oldElemID.GT.PP_nElems) IPWRITE(*,*) 'second too large'
  !IF(oldElemID.LT.1)         IPWRITE(*,*) 'second too small'
  !IF(oldElemID.EQ.-888)      IPWRITE(*,*) 'second not set'
  ! get background mesh cell of particle
  CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
  !CellX = MIN(GEO%FIBGMimax,CellX)
  !CellX = MAX(MIN(GEO%FIBGMimax,CellX),GEO%FIBGMimin)
  CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
  !CellY = MIN(GEO%FIBGMjmax,CellY)
  !CellY = MAX(MIN(GEO%FIBGMjmax,CellY),GEO%FIBGMjmin)
  CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
  !CellZ = MAX(MIN(GEO%FIBGMkmax,CellZ),GEO%FIBGMkmin)
  !CellZ = MIN(GEO%FIBGMkmax,CellZ)

        
  ! check all cells associated with this beckground mesh cell
  nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem
  ALLOCATE( Distance(1:nBGMElems) &
          , ListDistance(1:nBGMElems) )

  ! get closest element barycenter
  Distance=0.
  ListDistance=0
  DO iBGMElem = 1, nBGMElems
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    ListDistance(iBGMElem)=ElemID
    IF(ElemID.EQ.-1)CYCLE
    IF(ElemID.EQ.OldElemID)THEN
      Distance(iBGMElem)=-1.0
    ELSE
      Distance(iBGMElem)=SQRT((PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))  &
                             +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
                             +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID)) )
      IF(Distance(iBGMElem).GT.ElemRadiusNGeo(ElemID))THEN
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
    CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),ElemID)
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.ClipHit) THEN ! particle inside
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.epsOneCell) THEN ! particle inside
      PEM%Element(iPart) = ElemID
      ParticleFound(iPart)=.TRUE.
      EXIT
    END IF
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
      newXi=PartPosRef(1:3,iPart)
      newElemID=ElemID
    END IF
  END DO ! iBGMElem
  !IF(iPart.EQ.1) print*,'pos,elem',PartPosRef(:,ipart),ElemID
  !IF((iPart.EQ.85).AND.(iter.GE.30))THEN
  !  WRITE(*,*) 'Re-Found?', ParticleFound(iPart)
  !  WRITE(*,*) 'OldXi',PartPosRef(1:3,iPart)
  !  WRITE(*,*) 'best-xi',NewXi
  !END IF

  IF(.NOT.ParticleFound(iPart))THEN
    ! use best xi
    !IPWRITE(UNIT_stdOut,*) ' recover particle', iPart
    IF(MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi)))THEN
      PartPosRef(1:3,iPart)=OldXi
      PEM%Element(iPart)=oldElemID
    ELSE
      PartPosRef(1:3,iPart)=NewXi
      PEM%Element(iPart)=NewElemID
    END IF

    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.1.1) THEN
      ! check in which dimension the particle left the cell
      ! backtraze
      oldXi=ABS(PartPosRef(1:3,iPart))
      locSideID=-1
      DO ilocSide=1,6
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,PEM%Element(iPart))
        BCSideID=PartBCSideList(SideID)
        IF((BCSideID.GE.1).AND.(BCSideID.LE.nTotalBCSides))THEN
          SELECT CASE(ilocSide)
          CASE(XI_MINUS)
            IF(PartPosRef(1,iPart).LT.-1.0) THEN
              locSideID=ilocSide
              !BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
              Xi =PartPosRef(3,iPart)
              Eta=PartPosRef(2,iPart)
              PartPosRef(1,iPart)=-0.9999
              EXIT
            END IF
          CASE(XI_PLUS)
            IF(PartPosRef(1,iPart).GT.1.0) THEN
              locSideID=ilocSide
              Xi =PartPosRef(2,iPart)
              Eta=PartPosRef(3,iPart)
              PartPosRef(1,iPart)= 0.9999
              EXIT
            END IF
          CASE(ETA_MINUS)
            IF(PartPosRef(2,iPart).LT.-1.0) THEN
              locSideID=ilocSide
               Xi =PartPosRef(1,iPart)
               Eta=PartPosRef(3,iPart)
               PartPosRef(2,iPart)=-0.9999
              EXIT
            END IF
          CASE(ETA_PLUS)
            IF(PartPosRef(2,iPart).GT.1.0) THEN
              locSideID=ilocSide
              Xi =-PartPosRef(1,iPart)
              Eta=PartPosRef(3 ,iPart)
              PartPosRef(2,iPart)= 0.9999
              EXIT
            END IF
          CASE(ZETA_MINUS)
            IF(PartPosRef(3,iPart).LT.-1.0) THEN
              locSideID=ilocSide
              Xi =PartPosRef(2,iPart)
              Eta=PartPosRef(1,iPart)
              PartPosRef(3,iPart)=-0.9999
              EXIT
            END IF
          CASE(ZETA_PLUS)
            IF(PartPosRef(3,iPart).GT.1.0) THEN
              locSideID=ilocSide
              Xi =PartPosRef(1,iPart)
              Eta=PartPosRef(2,iPart)
              PartPosRef(3,iPart)= 0.9999
              EXIT
            END IF
          END SELECT
        END IF ! BCSideID
      END DO ! ilocSide

      !IF( oldXI(1).GT.oldXI(2) .AND. OldXi(1) .GT. oldXI(3))THEN
      !  tmplocSideID=1
      !ELSEIF( oldXI(2).GT.oldXI(1) .AND. OldXi(2) .GT. oldXI(3))THEN
      !  tmplocSideID=2
      !ELSE
      !  tmplocSideID=3
      !END IF
      !tmpLocSideID=MAXLOC(oldXi)
      !SELECT CASE(tmpLocSideID)
      !CASE(1)
      !  IF(PartPosRef(tmpLocSideID,iPart).GT.0)THEN
      !    LocSideID=XI_PLUS
      !  ELSE
      !    LocSideID=XI_MINUS
      !  END IF
      !CASE(2)
      !  IF(PartPosRef(tmpLocSideID,iPart).GT.0)THEN
      !    LocSideID=ETA_PLUS
      !  ELSE
      !    LocSideID=ETA_MINUS
      !  END IF
      !CASE(3)
      !  IF(PartPosRef(tmpLocSideID,iPart).GT.0)THEN
      !    LocSideID=ZETA_PLUS
      !  ELSE
      !    LocSideID=ZETA_MINUS
      !  END IF
      !END SELECT

      !SideID=PartElemToSide(E2S_SIDE_ID,locSideID,PEM%Element(iPart))
      !BCSideID=PartBCSideList(SideID)

      !IF((BCSideID.GE.1).AND.(BCSideID.LE.nTotalBCSides))THEN
      IF(locSideID.GT.0)THEN
        CALL ReComputeParticleBCInteraction(xi,eta,LocSideID,SideID,BCSideID,iPart) 
      ELSE ! If no BC SideID
        IPWRITE(UNIT_stdOut,*) ' xi          ', PartPosRef(1:3,iPart)
        IPWRITE(UNIT_stdOut,*) ' ParticlePos ', PartState(iPart,1:3)
#ifdef MPI
        InElem=PEM%Element(iPart)
        IF(InElem.LE.PP_nElems)THEN
        IPWRITE(UNIT_stdOut,*) ' ElemID       ', InElem+offSetElem
        ELSE
          IPWRITE(UNIT_stdOut,*) ' ElemID       ', offSetElemMPI(PartHaloToProc(NATIVE_PROC_ID,InElem)) &
                                                 + PartHaloToProc(NATIVE_ELEM_ID,InElem)
        END IF
#else
        IPWRITE(UNIT_stdOut,*) ' ElemID       ', PEM%Element(iPart)+offSetElem
#endif
        CALL abort(&
            __STAMP__, &
            'Particle Not inSide of Element, iPart',iPart)
      END IF ! BCSideID
    END IF

    ParticleFound(iPart)=.TRUE.
    !WRITE(*,*) ' iPart  :   ', iPart
    !WRITE(*,*) ' PartPos:   ', PartState(iPart,1:3)
    !WRITE(*,*) ' oldxi:     ', oldXi
    !WRITE(*,*) ' bestxi:    ', newXi
    !WRITE(*,*) ' oldelemid: ', oldElemID
    !WRITE(*,*) ' bestelemid ', NewElemID
!    CALL abort(__STAMP__,&
!      ' particle lost !! ')
  END IF
  DEALLOCATE( Distance)
  DEALLOCATE( ListDistance)

END DO ! iPart


END SUBROUTINE ParticleRefTracking


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
USE MOD_Particle_Mesh_Vars,          ONLY:PartElemToSide,PartSideToElem,nTotalSides,nTotalElems
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
        IF(SideType(SideID).EQ.PLANAR) THEN
          dolocSide(hitlocSide)=.FALSE.
        END IF
      ELSE
        ! inner side
        dolocSide=.TRUE.
        dolocSide(PartElemToElem(E2E_NB_LOC_SIDE_ID,hitlocSide,ElemID))=.FALSE.
        ElemID=PartElemToElem(E2E_NB_ELEM_ID,hitlocSide,ElemID)
        IF(ElemID.LE.0) CALL abort(&
            __STAMP__,&
           ' HaloRegion too small or critical error during halo region reconstruction!')
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


SUBROUTINE PeriodicMovement(PartID) 
!----------------------------------------------------------------------------------------------------------------------------------!
! move particle in the periodic direction, if particle is outside of the box
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Mesh_Vars,          ONLY:GEO
#ifdef MPI
USE MOD_Particle_MPI_Vars,           ONLY:PartShiftVector
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)              :: PartID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPV
!===================================================================================================================================

#ifdef MPI
PartShiftVector(1:3,PartID)=PartState(PartID,1:3)
#endif /*MPI*/
! x direction
IF(GEO%directions(1)) THEN
  IF(PartState(PartID,1).GT.GEO%xmaxglob) THEN
    DO iPV=1,GEO%nPeriodicVectors
      IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
    END DO
    IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
      PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    END IF
  END IF
  IF(PartState(PartID,1).LT.GEO%xminglob) THEN
    DO iPV=1,GEO%nPeriodicVectors
      IF(GEO%DirPeriodicVectors(iPV).EQ.1) EXIT
    END DO
    IF(GEO%PeriodicVectors(1,iPV).GT.0)THEN
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
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
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    END IF
  END IF
  IF(PartState(PartID,2).LT.GEO%yminglob) THEN
    DO iPV=1,GEO%nPeriodicVectors
      IF(GEO%DirPeriodicVectors(iPV).EQ.2) EXIT
    END DO
    IF(GEO%PeriodicVectors(2,iPV).GT.0)THEN
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
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
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    END IF
  END IF
  IF(PartState(PartID,3).LT.GEO%zminglob) THEN
    DO iPV=1,GEO%nPeriodicVectors
      IF(GEO%DirPeriodicVectors(iPV).EQ.3) EXIT
    END DO
    IF(GEO%PeriodicVectors(3,iPV).GT.0)THEN
      PartState(PartID,1:3)  =PartState(PartID,1:3)  +GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)+GEO%PeriodicVectors(1:3,iPV)
    ELSE
      PartState(PartID,1:3)  =PartState(PartID,1:3)  -GEO%PeriodicVectors(1:3,iPV)
      LastPartPos(PartID,1:3)=LastPartPos(PartID,1:3)-GEO%PeriodicVectors(1:3,iPV)
    END IF
  END IF
END IF

#ifdef MPI
PartShiftVector(1:3,PartID)=PartState(PartID,1:3)-PartShiftvector(1:3,PartID)
#endif /*MPI*/

END SUBROUTINE PeriodicMovement


SUBROUTINE ReComputeParticleBCInteraction(xi,eta,locSideID,SideID,BCSideID,PartID) 
!----------------------------------------------------------------------------------------------------------------------------------!
! The particle BC intersection is ignored. therefore, the particle is mapped onto the BC at the lost position and a wrong
! particle BC interaction is performed. FALLBACK!!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionRef
USE MOD_Particle_Vars,           ONLY:PDM,PEM,PartState,PartPosRef,lastpartpos
USE MOD_Eval_xyz,                ONLY:Eval_XYZ_Poly
USE MOD_Mesh_Vars,               ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Particle_Mesh_Vars,      ONLY:Geo,IsBCElem,epsInCell,epsOneCell,ElemRadiusNGeo,ElemBaryNGeo,BCElem
USE MOD_Particle_Mesh_Vars,      ONLY:GEO
USE MOD_Utils,                   ONLY:BubbleSortID,InsertionSort
USE MOD_Eval_xyz,                ONLY:eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars,  ONLY:nTracks
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloToProc
#endif
USE MOD_Mesh_Vars,               ONLY:OffSetElem
#ifdef MPI
USE MOD_MPI_Vars,                ONLY:offsetElemMPI
USE MOD_Particle_MPI_Vars,       ONLY:PartHaloToProc
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
REAL                        :: oldXi(3),newXi(3), LastPos(3),epsLowOne
INTEGER,ALLOCATABLE         :: ListDistance(:)
#ifdef MPI
INTEGER                     :: inElem
#endif /*MPI*/
!===================================================================================================================================

! SELECT CASE(locSideID)
! CASE(XI_MINUS)
!   !BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
!   Xi =PartPosRef(3,PartID)
!   Eta=PartPosRef(2,PartID)
!   PartPosRef(1,PartID)=-0.9999
! CASE(XI_PLUS)
!   Xi =PartPosRef(2,PartID)
!   Eta=PartPosRef(3,PartID)
!   PartPosRef(1,PartID)= 0.9999
! CASE(ETA_MINUS)
!   Xi =PartPosRef(1,PartID)
!   Eta=PartPosRef(3,PartID)
!   PartPosRef(2,PartID)=-0.9999
! CASE(ETA_PLUS)
!   !BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,NGeo-p,q)
!   ! hopefully correct
!   Xi =-PartPosRef(1,PartID)
!   Eta=PartPosRef(3,PartID)
!   PartPosRef(2,PartID)= 0.9999
! CASE(ZETA_MINUS)
!   Xi =PartPosRef(2,PartID)
!   Eta=PartPosRef(1,PartID)
!   PartPosRef(3,PartID)=-0.9999
! CASE(ZETA_PLUS)
!   Xi =PartPosRef(1,PartID)
!   Eta=PartPosRef(2,PartID)
!   PartPosRef(3,PartID)= 0.9999
! END SELECT

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory

! compute temporary last particle position
ElemID=PEM%Element(PartID)
CALL Eval_xyz_Poly(PartPosRef(:,PartID),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),LastPartPos(PartID,1:3))
PartState(PartID,1:3)=LastPartPos(PartID,1:3)+PartTrajectory*lengthPartTrajectory

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
IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).LE.epsOneCell) THEN ! particle inside
  PEM%Element(PartID)  = ElemID
  RETURN
ELSE
  ! particle has to to located in its final cell
  oldElemID=ElemID
  PDM%ParticleInside(PartID)=.FALSE.
  ! BGM cell
  CellX = CEILING((PartState(PartID,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
  !CellX = MIN(GEO%FIBGMimax,CellX)
  !CellX = MAX(MIN(GEO%FIBGMimax,CellX),GEO%FIBGMimin)
  CellY = CEILING((PartState(PartID,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
  !CellY = MIN(GEO%FIBGMjmax,CellY)
  !CellY = MAX(MIN(GEO%FIBGMjmax,CellY),GEO%FIBGMjmin)
  CellZ = CEILING((PartState(PartID,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
  !CellZ = MAX(MIN(GEO%FIBGMkmax,CellZ),GEO%FIBGMkmin)
  !CellZ = MIN(GEO%FIBGMkmax,CellZ)

  ! check all cells associated with this beckground mesh cell
  nBGMElems=GEO%FIBGM(CellX,CellY,CellZ)%nElem
  ALLOCATE( Distance(1:nBGMElems) &
          , ListDistance(1:nBGMElems) )

  ! get closest element barycenter
  Distance=0.
  ListDistance=0
  DO iBGMElem = 1, nBGMElems
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    ListDistance(iBGMElem)=ElemID
    IF(ElemID.EQ.-1)CYCLE
    IF(ElemID.EQ.OldElemID)THEN
      Distance(iBGMElem)=-1.0
    ELSE
      Distance(iBGMElem)=SQRT((PartState(PartID,1)-ElemBaryNGeo(1,ElemID))*(PartState(PartID,1)-ElemBaryNGeo(1,ElemID))  &
                             +(PartState(PartID,2)-ElemBaryNGeo(2,ElemID))*(PartState(PartID,2)-ElemBaryNGeo(2,ElemID)) &
                             +(PartState(PartID,3)-ElemBaryNGeo(3,ElemID))*(PartState(PartID,3)-ElemBaryNGeo(3,ElemID)) )
      IF(Distance(iBGMElem).GT.ElemRadiusNGeo(ElemID))THEN
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
    !IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LE.ClipHit) THEN ! particle inside
    IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).LE.epsOneCell) THEN ! particle inside
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
    IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).GT.1.1) THEN
      IPWRITE(UNIT_stdOut,*) ' xi          ', PartPosRef(1:3,PartID)
      IPWRITE(UNIT_stdOut,*) ' ParticlePos ', PartState(PartID,1:3)
#ifdef MPI
      InElem=PEM%Element(PartID)
      IF(InElem.LE.PP_nElems)THEN
      IPWRITE(UNIT_stdOut,*) ' ElemID       ', InElem+offSetElem
      ELSE
        IPWRITE(UNIT_stdOut,*) ' ElemID       ', offSetElemMPI(PartHaloToProc(NATIVE_PROC_ID,InElem)) &
                                               + PartHaloToProc(NATIVE_ELEM_ID,InElem)
      END IF
#else
      IPWRITE(UNIT_stdOut,*) ' ElemID       ', PEM%Element(PartID)+offSetElem
#endif
      CALL abort(&
          __STAMP__, &
          'Particle Not inSide of Element, PartID',PartID)

    END IF
  END IF
END IF

END SUBROUTINE ReComputeParticleBCInteraction

END MODULE MOD_Particle_Tracking
