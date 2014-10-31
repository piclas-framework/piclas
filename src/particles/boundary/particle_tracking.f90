#include "boltzplatz.h"

MODULE MOD_Particle_Tracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ParticleTracking
  MODULE PROCEDURE ParticleTracking
END INTERFACE

INTERFACE ParticleTrackingCurved
  MODULE PROCEDURE ParticleTrackingCurved
END INTERFACE

PUBLIC::ParticleTracking,ParticleTrackingCurved
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTrackingCurved()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
!USE MOD_Preproc
USE MOD_Globals,                     ONLY:abort
USE MOD_Mesh_Vars,                   ONLY:ElemToSide,nBCSides,NGeo!,NormVec
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:epsilontol,SideType,epsilonOne,neighborElemID,neighborlocSideID,epsilonbilinear
USE MOD_Particle_Surfaces_Vars,      ONLY:nPartCurved,BezierControlPoints3D,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,      ONLY:SuperSampledNodes,nQuads
USE MOD_TimeDisc_Vars,               ONLY:iter
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Particle_Vars,           ONLY:time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID!,LastElemID
INTEGER                       :: ilocSide,SideID,iIntersect
LOGICAL                       :: PartisDone,dolocSide(1:6)
REAL                          :: alpha,xi,eta
INTEGER                       :: lastlocSide
REAL                          :: oldXIntersection(1:3),distance,helpVec(1:3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
!===================================================================================================================================
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
    !LastElemID=ElemID
    ! DEBUUGGG
    !!!PartState(iPart,3)=LastPartPos(iPart,3)+2.0
    IF(DOT_PRODUCT(LastPartPos(iPart,1:2),LastPartPos(iPart,1:2)).GT.0.52)THEN
      print*,'wuuumm'
      STOP
    END IF
    IF(ABS(LastPartPos(ipart,3)).GT.0.2501)THEN
      print*,'out of space'
      STOP
    END IF

!    IF(time.GT.104.0)THEN
!      print*,'ElemID',ElemID
!    END IF
    !!! DEBUGG ENDE
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )
    PartTrajectory=PartTrajectory/lengthPartTrajectory
    lengthPartTrajectory=lengthPartTrajectory+epsilontol
    lastlocSide=-1
    oldXInterSection=HUGE(1.0)
    !print*,'partTrajectory',PartTrajectory


!    IF(iter.EQ.326)THEN
!      print*,'ElemID',ElemID
!      print*,'pos',LastPartPos(iPart,1:3)
!      print*,'Trajectory',PartTrajectory
!      read*
!    END IF
!    IF(iter.GE.420)THEN
!      print*,'ElemID',ElemID
!      print*,'pos',LastPartPos(iPart,1:3)
!      print*,'Trajectory',PartTrajectory
!      read*
!    END IF


    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    DO WHILE (.NOT.PartisDone)
      DO ilocSide=1,6
        alpha=-1.0
        IF(.NOT.dolocSide(ilocSide)) CYCLE
        SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
!        print*,'SideID',SideID
!        print*,'normvec',NormVec(:,0,0,SideID)
!        print*,'normvec',NormVec(:,PP_N,PP_N,SideID)
        ! find intersection point with surface: use NGeo for BezierControlPoints
        SELECT CASE(SideType(SideID))
        CASE(PLANAR)
          !xNodes(:,1)=BezierControlPoints3D(1:3,0   ,0   ,SideID)
          !xNodes(:,2)=BezierControlPoints3D(1:3,NGeo,0   ,SideID)
          !xNodes(:,3)=BezierControlPoints3D(1:3,NGeo,NGeo,SideID)
          !xNodes(:,4)=BezierControlPoints3D(1:3,0   ,NGeo,SideID)
!          CALL ComputePlanarIntersectionSuperSampled2([BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
!                                                      ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
!                                                      ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
!                                                      ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
!                                                      ,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart)
          CALL ComputePlanarIntersectionBezier(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
        CASE(BILINEAR)
          CALL ComputeBiLinearIntersectionSuperSampled2([BezierControlPoints3D(1:3,0   ,0   ,SideID)  &
                                                        ,BezierControlPoints3D(1:3,NGeo,0   ,SideID)  &
                                                        ,BezierControlPoints3D(1:3,NGeo,NGeo,SideID)  &
                                                        ,BezierControlPoints3D(1:3,0   ,NGeo,SideID)] &
                                                        ,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
        CASE(CURVED)
          CALL ComputeBezierIntersection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
        END SELECT
    !    print*,'alpha',alpha
    !    print*,'xi',xi
    !    print*,'eta',eta
        if(alpha.gt.0)THEN
    !    print*,'interpoint-x',LastPartPos(iPart,1)+alpha*PartTrajectory(1)
    !    print*,'interpoint-y',LastPartPos(iPart,2)+alpha*PartTrajectory(2)
    !    print*,'interpoint-z',LastPartPos(iPart,3)+alpha*PartTrajectory(3)
       endif
        !read*

        ! get correct intersection
        IF(alpha.EQ.-1.0) CYCLE
        IF(SideID.LE.nBCSides)THEN
          IF(SideType(SideID).NE.PLANAR)THEN
            IF(lastlocSide.EQ.ilocSide)THEN
              helpVec=oldXInterSection-LastPartPos(iPart,1:3)-alpha*PartTrajectory
              distance=SQRT(DOT_PRODUCT(helpVec,helpVec))
              IF(distance.LT.1e-5)THEN
                alpha=-1.0
                CYCLE
              END IF
            END IF
          END IF
          ! check if interesction is possible and take first intersection
          CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID)
          dolocSide=.TRUE.
          IF(SideType(SideID).NE.PLANAR) THEN
            lastlocSide=ilocSide
            oldXInterSection=LastPartPos(iPart,1:3)
          END IF
          dolocSide(ilocSide)=.FALSE.
          !IF(SideType(SideID).EQ.PLANAR) dolocSide(ilocSide)=.FALSE.
          EXIT
        ELSE ! no BC Side
          ! check if the found alpha statisfy the selection condition
          iInterSect=INT((ABS(xi)-2*epsilontol)/1.0)+INT((ABS(eta)-2*epsilontol)/1.0)
          IF(iInterSect.GT.0)THEN
            CALL abort(__STAMP__,&
                ' Particle went through edge or node. Not implemented yet.',999,999.)
          ELSE
            dolocSide=.TRUE.
            dolocSide(neighborlocSideID(ilocSide,ElemID))=.FALSE.
            ElemID=neighborElemID(ilocSide,ElemID)
            lastlocSide=-1
            EXIT
          END IF ! possible intersect
        END IF ! SideID.LT.nBCSides
      END DO ! ilocSide
      ! no intersection found
      IF(alpha.EQ.-1.0)THEN
        PEM%Element(iPart) = ElemID
        PartisDone=.TRUE.
      END IF
    END DO ! PartisDone=.FALSE.
  END IF ! Part inside
END DO ! iPart

END SUBROUTINE ParticleTrackingCurved

SUBROUTINE ComputeBezierIntersection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals,                 ONLY:Cross,abort
USE MOD_Mesh_Vars,               ONLY:NGeo,nBCSides
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne,SideDistance
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,ClipTolerance,ClipMaxInter,ClipMaxIter
USE MOD_Particle_Surfaces_Vars,  ONLY:locXi,locEta,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:arrayNchooseK,BoundingBoxIsEmpty
USE MOD_Utils,                   ONLY:BubbleSortID
USE MOD_Particle_Vars,           ONLY:time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)           :: PartTrajectory
REAL,INTENT(IN)                          :: lengthPartTrajectory
INTEGER,INTENT(IN)                       :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                         :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,iInter,p,q
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
REAL                                     :: coeffA,locSideDistance
INTEGER,ALLOCATABLE,DIMENSION(:)         :: locID
LOGICAL                                  :: foundInter,firstClip
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
alpha=-1.0
Xi   = 2.0
Eta  = 2.0

!!!! DEBUGGG fix side ID
!print*,'sideid',sideid
! If side is flat, than check if particle vector is perpenticular to side. if true, then particle moves parallel to or in side
IF(BoundingBoxIsEmpty(SideID))THEN
  IF(ABS(DOT_PRODUCT(PartTrajectory,SideNormVec(1:3,SideID))).LT.epsilontol) RETURN
!  coeffA=DOT_PRODUCT(SideNormVec(:,SideID),PartTrajectory)
!  IF(ABS(coeffA).LT.epsilontol) RETURN
!  !! interaction should be computed in last step
!  locSideDistance=SideDistance(SideID)-DOT_PRODUCT(LastPartPos(iPart,1:3),SideNormVec(:,SideID))
!  alpha=locSideDistance/coeffA
!  IF((alpha.GT.lengthPartTrajectory).OR.(alpha.LT.-epsilontol))THEN
!    alpha=-1.0
!    RETURN
!  END IF
!
!ELSE 
END IF ! BoundingBoxIsEmpty
  
! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
IF(.NOT.InsideBoundingBox(LastPartPos(iPart,1:3),SideID))THEN ! the old particle position is not inside the bounding box
  IF(.NOT.InsideBoundingBox(PartState(iPart,1:3),SideID))THEN ! the new particle position is not inside the bounding box
    IF(.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)) RETURN ! the particle does not intersect the 
                                                                                              ! bounding box
  END IF
END IF

!IF(time.GT.28.07)THEN
!  print*,'boundingbox hit- do clipping'
!  print*,'SideID',SideID
!  read*
!END IF
!
! 2.) Bezier intersection: transformation of bezier patch 3D->2D
IF(ABS(PartTrajectory(3)).LT.epsilontol)THEN
  n1=(/ -PartTrajectory(2)-PartTrajectory(3)  , PartTrajectory(1) ,PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3) , PartTrajectory(3) , -PartTrajectory(1)-PartTrajectory(2) /)
END IF
n1=n1/SQRT(DOT_PRODUCT(n1,n1))
n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
       , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
       , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
n2=n2/SQRT(DOT_PRODUCT(n2,n2))

!WRITE(*,'(A,2x,F12.8,2x,F12.8,2x,F12.8)') 'plane-n1',n1
!WRITE(*,'(A,2x,F12.8,2x,F12.8,2x,F12.8)') 'plane-n2',n2
!print*,'test n1*PartTrajectory',DOT_PRODUCT(PartTrajectory,n1)
!print*,'test n2*PartTrajectory',DOT_PRODUCT(PartTrajectory,n2)
!print*,'test n1*n2',DOT_PRODUCT(n1,n2)
!read*!CHANGETAG
DO q=0,NGeo
  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8,2x,F12.8)') &
!    'Bezier3D',p,q,BezierControlPoints3D(1,p,q,SideID),BezierControlPoints3d(2,p,q,SideID),BezierControlPoints3d(3,p,q,SideID)
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n1)
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n2)
!j    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D(1,p,q),BezierControlPoints2d(2,p,q)
!j    read*
  END DO
END DO

!print*,'init 3d patch'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8,2x,F12.8)') &
!    'Bezier3D',p,q,BezierControlPoints3D(1,p,q,SideID),BezierControlPoints3d(2,p,q,SideID),BezierControlPoints3d(3,p,q,SideID)
!  END DO
!END DO

!print*,'init 2d patch'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D(1,p,q),BezierControlPoints2d(2,p,q)
!  END DO
!END DO





!  this part in a new function or subroutine
locAlpha=-1.0
iClipIter=1
nXiClip=0
nEtaClip=0
nInterSections=0
firstClip=.TRUE.
CALL BezierClip(firstClip,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory&
               ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)

!IF(SideID.EQ.2)THEN
!  print*,'nInterSections',nInterSections
!  print*,'locAlpha',locAlpha
!  print*,'lengthPartTrajectory',lengthPartTrajectory
!  read*
!END IF
!
!
!IF(time.GT.103.0)THEN
!  print*, 'SideID',SideID
!  print*,'nInterSections',nInterSections
!  IF(SideID.EQ.2)THEN
!   print*,'LastPartPos',LastPartPos(iPart,1:3)
!   print*,'PartTrajectory',PartTrajectory
!   stop
!  END IF
!  IF(nInterSections.GE.1) read*
!  print*,'locAlpha',locAlpha
!  IF(SideID.EQ.8)THEN
!    print*,'nInterSections',nInterSections
!    print*,'locAlpha',locAlpha
!    read*
!  END IF
!END IF

! old and oobsolet
!IF(nInterSections.GT.1)THEN
!  ALLOCATE(locID(nInterSections))
!  DO iInter=1,nInterSections
!    locID(iInter)=iInter
!  END DO ! iInter
!END IF ! nInterSections

SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha=locAlpha(nInterSections)
  xi =locXi (1)
  eta=loceta(1)
CASE DEFAULT
  ! more than one intersection
  !ALLOCATE(locID(nInterSections))
  ALLOCATE(locID(nInterSections))
  DO iInter=1,nInterSections
    locID(iInter)=iInter
  END DO ! iInter
  ! sort intersection distance
  CALL BubbleSortID(locAlpha,locID,nIntersections)
  
  IF(SideID.LE.nBCSides)THEN
    ! requires first hit with BC
    DO iInter=1,nInterSections 
      IF(locAlpha(iInter).GT.-1.0)THEN
        alpha=locAlpha(iInter)
        xi =locXi (locID(iInter))
        eta=loceta(locID(iInter))
        SDEALLOCATE(locID)
        RETURN 
      END IF
    END DO ! iInter
  ELSE
    IF(MOD(nInterSections,2).EQ.0) THEN
      SDEALLOCATE(locID)
      RETURN ! leave and enter a cell multiple times
    END IF
    alpha=locAlpha(nInterSections)
    xi =locXi (locID(nInterSections))
    eta=loceta(locID(nInterSections))
  END IF
  DEALLOCATE(locID)
END SELECT


END SUBROUTINE ComputeBezierIntersection

RECURSIVE SUBROUTINE BezierClip(firstClip,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,iClipIter,nXiClip,nEtaClip&
                               ,nInterSections,iPart,SideID)
!================================================================================================================================
! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:XiArray,EtaArray,locAlpha,locXi,locEta
USE MOD_Particle_Surfaces_Vars,  ONLY:ClipTolerance,ClipMaxIter,ArrayNchooseK,FacNchooseK,ClipMaxInter
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,mEpsilontol
USE MOD_Particle_Vars,           ONLY:LastPartPos,Time
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
! INTEGER,INTENT(INOUT),DIMENSION(:)   :: locXi,locEta,locID
INTEGER,INTENT(INOUT)                :: iClipIter,nXiClip,nEtaClip,nInterSections
LOGICAL,INTENT(INOUT)                :: firstClip
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3,0:NGeo,0:NGeo)      :: ReducedBezierControlPoints
REAL,DIMENSION(0:NGeo,0:NGeo)        :: BezierControlPoints1D
REAL,DIMENSION(3)                    :: n1,n2,IntersectionVector
REAL,DIMENSION(2)                    :: LineNormVec
REAL                                 :: PatchDOF2D
REAL                                 :: minmax(1:2,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp(2,0:NGeo,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp2(2,0:NGeo,0:NGeo)
INTEGER                              :: p,q,idir,l,iDeCasteljau
REAL                                 :: Xi,Eta,XiMin,EtaMin,XiMax,EtaMax,XiSplit,EtaSplit,x,y,alpha
REAL                                 :: ZeroDistance,ClipTolerance2
LOGICAL                              :: DoXiClip,DoEtaClip,DoCheck
INTEGER                              :: iClip
REAL                                 :: PlusXi,MinusXi,PlusEta,MinusEta,tmpXi,tmpEta
INTEGER                              :: tmpnClip,tmpnXi,tmpnEta
!================================================================================================================================

PatchDOF2D=1.0/REAL((NGeo+1)*(NGeo+1))
ClipTolerance2=ClipTolerance*ClipTolerance

! 3.) Bezier intersection: solution Newton's method or Bezier clipping
! outcome: no intersection, single intersection, multiple intersection with patch
! BEZIER CLIPPING: xi- and eta-direction
IF(.NOT.FirstClip)THEN
  DoXiClip=.FALSE.
ELSE
  DoXiClip=.TRUE.
END IF
DoEtaClip=.TRUE.
DoCheck=.TRUE.
!nXiClip=0
!nEtaClip=0
!print*,'initer',iClipIter
DO iClipIter=iClipIter,ClipMaxIter
  !print*,'iClip',iClipIter
  !read*
  ! a) xi-direction
  IF(DoXiClip)THEN
    CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoCheck)
    IF(.NOT.DoCheck) EXIT
 !   print*,'LineVec',LineNormVec
    DO q=0,NGeo 
      DO p=0,NGeo
        BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec)
!        WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8)') 'Bezier1d',p,q,BezierControlPoints1D(p,q)
      END DO
    END DO
    DO l=0,NGeo
      minmax(2,l)=MAXVAL(BezierControlPoints1D(l,:)) 
      minmax(1,l)=MINVAL(BezierControlPoints1D(l,:)) 
    END DO ! l
    ! calc Smin and Smax and check boundaries
    CALL CalcSminSmax(minmax,XiMin,XiMax)
!IF(time.GT.28.07)THEN
!  IF(SideID.EQ.8)THEN
!    print*,'XiMinm,XiMax',XiMin,XiMax
!    read*
!  END IF
!END IF


    IF((XiMin.EQ.1.5).OR.(XiMax.EQ.-1.5))RETURN
    nXiClip=nXiClip+1
    ! 1.) CLIPPING xi
    IF((XiMax-XiMin).GT.1.2)THEN ! two possible intersections
      XiSplit=0.5*(XiMax+XiMin)
 !     print*,'XiSplit',Xisplit
      ! first split
      XiArray(:,iClipIter)=(/XiSplit,XiMax/)
 !     print*,'xisplit,ximax',xisplit,ximax
      !read*
      BezierControlPoints2D_temp=0.
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      IF(XiMax.NE.1.0)THEN
        PlusXi=1.0+XiMax
        MinusXi=1.0-XiMax
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
              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
                                               *(PlusXi**l)*(MinusXi**(p-l))

            END DO
          END DO
        END DO
      ELSE
        BezierControlPoints2D_temp=BezierControlPoints2D
      END IF ! XiMax

!print*,'upper'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D_temp(1,p,q)&
!        ,BezierControlPoints2d_temp(2,p,q)
!  END DO
!END DO


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
            BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
                                             +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
                                             *(PlusXi**l)*(MinusXi**(p-l))
          END DO
        END DO
      END DO

!print*,'split'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D_temp2(1,p,q)&
!        ,BezierControlPoints2d_temp2(2,p,q)
!  END DO
!END DO
      tmpnClip=iClipIter+1
      tmpnXi   =nXiClip
      tmpnEta  =nEtaClip
      firstClip=.FALSE.
      CALL BezierClip(firstClip,BezierControlPoints2D_temp2,PartTrajectory,lengthPartTrajectory &
                     ,tmpnClip,tmpnXi,tmpnEta,nInterSections,iPart,SideID)

      ! second split
      XiArray(:,iClipIter)=(/XiMin,XiSplit/)
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      BezierControlPoints2D_temp=0.
      PlusXi=1.0+XiSplit
      MinusXi=1.0-XiSplit
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
            BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                             +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
                                             *(PlusXi**l)*(MinusXi**(p-l))

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

              BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
                                               +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
                                               *(PlusXi**l)*(MinusXi**(p-l))
            END DO
          END DO
        END DO
      ELSE
        BezierControlPoints2D_temp2=BezierControlPoints2D_temp
      END IF
      tmpnClip=iClipIter+1
      tmpnXi   =nXiClip
      tmpnEta  =nEtaClip
      firstClip=.FALSE.
      CALL BezierClip(firstClip,BezierControlPoints2D_temp2,PartTrajectory,lengthPartTrajectory&
                     ,tmpnClip,tmpnXi,tmpnEta,nInterSections,iPart,SideID)
      DoCheck=.FALSE.
      EXIT
    ELSE ! only one possible intersection
      XiArray(:,iClipIter)=(/XiMin,XiMax/)
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      IF(XiMax.NE.1.0)THEN
        !print*,'do it 1'
        BezierControlPoints2D_temp=0.
        PlusXi=1.0+XiMax
        MinusXi=1.0-XiMax
        DO q=0,NGeo
          DO p=0,NGeo
            DO l=0,p
              !BezierControlPoints2D_temp(:,p,q)=&
              !BezierControlPoints2D_temp(:,p,q)+&
              !!BezierControlPoints2D(:,l,q)*B(p,l,Smax)
              !BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
              !                                 *arrayNchooseK(p,l) &
              !                                 *(1+XiMax)**l        &
              !                                 *(1-XiMax)**(p-l)
              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
                                               *(PlusXi**l)*(MinusXi**(p-l))
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
             BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)             &
                                              +BezierControlPoints2D  (:,NGeo-l,q)*FacNchooseK(p,l) &
                                              *(PlusXi**l)*(MinusXi**(p-l))
            END DO
          END DO
        END DO
        BezierControlPoints2D=BezierControlPoints2D_temp
      END IF

!print*,'new patch in xi'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D(1,p,q),BezierControlPoints2d(2,p,q)
!  END DO
!END DO
!
!
!      print*,'check xi-tolereance'  
      ! c) check Tolerance
      !IF(DoCheck)THEN
       ! check via mean value
      !x=SUM(BezierControlPoints2D(1,:,:))*PatchDOF2D
      !y=SUM(BezierControlPoints2D(2,:,:))*PatchDOF2D
      !IF(SQRT(x*x+y*y).LT.ClipTolerance)EXIT
      ! check via distance
      ZeroDistance=0.
      DO q=0,NGeo                                                                                   
        DO p=0,NGeo
          ZeroDistance=ZeroDistance+BezierControlPoints2D(1,p,q)*BezierControlPoints2d(1,p,q) &
                                   +BezierControlPoints2D(2,p,q)*BezierControlPoints2d(2,p,q)
        END DO
      END DO
      ZeroDistance=ZeroDistance*PatchDOF2D
      IF(ZeroDistance.LT.ClipTolerance2) EXIT

      !END IF ! DoCheck
      IF(ABS(XiMax-XiMin).LT.ClipTolerance) DoXiClip=.FALSE.




    END IF ! check clip size
  END IF!DoXiClip

  ! b) eta-direction
  IF(DoEtaClip)THEN
    IF(.NOT.FirstClip)THEN
      DoXiClip=.TRUE.
      FirstClip=.TRUE.
    END IF
    CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,0,NGeo,DoCheck)
    IF(.NOT.DoCheck) EXIT
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec)
      END DO
    END DO
    DO l=0,NGeo
      minmax(2,l)=MAXVAL(BezierControlPoints1D(:,l)) 
      minmax(1,l)=MINVAL(BezierControlPoints1D(:,l)) 
    END DO ! l
    ! calc Smin and Smax and check boundaries
    CALL CalcSminSmax(minmax,Etamin,Etamax)
!    print*,''
!    print*,'EtaMinm,EtaMax',EtaMin,EtaMax
!IF(time.GT.28.07)THEN
!  IF(SideID.EQ.8)THEN
!    print*,'EtaMinm,EtaMax',EtaMin,EtaMax
!    read*
!  END IF
!END IF

    IF((EtaMin.EQ.1.5).OR.(EtaMax.EQ.-1.5))RETURN
    nEtaClip=nEtaClip+1
    ! 2.) CLIPPING eta
    IF((EtaMax-EtaMin).GT.1.2)THEN ! two possible intersections
!      print*,'eta split'
      EtaSplit=0.5*(EtaMax+EtaMin)
      ! first clip
      EtaArray(:,iClipIter)=(/EtaSplit,EtaMax/)
!      print*,'etasplit,etamax',etasplit,etamax
      !read*
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      IF(EtaMax.NE.1.0)THEN
        BezierControlPoints2D_temp=0.
        PlusEta=1.0+EtaMax
        MinusEta=1.0-EtaMax
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
              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
                                               *(PlusEta**l)*(MinusEta**(p-l))

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
            BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                             +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
                                             *(PlusEta**l)*(MinusEta**(p-l))
          END DO
        END DO
      END DO
      ! new bezier-clip
      tmpnClip=iClipIter+1
      tmpnXi   =nXiClip
      tmpnEta  =nEtaClip
      firstClip=.TRUE.
      CALL BezierClip(firstClip,BezierControlPoints2D_temp2,PartTrajectory,lengthPartTrajectory &
                     ,tmpnClip,tmpnXi,tmpnEta,nInterSections,iPart,SideID)
      ! second split
!      print*,'etamin,etasplit',etamin,etasplit
      EtaArray(:,iClipIter)=(/EtaMin,EtaSplit/)
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      BezierControlPoints2D_temp=0.
      PlusEta=1.0+EtaSplit
      MinusEta=1.0-EtaSplit
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

            BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                             +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
                                             *(PlusEta**l)*(MinusEta**(p-l))
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

              BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                               +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
                                               *(PlusEta**l)*(MinusEta**(p-l))
            END DO
          END DO
        END DO
      ELSE
        BezierControlPoints2D_temp2=BezierControlPoints2D_temp
      END IF
      ! new bezier-clip
      tmpnClip=iClipIter+1
      tmpnXi   =nXiClip
      tmpnEta  =nEtaClip
      firstClip=.TRUE.
      CALL BezierClip(firstClip,BezierControlPoints2D_temp2,PartTrajectory,lengthPartTrajectory &
                     ,tmpnClip,tmpnXi,tmpnEta,nInterSections,iPart,SideID)
      DoCheck=.FALSE.
      EXIT
    ELSE ! only one possible clip in eta direction
      EtaArray(:,iClipIter)=(/EtaMin,EtaMax/)
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      IF(EtaMax.NE.1.0)THEN
        BezierControlPoints2D_temp=0.
        PlusEta=1.0+EtaMax
        MinusEta=1.0-EtaMax
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

              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
                                               *(PlusEta**l)*(MinusEta**(p-l))
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

              BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)        &
                                               +BezierControlPoints2D(:,q,NGeo-l)*FacNchooseK(p,l) &
                                               *(PlusEta**l)*(MinusEta**(p-l))
            END DO
          END DO
        END DO
        BezierControlPoints2D=BezierControlPoints2D_temp

!print*,'new patch in eta'
!DO q=0,NGeo
!  DO p=0,NGeo
!    WRITE(*,'(A,2x,I2.2,2x,I2.2,2x,F12.8,2x,F12.8)') 'Bezier2D',p,q,BezierControlPoints2D(1,p,q),BezierControlPoints2d(2,p,q)
!  END DO
!END DO

       !print*,'check eta-tolereance'  
       ! c) check Tolerance
       !IF(DoCheck)THEN

       ! check via mean value
       !x=SUM(BezierControlPoints2D(1,:,:))*PatchDOF2D
       !y=SUM(BezierControlPoints2D(2,:,:))*PatchDOF2D
       !IF(SQRT(x*x+y*y).LT.ClipTolerance)EXIT
       ! check via distance
       ZeroDistance=0.
       DO q=0,NGeo
         DO p=0,NGeo
           ZeroDistance=ZeroDistance+BezierControlPoints2D(1,p,q)*BezierControlPoints2d(1,p,q) &
                                    +BezierControlPoints2D(2,p,q)*BezierControlPoints2d(2,p,q)
         END DO
       END DO
       ZeroDistance=ZeroDistance*PatchDOF2D
       IF(ZeroDistance.LT.ClipTolerance2) EXIT


       
       !END IF ! DoCheck
       IF(ABS(EtaMax-EtaMin).LT.ClipTolerance) DoEtaClip=.FALSE.

      END IF
    END IF ! check clip size
  END IF ! DoEtaClip 
END DO

IF(DoCheck)THEN
  ! back transformation of sub-level clipping values to original bezier surface: ximean, etamean
  !   xi-direction
  !print*,'nXiClip',nXiClip
  !print*,'nEtaClip',nEtaClip
  IF(nXiClip.EQ.0)THEN
    Xi=0.
  ELSE
    Xi=0.5*SUM(XiArray(:,nXiClip))
    !print*,'deepest xi',xi
    DO iClip=nXiClip-1,1,-1
      Xi=XiArray(1,iClip)+0.5*(Xi+1)*(XiArray(2,iClip)-XiArray(1,iClip))
      !print*,'xi',xi
    END DO
  END IF ! nXIClip
  !   eta-direction
  IF(nEtaClip.EQ.0)THEN
    Eta=0.
  ELSE
    Eta=0.5*SUM(EtaArray(:,nEtaClip))
    !print*,'deepesst eta',eta
    DO iClip=nEtaClip-1,1,-1
      Eta=EtaArray(1,iClip)+0.5*(Eta+1)*(EtaArray(2,iClip)-EtaArray(1,iClip))
      !print*,eta
    END DO
  END IF ! nEtaclip
  !print*,'xi,eta',xi,eta
  ! Calculate intersection value in 3D (De Casteljau)
  tmpXi=XI
  tmpEta=Eta
  Xi=0.5*(Xi+1)
  Eta=0.5*(Eta+1)
  MinusXi =1.0-Xi
  MinusEta=1.0-Eta
  
  ReducedBezierControlPoints=BezierControlPoints3D(:,:,:,SideID)
  l=NGeo-1
  DO iDeCasteljau=1,NGeo
    !print*,'l',l
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
  
  
  !print*,'FoundPoint',ReducedBezierControlPoints(:,0,0)
  ! resulting point is ReducedBezierControlPoints(:,1,1)
  IntersectionVector=ReducedBezierControlPoints(:,0,0)-LastPartPos(iPart,1:3)
  alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)
  !print*,alpha
!IF(time.GT.28.07)THEN
!  print*,'here'
!  IF(SideID.EQ.8)THEN
!    print*,'Alpha',alpha
!    read*
!  END IF
!END IF

  
  
  !IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
  IF((alpha.LT.lengthPartTrajectory).AND.(alpha.GT.Mepsilontol))THEN
    ! found additional intersection point
    nInterSections=nIntersections+1
    IF(nInterSections.GT.ClipMaxInter)THEN
      nInterSections=0
      locAlpha=-1
      RETURN
    END IF
    locAlpha(nInterSections)=alpha
    locXi (nInterSections)=tmpXi
    locEta(nInterSections)=tmpEta
    ! defined in [0,1]
    !locXi (nInterSections)=Xi
    !locEta(nInterSections)=Eta
  END IF

END IF ! docheck

END SUBROUTINE BezierClip

SUBROUTINE calcLineNormVec(BezierControlPoints2D,LineNormVec,a,b,isParallel)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
!================================================================================================================================
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
LOGICAL,INTENT(INOUT)                :: isParallel
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length
!================================================================================================================================
LineNormVec=(BezierControlPoints2D(:,a,b)-BezierControlPoints2D(:,0,0))+&
            (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,b,a))
Length=SQRT(DOT_PRODUCT(LineNormVec,LineNormVec))
IF(Length.EQ.0)THEN
  isParallel=.FALSE.
  RETURN
END IF
LineNormVec=LineNormVec/Length
END SUBROUTINE calcLineNormVec

SUBROUTINE CalcSminSmax(minmax,Smin,Smax)
!================================================================================================================================
! find upper and lower intersection with convex hull (or no intersection)
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol
USE MOD_Mesh_Vars,               ONLY:NGeo,Xi_NGeo,DeltaXi_NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: minmax(1:2,0:NGeo)
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
    IF(minmax(2,l)*minmax(2,l+1).LE.0.)THEN
      m    = (minmax(2,l+1)-minmax(2,l))/DeltaXi_NGeo
      tmp  = Xi_NGeo(l)-minmax(2,l)/m
      Smin = MIN(tmp,Smin)
    END IF
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
END SUBROUTINE calcSminSmax

FUNCTION InsideBoundingBox(ParticlePosition,SideID)
!================================================================================================================================
! check is the particles is inside the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,BiLinearCoeff
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SlabNormals,SlabIntervalls,BezierControlPoints3D
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
x=DOT_PRODUCT(P,SlabNormals(:,1,SideID))
y=DOT_PRODUCT(P,SlabNormals(:,2,SideID))
z=DOT_PRODUCT(P,SlabNormals(:,3,SideID))
IF((x.LT.SlabIntervalls(1,SideID)-epsilontol).OR.(x.GT.SlabIntervalls(2,SideID)+epsilontol))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
IF((y.LT.SlabIntervalls(3,SideID)-epsilontol).OR.(y.GT.SlabIntervalls(4,SideID)+epsilontol))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
IF((z.LT.SlabIntervalls(5,SideID)-epsilontol).OR.(z.GT.SlabIntervalls(6,SideID)+epsilontol))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
InsideBoundingBox=.TRUE.
END FUNCTION InsideBoundingBox

FUNCTION BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol,BiLinearCoeff
USE MOD_Particle_Vars,            ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:SlabNormals,SlabIntervalls,BezierControlPoints3D
USE MOD_TimeDisc_Vars,               ONLY:iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: lengthPartTrajectory
!REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: iPart,SideID
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
! 1.) Calculate the projection of the PartTrajectory onto the SlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
DO i=1,3!x,y,z direction
  !dnk=DOT_PRODUCT(PartTrajectory,SlabNormals(i,:,SideID))
  dnk=DOT_PRODUCT(PartTrajectory,SlabNormals(:,i,SideID))
  IF(ABS(dnk).LT.epsilontol)THEN
    dnk=epsilontol ! ÜBERPRÜFEN OB SIGN sinn macht
  END IF
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SlabNormals(:,i,SideID))&
                                                                              +SlabIntervalls(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SlabNormals(:,i,SideID))&
                                                                              +SlabIntervalls(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SlabNormals(:,i,SideID))&
                                                                              +SlabIntervalls(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SlabNormals(:,i,SideID))&
                                                                              +SlabIntervalls(2*i  ,SideID) )/dnk!t_max
  END IF
END DO!i
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAXVAL(alpha(1,:))
minvalue=MINVAL(alpha(2,:))


!IF(iter.GT.430)THEN
!  !print*,'SideID',SideID
!  IF(SideID.eq.6)THEN
!    print*,'sideid',sideid
!    print*,'max',maxvalue
!    print*,'min',minvalue
!!    print*,'SlabNormals(:,1,SideID)',SlabNormals(:,1,SideID)
!!    print*,'SlabNormals(:,2,SideID)',SlabNormals(:,2,SideID)
!!    print*,'SlabNormals(:,3,SideID)',SlabNormals(:,3,SideID)
!!    print*,'SlabInt1',SlabIntervalls(1:2,SideID)
!!    print*,'SlabInt2',SlabIntervalls(3:4,SideID)
!!    print*,'SlabInt3',SlabIntervalls(5:6,SideID)
!    print*,'length', LengthPartTrajectory
!    print*,'Lastpartpos',LastPartPos(iPart,:)
!    print*,'traj',Parttrajectory
!    read*
!  END IF
!END IF

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
  IF((maxvalue.LT.lengthPartTrajectory+epsilontol).AND.(maxvalue+epsilontol.GT.0.))THEN
  !the first intersection is less than lengthPartTrajectory and greater 0
    BoundingBoxIntersection=.TRUE.
  ELSE
    BoundingBoxIntersection=.FALSE.
  END IF
ELSE
  BoundingBoxIntersection=.FALSE.
END IF
END FUNCTION BoundingBoxIntersection

SUBROUTINE ParticleTracking()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals,                     ONLY:abort
USE MOD_Mesh_Vars,                   ONLY:ElemToSide,nBCSides
USE MOD_Particle_Vars,               ONLY:PEM,PDM
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,      ONLY:epsilontol,SideType,epsilonOne,neighborElemID,neighborlocSideID,epsilonbilinear
USE MOD_Particle_Surfaces_Vars,      ONLY:nPartCurved, SuperSampledNodes,nQuads
USE MOD_TimeDisc_Vars,               ONLY:iter
!USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteractionSuperSampled
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID
INTEGER                       :: ilocSide,SideID,flip
INTEGER                       :: iInterSect,nInter
INTEGER                       :: p,q,QuadID,iQuad,minQuadID,maxQuadID
LOGICAL                       :: PartisDone,dolocSide(1:6)
REAL                          :: alpha,xi,eta
REAL                          :: alpha_loc(1:nQuads),xi_loc(1:nQuads),eta_loc(1:nQuads)
REAL                          :: xNodes(1:3,4),Displacement,xdisplace(1:3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
!===================================================================================================================================

!IF(iter.EQ.153) read*

!print*,'ici'
!read*
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             +PartTrajectory(2)*PartTrajectory(2) &
                             +PartTrajectory(3)*PartTrajectory(3) )
    PartTrajectory=PartTrajectory/lengthPartTrajectory
    lengthPartTrajectory=lengthPartTrajectory+epsilontol
    ! print*,'ElemID','new RK',ElemID
    ! print*,'partpos',PartState(iPart,1:3)
    ! print*,'par vel',PartState(iPart,4:6)
    ! print*,'lastpos',LastPartPos(iPart,1:3)
    ! print*,'PartTrajectory',PartTrajectory
    ! read*
    ! track particle vector until the final particle position is achieved
    dolocSide=.TRUE.
    DO WHILE (.NOT.PartisDone)
      ! debugggg
      DO ilocSide=1,6
      !ilocSide=6
        alpha_loc=-1.0
        IF(.NOT.dolocSide(ilocSide)) CYCLE
        SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
        flip  =ElemToSide(E2S_FLIP,ilocSide,ElemID)
        QuadID=0
        ! supersampling of each side
        DO q=0,NPartCurved-1
          DO p=0,NPartCurved-1
            QuadID=QuadID+1
            xNodes(:,1)=SuperSampledNodes(1:3,p  ,q  ,SideID)
            xNodes(:,2)=SuperSampledNodes(1:3,p+1,q  ,SideID)
            xNodes(:,3)=SuperSampledNodes(1:3,p+1,q+1,SideID)
            xNodes(:,4)=SuperSampledNodes(1:3,p  ,q+1,SideID)
            ! compute displacement || decision between planar or bi-linear plane 
            xdisplace(1:3) = xNodes(:,1)-xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
            Displacement = xdisplace(1)*xdisplace(1)+xdisplace(2)*xdisplace(2)+xdisplace(3)*xdisplace(3)
            !print*,displacement
            IF(Displacement.LT.epsilonbilinear)THEN
              CALL ComputePlanarIntersectionSuperSampled2(xNodes,PartTrajectory,lengthPartTrajectory &
                                                         ,alpha_loc(QuadID),xi_loc(QuadID),eta_loc(QuadID),iPart,SideID)
            ELSE
!           print*, CALL abort(__STAMP__,&
!                ' flip missing!!! ',999,999.)
              !print*,'locSideID,QuadID',ilocSide,QuadID
              !print*,'p,q',p,q
              !print*,xNodes(:,1)
              !print*,xNodes(:,2)
              !print*,xNodes(:,3)
              !print*,xNodes(:,4)
              !read*
!              CALL ComputeBiLinearIntersectionSuperSampled(xNodes,PartTrajectory &
!                                                          ,alpha_loc(QuadID),xi_loc(QuadID),eta_loc(QuadID),iPart,SideID)

              CALL ComputeBiLinearIntersectionSuperSampled2(xNodes,PartTrajectory,lengthPartTrajectory &
                                                          ,alpha_loc(QuadID),xi_loc(QuadID),eta_loc(QuadID),iPart,SideID)
              !print*,'alpha_loc',alpha_loc(QuadID)
              !print*,'xi_loc',xi_loc(QuadID)
              !print*,'eta_loc',eta_loc(QuadID)
              !read*
            END IF
          END DO ! p
        END DO ! q
        ! get correct intersection
        IF(SideID.LE.nBCSides)THEN
          alpha=HUGE(0.0)
          minQuadID=99999
          ! get smallest alpha
          DO iQuad=1,nQuads
            IF(alpha_loc(iQuad).GT.epsilontol)THEN
              IF(alpha.GT.alpha_loc(iQuad))THEN
                alpha=alpha_loc(iQuad)
                minQuadID=iQuad
              END IF ! alpha.GT.alpha_loc
            END IF ! alpha_loc.GT.espilontol
          END DO ! iQuad
          ! check if interesction is possible and take first intersection
          !print*,alpha,minQuadID
          !read*
          IF(alpha.GT.epsilontol.AND.alpha.LT.lengthPartTrajectory)THEN
            !print*,'alpha',alpha,xi_loc(minQUadID),eta_loc(minQuadID)
!            xi=xi_loc(minQuadID)
!            eta=eta_loc(minQuadID)
!            QuadID=minQuadID
            !print*,'Boundary interaction implemented for new method'
            !print*,'Side',SideID
            !print*,'oldstate',PartState(iPart,1:3)
            ! scaling alpha and Part Trajectory
!            alpha=alpha/LengthPartTrajectory
!            PartTrajectory=PartTrajectory*lengthPartTrajectory
            CALL GetBoundaryInteractionSuperSampled(PartTrajectory,lengthPartTrajectory,alpha,xi_loc(minQuadID),eta_loc(minQuadID),&
                                                                                            iPart,minQuadID,SideID,ElemID)
!            CALL abort(__STAMP__,&
!                ' Boundary interaction not implemented for new method.',999,999.)
            !print*,'newState',PartState(iPart,1:3)
            !print*,'newTrajectory',PartTrajectory
!            read*
             EXIT
          ELSE ! no intersection
            alpha=-1.0
          END IF
        ELSE ! no BC Side
!          print*,'alpha_loc',alpha_loc
!          print*,'xi_loc',xi_loc
!          print*,'eta_loc',eta_loc
          ! search max alpha
          alpha=-1.0
          maxQuadID=-1
          nInter=0
          ! get largest possible intersection
          DO iQuad=1,nQuads
            !print*,'alpha_loc',alpha_loc(iQuad)
            IF(alpha_loc(iQuad).GT.alpha)THEN
              IF(alpha_loc(iQuad)*alpha.LT.0) nInter=nInter+1
              IF(ABS(alpha_loc(iQuad)/alpha).GT.epsilonOne) nInter=nInter+1
              alpha=alpha_loc(iQuad)
              maxQuadID=iQuad
            END IF
          END DO ! iQuad
          !print*,'nInter',nInter
          IF(MOD(nInter,2).EQ.0) alpha=-1.0
          IF(alpha.GT.epsilontol)THEN
          !  print*,'alpha',alpha
          !  read*
!             print*,'next elem'
!             print*,'alpha',alpha
             xi=xi_loc(maxQuadID) 
             eta=eta_loc(maxQuadID)
             ! check if the found alpha statisfy the selection condition
             iInterSect=INT((ABS(xi)-2*epsilontol)/1.0)+INT((ABS(eta)-2*epsilontol)/1.0)
             IF(iInterSect.GT.0)THEN
               CALL abort(__STAMP__,&
                   ' Particle went through edge or node. Not implemented yet.',999,999.)
             ELSE
               dolocSide=.TRUE.
               dolocSide(neighborlocSideID(ilocSide,ElemID))=.FALSE.
          !     print*,'old elem'
               ElemID=neighborElemID(ilocSide,ElemID)
          !     print*,'new elem'
!               print*,'new elem id',ElemID
               !print*,'new particle positon',ElemID
  !             CALL abort(__STAMP__,&
  !                 ' Particle mapping to neighbor elem not verified!',999,999.)
               EXIT
             END IF ! possible intersect
           ELSE
             alpha=-1.0
           END IF ! alpha.GT.epsilontol
        END IF ! SideID.LT.nBCSides
      ! debugg 
      END DO ! ilocSide
      ! no intersection found
      IF(alpha.EQ.-1.0)THEN
        PEM%Element(iPart) = ElemID
        PartisDone=.TRUE.
      END IF
    END DO ! PartisDone=.FALSE.
  END IF ! Part inside
END DO ! iPart

END SUBROUTINE ParticleTracking


!SUBROUTINE ParticleTrackinglin()
!===================================================================================================================================
!! read required parameters
!===================================================================================================================================
!! MODULES
!USE MOD_Globals,                     ONLY:abort
!USE MOD_Mesh_Vars,                   ONLY:ElemToSide,nBCSides
!USE MOD_Particle_Vars,               ONLY:PEM,PDM
!USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
!USE MOD_Particle_Surfaces_Vars,      ONLY:epsilontol,SideIsPlanar,epsilonOne,neighborElemID,neighborlocSideID
!USE MOD_Particle_Boundary_Condition, ONLY:GetBoundaryInteraction
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                       :: iPart,ElemID
!INTEGER                       :: ilocSide,SideID
!INTEGER                       :: iInterSect
!LOGICAL                       :: PartisDone,dolocSide(1:6)
!!REAL                          :: alpha(1:6),xietaIntersect(1:2,1:6)
!REAL                          :: alpha,xi,eta!xietaIntersect(1:2,1:6)
!REAL                          :: PartTrajectory(1:3)
!===================================================================================================================================
!
!DO iPart=1,PDM%ParticleVecLength
!  IF(PDM%ParticleInside(iPart))THEN
!    PartisDone=.FALSE.
!    ElemID = PEM%lastElement(iPart)
!   !Element = PEM%lastElement(i)
!    print*,ElemID
!    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
!    ! track particle vector until the final particle position is achieved
!    alpha=-1.
!    dolocSide=.TRUE.
!    DO WHILE (.NOT.PartisDone)
!      DO ilocSide=1,6
!        IF(.NOT.dolocSide(ilocSide)) CYCLE
!        SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
!        IF(SideIsPlanar(SideID))THEN
!          !CALL ComputePlanarIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xietaIntersect(1,ilocSide) &
!          !                              ,XiEtaIntersect(2,ilocSide))
!          CALL ComputePlanarIntersection(PartTrajectory,alpha,xi,eta,iPart,SideID)
!        ELSE
!          CALL ComputeBiLinearIntersection(PartTrajectory,alpha,xi,eta,iPart,SideID)
!          !CALL ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xietaIntersect(1,ilocSide) &
!          !                              ,XiEtaIntersect(2,ilocSide))
!        END IF
!        !print*,ilocSide,alpha
!        !print*,'ilocSide,alpha,xi,eta',ilocSide,alpha,xi,eta
!        !print*,'neighborElemID',neighborElemID(ilocSide,ElemID)
!        ! check after each side if particle went through checked side
!        IF(alpha.GT.epsilontol)THEN ! or minus epsilontol
!          !IF(alpha+epsilontol.GE.epsilonOne) PartisDone=.TRUE.
!          IF(SideID.LE.nBCSides)THEN
!            print*,'Boundary interaction implemented for new method'
!            CALL GetBoundaryInteraction(PartTrajectory,alpha,xi,eta,iPart,SideID,ElemID)
!            !CALL abort(__STAMP__,&
!                !' Boundary interaction not implemented for new method.',999,999.)
!          END IF
!          iInterSect=INT((ABS(xi)-epsilontol)/1.0)+INT((ABS(eta)-epsilontol)/1.0)
!          IF(iInterSect.GT.0)THEN
!            CALL abort(__STAMP__,&
!                ' Particle went through edge or node. Not implemented yet.',999,999.)
!          ELSE
!            dolocSide=.TRUE.
!            dolocSide(neighborlocSideID(ilocSide,ElemID))=.FALSE.
!            ElemID=neighborElemID(ilocSide,ElemID)
!            !print*,'new particle positon',ElemID
!!            CALL abort(__STAMP__,&
!!                ' Particle mapping to neighbor elem not verified!',999,999.)
!            EXIT
!          END IF ! iInteSect
!        END IF
!      END DO ! ilocSide
!      !sop
!      !read*
!      ! no intersection found
!      IF(alpha.EQ.-1.0)THEN
!        PEM%Element(iPart) = ElemID
!        PartisDone=.TRUE.
!      END IF
!    END DO ! PartisDone=.FALSE.
!  END IF ! Part inside
!END DO ! iPart
!
!END SUBROUTINE ParticleTracking

!SUBROUTINE ComputePlanarIntersection(PartTrajectory,alpha,xi,eta,iPart,SideID)
!==================================================================================================================================
!! Compute the Intersection with planar surface
!==================================================================================================================================
!! MODULES
!USE MOD_Particle_Vars,           ONLY:LastPartPos
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,SideDistance,epsilonOne
!!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
!INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)                  :: alpha,xi,eta
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                              :: coeffA,coeffB,xInter(3)
!REAL                              :: Axz, Bxz, Cxz
!REAL                              :: Ayz, Byz, Cyz
!==================================================================================================================================
!
!! set alpha to minus 1, asume no intersection
!alpha=-1.0
!xi=-2.
!eta=-2.
!
!! check if the particle can intersect with the planar plane
!! if the normVec point in the opposite direction, cycle
!coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!!print*,'coeffA',coeffA
!IF(ABS(coeffA).LT.+epsilontol)THEN
!  ! particle tangential to surface ==> no interesection, particle remains in element
!  RETURN
!END IF
!
!! distance of plane fromn origion minus trajectory start point times normal vector of side
!coeffB=SideDistance(SideID)-DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(iPart,1:3))
!
!alpha=coeffB/coeffA
!!print*,'coeffB',coeffB
!!print*,'alpha',alpha
!!read*
!
!IF((alpha.GT.epsilonOne).OR.(alpha.LT.epsilontol))THEN
!  alpha=-1.0
!  RETURN
!END IF
!
!! compute intersection
!xInter(1:3) =LastPartPos(iPart,1:3)+alpha*PartTrajectory(1:3)
!
!!! theoretically, can be computed in advance
!Axz = BiLinearCoeff(1,2,SideID) - BiLinearCoeff(3,2,SideID)
!Bxz = BiLinearCoeff(1,3,SideID) - BiLinearCoeff(3,3,SideID)
!Cxz = xInter(1) - BiLinearCoeff(1,4,SideID) - xInter(3) + BiLinearCoeff(3,4,SideID)
!
!Ayz = BiLinearCoeff(2,2,SideID) - BiLinearCoeff(3,2,SideID)
!Byz = BiLinearCoeff(2,3,SideID) - BiLinearCoeff(3,3,SideID)
!Cyz = xInter(2) - BiLinearCoeff(2,4,SideID) - xInter(3) + BiLinearCoeff(3,4,SideID)
!
!print*,'Bxz,Byz',Bxz,Byz
!
!IF(ABS(Bxz).LT.epsilontol)THEN
!  xi = Axz + Bxz*Ayz/Byz
!  ! check denominator
!  xi = (Cxz - Bxz*Ayz/Byz*Cyz)/xi
!ELSE
!  xi = Ayz + Byz*Axz/Bxz
!  ! check denominator
!  xi = (Cyz - Byz*Axz/Bxz*Cxz)/xi
!END IF
!
!IF(ABS(xi).GT.epsilonOne) THEN 
!  ! xi outside of possible range
!  alpha=-1.0
!  RETURN
!END IF
!
!eta = Bxz+Byz
!eta = (Cxz+Cyz - (Axz+Ayz)*xi) / eta
!
!!print*,'xi,eta',xi,eta
!IF(ABS(eta).GT.epsilonOne) THEN 
!  ! eta outside of possible range
!  alpha=-1.0
!  RETURN
!END IF
!
!! here, eta,xi,alpha are computed
!
!END SUBROUTINE ComputePlanarIntersection

SUBROUTINE ComputePlanarIntersectionBezier(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals,                 ONLY:Cross,abort
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne,SideDistance
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
!USE MOD_Equations_Vars,          ONLY:epsMach
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,nlength,locSideDistance
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.

coeffA=DOT_PRODUCT(SideNormVec(:,SideID),PartTrajectory)

!! corresponding to particle starting in plane
!! interaction should be computed in last step
IF(ABS(coeffA).LT.+epsilontol)THEN 
  RETURN
END IF

locSideDistance=SideDistance(SideID)-DOT_PRODUCT(LastPartPos(iPart,1:3),SideNormVec(:,SideID))
alpha=locSideDistance/coeffA

!IF(SideID.EQ.8)THEN
!  print*,'nvec',SideNormVec(:,SideID)
!  print*,'LastPartPos,SideDistance',DOT_PRODUCT(LastPartPos(iPart,1:3),SideNormVec(:,sideID)),SideDistance(SideID)
!  print*,'alppha',alpha
!  print*,'len',lengthPartTrajectory
!  IF(alpha.LT.0) STOP
!END IF

!IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
IF((alpha.GT.lengthPartTrajectory).OR.(alpha.LT.-epsilontol))THEN
  alpha=-1.0
  RETURN
END IF

P1=0.25*(-BezierControlPoints3D(:,0,0   ,SideID)+BezierControlPoints3D(:,NGeo,0,SideID)      &
         -BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )

P2=0.25*(-BezierControlPoints3D(:,0,0,SideID)   -BezierControlPoints3D(:,NGeo,0,SideID)        &
         +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) )

P0=0.25*(BezierControlPoints3D(:,0,0,SideID)+BezierControlPoints3D(:,NGeo,0,SideID)         &
        +BezierControlPoints3D(:,0,NGeo,SideID)+BezierControlPoints3D(:,NGeo,NGeo,SideID) ) &
        -LastPartPos(iPart,1:3)-alpha*PartTrajectory


A1=P1(1)+P1(3)
B1=P2(1)+P2(3)
C1=P0(1)+P0(3)

A2=P1(2)+P1(3)
B2=P2(2)+P2(3)
C2=P0(2)+P0(3)

!IF(SideID.EQ.8)THEN
!  print*,'alppha',alpha
!  print*,'A1',A1
!  print*,'B1',B1
!  print*,'C1',C1
!  
!  print*,'A2',A2
!  print*,'B2',B2
!  print*,'C2',C2
!
!  print*,"A2-B2/B2*A2", A2-B2/B1*A1
!  print*,'B2/B1*C1-C2', B2/B1*C1-C2
!  print*,'gesamt',(B2/B1*C1-C2)/(A2-B2/B1*A1)
!END IF


IF(ABS(B1).GT.epsilontol)THEN
  xi = A2-B2/B1*A1
  IF(ABS(xi).LT.epsilontol)THEN
    xi=0.
  ELSE
    xi = (B2/B1*C1-C2)/xi
  END IF
ELSE
  xi = A1-B1/B2*A2
  IF(ABS(xi).LT.epsilontol)THEN
    xi=0.
  ELSE
    xi = (B1/B2*C2-C1)/xi
  END IF
END IF

!IF(SideID.EQ.8)THEN
!  print*,'xi',xi
!END IF

IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

eta=-((A1+A2)*xi+C1+C2)/(B1+B2)

!IF(SideID.EQ.8)THEN
!  print*,'eta',eta
!END IF

IF(ABS(eta).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

END SUBROUTINE ComputePlanarIntersectionBezier


!SUBROUTINE ComputePlanarIntersectionSuperSampled(xNodes,PartTrajectory,alpha,xi,eta,iPart)
SUBROUTINE ComputePlanarIntersectionSuperSampled2(xNodes,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals,                 ONLY:Cross,abort
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne,SideType
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2,nVec
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,coeffB,nlength
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.

! compute basis vectors of plane
! first vector of plane
P1 = -xNodes(:,1)+xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
! second vector
P2 = -xNodes(:,1)-xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
! base point
P0 = xNodes(:,1)+xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
P1=0.25*P1
P2=0.25*P2
! distance of plane fromn origion minus trajectory start point times normal vector of side
!P0=P0-LastPartPos(iPart,1:3)
P0=0.25*P0-LastPartPos(iPart,1:3)
! planar plane
! P1*xi + P2*eta+P0

nVec=CROSS(P1,P2)
nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2) +nVec(3)*nVec(3) 
nlength=SQRT(nlength)
nVec=nVec/nlength
!IF(flip.NE.0)THEN
!  nVec=-nVec
!END IF

 
!! compute distance along trajectory
coeffA=DOT_PRODUCT(nVec,PartTrajectory)
!print*,coeffA
!IF(flip.EQ.0)THEN ! master side ! is in flip for normVec
!  IF(coeffA.LT.epsilontol)RETURN
!ELSE ! slave sides
!  IF(coeffA.GT.-epsilontol)RETURN
!END IF ! flip

!! corresponding to particle starting in plane
!! interaction should be computed in last step
IF(ABS(coeffA).LT.+epsilontol)THEN 
  RETURN
END IF
coeffB=DOT_PRODUCT(P0,nVec)

alpha=coeffB/coeffA
!print*,'coeffA',coeffA
!print*,'coeffB',coeffB print*,'alpha',alpha
!read*
!read*

!IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
IF((alpha.GT.lengthPartTrajectory).OR.(alpha.LT.-epsilontol))THEN
  alpha=-1.0
  RETURN
END IF

! new computation
! compute xi and eta of intersection
P0=P0-alpha*PartTrajectory


A1=P1(1)+P1(3)
B1=P2(1)+P2(3)
C1=P0(1)+P0(3)

A2=P1(2)+P1(3)
B2=P2(2)+P2(3)
C2=P0(2)+P0(3)


IF(ABS(B1).GT.epsilontol)THEN
  xi = A2-B2/B1*A1
  xi = (B2/B1*C1-C2)/xi
ELSE
  xi = A1-B1/B2*A2
  xi = (B1/B2*C2-C1)/xi
END IF

IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

eta=-((A1+A2)*xi+C1+C2)/(B1+B2)
IF(ABS(eta).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF


!!  a1(2)= P1(1)*PartTrajectory(3)-P1(3)*PartTrajectory(1)
!!  a1(3)= P2(1)*PartTrajectory(3)-P2(3)*PartTrajectory(1)
!!  a1(4)= P0(1)*PartTrajectory(3) -P0(3)*PartTrajectory(1)
!!  
!!  a2(2)= P1(2)*PartTrajectory(3)-P1(3)*PartTrajectory(2)
!!  a2(3)= P2(2)*PartTrajectory(3)-P2(3)*PartTrajectory(2)
!!  a2(4)= P0(2)*PartTrajectory(3) &
!!        -P0(3)*PartTrajectory(2)
!!  
!!  !print*,'a12',a1(2)
!!  !print*,'a13',a1(3)
!!  !print*,'a14',a1(4)
!!  
!!  
!!  !print*,'a22',a2(2)
!!  !print*,'a23',a2(3)
!!  !print*,'a24',a2(4)
!!  !print*,'a23,a13',a2(3),a1(3)
!!  !print*,'a22,a12',a2(2),a1(2)
!!  
!!  ! old one               ! working in not all cases
!!  !! caution with accuracy
!!  !IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
!!  !  eta=0.
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !  ELSE
!!  !    ! compute xi
!!  !    xi=a1(2)-a2(2)
!!  !    xi=1.0/xi
!!  !    xi=(a2(4)-a1(4))*xi
!!  !  END IF
!!  !!  IF(ABS(xi).GT.epsilonOne)THEN
!!  !!    RETURN
!!  !!  END IF
!!  !ELSE ! a2(3) not zero
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !    eta=a1(3)-a2(3)
!!  !    eta=1.0/eta
!!  !    eta=(a2(4)-a1(4))*eta
!!  !  ELSE
!!  !    xi = a1(2) - a1(3)*a2(2)/a2(3)
!!  !    xi = 1.0/xi
!!  !    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !    ! check distance of xi 
!!  !  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !  !    RETURN
!!  !  !  END IF
!!  !    ! compute eta
!!  !    eta=a1(3)-a2(3)
!!  !    eta=1.0/eta
!!  !    eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !  END IF
!!  !END IF
!!  
!!  IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
!!    ! solution independent of eta
!!    !eta=0.
!!    IF(ABS(a2(2)).LT.epsilontol)THEN
!!      ! and independent of xi
!!      ! particle in plane
!!      !xi=99.
!!      !eta=99.
!!      print*,'here'
!!      xi=0.
!!    ELSE
!!      ! compute xi
!!      xi=a1(2)-a2(2)
!!      xi=1.0/xi
!!      xi=(a2(4)-a1(4))*xi
!!      ! compute eta
!!      IF(ABS(P2(1)).GT.epsilontol)THEN
!!        eta=(P0(1)-xi*P1(1))/P2(1)
!!      ELSE IF(ABS(P2(2)).GT.epsilontol)THEN
!!        eta=(P0(2)-xi*P1(3))/P2(2)
!!      ELSE IF(ABS(P2(3)).GT.epsilontol)THEN
!!        eta=(P0(2)-xi*P1(3))/P2(3)
!!      ELSE
!!        CALL abort(__STAMP__,&
!!                  ' error in computation of xi! iPart,eta ',iPart,eta)
!!      END IF
!!    END IF
!!  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !    RETURN
!!  !  END IF
!!  END IF
!!  IF(ABS(a2(2)).LT.epsilontol)THEN
!!    !xi=0.
!!    IF(ABS(a2(3)).LT.epsilontol)THEN
!!      ! and independent of xi
!!      ! particle in plane
!!      !xi=99.
!!      !eta=99.
!!      print*,'here'
!!      eta=0.
!!    ELSE
!!      eta=a1(3)-a2(3)
!!      eta=1.0/eta
!!      eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!      ! recompute xi
!!      IF(ABS(P1(1)).GT.epsilontol)THEN
!!        xi=(P0(1)-eta*P2(1))/P1(1)
!!      ELSE IF(ABS(P1(2)).GT.epsilontol)THEN
!!        xi=(P0(2)-eta*P2(3))/P1(2)
!!      ELSE IF(ABS(P1(3)).GT.epsilontol)THEN
!!        xi=(P0(2)-eta*P2(3))/P1(3)
!!      ELSE
!!        CALL abort(__STAMP__,&
!!                  ' error in computation of xi! iPart,eta ',iPart,eta)
!!      END IF
!!     ! ! compute xi
!!     ! xi=a1(2)-a2(2)
!!     ! xi=1.0/xi
!!     ! xi=(a2(4)-a1(4))*xi
!!    END IF
!!  END IF
!!  
!!  !ELSE ! a2(3) not zero ==> eta not zero?
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !    eta=a1(3)-a2(3)
!!  !    IF(ABS(eta).LT.epsilontol)THEN
!!  !      eta=0. ! here not 99??
!!  !    ELSE
!!  !      eta=1.0/eta
!!  !      eta=(a2(4)-a1(4))*eta
!!  !    END IF
!!  !  ELSE
!!  !    xi = a1(2) - a1(3)*a2(2)/a2(3)
!!  !    xi = 1.0/xi
!!  !    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !      ! check distance of xi 
!!  !  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !  !    RETURN
!!  !  !  END IF
!!  !    ! compute eta
!!  !    eta=a1(3)-a2(3)
!!  !    IF(ABS(eta).LT.epsilontol)THEN
!!  !      eta=0. ! here eta 99 ???
!!  !    ELSE ! eta not zero
!!  !     eta=1.0/eta
!!  !     eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !    END IF ! eta .LT.epsilontol
!!  !  END IF
!!  !END IF
!!  
!!  
!!  print*,'alpha,xi,eta',alpha,xi,eta
!!  
!!  !xi = a1(2) - a1(3)*a2(2)/a2(3) !xi = 1.0/xi
!!  !xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !! check distance of xi 
!!  IF(ABS(xi).GT.epsilonOne)THEN
!!    alpha=-1.0
!!    RETURN
!!  END IF
!!  !! compute eta
!!  !eta=a1(3)-a2(3)
!!  !eta=1.0/eta
!!  !eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !
!!  IF(ABS(eta).GT.epsilonOne)THEN
!!    alpha=-1.0
!!    RETURN
!!  END IF
!!  
!!  !! compute distance with intersection
!!  !IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!!  !  alpha =xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) -lastPartPos(iPart,1)
!!  !  alpha = alpha/ PartTrajectory(1)
!!  !ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!!  !  alpha =xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) -lastPartPos(iPart,2)
!!  !  alpha = alpha/ PartTrajectory(2)
!!  !ELSE
!!  !  alpha =xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) -lastPartPos(iPart,3)
!!  !  alpha = alpha/ PartTrajectory(3)
!!  !END IF
!!  !
!!  !IF((alpha.LT.epsilontol).OR.(alpha.GT.epsilonOne)) alpha=-1.0

END SUBROUTINE ComputePlanarIntersectionSuperSampled2

!SUBROUTINE ComputePlanarIntersectionSuperSampled(xNodes,PartTrajectory,alpha,xi,eta,iPart)
SUBROUTINE ComputePlanarIntersectionSuperSampled(xNodes,PartTrajectory,alpha,xi,eta,iPart)
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals,                 ONLY:Cross,abort
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
!REAL,INTENT(IN)                   :: lengthPartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart!,SideID!,ElemID
!INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2,nVec
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,coeffB,nlength
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.

! compute basis vectors of plane
! first vector of plane
P1 = -xNodes(:,1)+xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
! second vector
P2 = -xNodes(:,1)-xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
! base point
P0 = xNodes(:,1)+xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
P1=0.25*P1
P2=0.25*P2
! distance of plane fromn origion minus trajectory start point times normal vector of side
!P0=P0-LastPartPos(iPart,1:3)
P0=0.25*P0-LastPartPos(iPart,1:3)
! planar plane
! P1*xi + P2*eta+P0

nVec=CROSS(P1,P2)
nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2) +nVec(3)*nVec(3) 
nlength=SQRT(nlength)
nVec=nVec/nlength
!IF(flip.NE.0)THEN
!  nVec=-nVec
!END IF
 
!! compute distance along trajectory
coeffA=DOT_PRODUCT(nVec,PartTrajectory)
!print*,coeffA
!IF(flip.EQ.0)THEN ! master side ! is in flip for normVec
!  IF(coeffA.LT.epsilontol)RETURN
!ELSE ! slave sides
!  IF(coeffA.GT.-epsilontol)RETURN
!END IF ! flip

!! corresponding to particle starting in plane
!! interaction should be computed in last step
IF(ABS(coeffA).LT.+epsilontol)THEN 
  RETURN
END IF
coeffB=DOT_PRODUCT(P0,nVec)

alpha=coeffB/coeffA
!print*,'coeffA',coeffA
!print*,'coeffB',coeffB
!print*,'alpha',alpha
!read*

IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
  alpha=-1.0
  RETURN
END IF

! new computation
! compute xi and eta of intersection
P0=P0-alpha*PartTrajectory


A1=P1(1)+P1(3)
B1=P2(1)+P2(3)
C1=P0(1)+P0(3)

A2=P1(2)+P1(3)
B2=P2(2)+P2(3)
C2=P0(2)+P0(3)

IF(ABS(B1).GT.epsilontol)THEN
  xi = A2-B2/B1*A1
  xi = (B2/B1*C1-C2)/xi
ELSE
  xi = A1-B1/B2*A2
  xi = (B1/B2*C2-C1)/xi
END IF

IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

eta=-((A1+A2)*xi+C1+C2)/(B1+B2)
IF(ABS(eta).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF


!!  a1(2)= P1(1)*PartTrajectory(3)-P1(3)*PartTrajectory(1)
!!  a1(3)= P2(1)*PartTrajectory(3)-P2(3)*PartTrajectory(1)
!!  a1(4)= P0(1)*PartTrajectory(3) -P0(3)*PartTrajectory(1)
!!  
!!  a2(2)= P1(2)*PartTrajectory(3)-P1(3)*PartTrajectory(2)
!!  a2(3)= P2(2)*PartTrajectory(3)-P2(3)*PartTrajectory(2)
!!  a2(4)= P0(2)*PartTrajectory(3) &
!!        -P0(3)*PartTrajectory(2)
!!  
!!  !print*,'a12',a1(2)
!!  !print*,'a13',a1(3)
!!  !print*,'a14',a1(4)
!!  
!!  
!!  !print*,'a22',a2(2)
!!  !print*,'a23',a2(3)
!!  !print*,'a24',a2(4)
!!  !print*,'a23,a13',a2(3),a1(3)
!!  !print*,'a22,a12',a2(2),a1(2)
!!  
!!  ! old one               ! working in not all cases
!!  !! caution with accuracy
!!  !IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
!!  !  eta=0.
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !  ELSE
!!  !    ! compute xi
!!  !    xi=a1(2)-a2(2)
!!  !    xi=1.0/xi
!!  !    xi=(a2(4)-a1(4))*xi
!!  !  END IF
!!  !!  IF(ABS(xi).GT.epsilonOne)THEN
!!  !!    RETURN
!!  !!  END IF
!!  !ELSE ! a2(3) not zero
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !    eta=a1(3)-a2(3)
!!  !    eta=1.0/eta
!!  !    eta=(a2(4)-a1(4))*eta
!!  !  ELSE
!!  !    xi = a1(2) - a1(3)*a2(2)/a2(3)
!!  !    xi = 1.0/xi
!!  !    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !    ! check distance of xi 
!!  !  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !  !    RETURN
!!  !  !  END IF
!!  !    ! compute eta
!!  !    eta=a1(3)-a2(3)
!!  !    eta=1.0/eta
!!  !    eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !  END IF
!!  !END IF
!!  
!!  IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
!!    ! solution independent of eta
!!    !eta=0.
!!    IF(ABS(a2(2)).LT.epsilontol)THEN
!!      ! and independent of xi
!!      ! particle in plane
!!      !xi=99.
!!      !eta=99.
!!      print*,'here'
!!      xi=0.
!!    ELSE
!!      ! compute xi
!!      xi=a1(2)-a2(2)
!!      xi=1.0/xi
!!      xi=(a2(4)-a1(4))*xi
!!      ! compute eta
!!      IF(ABS(P2(1)).GT.epsilontol)THEN
!!        eta=(P0(1)-xi*P1(1))/P2(1)
!!      ELSE IF(ABS(P2(2)).GT.epsilontol)THEN
!!        eta=(P0(2)-xi*P1(3))/P2(2)
!!      ELSE IF(ABS(P2(3)).GT.epsilontol)THEN
!!        eta=(P0(2)-xi*P1(3))/P2(3)
!!      ELSE
!!        CALL abort(__STAMP__,&
!!                  ' error in computation of xi! iPart,eta ',iPart,eta)
!!      END IF
!!    END IF
!!  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !    RETURN
!!  !  END IF
!!  END IF
!!  IF(ABS(a2(2)).LT.epsilontol)THEN
!!    !xi=0.
!!    IF(ABS(a2(3)).LT.epsilontol)THEN
!!      ! and independent of xi
!!      ! particle in plane
!!      !xi=99.
!!      !eta=99.
!!      print*,'here'
!!      eta=0.
!!    ELSE
!!      eta=a1(3)-a2(3)
!!      eta=1.0/eta
!!      eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!      ! recompute xi
!!      IF(ABS(P1(1)).GT.epsilontol)THEN
!!        xi=(P0(1)-eta*P2(1))/P1(1)
!!      ELSE IF(ABS(P1(2)).GT.epsilontol)THEN
!!        xi=(P0(2)-eta*P2(3))/P1(2)
!!      ELSE IF(ABS(P1(3)).GT.epsilontol)THEN
!!        xi=(P0(2)-eta*P2(3))/P1(3)
!!      ELSE
!!        CALL abort(__STAMP__,&
!!                  ' error in computation of xi! iPart,eta ',iPart,eta)
!!      END IF
!!     ! ! compute xi
!!     ! xi=a1(2)-a2(2)
!!     ! xi=1.0/xi
!!     ! xi=(a2(4)-a1(4))*xi
!!    END IF
!!  END IF
!!  
!!  !ELSE ! a2(3) not zero ==> eta not zero?
!!  !  IF(ABS(a2(2)).LT.epsilontol)THEN
!!  !    xi=0.
!!  !    eta=a1(3)-a2(3)
!!  !    IF(ABS(eta).LT.epsilontol)THEN
!!  !      eta=0. ! here not 99??
!!  !    ELSE
!!  !      eta=1.0/eta
!!  !      eta=(a2(4)-a1(4))*eta
!!  !    END IF
!!  !  ELSE
!!  !    xi = a1(2) - a1(3)*a2(2)/a2(3)
!!  !    xi = 1.0/xi
!!  !    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !      ! check distance of xi 
!!  !  !  IF(ABS(xi).GT.epsilonOne)THEN
!!  !  !    RETURN
!!  !  !  END IF
!!  !    ! compute eta
!!  !    eta=a1(3)-a2(3)
!!  !    IF(ABS(eta).LT.epsilontol)THEN
!!  !      eta=0. ! here eta 99 ???
!!  !    ELSE ! eta not zero
!!  !     eta=1.0/eta
!!  !     eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !    END IF ! eta .LT.epsilontol
!!  !  END IF
!!  !END IF
!!  
!!  
!!  print*,'alpha,xi,eta',alpha,xi,eta
!!  
!!  !xi = a1(2) - a1(3)*a2(2)/a2(3) !xi = 1.0/xi
!!  !xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!!  !! check distance of xi 
!!  IF(ABS(xi).GT.epsilonOne)THEN
!!    alpha=-1.0
!!    RETURN
!!  END IF
!!  !! compute eta
!!  !eta=a1(3)-a2(3)
!!  !eta=1.0/eta
!!  !eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!!  !
!!  IF(ABS(eta).GT.epsilonOne)THEN
!!    alpha=-1.0
!!    RETURN
!!  END IF
!!  
!!  !! compute distance with intersection
!!  !IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!!  !  alpha =xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) -lastPartPos(iPart,1)
!!  !  alpha = alpha/ PartTrajectory(1)
!!  !ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!!  !  alpha =xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) -lastPartPos(iPart,2)
!!  !  alpha = alpha/ PartTrajectory(2)
!!  !ELSE
!!  !  alpha =xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) -lastPartPos(iPart,3)
!!  !  alpha = alpha/ PartTrajectory(3)
!!  !END IF
!!  !
!!  !IF((alpha.LT.epsilontol).OR.(alpha.GT.epsilonOne)) alpha=-1.0

END SUBROUTINE ComputePlanarIntersectionSuperSampled

SUBROUTINE ComputePlanarIntersection(PartTrajectory,alpha,xi,eta,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,SideDistance,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(2:4)                 :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                                :: coeffA,coeffB
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.

! compute distance of lastPartPos with planar plane

coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!print*,'coeffA',coeffA
!read*
! corresponding to particle starting in plane
! interaction should be computed in last step
IF(ABS(coeffA).LT.+epsilontol)THEN 
  RETURN
END IF
! distance of plane fromn origion minus trajectory start point times normal vector of side
coeffB=SideDistance(SideID)-DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(iPart,1:3))

alpha=coeffB/coeffA
!!print*,'coeffB',coeffB
!!print*,'alpha',alpha
!!read*

IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
  alpha=-1.0
  RETURN
END IF

a1(2)= BilinearCoeff(1,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(1)
a1(4)= (BilinearCoeff(1,4,SideID)-LastPartPos(iPart,1))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(2)= BilinearCoeff(2,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(2)
a2(4)= (BilinearCoeff(2,4,SideID)-LastPartPos(iPart,2))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(2)

!print*,'a23,a13',a2(3),a1(3)
!print*,'a22,a12',a2(2),a1(2)

!! caution with accuracy
IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
  eta=0.
  IF(ABS(a2(2)).LT.epsilontol)THEN
    xi=0.
  ELSE
    ! compute xi
    xi=a1(2)-a2(2)
    xi=1.0/xi
    xi=(a2(4)-a1(4))*xi
  END IF
!  IF(ABS(xi).GT.epsilonOne)THEN
!    RETURN
!  END IF
ELSE ! a2(3) not zero
  IF(ABS(a2(2)).LT.epsilontol)THEN
    xi=0.
    eta=a1(3)-a2(3)
    eta=1.0/eta
    eta=(a2(4)-a1(4))*eta
  ELSE
    xi = a1(2) - a1(3)*a2(2)/a2(3)
    xi = 1.0/xi
    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
    ! check distance of xi 
  !  IF(ABS(xi).GT.epsilonOne)THEN
  !    RETURN
  !  END IF
    ! compute eta
    eta=a1(3)-a2(3)
    eta=1.0/eta
    eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
  END IF
END IF

!xi = a1(2) - a1(3)*a2(2)/a2(3)
!xi = 1.0/xi
!xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!! check distance of xi 
IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF
!! compute eta
!eta=a1(3)-a2(3)
!eta=1.0/eta
!eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!
IF(ABS(eta).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

!! compute distance with intersection
!IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!  alpha =xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) -lastPartPos(iPart,1)
!  alpha = alpha/ PartTrajectory(1)
!ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!  alpha =xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) -lastPartPos(iPart,2)
!  alpha = alpha/ PartTrajectory(2)
!ELSE
!  alpha =xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) -lastPartPos(iPart,3)
!  alpha = alpha/ PartTrajectory(3)
!END IF
!
!IF((alpha.LT.epsilontol).OR.(alpha.GT.epsilonOne)) alpha=-1.0

END SUBROUTINE ComputePlanarIntersection

SUBROUTINE ComputeBiLinearIntersection(PartTrajectory,alpha,xitild,etatild,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, epsilontol,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2), q1(3)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0

a1(1)= BilinearCoeff(1,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(1)
a1(4)= (BilinearCoeff(1,4,SideID)-LastPartPos(iPart,1))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(2)
a2(4)= (BilinearCoeff(2,4,SideID)-LastPartPos(iPart,2))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(2)

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)
!print*,'A,B,C', A,B,C
CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
!print*,nRoot,Eta
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      !q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(1)=ComputeSurfaceDistance(xi(1),eta(1),PartTrajectory,iPart,SideID)
      IF((t(1).GE.+epsilontol).AND.(t(1).LE.epsilonOne))THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. epsilonOne
  ELSE ! eta not in reange
    RETURN 
  END IF ! eta .lt. epsilonOne
ELSE 
  nInter=0
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      ! q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !  WRITE(*,*) ' t ', t(2)
      !  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      t(1)=ComputeSurfaceDistance(xi(1),eta(1),PartTrajectory,iPart,SideID)
      IF((t(1).LT.epsilontol).OR.(t(1).GT.epsilonOne))THEN
        t(1)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(1).LT.epsilontol).AND.(t(1).GT.epsilonOne))THEN
!        t(1)=-1
!      END IF
    END IF
  END IF
  IF(ABS(eta(2)).LT.epsilonOne)THEN
    xi(2)=eta(2)*a2(1)-eta(2)*a1(1)+a2(2)-a1(2)
    xi(2)=1.0/xi(2)
    xi(2)=(eta(2)*a1(3)-eta(2)*a2(3)+a1(4)-a2(4))*xi(2)
    IF(ABS(xi(2)).LT.epsilonOne)THEN
      ! q1=xi(2)*eta(2)*BilinearCoeff(:,1)+xi(2)*BilinearCoeff(:,2)+eta(2)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(2)=ComputeSurfaceDistance(xi(2),eta(2),PartTrajectory,iPart,SideID)
      IF((t(2).LT.epsilontol).OR.(t(2).GT.epsilonOne))THEN
        t(2)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(2).LT.epsilontol).AND.(t(2).GT.epsilonOne))THEN
!        t(2)=-1
!      END IF
      !IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
      !!  WRITE(*,*) ' Second Intersection'
      !!  WRITE(*,*) ' t ', t(2)
      !!  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      !END IF 
    END IF
  END IF
  ! if no intersection, return
  IF(nInter.EQ.0) RETURN
  IF(SideID.LE.nBCSides)THEN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  ELSE ! no BC Side
    ! if two intersections, return, particle re-enters element
    IF(nInter.EQ.2) RETURN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  END IF ! SideID.LT.nCBSides
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersection


SUBROUTINE ComputeBiLinearIntersectionSuperSampled2(xNodes,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
!USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2), q1(3)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0

! compute initial vectors
BiLinearCoeff(:,1) = xNodes(:,1)-xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
BiLinearCoeff(:,2) =-xNodes(:,1)+xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
BiLinearCoeff(:,3) =-xNodes(:,1)-xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
BiLinearCoeff(:,4) = xNodes(:,1)+xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
BiLinearCoeff= 0.25*BiLinearCoeff

! compute product with particle trajectory
a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
a1(4)=(BilinearCoeff(1,4)-LastPartPos(iPart,1))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
a2(4)=(BilinearCoeff(2,4)-LastPartPos(iPart,2))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(iPart,3))*PartTrajectory(2)

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)
!print*,'A,B,C', A,B,C
CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
!read*
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF
!SELECT CASE(iter)
!CASE(162)
!  print*,nRoot,eta
!  print*,'PartTrajectory',PartTrajectory
!  print*,'length',LengthPartTrajectory
!  read*
!END SELECT

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      !q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      !IF((t(1).GE.-epsilontol).AND.(t(1).LE.epsilonOne))THEN
      IF((t(1).GE.-epsilontol).AND.(t(1).LE.lengthPartTrajectory))THEN
        alpha=t(1)!/LengthPartTrajectory
        xitild=xi(1)
        etatild=eta(1)
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. epsilonOne
  ELSE ! eta not in reange
    RETURN 
  END IF ! eta .lt. epsilonOne
ELSE 
  nInter=0
  t=-2.
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
     ! q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !  WRITE(*,*) ' t ', t(2)
      !  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      !t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      !IF((t(1).LT.-epsilontol).OR.(t(1).GT.epsilonOne))THEN
      IF((t(1).LT.-epsilontol).AND.(t(1).GT.lengthPartTrajectory))THEN
        t(1)=-2.0
      ELSE
        nInter=nInter+1
        t(1)=t(1)!/lengthPartTrajectory
      END IF
!      IF((t(1).LT.epsilontol).AND.(t(1).GT.epsilonOne))THEN
!        t(1)=-1
!      END IF
    END IF
  END IF
  IF(ABS(eta(2)).LT.epsilonOne)THEN
    xi(2)=eta(2)*a2(1)-eta(2)*a1(1)+a2(2)-a1(2)
    xi(2)=1.0/xi(2)
    xi(2)=(eta(2)*a1(3)-eta(2)*a2(3)+a1(4)-a2(4))*xi(2)
!    SELECT CASE(iter)
!    CASE(162)
!      print*,'xi',xi(2)
!      read*
!    END SELECT
    IF(ABS(xi(2)).LT.epsilonOne)THEN
      ! q1=xi(2)*eta(2)*BilinearCoeff(:,1)+xi(2)*BilinearCoeff(:,2)+eta(2)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
      t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
     ! IF((t(2).LT.-epsilontol).OR.(t(2).GT.epsilonOne))THEN
!      SELECT CASE(iter)
!      CASE(162)
!        print*,'t',t(2)
!        read*
!      END SELECT
      IF((t(2).LT.-epsilontol).AND.(t(2).GT.lengthPartTrajectory))THEN
        t(2)=-2.0
      ELSE
        t(2)=t(2)!/lengthPartTrajectory
        nInter=nInter+1
      END IF
!      IF((t(2).LT.epsilontol).AND.(t(2).GT.epsilonOne))THEN
!        t(2)=-1
!      END IF
      !IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
      !!  WRITE(*,*) ' Second Intersection'
      !!  WRITE(*,*) ' t ', t(2)
      !!  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      !END IF 
    END IF
  END IF
  ! if no intersection, return
!  print*,nInter
!  print*,'xi,eta,t',xi(2),eta(2),t(2)
!  print*,'t(1),t(2)',t(1),t(2)
  IF(nInter.EQ.0) RETURN
  IF(SideID.LE.nBCSides)THEN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  ELSE ! no BC Side
    ! if two intersections, return, particle re-enters element
    IF(nInter.EQ.2) RETURN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  END IF ! SideID.LT.nCBSides
 !print*,'xi,eta,t',xitild,etatild,t(2)
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersectionSuperSampled2

SUBROUTINE ComputeBiLinearIntersectionSuperSampled(xNodes,PartTrajectory,alpha,xitild,etatild,iPart,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2), q1(3)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0

! compute initial vectors
BiLinearCoeff(:,1) = xNodes(:,1)-xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
BiLinearCoeff(:,2) =-xNodes(:,1)+xNodes(:,2)+xNodes(:,3)-xNodes(:,4)
BiLinearCoeff(:,3) =-xNodes(:,1)-xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
BiLinearCoeff(:,4) = xNodes(:,1)+xNodes(:,2)+xNodes(:,3)+xNodes(:,4)
BiLinearCoeff= 0.25*BiLinearCoeff

! compute product with particle trajectory
a1(1)= BilinearCoeff(1,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(1)
a1(4)=(BilinearCoeff(1,4)-LastPartPos(iPart,1))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1)*PartTrajectory(3) - BilinearCoeff(3,1)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2)*PartTrajectory(3) - BilinearCoeff(3,2)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3)*PartTrajectory(3) - BilinearCoeff(3,3)*PartTrajectory(2)
a2(4)=(BilinearCoeff(2,4)-LastPartPos(iPart,2))*PartTrajectory(3) &
     -(BilinearCoeff(3,4)-LastPartPos(iPart,3))*PartTrajectory(2)

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)
!print*,'A,B,C', A,B,C
CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
!read*
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      !q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      IF((t(1).GE.-epsilontol).AND.(t(1).LE.epsilonOne))THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. epsilonOne
  ELSE ! eta not in reange
    RETURN 
  END IF ! eta .lt. epsilonOne
ELSE 
  nInter=0
  t=-2.
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
     ! q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !  WRITE(*,*) ' t ', t(2)
      !  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      IF((t(1).LT.-epsilontol).OR.(t(1).GT.epsilonOne))THEN
        t(1)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(1).LT.epsilontol).AND.(t(1).GT.epsilonOne))THEN
!        t(1)=-1
!      END IF
    END IF
  END IF
  IF(ABS(eta(2)).LT.epsilonOne)THEN
    xi(2)=eta(2)*a2(1)-eta(2)*a1(1)+a2(2)-a1(2)
    xi(2)=1.0/xi(2)
    xi(2)=(eta(2)*a1(3)-eta(2)*a2(3)+a1(4)-a2(4))*xi(2)
    IF(ABS(xi(2)).LT.epsilonOne)THEN
      ! q1=xi(2)*eta(2)*BilinearCoeff(:,1)+xi(2)*BilinearCoeff(:,2)+eta(2)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
      IF((t(2).LT.-epsilontol).OR.(t(2).GT.epsilonOne))THEN
        t(2)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(2).LT.epsilontol).AND.(t(2).GT.epsilonOne))THEN
!        t(2)=-1
!      END IF
      !IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
      !!  WRITE(*,*) ' Second Intersection'
      !!  WRITE(*,*) ' t ', t(2)
      !!  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      !END IF 
    END IF
  END IF
  ! if no intersection, return
!  print*,nInter
!  print*,'xi,eta,t',xi(2),eta(2),t(2)
!  print*,'t(1),t(2)',t(1),t(2)
  IF(nInter.EQ.0) RETURN
  IF(SideID.LE.nBCSides)THEN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  ELSE ! no BC Side
    ! if two intersections, return, particle re-enters element
    IF(nInter.EQ.2) RETURN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  END IF ! SideID.LT.nCBSides
 !print*,'xi,eta,t',xitild,etatild,t(2)
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersectionSuperSampled

SUBROUTINE QuatricSolver(A,B,C,nRoot,r1,r2)
!================================================================================================================================
! subroutine to compute the modified a,b,c equation, parameter already mapped in final version
!================================================================================================================================
USE MOD_Equation_Vars,       ONLY:epsMach
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: A,B,C
!--------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(OUT)     :: nRoot
REAL,INTENT(OUT)        :: R1,R2
!--------------------------------------------------------------------------------------------------------------------------------
! local variables
REAL                    :: radicant
!================================================================================================================================

radicant = B*B-4.0*A*C
IF(ABS(a).LT.epsMach)THEN
  IF(ABS(b).LT.epsMach)THEN
    nRoot=0
    R1=0.
    R2=0.
  ELSE
    nRoot=1
    R1=-c/b
    R2=0.
  END IF
ELSE
  IF(radicant.LT.-epsMach) THEN
    nRoot = 0
    R1=0.
    R2=0.
  ELSE IF (ABS(radicant).LT.epsMach)THEN
    nRoot =1
    R1 = -0.5*B/A
    R2 = 0.
  ELSE
    nRoot=2
    R1 = SQRT(B*B-4.0*A*C)
    R2 = -R1
    R1 = -B+R1
    R1 = 0.5*R1/A
    R2 = -B+R2
    R2 = 0.5*R2/A
  END IF
END IF

END SUBROUTINE QuatricSolver

FUNCTION ComputeSurfaceDistance(xi,eta,PartTrajectory,iPart,SideID)
!================================================================================================================================
! compute the required vector length to intersection
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol,BiLinearCoeff
USE MOD_Particle_Vars,            ONLY:PartState,LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: iPart,SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeSurfaceDistance
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
!================================================================================================================================

IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
  t =xi*eta*BiLinearCoeff(1,1,SideID)+xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) &
             -lastPartPos(iPart,1)
  t = t/ PartTrajectory(1)-epsilontol 
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1,SideID)+xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) &
             -lastPartPos(iPart,2)
  t = t/ PartTrajectory(2)-epsilontol 
ELSE
  t =xi*eta*BilinearCoeff(3,1,SideID)+xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) &
             -lastPartPos(iPart,3)
  t = t/ PartTrajectory(3)-epsilontol 
END IF

ComputeSurfaceDistance=t

END FUNCTION ComputeSurfaceDistance


FUNCTION ComputeSurfaceDistance2(BiLinearCoeff,xi,eta,PartTrajectory,iPart)
!================================================================================================================================
! compute the required vector length to intersection
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol
USE MOD_Particle_Vars,            ONLY:PartState,LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,DIMENSION(3),INTENT(IN)         :: BiLinearCoeff(1:3,4)
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: iPart
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeSurfaceDistance2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
!================================================================================================================================

IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -lastPartPos(iPart,1)
  t = t/ PartTrajectory(1)-epsilontol 
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -lastPartPos(iPart,2)
  t = t/ PartTrajectory(2)-epsilontol 
ELSE
  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -lastPartPos(iPart,3)
  t = t/ PartTrajectory(3)-epsilontol 
END IF

ComputeSurfaceDistance2=t

END FUNCTION ComputeSurfaceDistance2


END MODULE MOD_Particle_Tracking
