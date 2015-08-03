#include "boltzplatz.h"

MODULE MOD_Particle_InterSection
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ComputeBezierIntersection
  MODULE PROCEDURE ComputeBezierIntersection
END INTERFACE


INTERFACE PartInElemCheck
  MODULE PROCEDURE PartInElemCheck
END INTERFACE

INTERFACE ComputeBilinearIntersectionSuperSampled2
  MODULE PROCEDURE ComputeBilinearIntersectionSuperSampled2
END INTERFACE

INTERFACE ComputePlanarInterSectionBezier
  MODULE PROCEDURE ComputePlanarInterSectionBezier
END INTERFACE

INTERFACE ComputePlanarInterSectionBezierRobust
  MODULE PROCEDURE ComputePlanarInterSectionBezierRobust
END INTERFACE

INTERFACE ComputeBilinearIntersectionRobust
  MODULE PROCEDURE ComputeBilinearIntersectionRobust
END INTERFACE


PUBLIC::ComputeBezierIntersection
PUBLIC::ComputeBilinearIntersectionSuperSampled2
PUBLIC::ComputeBilinearIntersectionRobust
PUBLIC::ComputePlanarInterSectionBezier
PUBLIC::ComputePlanarInterSectionBezierRobust
PUBLIC::PartInElemCheck
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


SUBROUTINE PartInElemCheck(PartID,ElemID,Check)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,OneMepsilon,epsilonOne,BezierControlPoints3D,SideType
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: ElemID,PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                      :: Check
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,flip,SideID
REAL                                     :: PartTrajectory(1:3)
REAL                                     :: lengthPartTrajectory,tmpPos(3),xNodes(1:3,1:4),tmpLastPartPos(1:3)
LOGICAL                                  :: isHit
REAL                                     :: alpha,eta,xi
!===================================================================================================================================

! backup positions
tmpPos=PartState(PartID,1:3)
tmpLastPartPos(1:3)=LastPartPos(PartID,1:3)
! virtual move to element barycenter
LastPartPos(PartID,1:3)=PartState(PartID,1:3)
PartState(PartID,1:3)=ElemBaryNGeo(:,ElemID)
PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory
isHit=.FALSE.
DO ilocSide=1,6
  !SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
  SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
  flip  = PartElemToSide(E2S_FLIP,ilocSide,ElemID)
  SELECT CASE(SideType(SideID))
  CASE(PLANAR)
    CALL ComputePlanarIntersectionBezier(ishit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta  &
                                                                            ,PartID,flip,SideID)
  CASE(BILINEAR)
    xNodes(1:3,1)=BezierControlPoints3D(1:3,0   ,0   ,SideID)
    xNodes(1:3,2)=BezierControlPoints3D(1:3,NGeo,0   ,SideID)
    xNodes(1:3,3)=BezierControlPoints3D(1:3,NGeo,NGeo,SideID)
    xNodes(1:3,4)=BezierControlPoints3D(1:3,0   ,NGeo,SideID)
    CALL ComputeBiLinearIntersectionSuperSampled2(ishit,xNodes &
                                                        ,PartTrajectory,lengthPartTrajectory,Alpha &
                                                                                      ,xi                       &
                                                                                      ,eta                ,PartID,flip,SideID)
  CASE(CURVED)
    CALL ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha &
                                                                            ,xi                 &
                                                                            ,eta                ,PartID,SideID)
  END SELECT
  IF(alpha.GT.-1.0) THEN
    IF((ABS(xi).GT.1.0).OR.(ABS(eta).GT.1.0)) THEN
      isHit=.FALSE.
    END IF
  END IF
  IF(isHit) EXIT
END DO ! ilocSide
Check=.TRUE.
IF(isHit) Check=.FALSE.
PartState(PartID,1:3)   = tmpPos
LastPartPos(PartID,1:3) = tmpLastPartPos

END SUBROUTINE PartInElemCheck


SUBROUTINE ComputeBezierIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID)
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
USE MOD_TimeDisc_Vars,           ONLY:iter
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
LOGICAL,INTENT(OUT)                      :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,iInter,p,q
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,ALLOCATABLE,DIMENSION(:)         :: locID!,realInterID
LOGICAL                                  :: firstClip
INTEGER                                  :: realnInter
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
alpha=-1.0
Xi   = 2.0
Eta  = 2.0
isHit=.FALSE.



!!!! DEBUGGG fix side ID
!print*,'sideid',sideid
! If side is flat, than check if particle vector is perpenticular to side. if true, then particle moves parallel to or in side
IF(BoundingBoxIsEmpty(SideID))THEN
  IF(ABS(DOT_PRODUCT(PartTrajectory,SideNormVec(1:3,SideID))).LT.epsilontol) RETURN
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

!IF(ipart.EQ.223)THEN
!  print*,'hello '
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

!IF(ipart.EQ.223)THEN
!  print*,'do check intersection'
!END IF




!  this part in a new function or subroutine
locAlpha=-1.0
iClipIter=1
nXiClip=0
nEtaClip=0
nInterSections=0
firstClip=.TRUE.
CALL BezierClip(firstClip,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory&
               ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)


!IF(ipart.EQ.1899)THEN
!  IF(nInterSections.GT.0)THEN
!    print*,'nintersections',nintersections
!    print*,'alpha ',localpha
!  END IF
!   !read*
!END IF
!IF(ipart.EQ.214)THEN
!  print*,'ninter ',nintersections
!  read*
!END IF



SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha=locAlpha(nInterSections)
  xi =locXi (1)
  eta=loceta(1)
  isHit=.TRUE.
CASE DEFAULT
  isHit=.TRUE.
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
    realnInter=1
    !ALLOCATE(RealInterID(1:nInterSections))
    !print*,MOD(nInterSections,2)
    DO iInter=2,nInterSections
      IF(  (locAlpha(1)/locAlpha(iInter).LT.0.998) &
      .AND.(locAlpha(1)/locAlpha(iInter).GT.1.002))THEN
          realNInter=realnInter+1
      END IF
    END DO
    IF(MOD(realNInter,2).EQ.0) THEN
      DEALLOCATE(locID)
      !print*,'fuck here'
      RETURN ! leave and enter a cell multiple times
    ELSE
      !print*,'inter hit'
      alpha=locAlpha(nInterSections)
      xi =locXi (locID(nInterSections))
      eta=loceta(locID(nInterSections))
      DEALLOCATE(locID)
      RETURN
    END IF
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
USE MOD_Particle_Surfaces_Vars,  ONLY:ClipTolerance,ClipMaxIter,ArrayNchooseK,FacNchooseK,ClipMaxInter,SplitLimit,ClipHit
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,mEpsilontol
USE MOD_Particle_Vars,           ONLY:LastPartPos,Time
USE MOD_TimeDisc_Vars,               ONLY:iter
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
REAL,DIMENSION(3)                    :: IntersectionVector
REAL,DIMENSION(2)                    :: LineNormVec
REAL                                 :: PatchDOF2D
REAL                                 :: minmax(1:2,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp(2,0:NGeo,0:NGeo)
REAL                                 :: BezierControlPoints2D_temp2(2,0:NGeo,0:NGeo)
INTEGER                              :: p,q,l,iDeCasteljau
REAL                                 :: Xi,Eta,XiMin,EtaMin,XiMax,EtaMax,XiSplit,EtaSplit,alpha
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

DO iClipIter=iClipIter,ClipMaxIter
  ! a) xi-direction
  IF(DoXiClip)THEN
    CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoCheck)
    IF(.NOT.DoCheck) EXIT
    DO q=0,NGeo 
      DO p=0,NGeo
        BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec)
      END DO
    END DO
    DO l=0,NGeo
      minmax(2,l)=MAXVAL(BezierControlPoints1D(l,:)) 
      minmax(1,l)=MINVAL(BezierControlPoints1D(l,:)) 
    END DO ! l
    ! calc Smin and Smax and check boundaries
    CALL CalcSminSmax(minmax,XiMin,XiMax)

    IF((XiMin.EQ.1.5).OR.(XiMax.EQ.-1.5))RETURN
    nXiClip=nXiClip+1
    ! 1.) CLIPPING xi
    IF((XiMax-XiMin).GT.SplitLimit)THEN ! two possible intersections
      XiSplit=0.5*(XiMax+XiMin)
      ! first split
      XiArray(:,nXiClip)=(/XiSplit,XiMax/)
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

      tmpnClip=iClipIter+1
      tmpnXi   =nXiClip
      tmpnEta  =nEtaClip
      firstClip=.FALSE.
      CALL BezierClip(firstClip,BezierControlPoints2D_temp2,PartTrajectory,lengthPartTrajectory &
                     ,tmpnClip,tmpnXi,tmpnEta,nInterSections,iPart,SideID)

      ! second split
      XiArray(:,nXiClip)=(/XiMin,XiSplit/)
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
      XiArray(:,nXiClip)=(/XiMin,XiMax/)
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

      ! c) check Tolerance
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
      !IF(ZeroDistance.LT.ClipTolerance2) EXIT
      IF(SQRT(ZeroDistance).LT.ClipTolerance) EXIT

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
    IF((EtaMin.EQ.1.5).OR.(EtaMax.EQ.-1.5))RETURN
    nEtaClip=nEtaClip+1
    ! 2.) CLIPPING eta
    IF((EtaMax-EtaMin).GT.SplitLimit)THEN ! two possible intersections
      EtaSplit=0.5*(EtaMax+EtaMin)
      ! first clip
      EtaArray(:,nEtaClip)=(/EtaSplit,EtaMax/)
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
      EtaArray(:,nEtaClip)=(/EtaMin,EtaSplit/)
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
      EtaArray(:,nEtaClip)=(/EtaMin,EtaMax/)
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
       !IF(ZeroDistance.LT.ClipTolerance2) EXIT
       IF(SQRT(ZeroDistance).LT.ClipTolerance) EXIT

       IF(ABS(EtaMax-EtaMin).LT.ClipTolerance) DoEtaClip=.FALSE.

      END IF
    END IF ! check clip size
  END IF ! DoEtaClip 
END DO

IF(iClipIter.GE.ClipMaxIter)THEN
  WRITE(*,*) 'Bezier Clipping not converged!'
END IF

IF(DoCheck)THEN
  ! back transformation of sub-level clipping values to original bezier surface: ximean, etamean
  !   xi-direction
  IF(nXiClip.EQ.0)THEN
    Xi=0.
  ELSE
    Xi=0.5*SUM(XiArray(:,nXiClip))
    DO iClip=nXiClip-1,1,-1
      Xi=XiArray(1,iClip)+0.5*(Xi+1)*(XiArray(2,iClip)-XiArray(1,iClip))
    END DO
  END IF ! nXIClip
  !   eta-direction
  IF(nEtaClip.EQ.0)THEN
    Eta=0.
  ELSE
    Eta=0.5*SUM(EtaArray(:,nEtaClip))
    DO iClip=nEtaClip-1,1,-1
      Eta=EtaArray(1,iClip)+0.5*(Eta+1)*(EtaArray(2,iClip)-EtaArray(1,iClip))
    END DO
  END IF ! nEtaclip
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
  
  ! resulting point is ReducedBezierControlPoints(:,1,1)
  IntersectionVector=ReducedBezierControlPoints(:,0,0)-LastPartPos(iPart,1:3)
  alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)

  ! funny hard coded tolerance :), obtained by numerical experiments
  !IF((alpha/lengthPartTrajectory.LE.1.0000464802767983).AND.(alpha.GT.Mepsilontol))THEN
  IF((alpha/lengthPartTrajectory.LE.ClipHit).AND.(alpha.GT.Mepsilontol))THEN
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


SUBROUTINE ComputePlanarIntersectionBezier(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,flip,SideID)!,doTest)
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Globals_Vars,            ONLY:epsMach
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Vars,           ONLY:LastPartPos,PartState
USE MOD_Particle_Mesh_Vars,      ONLY:GEO,PartElemToSide,SidePeriodicDisplacement,SidePeriodicType
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne,SideDistance,ClipHit
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
!USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElementNoMap
!USE MOD_Equations_Vars,          ONLY:epsMach
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID,locSideID
INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2
REAL                              :: Inter1(1:3)
REAL                              :: locBezierControlPoints3D(1:3,0:1,0:1)
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,locSideDistance,SideBasePoint(1:3)
!INTEGER                           :: flip
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.
isHit=.FALSE.

!flip   = PartElemToSide(E2S_FLIP,LocSideID,ElemID)
! new
!IF(flip.EQ.0)THEN
!  NormVec  =SideNormVec(1:3,SideID)
!  locDistance=SideDistance(SideID)
!ELSE
!  NormVec  =-SideNormVec(1:3,SideID)
!  locDistance=-SideDistance(SideID)
!END IF
!coeffA=DOT_PRODUCT(NormVec,PartTrajectory)
!IF(coeffA.LE.0.0) RETURN

! old
coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!! corresponding to particle starting in plane
!! interaction should be computed in last step
!IF(iPart.EQ.257) IPWRITE(UNIT_stdOut,*) coeffA
!IF(ABS(coeffA).LT.+epsilontol)  RETURN
IF(ABS(coeffA).EQ.0.)  RETURN

! extension for periodic sides
IF(SidePeriodicType(SideID).EQ.0)THEN
  locSideDistance=SideDistance(SideID)-DOT_PRODUCT(LastPartPos(iPart,1:3),SideNormVec(:,SideID))
  !locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(iPart,1:3),NormVec)
  alpha=locSideDistance/coeffA
  locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
  locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
  locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
  locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
ELSE
  SideBasePoint=0.25*(locBezierControlPoints3D(:,0,0) &
                     +locBezierControlPoints3D(:,0,1) & 
                     +locBezierControlPoints3D(:,1,0) &
                     +locBezierControlPoints3D(:,1,1) )
  !flip   = PartElemToSide(E2S_FLIP,LocSideID,ElemID)
  ! nothing to do for master side, Side of elem to which owns the BezierPoints
  IF(flip.EQ.0)THEN
    locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
    locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
    locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
    locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
  ELSE
    locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)      -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)   -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)   -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)-SidePeriodicDisplacement(:,SidePeriodicType(SideID))
    ! caution, displacement is for master side computed, therefore, slave side requires negative displacement
    SideBasePoint=SideBasePoint-SidePeriodicDisplacement(:,SidePeriodicType(SideID))
  END IF
  !locSideDistance=DOT_PRODUCT(SideBasePoint-LastPartPos(iPart,1:3),NormVec)
  locSideDistance=DOT_PRODUCT(SideBasePoint-LastPartPos(iPart,1:3),SideNormVec(:,SideID))
  alpha=locSideDistance/coeffA
END IF ! SidePeriodicType


IF((alpha.GT.lengthPartTrajectory) .OR.(alpha.LT.-100*epsMach*coeffA))THEN
  alpha=-1.0
  RETURN
END IF



P1=(-locBezierControlPoints3D(:,0,0)+locBezierControlPoints3D(:,1,0)   &
    -locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) )

P2=(-locBezierControlPoints3D(:,0,0)-locBezierControlPoints3D(:,1,0)   &
    +locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) )

P0=(locBezierControlPoints3D(:,0,0)+locBezierControlPoints3D(:,1,0)    &
   +locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) ) 

P1=0.25*P1
P2=0.25*P2
Inter1=LastPartPos(iPart,1:3)+alpha*PartTrajectory
P0=0.25*P0-Inter1


A1=P1(1)+P1(3)
B1=P2(1)+P2(3)
C1=P0(1)+P0(3)

A2=P1(2)+P1(3)
B2=P2(2)+P2(3)
C2=P0(2)+P0(3)


IF(ABS(B1).GE.ABS(B2))THEN
  xi = A2-B2/B1*A1
  IF(ABS(xi).LT.epsilontol)THEN
    print*,'blabla'
    STOP
  END IF
  xi = (B2/B1*C1-C2)/xi
ELSE
  xi = A1-B1/B2*A2
  IF(ABS(xi).LT.epsilontol)THEN
    print*,'blabla'
    STOP
  END IF
  xi = (B1/B2*C2-C1)/xi
END IF
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

IF(ABS(xi).GT.ClipHit)THEN
!IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

eta=-((A1+A2)*xi+C1+C2)/(B1+B2)
IF(ABS(eta).GT.ClipHit)THEN
  alpha=-1.0
  RETURN
END IF
isHit=.TRUE.

END SUBROUTINE ComputePlanarIntersectionBezier


SUBROUTINE ComputePlanarIntersectionBezierRobust(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,flip,SideID)!,doTest)
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Globals_Vars,            ONLY:epsMach
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Vars,           ONLY:LastPartPos,PartState
USE MOD_Particle_Mesh_Vars,      ONLY:GEO,PartElemToSide,SidePeriodicDisplacement,SidePeriodicType
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,epsilonOne,SideDistance,ClipHit
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
!USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElementNoMap
!USE MOD_Equations_Vars,          ONLY:epsMach
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID,locSideID
INTEGER,INTENT(IN)                :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2
REAL                              :: NormVec(1:3),locDistance,Inter1(1:3)
REAL                              :: locBezierControlPoints3D(1:3,0:1,0:1)
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,locSideDistance,SideBasePoint(1:3)
!INTEGER                           :: flip
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.
isHit=.FALSE.

! new with flip
IF(flip.EQ.0)THEN
  NormVec  =SideNormVec(1:3,SideID)
  locDistance=SideDistance(SideID)
ELSE
  NormVec  =-SideNormVec(1:3,SideID)
  locDistance=-SideDistance(SideID)
END IF
coeffA=DOT_PRODUCT(NormVec,PartTrajectory)

!! corresponding to particle starting in plane
!! interaction should be computed in last step
IF(ABS(coeffA).EQ.0.)  RETURN

! extension for periodic sides
IF(.NOT.DoRefMapping)THEN
 IF(SidePeriodicType(SideID).EQ.0)THEN
   locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(iPart,1:3),NormVec)
   alpha=locSideDistance/coeffA
   locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
   locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
   locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
   locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
 ELSE
   locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
   locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
   locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
   locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
   SideBasePoint=0.25*(locBezierControlPoints3D(:,0,0) &
                      +locBezierControlPoints3D(:,0,1) & 
                      +locBezierControlPoints3D(:,1,0) &
                      +locBezierControlPoints3D(:,1,1) )
   !flip   = PartElemToSide(E2S_FLIP,LocSideID,ElemID)
   ! nothing to do for master side, Side of elem to which owns the BezierPoints
   IF(flip.EQ.0)THEN
     locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
     locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
     locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
     locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
   ELSE
     locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)      -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
     locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)   -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
     locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)   -SidePeriodicDisplacement(:,SidePeriodicType(SideID))
     locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)-SidePeriodicDisplacement(:,SidePeriodicType(SideID))
     ! caution, displacement is for master side computed, therefore, slave side requires negative displacement
     SideBasePoint=SideBasePoint-SidePeriodicDisplacement(:,SidePeriodicType(SideID))
   END IF
   !locSideDistance=DOT_PRODUCT(SideBasePoint-LastPartPos(iPart,1:3),NormVec)
   locSideDistance=DOT_PRODUCT(SideBasePoint-LastPartPos(iPart,1:3),SideNormVec(:,SideID))
   alpha=locSideDistance/coeffA
 END IF ! SidePeriodicType
ELSE
 locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(iPart,1:3),NormVec)
 alpha=locSideDistance/coeffA
 locBezierControlPoints3D(:,0,0)=BezierControlPoints3D(:,0,0,SideID)
 locBezierControlPoints3D(:,0,1)=BezierControlPoints3D(:,0,NGeo,SideID)
 locBezierControlPoints3D(:,1,0)=BezierControlPoints3D(:,NGeo,0,SideID)
 locBezierControlPoints3D(:,1,1)=BezierControlPoints3D(:,NGeo,NGeo,SideID)
END IF

IF(locSideDistance.LT.0)THEN
  ! particle is located outside of element, THEREFORE, an intersection were not detected
  alpha=0.
  isHit=.TRUE.
  RETURN
  ! do I have to compute the xi and eta value? first try: do not re-check new element!
END IF

!IF(iPart.EQ.238.AND.iter.GE.182) IPWRITE(UNIT_stdOut,*) 'a/l',alpha/lengthPartTrajectory
!IF(alpha.GT.lengthPartTrajectory) THEN !.OR.(alpha.LT.-epsilontol))THEN
!IF((alpha.GT.lengthPartTrajectory+epsilontol) .OR.(alpha.LT.-epsilontol))THEN
IF((alpha.GT.lengthPartTrajectory) .OR.(alpha.LT.-epsilontol))THEN
!IF((alpha.GT.lengthPartTrajectory) .OR.(alpha.LT.-100*epsMach*coeffA))THEN
  alpha=-1.0
  RETURN
END IF

P1=(-locBezierControlPoints3D(:,0,0)+locBezierControlPoints3D(:,1,0)   &
    -locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) )

P2=(-locBezierControlPoints3D(:,0,0)-locBezierControlPoints3D(:,1,0)   &
    +locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) )

P0=(locBezierControlPoints3D(:,0,0)+locBezierControlPoints3D(:,1,0)    &
   +locBezierControlPoints3D(:,0,1)+locBezierControlPoints3D(:,1,1) ) 

P1=0.25*P1
P2=0.25*P2
Inter1=LastPartPos(iPart,1:3)+alpha*PartTrajectory
P0=0.25*P0-Inter1

A1=P1(1)+P1(3)
B1=P2(1)+P2(3)
C1=P0(1)+P0(3)

A2=P1(2)+P1(3)
B2=P2(2)+P2(3)
C2=P0(2)+P0(3)

!IF((ABS(P2(1)).GE.ABS(P2(2))).AND.(ABS(P2(1))).GE.ABS(P2(3)))THEN
!  xi1= P1(2)-P2(2)/P2(1)*P1(1)
!  xi2= P1(3)-P2(3)/P2(1)*P1(1)
!  IF(ABS(xi1).GT.ABS(xi2))THEN
!    xi=-(P0(2)-P2(2)/P2(1)*P0(1))/xi1
!  ELSE
!    xi=-(P0(3)-P2(3)/P2(1)*P0(1))/xi2
!  END IF
!ELSE IF(ABS(P2(2)).GE.ABS(P2(3)))THEN
!  xi1= P1(1)-P2(1)/P2(2)*P1(2)
!  xi2= P1(3)-P2(3)/P2(2)*P1(2)
!  IF(ABS(xi1).GT.ABS(xi2))THEN
!    xi=-(P0(1)-P2(1)/P2(2)*P0(2))/xi1
!  ELSE
!    xi=-(P0(3)-P2(3)/P2(2)*P0(2))/xi2
!  END IF
!ELSE
!  xi1= P1(1)-P2(1)/P2(3)*P1(3)
!  xi2= P1(2)-P2(2)/P2(3)*P1(3)
!  IF(ABS(xi1).GT.ABS(xi2))THEN
!    xi=-(P0(1)-P2(1)/P2(3)*P0(3))/xi1
!  ELSE
!    xi=-(P0(2)-P2(2)/P2(3)*P0(3))/xi2
!  END IF
!END IF


IF(ABS(B1).GE.ABS(B2))THEN
  xi = A2-B2/B1*A1
  IF(ABS(xi).LT.epsilontol)THEN
    print*,'blabla'
    STOP
  END IF
  xi = (B2/B1*C1-C2)/xi
ELSE
  xi = A1-B1/B2*A2
  IF(ABS(xi).LT.epsilontol)THEN
    print*,'blabla'
    STOP
  END IF
  xi = (B1/B2*C2-C1)/xi
END IF
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

IF(ABS(xi).GT.ClipHit)THEN
!IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

eta=-((A1+A2)*xi+C1+C2)/(B1+B2)
IF(ABS(eta).GT.ClipHit)THEN
  alpha=-1.0
  RETURN
END IF
isHit=.TRUE.

END SUBROUTINE ComputePlanarIntersectionBezierRobust


SUBROUTINE ComputeBiLinearIntersectionSuperSampled2(isHit,xNodes,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                                   ,iPart,flip,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,epsilonOne,hitepsbi,cliphit
USE MOD_Particle_Vars,ONLY:PartState
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart,SideID,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0
isHit=.FALSE.

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
!IF((iPart.EQ.238).AND.(iter.GT.78))THEN
!  IPWRITE(UNIT_stdOut,*) 'a,b,c',a,b,c
!END IF

CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))

!IF((iPart.EQ.238).AND.(iter.GE.182))  IPWRITE(UNIT_stdOut,*) 'nroots',nroot
IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN

!  IF((iPart.EQ.238).AND.(iter.GE.182))THEN
!    IPWRITE(UNIT_stdOut,*) 'radicant', B*B-4.0*A*C
!    IPWRITE(UNIT_stdOut,*) 'partpos',PartState(iPart,1:3)
!    IPWRITE(UNIT_stdOut,*) 'lastpartpos',LastPartPos(iPart,1:3)
!    IPWRITE(UNIT_stdOut,*) 'trajectory',PartTrajectory
!    IPWRITE(UNIT_stdOut,*) 'length',lengthPartTrajectory
!    IPWRITE(UNIT_stdOut,*) 'eta',eta(1)
!    xi(1)=ComputeXi(a1,a2,eta(1))
!    IPWRITE(UNIT_stdOut,*) 'xi',xi(1)
!    IPWRITE(UNIT_stdOut,*) 'intetrsect', xi(1)*eta(1)*BiLinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)&
!                               +eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)
!    t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
!    IPWRITE(UNIT_stdOut,*) 'a/l',t(1)/lengthPartTrajectory
!  END IF

  !IF(ABS(eta(1)).LT.epsilonOne)THEN
  IF(ABS(eta(1)).LT.ClipHit)THEN
  !IF(ABS(eta(1)).LT.hitepsbi)THEN
    ! as paper ramsay
    xi(1)=ComputeXi(a1,a2,eta(1))

   ! IF(ABS(xi(1)).LT.hitepsbi)THEN
    !IF(ABS(xi(1)).LT.epsilonOne)THEN
    IF(ABS(xi(1)).LT.ClipHit)THEN
      !q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      !IF((t(1).GE.-epsilontol).AND.(t(1).LE.epsilonOne))THEN
      IF((t(1).GE.-epsilontol).AND.(t(1).LE.lengthPartTrajectory))THEN
        alpha=t(1)!/LengthPartTrajectory
        xitild=xi(1)
        etatild=eta(1)
        isHit=.TRUE.
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
  t=-1.
!  IF((iPart.EQ.238).AND.(iter.GE.182))THEN
!    IPWRITE(UNIT_stdOut,*) 'radicant', B*B-4.0*A*C
!    IPWRITE(UNIT_stdOut,*) 'partpos',PartState(iPart,1:3)
!    IPWRITE(UNIT_stdOut,*) 'eta',eta(1)
!    xi(1)=ComputeXi(a1,a2,eta(1))
!    IPWRITE(UNIT_stdOut,*) 'xi',xi(1)
!    t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
!    IPWRITE(UNIT_stdOut,*) 'a/l',t(1)/lengthPartTrajectory
!  END IF

  !IF(ABS(eta(1)).LT.hitepsbi)THEN
  !IF(ABS(eta(1)).LT.epsilonOne)THEN
  IF(ABS(eta(1)).LT.ClipHit)THEN
    ! as paper ramsay
    xi(1)=ComputeXi(a1,a2,eta(1))

    !IF(ABS(xi(1)).LT.hitepsbi)THEN
    !IF(ABS(xi(1)).LT.epsilonOne)THEN
    IF(ABS(xi(1)).LT.Cliphit)THEN
      t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
      IF((t(1).LT.-epsilontol).OR.(t(1).GT.lengthPartTrajectory))THEN
        t(1)=-1.0
      ELSE
        nInter=nInter+1
        t(1)=t(1)!/lengthPartTrajectory
      END IF
    END IF
  END IF ! eta(1)

!  IF((iPart.EQ.238).AND.(iter.GE.182))THEN
!    IPWRITE(UNIT_stdOut,*) 'eta',eta(2)
!    xi(2)=ComputeXi(a1,a2,eta(2))
!    IPWRITE(UNIT_stdOut,*) 'xi',xi(2)
!    t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
!    IPWRITE(UNIT_stdOut,*) 'a/l',t(2)/lengthPartTrajectory
!  END IF

 !IF(ABS(eta(2)).LT.hitepsbi)THEN
 !IF(ABS(eta(2)).LT.epsilonOne)THEN
 IF(ABS(eta(2)).LT.ClipHit)THEN
    ! as paper ramsay
    xi(2)=ComputeXi(a1,a2,eta(2))

    !IF(ABS(xi(2)).LT.hitepsbi)THEN
    !IF(ABS(xi(2)).LT.epsilonOne)THEN
    IF(ABS(xi(2)).LT.ClipHit)THEN
      t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
      IF((t(2).LT.-epsilontol).OR.(t(2).GT.lengthPartTrajectory))THEN
        !print*,'here'
        t(2)=-1.0
      ELSE
        !print*,'why'
        t(2)=t(2)!/lengthPartTrajectory
        nInter=nInter+1
      END IF
    END IF
  END IF
  IF(nInter.EQ.0) RETURN
  isHit=.TRUE.
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


SUBROUTINE ComputeBiLinearIntersectionRobust(isHit,xNodes,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                                   ,iPart,flip,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
! robust version
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,epsilonOne,hitepsbi,cliphit
USE MOD_Particle_Vars,ONLY:PartState
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVecBezier
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
USE MOD_Timedisc_vars,           ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
REAL,INTENT(IN),DIMENSION(1:3,4)  :: xNodes
INTEGER,INTENT(IN)                :: iPart,SideID,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2), normVec(3)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0
isHit=.FALSE.

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
!IF((iPart.EQ.238).AND.(iter.GT.78))THEN
!  IPWRITE(UNIT_stdOut,*) 'a,b,c',a,b,c
!END IF

CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))

!IF((iPart.EQ.238).AND.(iter.GE.182))  IPWRITE(UNIT_stdOut,*) 'nroots',nroot
IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
!  IF((iPart.EQ.238).AND.(iter.GE.182))THEN
!    IPWRITE(UNIT_stdOut,*) 'eta',eta(1)
!    IPWRITE(UNIT_stdOut,*) 'xi',xi(1)
!    IPWRITE(UNIT_stdOut,*) 'intetrsect', xi(1)*eta(1)*BiLinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)&
!                               +eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)
!    IPWRITE(UNIT_stdOut,*) 'a/l',t(1)/lengthPartTrajectory
!  END IF
  xi(1)=ComputeXi(a1,a2,eta(1))
  t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
  IF(flip.EQ.0)THEN
    normVec=CalcBiLinearNormVecBezier(xi(1),eta(1),SideID)
  ELSE
    normVec=-CalcBiLinearNormVecBezier(xi(1),eta(1),SideID)
  END IF
  !xInter=LastpartPos(iPart,1:3)+t(1)*PartTrajactory
  ! check if lastpartpos is inside of the element
  IF(DOT_PRODUCT(t(1)*PartTrajectory,NormVec).LT.0)THEN
    ! or not null??
    alpha=0.
    isHit=.TRUE.
    RETURN
  END IF
  !IF(ABS(eta(1)).LT.epsilonOne)THEN
  IF(ABS(eta(1)).LT.ClipHit)THEN
  !IF(ABS(eta(1)).LT.hitepsbi)THEN
    ! as paper ramsay
    xi(1)=ComputeXi(a1,a2,eta(1))
    IF(ABS(xi(1)).LT.ClipHit)THEN
      IF((t(1).GE.-epsilontol).AND.(t(1).LE.lengthPartTrajectory))THEN
        alpha=t(1)!/LengthPartTrajectory
        xitild=xi(1)
        etatild=eta(1)
        isHit=.TRUE.
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
  t=-1.
!  IF((iPart.EQ.238).AND.(iter.GE.182))THEN
!    IPWRITE(UNIT_stdOut,*) 'radicant', B*B-4.0*A*C
!    IPWRITE(UNIT_stdOut,*) 'partpos',PartState(iPart,1:3)
!    IPWRITE(UNIT_stdOut,*) 'eta',eta(1)
!    xi(1)=ComputeXi(a1,a2,eta(1))
!    IPWRITE(UNIT_stdOut,*) 'xi',xi(1)
!    t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
!    IPWRITE(UNIT_stdOut,*) 'a/l',t(1)/lengthPartTrajectory
!  END IF

  !IF(ABS(eta(1)).LT.hitepsbi)THEN
  !IF(ABS(eta(1)).LT.epsilonOne)THEN
  xi(1)=ComputeXi(a1,a2,eta(1))
  t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)
  IF(flip.EQ.0)THEN
    normVec=CalcBiLinearNormVecBezier(xi(1),eta(1),SideID)
  ELSE
    normVec=-CalcBiLinearNormVecBezier(xi(1),eta(1),SideID)
  END IF
  !xInter=LastpartPos(iPart,1:3)+t(1)*PartTrajactory
  ! check if lastpartpos is inside of the element
  IF(DOT_PRODUCT(t(1)*PartTrajectory,NormVec).LT.0)THEN
    ! or not null??
    alpha=0.
    isHit=.TRUE.
    RETURN
  END IF

  IF(ABS(eta(1)).LT.ClipHit)THEN
    ! as paper ramsay
    IF(ABS(xi(1)).LT.Cliphit)THEN
      IF((t(1).LT.-epsilontol).OR.(t(1).GT.lengthPartTrajectory))THEN
        t(1)=-1.0
      ELSE
        nInter=nInter+1
        t(1)=t(1)!/lengthPartTrajectory
      END IF
    END IF
  END IF ! eta(1)

  xi(2)=ComputeXi(a1,a2,eta(2))
  t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)
  IF(flip.EQ.0)THEN
    normVec=CalcBiLinearNormVecBezier(xi(2),eta(2),SideID)
  ELSE
    normVec=-CalcBiLinearNormVecBezier(xi(2),eta(2),SideID)
  END IF
  !xInter=LastpartPos(iPart,1:3)+t(1)*PartTrajactory
  ! check if lastpartpos is inside of the element
  IF(DOT_PRODUCT(t(2)*PartTrajectory,NormVec).LT.0)THEN
    ! or not null??
    alpha=0.
    isHit=.TRUE.
    RETURN
  END IF

 !IF(ABS(eta(2)).LT.epsilonOne)THEN
 IF(ABS(eta(2)).LT.ClipHit)THEN
    !IF(ABS(xi(2)).LT.epsilonOne)THEN
    IF(ABS(xi(2)).LT.ClipHit)THEN
      IF((t(2).LT.-epsilontol).OR.(t(2).GT.lengthPartTrajectory))THEN
        !print*,'here'
        t(2)=-1.0
      ELSE
        !print*,'why'
        t(2)=t(2)!/lengthPartTrajectory
        nInter=nInter+1
      END IF
    END IF
  END IF
  IF(nInter.EQ.0) RETURN
  isHit=.TRUE.
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

END SUBROUTINE ComputeBiLinearIntersectionRobust


SUBROUTINE QuatricSolver(A,B,C,nRoot,r1,r2)
!================================================================================================================================
! subroutine to compute the modified a,b,c equation, parameter already mapped in final version
!================================================================================================================================
USE MOD_Globals_Vars,       ONLY:epsMach
USE MOD_Globals,            ONLY:AlmostZero
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

!IF(radicant.LT.-epsMach)THEN
IF(ALMOSTZERO(a))THEN
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
  radicant = B*B-4.0*A*C
  IF(radicant.LT.0.)THEN
    nRoot=0
    R1=0.
    R2=0.
  !ELSE IF (radicant.LT.epsMach)THEN
  ELSE IF (radicant.EQ.0.)THEN
    nRoot=1
    R1=-0.5*B/A
    R2=0.
  ELSE 
    nRoot=2
    R1=SQRT(radicant)
    R2=-R1
    write(*,*) 'a',a
    R1=0.5*(-B+R1)/A
    R2=0.5*(-B+R2)/A ! sign above
  END IF
END IF
!IF(ABS(a).LT.epsMach)THEN
!  IF(ABS(b).LT.epsMach)THEN
!    nRoot=0
!    R1=0.
!    R2=0.
!  ELSE
!    nRoot=1
!    R1=-c/b
!    R2=0.
!  END IF
!ELSE
!  IF(radicant.LT.-epsMach) THEN
!    nRoot = 0
!    R1=0.
!    R2=0.
!  ELSE IF (ABS(radicant).LT.epsMach)THEN
!    nRoot =1
!    R1 = -0.5*B/A
!    R2 = 0.
!  ELSE
!    nRoot=2
!    R1 = SQRT(B*B-4.0*A*C)
!    R2 = -R1
!    R1 = -B+R1
!    R1 = 0.5*R1/A
!    R2 = -B+R2
!    R2 = 0.5*R2/A
!  END IF
!END IF

END SUBROUTINE QuatricSolver


FUNCTION ComputeSurfaceDistance2(BiLinearCoeff,xi,eta,PartTrajectory,iPart)
!================================================================================================================================
! compute the required vector length to intersection
! ramsey paper algorithm 3.4
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


IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(3))))THEN
  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -lastPartPos(iPart,1)
  t = t/ PartTrajectory(1)-epsilontol 
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -lastPartPos(iPart,2)
  t = t/ PartTrajectory(2)-epsilontol 
ELSE
  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -lastPartPos(iPart,3)
  t = t/ PartTrajectory(3)-epsilontol 
END IF

!IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!  t =xi*eta*BiLinearCoeff(1,1)+xi*BilinearCoeff(1,2)+eta*BilinearCoeff(1,3)+BilinearCoeff(1,4) -lastPartPos(iPart,1)
!  t = t/ PartTrajectory(1)!-epsilontol 
!ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!  t =xi*eta*BilinearCoeff(2,1)+xi*BilinearCoeff(2,2)+eta*BilinearCoeff(2,3)+BilinearCoeff(2,4) -lastPartPos(iPart,2)
!  t = t/ PartTrajectory(2)!-epsilontol 
!ELSE
!  t =xi*eta*BilinearCoeff(3,1)+xi*BilinearCoeff(3,2)+eta*BilinearCoeff(3,3)+BilinearCoeff(3,4) -lastPartPos(iPart,3)
!  t = t/ PartTrajectory(3)!-epsilontol 
!END IF

ComputeSurfaceDistance2=t

END FUNCTION ComputeSurfaceDistance2


FUNCTION ComputeXi(A1,A2,eta)
!================================================================================================================================
! compute the xi value with algorithm 3.3 of Ramsey paper
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(4),INTENT(IN)         :: A1,A2
REAL,INTENT(IN)                      :: eta
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeXi
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: a,b
!================================================================================================================================

a=eta*A2(1)+A2(2)
b=eta*(A2(1)-A1(1))+A2(2)-A1(2)

IF(ABS(B).GE.ABS(A))THEN
  ComputeXi=(eta*(A1(3)-A2(3))+A1(4)-A2(4))/b
ELSE
  ComputeXi=(-eta*A2(3)-A2(4))/a
END IF

END FUNCTION ComputeXi


END MODULE MOD_Particle_Intersection
