#include "boltzplatz.h"

MODULE MOD_Particle_InterSection
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ComputePlanarCurvedIntersection
  MODULE PROCEDURE ComputePlanarCurvedIntersection
END INTERFACE

INTERFACE ComputePlanarRectInterSection
  MODULE PROCEDURE ComputePlanarRectInterSection
END INTERFACE

INTERFACE ComputeBilinearIntersection
  MODULE PROCEDURE ComputeBilinearIntersection
END INTERFACE

INTERFACE ComputeCurvedIntersection
  MODULE PROCEDURE ComputeCurvedIntersection
END INTERFACE

PUBLIC::ComputePlanarRectInterSection
PUBLIC::ComputePlanarCurvedIntersection
PUBLIC::ComputeBilinearIntersection
PUBLIC::ComputeCurvedIntersection
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


SUBROUTINE ComputePlanarRectIntersection(isHit                       &
                                        ,PartTrajectory              &
                                        ,lengthPartTrajectory        &
                                        ,alpha                       &
                                        ,xi                          &
                                        ,eta                         &
                                        ,iPart                       &
                                        ,flip                        &
                                        ,SideID                      &
                                        ,opt_CriticalParllelInSide   )  
!===================================================================================================================================
! Compute the Intersection with planar surface
! equation of plane: P1*xi + P2*eta+P0
! equation to solve intersection point with plane
! P1*xi+P2*eta+P0-LastPartPos-alpha*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                 ONLY:Cross,abort
USE MOD_Globals_Vars,            ONLY:epsMach
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,epsilontol,OnePlusEps,SideDistance
USE MOD_Particle_Surfaces_Vars,  ONLY:BaseVectors0,BaseVectors1,BaseVectors2
USE MOD_Particle_Surfaces_Vars,  ONLY:BaseVectors0flip,BaseVectors1flip,BaseVectors2flip
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
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
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:3)               :: P0,P1,P2
REAL                              :: NormVec(1:3),locDistance,Inter1(1:3), alphaNorm
REAL                              :: locBezierControlPoints3D(1:3,0:1,0:1)
!REAL,DIMENSION(2:4)               :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                              :: a1,a2,b1,b2,c1,c2
REAL                              :: coeffA,locSideDistance,SideBasePoint(1:3)
REAL                              :: sdet
REAL                              :: epsLoc
LOGICAL                           :: CriticalParallelInSide
!INTEGER                           :: flip
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
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
CriticalParallelInSide=.FALSE.
IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.

! extension for periodic sides
IF(.NOT.DoRefMapping)THEN
  locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(iPart,1:3),NormVec)
ELSE
  locSideDistance=locDistance-DOT_PRODUCT(LastPartPos(iPart,1:3),NormVec)
END IF
  
IF(CriticalParallelInSide)THEN ! particle parallel to side
  IF(ALMOSTZERO(locSideDistance))THEN ! particle on/in side
    IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
    ! move particle eps into interior 
    alpha=-1.
    RETURN
  END IF
  IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
  alpha=-1.
  RETURN
ELSE
  IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
  alpha=locSideDistance/coeffA
END IF

IF(locSideDistance.LT.-100*epsMach)THEN
  ! particle is located outside of element, THEREFORE, an intersection were not detected
  alpha=-1. ! here, alpha was set to zero? why?
  isHit=.FALSE.
  RETURN
  ! do I have to compute the xi and eta value? first try: do not re-check new element!
END IF

alphaNorm=alpha/lengthPartTrajectory

IF((alphaNorm.GT.OnePlusEps) .OR.(alphaNorm.LT.-epsilontol))THEN
  ishit=.FALSE.
  alpha=-1.0
  RETURN
END IF

IF(.NOT.DoRefMapping)THEN
  ! iSide_temp = SideID2PlanarSideID(SideID)
  Inter1=LastPartPos(iPart,1:3)+alpha*PartTrajectory
  P0 =-0.25*BaseVectors0(:,SideID)+Inter1
  P1 = 0.25*BaseVectors1(:,SideID)
  P2 = 0.25*BaseVectors2(:,SideID)
ELSE
  ! iSide_temp = SideID2PlanarSideID(SideID)
  Inter1=LastPartPos(iPart,1:3)+alpha*PartTrajectory
  P0 =-0.25*BaseVectors0(:,SideID)+Inter1
  P1 = 0.25*BaseVectors1(:,SideID)
  P2 = 0.25*BaseVectors2(:,SideID)
END IF

A1=P1(1)*P1(1)+P1(2)*P1(2)+P1(3)*P1(3)
B1=P2(1)*P1(1)+P2(2)*P1(2)+P2(3)*P1(3)
C1=P1(1)*P0(1)+P1(2)*P0(2)+P1(3)*P0(3)

A2=B1
B2=P2(1)*P2(1)+P2(2)*P2(2)+P2(3)*P2(3)
C2=P2(1)*P0(1)+P2(2)*P0(2)+P2(3)*P0(3)

sdet=A1*B2-A2*B1
IF(ABS(sdet).EQ.0)THEN
  CALL abort(& 
  __STAMP__& 
  ,' ABS(sdet).EQ.0!')
END IF
sdet=1.0/sdet
epsLoc=1.0+100.*epsMach


xi=(B2*C1-B1*C2)*sdet

IF(ABS(xi).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF

eta=(-A2*C1+A1*C2)*sdet
IF(ABS(eta).GT.epsLoc)THEN
  alpha=-1.0
  RETURN
END IF
isHit=.TRUE.

END SUBROUTINE ComputePlanarRectIntersection


SUBROUTINE ComputePlanarCurvedIntersection(isHit                       &
                                           ,PartTrajectory              &
                                           ,lengthPartTrajectory        &
                                           ,alpha                       &
                                           ,xi                          &
                                           ,eta                         &
                                           ,iPart                       &
                                           ,flip                        &
                                           ,SideID                      &
                                           ,opt_CriticalParllelInSide)
!===================================================================================================================================
! Compute the intersection with a planar non rectangular face
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,            ONLY:PI
USE MOD_Globals,                 ONLY:Cross,abort,UNIT_stdOut,AlmostZero
USE MOD_Mesh_Vars,               ONLY:NGeo,nBCSides,nSides,BC
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SideType
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,epsilontol,BezierNewtonAngle
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,  ONLY:locXi,locEta,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals,epsilonTol
USE MOD_Utils,                   ONLY:InsertionSort !BubbleSortID
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipMaxIntersec,BezierClipMaxIter
USE MOD_Globals,                 ONLY:myrank
USE MOD_Particle_Surfaces_Vars,  ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)           :: PartTrajectory
REAL,INTENT(IN)                          :: lengthPartTrajectory
INTEGER,INTENT(IN)                       :: iPart,SideID
INTEGER,INTENT(IN)                       :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                         :: alpha,xi,eta
LOGICAL,INTENT(OUT)                      :: isHit
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,iInter,p,q
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
#ifdef CODE_ANALYZE
REAL                                     :: BezierControlPoints2D_tmp(2,0:NGeo,0:NGeo)
#endif /*CODE_ANALYZE*/
INTEGER,ALLOCATABLE,DIMENSION(:)         :: locID!,realInterID
LOGICAL                                  :: CriticalParallelInSide
INTEGER                                  :: realnInter,isInter
REAL                                     :: XiNewton(2)
REAL                                     :: PartFaceAngle,dXi,dEta
REAL                                     :: coeffA
!REAL                                     :: Interval1D,dInterVal1D
!===================================================================================================================================
!PartTrajectory = PartTrajectory
! set alpha to minus 1, asume no intersection
alpha=-1.0
Xi   = 2.0
Eta  = 2.0
isHit=.FALSE.

#ifdef CODE_ANALYZE
rBoundingBoxChecks=rBoundingBoxChecks+1.
#endif /*CODE_ANALYZE*/

CriticalParallelInSide=.FALSE.

IF(DoRefMapping)THEN
  IF(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory).LT.0.)RETURN
ELSE
  coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory) 
  IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.
  IF(flip.EQ.0)THEN
    IF(coeffA.LT.0.)RETURN
  ELSE
    IF(coeffA.GT.0.)RETURN
  END IF
END IF
IF(.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)) RETURN ! the particle does not intersect the 

!IF(DoRefMapping)THEN
!  IF(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory).LT.0.)RETURN
!ELSE
!  ! dependend on master/slave flip
!  coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory) 
!  IF(ALMOSTZERO(coeffA)) CriticalParallelInSide=.TRUE.
!  IF(flip.EQ.0)THEN
!    IF(coeffA.LT.0.)RETURN
!  ELSE
!    IF(coeffA.GT.0.)RETURN
!  END IF
!END IF
!! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
!IF(.NOT.InsideBoundingBox(LastPartPos(iPart,1:3),SideID))THEN ! the old particle position is not inside the bounding box
!  IF(.NOT.InsideBoundingBox(PartState(iPart,1:3),SideID))THEN ! the new particle position is not inside the bounding box
!    IF(.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)) RETURN ! the particle does not intersect the 
!                                                                                              ! bounding box
!  END IF
!END IF

! 2.) Bezier intersection: transformation of bezier patch 3D->2D
!PartTrajectory = PartTrajectoryOrig + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF(ABS(PartTrajectory(3)).LT.0.)THEN
  n1=(/ -PartTrajectory(2)-PartTrajectory(3)  , PartTrajectory(1) ,PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3) , PartTrajectory(3) , -PartTrajectory(1)-PartTrajectory(2) /)
END IF

n1=n1/SQRT(DOT_PRODUCT(n1,n1))
n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
       , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
       , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
n2=n2/SQRT(DOT_PRODUCT(n2,n2))

DO q=0,NGeo
  DO p=0,NGeo
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n1)
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n2)
  END DO
END DO

XiNewton=0.
CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,iPart,SideID)
nInterSections=0
IF(locAlpha(1).GT.-1) nInterSections=1
locXi (1)=XiNewton(1)
locEta(1)=XiNewton(2)

IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
IF(CriticalParallelInSide)THEN
  IF(ALMOSTZERO(locAlpha(1))) THEN
    IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
  END IF
END IF

SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha=locAlpha(1)
  xi =locXi (1)
  eta=loceta(1)
  isHit=.TRUE.
  RETURN
CASE DEFAULT
  CALL abort(&
__STAMP__&
,' The code should never go here')
END SELECT


END SUBROUTINE ComputePlanarCurvedIntersection


SUBROUTINE ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xitild,etatild &
                                                   ,iPart,flip,SideID)
!===================================================================================================================================
! Compute the Intersection with planar surface
! robust version
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides,nSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilontol,OnePlusEps,Beziercliphit,BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3
USE MOD_Particle_Mesh_Vars,          ONLY:PartBCSideList
!USE MOD_Particle_Surfaces_Vars,  ONLY:OnePlusEps,SideIsPlanar,BiLinearCoeff,SideNormVec
#ifdef MPI
USE MOD_Mesh_Vars,               ONLY:BC
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff
REAL                              :: A,B,C,alphaNorm
REAL                              :: xi(2),eta(2),t(2)!, normVec(3)
INTEGER                           :: nInter,nRoot!,BCSideID
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0
isHit=.FALSE.

! compute initial vectors
BiLinearCoeff(:,1) = 0.25*BaseVectors3(:,SideID)
BiLinearCoeff(:,2) = 0.25*BaseVectors1(:,SideID)
BiLinearCoeff(:,3) = 0.25*BaseVectors2(:,SideID)
BiLinearCoeff(:,4) = 0.25*BaseVectors0(:,SideID)

! a1(1)= 0.25 * (BaseVectors3(1,SideID)*PartTrajectory(3) - BaseVectors3(3,SideID)*PartTrajectory(1))
! a1(2)= 0.25 * (BaseVectors1(1,SideID)*PartTrajectory(3) - BaseVectors1(3,SideID)*PartTrajectory(1))
! a1(3)= 0.25 * (BaseVectors2(1,SideID)*PartTrajectory(3) - BaseVectors2(3,SideID)*PartTrajectory(1))
! a1(4)=(0.25*BaseVectors0(1,SideID)-LastPartPos(iPart,1))*PartTrajectory(3) &
!      -(0.25*BaseVectors0(3,SideID)-LastPartPos(iPart,3))*PartTrajectory(1)
!      
! a2(1)= 0.25 * (BaseVectors3(2,SideID)*PartTrajectory(3) - BaseVectors3(3,SideID)*PartTrajectory(2))
! a2(2)= 0.25 * (BaseVectors1(2,SideID)*PartTrajectory(3) - BaseVectors1(3,SideID)*PartTrajectory(2))
! a2(3)= 0.25 * (BaseVectors2(2,SideID)*PartTrajectory(3) - BaseVectors2(3,SideID)*PartTrajectory(2))
! a2(4)=(0.25*BaseVectors0(2,SideID)-LastPartPos(iPart,2))*PartTrajectory(3) &
!      -(0.25*BaseVectors0(3,SideID)-LastPartPos(iPart,3))*PartTrajectory(2)

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

CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  xi(1)=ComputeXi(a1,a2,eta(1))
  t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)

  IF(ABS(eta(1)).LT.BezierClipHit)THEN
    IF(ABS(xi(1)).LT.BezierClipHit)THEN
      alphaNorm=t(1)/lengthPartTrajectory
      IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
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
    END IF ! xi .lt. OnePlusEps
  ELSE ! eta not in reange
    RETURN 
  END IF ! eta .lt. OnePlusEps
ELSE 
  nInter=0
  t=-1.

  !IF(ABS(eta(1)).LT.BezierHitEpsBi)THEN
  !IF(ABS(eta(1)).LT.OnePlusEps)THEN
  xi(1)=ComputeXi(a1,a2,eta(1))
  t(1)=ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PartTrajectory,iPart)

  IF(ABS(eta(1)).LT.BezierClipHit)THEN
    ! as paper ramsay
    IF(ABS(xi(1)).LT.BezierCliphit)THEN
      alphaNorm=t(1)/lengthPartTrajectory
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GE.0.))THEN
      IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
        nInter=nInter+1
        isHit=.TRUE.
        t(1)=t(1)
      ELSE
        t(1)=-1.0
      END IF
    END IF
  END IF ! eta(1)


  xi(2)=ComputeXi(a1,a2,eta(2))
  t(2)=ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PartTrajectory,iPart)

 !IF(ABS(eta(2)).LT.OnePlusEps)THEN
 IF(ABS(eta(2)).LT.BezierClipHit)THEN
    !IF(ABS(xi(2)).LT.OnePlusEps)THEN
    IF(ABS(xi(2)).LT.BezierClipHit)THEN
      alphaNorm=t(2)/lengthPartTrajectory
      !IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GE.0.))THEN
      IF((alphaNorm.LT.OnePlusEps) .AND.(alphaNorm.GT.-epsilontol))THEN
        t(2)=t(2)!/lengthPartTrajectory
        isHit=.TRUE.
        nInter=nInter+2
      ELSE
        t(2)=-1.0
      END IF
    END IF
  END IF

  IF(nInter.EQ.0) RETURN
  isHit=.TRUE.
  IF(ALLOCATED(PartBCSideList))THEN ! correspond to DoRefMapping
    !BCSideID=PartBCSideList(SideID)
    !IF((BCSideID.GE.1).AND.(BCSideID.LE.nTotalBCSides))THEN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
    RETURN
  END IF
  ! no refmapping 
  IF(SideID.LE.nSides)THEN
    IF(SideID.LE.nBCSides)THEN
      ! take closest
      SELECT CASE(nInter)
      CASE(1)
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      CASE(2)
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      CASE(3)
        IF(ABS(t(1)).LT.ABS(t(2)))THEN
          alpha=t(1)
          xitild=xi(1)
          etatild=eta(1)
        ELSE
          alpha=t(2)
          xitild=xi(2)
          etatild=eta(2)
        END IF
      END SELECT
    ELSE
      SELECT CASE(nInter)
      CASE(1)
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      CASE(2)
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      CASE(3) ! double intersection leaves and entries element
        alpha=-1.0
        xitild=0.
        etatild=0.
      END SELECT
    END IF
#ifdef MPI
  ELSE
    ! halo side
    IF(BC(SideID).GT.0)THEN ! BC Sides
      ! take closest
      SELECT CASE(nInter)
      CASE(1)
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      CASE(2)
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      CASE(3)
        IF(ABS(t(1)).LT.ABS(t(2)))THEN
          alpha=t(1)
          xitild=xi(1)
          etatild=eta(1)
        ELSE
          alpha=t(2)
          xitild=xi(2)
          etatild=eta(2)
        END IF
      END SELECT
    ELSE
      SELECT CASE(nInter)
      CASE(1)
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      CASE(2)
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      CASE(3) ! double intersection leaves and entries element
        alpha=-1.0
        xitild=0.
        etatild=0.
      END SELECT
    END IF
#endif /*MPI*/
  END IF
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersection


SUBROUTINE ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,opt_CriticalParllelInSide)
!===================================================================================================================================
! Compute the intersection with a Bezier surface
! particle path = LastPartPos+lengthPartTrajectory*PartTrajectory
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,            ONLY:PI
USE MOD_Globals,                 ONLY:Cross,abort,UNIT_stdOut,AlmostZero
USE MOD_Mesh_Vars,               ONLY:NGeo,nBCSides,nSides,BC
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:SideType
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec,epsilontol,BezierNewtonAngle
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,  ONLY:locXi,locEta,locAlpha
USE MOD_Particle_Surfaces_Vars,  ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals,epsilonTol
USE MOD_Utils,                   ONLY:InsertionSort !BubbleSortID
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipMaxIntersec,BezierClipMaxIter
USE MOD_Globals,                 ONLY:myrank
USE MOD_Particle_Surfaces_Vars,  ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
#endif /*CODE_ANALYZE*/
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
LOGICAL,INTENT(OUT),OPTIONAL             :: opt_CriticalParllelInSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: n1(3),n2(3)
INTEGER                                  :: nInterSections,iInter,p,q
INTEGER                                  :: iClipIter,nXiClip,nEtaClip
REAL                                     :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
REAL                                     :: PartTrajectoryDif(3)
#ifdef CODE_ANALYZE
REAL                                     :: BezierControlPoints2D_tmp(2,0:NGeo,0:NGeo)
#endif /*CODE_ANALYZE*/
INTEGER,ALLOCATABLE,DIMENSION(:)         :: locID!,realInterID
LOGICAL                                  :: firstClip
INTEGER                                  :: realnInter,isInter
REAL                                     :: XiNewton(2)
REAL                                     :: PartFaceAngle,dXi,dEta
LOGICAL                                  :: CriticalParallelInSide
!REAL                                     :: Interval1D,dInterVal1D
!===================================================================================================================================
! set alpha to minus 1, asume no intersection
alpha=-1.0
Xi   = 2.0
Eta  = 2.0
isHit=.FALSE.

#ifdef CODE_ANALYZE
rBoundingBoxChecks=rBoundingBoxChecks+1.
#endif /*CODE_ANALYZE*/

CriticalParallelInSide=.FALSE.
IF(BoundingBoxIsEmpty(SideID))THEN
  IF(DoRefMapping)THEN
    IF(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory).LT.0.)RETURN
  ELSE
    IF(ALMOSTZERO(DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory))) CriticalParallelInSide=.TRUE.
  END IF
  IF(.NOT.FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)) RETURN ! the particle does not intersect the 
                                                                                                ! bounding box
ELSE
  ! 1.) Check if LastPartPos or PartState are within the bounding box. If yes then compute a Bezier intersection problem
  IF(.NOT.InsideBoundingBox(LastPartPos(iPart,1:3),SideID))THEN ! the old particle position is not inside the bounding box
    IF(.NOT.InsideBoundingBox(PartState(iPart,1:3),SideID))THEN ! the new particle position is not inside the bounding box
      IF(.NOT.BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)) RETURN ! the particle does not intersect the 
                                                                                                ! bounding box
    END IF
  END IF
END IF
! 2.) Bezier intersection: transformation of bezier patch 3D->2D
!PartTrajectoryOrig=PartTrajectory
!PartTrajectory = PartTrajectory + epsilontol !move minimal in arb. dir. for preventing collapsing BezierControlPoints2D
IF(ABS(PartTrajectory(3)).LT.0.)THEN
  n1=(/ -PartTrajectory(2)-PartTrajectory(3)  , PartTrajectory(1) ,PartTrajectory(1) /)
ELSE
  n1=(/ PartTrajectory(3) , PartTrajectory(3) , -PartTrajectory(1)-PartTrajectory(2) /)
END IF

! check angle to boundingbox (height normal vector)
PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,2,SideID))))
IF(ALMOSTZERO(PartFaceAngle*180/ACOS(-1.)))THEN
  n1=n1! +epsilontol
END IF

n1=n1/SQRT(DOT_PRODUCT(n1,n1))
n2(:)=(/ PartTrajectory(2)*n1(3)-PartTrajectory(3)*n1(2) &
       , PartTrajectory(3)*n1(1)-PartTrajectory(1)*n1(3) &
       , PartTrajectory(1)*n1(2)-PartTrajectory(2)*n1(1) /)
n2=n2/SQRT(DOT_PRODUCT(n2,n2))
!PartTrajectory = PartTrajectoryOrig !set back for preventing angles > 90 deg (0.5pi+eps)

DO q=0,NGeo
  DO p=0,NGeo
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n1)
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-LastPartPos(iPart,1:3),n2)
  END DO
END DO

! calculate angle between particle path and slab normal plane of face
! angle2=abs(90-RadianToDegree*acos(scalar_product/(VECTOR_LENGTH(t)*VECTOR_LENGTH(n2))))
!IF(BoundingBoxIsEmpty(SideID))THEN
!  PartFaceAngle=BezierNewtonAngle+1.0
!ELSE
!END IF
!IF(.NOT.BezierNewtonAngle)THEN
IF((PartFaceAngle.LT.BezierNewtonAngle))THEN ! 1° = 0.01745rad: critical side at the moment need: 0.57° angle
#ifdef CODE_ANALYZE
rPerformBezierClip=rPerformBezierClip+1.
#endif /*CODE_ANALYZE*/
  !print*,"bezier"
  !  this part in a new function or subroutine
  locAlpha=-1.0
  iClipIter=1
  nXiClip=0
  nEtaClip=0
  nInterSections=0
  firstClip=.TRUE.
  dXi =MAXVAL(BezierControlPoints2D(1,:,:))-MINVAL(BezierControlPoints2D(1,:,:))
  dEta=MAXVAL(BezierControlPoints2D(2,:,:))-MINVAL(BezierControlPoints2D(2,:,:))
  IF(dXi.LT.dEta) firstClip=.FALSE.
  ! CALL recursive Bezier clipping algorithm
  CALL BezierClip(firstClip,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory&
                ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
ELSE!BezierNewtonAngle
#ifdef CODE_ANALYZE
rPerformBezierNewton=rPerformBezierNewton+1.
BezierControlPoints2D_tmp=BezierControlPoints2D
locAlpha=-1.0
iClipIter=1
nXiClip=0
nEtaClip=0
nInterSections=0
firstClip=.TRUE.
CALL BezierClip(firstClip,BezierControlPoints2D_tmp,PartTrajectory,lengthPartTrajectory&
              ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
IF(nInterSections.GT.1)THEN
  CALL abort(&
  __STAMP__&
  ,' More then one intersection! Cannot use Newton!' ,nInterSections)
END IF
#endif /*CODE_ANALYZE*/
  !dInterVal1D =MINVAL(BezierControlPoints2D(1,:,:))
  !InterVal1D  =MAXVAL(BezierControlPoints2D(1,:,:))-dInterVal1D
  !XiNewton(1) =MIN(-1.0+ABS(dInterVal1D)/InterVal1D,1.0)
  !InterVal1D  =MINVAL(BezierControlPoints2D(2,:,:))
  !InterVal1D  =MAXVAL(BezierControlPoints2D(2,:,:))-dInterVal1D
  !XiNewton(2) =MIN(-1.0+ABS(dInterVal1D)/InterVal1D,1.0)
  XiNewton=0.
  !CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D_tmp,PartTrajectory,lengthPartTrajectory,iPart,SideID)
  CALL BezierNewton(locAlpha(1),XiNewton,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,iPart,SideID)
  nInterSections=0
  IF(locAlpha(1).GT.-1) nInterSections=1
#ifdef CODE_ANALYZE
  IF(nInterSections.EQ.1)THEN
    dXi =ABS(locXi(1)-XiNewton(1)) !/(400.*BezierClipTolerance)
    dEta=ABS(locEta(1)-XiNewton(2))!/(400.*BezierClipTolerance)
    dXi =dXi*dXi+dEta*dEta
    dXi =SQRT(dXi)/(400.*BezierClipTolerance)
    IF(dXi.GT.1.0) THEN
      IPWRITE(UNIT_stdout,*) ': Difference between Intersections > Tolerance'
      IPWRITE(UNIT_stdout,*) ': xi-clip,   xi-newton', locXi(1), XiNewton(1)
      IPWRITE(UNIT_stdout,*) ': eta-clip, eta-newton', loceta(1), XiNewton(2)
      !CALL abort(__STAMP__ &
      ! ' Wrong intersection in Xi! Clip/Newton=',nInterSections,dXi)
    END IF
    !IF(dXi.GT.1.0)THEN
    !  IPWRITE(UNIT_stdout,*) ' eta-clip, eta-newton', loceta(1), XiNewton(2)
    !  CALL abort(__STAMP__ &
    !   ' Wrong intersection in Eta! Clip/Newton=',nInterSections, dXi)
    !END IF
  END IF
#endif /*CODE_ANALYZE*/
  locXi (1)=XiNewton(1)
  locEta(1)=XiNewton(2)
END IF

SELECT CASE(nInterSections)
CASE(0)
  RETURN
CASE(1)
  alpha=locAlpha(1)
  xi =locXi (1)
  eta=loceta(1)
  isHit=.TRUE.
  IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
  IF(CriticalParallelInSide)THEN
    IF(ALMOSTZERO(alpha))THEN
      IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
    END IF
  END IF
  RETURN
CASE DEFAULT
  ! more than one intersection
  !ALLOCATE(locID(nInterSections))
  ALLOCATE(locID(nInterSections))
  DO iInter=1,nInterSections
    locID(iInter)=iInter
  END DO ! iInter
  ! sort intersection distance
!  CALL BubbleSortID(locAlpha,locID,nIntersections)
  CALL InsertionSort(locAlpha(1:nIntersections),locID,nIntersections)
  
  IF(DoRefMapping)THEN
    DO iInter=1,nInterSections 
      IF(locAlpha(iInter).GT.-1.0)THEN
        alpha=locAlpha(iInter)
        xi =locXi (locID(iInter))
        eta=loceta(locID(iInter))
        DEALLOCATE(locID)
        isHit=.TRUE.
        RETURN 
      END IF
    END DO ! iInter
  ELSE
    ! no ref mapping
    IF(SideID.LE.nSides)THEN
      IF(SideID.LE.nBCSides)THEN
        ! requires first hit with BC
        ! take closest
        DO iInter=1,nInterSections 
          IF(locAlpha(iInter).GT.-1.0)THEN
            alpha=locAlpha(iInter)
            xi =locXi (locID(iInter))
            eta=loceta(locID(iInter))
            DEALLOCATE(locID)
            isHit=.TRUE.
            IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
            IF(CriticalParallelInSide)THEN
              IF(ALMOSTZERO(alpha))THEN
                IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
              END IF
            END IF
            RETURN 
          END IF
        END DO ! iInter
      ELSE ! inner side
        realnInter=1
        isInter=1
        DO iInter=2,nInterSections
          IF(  (locAlpha(1)/locAlpha(iInter).LT.0.998) &
          .AND.(locAlpha(1)/locAlpha(iInter).GT.1.002))THEN
              realNInter=realnInter+1
              isInter=iInter
          END IF
        END DO
        IF(MOD(realNInter,2).EQ.0) THEN
          DEALLOCATE(locID)
          alpha=-1.0
          isHit=.FALSE.
          IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
          IF(CriticalParallelInSide)THEN
            IF(ALMOSTZERO(alpha))THEN
              IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
            END IF
          END IF
          RETURN ! leave and enter a cell multiple times
        ELSE
          alpha=locAlpha(isInter)
          xi =locXi (locID(isInter))
          eta=loceta(locID(isInter))
          isHit=.TRUE.
          DEALLOCATE(locID)
          IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
          IF(CriticalParallelInSide)THEN
            IF(ALMOSTZERO(alpha))THEN
              IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
            END IF
          END IF
          RETURN
        END IF
      END IF ! inner Side
    END IF ! SideID.LE.nSides
#ifdef MPI
    ! halo side
    IF(BC(SideID).GT.0)THEN ! BC Sides
      ! take closest
      DO iInter=1,nInterSections 
        IF(locAlpha(iInter).GT.-1.0)THEN
          alpha=locAlpha(iInter)
          xi =locXi (locID(iInter))
          eta=loceta(locID(iInter))
          DEALLOCATE(locID)
          isHit=.TRUE.
          IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
          IF(CriticalParallelInSide)THEN
            IF(ALMOSTZERO(alpha))THEN
              IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
            END IF
          END IF
          RETURN 
        END IF
      END DO ! iInter
    ELSE ! no BC Side
      realnInter=1
      isInter=1
      DO iInter=2,nInterSections
        IF(  (locAlpha(1)/locAlpha(iInter).LT.0.998) &
        .AND.(locAlpha(1)/locAlpha(iInter).GT.1.002))THEN
            realNInter=realnInter+1
            isInter=iInter
        END IF
      END DO
      IF(MOD(realNInter,2).EQ.0) THEN
        DEALLOCATE(locID)
        alpha=-1.0
        isHit=.FALSE.
        IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
        IF(CriticalParallelInSide)THEN
          IF(ALMOSTZERO(alpha))THEN
            IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
          END IF
        END IF
        RETURN ! leave and enter a cell multiple times
      ELSE
        alpha=locAlpha(isInter)
        xi =locXi (locID(isInter))
        eta=loceta(locID(isInter))
        isHit=.TRUE.
        IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.FALSE.
        IF(CriticalParallelInSide)THEN
          IF(ALMOSTZERO(alpha))THEN
            IF(PRESENT(opt_CriticalParllelInSide)) opt_CriticalParllelInSide=.TRUE.
          END IF
        END IF
        DEALLOCATE(locID)
        RETURN
      END IF
    END IF ! inner Side
#endif /*MPI*/
  END IF
  SDEALLOCATE(locID)
END SELECT

CALL abort(&
__STAMP__&
,' The code should never go here')

END SUBROUTINE ComputeCurvedIntersection


RECURSIVE SUBROUTINE BezierClip(firstClip,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,iClipIter,nXiClip,nEtaClip&
                               ,nInterSections,iPart,SideID)
!================================================================================================================================
! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
! original article:
!   author = {Nishita, Tomoyuki and Sederberg, Thomas W. and Kakimoto, Masanori},                            
!   title = {Ray Tracing Trimmed Rational Surface Patches},                                                  
!   year = {1990},
! book:
!   author = {Farin, Gerald},
!   title = {Curves and Surfaces for CAGD: A Practical Guide},
!   year = {2002},
!================================================================================================================================
!USE MOD_Globals,                 ONLY:MyRank
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:XiArray,EtaArray,locAlpha,locXi,locEta
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipMaxIter,FacNchooseK,BezierClipMaxIntersec
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,epsilontol,BezierClipHit,BezierSplitLimit
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
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
REAL                                 :: ZeroDistance,BezierClipTolerance2
LOGICAL                              :: DoXiClip,DoEtaClip,DoCheck
INTEGER                              :: iClip
REAL                                 :: alphaNorm
REAL                                 :: PlusXi,MinusXi,PlusEta,MinusEta,tmpXi,tmpEta
INTEGER                              :: tmpnClip,tmpnXi,tmpnEta
REAL                                 :: xiup(0:NGeo),etaup(0:NGeo),xidown(0:NGeo),etadown(0:NGeo)
REAL                                 :: XiBuf(0:NGeo,0:NGeo),EtaBuf(0:NGeo,0:NGeo)
!================================================================================================================================

PatchDOF2D=1.0/REAL((NGeo+1)*(NGeo+1))
BezierClipTolerance2=BezierClipTolerance*BezierClipTolerance

! 3.) Bezier intersection: solution Newton's method or Bezier clipping
! outcome: no intersection, single intersection, multiple intersection with patch
! BEZIER CLIPPING: xi- and eta-direction
IF(FirstClip)THEN
  DoXiClip=.TRUE.
ELSE
  DoXiClip=.FALSE.
END IF
DoEtaClip=.TRUE.
DoCheck=.TRUE.

DO iClipIter=iClipIter,BezierClipMaxIter
  ! a) xi-direction
  !IF(iClipIter.EQ.1)THEN
  !  CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoXiClip)
  !END IF
  IF(DoXiClip)THEN
    !CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoCheck)
    CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoCheck,Mode=1)
    IF(.NOT.DoCheck) EXIT
    DO q=0,NGeo 
      DO p=0,NGeo
        BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec)
      END DO
    END DO
    DO l=0,NGeo
      minmax(1,l)=MINVAL(BezierControlPoints1D(l,:)) 
      minmax(2,l)=MAXVAL(BezierControlPoints1D(l,:)) 
    END DO ! l
    ! DEBUG: additional abort criterion from
    !      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                  
    !      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},               
    !      YEAR = {2005},
    ! IF(MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipHit)THEN ! which tolerance is best used?
    !   DoEtaClip=.FALSE.
    !   BREAK
    ! END IF

    ! calc Smin and Smax and check boundaries
    CALL CalcSminSmax(minmax,XiMin,XiMax)
    IF(nXiClip.EQ.0)THEN
      XiMin=MIN(-1.0,XiMin)
      XiMax=Max( 1.0,XiMax)
    END IF
    IF((XiMin.EQ.1.5).OR.(XiMax.EQ.-1.5))RETURN
    IF(XiMin.GT.XiMax)THEN
      tmpXi=XiMax
      XiMax=XiMin
      XiMin=tmpXi
    END IF

    nXiClip=nXiClip+1
    ! 1.) CLIPPING xi
    IF((XiMax-XiMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
      XiSplit=0.5*(XiMax+XiMin)
      ! first split
      XiArray(:,nXiClip)=(/XiSplit,XiMax/)
      BezierControlPoints2D_temp=0.
      ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
      IF(XiMax.NE.1.0)THEN
        PlusXi=1.0+XiMax
        MinusXi=1.0-XiMax
        ! compute the required stuff || pseudo Horner or precomputation
        xiup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        xidown(0)=1.0
        DO l=1,NGeo
          xiup   (l)     =xiup   (l-1)*PlusXi
          xidown (l)     =xidown (l-1)*MinusXi
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
              !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
              !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
              !                                 *(PlusXi**l)*(MinusXi**(p-l))
              ! not Horner!
              !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
              !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
              !                                 *xiup(l)*XiDown(p-l)

              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                               +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
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
      ! compute the required stuff || pseudo Horner or precomputation
      xiup(0)=1.0
      ! caution, here the indicies are switched from n-j to j  for **down
      xidown(0)=1.0
      DO l=1,NGeo
        xiup   (l)     =xiup   (l-1)*PlusXi
        xidown (l)     =xidown (l-1)*MinusXi
      END DO ! l=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
        END DO ! l=0,p
      END DO ! p=0,NGeo

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
            !DEBUG: optimize this !
            !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
            !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
            !                                 *(PlusXi**l)*(MinusXi**(p-l))
            BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
                                                   +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
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
      ! compute the required stuff || pseudo Horner or precomputation
      xiup(0)=1.0
      ! caution, here the indicies are switched from n-j to j  for **down
      xidown(0)=1.0
      DO l=1,NGeo
        xiup   (l)     =xiup   (l-1)*PlusXi
        xidown (l)     =xidown (l-1)*MinusXi
      END DO ! l=0,NGeo
      DO p=0,NGeo
        DO l=0,p
          XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
        END DO ! l=0,p
      END DO ! p=0,NGeo

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
            !DEBUG: optimize this !
            !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
            !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
            !                                 *(PlusXi**l)*(MinusXi**(p-l))

            BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                             +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
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
        ! compute the required stuff || pseudo Horner or precomputation
        xiup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        xidown(0)=1.0
        DO l=1,NGeo
          xiup   (l)     =xiup   (l-1)*PlusXi
          xidown (l)     =xidown (l-1)*MinusXi
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
              !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
              !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
              !                                 *(PlusXi**l)*(MinusXi**(p-l))
              BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
                                                     +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
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
        ! compute the required stuff || pseudo Horner or precomputation
        xiup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        xidown(0)=1.0
        DO l=1,NGeo
          xiup   (l)     =xiup   (l-1)*PlusXi
          xidown (l)     =xidown (l-1)*MinusXi
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
              !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
              !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
              !                                 *(PlusXi**l)*(MinusXi**(p-l))
              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
                                               +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
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
        ! compute the required stuff || pseudo Horner or precomputation
        xiup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        xidown(0)=1.0
        DO l=1,NGeo
          xiup   (l)     =xiup   (l-1)*PlusXi
          xidown (l)     =xidown (l-1)*MinusXi
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
              !BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)             &
              !                                +BezierControlPoints2D  (:,NGeo-l,q)*FacNchooseK(p,l) &
              !                                *(PlusXi**l)*(MinusXi**(p-l))
              BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)                  &
                                                    +BezierControlPoints2D(:,NGeo-l,q)*XiBuf(p,l)
            END DO
          END DO
        END DO
        BezierControlPoints2D=BezierControlPoints2D_temp
      END IF

      ! c) check Tolerance
       ! check via mean value
      !x=SUM(BezierControlPoints2D(1,:,:))*PatchDOF2D
      !y=SUM(BezierControlPoints2D(2,:,:))*PatchDOF2D
      !IF(SQRT(x*x+y*y).LT.BezierClipTolerance)EXIT
      ! check via distance
      ZeroDistance=0.
      DO q=0,NGeo                                                                                   
        DO p=0,NGeo
          ZeroDistance=ZeroDistance+BezierControlPoints2D(1,p,q)*BezierControlPoints2d(1,p,q) &
                                   +BezierControlPoints2D(2,p,q)*BezierControlPoints2d(2,p,q)
        END DO
      END DO
      ZeroDistance=ZeroDistance*PatchDOF2D ! divide by number of points
      !IF(ZeroDistance.LT.BezierClipTolerance2) EXIT
      IF(SQRT(ZeroDistance).LT.BezierClipTolerance) EXIT

      IF(ABS(XiMax-XiMin).LT.BezierClipTolerance) DoXiClip=.FALSE.


    END IF ! check clip size
  END IF!DoXiClip

  ! b) eta-direction
  !IF(iClipIter.EQ.1)THEN
  !  CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,0,NGeo,DoEtaClip)
  !END IF
  IF(DoEtaClip)THEN
    IF(.NOT.FirstClip)THEN
      DoXiClip=.TRUE.
      FirstClip=.TRUE.
    END IF
    !CALL CalcLineNormVec(BezierControlPoints2D(:,:,:),LineNormVec,0,NGeo,DoCheck)
    CALL CalcLineNormVec2(BezierControlPoints2D(:,:,:),LineNormVec,NGeo,0,DoCheck,Mode=2)
    IF(.NOT.DoCheck) EXIT
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec)
      END DO
    END DO
    DO l=0,NGeo
      minmax(1,l)=MINVAL(BezierControlPoints1D(:,l)) 
      minmax(2,l)=MAXVAL(BezierControlPoints1D(:,l)) 
    END DO ! l
    ! DEBUG: additional abort criterion from
    !      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                 
    !      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},              
    !      YEAR = {2005},
!    IF(MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipTolerance)THEN ! which tolerance is best used?
!      DoEtaClip=.FALSE. ! stop considering this direction
      !print*,"IF(MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipTolerance)THEN"
      !print*,MAXVAL(minmax(1,:))-MINVAL(minmax(2,:))
      !print*,"BezierClipTolerance",BezierClipTolerance
      !read*
      !print*,"1"
!    ELSE
      ! calc Smin and Smax and check boundaries
      CALL CalcSminSmax(minmax,Etamin,Etamax)
      IF(nEtaClip.EQ.0)THEN
        EtaMin=MIN(-1.0,EtaMin)
        EtaMax=Max( 1.0,EtaMax)
      END IF
      IF((EtaMin.EQ.1.5).OR.(EtaMax.EQ.-1.5))RETURN
      IF(EtaMin.GT.EtaMax)THEN
        print*,'swwwaaaaaaaap etta'
      END IF
      nEtaClip=nEtaClip+1
      ! 2.) CLIPPING eta
      IF((EtaMax-EtaMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
        EtaSplit=0.5*(EtaMax+EtaMin)
        ! first clip
        EtaArray(:,nEtaClip)=(/EtaSplit,EtaMax/)
        ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
        IF(EtaMax.NE.1.0)THEN
          BezierControlPoints2D_temp=0.
          PlusEta=1.0+EtaMax
          MinusEta=1.0-EtaMax
          ! compute the required stuff || pseudo Horner or precomputation
          Etaup(0)=1.0
          ! caution, here the indicies are switched from n-j to j  for **down
          Etadown(0)=1.0
          DO l=1,NGeo
            Etaup   (l)     =Etaup   (l-1)*PlusEta
            Etadown (l)     =Etadown (l-1)*MinusEta
          END DO ! l=0,NGeo
          DO p=0,NGeo
            DO l=0,p
              EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
            END DO ! l=0,p
          END DO ! p=0,NGeo

          DO q=0,NGeo
            DO p=0,NGeo
              DO l=0,p
!                BezierControlPoints2D_temp(:,q,p)=&
!                BezierControlPoints2D_temp(:,q,p)+&
!                !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!                BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!                                                 *arrayNchooseK(p,l) &
!                                                 *(1.+Etamax)**l       &
!                                                 *(1.-Etamax)**(p-l)
                !DEBUG: optimize this !
!                BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                                 +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                                 *(PlusEta**l)*(MinusEta**(p-l))

                BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                                 +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
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
        ! compute the required stuff || pseudo Horner or precomputation
        Etaup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        Etadown(0)=1.0
        DO l=1,NGeo
          Etaup   (l)     =Etaup   (l-1)*PlusEta
          Etadown (l)     =Etadown (l-1)*MinusEta
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !  !DEBUG: optimize this !
              !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
              !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
              !                                 *(PlusEta**l)*(MinusEta**(p-l))
              BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                                     +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
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
        ! compute the required stuff || pseudo Horner or precomputation
        Etaup(0)=1.0
        ! caution, here the indicies are switched from n-j to j  for **down
        Etadown(0)=1.0
        DO l=1,NGeo
          Etaup   (l)     =Etaup   (l-1)*PlusEta
          Etadown (l)     =Etadown (l-1)*MinusEta
        END DO ! l=0,NGeo
        DO p=0,NGeo
          DO l=0,p
            EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
          END DO ! l=0,p
        END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
!              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                               *(PlusEta**l)*(MinusEta**(p-l))
!
              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                               +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
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
          ! compute the required stuff || pseudo Horner or precomputation
          Etaup(0)=1.0
          ! caution, here the indicies are switched from n-j to j  for **down
          Etadown(0)=1.0
          DO l=1,NGeo
            Etaup   (l)     =Etaup   (l-1)*PlusEta
            Etadown (l)     =Etadown (l-1)*MinusEta
          END DO ! l=0,NGeo
          DO p=0,NGeo
            DO l=0,p
              EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
            END DO ! l=0,p
          END DO ! p=0,NGeo

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
              !DEBUG: optimize this !
                !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
                !                                 *(PlusEta**l)*(MinusEta**(p-l))
                BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
                                                       +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
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
          ! compute the required stuff || pseudo Horner or precomputation
          Etaup(0)=1.0
          ! caution, here the indicies are switched from n-j to j  for **down
          Etadown(0)=1.0
          DO l=1,NGeo
            Etaup   (l)     =Etaup   (l-1)*PlusEta
            Etadown (l)     =Etadown (l-1)*MinusEta
          END DO ! l=0,NGeo
          DO p=0,NGeo
            DO l=0,p
              EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
            END DO ! l=0,p
          END DO ! p=0,NGeo

          DO q=0,NGeo
            DO p=0,NGeo
              DO l=0,p
!                BezierControlPoints2D_temp(:,q,p)=&
!                BezierControlPoints2D_temp(:,q,p)+&
!                !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!                BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!                                                 *arrayNchooseK(p,l) &
!                                                 *(1.+Etamax)**l       &
!                                                 *(1.-Etamax)**(p-l)
               !DEBUG: optimize this !
!                BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                                 +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!                                                 *(PlusEta**l)*(MinusEta**(p-l))
!
                BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
                                                 +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
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
          ! compute the required stuff || pseudo Horner or precomputation
          Etaup(0)=1.0
          ! caution, here the indicies are switched from n-j to j  for **down
          Etadown(0)=1.0
          DO l=1,NGeo
            Etaup   (l)     =Etaup   (l-1)*PlusEta
            Etadown (l)     =Etadown (l-1)*MinusEta
          END DO ! l=0,NGeo
          DO p=0,NGeo
            DO l=0,p
              EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
            END DO ! l=0,p
          END DO ! p=0,NGeo

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
                !DEBUG: optimize this !
                !BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)        &
                !                                 +BezierControlPoints2D(:,q,NGeo-l)*FacNchooseK(p,l) &
                !                                 *(PlusEta**l)*(MinusEta**(p-l))
                BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)             &
                                                       +BezierControlPoints2D     (:,q,NGeo-l)*EtaBuf(p,l)
              END DO
            END DO
          END DO
          BezierControlPoints2D=BezierControlPoints2D_temp

         ! check via mean value
         !x=SUM(BezierControlPoints2D(1,:,:))*PatchDOF2D
         !y=SUM(BezierControlPoints2D(2,:,:))*PatchDOF2D
         !IF(SQRT(x*x+y*y).LT.BezierClipTolerance)EXIT
         ! check via distance
         ZeroDistance=0.
         DO q=0,NGeo
           DO p=0,NGeo
             ZeroDistance=ZeroDistance+BezierControlPoints2D(1,p,q)*BezierControlPoints2d(1,p,q) &
                                      +BezierControlPoints2D(2,p,q)*BezierControlPoints2d(2,p,q)
           END DO
         END DO
         ZeroDistance=ZeroDistance*PatchDOF2D ! divide by number of points
         !IF(ZeroDistance.LT.BezierClipTolerance2) EXIT
         IF(SQRT(ZeroDistance).LT.BezierClipTolerance) EXIT

         IF(ABS(EtaMax-EtaMin).LT.BezierClipTolerance) DoEtaClip=.FALSE.

        END IF ! EtaMin.NE.-1.0
      END IF ! (EtaMax-EtaMin).GT.BezierSplitLimit 
!    END IF ! MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipTolerance
  END IF ! DoEtaClip 
  !IF(MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipTolerance)THEN ! which tolerance is best used?
    !print*,"2"
  !END IF
END DO ! iClipIter=iClipIter,BezierClipMaxIter
  !IF(MAXVAL(minmax(1,:))-MINVAL(minmax(2,:)).LT.BezierClipTolerance)THEN ! which tolerance is best used?
    !print*,"3"
  !END IF

IF(iClipIter.GE.BezierClipMaxIter)THEN
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

  IF((ABS(eta).GT.BezierClipHit).OR.(ABS(xi).GT.BezierClipHit))THEN
    RETURN
  END IF


  MinusXi =1.0-Xi
  MinusEta=1.0-Eta
  
  ! BEGIN DECASTELJAU ------------------------------------
  ! Wikipedia: "Although the algorithm is slower for most architectures 
  !             when compared with the direct approach, it is more numerically stable."
  ! DEBUG: keep decastejau or implement horner for direct evaluation
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
  ! END DECASTELJAU ------------------------------------
  alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)

  ! funny hard coded tolerance :), obtained by numerical experiments
  !IF((alpha/lengthPartTrajectory.LE.1.0000464802767983).AND.(alpha.GT.MinusEps))THEN
  alphaNorm=alpha/lengthPartTrajectory

  IF((alphaNorm.LE.BezierClipHit).AND.(alphaNorm.GT.-epsilontol))THEN
    ! found additional intersection point
    IF(nInterSections.GE.BezierClipMaxIntersec)THEN
      !nInterSections=nInterSections-1
      !locAlpha=-1
      RETURN
    END IF
    nInterSections=nIntersections+1
    locAlpha(nInterSections)=alpha
    locXi (nInterSections)=tmpXi
    locEta(nInterSections)=tmpEta
  END IF

END IF ! docheck

END SUBROUTINE BezierClip


SUBROUTINE BezierNewton(alpha,Xi,BezierControlPoints2D,PartTrajectory,lengthPartTrajectory,iPart,SideID)
!===================================================================================================================================
! Newton to find root in projected plane for curved tracking
! whole Newton operates in [-1,1]
! output: [-1,1]
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars,               ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipHit,BezierClipMaxIter,BezierControlPoints3D,epsilontol
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars,  ONLY:D_Bezier
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)    :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
REAL,INTENT(IN)    :: PartTrajectory(3),lengthPartTrajectory
INTEGER,INTENT(IN) :: iPart
INTEGER,INTENT(IN) :: SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Xi(2)
REAL,INTENT(OUT)   :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: dXi(2),sdet,dXi2,alphaNorm
REAL               :: dBezierControlPoints2D(2,2,0:NGeo,0:NGeo)
REAL               :: P(2),gradXi(2,2),InterP(3)
REAL               :: BezierClipTolerance2
REAL               :: IntersectionVector(3)
INTEGER            :: nIter
INTEGER            :: dd,nn,l,i,j
INTEGER            :: CycIJ(2),Cyc(2)
!===================================================================================================================================

BezierClipTolerance2=BezierClipTolerance**2
nIter=1
alpha=-1.0
dXi2=1.0
! compute gradient at each control point and use D-matrix
dBezierControlPoints2D=0.
DO nn=1,2
  DO dd=1,2 !iSize
    DO j=0,NGeo
      CycIJ(2)=j
      DO i=0,NGeo
        CycIJ(1)=i
        ! Matrix-vector multiplication
        Cyc=CycIJ
        DO l=0,NGeo
          Cyc(dd)=l  ! d/dxi_dd
          dBezierControlPoints2D(dd,nn,i,j) = dBezierControlPoints2D(dd,nn,i,j) + &
                                              D_Bezier(CycIJ(dd),l)*BezierControlPoints2D(nn,Cyc(1),Cyc(2))
        END DO ! l=0,NGeo
      END DO ! j=0,NGeo
    END DO ! i=0,NGeo
  END DO ! dd=1,2
END DO ! nn=1,2

DO WHILE((dXi2.GT.BezierClipTolerance2).AND.(nIter.LE.BezierClipMaxIter))
  ! compute f(xi) and df(xi)/dxi
  CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,2, BezierControlPoints2D(1:2,    0:NGeo,0:NGeo) &
                                                    ,dBezierControlPoints2D(1:2,1:2,0:NGeo,0:NGeo),Point=P,Gradient=gradXi)
  ! caution with index
  sdet=gradXi(1,1)*gradXi(2,2)-gradXi(1,2)*gradXi(2,1)
  !CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,2, BezierControlPoints2D(1:2,    0:NGeo,0:NGeo) &
  !                                                  , Point=P,Gradient=gradXi)
  ! caution with index

  IF(ABS(sdet).GT.epsilon(0.)) THEN
    sdet=1./sdet
  ELSE !shit
    alpha=-1.0
    Xi=1.5
    EXIT
    !CALL abort(__STAMP__ &
    !   'Bezier-Netwton singular. iter,sdetJac',nIter,sDet)
  END IF

  ! build 2x2 inverse and multiply by vector
  dXi(1) = gradXi(2,2)*P(1)-gradXi(2,1)*P(2)
  dXi(2) =-gradXi(1,2)*P(1)+gradXi(1,1)*P(2)

  dXi    = sdet*dXi
  dXi2   = dXi(1)*dXi(1)+dXi(2)*dXi(2)  

  ! update to new position
  Xi = Xi - dXi

  IF(ANY(ABS(Xi).GT.1.5)) THEN
    ! no intersection of ray and bezier patch
    Xi=1.5
    EXIT
  END IF
  nIter=nIter+1
END DO

IF(nIter.GT.BezierClipMaxIter) CALL abort(&
    __STAMP__&
    ,' Bezier-Newton does not yield root! ')

! check if found Xi,Eta are in parameter range
IF(ABS(xi(1)).GT.BezierClipHit) RETURN
IF(ABS(xi(2)).GT.BezierClipHit) RETURN

! compute 3D intersection
CALL EvaluateBezierPolynomialAndGradient(Xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=InterP) 
IntersectionVector=InterP-LastPartPos(iPart,1:3)

alpha=DOT_PRODUCT(IntersectionVector,PartTrajectory)

alphaNorm=alpha/lengthPartTrajectory

IF((alphaNorm.LE.BezierClipHit).AND.(alphaNorm.GT.-epsilontol)) RETURN
alpha=-1.0
 
END SUBROUTINE BezierNewton


SUBROUTINE calcLineNormVec2(BezierControlPoints2D,LineNormVec,a,b,DoCheck,Mode)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
!================================================================================================================================
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Mesh_Vars,               ONLY:NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
INTEGER,INTENT(IN)                   :: a,b,Mode
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: LineNormVec(1:2)
LOGICAL,INTENT(INOUT)                :: DoCheck
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length,alpha(2),dalpha, doPro,dcorr
REAL,DIMENSION(2)                    :: LXi, Leta
!================================================================================================================================

LXi=(BezierControlPoints2D(:,a,b)-BezierControlPoints2D(:,0,0))+&
    (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,b,a))
Length=SQRT(DOT_PRODUCT(LXi,LXi))

IF(Length.EQ.0)THEN
  DoCheck=.FALSE.
  ! DEBUG: is the complete IF statement dispensable?
  CALL abort(&
  __STAMP__&
  ,'Bezier Clipping -> LineNormVec is Null vector!')
  RETURN
END IF
LXi=LXi/Length

Leta=(BezierControlPoints2D(:,b,a)-BezierControlPoints2D(:,0,0))+&
     (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,a,b))
Length=SQRT(DOT_PRODUCT(Leta,Leta))
IF(Length.EQ.0)THEN
  DoCheck=.FALSE.
  ! DEBUG: is the complete IF statement dispensable?
  CALL abort(&
  __STAMP__&
  ,'Bezier Clipping -> LineNormVec is Null vector!')
  RETURN
END IF
Leta=Leta/Length


! if Lxi and Leta are orientated in opposite directions
doPro=DOT_PRODUCT(Lxi,Leta) ! can be negative

IF(ABS(doPro).GT.0.5)THEN
  ! only required here
  alpha(1)=atan2(Lxi (2),Lxi (1))
  alpha(2)=atan2(Leta(2),Leta(1))
  
  IF(alpha(1).LT.0) alpha(1)=alpha(1)+2.*PI
  IF(alpha(2).LT.0) alpha(2)=alpha(2)+2.*PI
  
  dalpha=alpha(1)-alpha(2)

  IF(dalpha.LT.0)THEN
    ! alpha2.GT.alpha1
    IF(dalpha.LT.-PI)THEN
      IF(dalpha.LT.-1.5*PI)THEN
        ! dcorr=2*PI+dalpha
        dcorr=2*PI+dalpha
      ELSE
        dcorr=dalpha+PI
      END IF
    ELSE
      IF(dalpha.LT.-0.5*PI)THEN
        ! dcorr=2*PI+dalpha
        dcorr=PI+dalpha
      ELSE
        dcorr=dalpha !+0.5*PI
      END IF
    END IF
  ELSE
    ! angle angle2.LT.angle1
    IF(dalpha.GT.PI)THEN
      IF(dalpha.GT.1.5*PI)THEN
        dcorr=-2*PI+dalpha
      ELSE
        dcorr=-PI+dalpha
      END IF
    ELSE
      IF(dalpha.GT.0.5*PI)THEN
        dcorr=-PI+dalpha
      ELSE
        dcorr=dalpha
      END IF
    END IF
  END IF 
  dcorr=dcorr*0.5
  IF(Mode.EQ.1)THEN
    alpha(1)=alpha(1)+dcorr
    LineNormVec(1)=COS(alpha(1))
    LineNormVec(2)=SIN(alpha(1))
  ELSE
    alpha(2)=alpha(2)-dcorr
    LineNormVec(1)=COS(alpha(2))
    LineNormVec(2)=SIN(alpha(2))
  END IF
ELSE
  IF(Mode.EQ.1)THEN
    LineNormVec=Lxi
  ELSE
    LineNormVec=Leta
  END IF
END IF

! DEBUG: fix from (could become zero)
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                   
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},                
!      YEAR = {2005},
END SUBROUTINE calcLineNormVec2



SUBROUTINE calcLineNormVec(BezierControlPoints2D,LineNormVec,a,b,DoCheck)
!================================================================================================================================
! Calculate the normal vector for the line Ls (with which the distance of a point to the line Ls is determined)
!================================================================================================================================
USE MOD_Globals
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
LOGICAL,INTENT(INOUT)                :: DoCheck
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Length
!================================================================================================================================
LineNormVec=(BezierControlPoints2D(:,a,b)-BezierControlPoints2D(:,0,0))+&
            (BezierControlPoints2D(:,NGeo,NGeo)-BezierControlPoints2D(:,b,a))
Length=SQRT(DOT_PRODUCT(LineNormVec,LineNormVec))
! DEBUG: fix from (could become zero)
!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                   
!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},                
!      YEAR = {2005},
IF(Length.EQ.0)THEN
  DoCheck=.FALSE.
  ! DEBUG: is the complete IF statement dispensable?
  CALL abort(&
  __STAMP__&
  ,'Bezier Clipping -> LineNormVec is Null vector!')
  RETURN
END IF
LineNormVec=LineNormVec/Length
END SUBROUTINE calcLineNormVec


SUBROUTINE CalcSminSmax(minmax,Smin,Smax)
!================================================================================================================================
! find upper and lower intersection with convex hull (or no intersection)
! find the largest and smallest roots of the convex hull, pre-sorted values minmax(:,:) are required
!================================================================================================================================
USE MOD_Mesh_Vars,               ONLY:NGeo,Xi_NGeo,DeltaXi_NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipTolerance,BezierClipHit
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
  
  ! adjust Smin and Smax to increase the current range
  ! adapted from: 2005, A. Efremov, Robust and numerically stable bezier clipping method for ray tracing nurbs surfaces
  IF(Smax.GT.-1.5)THEN
    !Smax=MIN(Smax+20.*BezierClipTolerance,1.0)
    Smax=MIN(Smax+100.*BezierClipTolerance,BezierClipHit)
  END IF
  IF(Smin.LT.1.5)THEN
    !Smin=MAX(Smin-20.*BezierClipTolerance,-1.0)
    Smin=MAX(Smin-100.*BezierClipTolerance,-BezierClipHit)
  END IF

END SUBROUTINE calcSminSmax


FUNCTION InsideBoundingBox(ParticlePosition,SideID)
!================================================================================================================================
! check is the particles is inside the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Particle_Surfaces_Vars,  ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
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
! y is perpendicular to xi & eta directions --> check first, smallest intervall
y=DOT_PRODUCT(P,SideSlabNormals(:,2,SideID))
!IF((y.LT.SideSlabIntervals(3,SideID)-epsilontol).OR.(y.GT.SideSlabIntervals(4,SideID)+epsilontol))THEN
IF((y.LT.SideSlabIntervals(3,SideID)-100.*epsMach).OR.(y.GT.SideSlabIntervals(4,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
! than xi
x=DOT_PRODUCT(P,SideSlabNormals(:,1,SideID))
!IF((x.LT.SideSlabIntervals(1,SideID)-epsilontol).OR.(x.GT.SideSlabIntervals(2,SideID)+epsilontol))THEN
IF((x.LT.SideSlabIntervals(1,SideID)-100.*epsMach).OR.(x.GT.SideSlabIntervals(2,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
! than eta
z=DOT_PRODUCT(P,SideSlabNormals(:,3,SideID))
!IF((z.LT.SideSlabIntervals(5,SideID)-epsilontol).OR.(z.GT.SideSlabIntervals(6,SideID)+epsilontol))THEN
IF((z.LT.SideSlabIntervals(5,SideID)-100.*epsMach).OR.(z.GT.SideSlabIntervals(6,SideID)+100.*epsMach))THEN
  InsideBoundingBox=.FALSE.
  RETURN
END IF
InsideBoundingBox=.TRUE.
END FUNCTION InsideBoundingBox


FUNCTION BoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Globals,                  ONLY:abort
USE MOD_Particle_Vars,            ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
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
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
DO i=1,3!x,y,z direction
  !dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(i,:,SideID))
  dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))
  !IF(ABS(dnk).LT.epsilontol)THEN
  IF(ABS(dnk).LT.100.*epsMach)THEN
    dnk=100.*epsMach ! ÜBERPRÜFEN OB SIGN sinn macht
  END IF
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
  END IF
END DO!i
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAXVAL(alpha(1,:))
minvalue=MINVAL(alpha(2,:))

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
!  IF((maxvalue.LT.0).AND.(minvalue.GT.0))THEN
!    CALL abort(&
!  __STAMP__& 
!  ,' BoundingBox check failed!')
!  ELSE
    IF((maxvalue.LT.lengthPartTrajectory+100*epsMach).AND.(maxvalue+100*epsMach.GT.0.))THEN
      !the first intersection is less than lengthPartTrajectory and greater 0
      BoundingBoxIntersection=.TRUE.
    ELSE
      BoundingBoxIntersection=.FALSE.
    END IF
  !END IF
ELSE
  BoundingBoxIntersection=.FALSE.
END IF
END FUNCTION BoundingBoxIntersection


FUNCTION FlatBoundingBoxIntersection(PartTrajectory,lengthPartTrajectory,iPart,SideID)
!================================================================================================================================
! check if the particle trajectory penetrates the bounding box, return TRUE/FALSE
!================================================================================================================================
USE MOD_Globals_Vars
USE MOD_Particle_Vars,            ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierClipTolerance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: lengthPartTrajectory
INTEGER,INTENT(IN)                   :: iPart,SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: FlatBoundingBoxIntersection
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: dnk,alpha(2,3)! alpha(2,2): dummy because we are lazy
REAL                                 :: maxvalue,minvalue,comparevalue
INTEGER                              :: i
!================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) Calculate the projection of the PartTrajectory onto the SideSlabNormals and sort accoring to the sign of T*n
!-----------------------------------------------------------------------------------------------------------------------------------
i=1
dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))

!IF(ABS(dnk).LT.epsilontol)THEN
IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0. ! ÜBERPRÜFEN OB SIGN sinn macht
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
  END IF
END IF
i=3
dnk=DOT_PRODUCT(PartTrajectory,SideSlabNormals(:,i,SideID))
!IF(ABS(dnk).LT.epsilontol)THEN
IF(ABS(dnk).LT.100.*epsMach)THEN
  dnk=0.
ELSE
  IF(dnk.LT.0.)THEN
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
  ELSE
    alpha(1,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i-1,SideID) )/dnk!t_min
    alpha(2,i)=( DOT_PRODUCT(BezierControlPoints3D(:,0,0,SideID)-LastPartPos(iPart,:),SideSlabNormals(:,i,SideID))&
                                                                              +SideSlabIntervals(2*i  ,SideID) )/dnk!t_max
  END IF
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) Get smallest subspace interval
!-----------------------------------------------------------------------------------------------------------------------------------

maxvalue=MAX(alpha(1,1),alpha(1,3))!only x,z directions due to flat surface
minvalue=MIN(alpha(2,1),alpha(2,3))!only x,z directions due to flat surface

IF(maxvalue.LE.minvalue)THEN!smallest interval exists with atleast one point
  IF((maxvalue.LT.0).AND.(minvalue.GT.0))THEN
    FlatBoundingBoxIntersection=.TRUE.
  ELSE
    IF((maxvalue.LT.lengthPartTrajectory+100.*epsMach).AND.(maxvalue+100.*epsMach.GT.0.))THEN
    !the first intersection is less than lengthPartTrajectory and greater 0
      FlatBoundingBoxIntersection=.TRUE.
    ELSE
      FlatBoundingBoxIntersection=.FALSE.
    END IF
  END IF
ELSE
  FlatBoundingBoxIntersection=.FALSE.
END IF

END FUNCTION FlatBoundingBoxIntersection



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
USE MOD_Particle_Vars,            ONLY:LastPartPos
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
