!#include "boltzplatz.h"
!
!MODULE MOD_XiClip
!!===================================================================================================================================
!! Contains global variables provided by the particle surfaces routines
!!===================================================================================================================================
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!PUBLIC
!
!INTERFACE PerformXiClip
!  MODULE PROCEDURE PerformXiClip
!END INTERFACE
!
!PUBLIC::PerformXiClip
!!===================================================================================================================================
!
!CONTAINS
!
!SUBROUTINE PerformXiClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
!                               ,iClipIter,nXiClip,nEtaClip&
!                               ,nInterSections,iPart,SideID)
!!================================================================================================================================
!! Performes the de-Casteljau alogrithm with Clipping to find the intersection between trajectory and surface
!! original article:
!!   author = {Nishita, Tomoyuki and Sederberg, Thomas W. and Kakimoto, Masanori},                            
!!   title = {Ray Tracing Trimmed Rational Surface Patches},                                                  
!!   year = {1990},
!! book:
!!   author = {Farin, Gerald},
!!   title = {Curves and Surfaces for CAGD: A Practical Guide},
!!   year = {2002},
!!================================================================================================================================
!USE MOD_Mesh_Vars,               ONLY:NGeo
!USE MOD_Particle_Surfaces_Vars,  ONLY:XiArray,EtaArray,locAlpha,locXi,locEta
!USE MOD_Particle_Surfaces_Vars,  ONLY:BezierClipLocalTol,BezierClipMaxIter,FacNchooseK,BezierClipMaxIntersec,BezierClipTolerance
!USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,epsilontol,BezierClipHit,BezierSplitLimit
!USE MOD_Particle_Vars,           ONLY:LastPartPos
!USE MOD_Particle_Surfaces,       ONLY:EvaluateBezierPolynomialAndGradient
!USE MOD_Globals,                 ONLY:MyRank,UNIT_stdOut
!#ifdef CODE_ANALYZE
!USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
!USE MOD_Particle_Surfaces,       ONLY:OutputBezierControlPoints
!#endif /*CODE_ANALYZE*/
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!--------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)                      :: lengthPartTrajectory
!REAL,INTENT(INOUT)                   :: BezierControlPoints2D(2,0:NGeo,0:NGeo)
!INTEGER,INTENT(IN)                   :: SideID,iPart
!REAL,INTENT(IN),DIMENSION(1:3)       :: PartTrajectory
!!--------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!! REAL,INTENT(INOUT),DIMENSION(:)      :: locAlpha
!! INTEGER,INTENT(INOUT),DIMENSION(:)   :: locXi,locEta,locID
!INTEGER,INTENT(INOUT)                  :: iClipIter,nXiClip,nEtaClip,nInterSections
!INTEGER(KIND=2),INTENT(INOUT)          :: ClipMode
!REAL,DIMENSION(2,2),INTENT(INOUT)      :: LineNormVec
!!--------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL,DIMENSION(3,0:NGeo,0:NGeo)      :: ReducedBezierControlPoints
!REAL,DIMENSION(0:NGeo,0:NGeo)        :: BezierControlPoints1D
!REAL,DIMENSION(3)                    :: IntersectionVector
!REAL                                 :: PatchDOF2D
!REAL                                 :: minmax(1:2,0:NGeo)
!REAL                                 :: BezierControlPoints2D_temp(2,0:NGeo,0:NGeo)
!REAL                                 :: BezierControlPoints2D_temp2(2,0:NGeo,0:NGeo)
!INTEGER                              :: p,q,l,iDeCasteljau
!REAL                                 :: Xi,Eta,XiMin,EtaMin,XiMax,EtaMax,XiSplit,EtaSplit,alpha
!!REAL                                 :: ZeroDistance,BezierClipTolerance2
!LOGICAL                              :: isNewIntersection
!INTEGER                              :: iClip,iInter
!REAL                                 :: alphaNorm
!REAL                                 :: PlusXi,MinusXi,PlusEta,MinusEta,tmpXi,tmpEta
!INTEGER                              :: tmpnClip,tmpnXi,tmpnEta
!REAL                                 :: xiup(0:NGeo),etaup(0:NGeo),xidown(0:NGeo),etadown(0:NGeo)
!REAL                                 :: XiBuf(0:NGeo,0:NGeo),EtaBuf(0:NGeo,0:NGeo)
!REAL                                 :: dmin,dmax
!INTEGER(KIND=2)                      :: tmpClipMode
!REAL,DIMENSION(2,2)                  :: tmpLineNormVec
!!================================================================================================================================
!
!! Bezier Clip (and Split) in xi
!DO q=0,NGeo
!  DO p=0,NGeo
!    BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec(:,1))
!  END DO
!END DO
!DO l=0,NGeo
!  minmax(1,l)=MINVAL(BezierControlPoints1D(l,:)) 
!  minmax(2,l)=MAXVAL(BezierControlPoints1D(l,:)) 
!END DO ! l
!dmin=MINVAL(minmax(1,:))
!dmax=MAXVAL(minmax(2,:))
!#ifdef CODE_ANALYZE
! IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!   IF(iPart.EQ.PARTOUT)THEN
!     IPWRITE(UNIT_stdout,*) ' dmax-dmin-xi ',ABS(dmax-dmin)
!   END IF
! END IF
!#endif /*CODE_ANALYZE*/
!
!!print*,'dmax,min,...',dmax,dmin,ABS(dmax-dmin),BezierClipTolerance
!! 1D abort criterion from
!!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                  
!!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},               
!!      YEAR = {2005},
!IF(ABS(dmax-dmin).LT.BezierClipLocalTol)THEN ! current patch is converged in xi, then skip xi
!  IF(ClipMode.EQ.3) THEN  ! eta has already converged
!    ClipMode=5            ! no more clipping, we have converged
!    RETURN                ! or stop
!  ELSE
!    ClipMode=4            ! xi is converged, but not eta
!  END IF
!ELSE ! xi not converged, next clip should be in eta
!  ClipMode=2 ! after clipping in xi, clip in eta
!END IF
!
!! we perform the clip  (with XiMin and XiMax) in Xi direction 
!!             or split (if (XiMax-XiMin).GT.BezierSplitLimit) in Xi direction 
!
!! calc Smin and Smax and check boundaries
!CALL CalcSminSmax(minmax,XiMin,XiMax)
!#ifdef CODE_ANALYZE
! IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!   IF(iPart.EQ.PARTOUT)THEN
!     IPWRITE(UNIT_stdout,*) ' XiMin,XiMax ',XiMin,XiMax
!   END IF
! END IF
!#endif /*CODE_ANALYZE*/
!
!IF(nXiClip.EQ.0)THEN
!  XiMin=MIN(-1.0,XiMin)
!  XiMax=Max( 1.0,XiMax)
!END IF
!IF((XiMin.EQ.1.5).OR.(XiMax.EQ.-1.5))RETURN
!
!IF(XiMin.GT.XiMax)THEN
!    print*,'swwwaaaaaaaap xi',XiMin,XiMax
!END IF
!
!nXiClip=nXiClip+1
!! 1.) CLIPPING xi
!IF((XiMax-XiMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
!  XiSplit=0.5*(XiMax+XiMin)
!  ! first split: xi upper
!  ! set mapping array
!  XiArray(:,nXiClip)=(/XiSplit,XiMax/)
!  BezierControlPoints2D_temp=0.
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  IF(XiMax.NE.1.0)THEN
!    PlusXi=1.0+XiMax
!    MinusXi=1.0-XiMax
!    ! compute the required stuff || pseudo Horner or precomputation
!    xiup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    xidown(0)=1.0
!    DO l=1,NGeo
!      xiup   (l)     =xiup   (l-1)*PlusXi
!      xidown (l)     =xidown (l-1)*MinusXi
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!!              BezierControlPoints2D_temp(:,p,q)=&
!!              BezierControlPoints2D_temp(:,p,q)+&
!!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!!              BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
!!                                               *arrayNchooseK(p,l) &
!!                                               *PlusXi**l        &
!!                                               *MinusXi**(p-l)
!!              DEBUG: optimize this !
!!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!!                                               *(PlusXi**l)*(MinusXi**(p-l))
!!               not Horner!
!!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!!                                               *xiup(l)*XiDown(p-l)
!
!          BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                           +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
!        END DO
!      END DO
!    END DO
!
!  ELSE
!    BezierControlPoints2D_temp=BezierControlPoints2D
!  END IF
!
!  ! BOTTOM (mirrored Bernstein Basis evaluation)
!  ! s = (smin+1)/(smax+1) for [-1, +1]
!  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
!  BezierControlPoints2D_temp2=0.
!  PlusXi=(XiSplit+1.0)/(XiMax+1.0)
!  ! MinusXi= 2.0*(1-PlusXi)-1
!  ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!  MinusXi=2.0*PlusXi
!  ! PlusXi=1+ 2.0*(1-s)-1
!  PlusXi=2.0-2.0*PlusXi
!  ! compute the required stuff || pseudo Horner or precomputation
!  xiup(0)=1.0
!  ! caution, here the indicies are switched from n-j to j  for **down
!  xidown(0)=1.0
!  DO l=1,NGeo
!    xiup   (l)     =xiup   (l-1)*PlusXi
!    xidown (l)     =xidown (l-1)*MinusXi
!  END DO ! l=0,NGeo
!  DO p=0,NGeo
!    DO l=0,p
!      XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!    END DO ! l=0,p
!  END DO ! p=0,NGeo
!
!  DO q=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        !BezierControlPoints2D_temp2(:,NGeo-p,q)=&
!        !BezierControlPoints2D_temp2(:,NGeo-p,q)+&
!        !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!        !BezierControlPoints2D_temp  (:,NGeo-l,q)*(1./(2.**p))                        &
!        !                                      *arrayNchooseK(p,l)                  &
!        !                                      *PlusXi**l &
!        !                                      *MinusXi**(p-l)
!        !DEBUG: optimize this !
!        !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
!        !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
!        !                                 *(PlusXi**l)*(MinusXi**(p-l))
!        BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
!                                               +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
!      END DO
!    END DO
!  END DO
!
!  ! Bezier Split
!  tmpnClip =iClipIter
!  ! backup current split level to compute correct intersection, required for back-trafo of intervals
!  tmpnXi   =nXiClip
!  tmpnEta  =nEtaClip
!tmpLineNormVec=LineNormVec
!  tmpClipMode   =ClipMode
!  ! MAYBE set ClipMode for NEXT clip
!  ! HERE, ClipMode currently set above
!#ifdef CODE_ANALYZE
!!      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!!        IF(iPart.EQ.PARTOUT)THEN
!!          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
!!          IPWRITE(UNIT_stdout,*) ' split xi-upper '
!!          CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
!!        END IF
!!      END IF
!#endif /*CODE_ANALYZE*/
!  ! HERE, ClipMode currently set above
!  ! Perform split xi-upper
!  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
!                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
!
!
!  ! second split: xi lower
!  ! restore values to allow for correct back-trafo of intervals (required for intersectionpoint)
!  nXiClip     =tmpnXi
!  nEtaClip    =tmpnEta
!  LineNormVec =tmpLineNormVec
!  ClipMode    =tmpClipMode
!
!  ! set mapping array
!  XiArray(:,nXiClip)=(/XiMin,XiSplit/)
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  BezierControlPoints2D_temp=0.
!  PlusXi=1.0+XiSplit
!  MinusXi=1.0-XiSplit
!  ! compute the required stuff || pseudo Horner or precomputation
!  xiup(0)=1.0
!  ! caution, here the indicies are switched from n-j to j  for **down
!  xidown(0)=1.0
!  DO l=1,NGeo
!    xiup   (l)     =xiup   (l-1)*PlusXi
!    xidown (l)     =xidown (l-1)*MinusXi
!  END DO ! l=0,NGeo
!  DO p=0,NGeo
!    DO l=0,p
!      XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!    END DO ! l=0,p
!  END DO ! p=0,NGeo
!
!  DO q=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        !BezierControlPoints2D_temp(:,p,q)=&
!        !BezierControlPoints2D_temp(:,p,q)+&
!        !!BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!        !BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
!        !                                 *arrayNchooseK(p,l) &
!        !                                 *(1+XiSplit)**l        &
!        !                                 *(1-XiSplit)**(p-l)
!        !DEBUG: optimize this !
!        !BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!        !                                 +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!        !                                 *(PlusXi**l)*(MinusXi**(p-l))
!
!        BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                         +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
!      END DO
!    END DO
!  END DO
!  ! BOTTOM (mirrored Bernstein Basis evaluation)
!  ! s = (smin+1)/(smax+1) for [-1, +1]
!  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
!  IF(XiMin.NE.-1.0)THEN
!    BezierControlPoints2D_temp2=0.
!    PlusXi=(XiMin+1.0)/(XiSplit+1.0)
!    ! MinusXi= 2.0*(1-PlusXi)-1
!    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!    MinusXi=2.0*PlusXi
!    ! PlusXi=1+ 2.0*(1-s)-1
!    PlusXi=2.0-2.0*PlusXi
!    ! compute the required stuff || pseudo Horner or precomputation
!    xiup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    xidown(0)=1.0
!    DO l=1,NGeo
!      xiup   (l)     =xiup   (l-1)*PlusXi
!      xidown (l)     =xidown (l-1)*MinusXi
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!    !      BezierControlPoints2D_temp2(:,NGeo-p,q)=&
!    !      BezierControlPoints2D_temp2(:,NGeo-p,q)+&
!    !      !BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!    !      BezierControlPoints2D_temp (:,NGeo-l,q)*(1./(2.**p))                        &
!    !                                            *arrayNchooseK(p,l)                  &
!    !                                            *(1.+2*((XiMin+1.)/(XiSplit+1.)))**(l-1) &
!    !                                            *(1.-2*((XiMin+1.)/(XiSplit+1.)))**(p-l)
!          !DEBUG: optimize this !
!          !BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)             &
!          !                                 +BezierControlPoints2D_temp  (:,NGeo-l,q)*FacNchooseK(p,l) &
!          !                                 *(PlusXi**l)*(MinusXi**(p-l))
!          BezierControlPoints2D_temp2(:,NGeo-p,q)=BezierControlPoints2D_temp2(:,NGeo-p,q)                  &
!                                                 +BezierControlPoints2D_temp(:,NGeo-l,q)*XiBuf(p,l)
!        END DO
!      END DO
!    END DO
!  ELSE
!    BezierControlPoints2D_temp2=BezierControlPoints2D_temp
!  END IF
!  tmpnClip      =iClipIter
!  tmpnXi        =nXiClip
!  tmpnEta       =nEtaClip
!tmpLineNormVec=LineNormVec
!  tmpClipMode   =ClipMode
!  ! MAYBE set ClipMode for NEXT clip
!  ! HERE, ClipMode currently set above
!#ifdef CODE_ANALYZE
!!      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!!        IF(iPart.EQ.PARTOUT)THEN
!!          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
!!          IPWRITE(UNIT_stdout,*) ' split xi-lower '
!!          CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
!!        END IF
!!      END IF
!#endif /*CODE_ANALYZE*/
!  ! HERE, ClipMode currently set above
!  ! Perform split xi-lower
!  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
!                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
!
!  ! after recursive steps, we are done!
!ELSE  ! no split necessary, only a clip
!
!  ! set mapping array
!  XiArray(:,nXiClip)=(/XiMin,XiMax/)
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  IF(XiMax.NE.1.0)THEN
!    BezierControlPoints2D_temp=0.
!    PlusXi=1.0+XiMax
!    MinusXi=1.0-XiMax
!    ! compute the required stuff || pseudo Horner or precomputation
!    xiup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    xidown(0)=1.0
!    DO l=1,NGeo
!      xiup   (l)     =xiup   (l-1)*PlusXi
!      xidown (l)     =xidown (l-1)*MinusXi
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!!              BezierControlPoints2D_temp(:,p,q)=&
!!              BezierControlPoints2D_temp(:,p,q)+&
!!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!!              BezierControlPoints2D     (:,l,q)*(1./(2.**p))       &
!!                                               *arrayNchooseK(p,l) &
!!                                               *(1+XiMax)**l        &
!!                                               *(1-XiMax)**(p-l)
!!             !DEBUG: optimize this !
!!              BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!!                                               +BezierControlPoints2D     (:,l,q)*FacNchooseK(p,l) &
!!                                               *(PlusXi**l)*(MinusXi**(p-l))
!!
!          BezierControlPoints2D_temp(:,p,q)=BezierControlPoints2D_temp(:,p,q)                  &
!                                           +BezierControlPoints2D     (:,l,q)*XiBuf(p,l)
!        END DO
!      END DO
!    END DO
!    BezierControlPoints2D=BezierControlPoints2D_temp
!  END IF
!  ! BOTTOM (mirrored Bernstein Basis evaluation)
!  ! s = (smin+1)/(smax+1) for [-1, +1]
!  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
!  IF(XiMin.NE.-1.0)THEN
!    BezierControlPoints2D_temp=0.
!    PlusXi=(XiMin+1.0)/(XiMax+1.0)
!    ! MinusXi= 2.0*(1-PlusXi)-1
!    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!    MinusXi=2.0*PlusXi
!    ! PlusXi=1+ 2.0*(1-s)-1
!    PlusXi=2.0-2.0*PlusXi
!    ! compute the required stuff || pseudo Horner or precomputation
!    xiup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    xidown(0)=1.0
!    DO l=1,NGeo
!      xiup   (l)     =xiup   (l-1)*PlusXi
!      xidown (l)     =xidown (l-1)*MinusXi
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        XiBuf(p,l)=XiUp(l)*XiDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!          !BezierControlPoints2D_temp(:,NGeo-p,q)=&
!          !BezierControlPoints2D_temp(:,NGeo-p,q)+&
!          !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!          !BezierControlPoints2D     (:,NGeo-l,q)*(1./(2.**p))                        &
!          !                                      *arrayNchooseK(p,l)                  &
!          !                                      *(1.+2*((XiMin+1.)/(XiMax+1.)))**(l-1) &
!          !                                      *(1.-2*((XiMin+1.)/(XiMax+1.)))**(p-l)
!          !DEBUG: optimize this !
!          !BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)             &
!          !                                +BezierControlPoints2D  (:,NGeo-l,q)*FacNchooseK(p,l) &
!          !                                *(PlusXi**l)*(MinusXi**(p-l))
!          BezierControlPoints2D_temp(:,NGeo-p,q)=BezierControlPoints2D_temp(:,NGeo-p,q)                  &
!                                                +BezierControlPoints2D(:,NGeo-l,q)*XiBuf(p,l)
!        END DO
!      END DO
!    END DO
!    BezierControlPoints2D=BezierControlPoints2D_temp
!  END IF
!  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
!                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
!
!  ! after recursive steps, we are done!
!END IF ! decision between Clip or Split
!
!END SUBROUTINE PerformXiClip
!
!END MODULE MOD_XiClip
