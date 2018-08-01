!#include "boltzplatz.h"
!
!MODULE MOD_EtaClip
!!===================================================================================================================================
!! Contains global variables provided by the particle surfaces routines
!!===================================================================================================================================
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!PUBLIC
!
!INTERFACE PerformEtaClip
!  MODULE PROCEDURE PerformEtaClip
!END INTERFACE
!
!PUBLIC::PerformEtaClip
!!===================================================================================================================================
!
!CONTAINS
!
!SUBROUTINE PerformEtaClip(ClipMode,BezierControlPoints2D,LineNormVec,PartTrajectory,lengthPartTrajectory &
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
!! Bezier Clip (and Split) in eta
!DO q=0,NGeo
!  DO p=0,NGeo
!    BezierControlPoints1D(p,q)=DOT_PRODUCT(BezierControlPoints2D(:,p,q),LineNormVec(:,2))
!  END DO
!END DO
!DO l=0,NGeo
!  minmax(1,l)=MINVAL(BezierControlPoints1D(:,l)) 
!  minmax(2,l)=MAXVAL(BezierControlPoints1D(:,l)) 
!END DO ! l
!dmin=MINVAL(minmax(1,:))
!dmax=MAXVAL(minmax(2,:))
!#ifdef CODE_ANALYZE
! IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!   IF(iPart.EQ.PARTOUT)THEN
!     IPWRITE(UNIT_stdout,*) ' dmax-dmin-eta ',ABS(dmax-dmin)
!   END IF
! END IF
!#endif /*CODE_ANALYZE*/
!
!!print*,'dmax,min,...',dmax,dmin,ABS(dmax-dmin),BezierClipTolerance
!! 1D abort criterion from
!!      AUTHOR = {Efremov, Alexander and Havran, Vlastimil and Seidel, Hans-Peter},                                  
!!      TITLE = {Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces},               
!!      YEAR = {2005},
!IF(ABS(dmax-dmin).LT.BezierClipLocalTol)THEN ! current patch is converged in eta, then skip eta
!  IF(ClipMode.EQ.4) THEN  ! xi has already converged
!    ClipMode=5            ! no more clipping, we have converged
!    RETURN                ! or stop
!  ELSE
!    ClipMode=3            ! eta is converged, but not xi
!  END IF
!ELSE ! eta not converged, next clip should be in xi
!  ClipMode=1 ! after clipping in eta, clip in xi
!END IF
!
!! we perform the clip  (with EtaMin and EtaMax) in Eta direction 
!!             or split (if (EtaMax-EtaMin).GT.BezierSplitLimit) in Eta direction 
!
!! calc Smin and Smax and check boundaries
!CALL CalcSminSmax(minmax,Etamin,Etamax)
!#ifdef CODE_ANALYZE
! IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!   IF(iPart.EQ.PARTOUT)THEN
!     IPWRITE(UNIT_stdout,*) ' EtaMin,EtaMax ',EtaMin,EtaMax
!   END IF
! END IF
!#endif /*CODE_ANALYZE*/
!
!IF(nEtaClip.EQ.0)THEN
!  EtaMin=MIN(-1.0,EtaMin)
!  EtaMax=Max( 1.0,EtaMax)
!END IF
!IF((EtaMin.EQ.1.5).OR.(EtaMax.EQ.-1.5))RETURN
!
!IF(EtaMin.GT.EtaMax)THEN
!    print*,'swwwaaaaaaaap etta',etamin,etamax
!END IF
!
!nEtaClip=nEtaClip+1
!! 2.) CLIPPING eta
!IF((EtaMax-EtaMin).GT.BezierSplitLimit)THEN ! two possible intersections: split the clipped patch at 50%
!  EtaSplit=0.5*(EtaMax+EtaMin)
!  ! first split: eta upper
!  ! set mapping array
!  EtaArray(:,nEtaClip)=(/EtaSplit,EtaMax/)
!  BezierControlPoints2D_temp=0.
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  IF(EtaMax.NE.1.0)THEN
!    PlusEta=1.0+EtaMax
!    MinusEta=1.0-EtaMax
!    ! compute the required stuff || pseudo Horner or precomputation
!    Etaup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    Etadown(0)=1.0
!    DO l=1,NGeo
!      Etaup   (l)     =Etaup   (l-1)*PlusEta
!      Etadown (l)     =Etadown (l-1)*MinusEta
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!!              BezierControlPoints2D_temp(:,q,p)=&
!!              BezierControlPoints2D_temp(:,q,p)+&
!!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!!              BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!!                                               *arrayNchooseK(p,l) &
!!                                               *(1.+Etamax)**l       &
!!                                               *(1.-Etamax)**(p-l)
!!              DEBUG: optimize this !
!!              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!!                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!!                                               *(PlusEta**l)*(MinusEta**(p-l))
!!              
!!              
!!              
!!              
!
!          BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                           +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
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
!  PlusEta=(EtaSplit+1.0)/(EtaMax+1.0)
!  ! MinusXi= 2.0*(1-PlusXi)-1
!  ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!  MinusEta=2.0*PlusEta
!  ! PlusXi=1+ 2.0*(1-s)-1
!  PlusEta=2.0-2.0*PlusEta
!  ! compute the required stuff || pseudo Horner or precomputation
!  Etaup(0)=1.0
!  ! caution, here the indicies are switched from n-j to j  for **down
!  Etadown(0)=1.0
!  DO l=1,NGeo
!    Etaup   (l)     =Etaup   (l-1)*PlusEta
!    Etadown (l)     =Etadown (l-1)*MinusEta
!  END DO ! l=0,NGeo
!  DO p=0,NGeo
!    DO l=0,p
!      EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!    END DO ! l=0,p
!  END DO ! p=0,NGeo
!
!  DO q=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        !BezierControlPoints2D_temp2(:,q,p)=&
!        !BezierControlPoints2D_temp2(:,q,p)+&
!        !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!        !BezierControlPoints2D_temp(:,q,NGeo-l)*(1./(2.**p))                     &
!        !                                      *arrayNchooseK(p,l)               &
!        !                                      *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
!        !                                      *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
!        !  !DEBUG: optimize this !
!        !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
!        !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
!        !                                 *(PlusEta**l)*(MinusEta**(p-l))
!        BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
!                                               +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
!      END DO
!    END DO
!  END DO
!
!  ! Bezier Split
!  tmpnClip =iClipIter
!  ! backup current split level to compute correct intersection, required for back-trafo of intervals
!  tmpnXi   =nXiClip
!  tmpnEta  =nEtaClip
!  tmpLineNormVec=LineNormVec
!  tmpClipMode   =ClipMode
!  ! MAYBE set ClipMode for NEXT clip
!  ! HERE, ClipMode currently set above
!#ifdef CODE_ANALYZE
!!      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
!!        IF(iPart.EQ.PARTOUT)THEN
!!          IPWRITE(UNIT_stdout,*) ' --------------------------------------- '
!!          IPWRITE(UNIT_stdout,*) ' split eta-upper '
!!          CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
!!        END IF
!!      END IF
!#endif /*CODE_ANALYZE*/
!  ! HERE, ClipMode currently set above
!  ! Perform split eta-upper
!  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
!                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
!
!  ! second split: eta lower
!  ! restore values to allow for correct back-trafo of intervals (required for intersectionpoint)
!  nXiClip     =tmpnXi
!  nEtaClip    =tmpnEta
!  LineNormVec =tmpLineNormVec
!  ClipMode    =tmpClipMode
!
!  ! set mapping array
!  EtaArray(:,nEtaClip)=(/EtaMin,EtaSplit/)
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  BezierControlPoints2D_temp=0.
!  PlusEta=1.0+EtaSplit
!  MinusEta=1.0-EtaSplit
!  ! compute the required stuff || pseudo Horner or precomputation
!  Etaup(0)=1.0
!  ! caution, here the indicies are switched from n-j to j  for **down
!  Etadown(0)=1.0
!  DO l=1,NGeo
!    Etaup   (l)     =Etaup   (l-1)*PlusEta
!    Etadown (l)     =Etadown (l-1)*MinusEta
!  END DO ! l=0,NGeo
!  DO p=0,NGeo
!    DO l=0,p
!      EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!    END DO ! l=0,p
!  END DO ! p=0,NGeo
!
!  DO q=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        !BezierControlPoints2D_temp(:,q,p)=&
!        !BezierControlPoints2D_temp(:,q,p)+&
!        !!BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!        !BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!        !                                 *arrayNchooseK(p,l) &
!        !                                 *(1.+EtaSplit)**l       &
!        !                                 *(1.-EtaSplit)**(p-l)
!        !DEBUG: optimize this !
!!            BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!!                                             +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!!                                             *(PlusEta**l)*(MinusEta**(p-l))
!!
!        BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                         +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
!      END DO
!    END DO
!  END DO
!  ! BOTTOM (mirrored Bernstein Basis evaluation)
!  ! s = (smin+1)/(smax+1) for [-1, +1]
!  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
!  IF(EtaMin.NE.-1.0)THEN
!    BezierControlPoints2D_temp2=0.
!    PlusEta=(EtaMin+1.0)/(EtaSplit+1.0)
!    ! MinusXi= 2.0*(1-PlusXi)-1
!    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!    MinusEta=2.0*PlusEta
!    ! PlusXi=1+ 2.0*(1-s)-1
!    PlusEta=2.0-2.0*PlusEta
!    ! compute the required stuff || pseudo Horner or precomputation
!    Etaup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    Etadown(0)=1.0
!    DO l=1,NGeo
!      Etaup   (l)     =Etaup   (l-1)*PlusEta
!      Etadown (l)     =Etadown (l-1)*MinusEta
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!        !  BezierControlPoints2D_temp2(:,q,p)=&
!        !  BezierControlPoints2D_temp2(:,q,p)+&
!        !  !BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!        !  BezierControlPoints2D_temp(:,q,NGeo-l)*(1./(2.**p))                     &
!        !                                        *arrayNchooseK(p,l)               &
!        !                                        *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
!        !                                        *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
!        !DEBUG: optimize this !
!          !BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
!          !                                 +BezierControlPoints2D_temp  (:,q,NGeo-l)*FacNchooseK(p,l) &
!          !                                 *(PlusEta**l)*(MinusEta**(p-l))
!          BezierControlPoints2D_temp2(:,q,NGeo-p)=BezierControlPoints2D_temp2(:,q,NGeo-p)             &
!                                                 +BezierControlPoints2D_temp(:,q,NGeo-l)*EtaBuf(p,l)
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
!!          IPWRITE(UNIT_stdout,*) ' split eta-lower '
!!          CALL OutputBezierControlPoints(BezierControlPoints2D_in=BezierControlPoints2D_temp2)
!!        END IF
!!      END IF
!#endif /*CODE_ANALYZE*/
!  ! HERE, ClipMode currently set above
!  ! Perform split eta-lower
!  CALL BezierClipRecursive(ClipMode,BezierControlPoints2D_temp2,LineNormVec,PartTrajectory,lengthPartTrajectory&
!                 ,iClipIter,nXiClip,nEtaClip,nInterSections,iPart,SideID)
!
!  ! after recursive steps, we are done!
!ELSE  ! no split necessary, only a clip
!
!  ! set mapping array
!  EtaArray(:,nEtaClip)=(/EtaMin,EtaMax/)
!  ! TOP, Bernstein polynomial B(n,k,x) = (1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k)
!  IF(EtaMax.NE.1.0)THEN
!    BezierControlPoints2D_temp=0.
!    PlusEta=1.0+EtaMax
!    MinusEta=1.0-EtaMax
!    ! compute the required stuff || pseudo Horner or precomputation
!    Etaup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    Etadown(0)=1.0
!    DO l=1,NGeo
!      Etaup   (l)     =Etaup   (l-1)*PlusEta
!      Etadown (l)     =Etadown (l-1)*MinusEta
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!!              BezierControlPoints2D_temp(:,q,p)=&
!!              BezierControlPoints2D_temp(:,q,p)+&
!!              !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
!!              BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
!!                                               *arrayNchooseK(p,l) &
!!                                               *(1.+Etamax)**l       &
!!                                               *(1.-Etamax)**(p-l)
!!             !DEBUG: optimize this !
!!              BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!!                                               +BezierControlPoints2D     (:,q,l)*FacNchooseK(p,l) &
!!                                               *(PlusEta**l)*(MinusEta**(p-l))
!!
!          BezierControlPoints2D_temp(:,q,p)=BezierControlPoints2D_temp(:,q,p)                  &
!                                           +BezierControlPoints2D     (:,q,l)*EtaBuf(p,l)
!        END DO
!      END DO
!    END DO
!    BezierControlPoints2D=BezierControlPoints2D_temp
!  END IF
!  ! BOTTOM (mirrored Bernstein Basis evaluation)
!  ! s = (smin+1)/(smax+1) for [-1, +1]
!  ! s = 2*(1-s)-1         for mirror input for bernstein (1-) and trafo [-1, +1] to [0, 1]
!  IF(EtaMin.NE.-1.0)THEN
!    BezierControlPoints2D_temp=0.
!    PlusEta=(EtaMin+1.0)/(EtaMax+1.0)
!    ! MinusXi= 2.0*(1-PlusXi)-1
!    ! MinusXi= 1.0-2.0*(1.0-s)+1.0
!    MinusEta=2.0*PlusEta
!    ! PlusXi=1+ 2.0*(1-s)-1
!    PlusEta=2.0-2.0*PlusEta
!    ! compute the required stuff || pseudo Horner or precomputation
!    Etaup(0)=1.0
!    ! caution, here the indicies are switched from n-j to j  for **down
!    Etadown(0)=1.0
!    DO l=1,NGeo
!      Etaup   (l)     =Etaup   (l-1)*PlusEta
!      Etadown (l)     =Etadown (l-1)*MinusEta
!    END DO ! l=0,NGeo
!    DO p=0,NGeo
!      DO l=0,p
!        EtaBuf(p,l)=EtaUp(l)*EtaDown(p-l)*FacNchooseK(p,l)
!      END DO ! l=0,p
!    END DO ! p=0,NGeo
!
!    DO q=0,NGeo
!      DO p=0,NGeo
!        DO l=0,p
!          !BezierControlPoints2D_temp(:,q,p)=&
!          !BezierControlPoints2D_temp(:,q,p)+&
!          !!BezierControlPoints2D(:,NGeo-l)*B(p-1,l-1,1-2*((Smin+1)/(Smax+1)))
!          !BezierControlPoints2D     (:,q,NGeo-l)*(1./(2.**p))                     &
!          !                                      *arrayNchooseK(p,l)               &
!          !                                      *(1+2*((EtaMin+1)/(EtaMax+1)))**(l-1) &
!          !                                      *(1-2*((EtaMin+1)/(EtaMax+1)))**(p-l)
!          !DEBUG: optimize this !
!          !BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)        &
!          !                                 +BezierControlPoints2D(:,q,NGeo-l)*FacNchooseK(p,l) &
!          !                                 *(PlusEta**l)*(MinusEta**(p-l))
!          BezierControlPoints2D_temp(:,q,NGeo-p)=BezierControlPoints2D_temp(:,q,NGeo-p)             &
!                                                 +BezierControlPoints2D     (:,q,NGeo-l)*EtaBuf(p,l)
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
!END SUBROUTINE PerformEtaClip
!
!END MODULE MOD_EtaClip
