#include "boltzplatz.h"

MODULE MOD_Particle_Surfaces
!===================================================================================================================================
! Contains subroutines to build the requiered data to track particles on (curviilinear) meshes, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleSurfaces
  MODULE PROCEDURE InitParticleSurfaces
END INTERFACE

INTERFACE FinalizeParticleSurfaces
  MODULE PROCEDURE FinalizeParticleSurfaces
END INTERFACE

INTERFACE GetBezierControlPoints3D
  MODULE PROCEDURE GetBezierControlPoints3D
END INTERFACE

INTERFACE CalcNormAndTangBilinear
  MODULE PROCEDURE CalcNormAndTangBilinear
END INTERFACE

INTERFACE CalcNormAndTangBezier
  MODULE PROCEDURE CalcNormAndTangBezier
END INTERFACE

INTERFACE GetSideSlabNormalsAndIntervals
  MODULE PROCEDURE GetSideSlabNormalsAndIntervals
END INTERFACE

INTERFACE GetElemSlabNormalsAndIntervals
  MODULE PROCEDURE GetElemSlabNormalsAndIntervals
END INTERFACE

INTERFACE GetBezierSampledAreas
  MODULE PROCEDURE GetBezierSampledAreas
END INTERFACE

INTERFACE EvaluateBezierPolynomialAndGradient
  MODULE PROCEDURE EvaluateBezierPolynomialAndGradient
END INTERFACE
  

PUBLIC::InitParticleSurfaces, FinalizeParticleSurfaces, GetBezierControlPoints3D, GetSideSlabNormalsAndIntervals, &
        GetElemSlabNormalsAndIntervals,GetBezierSampledAreas,EvaluateBezierPolynomialAndGradient

PUBLIC::CalcNormAndTangBilinear, CalcNormAndTangBezier

!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_vars
USE MOD_Preproc
USE MOD_Globals_Vars,               ONLY:epsMach
USE MOD_Mesh_Vars,                  ONLY:nSides,ElemToSide,NGeo,nBCSides,nSides
USE MOD_ReadInTools,                ONLY:GETREAL,GETINT,GETLOGICAL
USE MOD_Particle_Mesh_Vars,         ONLY:PartBCSideList
USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,     ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
#endif /*CODE_ANALYZE*/
!USE MOD_Particle_SFC_Vars,          ONLY:whichBoundBox
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: tmp,iSide!,iBCSide
CHARACTER(LEN=2)                :: dummy                         
!===================================================================================================================================

IF(ParticleSurfaceInitIsDone) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES ...!'

BezierNewtonAngle     = GETREAL('BezierNewtonAngle','1.570796326')! 1°=0.01754 (in rad)

!BezierHitEpsBi        = GETREAL('BezierHitEpsBi','1e-12')
!BezierHitEpsBi        = 1.0+BezierHitEpsBi
!BezierHitEpsBi=1.0+SQRT(EPSILON(0.0))
!BezierHitEpsBi=1.000800
BezierClipTolerance   = GETREAL('BezierClipTolerance','1e-8')
BezierSplitLimit      = GETREAL('BezierSplitLimit','0.6')
BezierSplitLimit      = 2.*BezierSplitLimit
BezierClipMaxIter     = GETINT('BezierClipMaxIter','100')

epsilontol            = GETREAL('epsilontol','0.')
! if nothing is entered, than a default value is used
! for tolerance issuses see, e.g. Haselbxxx PIC Tracking Paper
! epsilon approx 100*tolerance of the algorithm
IF(ALMOSTZERO(epsilontol)) epsilontol=100.*epsMach
MinusEps              = -epsilontol
OnePlusEps            = 1.0 + 100.*epsilontol
OneMinusEps           = 1.0 - epsilontol

BezierClipHit         = GETREAL('BezierClipHit','0.')
IF(ALMOSTZERO(BezierClipHit)) BezierClipHit=100.*BezierClipTolerance
BezierClipHit         = 1.0+BezierClipHit
tmp=2*(NGeo+1)
WRITE(dummy,'(I2.2)') tmp
BezierClipMaxIntersec = GETINT('BezierClipMaxIntersec',dummy)

IF(DoRefMapping)THEN
  !MultipleBCs    = GETLOGICAL('MultibleBCs',".FALSE.")
  ALLOCATE(PartBCSideList(1:nSides))
  PartBCSideList(:) = -1
  DO iSide=1,nBCSides
    PartBCSideList(iSide)=iSide
  END DO 
  !nTotalBCSides=nBCSides
 ! iBCSide=nBCSides
 ! DO iSide=nBCSides+1,nSides
 !   IF(BC(iSide).EQ.1) THEN
 !     iBCSide=iBCSide+1
 !     PartBCSideList(iSide)=iBCSide
 !   END IF
 ! END DO 
 ! nTotalBCSides=iBCSide
END IF

#ifdef CODE_ANALYZE
rBoundingBoxChecks=0.
rPerformBezierClip=0.
rPerformBezierNewton=0.
rTotalBBChecks    =0.
rTotalBezierClips =0.
rTotalBezierNewton=0.
#endif /*CODE_ANALYZE*/
!! ElemBaryNGeo are required for particle mapping| SingleParticleToExactElem
!IF(.NOT.DoRefMapping)THEN
!!   ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:PP_nElems) &
!!           ,slenXiEtaZetaBasis(1:6,1:PP_nElems) &
!!           ,ElemRadiusNGeo(1:PP_nElems)         &
!!           ,ElemBaryNGeo(1:3,1:PP_nElems)       )
!!   CALL BuildElementBasis()
!ELSE
!  !whichBoundBox = GETINT('PartSFC-BoundBox','1')
!END IF

ALLOCATE( locAlpha(1:BezierClipMaxIntersec) &
        , locXi   (1:BezierClipMaxIntersec) &
        , locEta  (1:BezierClipMaxIntersec) )
ALLOCATE( XiArray (1:2,1:BezierClipMaxIter) &
        , EtaArray(1:2,1:BezierClipMaxIter) )


! moved into mesh init
! construct connections to neighbor elems
!ALLOCATE( neighborElemID    (1:6,1:PP_nElems) &
!        , neighborlocSideID (1:6,1:PP_nElems) )
!neighborElemID=-1
!neighborlocSideID=-1

!DO iElem=1,PP_nElems
!  DO ilocSide=1,6
!    flip = ElemToSide(E2S_FLIP,ilocSide,iElem)
!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(flip.EQ.0)THEN
!      ! SideID of slave
!      neighborlocSideID(ilocSide,iElem)=SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
!      neighborElemID   (ilocSide,iElem)=SideToElem(S2E_NB_ELEM_ID,SideID)
!    ELSE
!      ! SideID of master
!      neighborlocSideID(ilocSide,iElem)=SideToElem(S2E_LOC_SIDE_ID,SideID)
!      neighborElemID   (ilocSide,iElem)=SideToElem(S2E_ELEM_ID,SideID)
!    END IF
!  END DO ! ilocSide
!END DO ! Elem

ParticleSurfaceInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleSurfaces

SUBROUTINE FinalizeParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(SideType)
!SDEALLOCATE(BiLinearCoeff)
SDEALLOCATE(SideNormVec)
SDEALLOCATE(SideDistance)
SDEALLOCATE(BezierControlPoints3D)
!SDEALLOCATE(SuperSampledBiLinearCoeff)
SDEALLOCATE(SideSlabNormals)
SDEALLOCATE(SideSlabIntervals)
SDEALLOCATE(ElemSlabNormals)
SDEALLOCATE(ElemSlabIntervals)
SDEALLOCATE(BoundingBoxIsEmpty)
SDEALLOCATE(ElevationMatrix)
SDEALLOCATE(BezierControlPoints3DElevated)
SDEALLOCATE(locAlpha)
SDEALLOCATE(locXi)
SDEALLOCATE(locEta)
SDEALLOCATE(XiArray)
SDEALLOCATE(EtaArray)
SDEALLOCATE(Vdm_Bezier)
SDEALLOCATE(sVdm_Bezier)
SDEALLOCATE(D_Bezier)
SDEALLOCATE(arrayNChooseK)
SDEALLOCATE(BezierControlPoints3DElevated)
SDEALLOCATE(FacNchooseK)
SDEALLOCATE(SideType)
SDEALLOCATE(BezierSampleXi)
SDEALLOCATE(SurfMeshSubSideData)
SDEALLOCATE(SurfMeshSideAreas)
SDEALLOCATE(BCdata_auxSF)
!SDEALLOCATE(gElemBCSide)
ParticleSurfaceInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleSurfaces


SUBROUTINE CalcNormAndTangBilinear(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Globals,                              ONLY:CROSSNORM,UNITVECTOR
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,INTENT(OUT)                       :: nVec(3)
REAL,INTENT(OUT),OPTIONAL              :: tang1(3), tang2(3)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: a,b
!================================================================================================================================


b=xi*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)  & 
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
   +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)   &
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)+BezierControlPoints3D(:,0   ,NGeo,SideID) )

a=eta*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
    +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) )

nVec=CROSSNORM(a,b)
IF(PRESENT(tang1)) tang1=UNITVECTOR(a)
IF(PRESENT(tang2)) tang2=UNITVECTOR(b)

END SUBROUTINE CalcNormAndTangBilinear


SUBROUTINE CalcNormAndTangBezier(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Globals,                              ONLY:CROSSNORM,UNITVECTOR
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,               ONLY:SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,DIMENSION(3),INTENT(OUT)          :: nVec
REAL,DIMENSION(3),INTENT(OUT),OPTIONAL :: tang1,tang2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(2,3)                    :: gradXiEta
!================================================================================================================================

! caution we require the formula in [0;1]
CALL EvaluateBezierPolynomialAndGradient((/xi,eta/),NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Gradient=gradXiEta)
nVec =CROSSNORM(gradXiEta(1,:),gradXiEta(2,:))
gradXiEta(2,:)=CROSSNORM(nVec,gradXiEta(1,:))
IF(PRESENT(tang1)) tang1=UNITVECTOR(gradXiEta(1,:))
IF(PRESENT(tang2)) tang2=gradXiEta(2,:)

END SUBROUTINE CalcNormAndTangBezier


SUBROUTINE EvaluateBezierPolynomialAndGradient(Xi,N_in,iSize,BezierControlPoints,dBezierControlPoints,Point,Gradient) 
!===================================================================================================================================
! evaluate the BezierPolynomial and/or gradient by means of direct evaluation and help of a pseudo-horner scheme
! faster for the evaluation of gradients
! results depend on mode
! mode=1 compute/interpolate point in patch at xi
! mode=2 compute gradient in xi at xi
! mode=3 point and gradient || both
! range: [-1,1]
! info for gradient
! grad(1,:) are all gradients in u direction,...
! book:
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! 1) computation of Bezier Polynomial is from Farin, Chap. 5, p. 57, EQ. 1
!    precompute t^i (1-t)^(n-i) by means of Horner Schema and store values direcly sorted 
! 2) point value
! 3) evaluate gradients
! optimized in terms of floating point operations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,               ONLY:facNchooseK,D_Bezier
USE MOD_TimeDisc_Vars,                        ONLY: iter
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)              :: Xi(2)                 ! value of the plane (u,v)
INTEGER,INTENT(IN)           :: N_in
INTEGER,INTENT(IN)           :: iSize                 ! depending if two or three dimensional
REAL,INTENT(IN)              :: BezierControlPoints(1:iSize,0:N_in,0:N_in)
REAL,INTENT(IN),OPTIONAL     :: dBezierControlPoints(1:2,1:iSize,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT),OPTIONAL    :: Point(1:iSize)        ! point value
REAL,INTENT(OUT),OPTIONAL    :: Gradient(1:2,1:iSize) ! gradient in u,v
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: Mode,iMode
INTEGER                      :: p,q,m
REAL                         :: MinusXi,PlusXI,MinusEta,PlusEta
REAL                         :: xiup(0:N_in),etaup(0:N_in),xidown(0:N_in),etadown(0:N_in)
REAL                         :: xitild(0:N_in), etatild(0:N_in)
!===================================================================================================================================

Mode=0
IF(PRESENT(Point)) Mode=1
iMode=0
IF(Present(Gradient)) iMode=2
Mode=Mode+iMode

M=N_in-1
MinusXi =1.0-xi(1)
PlusXi  =1.0+xi(1)
MinusEta=1.0-xi(2)
PlusEta =1.0+xi(2)

! compute the required stuff || pseudo Horner or precomputation
xiup(0)=1.0
etaup(0)=1.0
! caution, here the indicies are switched from n-j to j  for **down
xidown(N_in)=1.0
etadown(N_in)=1.0
DO p=1,N_in
  xiup   (p)     =xiup   (p-1)*PlusXi
  etaup  (p)     =etaup  (p-1)*PlusEta
  xidown (N_in-p)=xidown (N_in-p+1)*MinusXi
  etadown(N_in-p)=etadown(N_in-p+1)*MinusEta
END DO ! p

! only once
DO p=0,N_in
  xitild (p)=xiup (p)*xidown (p)
  etatild(p)=etaup(p)*etadown(p)
END DO ! p=0,N_in

IF((Mode.EQ.1).OR.(Mode.EQ.3))THEN
  Point=0.
  DO q=0,N_in
    DO p=0,N_in
      ! derivative in xi: d/dxi = u
      !Point(:)=Point(:)+BezierControlPoints(:,p,q)              &
      !                 *facNchooseK(N_In,p)*xiup(p) *xidown(p)  &
      !                 *facNChooseK(N_In,q)*etaup(q)*etadown(q)
      Point(:)=Point(:)+BezierControlPoints(:,p,q)     &
                       *facNchooseK(N_In,p)*xitild (p) &
                       *facNChooseK(N_In,q)*etatild(q)
    END DO ! p
  END DO ! q
  !Point=Point
END IF

IF(Mode.GE.2)THEN ! gradient
  gradient=0.
  IF(PRESENT(dBezierControlPoints))THEN
    DO q=0,N_in
      DO p=0,N_in
        gradient(:,:)=gradient(:,:)+dBezierControlPoints(:,:,p,q)  &
                                   *facNchooseK(N_In,p)*xitild (p) &
                                   *facNChooseK(N_In,q)*etatild(q)
      END DO ! p
    END DO ! q
  ELSE
    DO q=0,N_in
      DO p=0,M
        ! derivative in xi: d/dxi = u
        gradient(1,:)=gradient(1,:)+(BezierControlPoints(:,p+1,q)-BezierControlPoints(:,p,q)) &
                                   *facNchooseK(M,p)*xiup(p)*xidown(p+1)                      &
                                   *facNChooseK(N_in,q)*etatild(q)
        ! derivative in eta: d/deta = v ! caution - exchange indicies
        gradient(2,:)=gradient(2,:)+(BezierControlPoints(:,q,p+1)-BezierControlPoints(:,q,p)) &
                                   *facNchooseK(N_in,q)*xitild(q)                             &
                                   *facNchooseK(M,p)*etaup(p)*etadown(p+1) 
      END DO ! p
    END DO ! q
    ! caution, maybe one 1/(b-a) missing??
    gradient=0.5*REAL(N_in)*gradient
  END IF
  ! apply gradient
  !dBezierControlPoints=0.
  !DO nn=1,iSize
  !  DO dd=1,2 
  !    DO q=0,N_In
  !      CycPQ(2)=q
  !      DO p=0,N_In
  !        CycPQ(1)=p
  !        ! Matrix-vector multiplication
  !        Cyc=CycPQ
  !        DO l=0,N_in
  !          Cyc(dd)=l  ! d/dxi_dd
  !          dBezierControlPoints(dd,nn,p,q) = dBezierControlPoints(dd,nn,p,q) + &
  !                                            D_Bezier(CycPQ(dd),l)*BezierControlPoints(nn,Cyc(1),Cyc(2))
  !        END DO ! l=0,N_in
  !      END DO ! p=0,N_in
  !    END DO ! p=0,N_in
  !  END DO ! dd=1,iSize
  !END DO ! nn=1,iSize
END IF

END SUBROUTINE EvaluateBezierPolynomialAndGradient


SUBROUTINE GetBezierControlPoints3D(XCL_NGeo,iElem)
!===================================================================================================================================
! computes the nodes for Bezier Control Points for [P][I][C] [A]daptive [S]uper [S]ampled Surfaces [O]perations
! the control points (coeffs for bezier basis) are calculated using the change basis subroutine that interpolates the points 
! from the curved lagrange basis geometry (pre-computed inverse of vandermonde is required)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:ElemToSide,NGeo
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,sVdm_Bezier
USE MOD_Mesh_Vars,                ONLY:nBCSides,nInnerSides,nMPISides_MINE
USE MOD_ChangeBasis,              ONLY:ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: lastSideID,SideID
INTEGER                           :: p,q
REAL                              :: tmp(3,0:NGeo,0:NGeo)  

!===================================================================================================================================
!print*,"SIZE(BezierControlPoints)"
!print*,SIZE(BezierControlPoints)
!print*,"SHAPE(BezierControlPoints)"
!print*,SHAPE(BezierControlPoints)
! BCSides, InnerSides and MINE MPISides are filled
lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
!IF(DoRefMapping) lastSideID  = nBCSides
!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) XI_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,XI_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) XI_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,XI_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:),tmp)
    !print*,'ixi'
    BezierControlPoints3D(:,:,:,SideID)=tmp
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) ETA_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:),BezierControlPoints3D(1:3,:,:,sideID))
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
   END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 4.) ETA_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,NGeo-p,q)
      END DO !p
    END DO !q
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 5.) ZETA_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 6.) ZETA_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo),BezierControlPoints3D(1:3,:,:,sideID))
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
END SUBROUTINE GetBezierControlPoints3D

SUBROUTINE GetSideSlabNormalsAndIntervals(NGeo,SideID)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervalls)
! see article:  
!    author = {Shyue-wu Wang and Zen-chung Shih and Ruei-chuan Chang},                    
!    title = {An Efficient and Stable Ray Tracing Algorithm for Parametric Surfaces},     
!    year = {2001},
! original article: oriented-slab box (cartesian bounding box)
!   author = {Yen, Jonathan and Spach, Susan and Smith, Mark and Pulleyblank, Ron},       
!   title = {Parallel Boxing in B-Spline Intersection},                                   
!   issue_date = {January 1991},
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_Globals_Vars,    ONLY:EpsMach
USE MOD_Preproc
USE MOD_Particle_Surfaces_Vars,   ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,   ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,BezierElevation
USE MOD_Particle_Surfaces_Vars,   ONLY: BezierEpsilonBilinear
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,   ONLY: SideBoundingBoxVolume
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: SideID,NGeo
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                           :: lastSideID,flip,SideID
INTEGER            :: p,q
!REAL                              :: tmp(3,0:NGeo,0:NGeo)  
REAL               :: skalprod(3),dx,dy,dz,dMax,dMin,w,h,l
LOGICAL            :: SideIsCritical
!===================================================================================================================================


IF(BezierElevation.EQ.0)THEN
  BezierControlPoints3DElevated=BezierControlPoints3D
ELSE
  CALL GetBezierControlPoints3DElevated(NGeo,SideID)
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors
!-----------------------------------------------------------------------------------------------------------------------------------
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
SideSlabNormals(:,1,SideID)=BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideID)                      &
                           -BezierControlPoints3DElevated(:,0,0,SideID)                                         &
                           +BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideID)   &
                           -BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideID)
SideSlabNormals(:,1,SideID)=SideSlabNormals(:,1,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
SideSlabNormals(:,2,SideID)=BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideID)                      &
                           -BezierControlPoints3DElevated(:,0,0,SideID)                                         &
                           +BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideID)   &
                           -BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideID)

!print*,"length(SideSlabNormals(:,1,SideID))",DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID))
!fehlt das?
!SideSlabNormals(:,2,SideID)=SideSlabNormals(:,2,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
SideSlabNormals(:,2,SideID)=CROSSNORM(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID))
!print*,"length(SideSlabNormals(:,2,SideID))",DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID))

!SideSlabNormals(:,2,SideID)=SideSlabNormals(:,2,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
! n_3=n_1 x n_2
SideSlabNormals(:,3,SideID)=CROSSNORM(SideSlabNormals(:,2,SideID),SideSlabNormals(:,1,SideID))

! check vector length=1
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 2 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,3,SideID),SideSlabNormals(:,3,SideID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 3 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,3,SideID),SideSlabNormals(:,3,SideID)))

! check perpendicularity
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 and 2 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,3,SideID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,3,SideID))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,3,SideID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 2 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,3,SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervalls beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!SideSlabIntervals(x- x+ y- y+ z- z+, SideID)


! Intervall beta_1
!print*,"SideID",SideID
SideSlabIntervals(:, SideID)=0.

DO q=0,NGeo+BezierElevation
  DO p=0,NGeo+BezierElevation
    IF((p.EQ.0).AND.(q.EQ.0))CYCLE
    skalprod(1)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,1,SideID))
    skalprod(2)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,2,SideID))
    skalprod(3)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,3,SideID))
    IF    (skalprod(1).LT.0.)THEN
      SideSlabIntervals(1, SideID)=MIN(SideSlabIntervals(1,SideID),skalprod(1))
    ELSEIF(skalprod(1).GT.0.)THEN
      SideSlabIntervals(2, SideID)=MAX(SideSlabIntervals(2,SideID),skalprod(1))
    END IF
    IF    (skalprod(2).LT.0.)THEN
      SideSlabIntervals(3, SideID)=MIN(SideSlabIntervals(3,SideID),skalprod(2))
    ELSEIF(skalprod(2).GT.0.)THEN
      SideSlabIntervals(4, SideID)=MAX(SideSlabIntervals(4,SideID),skalprod(2))
    END IF
    IF    (skalprod(3).LT.0.)THEN
      SideSlabIntervals(5, SideID)=MIN(SideSlabIntervals(5,SideID),skalprod(3))
    ELSEIF(skalprod(3).GT.0.)THEN
      SideSlabIntervals(6, SideID)=MAX(SideSlabIntervals(6,SideID),skalprod(3))
    END IF
  END DO !p
END DO !q
!print*,"SideID",SideID," is planar:",SideIsPlanar(SideID)
dx=ABS(ABS(SideSlabIntervals(2, SideID))-ABS(SideSlabIntervals(1, SideID)))
dy=ABS(ABS(SideSlabIntervals(4, SideID))-ABS(SideSlabIntervals(3, SideID)))
dz=ABS(ABS(SideSlabIntervals(6, SideID))-ABS(SideSlabIntervals(5, SideID)))

!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) Is Side critical? (particle path parallel to the larger surface, therefore numerous intersections are possilbe)
! from Wang et al. 2001 An Efficient and Stable ray tracing algorithm for parametric surfaces
! critical sides: h << w <= l
! h: height
! w: width (short side)
! l: length (long side)
!-----------------------------------------------------------------------------------------------------------------------------------

h=dy
w=MIN(dx,dz)
l=MAX(dx,dz)

IF((h.LT.w*1.E-2))THEN!.AND.(w.LE.l))THEN ! critical case possible: second condition is obsolete
  SideIsCritical=.TRUE.
ELSE
  SideIsCritical=.FALSE.
END IF


!-----------------------------------------------------------------------------------------------------------------------------------
! 4.) determine is side is planar -> "BoundingBoxIsEmpty(SideID)"
!     this results also in the decision whether a side is also considered flat or bilinear! 
!-----------------------------------------------------------------------------------------------------------------------------------

dMax=MAX(dx,dy,dz)
dMin=MIN(dx,dy,dz)
IF(dx/dMax.LT.BezierEpsilonBilinear)THEN
  CALL Abort(&
__STAMP__&
,'Bezier side length is degenerated. dx/dMax.LT.BezierEpsilonBilinear ->',0,dx/dMax)
END IF
IF(dy/dMax.LT.BezierEpsilonBilinear)THEN
  SideSlabIntervals(3:4, SideID)=0.
  dy=0.
END IF
IF(dz/dMax.LT.BezierEpsilonBilinear)THEN
  CALL Abort(&
__STAMP__&
,'Bezier side length is degenerated. dz/dMax.LT.BezierEpsilonBilinear ->',0,dz/dMax)
END IF

IF(dx*dy*dz.LT.0) THEN
  IPWRITE(UNIT_stdOut,*) ' Warning, no bounding box'
  IF(dx*dy*dz.LT.0) CALL Abort(&
__STAMP__&
,'A bounding box (for sides) is negative!?. dx*dy*dz.LT.0 ->',0,(dx*dy*dz))
END IF

IF(ALMOSTZERO(dx*dy*dz))THEN ! bounding box volume is approx zeros
  BoundingBoxIsEmpty(SideID)=.TRUE.
ELSE
  BoundingBoxIsEmpty(SideID)=.FALSE.
END IF

#ifdef CODE_ANALYZE
SideBoundingBoxVolume(SideID)=dx*dy*dz
#endif /*CODE_ANALYZE*/
END SUBROUTINE GetSideSlabNormalsAndIntervals


SUBROUTINE GetElemSlabNormalsAndIntervals(NGeo,ElemID)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervalls)
! of each element. This routine must be called after GetSideSlabNormalsAndIntervals(...), because the elevation takes place there
! see article:  
!    author = {Shyue-wu Wang and Zen-chung Shih and Ruei-chuan Chang},                                                              
!    title = {An Efficient and Stable Ray Tracing Algorithm for Parametric Surfaces},                                               
!    year = {2001},
! original article: oriented-slab box (cartesian bounding box)
!   author = {Yen, Jonathan and Spach, Susan and Smith, Mark and Pulleyblank, Ron},                                                 
!   title = {Parallel Boxing in B-Spline Intersection},                                                                             
!   issue_date = {January 1991},
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,       ONLY:PartElemToSide,GEO,PartBCSideList,RefMappingEps 
USE MOD_Particle_Surfaces_Vars,   ONLY:ElemSlabNormals,ElemSlabIntervals,BezierControlPoints3DElevated,BezierElevation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID,NGeo
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,iLocSide,SideID,SideIDOrigin
REAL               :: skalprod(3),dx,dy,dz
!===================================================================================================================================

!BezierControlPoints(:,:,:,ElemID)
!ElemSlabNormals( x y z,1 2 3 , ElemID)
IF(GEO%nPeriodicVectors.GT.0)THEN
    CALL  Abort(&
  __STAMP__&
  ,' computation of wrong bounding box!')
  SWRITE(*,*) ' Computation of wrong bounding box in peridodic' 
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------
!ElemIsPlanar=.FALSE.

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors (use the first local element side)
!-----------------------------------------------------------------------------------------------------------------------------------
SideIDOrigin=PartElemToSide(E2S_SIDE_ID,1,ElemID)
ElemSlabNormals(:,0,ElemID)=SideIDOrigin
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
ElemSlabNormals(:,1,ElemID)= &
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideIDOrigin)- &
                BezierControlPoints3DElevated(:,0,0,SideIDOrigin)+&
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideIDOrigin)-&
                BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideIDOrigin)
ElemSlabNormals(:,1,ElemID)=ElemSlabNormals(:,1,ElemID)/SQRT(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
ElemSlabNormals(:,2,ElemID)= &
                BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideIDOrigin)   &
                -BezierControlPoints3DElevated(:,0,0,SideIDOrigin)+&
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideIDOrigin)&
                -BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideIDOrigin)
ElemSlabNormals(:,2,ElemID)=CROSSNORM(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID))
! n_3=n_1 x n_2
ElemSlabNormals(:,3,ElemID)=CROSSNORM(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,1,ElemID))

! check vector length=1
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 1 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID)))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,2,ElemID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 2 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,2,ElemID)))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,3,ElemID),ElemSlabNormals(:,3,ElemID))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 3 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,3,ElemID),ElemSlabNormals(:,3,ElemID)))

! check perpendicularity
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 1 and 2 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID))))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,3,ElemID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 1 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,3,ElemID))))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,3,ElemID)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Element slab normal 2 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,3,ElemID))))

!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervalls beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!ElemSlabIntervals(x- x+ y- y+ z- z+, ElemID)

ElemSlabIntervals(:,ElemID)=0.
DO iLocSide=1,6
  SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
  DO q=0,NGeo+BezierElevation
    DO p=0,NGeo+BezierElevation
      IF((p.EQ.0).AND.(q.EQ.0))CYCLE
      skalprod(1)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,1,ElemID))
      skalprod(2)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,2,ElemID))
      skalprod(3)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,3,ElemID))
      IF    (skalprod(1).LT.0.)THEN
        ElemSlabIntervals(1, ElemID)=MIN(ElemSlabIntervals(1,ElemID),skalprod(1))
      ELSEIF(skalprod(1).GT.0.)THEN
        ElemSlabIntervals(2, ElemID)=MAX(ElemSlabIntervals(2,ElemID),skalprod(1))
      END IF
      IF    (skalprod(2).LT.0.)THEN
        ElemSlabIntervals(3, ElemID)=MIN(ElemSlabIntervals(3,ElemID),skalprod(2))
      ELSEIF(skalprod(2).GT.0.)THEN
        ElemSlabIntervals(4, ElemID)=MAX(ElemSlabIntervals(4,ElemID),skalprod(2))
      END IF
      IF    (skalprod(3).LT.0.)THEN
        ElemSlabIntervals(5, ElemID)=MIN(ElemSlabIntervals(5,ElemID),skalprod(3))
      ELSEIF(skalprod(3).GT.0.)THEN
        ElemSlabIntervals(6, ElemID)=MAX(ElemSlabIntervals(6,ElemID),skalprod(3))
      END IF
    END DO !p
  END DO !q
END DO !iLocSide=1:6
dx=ABS(ABS(ElemSlabIntervals(2, ElemID))-ABS(ElemSlabIntervals(1, ElemID)))
dy=ABS(ABS(ElemSlabIntervals(4, ElemID))-ABS(ElemSlabIntervals(3, ElemID)))
dz=ABS(ABS(ElemSlabIntervals(6, ElemID))-ABS(ElemSlabIntervals(5, ElemID)))
ElemSlabIntervals(1,ElemID)=ElemSlabInterVals(1,ElemID)-RefMappingEps 
ElemSlabIntervals(2,ElemID)=ElemSlabInterVals(2,ElemID)+RefMappingEps 
ElemSlabIntervals(3,ElemID)=ElemSlabInterVals(3,ElemID)-RefMappingEps 
ElemSlabIntervals(4,ElemID)=ElemSlabInterVals(4,ElemID)+RefMappingEps 
ElemSlabIntervals(5,ElemID)=ElemSlabInterVals(5,ElemID)-RefMappingEps 
ElemSlabIntervals(6,ElemID)=ElemSlabInterVals(6,ElemID)+RefMappingEps 
IF(dx*dy*dz.LT.0) CALL Abort(&
__STAMP__&
,'A bounding box (for elements) is negative!?. dx*dy*dz.LT.0 ->',0,(dx*dy*dz))
!IF((dx*dy*dz).LT.GEO%Volume(ElemID))THEN
!  IPWRITE(*,*) 'Volume', dx*dy*dz
!  IPWRITE(*,*) 'DG-Volume', GEO%Volume(ElemID)
!  CALL Abort(&
!  __STAMP__&
!  'The bounding box is smaller than element! BoundingBox-volume:')
!END IF
END SUBROUTINE GetElemSlabNormalsAndIntervals


SUBROUTINE GetBezierControlPoints3DElevated(NGeo,SideID)
!===================================================================================================================================
! compute the elevated bezier control points (use pre-computed ElevationMatrix with binomial coefficicents)
! uses the tensor-product basis and combines two 1-D elevation steps
!===================================================================================================================================
! MODULES
!USE MOD_Globals
!USE MOD_Globals_Vars,    ONLY:EpsMach
!USE MOD_Preproc
!USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D,BoundingBoxIsEmpty
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierEpsilonBilinear
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,BezierControlPoints3DElevated,BezierElevation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: SideID,NGeo
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                           :: lastSideID,flip,SideID
INTEGER            :: p,q
!REAL                              :: tmp(3,0:NGeo,0:NGeo)  
!REAL               :: skalprod(3),dx,dy,dz,dMax,dMin,w,h,l
!LOGICAL            :: SideIsCritical
REAL               :: temp(1:3,0:NGeo+BezierElevation,0:NGeo)
!===================================================================================================================================
temp=0.
! p-direction
DO q=0,NGeo
  temp(:,:,q)=ElevateBezierPolynomial(NGeo,BezierControlPoints3D(:,:,q,SideID))
END DO
! q-direction
DO p=0,NGeo+BezierElevation
  BezierControlPoints3DElevated(:,p,:,SideID)=ElevateBezierPolynomial(NGeo,temp(:,p,:))
END DO
END SUBROUTINE GetBezierControlPoints3DElevated


SUBROUTINE GetBezierSampledAreas(SideID,BezierSampleN,ProjectionVector_opt,SurfMeshSubSideAreas,SurfMeshSideArea_opt &
                                ,SurfMeshSubSideVec_nOut_opt,SurfMeshSubSideVec_t1_opt,SurfMeshSubSideVec_t2_opt)
!===================================================================================================================================
! equidistanlty super-sampled bezier surface area and vector calculation
! --------------------------------------
! book: see also for general remarks
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! --------------------------------------
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,BezierControlPoints3D
USE MOD_Basis,                  ONLY:LegendreGaussNodesAndWeights
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Mesh_Vars,              ONLY:nBCSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: SideID,BezierSampleN
REAL,INTENT(IN),OPTIONAL            :: ProjectionVector_opt(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:BezierSampleN,1:BezierSampleN),INTENT(OUT)              :: SurfMeshSubSideAreas
REAL,INTENT(OUT),OPTIONAL                                                :: SurfMeshSideArea_opt
REAL,DIMENSION(1:3,1:BezierSampleN,1:BezierSampleN),INTENT(OUT),OPTIONAL :: SurfMeshSubSideVec_nOut_opt &
                                                                           ,SurfMeshSubSideVec_t1_opt &
                                                                           ,SurfMeshSubSideVec_t2_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q
INTEGER                                :: I,J,iSample,jSample
REAL                                   :: areaTotal,areaTotalAbs,area,deltaXi,tmp1,E,F,G,D
REAL                                   :: tmpI2,tmpJ2
REAL                                   :: BezierControlPoints2D(1:2,0:NGeo,0:NGeo)
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,DIMENSION(2,2)                    :: gradXiEta2D
REAL,DIMENSION(2)                      :: XiOut(2)
REAL,DIMENSION(2)                      :: Xi
REAL,DIMENSION(3)                      :: n1,n2
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
LOGICAL                                :: BezierSampleProjection
REAL,DIMENSION(3)                      :: ProjectionVector
REAL,DIMENSION(0:BezierSampleN)        :: BezierSampleXi
REAL,DIMENSION(1:3,1:BezierSampleN,1:BezierSampleN) :: SurfMeshSubSideVec_nOut,SurfMeshSubSideVec_t1,SurfMeshSubSideVec_t2
!===================================================================================================================================
IF (PRESENT(ProjectionVector_opt)) THEN
  BezierSampleProjection=.TRUE.
  ProjectionVector=ProjectionVector_opt
ELSE
  BezierSampleProjection=.FALSE.
  ProjectionVector=0. !dummy
END IF

Xi(1) =-SQRT(1./3.)
Xi(2) = SQRT(1./3.)
deltaXi=2.0/BezierSampleN
tmp1=deltaXi/2.0 !(b-a)/2
DO jSample=0,BezierSampleN
  BezierSampleXi(jSample)=-1.+deltaXi*jSample
END DO

!===================================================================================================================================

SurfMeshSubSideAreas=0.
SurfMeshSubSideVec_nOut=0.
SurfMeshSubSideVec_t1=0.
SurfMeshSubSideVec_t2=0.

IF(BezierSampleProjection)THEN
  ! transformation of bezier patch 3D->2D
  IF(ABS(ProjectionVector(3)).LT.epsilontol)THEN
    n1=(/ -ProjectionVector(2)-ProjectionVector(3)  , &
      ProjectionVector(1),ProjectionVector(1) /)
  ELSE
    n1=(/ ProjectionVector(3),ProjectionVector(3) ,&
      -ProjectionVector(1)-ProjectionVector(2) /)
  END IF
  n1=n1/SQRT(DOT_PRODUCT(n1,n1))
  n2(:)=CROSSNORM(ProjectionVector,n1)
  DO q=0,NGeo; DO p=0,NGeo
    BezierControlPoints2D(1,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID),n1)
    ! origin is (0,0,0)^T
    BezierControlPoints2D(2,p,q)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID),n2)
    ! origin is (0,0,0)^T
  END DO; END DO
END IF!(BezierSampleProjection)THEN

ALLOCATE(Xi_NGeo( 0:NGeo)  &
        ,wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

areaTotal=0.
areaTotalAbs=0.
DO jSample=1,BezierSampleN; DO iSample=1,BezierSampleN
  area=0.
  tmpI2=(BezierSampleXi(iSample-1)+BezierSampleXi(iSample))/2. ! (a+b)/2
  tmpJ2=(BezierSampleXi(jSample-1)+BezierSampleXi(jSample))/2. ! (a+b)/2
  ! ---------------------------------------
  ! calc integral
  DO I=0,NGeo; DO J=0,NGeo
    XiOut(1)=tmp1*Xi_NGeo(I)+tmpI2
    XiOut(2)=tmp1*Xi_NGeo(J)+tmpJ2
    IF(BezierSampleProjection)THEN
      ! get gradients
      CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,2,BezierControlPoints2D(1:2,0:NGeo,0:NGeo),Gradient=gradXiEta2D)
      ! calculate first fundamental form
      E=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
      F=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
      G=DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
    ELSE
      CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
                                              ,Gradient=gradXiEta3D)
      ! calculate first fundamental form
      E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
      F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
      G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
    END IF
    D=SQRT(E*G-F*F)
    area=area+tmp1*tmp1*D*wGP_NGeo(i)*wGP_NGeo(j)
  END DO; END DO
  ! ---------------------------------------
  IF (BezierSampleProjection .OR. &
      PRESENT(SurfMeshSubSideVec_nOut_opt) .OR. &
      PRESENT(SurfMeshSubSideVec_t1_opt) .OR. &
      PRESENT(SurfMeshSubSideVec_t2_opt) ) THEN
    CALL CalcNormAndTangBezier( nVec=SurfMeshSubSideVec_nOut(:,iSample,jSample) &
                              ,tang1=SurfMeshSubSideVec_t1(:,iSample,jSample) &
                              ,tang2=SurfMeshSubSideVec_t2(:,iSample,jSample) &
                              ,xi=tmpI2,eta=tmpJ2,SideID=SideID )
  END IF
  IF(BezierSampleProjection)THEN
    !!!IPWRITE(*,*) 'vec_nOut', SurfMeshSubSideVec_nOut(iSample,jSample)
    ! add facing sides and substract non-facing sides
    SurfMeshSubSideAreas(iSample,jSample) = SIGN(area,-DOT_PRODUCT(ProjectionVector,SurfMeshSubSideVec_nOut(:,iSample,jSample)))
    !!!IPWRITE(*,*)'sign-area', SurfMeshSubSideAreas(iSample,jSample)
  ELSE
    SurfMeshSubSideAreas(iSample,jSample) = area
  END IF
  areaTotal=areaTotal+SurfMeshSubSideAreas(iSample,jSample)
  areaTotalAbs=areaTotalAbs+area
  !!!IPWRITE(*,*)"areaTotal",areaTotal
  !!!IPWRITE(*,*)"areaTotalAbs",areaTotalAbs
END DO; END DO !jSample=1,BezierSampleN;iSample=1,BezierSampleN
  
DEALLOCATE(Xi_NGeo,wGP_NGeo)

IF (PRESENT(SurfMeshSideArea_opt)) SurfMeshSideArea_opt=areaTotal
IF (PRESENT(SurfMeshSubSideVec_nOut_opt)) SurfMeshSubSideVec_nOut_opt=SurfMeshSubSideVec_nOut
IF (PRESENT(SurfMeshSubSideVec_t1_opt)) SurfMeshSubSideVec_t1_opt=SurfMeshSubSideVec_t1
IF (PRESENT(SurfMeshSubSideVec_t2_opt)) SurfMeshSubSideVec_t2_opt=SurfMeshSubSideVec_t2

#ifdef CODE_ANALYZE
IPWRITE(*,*)" ===== SINGLE AREA ====="
IPWRITE(*,*)"SideID:",SideID,"areaTotal   =",areaTotal
IPWRITE(*,*)"SideID:",SideID,"areaTotalAbs=",areaTotalAbs
IPWRITE(*,*)" ============================"
#endif /*CODE_ANALYZE*/

END SUBROUTINE GetBezierSampledAreas


FUNCTION ElevateBezierPolynomial(NGeo,BezierPolynomial)
!===================================================================================================================================
! this function creates a new equidistantly distributed set of control points (bézier polynomial basis coefficients) based on the 
! control points "BezierPolynomial" with order "N" and elevates them to "ElevateBezierPolynomial" on "Xi_NGeo_elevated" with order 
! "p+BezierElevation"
! book: single step procedure (1D)
!   author = {Piegl, Les and Tiller, Wayne},                                                                       
!   title = {The NURBS Book (2Nd Ed.)},                                                                            
!   year = {1997},
! book: see also for general remarks
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
!===================================================================================================================================
! MODULES
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierElevation,ElevationMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo
REAL,INTENT(IN)    :: BezierPolynomial(1:3,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL               :: ElevateBezierPolynomial(1:3,0:NGeo+BezierElevation)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL            :: length
INTEGER            :: i,j,jStart,jEnd
!===================================================================================================================================
! algorithm originally from "The NURBS Book" by Les Piegl, Wayne Tiller (p.205)
ElevateBezierPolynomial = 0.
! the first and last points remain the same!
! edge points remain: P_0 and P_p
ElevateBezierPolynomial(:,0)                    = BezierPolynomial(:,0)
ElevateBezierPolynomial(:,NGeo+BezierElevation) = BezierPolynomial(:,NGeo)
! inner points change: P_1,...,P_p-1
DO i=1,NGeo+BezierElevation-1
  jStart = MAX(0,i-BezierElevation)
  jEnd   = MIN(NGeo,i)
  DO j=jStart,jEnd 
    ElevateBezierPolynomial(:,i)=ElevateBezierPolynomial(:,i)+ElevationMatrix(i,j)*BezierPolynomial(:,j)
  END DO
END DO
END FUNCTION ElevateBezierPolynomial


END MODULE MOD_Particle_Surfaces
