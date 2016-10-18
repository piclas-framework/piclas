#include "boltzplatz.h"

MODULE MOD_Metrics
!===================================================================================================================================
! calculate the Volume Metric terms
!           Metrics_fTilde(n=1:3,i,j,k,iElem) 
!           Metrics_gTilde(n=1:3,i,j,k,iElem) 
!           Metrics_hTilde(n=1:3,i,j,k,iElem) 
!       using a high precision mapping X_n(xi_i) of the Geometry using 
!   Chebyshev-Lobatto points with (NGeo+1) 1D points, saved in XCL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!   Per Element we do:
!   1.) calculate the gradient: compute the derivative of the mapping XCL_NGeo in (xi_1,xi_2,xi_3) direction,
!       using a Polynomial derivative Matrix
!   2.) for each direction n
!       a.) compute the nth vector and for each Chebyshev point (:,i,j,k)
!          (dXCL_n^1,dXCL_n^2,dXCL_n^3)^T=(X_l grad_xi (X_m) ) for n=1,2,3 and (n,m,l) cyclic
!       b.) interpolate the dXCL_n vector defined primarily on (NGeo+1)x(NGeo+1)x(Ngeo+1) Chebyshev-Lobatto points to
!             (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points and write to Ja_n(1:3,i,j,k) i,j,k=[0:N]
!       c.) compute the curl of vector Ja_n(1:3,i,j,k) using the derivative Matrix DCL_N [NxN]
!       d.) interpolate from (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points to  Gauss-Points (N+1)x(N+1)x(N+1) (exact!) 
!       e.) save Ja_n to the Metric terms
!             Metrics_fTilde(n,iElem)=Ja_n^1
!             Metrics_gTilde(n,iElem)=Ja_n^2
!             Metrics_hTilde(n,iElem)=Ja_n^3
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CalcMetrics
  MODULE PROCEDURE CalcMetrics
END INTERFACE

PUBLIC::CalcMetrics
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcMetrics(NodeCoords,XCL_NGeo_Out,dXCL_NGeo_out)
!===================================================================================================================================
! calculate the Volume Metric terms
!           Metrics_fTilde(n=1:3) 
!           Metrics_gTilde(n=1:3) 
!           Metrics_hTilde(n=1:3) 
!       using a high precision mapping X_n(xi_i) of the Geometry using 
!   Chebyshev-Lobatto points with (NGeo+1) 1D points, saved in XCL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!   Per Element we do:
!   1.) calculate the Jacobi Matrix: compute the derivative of the mapping XCL_NGeo in (xi_1,xi_2,xi_3) direction,
!       using a Polynomial derivative Matrix
!       b.) calculate the determinant of the Jacobi Matrix and save to sJ=1/detJac
!   2.) for each direction n
!       a.) compute the nth vector and for all (NGeo+1)x(NGeo+1)x(Ngeo+1) Chebyshev points, using cyclic indices
!          (dXCL_n^1,dXCL_n^2,dXCL_n^3)^T=(X_l grad_xi (X_m) ) for n=1,2,3 and (n,m,l) cyclic
!       b.) interpolate the dXCL_n vector defined primarily on (NGeo+1)x(NGeo+1)x(Ngeo+1) Chebyshev-Lobatto points to
!             (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points and write to Ja(1:3,1:3,i,j,k) i,j,k=[0:N]
!       c.) compute the curl of vector Ja(1:3,1:3,i,j,k) using the derivative Matrix DCL_N [NxN]
!       d.) interpolate from (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points to Gauss-Points (N+1)x(N+1)x(N+1) (exact!) 
!       e.) save Ja to the Metric terms
!             Metrics_fTilde(n,iElem)=Ja(1,n)
!             Metrics_gTilde(n,iElem)=Ja(2,n)
!             Metrics_hTilde(n,iElem)=Ja(3,n)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,               ONLY:NGeo,NGeoRef
USE MOD_Mesh_Vars,               ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,crossProductMetrics
USE MOD_Mesh_Vars,               ONLY:Face_xGP,normVec,surfElem,TangVec1,TangVec2
USE MOD_Mesh_Vars,               ONLY:nElems,sideID_minus_upper,nBCSides
USE MOD_Mesh_Vars,               ONLY:detJac_Ref
USE MOD_Mesh_Vars,               ONLY:crossProductMetrics
USE MOD_Mesh_Vars,               ONLY:nElems,sideID_minus_upper,nBCSides
USE MOD_Interpolation,           ONLY:GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_ChangeBasis,             ONLY:changeBasis3D,ChangeBasis3D_XYZ
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:NodeTypeG,NodeTypeGL,NodeTypeCL,NodeTypeVISU
#ifdef PARTICLES
USE MOD_Particle_Surfaces,       ONLY:GetBezierControlPoints3D
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#endif /*PARTICLES*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN)         :: PP_N
REAL,INTENT(IN)             :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT),OPTIONAL  :: XCL_Ngeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)      ! mapping X(xi) P\in Ngeo
REAL ,INTENT(INOUT),OPTIONAL :: dXCL_Ngeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)   ! jacobi matrix on CL Ngeo
!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,q,iElem
INTEGER :: dd,ee,ff
INTEGER :: nn,mm,ll
INTEGER :: Cyclic(5),CycIJK(3),Cyc(3)
! Jacobian on CL N and NGeoRef
REAL    :: DetJac_N( 1,0:PP_N,   0:PP_N,   0:PP_N)
REAL    :: tmp(      1,0:NgeoRef,0:NgeoRef,0:NgeoRef)
!REAL    :: tmp2(     1,0:Ngeo,0:Ngeo,0:Ngeo)
! interpolation points and derivatives on CL N
REAL    :: XCL_N(      3,  0:PP_N,0:PP_N,0:PP_N)          ! mapping X(xi) P\in N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: XCL_N_quad( 3,  0:PP_N,0:PP_N,0:PP_N)          ! mapping X(xi) P\in N
REAL    :: dXCL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)          ! jacobi matrix on CL N
REAL    :: dXCL_Ngeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on CL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

REAL    :: R_CL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! buffer for metric terms, uses XCL_N,dXCL_N
REAL    :: JaCL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! metric terms P\in N
REAL    :: JaCL_N_quad(3,3,0:PP_N,0:PP_N,0:PP_N)    ! metric terms P\in N
REAL    :: scaledJac(2)

! Polynomial derivativion matrices
REAL    :: DCL_NGeo(0:Ngeo,0:Ngeo)
REAL    :: DCL_N(   0:PP_N,0:PP_N)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNgeo_CLNgeo( 0:Ngeo   ,0:Ngeo)
REAL    :: Vdm_CLNGeo_NgeoRef(0:NgeoRef,0:Ngeo)
REAL    :: Vdm_NgeoRef_N(     0:PP_N   ,0:NgeoRef)
REAL    :: Vdm_CLNGeo_CLN(    0:PP_N   ,0:Ngeo)
REAL    :: Vdm_CLN_N(         0:PP_N   ,0:PP_N)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL,DIMENSION(0:NgeoRef,0:NgeoRef) :: Vdm_xi_Ref,Vdm_eta_Ref,Vdm_zeta_Ref
REAL,DIMENSION(0:PP_N   ,0:PP_N)    :: Vdm_xi_N  ,Vdm_eta_N  ,Vdm_zeta_N
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)
REAL    :: xiCL_N(0:PP_N)   ,wBaryCL_N(0:PP_N)
REAL    :: xi0(3),dxi(3),length(3)

#ifdef PARTICLES
INTEGER            :: iSide,lowerLimit
REAL               :: StartT2,BezierTime
#endif /*PARTICLES*/
REAL               :: StartT,EndT
!===================================================================================================================================


StartT=BOLTZPLATZTIME(MPI_COMM_WORLD)
#ifdef PARTICLES
BezierTime=0.
#endif

! Prerequisites
Metrics_fTilde=0.
Metrics_gTilde=0.
Metrics_hTilde=0.
! 
Cyclic=(/1,2,3,1,2/)

! Initialize Vandermonde and D matrices
! Only use modal Vandermonde for terms that need to be conserved as Jacobian if N_out>PP_N
! Always use interpolation for the rest!

! 1.a) NodeCoords: EQUI Ngeo to CLNgeo and CLN
CALL GetVandermonde(    Ngeo   , NodeTypeVISU, Ngeo    , NodeTypeCL, Vdm_EQNgeo_CLNgeo , modal=.FALSE.)

! 1.b) dXCL_Ngeo:
CALL GetDerivativeMatrix(Ngeo  , NodeTypeCL  , DCL_Ngeo)

! 1.c) Jacobian: CLNgeo to NgeoRef, CLNgeoRef to N
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNgeo_NgeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NgeoRef_N     , modal=.TRUE.)
CALL GetNodesAndWeights(NgeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)

! 1.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , PP_N    , NodeTypeCL, Vdm_CLNgeo_CLN    , modal=.FALSE.)
CALL GetDerivativeMatrix(PP_N  , NodeTypeCL  , DCL_N)

! 2.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
CALL GetVandermonde(    PP_N   , NodeTypeCL  , PP_N    , NodeType,   Vdm_CLN_N         , modal=.FALSE.)
CALL GetNodesAndWeights(PP_N   , NodeTypeCL  , xiCL_N  , wIPBary=wBaryCL_N)

! Outer loop over all elements
detJac_Ref=0.
DO iElem=1,nElems
  !1.a) Transform from EQUI_Ngeo to CL points on Ngeo and N
  !IF(interpolateFromTree)THEN
  !  xi0   =xiMinMax(:,1,iElem)
  !  length=xiMinMax(:,2,iElem)-xi0
  !  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,TreeCoords(:,:,:,:,ElemToTree(iElem)),XCL_Ngeo)
  !ELSE
    CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem)            ,XCL_Ngeo)
  !END IF
  CALL   ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,   XCL_Ngeo                             ,XCL_N)

  !1.b) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_NGeo(dd,nn,i,j,k))
  dXCL_NGeo=0.
  DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
    ! Matrix-vector multiplication
    DO ll=0,Ngeo
      dXCL_Ngeo(1,:,i,j,k)=dXCL_Ngeo(1,:,i,j,k) + DCL_Ngeo(i,ll)*XCL_Ngeo(:,ll,j,k)
      dXCL_Ngeo(2,:,i,j,k)=dXCL_Ngeo(2,:,i,j,k) + DCL_Ngeo(j,ll)*XCL_Ngeo(:,i,ll,k)
      dXCL_Ngeo(3,:,i,j,k)=dXCL_Ngeo(3,:,i,j,k) + DCL_Ngeo(k,ll)*XCL_Ngeo(:,i,j,ll)
    END DO !l=0,N
  END DO; END DO; END DO !i,j,k=0,Ngeo

  ! 1.c)Jacobians! grad(X_1) (grad(X_2) x grad(X_3))
  ! Compute Jacobian on NGeo and then interpolate:
  ! required to guarantee conservativity when restarting with N<NGeo
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,1,:,:,:),dX_NgeoRef(:,1,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,2,:,:,:),dX_NgeoRef(:,2,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,3,:,:,:),dX_NgeoRef(:,3,:,:,:))
  DO k=0,NgeoRef; DO j=0,NgeoRef; DO i=0,NgeoRef
    DO dd=1,3
      ee=Cyclic(dd+1) !cyclic!
      ff=Cyclic(dd+2) !cyclic!
      detJac_Ref(1,i,j,k,iElem)=detJac_Ref(1,i,j,k,iElem) + &
                             dX_NgeoRef(dd,1,i,j,k) * &
                           ( dX_NgeoRef(ee,2,i,j,k) * dX_NgeoRef(ff,3,i,j,k) - &
                             dX_NgeoRef(ff,2,i,j,k) * dX_NgeoRef(ee,3,i,j,k)   )
    END DO !nn=1,3
  END DO; END DO; END DO !i,j,k=0,NgeoRef

  !tmp2=0.
  !DO dd=1,3
  !  ee=Cyclic(dd+1) !cyclic!
  !  ff=Cyclic(dd+2) !cyclic!
  !  DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
  !        tmp2(1,i,j,k)=tmp2(1,i,j,k) + &
  !                          dXCL_Ngeo(dd,1,i,j,k) * &
  !                        ( dXCL_Ngeo(ee,2,i,j,k) * dXCL_Ngeo(ff,3,i,j,k) - &
  !                          dXCL_Ngeo(ff,2,i,j,k) * dXCL_Ngeo(ee,3,i,j,k)   )
  !  END DO; END DO; END DO !i,j,k=0,NgeoRef
  !END DO !nn=1,3
  !CALL ChangeBasis3D(1,Ngeo,PP_N,Vdm_CLNgeo_CLN,tmp2,XCL_N_quad)
  !CALL ChangeBasis3D(1,PP_N,PP_N,Vdm_CLN_N,XCL_N_quad,DetJac_N)

  !IF(interpolateFromTree)THEN
  !  !interpolate detJac to the GaussPoints
  !  DO i=0,NgeoRef
  !    dxi=0.5*(xiRef(i)+1.)*Length
  !    CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),NgeoRef,xiRef,wBaryRef,Vdm_xi_Ref(  i,:))
  !    CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),NgeoRef,xiRef,wBaryRef,Vdm_eta_Ref( i,:))
  !    CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),NgeoRef,xiRef,wBaryRef,Vdm_zeta_Ref(i,:))
  !  END DO
  !  tmp=DetJac_Ref(:,:,:,:,iElem)
  !  CALL ChangeBasis3D_XYZ(1,NgeoRef,NgeoRef,Vdm_xi_Ref,Vdm_eta_Ref,Vdm_zeta_Ref,&
  !                         tmp,DetJac_Ref(:,:,:,:,iElem))
  !END IF
  ! interpolate detJac_ref to the solution points
  CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N)

  ! assign to global Variable sJ
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    sJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
  END DO; END DO; END DO !i,j,k=0,PP_N

  ! check for negative Jacobians
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IF(detJac_N(1,i,j,k).LE.0.)&
      WRITE(Unit_StdOut,*) 'Negative Jacobian found on Gauss point. Coords:', Elem_xGP(:,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,N
  ! check scaled Jacobians
  scaledJac(2)=MINVAL(detJac_N(1,:,:,:))/MAXVAL(detJac_N(1,:,:,:))
  IF(scaledJac(2).LT.0.01) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac
    CALL abort(__STAMP__,'Scaled Jacobian lower then tolerance!',iElem)
  END IF

  !2.a) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_N(dd,nn,i,j,k))
  ! N>=Ngeo: interpolate from dXCL_Ngeo (default)
  ! N< Ngeo: directly derive XCL_N
  IF(PP_N.GE.NGeo)THEN !compute first derivative on Ngeo and then interpolate
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,1,:,:,:),dXCL_N(:,1,:,:,:))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,2,:,:,:),dXCL_N(:,2,:,:,:))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:),dXCL_N(:,3,:,:,:))
  ELSE  !N<Ngeo: first interpolate and then compute derivative (important if curved&periodic)
    dXCL_N=0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! Matrix-vector multiplication
      DO ll=0,PP_N
        dXCL_N(1,:,i,j,k)=dXCL_N(1,:,i,j,k) + DCL_N(i,ll)*XCL_N(:,ll,j,k)
        dXCL_N(2,:,i,j,k)=dXCL_N(2,:,i,j,k) + DCL_N(j,ll)*XCL_N(:,i,ll,k)
        dXCL_N(3,:,i,j,k)=dXCL_N(3,:,i,j,k) + DCL_N(k,ll)*XCL_N(:,i,j,ll)
      END DO !l=0,N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=Ngeo

  JaCL_N=0.
  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO nn=1,3
      mm=Cyclic(nn+1) !cyclic!
      ll=Cyclic(nn+2) !cyclic!
      DO dd=1,3
        ee=Cyclic(dd+1) !cyclic!
        ff=Cyclic(dd+2) !cyclic!
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
          ! exact (cross-product) form
          JaCL_N(dd,nn,i,j,k)=dXCL_N(ee,mm,i,j,k)*dXCL_N(ff,ll,i,j,k) - dXCL_N(ee,ll,i,j,k)*dXCL_N(ff,mm,i,j,k)
        END DO; END DO; END DO ! i,j,k
      END DO !dd=1,3
    END DO !nn=1,3
  ELSE ! curl metrics
    ! 2. b.) Calculate X_l grad (X_m)
    R_CL_N=0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      DO nn=1,3
        mm=Cyclic(nn+1) !cyclic!
        ll=Cyclic(nn+2) !cyclic!
        DO dd=1,3
          ! conservative curl form
          !R_CL_N(dd,nn,i,j,k)=XCL_N(ll,i,j,k)*dXCL_N(dd,mm,i,j,k)
          ! invariant curl form
          R_CL_N(dd,nn,i,j,k)=0.5*(XCL_N(ll,i,j,k)*dXCL_N(dd,mm,i,j,k) - &
                                   XCL_N(mm,i,j,k)*dXCL_N(dd,ll,i,j,k)   )
        END DO !dd=1,3
      END DO !nn=1,3
    END DO; END DO; END DO ! i,j,k

    DO nn=1,3; DO dd=1,3
      ee=Cyclic(dd+1) !cyclic!
      ff=Cyclic(dd+2) !cyclic!
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        CycIJK(1)=i; CycIJK(2)=j; CycIJK(3)=k
        ! the f***ing curl is needed for the metric terms
        !
        !Metrics_fTilde(nn)=d/dxi_m (R_CL(1))_l - d/dxi_l (R_CL(1))_m
        !Metrics_gTilde(nn)=d/dxi_m (R_CL(2))_l - d/dxi_l (R_CL(2))_m
        !Metrics_hTilde(nn)=d/dxi_m (R_CL(3))_l - d/dxi_l (R_CL(3))_m
        !
        ! first part of the curl with cyclic indices
        !
        Cyc=CycIJK
        DO q=0,PP_N
          Cyc(ee)=q    !d/dxi_m
          JaCL_N(dd,nn,i,j,k)=JaCL_N(dd,nn,i,j,k) - &
                              DCL_N(CycIJK(ee),q)*R_CL_N(ff,nn,Cyc(1),Cyc(2),Cyc(3))
        END DO!q=0,PP_N
      END DO; END DO; END DO ! i,j,k
    END DO; END DO ! nn,dd=1,3

    DO nn=1,3; DO dd=1,3
      ee=Cyclic(dd+1) !cyclic!
      ff=Cyclic(dd+2) !cyclic!
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        CycIJK(1)=i; CycIJK(2)=j; CycIJK(3)=k
        !
        ! second part of the curl with cyclic indices
        !
        Cyc=CycIJK
        DO q=0,PP_N
          Cyc(ff)=q    !d/dxi_l
          JaCL_N(dd,nn,i,j,k)=JaCL_N(dd,nn,i,j,k) + &
                              DCL_N(CycIJK(ff),q)*R_CL_N(ee,nn,Cyc(1),Cyc(2),Cyc(3))
        END DO!q=0,PP_N
      END DO; END DO; END DO ! i,j,k=0,PP_N
    END DO;END DO ! nn,dd=1,3
  END IF !crossProductMetrics


  !IF(interpolateFromTree)THEN
  !  ! interpolate Metrics from Cheb-Lobatto N on tree level onto GaussPoints N on quad level
  !  DO i=0,PP_N
  !    dxi=0.5*(xGP(i)+1.)*length
  !    CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),PP_N,xiCL_N,wBaryCL_N,Vdm_xi_N(  i,:))
  !    CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),PP_N,xiCL_N,wBaryCL_N,Vdm_eta_N( i,:))
  !    CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),PP_N,xiCL_N,wBaryCL_N,Vdm_zeta_N(i,:))
  !  END DO
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,XCL_N            ,Elem_xGP(      :,:,:,:,iElem))
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem))
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem))
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem))
  !  ! for the metrics and the jacobian, we have to take into account the level !!!!!
  !  Metrics_fTilde(:,:,:,:,iElem)=(length(1)/2.)**2*Metrics_fTilde(:,:,:,:,iElem)
  !  Metrics_gTilde(:,:,:,:,iElem)=(length(2)/2.)**2*Metrics_gTilde(:,:,:,:,iElem)
  !  Metrics_hTilde(:,:,:,:,iElem)=(length(3)/2.)**2*Metrics_hTilde(:,:,:,:,iElem)
  !  sJ(:,:,:,iElem)=(8./PRODUCT(length))*sJ(:,:,:,iElem) ! scale down sJ

  !  ! interpolate Metrics and grid to Cheb-Lobatto on quadrant level for Surface metrics
  !  DO i=0,PP_N
  !    dxi=0.5*(xiCL_N(i)+1.)*length
  !    CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),PP_N,xiCL_N,wBaryCL_N,Vdm_xi_N(  i,:))
  !    CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),PP_N,xiCL_N,wBaryCL_N,Vdm_eta_N( i,:))
  !    CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),PP_N,xiCL_N,wBaryCL_N,Vdm_zeta_N(i,:))
  !  END DO
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,XCL_N            ,XCL_N_quad            )
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(1,:,:,:,:),JaCL_N_quad(1,:,:,:,:))
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(2,:,:,:,:),JaCL_N_quad(2,:,:,:,:))
  !  CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(3,:,:,:,:),JaCL_N_quad(3,:,:,:,:))
  !  !TODO: scale Ja for anisotropic
  !  JaCL_N_quad(:,1,:,:,:)=(length(2)*length(3)/4.)*JaCL_N_quad(:,1,:,:,:)
  !  JaCL_N_quad(:,2,:,:,:)=(length(1)*length(3)/4.)*JaCL_N_quad(:,2,:,:,:)
  !  JaCL_N_quad(:,3,:,:,:)=(length(1)*length(2)/4.)*JaCL_N_quad(:,3,:,:,:)
  !  CALL CalcSurfMetrics(PP_N,JaCL_N_quad,XCL_N_quad,Vdm_CLN_N,iElem,&
  !                       NormVec,TangVec1,TangVec2,SurfElem,Face_xGP)
  !ELSE
    ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,XCL_N            ,Elem_xGP(      :,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem))
    CALL CalcSurfMetrics(PP_N,JaCL_N,XCL_N,Vdm_CLN_N,iElem,&
                         NormVec,TangVec1,TangVec2,SurfElem,Face_xGP)
  !END IF

  IF(PRESENT(XCL_Ngeo_Out))   XCL_Ngeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)= XCL_Ngeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo)
  IF(PRESENT(dXCL_Ngeo_Out)) dXCL_Ngeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)=dXCL_Ngeo(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo)
#ifdef PARTICLES
  StartT2=BOLTZPLATZTIME(MPI_COMM_WORLD)
  CALL GetBezierControlPoints3D(XCL_NGeo(:,:,:,:),iElem)
  endT=BOLTZPLATZTIME(MPI_COMM_WORLD)
  BezierTime=BezierTime+endT-StartT2
#endif /*PARTICLES*/
END DO !iElem=1,nElems

#ifdef PARTICLES
SWRITE(UNIT_stdOut,'(A)') ' '
SWRITE(UNIT_stdOut,'(A)') ' VALIDATION OF BEZIERCONTROLPOINTS ...'
lowerLimit=SideID_minus_upper
!IF(DoRefMapping) lowerLimit=nBCSides
DO iSide=1,lowerLimit
  IF(SUM(ABS(BezierControlPoints3D(:,:,:,iSide))).LT.1e-10)THEN
    IPWRITE(UNIT_stdOut,'(I6,A,I6)') ' Warning, BezierControlPoint is zero! SideID:', iSide
    IPWRITE(UNIT_stdOut,*) 'Points',BezierControlPoints3D(:,:,:,iSide)
  END IF
END DO 

SWRITE(UNIT_stdOut,'(A)') ' '
endt=BOLTZPLATZTIME(MPI_COMM_WORLD)
SWRITE(UNIT_stdOut,'(A,F8.3,A)',ADVANCE='YES')' Calculation of Bezier control points took [',BezierTime - EndT+StartT,'s]'
#else
endt=BOLTZPLATZTIME(MPI_COMM_WORLD)
SWRITE(UNIT_stdOut,'(A,F8.3,A)',ADVANCE='YES')' Calculation of metrics took               [',EndT-StartT,'s]'
#endif /*PARTICLES*/



END SUBROUTINE CalcMetrics 


SUBROUTINE CalcSurfMetrics(Nloc,JaCL_N,XCL_N,Vdm_CLN_N,iElem,NormVec,TangVec1,TangVec2,SurfElem,Face_xGP)
!===================================================================================================================================
! Compute normal and tangential vectors from element metrics. Input is JaCL_N, the 3D element metrics on Cebychev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,     ONLY:CROSS
USE MOD_Mesh_Vars,   ONLY:ElemToSide,nBCSides,nBCSides,nSides!,MortarType
USE MOD_Mesh_Vars,   ONLY:sideid_minus_lower,sideid_minus_upper,nBCSides
USE MOD_Mappings,    ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis, ONLY:ChangeBasis2D
!USE MOD_Mortar_Geo,  ONLY:MortarGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                !< (IN) polynomial degree
INTEGER,INTENT(IN) :: iElem                               !< (IN) element index
REAL,INTENT(IN)    :: JaCL_N(  3,3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) volume metrics of element
REAL,INTENT(IN)    :: XCL_N(     3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) element geo. interpolation points (CL)
REAL,INTENT(IN)    :: Vdm_CLN_N(   0:Nloc,0:Nloc)         !< (IN) Vandermonde matrix from Cheby-Lob on N to final nodeset on N
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   ::    NormVec(3,0:Nloc,0:Nloc,SideID_minus_Lower:SideID_minus_Upper) !< (OUT) element face normal vectors
REAL,INTENT(OUT)   ::   TangVec1(3,0:Nloc,0:Nloc,SideID_minus_Lower:SideID_minus_Upper) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   TangVec2(3,0:Nloc,0:Nloc,SideID_minus_Lower:SideID_minus_Upper) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   SurfElem(  0:Nloc,0:Nloc,SideID_minus_Lower:SideID_minus_Upper) !< (OUT) element face surface area
REAL,INTENT(OUT)   :: Face_xGP(3,0:Nloc,0:Nloc,1:nSides)                       !< (OUT) element face interpolation points
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID
INTEGER            :: NormalDirs(6),TangDirs(6)
INTEGER            :: NormalDir    ,TangDir
REAL               :: NormalSign,   NormalSigns(6)
REAL               :: Ja_Face(3,3,0:Nloc,0:Nloc)
REAL               :: tmp(      3,0:Nloc,0:Nloc)
REAL               :: tmp2(     3,0:Nloc,0:Nloc)
!===================================================================================================================================

NormalDirs(XI_MINUS)  =1; TangDirs(XI_MINUS)    =2; NormalSigns(XI_MINUS)  =-1.
NormalDirs(XI_PLUS)   =1; TangDirs(XI_PLUS)     =2; NormalSigns(XI_PLUS)   = 1.
NormalDirs(ETA_MINUS) =2; TangDirs(ETA_MINUS)   =3; NormalSigns(ETA_MINUS) =-1.
NormalDirs(ETA_PLUS)  =2; TangDirs(ETA_PLUS)    =3; NormalSigns(ETA_PLUS)  = 1.
NormalDirs(ZETA_MINUS)=3; TangDirs(ZETA_MINUS)  =1; NormalSigns(ZETA_MINUS)=-1.
NormalDirs(ZETA_PLUS) =3; TangDirs(ZETA_PLUS)   =1; NormalSigns(ZETA_PLUS) = 1.

DO iLocSide=1,6
  IF(ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) CYCLE ! only master sides with flip=0
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)

  SELECT CASE(iLocSide)
  CASE(XI_MINUS)
    tmp=XCL_N(1:3,0   ,:   ,:   )
  CASE(XI_PLUS)
    tmp=XCL_N(1:3,Nloc,:   ,:   )
  CASE(ETA_MINUS)
    tmp=XCL_N(1:3,:   ,0   ,:   )
  CASE(ETA_PLUS)
    tmp=XCL_N(1:3,:   ,Nloc,:   )
  CASE(ZETA_MINUS)
    tmp=XCL_N(1:3,:   ,:   ,0   )
  CASE(ZETA_PLUS)
    tmp=XCL_N(1:3,:   ,:   ,Nloc)
  END SELECT
  CALL ChangeBasis2D(3,Nloc,Nloc,Vdm_CLN_N,tmp,tmp2)
  ! turn into right hand system of side
  DO q=0,Nloc; DO p=0,Nloc
    pq=CGNS_SideToVol2(p,q,iLocSide)
    ! Compute Face_xGP for sides
    Face_xGP(1:3,p,q,sideID)=tmp2(:,pq(1),pq(2))
  END DO; END DO ! p,q

  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide);
  DO dd=1,3
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=JaCL_N(dd,1:3,0   ,:   ,:   )
    CASE(XI_PLUS)
      tmp=JaCL_N(dd,1:3,Nloc,:   ,:   )
    CASE(ETA_MINUS)
      tmp=JaCL_N(dd,1:3,:   ,0   ,:   )
    CASE(ETA_PLUS)
      tmp=JaCL_N(dd,1:3,:   ,Nloc,:   )
    CASE(ZETA_MINUS)
      tmp=JaCL_N(dd,1:3,:   ,:   ,0   )
    CASE(ZETA_PLUS)
      tmp=JaCL_N(dd,1:3,:   ,:   ,Nloc)
    END SELECT
    CALL ChangeBasis2D(3,Nloc,Nloc,Vdm_CLN_N,tmp,tmp2)
    ! turn into right hand system of side
    DO q=0,Nloc; DO p=0,Nloc
      pq=CGNS_SideToVol2(p,q,iLocSide)
      Ja_Face(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd
  DO q=0,Nloc; DO p=0,Nloc
    SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(NormalDir,:,p,q)**2))
    NormVec( :,p,q,SideID) = NormalSign*Ja_Face(NormalDir,:,p,q)/SurfElem(p,q,SideID)
    TangVec1(:,p,q,SideID) = Ja_Face(TangDir,:,p,q) &
        -SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
    TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
    TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
  END DO; END DO ! p,q
!  IF(MortarType(SideID).GT.0) CALL mortarGeo(SideID,iLocSide,Ja_Face) !set Normal,tangential and surfelem on mortars
END DO

END SUBROUTINE CalcSurfMetrics




END MODULE MOD_Metrics
