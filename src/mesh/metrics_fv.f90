!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Metrics_FV
!===================================================================================================================================
!> This module contains routines for computing the geometries volume and surface metric terms for the finite volume solver,
!> using a temporal polynomial degree of 1 (PP_1), potentially different from the one used for HDG (PP_N).
!> The PP_1 metrics are then averaged into _FV metrics corresponding to a polynomial degree of 0.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_FV
PUBLIC::CalcMetrics_PP_1
PUBLIC::CalcSurfMetrics_PP_1
!==================================================================================================================================

CONTAINS

SUBROUTINE CalcMetrics_PP_1(XCL_NGeo_Out,dXCL_NGeo_out)
!===================================================================================================================================
!> This routine computes the geometries volume metric terms.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,NGeoRef,nGlobalElems,xyzMinMax,GetMeshMinMaxBoundariesIsDone,sJ,crossProductMetrics
USE MOD_Mesh_Vars_FV       ,ONLY: Metrics_fTilde_FV,Metrics_gTilde_FV,Metrics_hTilde_FV,JaCL_N_PP_1
USE MOD_Mesh_Vars_FV       ,ONLY: dXCL_N_PP_1
USE MOD_Mesh_Vars          ,ONLY: DetJac_Ref
USE MOD_Mesh_Vars          ,ONLY: crossProductMetrics
USE MOD_Mesh_Vars          ,ONLY: NodeCoords,Elem_xGP
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem
USE MOD_Mesh_Vars_FV       ,ONLY: XCL_N_PP_1,Vdm_CLN_N_PP_1
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_ChangeBasis        ,ONLY: changeBasis3D,ChangeBasis3D_XYZ
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType
USE MOD_ReadInTools        ,ONLY: GETLOGICAL
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT),OPTIONAL  :: XCL_Ngeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)      ! mapping X(xi) P\in Ngeo
REAL ,INTENT(INOUT),OPTIONAL :: dXCL_Ngeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)   ! jacobi matrix on CL Ngeo
!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,q,iElem
INTEGER :: ll
! Jacobian on CL N and NGeoRef
REAL    :: DetJac_N_PP_1( 1,0:PP_1,   0:PP_1,   0:PP_1)
! interpolation points and derivatives on CL N
! REAL    :: XCL_N(      3,  0:PP_1,0:PP_1,0:PP_1)          ! mapping X(xi) P\in N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: dXCL_Ngeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on CL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

REAL    :: R_CL_N(     3,3,0:PP_1,0:PP_1,0:PP_1)    ! buffer for metric terms, uses XCL_N_PP_1,dXCL_N_PP_1
REAL    :: scaledJac(2)

! Polynomial derivativion matrices
REAL    :: DCL_NGeo(0:Ngeo,0:Ngeo)
REAL    :: DCL_N(   0:PP_1,0:PP_1)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNGeo_CLNgeo( 0:Ngeo   ,0:Ngeo)
REAL    :: Vdm_CLNGeo_NgeoRef(0:NgeoRef,0:Ngeo)
REAL    :: Vdm_NgeoRef_N(     0:PP_1   ,0:NgeoRef)
REAL    :: Vdm_CLNGeo_CLN(    0:PP_1   ,0:Ngeo)
! REAL    :: Vdm_CLN_N(         0:PP_1   ,0:PP_1)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)
REAL    :: xiCL_N(0:PP_1)   ,wBaryCL_N(0:PP_1)
REAL    :: xiCL_NGeo(0:NGeo)   ,wBaryCL_NGeo(0:NGeo)

REAL               :: StartT,EndT
LOGICAL            :: meshCheckRef
REAL,ALLOCATABLE   :: scaledJacRef(:,:,:)
REAL               :: SmallestscaledJacRef
REAL,PARAMETER     :: scaledJacRefTol=0.01

! REAL :: Face_xGP_PP_1      (3,0:PP_1,0:PP_1,1:nSides)
! REAL :: NormVec_PP_1       (3,0:PP_1,0:PP_1,1:nSides)
! REAL :: TangVec1_PP_1      (3,0:PP_1,0:PP_1,1:nSides)
! REAL :: TangVec2_PP_1      (3,0:PP_1,0:PP_1,1:nSides)
! REAL :: SurfElem_PP_1      (0:PP_1,0:PP_1,1:nSides)
! REAL :: Ja_Face_PP_1       (3,3,0:PP_1,0:PP_1,1:nSides)
! REAL :: dXCL_N_PP_1        (3,3,0:PP_1,0:PP_1,0:PP_1,nElems)
REAL :: Metrics_fTilde_PP_1(3,0:PP_1,0:PP_1,0:PP_1,nElems)
REAL :: Metrics_gTilde_PP_1(3,0:PP_1,0:PP_1,0:PP_1,nElems)
REAL :: Metrics_hTilde_PP_1(3,0:PP_1,0:PP_1,0:PP_1,nElems)
!===================================================================================================================================

StartT=PICLASTIME()

! Prerequisites
Metrics_fTilde_FV=0.
Metrics_gTilde_FV=0.
Metrics_hTilde_FV=0.
Metrics_fTilde_PP_1 = 0.
Metrics_gTilde_PP_1 = 0.
Metrics_hTilde_PP_1 = 0.
! Face_xGP_PP_1       = 0.
! NormVec_PP_1        = 0.
! TangVec1_PP_1       = 0.
! TangVec2_PP_1       = 0.
! SurfElem_PP_1       = 0.
! Ja_Face_PP_1        = 0.

! Check Jacobians in Ref already (no good if we only go on because N doesn't catch misbehaving points)
meshCheckRef=GETLOGICAL('meshCheckRef-FV','.TRUE.')

ALLOCATE(    DetJac_Ref(1,0:NGeoRef,0:NGeoRef,0:NGeoRef,nElems))

! Initialize min/max coordinates on faces (but not during load balance restart)
IF(.NOT.GetMeshMinMaxBoundariesIsDone)THEN
  xyzMinMax(1) = HUGE(1.)
  xyzMinMax(2) = -HUGE(1.)
  xyzMinMax(3) = HUGE(1.)
  xyzMinMax(4) = -HUGE(1.)
  xyzMinMax(5) = HUGE(1.)
  xyzMinMax(6) = -HUGE(1.)
END IF

! Initialize Vandermonde and D matrices
! Only use modal Vandermonde for terms that need to be conserved as Jacobian if N_out>PP_N
! Always use interpolation for the rest!

! 1.a) NodeCoords: EQUI Ngeo to CLNgeo and CLN
CALL GetVandermonde(    Ngeo   , NodeTypeVISU, Ngeo    , NodeTypeCL, Vdm_EQNGeo_CLNgeo , modal=.FALSE.)

! 1.b) dXCL_Ngeo:
CALL GetDerivativeMatrix(Ngeo  , NodeTypeCL  , DCL_Ngeo)

! 1.c) Jacobian: CLNgeo to NgeoRef, CLNgeoRef to N
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNgeo_NgeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_1    , NodeType  , Vdm_NgeoRef_N     , modal=.TRUE.)
CALL GetNodesAndWeights(NgeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)

! 1.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
CALL GetVandermonde(    NGeo   , NodeTypeCL  , PP_1    , NodeTypeCL, Vdm_CLNGeo_CLN, modal=.FALSE.)
CALL GetDerivativeMatrix(PP_1  , NodeTypeCL  , DCL_N)

! 2.d) derivatives (dXCL) by projection or by direct derivation (D_CL):

CALL GetNodesAndWeights(NGeo   , NodeTypeCL  , XiCL_NGeo  , wIPBary=wBaryCL_NGeo)

CALL GetVandermonde(    PP_1   , NodeTypeCL  , PP_1    , NodeType,   Vdm_CLN_N_PP_1    , modal=.FALSE.)
Vdm_CLN_N_PP_1(:,:)=0.5
CALL GetNodesAndWeights(PP_1   , NodeTypeCL  , xiCL_N  , wIPBary=wBaryCL_N)


! Outer loop over all elements
DetJac_Ref=0.
dXCL_N_PP_1=0.
SmallestscaledJacRef=HUGE(1.)
DO iElem=1,nElems
  !1.a) Transform from EQUI_Ngeo to CL points on Ngeo and N
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem) ,XCL_Ngeo)
  CALL ChangeBasis3D(3,NGeo,PP_1,Vdm_CLNGeo_CLN,   XCL_Ngeo                  ,XCL_N_PP_1(:,:,:,:,iElem))


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
  ! required to guarantee conservation when restarting with N<NGeo
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,1,:,:,:),dX_NgeoRef(:,1,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,2,:,:,:),dX_NgeoRef(:,2,:,:,:))
  CALL ChangeBasis3D(3,Ngeo,NgeoRef,Vdm_CLNGeo_NgeoRef,dXCL_NGeo(:,3,:,:,:),dX_NgeoRef(:,3,:,:,:))
  DO k=0,NgeoRef; DO j=0,NgeoRef; DO i=0,NgeoRef
    DetJac_Ref(1,i,j,k,iElem)=DetJac_Ref(1,i,j,k,iElem) &
      + dX_NgeoRef(1,1,i,j,k)*(dX_NgeoRef(2,2,i,j,k)*dX_NgeoRef(3,3,i,j,k) - dX_NgeoRef(3,2,i,j,k)*dX_NgeoRef(2,3,i,j,k))  &
      + dX_NgeoRef(2,1,i,j,k)*(dX_NgeoRef(3,2,i,j,k)*dX_NgeoRef(1,3,i,j,k) - dX_NgeoRef(1,2,i,j,k)*dX_NgeoRef(3,3,i,j,k))  &
      + dX_NgeoRef(3,1,i,j,k)*(dX_NgeoRef(1,2,i,j,k)*dX_NgeoRef(2,3,i,j,k) - dX_NgeoRef(2,2,i,j,k)*dX_NgeoRef(1,3,i,j,k))
  END DO; END DO; END DO !i,j,k=0,NgeoRef

  ! Check Jacobians in Ref already (no good if we only go on because N doesn't catch misbehaving points)
  IF (meshCheckRef) THEN
    ALLOCATE(scaledJacRef(0:NGeoRef,0:NGeoRef,0:NGeoRef))
    DO k=0,NGeoRef; DO j=0,NGeoRef; DO i=0,NGeoRef
      scaledJacRef(i,j,k)=DetJac_Ref(1,i,j,k,iElem)/MAXVAL(DetJac_Ref(1,:,:,:,iElem))
      SmallestscaledJacRef=MIN(SmallestscaledJacRef,scaledJacRef(i,j,k))
      IF(scaledJacRef(i,j,k).LT.scaledJacRefTol) THEN
        WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJacRef(i,j,k)
        WRITE(Unit_StdOut,*) 'Coords near:', Elem_xGP(:,INT(PP_N/2),INT(PP_N/2),INT(PP_N/2),iElem)
        WRITE(Unit_StdOut,*) 'This check is optional. You can disable it by setting meshCheckRef = F'
        CALL abort(__STAMP__,&
          'Scaled Jacobian in reference system lower then tolerance in global element:',iElem+offsetElem)
      END IF
    END DO; END DO; END DO !i,j,k=0,N
    DEALLOCATE(scaledJacRef)
  END IF

  ! project detJac_ref onto the solution basis
  CALL ChangeBasis3D(1,NgeoRef,PP_1,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N_PP_1)

  ! assign to global Variable sJ
  sJ(0,0,0,iElem)=1./SUM(SUM(SUM(DetJac_N_PP_1(1,:,:,:),3),2),1)

  ! check for negative Jacobians
  IF(DetJac_N_PP_1(1,0,0,0).LE.0.)&
    WRITE(Unit_StdOut,*) 'Negative Jacobian found on Gauss point. Coords:', Elem_xGP(:,0,0,0,iElem)
  ! check scaled Jacobians
  scaledJac(2)=MINVAL(DetJac_N_PP_1(1,:,:,:))/MAXVAL(DetJac_N_PP_1(1,:,:,:))
  IF(scaledJac(2).LT.0.01) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac
    CALL abort(__STAMP__,&
      'Scaled Jacobian lower then tolerance in global element:',iElem+offsetElem)
  END IF

  !2.a) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_N_PP_1(dd,nn,i,j,k))
  ! N>=Ngeo: interpolate from dXCL_Ngeo (default)
  ! N< Ngeo: directly derive XCL_N_PP_1
  IF(PP_1.GE.NGeo)THEN !compute first derivative on Ngeo and then interpolate
    CALL ChangeBasis3D(3,NGeo,PP_1,Vdm_CLNGeo_CLN,dXCL_NGeo(:,1,:,:,:),dXCL_N_PP_1(:,1,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_1,Vdm_CLNGeo_CLN,dXCL_NGeo(:,2,:,:,:),dXCL_N_PP_1(:,2,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_1,Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:),dXCL_N_PP_1(:,3,:,:,:,iElem))
  ELSE  !N<Ngeo: first interpolate and then compute derivative (important if curved&periodic)
    DO k=0,PP_1; DO j=0,PP_1; DO i=0,PP_1
      ! Matrix-vector multiplication
      ASSOCIATE(dXCL => dXCL_N_PP_1(:,:,i,j,k,iElem))
      DO ll=0,PP_1
        dXCL(1,:)=dXCL(1,:) + DCL_N(i,ll)*XCL_N_PP_1(:,ll,j,k,iElem)
        dXCL(2,:)=dXCL(2,:) + DCL_N(j,ll)*XCL_N_PP_1(:,i,ll,k,iElem)
        dXCL(3,:)=dXCL(3,:) + DCL_N(k,ll)*XCL_N_PP_1(:,i,j,ll,iElem)
      END DO !l=0,N
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=Ngeo

  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO k=0,PP_1; DO j=0,PP_1; DO i=0,PP_1
      ASSOCIATE(dXCL => dXCL_N_PP_1(:,:,i,j,k,iElem))
      ! exact (cross-product) form
      ! Ja(:)^nn = ( d/dxi_(nn+1) XCL_N_PP_1(:) ) x (d/xi_(nn+2) XCL_N_PP_1(:))
      !
      ! JaCL_N_PP_1(dd,nn) = dXCL_N_PP_1(dd+1,nn+1)*dXCL_N_PP_1(dd+2,nn+2) -dXCL_N_PP_1(dd+1,nn+2)*dXCL_N_PP_1(dd+2,nn+1)
      JaCL_N_PP_1(1,1,i,j,k,iElem)=dXCL(2,2)*dXCL(3,3) - dXCL(2,3)*dXCL(3,2)
      JaCL_N_PP_1(2,1,i,j,k,iElem)=dXCL(3,2)*dXCL(1,3) - dXCL(3,3)*dXCL(1,2)
      JaCL_N_PP_1(3,1,i,j,k,iElem)=dXCL(1,2)*dXCL(2,3) - dXCL(1,3)*dXCL(2,2)
      JaCL_N_PP_1(1,2,i,j,k,iElem)=dXCL(2,3)*dXCL(3,1) - dXCL(2,1)*dXCL(3,3)
      JaCL_N_PP_1(2,2,i,j,k,iElem)=dXCL(3,3)*dXCL(1,1) - dXCL(3,1)*dXCL(1,3)
      JaCL_N_PP_1(3,2,i,j,k,iElem)=dXCL(1,3)*dXCL(2,1) - dXCL(1,1)*dXCL(2,3)
      JaCL_N_PP_1(1,3,i,j,k,iElem)=dXCL(2,1)*dXCL(3,2) - dXCL(2,2)*dXCL(3,1)
      JaCL_N_PP_1(2,3,i,j,k,iElem)=dXCL(3,1)*dXCL(1,2) - dXCL(3,2)*dXCL(1,1)
      JaCL_N_PP_1(3,3,i,j,k,iElem)=dXCL(1,1)*dXCL(2,2) - dXCL(1,2)*dXCL(2,1)
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  ELSE ! curl metrics
    ! invariant curl form, as cross product: R^dd = 1/2( XCL_N_PP_1(:) x (d/dxi_dd XCL_N_PP_1(:)))
    !
    !R_CL_N(dd,nn)=1/2*( XCL_N_PP_1(nn+2)* d/dxi_dd XCL_N_PP_1(nn+1) - XCL_N_PP_1(nn+1)* d/dxi_dd XCL_N_PP_1(nn+2))
    DO k=0,PP_1; DO j=0,PP_1; DO i=0,PP_1
      ASSOCIATE(dXCL => dXCL_N_PP_1(:,:,i,j,k,iElem))
      R_CL_N(:,1,i,j,k)=0.5*(XCL_N_PP_1(3,i,j,k,iElem)*dXCL(:,2) - XCL_N_PP_1(2,i,j,k,iElem)*dXCL(:,3) )
      R_CL_N(:,2,i,j,k)=0.5*(XCL_N_PP_1(1,i,j,k,iElem)*dXCL(:,3) - XCL_N_PP_1(3,i,j,k,iElem)*dXCL(:,1) )
      R_CL_N(:,3,i,j,k)=0.5*(XCL_N_PP_1(2,i,j,k,iElem)*dXCL(:,1) - XCL_N_PP_1(1,i,j,k,iElem)*dXCL(:,2) )
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
    ! Metrics are the curl of R:  Ja(:)^nn = -(curl R_CL(:,nn))
    ! JaCL_N_PP_1(dd,nn)= -[d/dxi_(dd+1) RCL(dd+2,nn) - d/dxi_(dd+2) RCL(dd+1,nn) ]
    !              =   d/dxi_(dd+2) RCL(dd+1,nn) - d/dxi_(dd+1) RCL(dd+2,nn)
    DO k=0,PP_1; DO j=0,PP_1; DO i=0,PP_1
      ASSOCIATE(JaCL_PP_1 => JaCL_N_PP_1(:,:,i,j,k,iElem))
      DO q=0,PP_1
        JaCL_PP_1(1,:)=JaCL_PP_1(1,:) - DCL_N(j,q)*R_CL_N(3,:,i,q,k)
        JaCL_PP_1(2,:)=JaCL_PP_1(2,:) - DCL_N(k,q)*R_CL_N(1,:,i,j,q)
        JaCL_PP_1(3,:)=JaCL_PP_1(3,:) - DCL_N(i,q)*R_CL_N(2,:,q,j,k)
      END DO!q=0,PP_N
      DO q=0,PP_1
        JaCL_PP_1(1,:)=JaCL_PP_1(1,:) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
        JaCL_PP_1(2,:)=JaCL_PP_1(2,:) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
        JaCL_PP_1(3,:)=JaCL_PP_1(3,:) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
      END DO!q=0,PP_N
      END ASSOCIATE
! same with only one loop, gives different roundoff ...
!      DO q=0,PP_N
!        JaCL_N_PP_1(1,:,i,j,k)=JaCL_N_PP_1(1,:,i,j,k) - DCL_N(j,q)*R_CL_N(3,:,i,q,k) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
!        JaCL_N_PP_1(2,:,i,j,k)=JaCL_N_PP_1(2,:,i,j,k) - DCL_N(k,q)*R_CL_N(1,:,i,j,q) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
!        JaCL_N_PP_1(3,:,i,j,k)=JaCL_N_PP_1(3,:,i,j,k) - DCL_N(i,q)*R_CL_N(2,:,q,j,k) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
!      END DO!q=0,PP_N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !crossProductMetrics

  ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
  CALL ChangeBasis3D(3,PP_1,PP_1,Vdm_CLN_N_PP_1,JaCL_N_PP_1(1,:,:,:,:,iElem),Metrics_fTilde_PP_1(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_1,PP_1,Vdm_CLN_N_PP_1,JaCL_N_PP_1(2,:,:,:,:,iElem),Metrics_gTilde_PP_1(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_1,PP_1,Vdm_CLN_N_PP_1,JaCL_N_PP_1(3,:,:,:,:,iElem),Metrics_hTilde_PP_1(:,:,:,:,iElem))
  ! CALL CalcSurfMetrics_PP_1(PP_1,JaCL_N_PP_1,XCL_N_PP_1,Vdm_CLN_N_PP_1,iElem,&
  !                       NormVec_PP_1,TangVec1_PP_1,TangVec2_PP_1,SurfElem_PP_1,Face_xGP_PP_1,Ja_Face_PP_1)

  ! particle mapping
  IF(PRESENT(XCL_Ngeo_Out))   XCL_Ngeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)= XCL_Ngeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo)
  IF(PRESENT(dXCL_ngeo_out)) dXCL_Ngeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)=dXCL_Ngeo(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo)
END DO !iElem=1,nElems

! PP_1 metrics to global ones
! Face_xGP_FV      (:,0,0,:)   = Face_xGP_PP_1      (:,0,0,:)
! NormVec_FV       (:,0,0,:)   = NormVec_PP_1       (:,0,0,:)
! TangVec1_FV      (:,0,0,:)   = TangVec1_PP_1      (:,0,0,:)
! TangVec2_FV      (:,0,0,:)   = TangVec2_PP_1      (:,0,0,:)
! SurfElem_FV      (0,0,:)     = SUM(SUM(SurfElem_PP_1(:,:,:),2),1)
! Ja_Face_FV       (:,:,0,0,:) = SUM(SUM(Ja_Face_PP_1(:,:,:,:,:),4),3)
Metrics_fTilde_FV(:,0,0,0,:) = SUM(SUM(Metrics_fTilde_PP_1(:,0,:,:,:),3),2)
Metrics_gTilde_FV(:,0,0,0,:) = SUM(SUM(Metrics_gTilde_PP_1(:,:,0,:,:),3),2)
Metrics_hTilde_FV(:,0,0,0,:) = SUM(SUM(Metrics_hTilde_PP_1(:,:,:,0,:),3),2)

! Communicate smallest ref. Jacobian and display
#if USE_MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , SmallestscaledJacRef , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
ELSE
  CALL MPI_REDUCE(SmallestscaledJacRef   , 0          , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
END IF
#endif /*USE_MPI*/
LBWRITE (*,'(A,ES18.10E3,A,I0,A,ES13.5E3)') " Smallest scaled Jacobian in reference system: ",SmallestscaledJacRef,&
    " (",nGlobalElems," global elements). Abort threshold is set to:", scaledJacRefTol

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'Calculation of PP_1 metrics took ', DisplayLine=.FALSE.)

END SUBROUTINE CalcMetrics_PP_1


SUBROUTINE CalcSurfMetrics_PP_1(JaCL_N_PP_1,iElem)
!===================================================================================================================================
! Compute normal and tangential vectors from element metrics. Input is JaCL_N_PP_1, the 3D element metrics on Cebychev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,        ONLY:CROSS
USE MOD_Mesh_Vars,      ONLY:ElemToSide,MortarType,xyzMinMax,GetMeshMinMaxBoundariesIsDone
USE MOD_Mesh_Vars,      ONLY:NormalDirs,TangDirs,NormalSigns
USE MOD_Mesh_Vars_FV,   ONLY:XCL_N_PP_1,Vdm_CLN_N_PP_1
USE MOD_Mesh_Vars_FV,   ONLY:NormVec_PP_1,TangVec1_PP_1,TangVec2_PP_1,SurfElem_PP_1,Face_xGP_PP_1,Ja_Face_PP_1
USE MOD_Mappings,       ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis,    ONLY:ChangeBasis2D
USE MOD_Mortar_Metrics, ONLY:Mortar_CalcSurfMetrics
USE MOD_Metrics,        ONLY:SurfMetricsFromJa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: JaCL_N_PP_1(1:3,1:3,0:PP_1,0:PP_1,0:PP_1) !< (IN) volume metrics of element
INTEGER,INTENT(IN) :: iElem                                !< (IN) element index
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID,SideID2,iMortar,nbSideIDs(4)
INTEGER            :: NormalDir,TangDir
REAL               :: NormalSign
REAL               :: Ja_Face_l_PP_1(3,3,0:PP_1,0:PP_1)
REAL               :: Mortar_Ja_PP_1(3,3,0:PP_1,0:PP_1,4)
REAL               :: Mortar_xGP_PP_1( 3,0:PP_1,0:PP_1,4)
REAL               :: tmp_PP_1(        3,0:PP_1,0:PP_1)
REAL               :: tmp2_PP_1(       3,0:PP_1,0:PP_1)
!==================================================================================================================================

#if PP_dim == 3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif
  IF(ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) CYCLE ! only master sides with flip=0
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)

  SELECT CASE(iLocSide)
  CASE(XI_MINUS)
    tmp_PP_1=XCL_N_PP_1(1:3,0   ,:   ,:   ,iElem)
  CASE(XI_PLUS)
    tmp_PP_1=XCL_N_PP_1(1:3,PP_1,:   ,:   ,iElem)
  CASE(ETA_MINUS)
    tmp_PP_1=XCL_N_PP_1(1:3,:   ,0   ,:   ,iElem)
  CASE(ETA_PLUS)
    tmp_PP_1=XCL_N_PP_1(1:3,:   ,PP_1,:   ,iElem)
  CASE(ZETA_MINUS)
    tmp_PP_1=XCL_N_PP_1(1:3,:   ,:   ,0   ,iElem)
  CASE(ZETA_PLUS)
    tmp_PP_1=XCL_N_PP_1(1:3,:   ,:   ,PP_1,iElem)
  END SELECT
  CALL ChangeBasis2D(3,PP_1,PP_1,Vdm_CLN_N_PP_1,tmp_PP_1,tmp2_PP_1)

  ! Get min/max coordinate from current face
  IF(.NOT.GetMeshMinMaxBoundariesIsDone)THEN
    xyzMinMax(1) = MIN(xyzMinMax(1), MINVAL(tmp_PP_1(1,:,:)))
    xyzMinMax(2) = MAX(xyzMinMax(2), MAXVAL(tmp_PP_1(1,:,:)))

    xyzMinMax(3) = MIN(xyzMinMax(3), MINVAL(tmp_PP_1(2,:,:)))
    xyzMinMax(4) = MAX(xyzMinMax(4), MAXVAL(tmp_PP_1(2,:,:)))

    xyzMinMax(5) = MIN(xyzMinMax(5), MINVAL(tmp_PP_1(3,:,:)))
    xyzMinMax(6) = MAX(xyzMinMax(6), MAXVAL(tmp_PP_1(3,:,:)))
  END IF ! .NOT.GetMeshMinMaxBoundariesIsDone

  ! turn into right hand system of side
  DO q=0,PP_1; DO p=0,PP_1
    pq=CGNS_SideToVol2(PP_1,p,q,iLocSide)
    ! Compute Face_xGP_PP_1 for sides
    Face_xGP_PP_1(1:3,p,q,SideID)=tmp2_PP_1(:,pq(1),pq(2))
  END DO; END DO ! p,q

  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide);
  DO dd=1,3
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,0   ,:   ,:   )
    CASE(XI_PLUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,PP_1,:   ,:   )
    CASE(ETA_MINUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,:   ,0   ,:   )
    CASE(ETA_PLUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,:   ,PP_1,:   )
    CASE(ZETA_MINUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,:   ,:   ,0   )
    CASE(ZETA_PLUS)
      tmp_PP_1=JaCL_N_PP_1(dd,1:3,:   ,:   ,PP_1)
    END SELECT
    CALL ChangeBasis2D(3,PP_1,PP_1,Vdm_CLN_N_PP_1,tmp_PP_1,tmp2_PP_1)
    ! turn into right hand system of side
    DO q=0,PP_1; DO p=0,PP_1
      pq=CGNS_SideToVol2(PP_1,p,q,iLocSide)
      Ja_Face_l_PP_1(dd,1:3,p,q)=tmp2_PP_1(:,pq(1),pq(2))
    ! DEBUG old version
      !Ja_Face_PP_1(dd,1:3,p,q)=tmp2_PP_1(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd
  Ja_Face_PP_1(:,:,:,:,SideID)=Ja_Face_l_PP_1


  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)
  CALL SurfMetricsFromJa(PP_1,NormalDir,TangDir,NormalSign,Ja_Face_l_PP_1,&
                         NormVec_PP_1(:,:,:,SideID),TangVec1_PP_1(:,:,:,SideID),&
                         TangVec2_PP_1(:,:,:,SideID),SurfElem_PP_1(:,:,SideID))

  !compute metrics for mortar faces, interpolate Ja_Face_PP_1 to small sides
  IF(MortarType(1,SideID).GT.0)THEN
    CALL Mortar_CalcSurfMetrics(SideID,PP_1,Ja_Face_l_PP_1,Face_xGP_PP_1(:,:,:,SideID),&
                                            Mortar_Ja_PP_1,Mortar_xGP_PP_1,nbSideIDs)
    DO iMortar=1,4
      SideID2=nbSideIDs(iMortar)
      IF(SideID2.LT.1) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
      Ja_Face_PP_1(:,:,:,:,SideID2)=Mortar_Ja_PP_1(:,:,:,:,iMortar)
      Face_xGP_PP_1(:,:,:,SideID2) = Mortar_xGP_PP_1(:,:,:,iMortar)
      CALL SurfMetricsFromJa(PP_1,NormalDir,TangDir,NormalSign,Mortar_Ja_PP_1(:,:,:,:,iMortar),&
                             NormVec_PP_1(:,:,:,SideID2),TangVec1_PP_1(:,:,:,SideID2),&
                             TangVec2_PP_1(:,:,:,SideID2),SurfElem_PP_1(:,:,SideID2))
    END DO

  END IF
END DO

END SUBROUTINE CalcSurfMetrics_PP_1

#endif /*USE_FV*/
END MODULE MOD_Metrics_FV