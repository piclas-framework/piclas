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

MODULE MOD_Metrics
!===================================================================================================================================
!> \brief This module contains routines for computing the geometries volume and surface metric terms.
!>
!> Compute the volume and surface metric terms:
!>     Metrics_fTilde(n=1:3,i,j,k,iElem)=Ja_n^1
!>     Metrics_gTilde(n=1:3,i,j,k,iElem)=Ja_n^2
!>     Metrics_hTilde(n=1:3,i,j,k,iElem)=Ja_n^3
!>
!>   Per Element we do:
!>   1.) a.) Preparation: the geometry (equidistant nodal basis, NGeo+1 points/dir) is interpolated to a high precision
!>           mapping X_n(xi_i) using a Chebyshev-Lobatto basis and stored in XCL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!>       b.) Computing the gradients: compute the derivative of the mapping XCL_NGeo in \f$ (xi_1,xi_2,xi_3) \f$ direction,
!>           using a polynomial derivative Matrix at degree NGeo.
!>       c.) Computing the Jacobian: compute Jacobian JRef at a degree of NGeoRef=3*NGeo (exact).
!>                                   For this gradients have to be interpolated to NGeoRef first.
!>                                   Then project JRef down to degree N. Finally check for negative Jacobians.
!>       d.) For computing Ja the gradients at degree N are required: if N>=NGeo directly interpolate dXCL_NGeo to dXCL_N,
!>                                                                    else compute dXCL_N from XCL_N directly.
!>
!>   2.) for each direction n
!>       a.) compute the nth vector and for each Chebyshev point (:,i,j,k)
!>          \f$(dXCL_n^1,dXCL_n^2,dXCL_n^3)^T=(X_l grad_xi (X_m) )\f$ for n=1,2,3 and (n,m,l) cyclic
!>       b.) interpolate the dXCL_n vector defined primarily on (NGeo+1)x(NGeo+1)x(Ngeo+1) Chebyshev-Lobatto points to
!>             (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points and write to Ja_n(1:3,i,j,k) i,j,k=[0:N]
!>       c.) compute the curl of vector Ja_n(1:3,i,j,k) using the derivative Matrix DCL_N [NxN]
!>       d.) interpolate from (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points to  Gauss-Points (N+1)x(N+1)x(N+1) (exact!)
!>       e.) store Ja_n in the Metrics arrays
!>
!>   3.) Compute the surface metrics (normal/tangential vectors, surface area) from volume metrics for each side.
!>
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
INTERFACE BuildCoords
  MODULE PROCEDURE BuildCoords
END INTERFACE

INTERFACE CalcMetrics
  MODULE PROCEDURE CalcMetrics
END INTERFACE

INTERFACE CalcSurfMetrics
  MODULE PROCEDURE CalcSurfMetrics
END INTERFACE

INTERFACE SurfMetricsFromJa
  MODULE PROCEDURE SurfMetricsFromJa
END INTERFACE

INTERFACE CalcMetricsErrorDiff
  MODULE PROCEDURE CalcMetricsErrorDiff
END INTERFACE

PUBLIC::BuildCoords
PUBLIC::CalcMetrics
PUBLIC::CalcSurfMetrics
PUBLIC::SurfMetricsFromJa
PUBLIC::CalcMetricsErrorDiff
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine takes the equidistant node coordinats of the mesh (on NGeo+1 points) and uses them to build the coordinates
!> of solution/interpolation points of type NodeType on polynomial degree Nloc (Nloc+1 points per direction).
!==================================================================================================================================
SUBROUTINE BuildCoords(NodeCoords,Nloc,VolumeCoords)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,nElems
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ, ChangeBasis3D
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems)         !< Equidistant mesh coordinates
INTEGER,INTENT(IN)            :: Nloc                                                  !< Convert to Nloc+1 points per direction
REAL,INTENT(OUT)              :: VolumeCoords(3,0:Nloc,0:Nloc,0:Nloc,nElems)       !< OUT: Coordinates of solution/interpolation
                                                                                       !< points
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
REAL                          :: Vdm_EQNGeo_CLNloc(0:Nloc ,0:Ngeo)
REAL                          :: Vdm_CLNloc_Nloc  (0:Nloc ,0:Nloc)
!==================================================================================================================================

CALL GetVandermonde(    NGeo, NodeTypeVISU, NLoc, NodeTypeCL, Vdm_EQNGeo_CLNloc,  modal=.FALSE.)
CALL GetVandermonde(    Nloc, NodeTypeCL  , Nloc, NodeType  , Vdm_CLNloc_Nloc,     modal=.FALSE.)

! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!       Important for curved meshes if NGeo<N, no effect for N>=NGeo

!1.a) Transform from EQUI_NGeo to solution points on Nloc
Vdm_EQNGeo_CLNloc=MATMUL(Vdm_CLNloc_Nloc,Vdm_EQNGeo_CLNloc)
DO iElem=1,nElems
  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),VolumeCoords(:,:,:,:,iElem))
END DO

END SUBROUTINE BuildCoords


SUBROUTINE CalcMetrics(XCL_NGeo_Out,dXCL_NGeo_out)
!===================================================================================================================================
!> This routine computes the geometries volume metric terms.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,NGeoRef,nGlobalElems,xyzMinMax,GetMeshMinMaxBoundariesIsDone
USE MOD_Mesh_Vars          ,ONLY: sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,crossProductMetrics
USE MOD_Mesh_Vars          ,ONLY: Face_xGP,normVec,surfElem,TangVec1,TangVec2
USE MOD_Mesh_Vars          ,ONLY: nElems,dXCL_N
USE MOD_Mesh_Vars          ,ONLY: DetJac_Ref,Ja_Face
USE MOD_Mesh_Vars          ,ONLY: crossProductMetrics
USE MOD_Mesh_Vars          ,ONLY: NodeCoords,Elem_xGP
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem
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
REAL    :: DetJac_N( 1,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL    :: XCL_N(      3,  0:PP_N,0:PP_N,0:PP_N)          ! mapping X(xi) P\in N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: dXCL_Ngeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on CL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

REAL    :: R_CL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! buffer for metric terms, uses XCL_N,dXCL_N
REAL    :: JaCL_N(     3,3,0:PP_N,0:PP_N,0:PP_N)    ! metric terms P\in N
REAL    :: scaledJac(2)

! Polynomial derivativion matrices
REAL    :: DCL_NGeo(0:Ngeo,0:Ngeo)
REAL    :: DCL_N(   0:PP_N,0:PP_N)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNGeo_CLNgeo( 0:Ngeo   ,0:Ngeo)
REAL    :: Vdm_CLNGeo_NgeoRef(0:NgeoRef,0:Ngeo)
REAL    :: Vdm_NgeoRef_N(     0:PP_N   ,0:NgeoRef)
REAL    :: Vdm_CLNGeo_CLN(    0:PP_N   ,0:Ngeo)
REAL    :: Vdm_CLN_N(         0:PP_N   ,0:PP_N)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)
REAL    :: xiCL_N(0:PP_N)   ,wBaryCL_N(0:PP_N)
REAL    :: xiCL_NGeo(0:NGeo)   ,wBaryCL_NGeo(0:NGeo)

REAL               :: StartT,EndT
LOGICAL            :: meshCheckRef
REAL,ALLOCATABLE   :: scaledJacRef(:,:,:)
REAL               :: SmallestscaledJacRef
REAL,PARAMETER     :: scaledJacRefTol=0.01
!===================================================================================================================================

StartT=PICLASTIME()

! Prerequisites
Metrics_fTilde=0.
Metrics_gTilde=0.
Metrics_hTilde=0.

! Check Jacobians in Ref already (no good if we only go on because N doesn't catch misbehaving points)
meshCheckRef=GETLOGICAL('meshCheckRef','.TRUE.')

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
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NgeoRef_N     , modal=.TRUE.)
CALL GetNodesAndWeights(NgeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)

! 1.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , PP_N    , NodeTypeCL, Vdm_CLNgeo_CLN    , modal=.FALSE.)
CALL GetDerivativeMatrix(PP_N  , NodeTypeCL  , DCL_N)

! 2.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
CALL GetVandermonde(    PP_N   , NodeTypeCL  , PP_N    , NodeType,   Vdm_CLN_N         , modal=.FALSE.)
CALL GetNodesAndWeights(PP_N   , NodeTypeCL  , xiCL_N  , wIPBary=wBaryCL_N)

CALL GetNodesAndWeights(NGeo   , NodeTypeCL  , XiCL_NGeo  , wIPBary=wBaryCL_NGeo)

! Outer loop over all elements
DetJac_Ref=0.
dXCL_N=0.
SmallestscaledJacRef=HUGE(1.)
DO iElem=1,nElems
  !1.a) Transform from EQUI_Ngeo to CL points on Ngeo and N
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem)            ,XCL_Ngeo)
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,   XCL_Ngeo                             ,XCL_N)

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
  CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N)

  ! assign to global Variable sJ
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    sJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
  END DO; END DO; END DO !i,j,k=0,PP_N

  ! check for negative Jacobians
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IF(DetJac_N(1,i,j,k).LE.0.)&
      WRITE(Unit_StdOut,*) 'Negative Jacobian found on Gauss point. Coords:', Elem_xGP(:,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,N
  ! check scaled Jacobians
  scaledJac(2)=MINVAL(DetJac_N(1,:,:,:))/MAXVAL(DetJac_N(1,:,:,:))
  IF(scaledJac(2).LT.0.01) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac
    CALL abort(__STAMP__,&
      'Scaled Jacobian lower then tolerance in global element:',iElem+offsetElem)
  END IF

  !2.a) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_N(dd,nn,i,j,k))
  ! N>=Ngeo: interpolate from dXCL_Ngeo (default)
  ! N< Ngeo: directly derive XCL_N
  IF(PP_N.GE.NGeo)THEN !compute first derivative on Ngeo and then interpolate
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,1,:,:,:),dXCL_N(:,1,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,2,:,:,:),dXCL_N(:,2,:,:,:,iElem))
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:),dXCL_N(:,3,:,:,:,iElem))
  ELSE  !N<Ngeo: first interpolate and then compute derivative (important if curved&periodic)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! Matrix-vector multiplication
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      DO ll=0,PP_N
        dXCL(1,:)=dXCL(1,:) + DCL_N(i,ll)*XCL_N(:,ll,j,k)
        dXCL(2,:)=dXCL(2,:) + DCL_N(j,ll)*XCL_N(:,i,ll,k)
        dXCL(3,:)=dXCL(3,:) + DCL_N(k,ll)*XCL_N(:,i,j,ll)
      END DO !l=0,N
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=Ngeo

  JaCL_N=0.
  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      ! exact (cross-product) form
      ! Ja(:)^nn = ( d/dxi_(nn+1) XCL_N(:) ) x (d/xi_(nn+2) XCL_N(:))
      !
      ! JaCL_N(dd,nn) = dXCL_N(dd+1,nn+1)*dXCL_N(dd+2,nn+2) -dXCL_N(dd+1,nn+2)*dXCL_N(dd+2,nn+1)
      JaCL_N(1,1,i,j,k)=dXCL(2,2)*dXCL(3,3) - dXCL(2,3)*dXCL(3,2)
      JaCL_N(2,1,i,j,k)=dXCL(3,2)*dXCL(1,3) - dXCL(3,3)*dXCL(1,2)
      JaCL_N(3,1,i,j,k)=dXCL(1,2)*dXCL(2,3) - dXCL(1,3)*dXCL(2,2)
      JaCL_N(1,2,i,j,k)=dXCL(2,3)*dXCL(3,1) - dXCL(2,1)*dXCL(3,3)
      JaCL_N(2,2,i,j,k)=dXCL(3,3)*dXCL(1,1) - dXCL(3,1)*dXCL(1,3)
      JaCL_N(3,2,i,j,k)=dXCL(1,3)*dXCL(2,1) - dXCL(1,1)*dXCL(2,3)
      JaCL_N(1,3,i,j,k)=dXCL(2,1)*dXCL(3,2) - dXCL(2,2)*dXCL(3,1)
      JaCL_N(2,3,i,j,k)=dXCL(3,1)*dXCL(1,2) - dXCL(3,2)*dXCL(1,1)
      JaCL_N(3,3,i,j,k)=dXCL(1,1)*dXCL(2,2) - dXCL(1,2)*dXCL(2,1)
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  ELSE ! curl metrics
    ! invariant curl form, as cross product: R^dd = 1/2( XCL_N(:) x (d/dxi_dd XCL_N(:)))
    !
    !R_CL_N(dd,nn)=1/2*( XCL_N(nn+2)* d/dxi_dd XCL_N(nn+1) - XCL_N(nn+1)* d/dxi_dd XCL_N(nn+2))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      R_CL_N(:,1,i,j,k)=0.5*(XCL_N(3,i,j,k)*dXCL(:,2) - XCL_N(2,i,j,k)*dXCL(:,3) )
      R_CL_N(:,2,i,j,k)=0.5*(XCL_N(1,i,j,k)*dXCL(:,3) - XCL_N(3,i,j,k)*dXCL(:,1) )
      R_CL_N(:,3,i,j,k)=0.5*(XCL_N(2,i,j,k)*dXCL(:,1) - XCL_N(1,i,j,k)*dXCL(:,2) )
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
    ! Metrics are the curl of R:  Ja(:)^nn = -(curl R_CL(:,nn))
    ! JaCL_N(dd,nn)= -[d/dxi_(dd+1) RCL(dd+2,nn) - d/dxi_(dd+2) RCL(dd+1,nn) ]
    !              =   d/dxi_(dd+2) RCL(dd+1,nn) - d/dxi_(dd+1) RCL(dd+2,nn)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(JaCL => JaCL_N(:,:,i,j,k))
      DO q=0,PP_N
        JaCL(1,:)=JaCL(1,:) - DCL_N(j,q)*R_CL_N(3,:,i,q,k)
        JaCL(2,:)=JaCL(2,:) - DCL_N(k,q)*R_CL_N(1,:,i,j,q)
        JaCL(3,:)=JaCL(3,:) - DCL_N(i,q)*R_CL_N(2,:,q,j,k)
      END DO!q=0,PP_N
      DO q=0,PP_N
        JaCL(1,:)=JaCL(1,:) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
        JaCL(2,:)=JaCL(2,:) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
        JaCL(3,:)=JaCL(3,:) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
      END DO!q=0,PP_N
      END ASSOCIATE
! same with only one loop, gives different roundoff ...
!      DO q=0,PP_N
!        JaCL_N(1,:,i,j,k)=JaCL_N(1,:,i,j,k) - DCL_N(j,q)*R_CL_N(3,:,i,q,k) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
!        JaCL_N(2,:,i,j,k)=JaCL_N(2,:,i,j,k) - DCL_N(k,q)*R_CL_N(1,:,i,j,q) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
!        JaCL_N(3,:,i,j,k)=JaCL_N(3,:,i,j,k) - DCL_N(i,q)*R_CL_N(2,:,q,j,k) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
!      END DO!q=0,PP_N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !crossProductMetrics

  ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem))
  CALL CalcSurfMetrics(PP_N,JaCL_N,XCL_N,Vdm_CLN_N,iElem,&
                        NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face)
#ifdef maxwell
#if defined(ROS) || defined(IMPA)
  CALL CalcElemLocalSurfMetrics(PP_N,JaCL_N,Vdm_CLN_N,iElem)
#endif /*ROS or IMPA*/
#endif /*maxwell*/

  ! particle mapping
  IF(PRESENT(XCL_Ngeo_Out))   XCL_Ngeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)= XCL_Ngeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo)
  IF(PRESENT(dXCL_ngeo_out)) dXCL_Ngeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)=dXCL_Ngeo(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo)
END DO !iElem=1,nElems

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
CALL DisplayMessageAndTime(EndT-StartT, 'Calculation of metrics took!', DisplayLine=.FALSE.)

END SUBROUTINE CalcMetrics


SUBROUTINE CalcSurfMetrics(Nloc,JaCL_N,XCL_N,Vdm_CLN_N,iElem,NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face)
!===================================================================================================================================
! Compute normal and tangential vectors from element metrics. Input is JaCL_N, the 3D element metrics on Cebychev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,     ONLY:CROSS
USE MOD_Mesh_Vars,   ONLY:ElemToSide,nSides,MortarType,xyzMinMax,GetMeshMinMaxBoundariesIsDone
USE MOD_Mesh_Vars,   ONLY:NormalDirs,TangDirs,NormalSigns
USE MOD_Mappings,    ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis, ONLY:ChangeBasis2D
USE MOD_Mortar_Metrics, ONLY:Mortar_CalcSurfMetrics
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                !< (IN) polynomial degree
INTEGER,INTENT(IN) :: iElem                               !< (IN) element index
REAL,INTENT(IN)    :: JaCL_N(1:3,1:3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) volume metrics of element
REAL,INTENT(IN)    :: XCL_N(     1:3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) element geo. interpolation points (CL)
REAL,INTENT(IN)    :: Vdm_CLN_N(   0:Nloc,0:Nloc)         !< (IN) Vandermonde matrix from Cheby-Lob on N to final nodeset on N
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   ::    NormVec(3,0:Nloc,0:Nloc,1:nSides) !< (OUT) element face normal vectors
REAL,INTENT(OUT)   ::   TangVec1(3,0:Nloc,0:Nloc,1:nSides) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   TangVec2(3,0:Nloc,0:Nloc,1:nSides) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   SurfElem(  0:Nloc,0:Nloc,1:nSides) !< (OUT) element face surface area
REAL,INTENT(OUT)   :: Face_xGP(1:3,0:Nloc,0:Nloc,1:nSides)                       !< (OUT) element face interpolation points
REAL,INTENT(OUT),OPTIONAL :: Ja_Face(3,3,0:Nloc,0:Nloc,1:nSides)  !< (OUT) surface metrics
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID,SideID2,iMortar,nbSideIDs(4)
INTEGER            :: NormalDir,TangDir
REAL               :: NormalSign
REAL               :: Ja_Face_l(3,3,0:Nloc,0:Nloc)
REAL               :: Mortar_Ja(3,3,0:Nloc,0:Nloc,4)
REAL               :: Mortar_xGP( 3,0:Nloc,0:Nloc,4)
REAL               :: tmp(        3,0:Nloc,0:Nloc)
REAL               :: tmp2(       3,0:Nloc,0:Nloc)
!==================================================================================================================================

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

  ! Get min/max coordinate from current face
  IF(.NOT.GetMeshMinMaxBoundariesIsDone)THEN
    xyzMinMax(1) = MIN(xyzMinMax(1), MINVAL(tmp(1,:,:)))
    xyzMinMax(2) = MAX(xyzMinMax(2), MAXVAL(tmp(1,:,:)))

    xyzMinMax(3) = MIN(xyzMinMax(3), MINVAL(tmp(2,:,:)))
    xyzMinMax(4) = MAX(xyzMinMax(4), MAXVAL(tmp(2,:,:)))

    xyzMinMax(5) = MIN(xyzMinMax(5), MINVAL(tmp(3,:,:)))
    xyzMinMax(6) = MAX(xyzMinMax(6), MAXVAL(tmp(3,:,:)))
  END IF ! .NOT.GetMeshMinMaxBoundariesIsDone

  ! turn into right hand system of side
  DO q=0,Nloc; DO p=0,Nloc
    pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
    ! Compute Face_xGP for sides
    Face_xGP(1:3,p,q,SideID)=tmp2(:,pq(1),pq(2))
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
      pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
      Ja_Face_l(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    ! DEBUG old version
      !Ja_Face(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd
  IF(PRESENT(Ja_Face)) Ja_Face(:,:,:,:,SideID)=Ja_Face_l


  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)
  CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face_l,&
                         NormVec(:,:,:,SideID),TangVec1(:,:,:,SideID),&
                         TangVec2(:,:,:,SideID),SurfElem(:,:,SideID))

  !compute metrics for mortar faces, interpolate Ja_Face to small sides
  IF(MortarType(1,SideID).GT.0)THEN
    CALL Mortar_CalcSurfMetrics(SideID,Nloc,Ja_Face_l,Face_xGP(:,:,:,SideID),&
                                            Mortar_Ja,Mortar_xGP,nbSideIDs)
    DO iMortar=1,4
      SideID2=nbSideIDs(iMortar)
      IF(SideID2.LT.1) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
      IF(PRESENT(Ja_Face)) Ja_Face(:,:,:,:,SideID2)=Mortar_Ja(:,:,:,:,iMortar)
      Face_xGP(:,:,:,SideID2) = Mortar_xGP(:,:,:,iMortar)
      CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Mortar_Ja(:,:,:,:,iMortar),&
                             NormVec(:,:,:,SideID2),TangVec1(:,:,:,SideID2),&
                             TangVec2(:,:,:,SideID2),SurfElem(:,:,SideID2))
    END DO

  END IF
END DO

END SUBROUTINE CalcSurfMetrics


!==================================================================================================================================
!> Computes surface normal and tangential vectors and surface area from surface metrics Ja_Face.
!==================================================================================================================================
SUBROUTINE SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face,NormVec,TangVec1,TangVec2,SurfElem)
! MODULES
USE MOD_Globals,     ONLY: CROSS
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                       !< polynomial degree
INTEGER,INTENT(IN) :: NormalDir                  !< direction of normal vector
INTEGER,INTENT(IN) :: TangDir                    !< direction of 1. tangential vector
REAL,INTENT(IN)    :: NormalSign                 !< sign of normal vector
REAL,INTENT(IN)    :: Ja_Face(3,3,0:Nloc,0:Nloc) !< face metrics
REAL,INTENT(OUT)   ::   NormVec(3,0:Nloc,0:Nloc) !< element face normal vectors
REAL,INTENT(OUT)   ::  TangVec1(3,0:Nloc,0:Nloc) !< element face tangential vectors
REAL,INTENT(OUT)   ::  TangVec2(3,0:Nloc,0:Nloc) !< element face tangential vectors
REAL,INTENT(OUT)   ::  SurfElem(  0:Nloc,0:Nloc) !< element face surface area
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,Nloc; DO p=0,Nloc
  SurfElem(  p,q) = SQRT(SUM(Ja_Face(NormalDir,:,p,q)**2))
  NormVec( :,p,q) = NormalSign*Ja_Face(NormalDir,:,p,q)/SurfElem(p,q)
  TangVec1(:,p,q) = Ja_Face(TangDir,:,p,q) - SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q)) &
                    *NormVec(:,p,q)
  TangVec1(:,p,q) = TangVec1(:,p,q)/SQRT(SUM(TangVec1(:,p,q)**2))
  TangVec2(:,p,q) = CROSS(NormVec(:,p,q),TangVec1(:,p,q))
END DO; END DO ! p,q
END SUBROUTINE SurfMetricsFromJa

#ifdef maxwell
#if defined(ROS) || defined(IMPA)
SUBROUTINE CalcElemLocalSurfMetrics(Nloc,JaCL_N,Vdm_CLN_N,iElem)
!===================================================================================================================================
! Compute the element-local normal vectors and SurfElem from element metrics. Input is JaCL_N, the 3D element metrics on
! Cebychev-Lobatto points. The orientation of the element-local vectors correspond to the volume-DOFs ijk
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,     ONLY:CROSS
USE MOD_Mesh_Vars,   ONLY:ElemToSide
USE MOD_Mesh_Vars,   ONLY:NormalDirs,TangDirs,NormalSigns
USE MOD_Mesh_Vars,   ONLY:nVecLoc,SurfLoc
USE MOD_Mappings,    ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis, ONLY:ChangeBasis2D
USE MOD_Mortar_Metrics, ONLY:Mortar_CalcSurfMetrics
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                  !< (IN) polynomial degree
INTEGER,INTENT(IN) :: iElem                                 !< (IN) element index
REAL,INTENT(IN)    :: JaCL_N(1:3,1:3,0:Nloc,0:Nloc,0:Nloc)  !< (IN) volume metrics of element
REAL,INTENT(IN)    :: Vdm_CLN_N(   0:Nloc,0:Nloc)           !< (IN) Vandermonde matrix from Cheby-Lob on N to final nodeset on N
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,dd,iLocSide,SideID
INTEGER            :: NormalDir,TangDir
REAL               :: NormalSign
REAL               :: Ja_Face_l(3,3,0:Nloc,0:Nloc)
REAL               :: tmp(        3,0:Nloc,0:Nloc)
REAL               :: tmp2(        3,0:Nloc,0:Nloc)
!==================================================================================================================================

DO iLocSide=1,6
 ! compute the local normVec and SurfElem form each element
 ! this is currently only required for the Maxwell case (HDG should be handled similar)
 ! this allows to build the preconditioner with the correct normVecs and SurfElems
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)

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
      Ja_Face_l(dd,1:3,p,q)=tmp2(:,p,q)
    END DO; END DO ! p,q
  END DO ! dd

  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)
 ! compute Surf and normal vector in ijk orientation of volume
  DO q=0,Nloc; DO p=0,Nloc
    Surfloc (  p,q,iLocSide,iElem) = SQRT(SUM(Ja_Face_l(NormalDir,:,p,q)**2))
    nVecLoc( :,p,q,iLocSide,iElem) = NormalSign*Ja_Face_l(NormalDir,:,p,q)/SurfLoc(p,q,iLocSide,iElem)
   !TangVec1(:,p,q) = Ja_Face(TangDir,:,p,q) - SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q)) &
   !                  *NormVec(:,p,q)
   !TangVec1(:,p,q) = TangVec1(:,p,q)/SQRT(SUM(TangVec1(:,p,q)**2))
   !TangVec2(:,p,q) = CROSS(NormVec(:,p,q),TangVec1(:,p,q))
  END DO; END DO ! p,q

END DO

END SUBROUTINE CalcElemLocalSurfMetrics
#endif /*ROS or IMPA*/
#endif /*maxwell*/

SUBROUTINE CalcMetricsErrorDiff()
!===================================================================================================================================
!> This routine computes the geometries volume metric terms (sJ only!)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,               ONLY:NGeo,NGeoRef
USE MOD_Mesh_Vars,               ONLY:sJ
USE MOD_Mesh_Vars,               ONLY:nElems
USE MOD_Mesh_Vars,               ONLY:NodeCoords
USE MOD_Interpolation,           ONLY:GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_ChangeBasis,             ONLY:changeBasis3D,ChangeBasis3D_XYZ
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:NodeTypeCL,NodeTypeVISU,NodeType
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
INTEGER :: ll
! Jacobian on CL N and NGeoRef
REAL    :: DetJac_N( 1,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL    :: XCL_N(      3,  0:PP_N,0:PP_N,0:PP_N)          ! mapping X(xi) P\in N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: dXCL_Ngeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on CL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

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
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)
REAL    :: xiCL_N(0:PP_N)   ,wBaryCL_N(0:PP_N)
REAL    :: xiCL_NGeo(0:NGeo)   ,wBaryCL_NGeo(0:NGeo)

REAL    :: DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems)      !< determinant of the mesh Jacobian for each Gauss point at degree 3*NGeo
!===================================================================================================================================

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

CALL GetNodesAndWeights(NGeo   , NodeTypeCL  , XiCL_NGeo  , wIPBary=wBaryCL_NGeo)

! Outer loop over all elements
DetJac_Ref=0.
DO iElem=1,nElems
  !1.a) Transform from EQUI_Ngeo to CL points on Ngeo and N
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem)            ,XCL_Ngeo)
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,   XCL_Ngeo                             ,XCL_N)

  !1.b) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_NGeo(dd,nn,i,j,k))
  dXCL_NGeo=0.
  DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
    ! Matrix-vector multiplication
    DO ll=0,Ngeo
      dXCL_NGeo(1,:,i,j,k)=dXCL_NGeo(1,:,i,j,k) + DCL_NGeo(i,ll)*XCL_NGeo(:,ll,j,k)
      dXCL_NGeo(2,:,i,j,k)=dXCL_NGeo(2,:,i,j,k) + DCL_NGeo(j,ll)*XCL_NGeo(:,i,ll,k)
      dXCL_NGeo(3,:,i,j,k)=dXCL_NGeo(3,:,i,j,k) + DCL_NGeo(k,ll)*XCL_NGeo(:,i,j,ll)
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

  ! interpolate DetJac_ref to the solution points
  CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N)

  ! assign to global Variable sJ
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    sJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
  END DO; END DO; END DO !i,j,k=0,PP_N
END DO !iElem=1,nElems


END SUBROUTINE CalcMetricsErrorDiff

END MODULE MOD_Metrics
