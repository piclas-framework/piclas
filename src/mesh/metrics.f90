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
PUBLIC::BuildElem_xGP
PUBLIC::CalcMetrics
PUBLIC::CalcSurfMetrics
PUBLIC::SurfMetricsFromJa
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine takes the equidistant node coordinates of the mesh (on NGeo+1 points) and uses them to build the coordinates
!> of solution/interpolation points of type NodeType on polynomial degree Nloc (Nloc+1 points per direction).
!==================================================================================================================================
SUBROUTINE BuildElem_xGP(NodeCoords)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,nElems,N_VolMesh, offSetElem
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType,Nmax
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ, ChangeBasis3D
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
#endif /*!(PP_TimeDiscMethod==700)*/
#if USE_FV
USE MOD_Mesh_Vars_FV,       ONLY: Elem_xGP_PP_1,Elem_xGP_FV
#else
USE MOD_Interpolation_Vars ,ONLY: Nmin
#endif /*USE_FV*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems)         !< Equidistant mesh coordinates
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,Nloc

TYPE VdmType
  REAL, ALLOCATABLE           :: Vdm_EQNGeo_CLNloc(:,:)
  REAL, ALLOCATABLE           :: Vdm_CLNloc_Nloc  (:,:)
END TYPE VdmType

TYPE(VdmType), DIMENSION(:), ALLOCATABLE :: Vdm
!==================================================================================================================================

#if USE_FV
! Element centers
SDEALLOCATE(Elem_xGP_FV)
ALLOCATE(Elem_xGP_FV   (3,0:0,0:0,0:0,nElems))!
! Output points
SDEALLOCATE(Elem_xGP_PP_1)
ALLOCATE(Elem_xGP_PP_1 (3,0:PP_1,0:PP_1,0:PP_1,nElems))
ASSOCIATE( Nmin => 0 , Nmax => MAX(Nmax,1) )
#endif /*USE_FV*/
! Build Vdm for every degree
ALLOCATE(Vdm(Nmin:Nmax))
DO Nloc = Nmin, Nmax
  ALLOCATE(Vdm(Nloc)%Vdm_EQNGeo_CLNloc(0:Nloc,0:NGeo))
  ALLOCATE(Vdm(Nloc)%Vdm_CLNloc_Nloc(0:Nloc,0:Nloc))
  CALL GetVandermonde(NGeo, NodeTypeVISU, NLoc, NodeTypeCL, Vdm(Nloc)%Vdm_EQNGeo_CLNloc,  modal=.FALSE.)
  CALL GetVandermonde(Nloc, NodeTypeCL  , Nloc, NodeType  , Vdm(Nloc)%Vdm_CLNloc_Nloc,     modal=.FALSE.)

! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!       Important for curved meshes if NGeo<N, no effect for N>=NGeo

!1.a) Transform from EQUI_NGeo to solution points on Nloc
  Vdm(Nloc)%Vdm_EQNGeo_CLNloc=MATMUL(Vdm(Nloc)%Vdm_CLNloc_Nloc, Vdm(Nloc)%Vdm_EQNGeo_CLNloc)
END DO ! Nloc = Nmin, Nmax
#if USE_FV
END ASSOCIATE
#endif /*USE_FV*/

! Set Elem_xGP for each element
DO iElem=1,nElems
#if !(PP_TimeDiscMethod==700)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
#else
  Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/
  ALLOCATE(N_VolMesh(iElem)%Elem_xGP(3,0:Nloc,0:Nloc,0:Nloc))
  !WRITE (*,*) "NodeCoords(:,:,:,:,iElem) =", NodeCoords(:,:,:,:,iElem)
  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm(Nloc)%Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),N_VolMesh(iElem)%Elem_xGP(:,:,:,:))
#if USE_FV
  ! Element centers
  Nloc = 0
  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm(Nloc)%Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),Elem_xGP_FV(:,:,:,:,iElem))
  ! Output points
  Nloc = PP_1
  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm(Nloc)%Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),Elem_xGP_PP_1(:,:,:,:,iElem))
#endif /*USE_FV*/
END DO

END SUBROUTINE BuildElem_xGP


SUBROUTINE CalcMetrics(XCL_NGeo_Out,dXCL_NGeo_out)
!===================================================================================================================================
!> This routine computes the geometries volume metric terms.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,NGeoRef,nGlobalElems,xyzMinMax,GetMeshMinMaxBoundariesIsDone
USE MOD_Mesh_Vars          ,ONLY: crossProductMetrics, offSetElem
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Mesh_Vars          ,ONLY: DetJac_Ref
USE MOD_Mesh_Vars          ,ONLY: crossProductMetrics
USE MOD_Mesh_Vars          ,ONLY: NodeCoords,N_VolMesh,N_VolMesh2
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_ChangeBasis        ,ONLY: changeBasis3D,ChangeBasis3D_XYZ
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType,Nmin,Nmax,NInfo
USE MOD_ReadInTools        ,ONLY: GETLOGICAL
USE MOD_Globals_Vars       ,ONLY: PI
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars          ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Basis              ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights
#if USE_HDG
USE MOD_Symmetry_Vars      ,ONLY: Symmetry
USE MOD_Globals_Vars       ,ONLY: PI
#endif /*USE_HDG*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT),OPTIONAL  :: XCL_NGeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)      ! mapping X(xi) P\in Ngeo
REAL ,INTENT(INOUT),OPTIONAL :: dXCL_NGeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)   ! jacobi matrix on CL Ngeo
!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,q,iElem
INTEGER :: ll
! interpolation points and derivatives on CL N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:Ngeo)          ! mapping X(xi) P\in Ngeo
REAL    :: dXCL_NGeo(  3,3,0:Ngeo,0:Ngeo,0:Ngeo)          ! jacobi matrix on CL Ngeo
REAL    :: dX_NgeoRef( 3,3,0:NgeoRef,0:NgeoRef,0:NgeoRef) ! jacobi matrix on SOL NgeoRef

REAL    :: scaledJac(2)

! Polynomial derivation matrices
REAL    :: DCL_NGeo(0:Ngeo,0:Ngeo)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNGeo_CLNgeo( 0:Ngeo   ,0:Ngeo)
REAL    :: Vdm_CLNGeo_NGeoRef(0:NgeoRef,0:Ngeo)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL    :: xiRef( 0:NgeoRef),wBaryRef( 0:NgeoRef)

REAL               :: StartT,EndT
LOGICAL            :: meshCheckRef
REAL,ALLOCATABLE   :: scaledJacRef(:,:,:)
REAL               :: SmallestscaledJacRef
REAL,PARAMETER     :: scaledJacRefTol=0.01

INTEGER           :: Nloc

!===================================================================================================================================
StartT=PICLASTIME()

ALLOCATE(NInfo(1:Nmax))

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
CALL GetVandermonde(    Ngeo   , NodeTypeVISU, Ngeo    , NodeTypeCL, Vdm_EQNGeo_CLNGeo , modal=.FALSE.)

! 1.b) dXCL_NGeo:
CALL GetDerivativeMatrix(Ngeo  , NodeTypeCL  , DCL_NGeo)

! 1.c) Jacobian: CLNgeo to NgeoRef, CLNgeoRef to N
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNgeo_NGeoRef, modal=.FALSE.)
CALL GetNodesAndWeights(NgeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)

! Initialize Vandermonde and D matrices
! Only use modal Vandermonde for terms that need to be conserved as Jacobian if N_out>PP_N
! Always use interpolation for the rest!
DO Nloc = Nmin, Nmax
  ! interpolation points and derivatives on CL N
  ALLOCATE( NInfo(Nloc)%XCL_N(      3,  0:Nloc,0:Nloc,0:Nloc) )  ! mapping X(xi) P\in N
  ALLOCATE( NInfo(Nloc)%R_CL_N(     3,3,0:Nloc,0:Nloc,0:Nloc) )  ! buffer for metric terms, uses XCL_N,dXCL_N
  ALLOCATE( NInfo(Nloc)%JaCL_N(     3,3,0:Nloc,0:Nloc,0:Nloc) )  ! metric terms P\in N
  ALLOCATE( NInfo(Nloc)%DetJac_N(     1,0:Nloc,0:Nloc,0:Nloc) )
  ALLOCATE( NInfo(Nloc)%DCL_N(          0:Nloc,0:Nloc)        )
  ALLOCATE( NInfo(Nloc)%Vdm_NgeoRef_N(  0:Nloc,0:NgeoRef)     )
  ALLOCATE( NInfo(Nloc)%Vdm_CLNGeo_CLN( 0:Nloc,0:Ngeo)        )
  ALLOCATE( NInfo(Nloc)%Vdm_CLN_N(      0:Nloc,0:Nloc)        )

  ! 1.c) Jacobian: CLNgeo to NgeoRef, CLNgeoRef to N
  CALL GetVandermonde(    NgeoRef, NodeType    , Nloc    , NodeType  , NInfo(Nloc)%Vdm_NgeoRef_N     , modal=.TRUE.)

! 1.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
  CALL GetVandermonde(    Ngeo   , NodeTypeCL  , Nloc    , NodeTypeCL, NInfo(Nloc)%Vdm_CLNGeo_CLN    , modal=.FALSE.)
  CALL GetDerivativeMatrix(Nloc  , NodeTypeCL  , NInfo(Nloc)%DCL_N)

! 2.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
  CALL GetVandermonde(    Nloc   , NodeTypeCL  , Nloc    , NodeType,   NInfo(Nloc)%Vdm_CLN_N         , modal=.FALSE.)

END DO ! Nloc = Nmin, Nmax

! Chebyshev-Lobatto NGeo: XiCL_NGeo, wBaryCL_NGeo
ALLOCATE(wBaryCL_NGeo(0:NGeo))
ALLOCATE(XiCL_NGeo(0:NGeo))
CALL ChebyGaussLobNodesAndWeights(NGeo,XiCL_NGeo)
CALL BarycentricWeights(NGeo,XiCL_NGeo,wBaryCL_NGeo)

! Outer loop over all elements
DetJac_Ref=0.
DO iElem=1,nElems
SmallestscaledJacRef=HUGE(1.)
  N_VolMesh2(iElem)%dXCL_N=0.
  ! Get N
#if !(PP_TimeDiscMethod==700)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
#else
  Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/

  ! Init
  N_VolMesh(iElem)%Metrics_fTilde=0.
  N_VolMesh(iElem)%Metrics_gTilde=0.
  N_VolMesh(iElem)%Metrics_hTilde=0.

  !1.a) Transform from EQUI_Ngeo to CL points on Ngeo and N
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem)            ,XCL_Ngeo)
  CALL ChangeBasis3D(3, NGeo, Nloc, NInfo(Nloc)%Vdm_CLNGeo_CLN, XCL_Ngeo                 , NInfo(Nloc)%XCL_N)
  ! Save XCL_N for LB communication
  N_VolMesh(iElem)%XCL_N = NInfo(Nloc)%XCL_N

  !1.b) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_NGeo(dd,nn,i,j,k))
  dXCL_NGeo=0.
  DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
    ! Matrix-vector multiplication
    DO ll=0,Ngeo
      dXCL_NGeo(1,:,i,j,k)=dXCL_NGeo(1,:,i,j,k) + DCL_NGeo(i,ll)*XCL_Ngeo(:,ll,j,k)
      dXCL_NGeo(2,:,i,j,k)=dXCL_NGeo(2,:,i,j,k) + DCL_NGeo(j,ll)*XCL_Ngeo(:,i,ll,k)
      dXCL_NGeo(3,:,i,j,k)=dXCL_NGeo(3,:,i,j,k) + DCL_NGeo(k,ll)*XCL_Ngeo(:,i,j,ll)
    END DO !l=0,N
#if USE_HDG
    ! AXISYMMETRIC HDG
    IF(Symmetry%Axisymmetric) dXCL_Ngeo(3,3,i,j,k)=PI*XCL_Ngeo(2,i,j,k)
#endif /*USE_HDG*/
  END DO; END DO; END DO !i,j,k=0,Ngeo

  ! 1.c)Jacobians! grad(X_1) (grad(X_2) x grad(X_3))
  ! Compute Jacobian on NGeo and then interpolate:
  ! required to guarantee conservation when restarting with N<NGeo
  CALL ChangeBasis3D(3,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(:,1,:,:,:),dX_NGeoRef(:,1,:,:,:))
  CALL ChangeBasis3D(3,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(:,2,:,:,:),dX_NGeoRef(:,2,:,:,:))
  CALL ChangeBasis3D(3,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(:,3,:,:,:),dX_NGeoRef(:,3,:,:,:))
  DO k=0,NgeoRef; DO j=0,NgeoRef; DO i=0,NgeoRef
    DetJac_Ref(1,i,j,k,iElem)=DetJac_Ref(1,i,j,k,iElem) &
      + dX_NGeoRef(1,1,i,j,k)*(dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(3,3,i,j,k) - dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(2,3,i,j,k))  &
      + dX_NGeoRef(2,1,i,j,k)*(dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(1,3,i,j,k) - dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(3,3,i,j,k))  &
      + dX_NGeoRef(3,1,i,j,k)*(dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(2,3,i,j,k) - dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(1,3,i,j,k))
  END DO; END DO; END DO !i,j,k=0,NgeoRef

  ! Check Jacobians in Ref already (no good if we only go on because N doesn't catch misbehaving points)
  IF (meshCheckRef) THEN
    ALLOCATE(scaledJacRef(0:NGeoRef,0:NGeoRef,0:NGeoRef))
    DO k=0,NGeoRef; DO j=0,NGeoRef; DO i=0,NGeoRef
      scaledJacRef(i,j,k)=DetJac_Ref(1,i,j,k,iElem)/MAXVAL(DetJac_Ref(1,:,:,:,iElem))
      SmallestscaledJacRef=MIN(SmallestscaledJacRef,scaledJacRef(i,j,k))
      IF(scaledJacRef(i,j,k).LT.scaledJacRefTol) THEN
        WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJacRef(i,j,k)
        WRITE(Unit_StdOut,*) 'Coords near:', N_VolMesh(iElem)%Elem_xGP(:,INT(Nloc/2),INT(Nloc/2),INT(Nloc/2))
        WRITE(Unit_StdOut,*) 'This check is optional. You can disable it by setting meshCheckRef = F'
        CALL abort(__STAMP__,'Scaled Jacobian in reference system lower then tolerance in global element:',iElem+offsetElem)
      END IF
    END DO; END DO; END DO !i,j,k=0,N
    DEALLOCATE(scaledJacRef)
  END IF

  ! project detJac_ref onto the solution basis
  CALL ChangeBasis3D(1,NgeoRef,Nloc,NInfo(Nloc)%Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:,iElem),NInfo(Nloc)%DetJac_N)

  ! assign to global Variable sJ
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    N_VolMesh(iElem)%sJ(i,j,k)=1./NInfo(Nloc)%DetJac_N(1,i,j,k)
  END DO; END DO; END DO !i,j,k=0,Nloc

  ! check for negative Jacobians
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    IF(NInfo(Nloc)%DetJac_N(1,i,j,k).LE.0.)&
      WRITE(Unit_StdOut,*) 'Negative Jacobian found on Gauss point. Coords:', N_VolMesh(iElem)%Elem_xGP(:,i,j,k)
  END DO; END DO; END DO !i,j,k=0,N
  ! check scaled Jacobians
  scaledJac(2)=MINVAL(NInfo(Nloc)%DetJac_N(1,:,:,:))/MAXVAL(NInfo(Nloc)%DetJac_N(1,:,:,:))
  IF(scaledJac(2).LT.0.01) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac
    CALL abort(__STAMP__,'Scaled Jacobian lower then tolerance in global element:',iElem+offsetElem)
  END IF

  !2.a) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_N(dd,nn,i,j,k))
  ! N>=Ngeo: interpolate from dXCL_NGeo (default)
  ! N< Ngeo: directly derive XCL_N
  IF(Nloc.GE.NGeo)THEN !compute first derivative on Ngeo and then interpolate
    CALL ChangeBasis3D(3,NGeo,Nloc,NInfo(Nloc)%Vdm_CLNGeo_CLN,dXCL_NGeo(:,1,:,:,:),N_VolMesh2(iElem)%dXCL_N(:,1,:,:,:))
    CALL ChangeBasis3D(3,NGeo,Nloc,NInfo(Nloc)%Vdm_CLNGeo_CLN,dXCL_NGeo(:,2,:,:,:),N_VolMesh2(iElem)%dXCL_N(:,2,:,:,:))
    CALL ChangeBasis3D(3,NGeo,Nloc,NInfo(Nloc)%Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:),N_VolMesh2(iElem)%dXCL_N(:,3,:,:,:))
  ELSE  !N<Ngeo: first interpolate and then compute derivative (important if curved&periodic)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      ! Matrix-vector multiplication
      ASSOCIATE(dXCL => N_VolMesh2(iElem)%dXCL_N(:,:,i,j,k))
      DO ll=0,Nloc
        dXCL(1,:)=dXCL(1,:) + NInfo(Nloc)%DCL_N(i,ll)*NInfo(Nloc)%XCL_N(:,ll,j,k)
        dXCL(2,:)=dXCL(2,:) + NInfo(Nloc)%DCL_N(j,ll)*NInfo(Nloc)%XCL_N(:,i,ll,k)
        dXCL(3,:)=dXCL(3,:) + NInfo(Nloc)%DCL_N(k,ll)*NInfo(Nloc)%XCL_N(:,i,j,ll)
      END DO !l=0,N
#if USE_HDG
      ! AXISYMMETRIC HDG
      IF(Symmetry%Axisymmetric) dXCL(:,3)=PI * NInfo(Nloc)%XCL_N(2,i,j,k)
#endif /*USE_HDG*/
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=Ngeo

  NInfo(Nloc)%JaCL_N=0.
  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      ASSOCIATE(dXCL => N_VolMesh2(iElem)%dXCL_N(:,:,i,j,k))
      ! exact (cross-product) form
      ! Ja(:)^nn = ( d/dxi_(nn+1) XCL_N(:) ) x (d/xi_(nn+2) XCL_N(:))
      !
      ! NInfo(Nloc)%JaCL_N(dd,nn) = dXCL_N(dd+1,nn+1)*dXCL_N(dd+2,nn+2) -dXCL_N(dd+1,nn+2)*dXCL_N(dd+2,nn+1)
      NInfo(Nloc)%JaCL_N(1,1,i,j,k)=dXCL(2,2)*dXCL(3,3) - dXCL(2,3)*dXCL(3,2)
      NInfo(Nloc)%JaCL_N(2,1,i,j,k)=dXCL(3,2)*dXCL(1,3) - dXCL(3,3)*dXCL(1,2)
      NInfo(Nloc)%JaCL_N(3,1,i,j,k)=dXCL(1,2)*dXCL(2,3) - dXCL(1,3)*dXCL(2,2)
      NInfo(Nloc)%JaCL_N(1,2,i,j,k)=dXCL(2,3)*dXCL(3,1) - dXCL(2,1)*dXCL(3,3)
      NInfo(Nloc)%JaCL_N(2,2,i,j,k)=dXCL(3,3)*dXCL(1,1) - dXCL(3,1)*dXCL(1,3)
      NInfo(Nloc)%JaCL_N(3,2,i,j,k)=dXCL(1,3)*dXCL(2,1) - dXCL(1,1)*dXCL(2,3)
      NInfo(Nloc)%JaCL_N(1,3,i,j,k)=dXCL(2,1)*dXCL(3,2) - dXCL(2,2)*dXCL(3,1)
      NInfo(Nloc)%JaCL_N(2,3,i,j,k)=dXCL(3,1)*dXCL(1,2) - dXCL(3,2)*dXCL(1,1)
      NInfo(Nloc)%JaCL_N(3,3,i,j,k)=dXCL(1,1)*dXCL(2,2) - dXCL(1,2)*dXCL(2,1)
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  ELSE ! curl metrics
    ! AXISYMMETRIC HDG: Does not work without some work :(
    ! invariant curl form, as cross product: R^dd = 1/2( XCL_N(:) x (d/dxi_dd XCL_N(:)))
    !
    !NInfo(Nloc)%R_CL_N(dd,nn)=1/2*( XCL_N(nn+2)* d/dxi_dd XCL_N(nn+1) - XCL_N(nn+1)* d/dxi_dd XCL_N(nn+2))
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      ASSOCIATE(dXCL => N_VolMesh2(iElem)%dXCL_N(:,:,i,j,k))
      NInfo(Nloc)%R_CL_N(:,1,i,j,k)=0.5*(NInfo(Nloc)%XCL_N(3,i,j,k)*dXCL(:,2) - NInfo(Nloc)%XCL_N(2,i,j,k)*dXCL(:,3) )
      NInfo(Nloc)%R_CL_N(:,2,i,j,k)=0.5*(NInfo(Nloc)%XCL_N(1,i,j,k)*dXCL(:,3) - NInfo(Nloc)%XCL_N(3,i,j,k)*dXCL(:,1) )
      NInfo(Nloc)%R_CL_N(:,3,i,j,k)=0.5*(NInfo(Nloc)%XCL_N(2,i,j,k)*dXCL(:,1) - NInfo(Nloc)%XCL_N(1,i,j,k)*dXCL(:,2) )
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
    ! Metrics are the curl of R:  Ja(:)^nn = -(curl R_CL(:,nn))
    ! JaCL_N(dd,nn)= -[d/dxi_(dd+1) RCL(dd+2,nn) - d/dxi_(dd+2) RCL(dd+1,nn) ]
    !              =   d/dxi_(dd+2) RCL(dd+1,nn) - d/dxi_(dd+1) RCL(dd+2,nn)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      ASSOCIATE(JaCL => NInfo(Nloc)%JaCL_N(:,:,i,j,k))
      DO q=0,Nloc
        JaCL(1,:)=JaCL(1,:) - NInfo(Nloc)%DCL_N(j,q)*NInfo(Nloc)%R_CL_N(3,:,i,q,k)
        JaCL(2,:)=JaCL(2,:) - NInfo(Nloc)%DCL_N(k,q)*NInfo(Nloc)%R_CL_N(1,:,i,j,q)
        JaCL(3,:)=JaCL(3,:) - NInfo(Nloc)%DCL_N(i,q)*NInfo(Nloc)%R_CL_N(2,:,q,j,k)
      END DO!q=0,Nloc
      DO q=0,Nloc
        JaCL(1,:)=JaCL(1,:) + NInfo(Nloc)%DCL_N(k,q)*NInfo(Nloc)%R_CL_N(2,:,i,j,q)
        JaCL(2,:)=JaCL(2,:) + NInfo(Nloc)%DCL_N(i,q)*NInfo(Nloc)%R_CL_N(3,:,q,j,k)
        JaCL(3,:)=JaCL(3,:) + NInfo(Nloc)%DCL_N(j,q)*NInfo(Nloc)%R_CL_N(1,:,i,q,k)
      END DO!q=0,Nloc
      END ASSOCIATE
! same with only one loop, gives different roundoff ...
!      DO q=0,Nloc
!        JaCL_N(1,:,i,j,k)=JaCL_N(1,:,i,j,k) - DCL_N(j,q)*R_CL_N(3,:,i,q,k) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
!        JaCL_N(2,:,i,j,k)=JaCL_N(2,:,i,j,k) - DCL_N(k,q)*R_CL_N(1,:,i,j,q) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
!        JaCL_N(3,:,i,j,k)=JaCL_N(3,:,i,j,k) - DCL_N(i,q)*R_CL_N(2,:,q,j,k) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
!      END DO!q=0,Nloc
    END DO; END DO; END DO !i,j,k=0,N
  END IF !crossProductMetrics

  ! Save JaCL_N for LB communication
  N_VolMesh2(iElem)%JaCL_N = NInfo(Nloc)%JaCL_N

  ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
  CALL ChangeBasis3D(3,Nloc,Nloc,NInfo(Nloc)%Vdm_CLN_N,NInfo(Nloc)%JaCL_N(1,:,:,:,:),N_VolMesh(iElem)%Metrics_fTilde(:,:,:,:))
  CALL ChangeBasis3D(3,Nloc,Nloc,NInfo(Nloc)%Vdm_CLN_N,NInfo(Nloc)%JaCL_N(2,:,:,:,:),N_VolMesh(iElem)%Metrics_gTilde(:,:,:,:))
  CALL ChangeBasis3D(3,Nloc,Nloc,NInfo(Nloc)%Vdm_CLN_N,NInfo(Nloc)%JaCL_N(3,:,:,:,:),N_VolMesh(iElem)%Metrics_hTilde(:,:,:,:))

  !CALL CalcSurfMetrics(iElem)
  ! particle mapping
  IF(PRESENT(XCL_Ngeo_Out))   XCL_NGeo_Out(1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem)     = XCL_Ngeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo)
  IF(PRESENT(dXCL_ngeo_out)) dXCL_NGeo_Out(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo,iElem) = dXCL_NGeo(1:3,1:3,0:Ngeo,0:Ngeo,0:Ngeo)
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
CALL DisplayMessageAndTime(EndT-StartT, 'Calculation of metrics took ', DisplayLine=.FALSE.)

END SUBROUTINE CalcMetrics


SUBROUTINE CalcSurfMetrics(iElem)
!===================================================================================================================================
! Compute normal and tangential vectors from element metrics. Input is JaCL_N, the 3D element metrics on Cebychev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,        ONLY:CROSS
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,MortarType,xyzMinMax,GetMeshMinMaxBoundariesIsDone!,nSides
USE MOD_Mesh_Vars          ,ONLY: NormalDirs,TangDirs,NormalSigns, N_SurfMesh
USE MOD_Mappings,       ONLY:CGNS_SideToVol2
USE MOD_ChangeBasis,    ONLY:ChangeBasis2D
USE MOD_Mortar_Metrics, ONLY:Mortar_CalcSurfMetrics
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_Interpolation_Vars ,ONLY: Nmax,NInfo!,PREF_VDM,N_Inter
USE MOD_Mesh_Vars,          ONLY: SideToElem, offSetElem,N_VolMesh,N_VolMesh2
USE MOD_Interpolation_Vars ,ONLY: PREF_VDM!,N_Inter
! #if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)) && defined(maxwell)
! USE MOD_Equation_Vars      ,ONLY: DoExactFlux ! Required for skipping cycle because NormVec is then not built for loc.LT.NSideMax sides
! #endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)) && defined(maxwell)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem                                !< (IN) element index

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID,SideID2,iMortar,nbSideIDs(4)
INTEGER            :: NormalDir,TangDir,NSideMax
REAL               :: NormalSign
REAL               :: Ja_Face_l(3,3,0:Nmax,0:Nmax)
REAL               :: Mortar_Ja(3,3,0:Nmax,0:Nmax,4), Mortar_Ja_loc(3,3,0:Nmax,0:Nmax)
REAL               :: Mortar_xGP( 3,0:Nmax,0:Nmax,4), Mortar_xGP_loc( 3,0:Nmax,0:Nmax)
REAL               :: tmp(        3,0:Nmax,0:Nmax)
REAL               :: tmp2(       3,0:Nmax,0:Nmax)
REAL               :: tmpflip(    3,0:Nmax,0:Nmax)
INTEGER            :: Nloc,flip,NSideMortar!,NSideMin,NSideMax
LOGICAL            :: flipSide
!#if USE_HDG
!INTEGER            :: iSide
!#endif /*USE_HDG*/
!==================================================================================================================================

#if PP_dim == 3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  flipSide=.FALSE.
#if !(PP_TimeDiscMethod==700)
  ! Use maximum polynomial degree of master/slave sides
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  NSideMax = MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID))
  !WRITE (*,*) "Nloc,NSideMax,ElemToSide(E2S_FLIP,iLocSide,iElem) =", Nloc,NSideMax,ElemToSide(E2S_FLIP,iLocSide,iElem)

  ! TODO: maybe this has to be done differently between HDG and Maxwell
  IF(Nloc.LT.NSideMax)THEN
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)) && defined(maxwell)
    ! Do not cycle here when exact flux is used because NormVec is then not built for loc.LT.NSideMax sides
!    IF(.NOT.DoExactFlux)&
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400)) && defined(maxwell)*/
!    CYCLE
  ELSE
    IF(DG_Elems_master(SideID).EQ.DG_Elems_slave(SideID).AND.(ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) ) CYCLE ! only master sides with flip=0
    IF(ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) flipSide=.TRUE.
  END IF
  !IF (ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0) CYCLE
  !Nloc = NSideMax
#else
  Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/

  SELECT CASE(iLocSide)
  CASE(XI_MINUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , 0      , 0:Nloc , 0:Nloc )
  CASE(XI_PLUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , Nloc   , 0:Nloc , 0:Nloc )
  CASE(ETA_MINUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , 0:Nloc , 0      , 0:Nloc )
  CASE(ETA_PLUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , 0:Nloc , Nloc   , 0:Nloc )
  CASE(ZETA_MINUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , 0:Nloc , 0:Nloc , 0      )
  CASE(ZETA_PLUS)
    tmp(1:3 , 0:Nloc , 0:Nloc)=N_VolMesh(iElem)%XCL_N(1:3 , 0:Nloc , 0:Nloc , Nloc   )
  END SELECT

  ! slave sides with higher polynomial degree
  IF(flipSide)THEN
    flip = SideToElem(S2E_FLIP,SideID)
    !WRITE (*,*) "SideID,ilocside,Nloc,MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID)) =", SideID,ilocside,Nloc,MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID))
    !WRITE (*,*) "flip =", flip
    !WRITE (*,*) "flip,ElemToSide(E2S_FLIP,iLocSide,iElem) =", flip,ElemToSide(E2S_FLIP,iLocSide,iElem)
    !read*
    tmpflip(1:3 , 0:Nloc , 0:Nloc) = tmp(1:3 , 0:Nloc , 0:Nloc)
    SELECT CASE(flip)
      !CASE(0) ! master side
      !  Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,q,p)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,Nloc-p,q)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,Nloc-q,Nloc-p)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,p,Nloc-q)
          END DO ! p
        END DO ! q
    END SELECT
  END IF ! flipSide


  CALL ChangeBasis2D(3, Nloc, Nloc, NInfo(Nloc)%Vdm_CLN_N, tmp(1:3,0:Nloc,0:Nloc), tmp2(1:3,0:Nloc,0:Nloc))

  ! Get min/max coordinate from current face
  IF(.NOT.GetMeshMinMaxBoundariesIsDone)THEN
    xyzMinMax(1) = MIN(xyzMinMax(1), MINVAL(tmp(1,0:Nloc,0:Nloc)))
    xyzMinMax(2) = MAX(xyzMinMax(2), MAXVAL(tmp(1,0:Nloc,0:Nloc)))

    xyzMinMax(3) = MIN(xyzMinMax(3), MINVAL(tmp(2,0:Nloc,0:Nloc)))
    xyzMinMax(4) = MAX(xyzMinMax(4), MAXVAL(tmp(2,0:Nloc,0:Nloc)))

    xyzMinMax(5) = MIN(xyzMinMax(5), MINVAL(tmp(3,0:Nloc,0:Nloc)))
    xyzMinMax(6) = MAX(xyzMinMax(6), MAXVAL(tmp(3,0:Nloc,0:Nloc)))
  END IF ! .NOT.GetMeshMinMaxBoundariesIsDone

  ! turn into right hand system of side
  DO q=0,Nloc; DO p=0,Nloc
    pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
    ! Compute Face_xGP for sides
    N_SurfMesh(SideID)%Face_xGP(1:3,p,q)=tmp2(:,pq(1),pq(2))
  END DO; END DO ! p,q

  DO dd=1,3
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,0   ,:   ,:   )
    CASE(XI_PLUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,Nloc,:   ,:   )
    CASE(ETA_MINUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,:   ,0   ,:   )
    CASE(ETA_PLUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,:   ,Nloc,:   )
    CASE(ZETA_MINUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,:   ,:   ,0   )
    CASE(ZETA_PLUS)
      tmp(1:3,0:Nloc,0:Nloc)=N_VolMesh2(iElem)%JaCL_N(dd,1:3,:   ,:   ,Nloc)
    END SELECT

  ! slave sides with higher polynomial degree
  IF(flipSide)THEN
    flip = SideToElem(S2E_FLIP,SideID)
    tmpflip(1:3 , 0:Nloc , 0:Nloc) = tmp(1:3 , 0:Nloc , 0:Nloc)
    SELECT CASE(flip)
      !CASE(0) ! master side
      !  Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,q,p)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,Nloc-p,q)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,Nloc-q,Nloc-p)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,Nloc
          DO p=0,Nloc
            tmp(:,p,q)=tmpflip(:,p,Nloc-q)
          END DO ! p
        END DO ! q
    END SELECT
  END IF ! flipSide

    CALL ChangeBasis2D(3, Nloc, Nloc, NInfo(Nloc)%Vdm_CLN_N, tmp(1:3,0:Nloc,0:Nloc), tmp2(1:3,0:Nloc,0:Nloc))
    ! turn into right hand system of side
    DO q=0,Nloc; DO p=0,Nloc
      pq=CGNS_SideToVol2(Nloc,p,q,iLocSide)
      Ja_Face_l(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    ! DEBUG old version
      !Ja_Face(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd

  NormalDir  = NormalDirs(iLocSide)
  TangDir    = TangDirs(iLocSide)
  NormalSign = NormalSigns(iLocSide)
  IF(flipSide)NormalSign=-NormalSign
  CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face_l(1:3,1:3,0:Nloc,0:Nloc),&
                         N_SurfMesh(SideID)%NormVec(:,:,:),N_SurfMesh(SideID)%TangVec1(:,:,:),&
                         N_SurfMesh(SideID)%TangVec2(:,:,:),N_SurfMesh(SideID)%SurfElem(:,:))

  !compute metrics for mortar faces, interpolate Ja_Face to small sides
  IF(MortarType(1,SideID).GT.0)THEN
    CALL Mortar_CalcSurfMetrics(SideID,Nloc,Ja_Face_l(1:3,1:3,0:Nloc,0:Nloc),N_SurfMesh(SideID)%Face_xGP(:,:,:),&
                                            Mortar_Ja(1:3,1:3,0:Nloc,0:Nloc,1:4),Mortar_xGP(1:3,0:Nloc,0:Nloc,1:4),nbSideIDs)
    DO iMortar=1,4
      SideID2=nbSideIDs(iMortar)
      IF(SideID2.LT.1) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
      NSideMortar = MAX(DG_Elems_slave(SideID2),DG_Elems_master(SideID2))
      IF (Nloc.EQ.NSideMortar) THEN
        N_SurfMesh(SideID2)%Face_xGP(:,:,:) = Mortar_xGP(1:3,0:Nloc,0:Nloc,iMortar)
        CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Mortar_Ja(1:3,1:3,0:Nloc,0:Nloc,iMortar),&
                               N_SurfMesh(SideID2)%NormVec(:,:,:),N_SurfMesh(SideID2)%TangVec1(:,:,:),&
                               N_SurfMesh(SideID2)%TangVec2(:,:,:),N_SurfMesh(SideID2)%SurfElem(:,:))
      ELSE IF (NSideMortar.GT.NLoc) THEN
        CALL ChangeBasis2D(3, NLoc, NSideMortar, PREF_VDM(NLoc, NSideMortar)%Vdm, &
                 Mortar_xGP(1:3,0:Nloc,0:Nloc,iMortar), Mortar_xGP_loc(1:3,0:NSideMortar,0:NSideMortar))
        N_SurfMesh(SideID2)%Face_xGP(:,:,:) = Mortar_xGP_loc(1:3,0:NSideMortar,0:NSideMortar)
        CALL ChangeBasis2D(3, NLoc, NSideMortar, PREF_VDM(NLoc, NSideMortar)%Vdm, &
                 Mortar_Ja(1:3,1,0:Nloc,0:Nloc,iMortar), Mortar_Ja_loc(1:3,1,0:NSideMortar,0:NSideMortar))
        CALL ChangeBasis2D(3, NLoc, NSideMortar, PREF_VDM(NLoc, NSideMortar)%Vdm, &
                 Mortar_Ja(1:3,2,0:Nloc,0:Nloc,iMortar), Mortar_Ja_loc(1:3,2,0:NSideMortar,0:NSideMortar))
        CALL ChangeBasis2D(3, NLoc, NSideMortar, PREF_VDM(NLoc, NSideMortar)%Vdm, &
                 Mortar_Ja(1:3,3,0:Nloc,0:Nloc,iMortar), Mortar_Ja_loc(1:3,3,0:NSideMortar,0:NSideMortar))
        CALL SurfMetricsFromJa(NSideMortar,NormalDir,TangDir,NormalSign,Mortar_Ja_loc(1:3,1:3,0:NSideMortar,0:NSideMortar),&
                               N_SurfMesh(SideID2)%NormVec(:,:,:),N_SurfMesh(SideID2)%TangVec1(:,:,:),&
                               N_SurfMesh(SideID2)%TangVec2(:,:,:),N_SurfMesh(SideID2)%SurfElem(:,:))
      ELSE IF (NSideMortar.LT.NLoc) THEN
!        CALL ChangeBasis2D(3, NLoc, NLoc, N_Inter(NLoc)%sVdm_Leg, &
!                Mortar_xGP(1:3,0:Nloc,0:Nloc,iMortar), Mortar_xGP(1:3,0:Nloc,0:Nloc,iMortar))
!        ! Switch back to nodal basis but cut-off the higher-order DOFs (go from N_master to N_master)
!        CALL ChangeBasis2D(3, NSideMortar, NSideMortar, N_Inter(NSideMortar)%Vdm_Leg, &
!                Mortar_xGP(1:3, 0:NSideMortar, 0:NSideMortar ,iMortar), Mortar_xGP_loc(1:3,0:NSideMortar,0:NSideMortar))
!        N_SurfMesh(SideID2)%Face_xGP(:,:,:) = Mortar_xGP_loc(1:3,0:NSideMortar,0:NSideMortar)
!        CALL ChangeBasis2D(3, NLoc, NLoc, N_Inter(NLoc)%sVdm_Leg, &
!                Mortar_Ja(1:3,1,0:Nloc,0:Nloc,iMortar), Mortar_Ja(1:3,1,0:Nloc,0:Nloc,iMortar))
!        CALL ChangeBasis2D(3, NLoc, NLoc, N_Inter(NLoc)%sVdm_Leg, &
!                Mortar_Ja(1:3,2,0:Nloc,0:Nloc,iMortar), Mortar_Ja(1:3,2,0:Nloc,0:Nloc,iMortar))
!        CALL ChangeBasis2D(3, NLoc, NLoc, N_Inter(NLoc)%sVdm_Leg, &
!                Mortar_Ja(1:3,3,0:Nloc,0:Nloc,iMortar), Mortar_Ja(1:3,3,0:Nloc,0:Nloc,iMortar))
!        CALL ChangeBasis2D(3, NSideMortar, NSideMortar, N_Inter(NSideMortar)%Vdm_Leg, &
!                Mortar_Ja(1:3,1, 0:NSideMortar, 0:NSideMortar ,iMortar), Mortar_Ja_loc(1:3,1,0:NSideMortar,0:NSideMortar))
!        CALL ChangeBasis2D(3, NSideMortar, NSideMortar, N_Inter(NSideMortar)%Vdm_Leg, &
!                Mortar_Ja(1:3,2, 0:NSideMortar, 0:NSideMortar ,iMortar), Mortar_Ja_loc(1:3,2,0:NSideMortar,0:NSideMortar))
!        CALL ChangeBasis2D(3, NSideMortar, NSideMortar, N_Inter(NSideMortar)%Vdm_Leg, &
!                Mortar_Ja(1:3,3, 0:NSideMortar, 0:NSideMortar ,iMortar), Mortar_Ja_loc(1:3,3,0:NSideMortar,0:NSideMortar))
!        CALL SurfMetricsFromJa(NSideMortar,NormalDir,TangDir,NormalSign,Mortar_Ja_loc(1:3,3,0:NSideMortar,0:NSideMortar),&
!                               N_SurfMesh(SideID2)%NormVec(:,:,:),N_SurfMesh(SideID2)%TangVec1(:,:,:),&
!                               N_SurfMesh(SideID2)%TangVec2(:,:,:),N_SurfMesh(SideID2)%SurfElem(:,:))
      END IF
    END DO ! iMortar=1,4

  END IF ! MortarType(1,SideID).GT.0

END DO ! iLocSide=1,6

#if USE_HDG
!! Build SurfElemMin for all sides (including Mortar sides)
!DO iSide = 1, nSides
!  ! Get SurfElemMin
!  NSideMax = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
!  NSideMin = N_SurfMesh(iSide)%NSideMin
!  IF(NSideMax.EQ.NSideMin)THEN
!    N_SurfMesh(iSide)%SurfElemMin(:,:) = N_SurfMesh(iSide)%SurfElem(:,:)
!  ELSE
!    ! From high to low
!    ! Transform the slave side to the same degree as the master: switch to Legendre basis
!    CALL ChangeBasis2D(1, NSideMax, NSideMax, N_Inter(NSideMax)%sVdm_Leg, N_SurfMesh(iSide)%SurfElem(0:NSideMax,0:NSideMax), tmp(1,0:NSideMax,0:NSideMax))
!     !Switch back to nodal basis
!    CALL ChangeBasis2D(1, NSideMin, NSideMin, N_Inter(NSideMin)%Vdm_Leg , tmp(1,0:NSideMax,0:NSideMax)                      , N_SurfMesh(iSide)%SurfElemMin(0:NSideMin,0:NSideMin))
!  END IF ! NSideMax.EQ.NSideMin
!END DO ! iSide = 1, nSides
#endif /*USE_HDG*/

END SUBROUTINE CalcSurfMetrics


!==================================================================================================================================
!> Computes surface normal and tangential vectors and surface area from surface metrics Ja_Face.
!==================================================================================================================================
SUBROUTINE SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face,NormVec,TangVec1,TangVec2,SurfElem)
! MODULES
USE MOD_Globals       ,ONLY: CROSS,abort
#if USE_HDG
USE MOD_Symmetry_Vars ,ONLY: Symmetry
#endif /*USE_HDG*/
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
CHARACTER(32)      :: hilf
!==================================================================================================================================
WRITE(UNIT=hilf,FMT='(I0)') Nloc
DO q=0,Nloc; DO p=0,Nloc
  SurfElem(p,q) = SUM(Ja_Face(NormalDir,:,p,q)**2)
  IF(SurfElem(p,q).LT.0.)THEN
    CALL abort(__STAMP__,'Nloc='//TRIM(hilf)//': SurfElem(p,q).LT.0.',RealInfoOpt=SurfElem(p,q))
#if USE_HDG
  ELSEIF(Symmetry%Axisymmetric.AND.(SurfElem(p,q).EQ.0.))THEN
    NormVec( :,p,q) = (/0.,-1., 0./)
    TangVec1(:,p,q) = (/1., 0., 0./)
    TangVec2(:,p,q) = (/0., 0., 1./)
#endif /*USE_HDG*/
  ELSE
    IF(Nloc.GT.0.AND.ABS(SurfElem(p,q)).LE.0.0) CALL abort(__STAMP__,'Nloc='//TRIM(hilf)//&
      ': SUM(Ja_Face(NormalDir,:,p,q)**2) <= 0',RealInfoOpt=SurfElem(p,q))
    SurfElem(  p,q) = SQRT(SurfElem(p,q))
    NormVec( :,p,q) = NormalSign*Ja_Face(NormalDir,:,p,q)/SurfElem(p,q)
    TangVec1(:,p,q) = Ja_Face(TangDir,:,p,q) - SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q)) * NormVec(:,p,q)
    TangVec1(:,p,q) = SUM(TangVec1(:,p,q)**2)
    IF(ANY(ABS(TangVec1(:,p,q)).LE.0.0)) CALL abort(__STAMP__,'SUM(TangVec1(:,p,q)**2) <= 0')
    TangVec1(:,p,q) = TangVec1(:,p,q)/SQRT(TangVec1(:,p,q))
    TangVec2(:,p,q) = CROSS(NormVec(:,p,q),TangVec1(:,p,q))
  END IF ! SurfElem(p,q).LT.0.
END DO; END DO ! p,q
END SUBROUTINE SurfMetricsFromJa


END MODULE MOD_Metrics