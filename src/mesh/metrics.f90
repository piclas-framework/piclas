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

SUBROUTINE CalcMetrics()!XCL_NGeo)
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
USE MOD_Mesh_Vars, ONLY:NGeo,dXCL_NGeo,XCL_NGeo
USE MOD_Mesh_Vars, ONLY:Vdm_CLNGeo_GaussN,Vdm_CLNGeo_CLN,Vdm_CLN_GaussN
USE MOD_Mesh_Vars, ONLY:DCL_NGeo,DCL_N
USE MOD_Mesh_Vars, ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,crossProductMetrics
USE MOD_Mesh_Vars, ONLY:nElems
#ifdef PARTICLES
USE MOD_Particle_Surfaces,  ONLY:GetSuperSampledSurface
#endif /*PARTICLES*/
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_ChangeBasis,        ONLY:changeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,q,iElem
INTEGER            :: dd,ee,ff
INTEGER            :: nn,mm,ll
INTEGER            :: iGeo,jGeo,kGeo,lGeo
INTEGER            :: Cyclic(5),CycIJK(3),Cyc(3)
REAL               :: XCL_NGeo_loc(3,0:NGeo,0:NGeo,0:NGeo)    !mapping X(xi) P\in NGeo
!REAL               :: dXCL_NGeo(3,3,0:NGeo,0:NGeo,0:NGeo) !jacobi matrix of the mapping P\in NGeo
REAL               :: DetJacCL_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL               :: DetJacGauss_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL               :: XCL_N(3,0:PP_N,0:PP_N,0:PP_N)       ! mapping X(xi) P\in N
REAL               :: dXCL_N(3,3,0:PP_N,0:PP_N,0:PP_N)    !jacobi matrix interpolated on P\in N
REAL               :: R_CL_N(3,3,0:PP_N,0:PP_N,0:PP_N)    ! buffer for metric terms, uses XCL_N,dXCL_N 
REAL               :: JaCL_N(3,3,0:PP_N,0:PP_N,0:PP_N)    ! metric terms P\in N
REAL               :: scaledJac(2)
!===================================================================================================================================
ALLOCATE(dXCL_NGeo(3,3,0:NGeo,0:NGeo,0:NGeo,1:PP_nElems)) !jacobi matrix of the mapping P\in NGeo
! null outside!!
dXCL_NGeo=0.
! Prerequisites

Metrics_fTilde=0.
Metrics_gTilde=0.
Metrics_hTilde=0.
! 
Cyclic=(/1,2,3,1,2/)
! Outer loop over all elements
DO iElem=1,nElems
  XCL_NGeo_loc(:,:,:,:)=XCL_NGeo(:,:,:,:,iElem)
  !1.a) Jacobi Matrix of d/dxi_dd(X_nn): dXCL_NGeo(dd,nn,i,j,k,iElem)) 
  R_CL_N=0.
  dXCL_N=0.
  detJacCL_N=0.
  JaCL_N=0.
  DO nn=1,3
    DO dd=1,3
      DO kGeo=0,NGeo
                      CycIJK(3)=kGeo
        DO jGeo=0,NGeo
                        CycIJK(2)=jGeo
          DO iGeo=0,NGeo
                          CycIJK(1)=iGeo
            ! Matrix-vector multiplication
            Cyc=CycIJK
            DO lGeo=0,NGeo
              Cyc(dd)=lGeo   !d/dxi_dd 
              dXCL_NGeo(dd,nn,iGeo,jGeo,kGeo,iElem)=dXCL_NGeo(dd,nn,iGeo,jGeo,kGeo,iElem) + &
                                               DCL_NGeo(CycIJK(dd),lGeo)*XCL_NGeo_loc( nn ,Cyc(1),Cyc(2),Cyc(3))
            END DO !lGeo=0,NGeo
          END DO !iGeo=0,NGeo
        END DO !jGeo=0,NGeo
      END DO !kGeo=0,NGeo
    END DO !dd=1,3
  END DO !nn=1,3
  ! Interpolate the gradient of the mapping (living in a NGeo world) to the N (Cheb Lob) world
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,1,:,:,:,iElem),dXCL_N(:,1,:,:,:))
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,2,:,:,:,iElem),dXCL_N(:,2,:,:,:))
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:,iElem),dXCL_N(:,3,:,:,:))
  ! 1.b)Jacobians! grad(X_1) (grad(X_2) x grad(X_3))
  DO dd=1,3
    ee=Cyclic(dd+1) !cyclic!
    ff=Cyclic(dd+2) !cyclic!
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          detJacCL_N(1,i,j,k)=detJacCL_N(1,i,j,k) + &
                              dXCL_N(dd,1,i,j,k) * &
                              ( dXCL_N(ee,2,i,j,k)*dXCL_N(ff,3,i,j,k) - &
                                dXCL_N(ff,2,i,j,k)*dXCL_N(ee,3,i,j,k)   )

        END DO !i=0,N
      END DO !j=0,N
    END DO !k=0,N
  END DO !nn=1,3
  ! Interpolate the coordinates of the mapping (living in a NGeo world) to the N (Cheb Lob) world
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_CLN,XCL_NGeo_loc(:,:,:,:),XCL_N(:,:,:,:))
  !interpolate detJac to the GaussPoints
  CALL ChangeBasis3D(1,PP_N,PP_N,Vdm_CLN_GaussN,detJacCL_N(:,:,:,:),DetJacGauss_N(:,:,:,:))
  ! check scaled Jacobians
  scaledJac(1)=MINVAL(detJacCL_N(1,:,:,:))/MAXVAL(detJacCL_N(1,:,:,:))
  scaledJac(2)=MINVAL(detJacGauss_N(1,:,:,:))/MAXVAL(detJacGauss_N(1,:,:,:))
  IF(ANY(scaledJac.LT.0.01)) THEN
    WRITE(Unit_StdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac 
    CALL abort(__STAMP__,'Scaled Jacobian lower then tolerance!',iElem)
  END IF
  ! check for negative Jacobians
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        IF(detJacCL_N(1,i,j,k).LE.0.)THEN
          WRITE(Unit_StdOut,*) 'Negative Jacobian found in element on CL points. Coords:', XCL_NGeo_loc(:,i,j,k)
          WRITE(Unit_StdOut,*) 'Jacobian is:', detJacCL_N(1,i,j,k)
          CALL abort(__STAMP__,'Negative Jacobian found! Elem:',iElem)
        END IF
        IF(detJacGauss_N(1,i,j,k).LE.0.)THEN
          WRITE(Unit_StdOut,*) 'Negative Jacobian found in element on Gauss points. Coords:', XCL_N(:,i,j,k)
          WRITE(Unit_StdOut,*) 'Jacobian is:', detJacGauss_N(1,i,j,k)
          CALL abort(__STAMP__,'Negative Jacobian found! Elem:', iElem)
        END IF
      END DO !i=0,N
    END DO !j=0,N
  END DO !k=0,N

  ! assign to global Variable sJ
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        sJ(i,j,k,iElem)=1./DetJacGauss_N(1,i,j,k)
      END DO !iGeo=0,NGeo
    END DO !jGeo=0,NGeo
  END DO !kGeo=0,NGeo


  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO nn=1,3
      mm=Cyclic(nn+1) !cyclic!
      ll=Cyclic(nn+2) !cyclic!
      DO dd=1,3
        ee=Cyclic(dd+1) !cyclic!
        ff=Cyclic(dd+2) !cyclic!
        DO k=0,PP_N
          DO j=0,PP_N
            DO i=0,PP_N
              ! exact (cross-product) form
              JaCL_N(dd,nn,i,j,k)=dXCL_N(ee,mm,i,j,k)*dXCL_N(ff,ll,i,j,k) - dXCL_N(ee,ll,i,j,k)*dXCL_N(ff,mm,i,j,k)  
            END DO !i=0,N
          END DO !j=0,N
        END DO !k=0,N
      END DO !dd=1,3
    END DO !nn=1,3
  ELSE ! curl metrics 
    ! 2a.) Calculate X_l grad (X_m) 
    DO nn=1,3
      mm=Cyclic(nn+1) !cyclic!
      ll=Cyclic(nn+2) !cyclic!
      DO dd=1,3
        DO k=0,PP_N
          DO j=0,PP_N
            DO i=0,PP_N
              ! conservative curl form
              !R_CL_N(dd,nn,i,j,k)=XCL_N(ll,i,j,k)*dXCL_N(dd,mm,i,j,k) 
              ! invariant curl form
              R_CL_N(dd,nn,i,j,k)=0.5*(XCL_N(ll,i,j,k)*dXCL_N(dd,mm,i,j,k) - &
                                         XCL_N(mm,i,j,k)*dXCL_N(dd,ll,i,j,k)   )
            END DO !i=0,N
          END DO !j=0,N
        END DO !k=0,N
      END DO !dd=1,3
    END DO !nn=1,3
    DO nn=1,3
      DO dd=1,3
        ee=Cyclic(dd+1) !cyclic!
        ff=Cyclic(dd+2) !cyclic!
        DO k=0,PP_N
                  CycIJK(3)=k
          DO j=0,PP_N
                    CycIJK(2)=j
            DO i=0,PP_N
                      CycIJK(1)=i
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
            END DO !i=0,PP_N
          END DO !j=0,PP_N
        END DO !k=0,PP_N
      END DO!dd=1,3
    END DO!nn=1,3
    DO nn=1,3
      DO dd=1,3
        ee=Cyclic(dd+1) !cyclic!
        ff=Cyclic(dd+2) !cyclic!
        DO k=0,PP_N
                  CycIJK(3)=k
          DO j=0,PP_N
                    CycIJK(2)=j
            DO i=0,PP_N
                      CycIJK(1)=i
              !
              ! second part of the curl with cyclic indices
              !
              Cyc=CycIJK
              DO q=0,PP_N
                Cyc(ff)=q    !d/dxi_l
                JaCL_N(dd,nn,i,j,k)=JaCL_N(dd,nn,i,j,k) + &
                                    DCL_N(CycIJK(ff),q)*R_CL_N(ee,nn,Cyc(1),Cyc(2),Cyc(3))
              END DO!q=0,PP_N
            END DO !i=0,PP_N
          END DO !j=0,PP_N
        END DO !k=0,PP_N
      END DO!dd=1,3
    END DO!nn=1,3
  END IF !crossProductMetrics
  ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(:,:,:,:),Elem_xGP(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem))
  CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem))
  CALL CalcSurfMetrics(JaCL_N,XCL_N,iElem)
  CALL GetSuperSampledSurface(XCL_NGeo(:,:,:,:,iElem),iElem)
END DO !iElem=1,nElems
END SUBROUTINE CalcMetrics 



SUBROUTINE CalcSurfMetrics(JaCL_N,XCL_N,iElem)
!===================================================================================================================================
! Compute normal and tangential vectors from element metrics. Input is JaCL_N, the 3D element metrics on Cebychev-Lobatto points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,     ONLY:CROSS
USE MOD_Mesh_Vars,   ONLY:NGeo
USE MOD_Mesh_Vars,   ONLY:Vdm_CLN_GaussN
USE MOD_Mesh_Vars,   ONLY:ElemToSide,BCFace_xGP,nInnerSides,nBCSides,Face_xGP
USE MOD_Mesh_Vars,   ONLY:NormVec,TangVec1,TangVec2,SurfElem
USE MOD_Analyze_Vars,ONLY:CalcPoyntingInt
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_ChangeBasis, ONLY:ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
REAL,INTENT(IN)    :: JaCL_N(3,3,0:PP_N,0:PP_N,0:PP_N)  !Metrics from Element iElem
REAL,INTENT(IN)    :: XCL_N(3,0:PP_N,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
INTEGER            :: dd
INTEGER            :: sideID
REAL               :: Ja_Face(3,3,0:PP_N,0:PP_N)
REAL               :: tmp(3,0:PP_N,0:PP_N)
!===================================================================================================================================

! interpolate to xi sides
IF(ElemToSide(E2S_FLIP,XI_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,0,:,:),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        BCFace_xGP(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,0,:,:),tmp)
      ! turn into right hand system of side
      DO q=0,PP_N
        DO p=0,PP_N
          Face_xGP(1:3,p,q,sideID)=tmp(:,q,p)
        END DO !p
      END DO !q
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,0,:,:),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        Ja_Face(dd,1:3,p,q)=tmp(:,q,p)
      END DO !p
    END DO !q
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(1,:,p,q)**2))
      NormVec( :,p,q,SideID) = -Ja_Face(1,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(2,:,p,q)-SUM(Ja_Face(2,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0

IF(ElemToSide(E2S_FLIP,XI_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,PP_N,:,:),BCFace_xGP(1:3,:,:,sideID))
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF ! CalcPoyntintIntegral
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,PP_N,:,:),Face_xGP(1:3,:,:,SideID))
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,PP_N,:,:),Ja_Face(dd,1:3,:,:))
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(1,:,p,q)**2))
      NormVec( :,p,q,SideID) = Ja_Face(1,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(2,:,p,q)-SUM(Ja_Face(2,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0

! interpolate to eta sides
IF(ElemToSide(E2S_FLIP,ETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,0,:),BCFace_xGP(1:3,:,:,sideID))
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF ! CalcPoyntingInt
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,0,:),Face_xGP(1:3,:,:,SideID))
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,:,0,:),Ja_Face(dd,1:3,:,:))
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(2,:,p,q)**2))
      NormVec( :,p,q,SideID) = -Ja_Face(2,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(3,:,p,q)-SUM(Ja_Face(3,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0
  
IF(ElemToSide(E2S_FLIP,ETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,PP_N,:),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        BCFace_xGP(1:3,p,q,sideID)=tmp(:,PP_N-p,q)
      END DO !p
    END DO !q
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,PP_N,:),tmp)
      ! turn into right hand system of side
      DO q=0,PP_N
        DO p=0,PP_N
          Face_xGP(1:3,p,q,sideID)=tmp(:,PP_N-p,q)
        END DO !p
      END DO !q
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,:,PP_N,:),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        Ja_Face(dd,1:3,p,q)=tmp(:,PP_N-p,q)
      END DO !p
    END DO !q
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(2,:,p,q)**2))
      NormVec( :,p,q,SideID) = Ja_Face(2,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(3,:,p,q)-SUM(Ja_Face(3,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0

! interpolate to zeta sides
IF(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,:,0),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        BCFace_xGP(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,:,0),tmp)
      ! turn into right hand system of side
      DO q=0,PP_N
        DO p=0,PP_N
          Face_xGP(1:3,p,q,sideID)=tmp(:,q,p)
        END DO !p
      END DO !q
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,:,:,0),tmp)
    ! turn into right hand system of side
    DO q=0,PP_N
      DO p=0,PP_N
        Ja_Face(dd,1:3,p,q)=tmp(:,q,p)
      END DO !p
    END DO !q
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(3,:,p,q)**2))
      NormVec( :,p,q,SideID) = -Ja_Face(3,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(1,:,p,q)-SUM(Ja_Face(1,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0

IF(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
  IF ((sideID.LE.nBCSides))THEN !BC
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,:,PP_N),BCFace_xGP(1:3,:,:,sideID))
    IF (CalcPoyntingInt) THEN
      Face_xGP(:,:,:,SideID) = BCFace_xGP(:,:,:,SideID)
    END IF ! PoyntingIntegral
  END IF !BC
  IF (CalcPoyntingInt) THEN
    IF((SideID.GT.nBCSides))THEN ! > BC
       CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,XCL_N(1:3,:,:,PP_N),Face_xGP(1:3,:,:,SideID))
    END IF ! BC
  END IF ! PoyntingIntegral
  DO dd=1,3
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_GaussN,JaCL_N(dd,1:3,:,:,PP_N),Ja_Face(dd,1:3,:,:))
  END DO
  DO q=0,PP_N
    DO p=0,PP_N
      SurfElem(  p,q,SideID) = SQRT(SUM(Ja_Face(3,:,p,q)**2))
      NormVec( :,p,q,SideID) = Ja_Face(3,:,p,q)/SurfElem(p,q,SideID)
      TangVec1(:,p,q,SideID) = Ja_Face(1,:,p,q)-SUM(Ja_Face(1,:,p,q)*NormVec(:,p,q,SideID))*NormVec(:,p,q,SideID)
      TangVec1(:,p,q,SideID) = TangVec1(:,p,q,SideID)/SQRT(SUM(TangVec1(:,p,q,SideID)**2))
      TangVec2(:,p,q,SideID) = CROSS(NormVec(:,p,q,SideID),TangVec1(:,p,q,SideID))
    END DO
  END DO
END IF !flip=0

END SUBROUTINE CalcSurfMetrics

END MODULE MOD_Metrics
