#include "boltzplatz.h"

MODULE MOD_Jac_Ex
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitJac_Ex
  MODULE PROCEDURE InitJac_Ex
END INTERFACE

INTERFACE Jac_Ex
  MODULE PROCEDURE Jac_Ex
END INTERFACE

INTERFACE Jac_Ex_Neighbor
  MODULE PROCEDURE Jac_Ex_Neighbor
END INTERFACE

INTERFACE Jac_Ex1D
  MODULE PROCEDURE Jac_Ex1D
END INTERFACE


INTERFACE FinalizeJac_Ex
  MODULE PROCEDURE FinalizeJac_Ex
END INTERFACE


PUBLIC::InitJac_Ex,Jac_Ex,FinalizeJac_Ex,Jac_Ex_Neighbor,Jac_Ex1D
!===================================================================================================================================

CONTAINS

SUBROUTINE InitJac_Ex()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Jac_Ex_Vars
USE MOD_Interpolation_Vars ,ONLY:L_minus,L_plus
USE MOD_DG_Vars            ,ONLY:L_Hatminus,L_Hatplus
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                      :: i,j
!===================================================================================================================================

IF(Jac_Ex_InitIsDone) RETURN

!SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT BLOCK JACOBIAN...'
ALLOCATE(LL_minus(0:PP_N,0:PP_N), LL_plus(0:PP_N,0:PP_N))

DO j=0,PP_N
  DO i=0,PP_N
  LL_minus(i,j) = L_Hatminus(i)*L_minus(j)
  LL_plus(i,j)  = L_Hatplus(i) *L_plus(j)
  END DO
END DO 

Jac_Ex_InitIsDone=.TRUE.
!SWRITE(UNIT_stdOut,'(A)')' EXACT BLOCK JACOBIAN DONE!'
END SUBROUTINE InitJac_Ex


SUBROUTINE Jac_Ex(iElem,BJ)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_LinearSolver_Vars ,ONLY:nDOFElem
USE MOD_JacSurfInt        ,ONLY:DGJacSurfInt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(1:nDOFelem,1:nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! element local preconditioner
CALL DGVolIntJac(BJ,iElem)  ! without sJ!
CALL DGJacSurfInt(BJ,iElem) ! without sJ!
CALL Apply_sJ(BJ,iElem)

END SUBROUTINE Jac_Ex

SUBROUTINE Jac_Ex_Neighbor(locSideID,iElem,BJ)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_LinearSolver_Vars ,ONLY:nDOFElem
USE MOD_JacSurfInt    ,ONLY:DGJacSurfInt_Neighbor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(1:nDOFelem,1:nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! element local preconditioner
CALL DGJacSurfInt_Neighbor(BJ,locSideID,iElem) 
!CALL Apply_sJ(BJ,iElem)

END SUBROUTINE Jac_Ex_Neighbor

SUBROUTINE  DGVolIntJac(BJ,iElem)
!===================================================================================================================================
! volume integral Jacobian of the euler flux and the viscous part dF^v/dU
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars           ,ONLY: D_hat
USE MOD_Mesh_Vars         ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_LinearSolver_Vars ,ONLY: nDOFelem
USE MOD_Jacobian          ,ONLY: EvalFluxJacobian, EvalFluxJacobianDielectric
USE MOD_Dielectric_Vars,   ONLY: DoDielectric,isDielectricElem, ElemToDielectric,DielectricConstant_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: BJ(1:nDOFelem,1:nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                                 :: i,mm,nn,oo
INTEGER                                                 :: s,r1,r2,r3,ll,vn1,vn2
REAL                                                    :: delta(0:PP_N,0:PP_N)
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJac,gJac,hJac
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacTilde,gJacTilde,hJacTilde
LOGICAL                                                 :: isDielectric
!===================================================================================================================================
  vn1=PP_nVar*(PP_N+1)
  vn2=vn1*(PP_N+1)
  delta=0.
  DO i=0,PP_N
    delta(i,i)=1.
  END DO

  ! call this always, independent if dielectric or not
  CALL EvalFluxJacobian(fJac,gJac,hJac) !fills fJac,gJac,hJac

  isDielectric=.FALSE.
  IF(DoDielectric)THEN
    IF(isDielectricElem(iElem)) THEN ! 1.) select dielectric elements
      isDielectric=.TRUE.
    END IF
  END IF


  !VERSION 00
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        fJacTilde=fJac(:,:,i,j,k)
!        gJacTilde=gJac(:,:,i,j,k)
!        hJacTilde=hJac(:,:,i,j,k)
!        ! Compute the transformed fluxes with the metric terms
!        ! Attention 1: we store the transformed fluxes in f,g,h again
!        fJac(:,:,i,j,k) = fJacTilde(:,:)*Metrics_fTilde(1,i,j,k,iElem) + &
!                          gJacTilde(:,:)*Metrics_fTilde(2,i,j,k,iElem) + &
!                          hJacTilde(:,:)*Metrics_fTilde(3,i,j,k,iElem) 
!        gJac(:,:,i,j,k) = fJacTilde(:,:)*Metrics_gTilde(1,i,j,k,iElem) + &
!                          gJacTilde(:,:)*Metrics_gTilde(2,i,j,k,iElem) + &
!                          hJacTilde(:,:)*Metrics_gTilde(3,i,j,k,iElem)
!        hJac(:,:,i,j,k) = fJacTilde(:,:)*Metrics_hTilde(1,i,j,k,iElem) + &
!                          gJacTilde(:,:)*Metrics_hTilde(2,i,j,k,iElem) + &
!                          hJacTilde(:,:)*Metrics_hTilde(3,i,j,k,iElem)
!      END DO ! i
!    END DO ! j
!  END DO ! k
!  s=0
!  DO oo=0,PP_N
!    DO nn=0,PP_N
!      DO mm=0,PP_N
!        r=0
!        DO k=0,PP_N
!          DO j=0,PP_N
!            DO i=0,PP_N
!             BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)                       &
!                                                                 + D_hat(i,mm)*delta(j,nn)*delta(k,oo)*fJac(:,:,mm,j,k) &
!                                                                 + delta(i,mm)*D_hat(j,nn)*delta(k,oo)*gJac(:,:,i,nn,k) &
!                                                                 + delta(i,mm)*delta(j,nn)*D_hat(k,oo)*hJac(:,:,i,j,oo)
!              r=r+PP_nVar
!            END DO !i
!          END DO !j
!        END DO !k 
!        s=s+PP_nVar
!      END DO !mm
!    END DO !nn
!  END DO !oo 

  IF(isDielectric)THEN
    s=0
    DO oo=0,PP_N
      DO nn=0,PP_N
        DO mm=0,PP_N
          !fills fJac,gJac,hJac
          CALL EvalFluxJacobianDielectric(DielectricConstant_Inv(mm,nn,oo,ElemToDielectric(iElem)),fJac,gJac,hJac) 
          fJacTilde(:,:) = fJac(:,:)*Metrics_fTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_fTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_fTilde(3,mm,nn,oo,iElem) 
          gJacTilde(:,:) = fJac(:,:)*Metrics_gTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_gTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_gTilde(3,mm,nn,oo,iElem)
          hJacTilde(:,:) = fJac(:,:)*Metrics_hTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_hTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_hTilde(3,mm,nn,oo,iElem)
          r1=           vn1*nn+vn2*oo
          r2=mm*PP_nVar      +vn2*oo
          r3=mm*PP_nVar+vn1*nn
          DO ll=0,PP_N
                BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar)                &
                                                                                  + D_hat(ll,mm)*fJacTilde(:,:) 
                BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar)&
                                                                                  + D_hat(ll,nn)*gJacTilde(:,:) 
                BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar)&
                                                                                  + D_hat(ll,oo)*hJacTilde(:,:) 
                r1=r1+PP_nVar
                r2=r2+vn1
                r3=r3+vn2
          END DO !k 
          s=s+PP_nVar
        END DO !mm
      END DO !nn
    END DO !oo 
  ELSE
    s=0
    DO oo=0,PP_N
      DO nn=0,PP_N
        DO mm=0,PP_N
          fJacTilde(:,:) = fJac(:,:)*Metrics_fTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_fTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_fTilde(3,mm,nn,oo,iElem) 
          gJacTilde(:,:) = fJac(:,:)*Metrics_gTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_gTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_gTilde(3,mm,nn,oo,iElem)
          hJacTilde(:,:) = fJac(:,:)*Metrics_hTilde(1,mm,nn,oo,iElem) + &
                           gJac(:,:)*Metrics_hTilde(2,mm,nn,oo,iElem) + &
                           hJac(:,:)*Metrics_hTilde(3,mm,nn,oo,iElem)
          r1=           vn1*nn+vn2*oo
          r2=mm*PP_nVar      +vn2*oo
          r3=mm*PP_nVar+vn1*nn
          DO ll=0,PP_N
                BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar) = BJ(r1+1:r1+PP_nVar,s+1:s+PP_nVar)                &
                                                                                  + D_hat(ll,mm)*fJacTilde(:,:) 
                BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar) = BJ(r2+1:r2+PP_nVar,s+1:s+PP_nVar)&
                                                                                  + D_hat(ll,nn)*gJacTilde(:,:) 
                BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar) = BJ(r3+1:r3+PP_nVar,s+1:s+PP_nVar)&
                                                                                  + D_hat(ll,oo)*hJacTilde(:,:) 
                r1=r1+PP_nVar
                r2=r2+vn1
                r3=r3+vn2
          END DO !k 
          s=s+PP_nVar
        END DO !mm
      END DO !nn
    END DO !oo 
  END IF

END SUBROUTINE DGVolIntJac

SUBROUTINE  DGVolIntJac1D(dRdXi,dRdEta,dRdZeta,iElem)
!===================================================================================================================================
! volume integral Jacobian of the euler flux and the viscous part dF^v/dU
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY: D_hat
USE MOD_Mesh_Vars     ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_LinearSolver_Vars ,ONLY: nDOFLine
USE MOD_Jacobian          ,ONLY: EvalFluxJacobian, EvalFluxJacobianDielectric
USE MOD_Dielectric_Vars,   ONLY: DoDielectric,isDielectricElem, ElemToDielectric,DielectricConstant_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: dRdXi  (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdEta (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdZeta(1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                                 :: i,j,k,l,vni,vnj,vnk,s
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJac,gJac,hJac
REAL,DIMENSION(PP_nVar,PP_nVar)                         :: fJacTilde,gJacTilde,hJacTilde
LOGICAL                                                 :: isDielectric
!===================================================================================================================================

  CALL EvalFluxJacobian(fJac,gJac,hJac) !fills fJac,gJac,hJac

  isDielectric=.FALSE.
  IF(DoDielectric)THEN
    IF(isDielectricElem(iElem)) THEN ! 1.) select dielectric elements
      isDielectric=.TRUE.
    END IF
  END IF

  IF(isDielectric)THEN
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          !fills fJac,gJac,hJac
          CALL EvalFluxJacobianDielectric(DielectricConstant_Inv(i,j,k,ElemToDielectric(iElem)),fJac,gJac,hJac) 

          fJacTilde(:,:) = fJac(:,:)*Metrics_fTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_fTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_fTilde(3,i,j,k,iElem) 
          gJacTilde(:,:) = fJac(:,:)*Metrics_gTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_gTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_gTilde(3,i,j,k,iElem)
          hJacTilde(:,:) = fJac(:,:)*Metrics_hTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_hTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_hTilde(3,i,j,k,iElem)
          vni=PP_nVar*i
          vnj=PP_nVar*j
          vnk=PP_nVar*k
          s=0
          DO l=0,PP_N
            !dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,j,k) = dRdXi  (vni+1:vni+1,s+1:s+PP_nVar,j,k) +D_hat(l,i)*fJacTilde(:,:) 
            !dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,i,k) = dRdEta (vnj+1:vnj+1,s+1:s+PP_nVar,i,k) +D_hat(l,j)*gJacTilde(:,:) 
            !dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,i,j) = dRdZeta(vnk+1:vnk+1,s+1:s+PP_nVar,i,j) +D_hat(l,k)*hJacTilde(:,:) 
            dRdXi  (s+1:s+PP_nVar,vni+1:vni+PP_nVar,j,k) = dRdXi  (s+1:s+PP_nVar,vni+1:vni+PP_nVar,j,k) +D_hat(l,i)*fJacTilde(:,:) 
            dRdEta (s+1:s+PP_nVar,vnj+1:vnj+PP_nVar,i,k) = dRdEta (s+1:s+PP_nVar,vnj+1:vnj+PP_nVar,i,k) +D_hat(l,j)*gJacTilde(:,:) 
            dRdZeta(s+1:s+PP_nVar,vnk+1:vnk+PP_nVar,i,j) = dRdZeta(s+1:s+PP_nVar,vnk+1:vnk+PP_nVar,i,j) +D_hat(l,k)*hJacTilde(:,:) 
            s=s+PP_nVar
          END DO ! l
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          fJacTilde(:,:) = fJac(:,:)*Metrics_fTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_fTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_fTilde(3,i,j,k,iElem) 
          gJacTilde(:,:) = fJac(:,:)*Metrics_gTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_gTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_gTilde(3,i,j,k,iElem)
          hJacTilde(:,:) = fJac(:,:)*Metrics_hTilde(1,i,j,k,iElem) + &
                           gJac(:,:)*Metrics_hTilde(2,i,j,k,iElem) + &
                           hJac(:,:)*Metrics_hTilde(3,i,j,k,iElem)
          vni=PP_nVar*i
          vnj=PP_nVar*j
          vnk=PP_nVar*k
          s=0
          DO l=0,PP_N
            !dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,j,k) = dRdXi  (vni+1:vni+1,s+1:s+PP_nVar,j,k) +D_hat(l,i)*fJacTilde(:,:) 
            !dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,i,k) = dRdEta (vnj+1:vnj+1,s+1:s+PP_nVar,i,k) +D_hat(l,j)*gJacTilde(:,:) 
            !dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,i,j) = dRdZeta(vnk+1:vnk+1,s+1:s+PP_nVar,i,j) +D_hat(l,k)*hJacTilde(:,:) 
            dRdXi  (s+1:s+PP_nVar,vni+1:vni+PP_nVar,j,k) = dRdXi  (s+1:s+PP_nVar,vni+1:vni+PP_nVar,j,k) +D_hat(l,i)*fJacTilde(:,:) 
            dRdEta (s+1:s+PP_nVar,vnj+1:vnj+PP_nVar,i,k) = dRdEta (s+1:s+PP_nVar,vnj+1:vnj+PP_nVar,i,k) +D_hat(l,j)*gJacTilde(:,:) 
            dRdZeta(s+1:s+PP_nVar,vnk+1:vnk+PP_nVar,i,j) = dRdZeta(s+1:s+PP_nVar,vnk+1:vnk+PP_nVar,i,j) +D_hat(l,k)*hJacTilde(:,:) 
            s=s+PP_nVar
          END DO ! l
        END DO ! i
      END DO ! j
    END DO ! k
  END IF

END SUBROUTINE DGVolIntJac1D

SUBROUTINE Apply_sJ(BJ,iElem)
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars ,ONLY:nDOFElem
USE MOD_Mesh_Vars     ,ONLY:sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(nDOFelem,nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: r,s,i,j,k
!===================================================================================================================================

DO s=0,nDOFelem-1,PP_nVar
  r=0
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = -sJ(i,j,k,iElem)*BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)
        r=r+PP_nVar
      END DO !i
    END DO !j
  END DO !k
END DO ! s

END SUBROUTINE Apply_sJ

SUBROUTINE Apply_sJ1D(dRdXi,dRdEta,dRdZeta,iElem)
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine,mass
USE MOD_Mesh_Vars         ,ONLY:sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: dRdXi  (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdEta (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdZeta(1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: p,q,l,s
!===================================================================================================================================

DO p=0,PP_N
  DO q=0,PP_N
      DO l=0,PP_N
      s=l*PP_nVar
      dRdXi  (s+1:s+PP_nVar,:,p,q) = -sJ(l,p,q,iElem)*mass(1,l,p,q,iElem)*dRdXi  (s+1:s+PP_nVar,:,p,q)
      dRdEta (s+1:s+PP_nVar,:,p,q) = -sJ(p,l,q,iElem)*mass(1,p,l,q,iElem)*dRdEta (s+1:s+PP_nVar,:,p,q)
      dRdZeta(s+1:s+PP_nVar,:,p,q) = -sJ(p,q,l,iElem)*mass(1,p,q,l,iElem)*dRdZeta(s+1:s+PP_nVar,:,p,q)
    END DO !l
  END DO !p
END DO ! q

END SUBROUTINE Apply_sJ1D

SUBROUTINE Jac_Ex1D(dRdXi,dRdEta,dRdZeta,iElem)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_JacSurfInt        ,ONLY:DGJacSurfInt1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: dRdXi  (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdEta (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT) :: dRdZeta(1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! element local preconditioner
CALL DGVolIntJac1D(dRdXi,dRdEta,dRdZeta,iElem)  ! without sJ!
CALL DGJacSurfInt1D(dRdXi,dRdEta,dRdZeta,iElem) ! without sJ!
CALL Apply_sJ1D(dRdXi,dRdEta,dRdZeta,iElem)

END SUBROUTINE Jac_Ex1D


SUBROUTINE FinalizeJac_Ex()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
!SDEALLOCATE()
END SUBROUTINE FinalizeJac_Ex

END MODULE MOD_Jac_Ex
