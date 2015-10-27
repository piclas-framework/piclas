#include "boltzplatz.h"

MODULE MOD_ApplyPreconditioner
!===================================================================================================================================
! Module for the Block-Jacobi Preconditioner  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE Preconditioner
  MODULE PROCEDURE Preconditioner
END INTERFACE

!PUBLIC :: ApplyPrecond_Elem
!PUBLIC :: ApplyPrecond_DOF
!PUBLIC :: ApplyPrecond_Jacobi
!PUBLIC :: ApplyILU
!PUBLIC :: ApplyPrecond_GMRES
PUBLIC :: Preconditioner
!===================================================================================================================================

CONTAINS

SUBROUTINE Preconditioner(coeff,V,Vprecond)
!===================================================================================================================================
! select preconditioner
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars,           ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)  :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL,INTENT(INOUT)  :: Vprecond(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL                :: coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: i
!===================================================================================================================================

SELECT CASE(PrecondType)
  CASE(1,2) ! Finite-Differences and analytical Preconditioner
    ! BJ per Element
    CALL ApplyPrecond_Elem(V,Vprecond)
  CASE(3) 
    ! BJ per DOF, only analytic
    CALL ApplyPrecond_DOF(V,Vprecond)
  CASE(4)
    !CALL abort(__STAMP__,'Preconditioner not implemented!',PrecondType,999.)
    ! iterative through 
    CALL ApplyPrecond_GMRES(coeff,V,Vprecond)
  CASE(5,6,7,8,9,10,11)
    ! Jaboci Preconditioner 
    CALL ApplyPrecond_Jacobi(V,Vprecond)
  CASE(22)
    CALL ApplyILU(V,Vprecond)
  CASE(25)
    CALL ApplyBILU0(V,Vprecond)
  CASE(26)
    CALL ApplyBILU0_BCSR(V,Vprecond)
  CASE(200,205) ! - 205 addaptive schwarz preconditioner
    CALL ApplyTensorBJ(V,Vprecond)
  CASE(203)
    CALL ApplyTensorBJopt(V,Vprecond)
  CASE(201,202)
    CALL ApplyTensorProductBJ(V,Vprecond)
  CASE(210,212) ! multiplicative schwarz preconditioner
    CALL ApplyMultiSchwarz(V,Vprecond)
  CASE(211) ! multiplicative schwarz preconditioner
    CALL ApplyMultiSchwarzDG(coeff,V,Vprecond)
  CASE DEFAULT
    Vprecond=V
END SELECT

END SUBROUTINE Preconditioner

SUBROUTINE ApplyPrecond_Elem(v,z)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
USE MOD_Precond_Vars      ,ONLY:invP,PrecondMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: v(nDOFelem,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: z(nDOFelem,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: r,s,iElem
!===================================================================================================================================

SELECT CASE(PrecondMethod)
  CASE(0) ! using a loop
    ! loop over all elements
    DO iElem=1,PP_nElems
      z(:,iElem) = 0.
      DO s=1,nDOFElem
        DO r=1,nDOFElem
          z(r,iElem) = z(r,iElem)+invP(r,s,iElem)*v(s,iElem)
        END DO ! r
      END DO ! s
    END DO ! iElem
  CASE(1) ! Matmul
    DO iElem=1,PP_nElems
      z(:,iElem) = MATMUL(invP(:,:,iElem),v(:,iElem))
    END DO !iElem
END SELECT

END SUBROUTINE  ApplyPrecond_Elem

SUBROUTINE ApplyPrecond_DOF(v,z)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per DOF
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
USE MOD_Precond_Vars      ,ONLY:invBJ,PrecondMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: v(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: z(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k,iVar,jVar
!===================================================================================================================================

SELECT CASE(PrecondMethod)
  CASE(0) ! using a loop
    ! loop over all DOFS
    z=0.
    DO iElem=1,PP_nElems
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO jVar=1,PP_nVar
              DO iVar=1,PP_nVar
                z(iVar,i,j,k,iElem) = z(iVar,i,j,k,iElem) + invBJ(iVar,jVar,i,j,k,iElem)*v(jVar,i,j,k,iElem)
              END DO ! iVar
            END DO ! jVar
          END DO ! i
        END DO ! j
      END DO ! k
    END DO ! Elem
  CASE(1) ! Matmul
    DO iElem=1,PP_nElems
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            z(:,i,j,k,iElem) = MATMUL(invBJ(:,:,i,j,k,iElem),v(:,i,j,k,iElem))
          END DO ! i
        END DO ! j
      END DO ! k
    END DO !iElem
END SELECT

END SUBROUTINE  ApplyPrecond_DOF

SUBROUTINE ApplyPrecond_Jacobi(v,z)
!===================================================================================================================================
! Apply Jacobi Preconditioner which is constructed per DOF
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars  ,ONLY:invJ,PrecondMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: v(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: z(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k,iVar
!===================================================================================================================================

 ! loop over all DOFS
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          z(iVar,i,j,k,iElem) =  invJ(iVar,i,j,k,iElem)*v(iVar,i,j,k,iElem)
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! Elem

END SUBROUTINE  ApplyPrecond_Jacobi

SUBROUTINE ApplyPrecond_GMRES(coeff,Vin,Vout)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess   ! why not Un??
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_JacDG,                ONLY: JacDG
USE MOD_LinearSolver_Vars,    ONLY: eps_LinearSolver
USE MOD_LinearSolver_Vars,    ONLY: nRestarts
USE MOD_LinearOperator,       ONLY: ElementVectorDotProduct
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: coeff, Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: Un, W,Z,R0
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: V
REAL,ALLOCATABLE                        :: Gam(:),C(:),S(:),H(:,:),Alp(:)
REAL                                    :: Norm_R0,Resu,Temp,Bet
REAL                                    :: AbortCrit
INTEGER                                 :: Restart
INTEGER                                 :: m,nn,o,IterPrecond
!REAL                                    :: tS,tE
INTEGER                                 :: nKDimPrecond
INTEGER                                 :: iElem
!===================================================================================================================================

nKDimPrecond=3

ALLOCATE( Gam(1:nKDimPrecond+1)           &
        , C  (1:nKDimPrecond)             &
        , S  (1:nKDimPrecond)             &
        , H  (1:nKDimPrecond+1,1:nKDIMPrecond+1) &
        , Alp(1:nKDimPrecond)             )

ALLOCATE( Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
        , W (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
        , Z (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
        , R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
        , V (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nKDimPrecond)   )

! time measurement
!CALL CPU_TIME(tS)
DO iElem=1,PP_nElems
  ! start GMRES
  Restart=0
  Un=Vin(:,:,:,:,iElem)
  IterPrecond=0
!  DO WHILE (Restart<nRestarts)
    ! start residuum berrechnen
    CALL JacDG(coeff,iElem,Un,R0)
    R0=Un-R0
    CALL ElementVectorDotProduct(R0,R0,Norm_R0)
    Norm_R0=SQRT(Norm_R0)
    !AbortCrit=Norm_R0*eps_LinearSolver
    AbortCrit=Norm_R0*1e-3
    ! GMRES(m)  inner loop
    V(:,:,:,:,1)=R0/Norm_R0
    Gam(1)=Norm_R0
    DO m=1,nKDimPrecond
      IterPrecond=IterPrecond+1
      ! matrix vector
      CALL JacDG(coeff,iElem,V(:,:,:,:,m),W)
      ! Gram-Schmidt
      DO nn=1,m
        CALL ElementVectorDotProduct(V(:,:,:,:,nn),W,H(nn,m))
        W=W-H(nn,m)*V(:,:,:,:,nn)
      END DO !nn
      CALL ElementVectorDotProduct(W,W,Resu)
      H(m+1,m)=SQRT(Resu)
      ! Givens Rotation
      DO nn=1,m-1
        Temp     =   C(nn)*H(nn,m) + S(nn)*H(nn+1,m)
        H(nn+1,m) = - S(nn)*H(nn,m) + C(nn)*H(nn+1,m)
        H(nn,m)   =   Temp
      END DO !nn
      Bet=SQRT(H(m,m)*H(m,m)+H(m+1,m)*H(m+1,m))
      S(m)=H(m+1,m)/Bet
      C(m)=H(m,m)/Bet 
      H(m,m)=Bet
      Gam(m+1)=-S(m)*Gam(m)
      Gam(m)=C(m)*Gam(m)
      IF ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDimPrecond)) THEN !converge or max Krylov reached
        DO nn=m,1,-1
           Alp(nn)=Gam(nn) 
           DO o=nn+1,m
             Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
           END DO !o
           Alp(nn)=Alp(nn)/H(nn,nn)
        END DO !nn
        DO nn=1,m
          Un=Un+Alp(nn)*V(:,:,:,:,nn)
        END DO !nn
        IF ((ABS(Gam(m+1)).LE.AbortCrit).OR.(M.EQ.nKDimPrecond)) THEN !converged
          Vout(:,:,:,:,iElem)=Un
          
  !        CALL CPU_TIME(tE)
  !        ! Debug Ausgabe, Anzahl der Iterationen...
  !        SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter Precond       : ',IterPrecond
  !        SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in P-GMRES    : ',tE-tS
  !        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Gam(1)
  !        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Gam(m+1)
          !SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R/Norm_R0       : ',Gam(m+1)/Gam(1)
          EXIT
        END IF  ! converged
      ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
        V(:,:,:,:,m+1)=W/H(m+1,m)
      END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
    END DO ! m 
    ! Restart needed
!    Restart=Restart+1
!  END DO ! Restart

END DO ! iElem
DEALLOCATE(Un,W,R0,Z,V)
DEALLOCATE(Gam,C,S,H,Alp)

RETURN

CALL abort(__STAMP__, &
     'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(IterPrecond))

END SUBROUTINE ApplyPrecond_GMRES


!SUBROUTINE Apply_J(BJ,iElem)
!===================================================================================================================================
!! Deallocate global variables
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Implicit_Vars ,ONLY:nDOFElem
!USE MOD_Mesh_Vars     ,ONLY:sJ
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: BJ(nDOFelem,nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!INTEGER :: r,s,i,j,k
!===================================================================================================================================
!DO s=0,nDOFelem-1,PP_nVar
!  r=0
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        !BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = -sJ(i,j,k,iElem)*BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)
!        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)/sJ(i,j,k,iElem)
!        r=r+PP_nVar
!      END DO !i
!    END DO !j
!  END DO !k
!END DO ! s
!END SUBROUTINE Apply_J
!
!SUBROUTINE Add_J(BJ,iElem)
!===================================================================================================================================
!! Deallocate global variables
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Implicit_Vars ,ONLY:nDOFElem
!USE MOD_Mesh_Vars     ,ONLY:sJ
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: BJ(nDOFelem,nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!INTEGER :: r,s,i,j,k,vn1,vn2
!===================================================================================================================================
!
!vn1 = PP_nVar * (PP_N + 1)
!vn2 = vn1 * (PP_N +1)
!
!
!DO k = 0,PP_N
!   DO j = 0,PP_N
!     DO i = 0,PP_N
!       r = vn2 * k + vn1 * j + PP_nVar * i
!       s = r
!       BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) +1./sJ(i,j,k,iElem)
!   END DO ! i
!  END DO ! j
!END DO ! k
!
!END SUBROUTINE Add_J

SUBROUTINE ApplyILU(Vin,Vout)
!==================================================================================================================================
! Application of block ILU0 preconditioner
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY: nDOFelem
USE MOD_CSR_Vars               ,ONLY: DE,IL,IU
USE MOD_CSR_Vars               ,ONLY: nUNonZeros,nLNonZeros,nMTriangle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                    :: Vin(1:nDOFElem,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                      :: Vout(1:nDOFElem,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                               :: Vcalc(1:nDOFElem,1:PP_nElems)
INTEGER                                               :: ii,k1,k2,jj,iEntry,iElem
REAL                                                  :: tS,tE
logical::test
!==================================================================================================================================

!test=.FALSE.
!DO iElem=1,PP_nElems
!  DO ii=1,nDOFElem
!    IF(ISNAN(Vin(ii,iElem)))THEN
!      !WRITE(*,*) 'NAN. ii,iElem',ii,iElem
!      TEST=.TRUE.
!    ENDIF
!  END DO
!END DO
!IF(TEST) THEN
!  print*,'buggy init'
!  read*
!END IF

Vcalc=Vin
DO iElem=1,PP_nElems
  DO ii=1,nMTriangle,1
    iEntry=ii+1
    k1=IL(iElem)%iEntry(ii)
    k2=IL(iElem)%iEntry(ii+1)-1
    DO jj=k1,k2
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IL(iElem)%Entry(jj)*Vcalc(IL(iELEM)%jEntry(jj),iElem)
    END DO ! jj
  END DO ! ii
  ! backward elimination
  ! init backward Gauss
  Vout(nDOFElem,iElem) = Vcalc(nDOFElem,iElem)/DE(nDOFElem,iElem)
  DO ii=nMTriangle,1,-1
    k1=IU(iElem)%iEntry(ii)
    k2=IU(iElem)%iEntry(ii+1)-1
    iEntry=ii
    DO jj=k1,k2
      !Vin(iEntry)=Vin(iEntry)-UE(jj)*Vout(JU(jj))
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IU(iElem)%Entry(jj)*Vout(IU(iElem)%jEntry(jj),iElem)
    END DO ! jj
    Vout(iEntry,iElem)=Vcalc(iEntry,iElem)/DE(iEntry,iElem)
  END DO ! ii
END DO


!DO iElem=1,PP_nElems
!  DO ii=1,nMTriangle(iElem),2
!    k1=IL(iElem)%iEntry(ii)
!    k2=IL(iElem)%iEntry(ii+1)
!    iEntry=0.5*(ii+1)+1
!    DO jj=k1,k2
!      !Vin(iEntry)=Vin(iEntry)-LE(jj)*Vin(JL(jj))
!      Vin(iEntry,iElem)=Vin(iEntry,iElem)-IL(iElem)%Entry(jj)*Vin(IL(iELEM)%jEntry(jj),iElem)
!    END DO ! jj
!  END DO ! ii
!  ! backward elimination
!  ! init backward Gauss
!  Vout(nDOFElem,iElem) = Vin(nDOFElem,iElem)/DE(nDOFElem,iElem)
!  DO ii=nMTriangle(iElem)-1,1,-2
!    k1=IU(iElem)%iEntry(ii)
!    k2=IU(iElem)%iEntry(ii+1)
!    iEntry=0.5*(ii+1)
!    DO jj=k1,k2
!      !Vin(iEntry)=Vin(iEntry)-UE(jj)*Vout(JU(jj))
!      Vin(iEntry,iElem)=Vin(iEntry,iElem)-IU(iElem)%Entry(jj)*Vout(IU(iElem)%jEntry(jj),iElem)
!    END DO ! jj
!    Vout(iEntry,iElem)=Vin(iEntry,iElem)/DE(iEntry,iElem)
!  END DO ! ii
!END DO


!test=.FALSE.
!DO iElem=1,PP_nElems
!  DO ii=1,nDOFElem
!    IF(ISNAN(Vout(ii,iElem)))THEN
!      TEST=.TRUE.
!    ENDIF
!  END DO
!END DO
!IF(TEST) THEN
!  print*,'Vout NAN'
!  read*
!END IF
END SUBROUTINE ApplyILU


SUBROUTINE ApplyBILU0(Vin,Vout)
!==================================================================================================================================
! Application of block ILU0 preconditioner
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY:nDOFelem,nDOFLine
USE MOD_Precond_Vars           ,ONLY:nBlockSize
USE MOD_ILU_Vars               ,ONLY:DiagBILU0,XiBILU0,EtaBILU0,ZetaBILU0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: Vcalc(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
INTEGER                      :: i,j,k,l,vni,vnj,vnk,vnl,vn1,vn2,s,r,nBlocks,iElem,iVar,iVar2
LOGICAL                      :: dojump
!==================================================================================================================================

vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)
DO iElem=1,PP_nElems
  ! nullify
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Vout(iVar,i,j,k,iElem)=0.
          Vcalc(iVar,i,j,k)     =Vin(iVar,i,j,k,iElem)
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! L part
  ! xi and eta and zeta direction
  DO k=0,PP_N
    vnk=k*PP_nVar
    DO j=0,PP_N
      vnj=j*PP_nVar
      DO i=0,PP_N
        vni=PP_nVar*i
        r = vn2 * k + vn1 * j + PP_nVar * i
        DO l=0,i-1!PP_N
          vnl=PP_nVar*l
          ! xi
          !s = vn2 * k + vn1 * j + PP_nVar * l
          !IF(s.LT.r)THEN
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                       -XiBILU0  (vni+iVar,vnl+iVar2,j,k,iElem)*Vcalc(iVar2,l,j,k)
                                       !-XiBILU0  (vni+iVar,vnl+iVar2,j,k,iElem)*Vin(iVar2,l,j,k,iElem)

              END DO ! iVar2
            END DO ! iVar
          !END IF
        END DO ! l
          ! eta
        DO l=0,j-1
          vnl=PP_nVar*l
          !s = vn2 * k + vn1 * l + PP_nVar * i
          !IF(s.LT.r)THEN
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                      -EtaBILU0 (vnj+iVar,vnl+iVar2,i,k,iElem)*Vcalc(iVar2,i,l,k)
                                      !-EtaBILU0 (vnj+iVar,vnl+iVar2,i,k,iElem)*Vin(iVar2,i,l,k,iElem)

              END DO ! iVar2
            END DO ! iVar
          !END IF
        END DO ! l
          ! zeta
        DO l=0,k-1
          vnl=PP_nVar*l
         ! s = vn2 * l + vn1 * j + PP_nVar * i
          !IF(s.LT.r)THEN
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                      -ZetaBILU0(vnk+iVar,vnl+iVar2,i,j,iElem)*Vcalc(iVar2,i,j,l)
                                      !-ZetaBILU0(vnk+iVar,vnl+iVar2,i,j,iElem)*Vin(iVar2,i,j,l,iElem)
              END DO ! iVar2
            END DO ! iVar
          !END IF
        END DO ! l
        ! done
      END DO ! i
    END DO ! j
  END DO ! k

!!    dojump=.FALSE.
!!    ! eta direction
!!    DO k=0,PP_N
!!      DO j=0,PP_N
!!        vnj=PP_nVar*j
!!        DO i=0,PP_N
!!          r = vn2 * k + vn1 * j + PP_nVar * i
!!          DO l=0,PP_N
!!            vnl=PP_nVar*l
!!            s = vn2 * k + vn1 * l + PP_nVar * i
!!          END DO ! l
!!          IF(dojump) EXIT
!!        END DO ! i
!!        IF(dojump) EXIT
!!      END DO ! j
!!      IF(dojump) EXIT
!!    END DO ! k
!!  
!!    dojump=.FALSE.
!!    ! zeta direction
!!    DO k=0,PP_N
!!      vnk=PP_nVar*k
!!      DO j=0,PP_N
!!        DO i=0,PP_N
!!          r = vn2 * k + vn1 * j + PP_nVar * i
!!          DO l=0,PP_N
!!            s = vn2 * l + vn1 * j + PP_nVar * i
!!            vnl=PP_nVar*l
!!            END DO ! l
!!          IF(dojump) EXIT
!!        END DO ! i
!!        IF(dojump) EXIT
!!      END DO ! j
!!      IF(dojump) EXIT
!!    END DO ! k

  !Vout(nDOFElem-PP_nVar+1:nDOFElem,iElem) = MATMUL(DiagBILU0(:,:,nBlockSize,iElem),Vcalc(nDOFElem-PP_nVar+1:nDOFElem))

  nBlocks=nBlockSize
  DO k=PP_N,0,-1
    vnk=PP_nVar*k
    DO j=PP_N,0,-1
      vnj=PP_nVar*j
      DO i=PP_N,0,-1
        r = vn2 * k + vn1 * j + PP_nVar * i
        vni=PP_nVar*i
        DO l=PP_N,i+1,-1
          vnl=PP_nVar*l
          ! xi
          !s = vn2 * k + vn1 * j + PP_nVar * l
          !IF(s.GT.r)THEN
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                       -XiBILU0  (vni+iVar,vnl+iVar2,j,k,iElem)*Vout(iVar2,l,j,k,iElem)
              END DO ! iVar2
            END DO ! iVar
          !END IF
        END DO ! l
        DO l=PP_N,j+1,-1
          vnl=PP_nVar*l
          !s = vn2 * k + vn1 * l + PP_nVar * i
          ! stuff for eta
          !IF(s.GT.r)THEN
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                      -EtaBILU0 (vnj+iVar,vnl+iVar2,i,k,iElem)*Vout(iVar2,i,l,k,iElem)

              END DO ! iVar2
            END DO ! iVar
          !END IF
        END DO ! l
        DO l=PP_N,k+1,-1
          vnl=PP_nVar*l
          !s = vn2 * l + vn1 * j + PP_nVar * i
          !IF(s.GT.r)THEN
            ! stuff for zeta
            DO iVar=1,PP_nVar
              DO iVar2=1,PP_nVar
                Vcalc(iVar,i,j,k)=Vcalc(iVar,i,j,k)                             &
                                      -ZetaBILU0(vnk+iVar,vnl+iVar2,i,j,iElem)*Vout(iVar2,i,j,l,iElem)
              END DO ! iVar2
            END DO ! iVar
       !   END IF
        END DO ! l
      !END DO ! l
        ! apply diagonal inverse
        DO iVar=1,PP_nVar
          DO iVar2=1,PP_nVar
            Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)+DiagBILU0(iVar,iVar2,nBlocks,iElem)*Vcalc(iVar2,i,j,k)
          END DO ! iVar2
        END DO ! iVar
        nBlocks=nBlocks-1
      END DO ! i
    END DO ! j
  END DO ! k

END DO ! iElem

END SUBROUTINE ApplyBILU0

SUBROUTINE ApplyBILU0_BCSR(Vin,Vout)
!==================================================================================================================================
! Application of block ILU0 preconditioner
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY:nDOFelem,nDOFLine
USE MOD_Precond_Vars           ,ONLY:nBlockSize
USE MOD_ILU_Vars               ,ONLY:BlockAA,BlockIA,BlockJA,nBDOF,BlockDiag
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,1:nBDOF,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,1:nBDOF,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: Vcalc(1:PP_nVar,1:nBDOF)
INTEGER                      :: ii,jj,k1,k2,iVar,iVar2,jrow,iElem
!==================================================================================================================================


DO iElem=1,PP_nElems
  Vout(:,:,iElem)=0.
  Vcalc(:,:) = Vin(:,:,iElem)
  ! lower sweep
  DO ii=2,nBDOF
    k1=BlockIA(ii)
    k2=BlockDiag(ii)-1
    !print*,'k1,k2',k1,k2
    DO jj=k1,k2
      jrow=BlockJA(jj)
      !print*,'jrwo',jrow
      !read*
      DO iVar=1,PP_nVar
        DO iVar2=1,PP_nVar
          Vcalc(iVar,ii) = Vcalc(iVar,ii)-BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,jrow)
        END DO ! iVar2
      END DO ! ivar 
    END DO ! jj
  END DO ! ii
  ! init
  jj=BlockDiag(nBDOF)
  DO iVar=1,PP_nVar
    DO iVar2=1,PP_nVar
      Vout(iVar,nBDOF,iElem)=Vout(iVar,nBDOF,iElem)+BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,nBDOF)
    END DO ! iVar2
  END DO ! iVar
  ! lower sweep
  DO ii=nBDOF-1,1,-1
    k1=BlockDiag(ii)+1
    k2=BlockIA(ii+1)-1
    !print*,'k1,k2',k1,k2
    DO jj=k1,k2
      jrow=BlockJA(jj)
      !print*,'jrwo',jrow
      DO iVar=1,PP_nVar
        DO iVar2=1,PP_nVar
          Vcalc(iVar,ii)=Vcalc(iVar,ii)-BlockAA(iVar,iVar2,jj,iElem)*Vout(iVar2,jrow,iElem)
        END DO ! iVar2
      END DO ! iVar
    END DO ! jj
    jj=BlockDiag(ii)
    jrow=BlockJA(jj)
    !print*,'jj',jj
    DO iVar=1,PP_nVar
      DO iVar2=1,PP_nVar
        Vout(iVar,ii,iElem)=Vout(iVar,ii,iElem)+BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,jrow)
      END DO ! iVar2
    END DO ! iVar
  END DO ! ii
END DO ! iElem


END SUBROUTINE ApplyBILU0_BCSR


SUBROUTINE ApplyTensorBJ(Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:DENSE_MATMUL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc2(1:PP_nVar,0:PP_N)
INTEGER                      :: iElem,i,j,k,iVar,l,vn1,vn2,vni,vnj,vnk
!===================================================================================================================================

DO iElem=1,PP_nElems
 ! xi direction
  DO j=0,PP_N
    DO k=0,PP_N
      CALL DENSE_MATMUL(nDOFLine,invXI(:,:,j,k,iElem),Vin(:,:,j,k,iElem),Vcalc2)
      Vout(:,:,j,k,iElem)=Vcalc2
    END DO ! j
  END DO ! k
 ! eta direction
  DO k=0,PP_N
    DO i=0,PP_N
      CALL DENSE_MATMUL(nDOFLine,invEta(:,:,i,k,iElem),Vin(:,i,:,k,iElem),Vcalc2)
      Vout(:,i,:,k,iElem)=Vout(:,i,:,k,iElem)+Vcalc2
    END DO ! i
  END DO ! k
  ! zeta direction
  DO j=0,PP_N
    DO i=0,PP_N
      !CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vout(:,i,j,:,iElem),Vout(:,i,j,:,iElem))
      CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vin(:,i,j,:,iElem),Vcalc2)
      Vout(:,i,j,:,iElem)=Vout(:,i,j,:,iElem)+Vcalc2
    END DO ! i
  END DO ! j
  ! missing?
  !Vout=1./3.*Vout -- not working

END DO ! iElem
!Vout=1./3.*Vout

END SUBROUTINE  ApplyTensorBJ

SUBROUTINE ApplyTensorBJOpt(Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:DENSE_MATMUL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc1(1:PP_nVar,0:PP_N)
REAL                         :: Vcalc2(1:PP_nVar,0:PP_N)
REAL                         :: Vcalc3(1:PP_nVar,0:PP_N)
INTEGER                      :: iElem,i,j,k,iVar,l,vnl,iVar2,vni,vnj,vnk
INTEGER                      :: p,q
!===================================================================================================================================



DO iElem=1,PP_nElems
  ! null
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Vout(iVar,i,j,k,iElem)=0.
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! xi & eta & zeta
  DO k=0,PP_N
    vnk=PP_nVar*k
    DO j=0,PP_N
      vnj=PP_nVar*j
      DO i=0,PP_N
        vni=PP_nVar*i
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                                       &
                                    +invXi  (vnl+iVar2,vni+iVar,j,k,iElem)*Vin(iVar2,l,j,k,iElem) &
                                    +invEta (vnl+iVar2,vnj+iVar,i,k,iElem)*Vin(iVar2,i,l,k,iElem) &
                                    +invZeta(vnl+iVar2,vnk+iVar,i,j,iElem)*Vin(iVar2,i,j,l,iElem)
                                ! old
                                    !+invXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vin(iVar2,l,j,k,iElem) &
                                    !+invEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vin(iVar2,i,l,k,iElem) &
                                    !+invZeta(vnk+iVar,vnl+iVar2,i,j,iElem)*Vin(iVar2,i,j,l,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem

!  ! xi direction
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        vni=PP_nVar*i
!        DO l=0,PP_N
!          vnl=PP_nVar*l
!          DO iVar=1,PP_nVar
!            DO iVar2=1,PP_nVar
!              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                             &
!                                    +invXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vin(iVar2,l,j,k,iElem)
!
!              !Vout(iVar2,l,j,k,iElem)=Vout(iVar2,l,j,k,iElem)                             &
!              !                      +invXi  (vnl+iVar2,vni+iVar,j,k,iElem)*Vin(iVar,i,j,k,iElem)
!            END DO ! iVar2
!          END DO ! iVar
!        END DO ! l
!      END DO ! i
!    END DO ! j
!  END DO ! k
!
!  ! eta direction
!  DO k=0,PP_N
!    DO j=0,PP_N
!      vnj=PP_nVar*j
!      DO i=0,PP_N
!        DO l=0,PP_N
!          vnl=PP_nVar*l
!          DO iVar=1,PP_nVar
!            DO iVar2=1,PP_nVar
!              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                             &
!                                    +invEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vin(iVar2,i,l,k,iElem)
!
!              !Vout(iVar2,i,l,k,iElem)=Vout(iVar2,i,l,k,iElem)                             &
!              !                      +invEta (vnl+iVar2,vnj+iVar,i,k,iElem)*Vin(iVar,i,j,k,iElem)
!            END DO ! iVar2
!          END DO ! iVar
!        END DO ! l
!      END DO ! i
!    END DO ! j
!  END DO ! k
!
!  ! zeta direction
!  DO k=0,PP_N
!    vnk=PP_nVar*k
!    DO j=0,PP_N
!      DO i=0,PP_N
!        DO l=0,PP_N
!          vnl=PP_nVar*l
!          DO iVar=1,PP_nVar
!            DO iVar2=1,PP_nVar
!              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                             &
!                                    +invZeta(vnk+iVar,vnl+iVar2,i,j,iElem)*Vin(iVar2,i,j,l,iElem)
!              !Vout(iVar2,i,j,l,iElem)=Vout(iVar2,i,j,l,iElem)                             &
!              !                      +invZeta(vnl+iVar2,vnk+iVar,i,j,iElem)*Vin(iVar,i,j,k,iElem)
!            END DO ! iVar2
!          END DO ! iVar
!        END DO ! l
!      END DO ! i
!    END DO ! j
!  END DO ! k


END SUBROUTINE  ApplyTensorBJOpt

SUBROUTINE ApplyMultiSchwarz(Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_Precond_Vars      ,ONLY:dRdXi,dRdEta,dRdZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:DENSE_MATMUL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!REAL                         :: Vcalc2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vtild1(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vtild2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
INTEGER                      :: iElem,i,j,k,iVar,l,vnl,iVar2,vni,vnj,vnk
INTEGER                      :: p,q
!===================================================================================================================================



DO iElem=1,PP_nElems
  ! null
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Vout(iVar,i,j,k,iElem)=0.
          Vtild1(iVar,i,j,k)     =Vin(iVar,i,j,k,iElem)
          Vtild2(iVar,i,j,k)     =Vin(iVar,i,j,k,iElem)
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! first, all in xi
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        vni=PP_nVar*i
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                                       &
                                    +invXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vin(iVar2,l,j,k,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! apply A
  DO k=0,PP_N
    vnk=PP_nVar*k
    DO j=0,PP_N
      vnj=PP_nVar*j
      DO i=0,PP_N
        vni=PP_nVar*i
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vtild1(iVar,i,j,k)=Vtild1(iVar,i,j,k)                                        &
                                     -dRdXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vout(iVar2,l,j,k,iElem) &
                                     -dRdEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vout(iVar2,i,l,k,iElem) &
                                     -dRdZeta(vnk+iVar,vnl+iVar2,i,j,iElem)*Vout(iVar2,i,j,l,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! eta
  DO k=0,PP_N
    DO j=0,PP_N
      vnj=PP_nVar*j
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                             &
                               +invEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vtild1(iVar2,i,l,k)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! apply A
  DO k=0,PP_N
    vnk=PP_nVar*k
    DO j=0,PP_N
      vnj=PP_nVar*j
      DO i=0,PP_N
        vni=PP_nVar*i
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vtild2(iVar,i,j,k)=Vtild2(iVar,i,j,k)                                     &
                                     -dRdXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vout(iVar2,l,j,k,iElem) &
                                     -dRdEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vout(iVar2,i,l,k,iElem) &
                                     -dRdZeta(vnk+iVar,vnl+iVar2,i,j,iElem)*Vout(iVar2,i,j,l,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! zeta
  DO k=0,PP_N
    vnk=PP_nVar*k
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                                       &
                               +invZeta (vnk+iVar,vnl+iVar2,i,j,iElem)*Vtild2(iVar2,i,j,l)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem

END SUBROUTINE ApplyMultiSchwarz

SUBROUTINE ApplyMultiSchwarzDG(coeff,Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_Precond_Vars      ,ONLY:dRdXi,dRdEta,dRdZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:MatrixVector
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL,INTENT(IN)              :: coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!REAL                         :: Vcalc2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vtild(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
INTEGER                      :: iElem,i,j,k,iVar,l,vnl,iVar2,vni,vnj,vnk
!===================================================================================================================================



DO iElem=1,PP_nElems
  ! null
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Vout(iVar,i,j,k,iElem)=0.
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
  ! first, all in xi
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        vni=PP_nVar*i
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                                       &
                                    +invXi  (vni+iVar,vnl+iVar2,j,k,iElem)*Vin(iVar2,l,j,k,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem

CALL MatrixVector(0.,coeff,Vout,Vtild)
Vtild=Vin-Vtild

DO iElem=1,PP_nElems
  ! eta
  DO k=0,PP_N
    DO j=0,PP_N
      vnj=PP_nVar*j
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                             &
                               +invEta (vnj+iVar,vnl+iVar2,i,k,iElem)*Vtild(iVar2,i,l,k,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem

CALL MatrixVector(0.,coeff,Vout,Vtild)
Vtild=Vin-Vtild

DO iElem=1,PP_nElems
  ! zeta
  DO k=0,PP_N
    vnk=PP_nVar*k
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          DO l=0,PP_N
            vnl=PP_nVar*l
            DO iVar2=1,PP_nVar
              Vout(iVar,i,j,k,iElem)=Vout(iVar,i,j,k,iElem)                                  &
                               +invZeta (vnk+iVar,vnl+iVar2,i,j,iElem)*Vtild(iVar2,i,j,l,iElem)

            END DO ! iVar2
          END DO ! l
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem

END SUBROUTINE ApplyMultiSchwarzDG

SUBROUTINE ApplyTensorProductBJ(Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:DENSE_MATMUL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc(1:PP_nVar)!,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc1(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc3(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
INTEGER                      :: iElem,i,j,k,iVar,l,vn1,vn2,vni,vnj,vnk,iVar2
!===================================================================================================================================


DO iElem=1,PP_nElems
  !! init nullify
  !DO k=0,PP_N
  !  DO j=0,PP_N
  !    DO i=0,PP_N
  !      DO iVar=1,PP_nVar
  !        Vcalc(iVar,i,j,k)=0.
  !      END DO ! iVar
  !    END DO ! i
  !  END DO ! j
  !END DO ! k
  ! xi direction
  DO k=0,PP_N
    DO j=0,PP_N
    !  Vout(:,:,j,k,iElem)=DENSE_MATMUL(nDOFLine,invXI(:,:,j,k,iElem),Vin(:,:,j,k,iElem))
      CALL DENSE_MATMUL(nDOFLine,invXI  (:,:,j,k,iElem),Vin(:,:,j,k,iElem),Vout(:,:,j,k,iElem))
      !CALL DENSE_MATMUL(nDOFLine,invXI  (:,:,j,k,iElem),Vin(:,:,j,k,iElem),Vcalc1)
    END DO ! j
  END DO ! k
 ! eta direction
  DO k=0,PP_N
    DO i=0,PP_N
      !CALL DENSE_MATMUL(nDOFLine,invEta (:,:,i,k,iElem),Vin(:,i,:,k,iElem),Vcalc2)
      CALL DENSE_MATMUL(nDOFLine,invEta(:,:,i,k,iElem),Vout(:,i,:,k,iElem),Vout(:,i,:,k,iElem))
    END DO ! i
  END DO ! k
  ! zeta direction
  DO j=0,PP_N
    DO i=0,PP_N
      !CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vin(:,i,j,:,iElem),Vcalc3)
      CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vout(:,i,j,:,iElem),Vout(:,i,j,:,iElem))
    END DO ! i
  END DO ! j
  

!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        DO iVar=1,PP_nVar
!          Vout(iVar,i,j,k,iElem)=Vcalc1(iVar,i,j,k)*Vcalc2(iVar,i,j,k)*Vcalc3(iVar,i,j,k)
!        END DO ! iVar
!      END DO ! i
!    END DO ! j
!  END DO ! k
END DO ! iElem

END SUBROUTINE  ApplyTensorProductBJ

END MODULE MOD_ApplyPreconditioner
