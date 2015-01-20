MODULE MOD_LinearSolver
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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

PUBLIC:: LinearSolver_CGS, LinearSolver, LinearSolver_StabBiCGSTAB
!===================================================================================================================================
CONTAINS

SUBROUTINE LinearSolver_CGS(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using CGS
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars,ONLY:eps_LinearSolver,maxIter_LinearSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iter
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Q(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Tvec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,sigma,alpha,beta,Norm_R
INTEGER          :: chance
!===================================================================================================================================
chance=0
! Save old time step
1001 CONTINUE
Un=U
Uold = U
! Compute A*U^n
CALL MatrixVector(t,Coeff)
! Compute first residual b-A*U^n 
! Equation System:
! U_n+1 = U_n + dt*(DG_Operator(U_n+1)+Sources_n+1)
! (I-dt*DG_Operator)U_n+1 = U_n+dt(Sources_n+1)
!          A          x   =        b
! ==>
! R0 = b - Ax_0
!    = U_n + dt*Sources_n+1 - (I-dt*DG_Operator)*x_0
!    = U_n + dt*Sources_n+1 - x_0 + dt*DG_Operator*x_0
! with x_0 = U_n and Ut from MatrixVectorSource (see above)
!    = U_n - Ut
R0 = Uold - Ut
! Init
P=R0
R=R0
Tvec=R0

CALL VectorDotProduct(R,R,Norm_R0)
Norm_R0=SQRT(Norm_R0)

DO iter=1,maxIter_LinearSolver
  U=P
  CALL MatrixVector(t,Coeff)
  V=Ut
  CALL VectorDotProduct(V,R0,alpha)
  CALL VectorDotProduct(R,R0,sigma)
  alpha=sigma/alpha
  Q=Tvec-alpha*V
  Un=Un+alpha*(Tvec+Q)
  U=Tvec+Q
  CALL MatrixVector(t,Coeff)
  R=R-alpha*Ut
  CALL VectorDotProduct(R,R0,beta) 
  beta=beta/sigma
  Tvec=R+beta*Q
  P=Tvec+beta*(Q+beta*P)
  CALL VectorDotProduct(R,R,Norm_R)
  Norm_R=SQRT(Norm_R)
! Debug Ausgabe, Anzahl der Iterationen...
  !SWRITE(*,*)'iter,Res: ',iter,Norm_R/Norm_R0
  IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
  !IF((Norm_R.LE.eps_LinearSolver*Norm_R0)) THEN
    U=Un
    ! Debug Ausgabe, Anzahl der Iterationen...
    !SWRITE(*,*)'iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    RETURN
  ENDIF
  IF((iter.GT.maxIter_LinearSolver/2).AND.(chance.LT.2)) THEN
    U=0.5*(Uold+Un)
    chance=chance+1
    WRITE(*,*)'SAUARSCH,iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    GOTO 1001
  END IF 
END DO !iter
WRITE(*,*)'No Convergence: maxIter,last,initial residual',iter,Norm_R,Norm_R0
WRITE(*,*)'flummi',Norm_R,eps_LinearSolver,Norm_R0
STOP
END SUBROUTINE LinearSolver_CGS



SUBROUTINE LinearSolver(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab 
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars,ONLY:eps_LinearSolver,maxIter_LinearSolver
USE MOD_Equation_Vars, ONLY : eps0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iter
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,sigma,alpha,Norm_T2,omega,beta,Norm_R,eps0inv
INTEGER          :: chance
!===================================================================================================================================
chance=0
eps0inv = 1./eps0
1002 CONTINUE 
! Save old time step
Un=U
UOld=U
! Compute A*U^n
!CALL MatrixVector(t,Coeff)
CALL MatrixVectorSource(t,Coeff)
! Compute first residual b-A*U^n 
! Equation System:
! U_n+1 = U_n + dt*(DG_Operator(U_n+1)+Sources_n+1)
! (I-dt*DG_Operator)U_n+1 = U_n+dt(Sources_n+1)
!          A          x   =        b
! ==>
! R0 = b - Ax_0
!    = U_n + dt*Sources_n+1 - (I-dt*DG_Operator)*x_0
!    = U_n + dt*Sources_n+1 - x_0 + dt*DG_Operator*x_0
! with x_0 = U_n and Ut from MatrixVectorSource (see above)
!    = U_n - Ut
R0 = UOld - Ut

! Init
P=R0
R=R0

CALL VectorDotProduct(R,R,Norm_R0)
Norm_R0=SQRT(Norm_R0)

DO iter=1,maxIter_LinearSolver
  U=P
  CALL MatrixVector(t,Coeff)
  V=Ut
  CALL VectorDotProduct(V,R0,sigma)
  CALL VectorDotProduct(R,R0,alpha)
  alpha=alpha/sigma
  S = R - alpha*V
  U=S
  CALL MatrixVector(t,Coeff)
  TVec=Ut
  CALL VectorDotProduct(TVec,TVec,Norm_T2)
  CALL VectorDotProduct(TVec,S,omega)
  omega=omega/Norm_T2
  Un=Un+alpha*P+omega*S
  R=S-omega*TVec
  CALL VectorDotProduct(R,R0,beta)
  beta=beta/(omega*sigma)
  P=R+beta*(P-omega*V)
  CALL VectorDotProduct(R,R,Norm_R)
  Norm_R=SQRT(Norm_R)
  !IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
  IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
    U=Un
    ! Debug Ausgabe, Anzahl der Iterationen...
    WRITE(*,*)'iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    RETURN
  ENDIF
  IF((iter.GT.maxIter_LinearSolver/2).AND.(chance.LT.2)) THEN
    U=0.5*(UOld+Un)
    chance=chance+1
    WRITE(*,*)'SAUARSCH,iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    GOTO 1002
  END IF 
END DO !iter
WRITE(*,*)'No Convergence: maxIter,last,initial residual',iter,Norm_R,Norm_R0
WRITE(*,*)'flummi',Norm_R,eps_LinearSolver,Norm_R0
STOP
END SUBROUTINE LinearSolver

SUBROUTINE LinearSolver_StabBiCGSTAB(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using stabilized BiCGStab 
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars,ONLY:eps_LinearSolver,maxIter_LinearSolver
USE MOD_Equation_Vars, ONLY : eps0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iter
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,sigma,alpha,Norm_T2,omega,beta,Norm_R,eps0inv, Norm_R02, Norm_V, Norm_S
INTEGER          :: chance
!===================================================================================================================================
chance=0
eps0inv = 1./eps0
1002 CONTINUE 
! Save old time step
Un=U
UOld=U
! Compute A*U^n
!CALL MatrixVector(t,Coeff)
CALL MatrixVectorSource(t,Coeff)
! Compute first residual b-A*U^n 
! Equation System:
! U_n+1 = U_n + dt*(DG_Operator(U_n+1)+Sources_n+1)
! (I-dt*DG_Operator)U_n+1 = U_n+dt(Sources_n+1)
!          A          x   =        b
! ==>
! R0 = b - Ax_0
!    = U_n + dt*Sources_n+1 - (I-dt*DG_Operator)*x_0
!    = U_n + dt*Sources_n+1 - x_0 + dt*DG_Operator*x_0
! with x_0 = U_n and Ut from MatrixVectorSource (see above)
!    = U_n - Ut
R0 = UOld - Ut

! Init
P=R0
R=R0

CALL VectorDotProduct(R,R,Norm_R0)
Norm_R0=SQRT(Norm_R0)

DO iter=1,maxIter_LinearSolver
  U=P
  CALL MatrixVector(t,Coeff)
  V=Ut
  CALL VectorDotProduct(V,R0,sigma)
  CALL VectorDotProduct(R0,R0,Norm_R02)
  CALL VectorDotProduct(V,V,Norm_V) 
 ! IF((sigma.GT.1E-4*Norm_R02*Norm_V).AND.(iter.GT.10)) THEN
    CALL VectorDotProduct(R,R0,alpha)
    alpha=alpha/sigma
    S = R - alpha*V
    CALL VectorDotProduct(S,S,Norm_S)
    IF(Norm_S.GT.1.E-12) THEN
      U=S
      CALL MatrixVector(t,Coeff)
      TVec=Ut
      CALL VectorDotProduct(TVec,TVec,Norm_T2)
      CALL VectorDotProduct(TVec,S,omega)
      omega=omega/Norm_T2
      Un=Un+alpha*P+omega*S
      R=S-omega*TVec
      CALL VectorDotProduct(R,R0,beta)
      beta=beta/(omega*sigma)
      P=R+beta*(P-omega*V)
      CALL VectorDotProduct(R,R,Norm_R)
      Norm_R=SQRT(Norm_R)
      !IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
      IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
        U=Un
        ! Debug Ausgabe, Anzahl der Iterationen...
        WRITE(*,*)'iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
        RETURN
      ENDIF
      IF((iter.GT.maxIter_LinearSolver/2).AND.(chance.LT.2)) THEN
        U=0.5*(UOld+Un)
        chance=chance+1
        WRITE(*,*)'SAUARSCH,iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
        GOTO 1002
      END IF 
    ELSE
      Un=Un+alpha*P
      R=S
      CALL VectorDotProduct(R,R,Norm_R)
      Norm_R=SQRT(Norm_R)
      !IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
      IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
        U=Un
        ! Debug Ausgabe, Anzahl der Iterationen...
        WRITE(*,*)'iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
        RETURN
      ENDIF
    END IF
 ! ELSE
 !   U=Un
 !   SWRITE(*,*)'Restart Stab,iter,t,sigma,Norm_R02*Norm_V',iter,t,sigma,Norm_R02*Norm_V
 !   GOTO 1002
 ! END IF
END DO !iter
WRITE(*,*)'No Convergence: maxIter,last,initial residual',iter,Norm_R,Norm_R0
WRITE(*,*)'flummi',Norm_R,eps_LinearSolver,Norm_R0
STOP
END SUBROUTINE LinearSolver_StabBiCGSTAB

SUBROUTINE MatrixVector(t,Coeff)
!===================================================================================================================================
! Computes Matrix Vector Product y=A*x for linear Equations only
! Matrix A = I - Coeff*R
! Attention: Vector x is always U 
! Attention: Result y is always stored in Ut
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut
USE MOD_PreProc
USE MOD_DG,ONLY:DGTimeDerivative_WoSource_weakForm
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL DGTimeDerivative_WoSource_weakForm(t,t,0)
! y = (I-Coeff*R)*x = x - Coeff*R*x 
Ut = U - Coeff*Ut
END SUBROUTINE MatrixVector

SUBROUTINE MatrixVectorSource(t,Coeff)
!===================================================================================================================================
! Computes Matrix Vector Product y=A*x for linear Equations only
! Matrix A = I - Coeff*R
! Attention: Vector x is always U 
! Attention: Result y is always stored in Ut
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:U,Ut
USE MOD_PreProc
USE MOD_DG,ONLY:DGTimeDerivative_weakForm
USE MOD_Equation,      ONLY: CalcSource
USE MOD_Equation,ONLY:DivCleaningDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL DGTimeDerivative_weakForm(t,t,0)
!CALL DivCleaningDamping()
! y = (I-Coeff*R)*x = x - Coeff*R*x 
Ut = U - Coeff*Ut

END SUBROUTINE MatrixVectorSource


SUBROUTINE VectorDotProduct(a,b,resu)
!===================================================================================================================================
! Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: a(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL,INTENT(IN)   :: b(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar,i,j,k,iElem
!===================================================================================================================================

resu=0.
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          resu=resu + a(iVar,i,j,k,iElem)*b(iVar,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO
END DO

#ifdef MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

END SUBROUTINE VectorDotProduct

END MODULE MOD_LinearSolver
