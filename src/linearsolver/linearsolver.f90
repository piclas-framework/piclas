#include "boltzplatz.h"

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


INTERFACE LinearSolver
  MODULE PROCEDURE LinearSolver
!  MODULE PROCEDURE LinearSolver_GMRES_P
!  MODULE PROCEDURE LinearSolver_StabBiCGSTAB  
!  MODULE PROCEDURE LinearSolver_StabBiCGSTAB_P
!  MODULE PROCEDURE LinearSolver_BiCGSTAB_P
!  MODULE PROCEDURE LinearSolver_BiCGSTAB_PM
!  MODULE PROCEDURE LinearSolver_BiCGSTAB
END INTERFACE

PUBLIC:: InitLinearSolver, FinalizeLinearSolver
PUBLIC:: LinearSolver!,LinearSolver_BiCGSTAB_PM,LinearSolver_GMRES_P
!===================================================================================================================================

CONTAINS

SUBROUTINE InitLinearSolver()
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LinearSolver_Vars
USE MOD_Interpolation_Vars,   ONLY:wGP
USE MOD_Mesh_Vars,            ONLY:MeshInitIsDone,sJ
USE MOD_Interpolation_Vars,   ONLY:InterpolationInitIsDone
USE MOD_ReadInTools,          ONLY:GETINT,GETREAL,GETLOGICAL
USE MOD_Precond,              ONLY:InitPrecond
USE MOD_Predictor,            ONLY:InitPredictor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER    :: i,j,k,iElem
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.LinearSolverInitIsDone)THEN
   CALL abort(&
__STAMP__&
,'InitImplicit not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LINEAR SOLVER...'

nGP2D=(PP_N+1)**2
nGP3D=nGP2D*(PP_N+1)
nDOFLine=PP_nVar*(PP_N+1)
nDOFside=PP_nVar*nGP2D
nDOFelem=PP_nVar*nGP3D
nDOFGlobal=nDOFelem*PP_nElems

ALLOCATE(ImplicitSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(LinSolverRHS  (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
!#if (PP_TimeDiscMethod==100)
!  ALLOCATE(FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1))
!#endif
!#if (PP_TimeDiscMethod==102)
!  ALLOCATE(FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:6))
!#endif 

#if (PP_TimeDiscMethod==104)
!ALLOCATE(Q(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)) LinSolverRHS
! the time derivative computed at the actual Newton iteration value "xk"
ALLOCATE(R_Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Eps2Newton =GETREAL('EpsNewton','0.001')
Eps2Newton =Eps2Newton**2
nNewtonIter=GETINT('nNewtonIter','20')
EisenstatWalker=GETLOGICAL('EisenstatWalker','.FALSE.')
gammaEW=GETREAL('gammaEW','0.9')
#endif

nDofGlobalMPI=nDofGlobal
#ifdef MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,nDofGlobalMPI,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

nRestarts             = GETINT('nRestarts','1')
eps_LinearSolver      = GETREAL('eps_LinearSolver','1e-3')
epsTilde_LinearSolver = eps_LinearSolver!GETREAL('epsTilde_LinearSolver')
eps2_LinearSolver     = eps_LinearSolver *eps_LinearSolver 
maxIter_LinearSolver  = GETINT('maxIter_LinearSolver','60')

nKDim=GETINT('nKDim','25')
nInnerIter=0
totalIterLinearSolver = 0

ALLOCATE(Mass(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        Mass(1:PP_nVar,i,j,k,iElem)=wGP(i)*wGP(j)*wGP(k) /sJ(i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iElem=1,PP_nElems
IF(.NOT.GETLOGICAL('withmass','F')) mass=1.




LinSolver= GETINT('LinSolver','2')
SELECT CASE(LinSolver)
CASE(1)
  SWRITE(*,*) ' Linear solver: CGS with right preconditioner'
CASE(2)
  SWRITE(*,*) ' Linear solver: BiCGSTAB with right preconditioner'
CASE(3)
  SWRITE(*,*) ' Linear solver: stabilized BiCGSTAB with right preconditioner'
CASE(4)
  SWRITE(*,*) ' Linear solver: GMRES with right preconditioner'
CASE(5)
  SWRITE(*,*) ' Linear solver: BiCGSTAB with left and right preconditioner'
CASE(6)
  SWRITE(*,*) ' Linear solver: BiCGSTAB with left preconditioner'
CASE(7)
  SWRITE(*,*) ' Linear solver: BiCGSTAB(l)'
  ldim  = GETINT('ldim','1')
  maxIter_LinearSolver=maxIter_LinearSolver/ldim+1
  SWRITE(*,'(A,I4)') ' New number of max. Iterations: ', maxIter_LinearSolver
CASE DEFAULT
  CALL abort(&
__STAMP__ &
,'WRONG TYPE OF LINEAR SOLVER:',LinSolver,999.)
END SELECT

LinearSolverInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LINEAR SOLVER DONE!'
! init predictor
CALL InitPredictor()
! init preconditoner
CALL InitPrecond()


END SUBROUTINE InitLinearSolver

SUBROUTINE LinearSolver(t,Coeff)
!==================================================================================================================================
! Selection between different linear solvers
!==================================================================================================================================
! MODULES
USE MOD_LinearSolver_Vars              ,ONLY: LinSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)  :: t,Coeff
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

SELECT CASE(LinSolver)
CASE(1)
 CALL LinearSolver_CGS(t,Coeff)
CASE(2)
 CALL LinearSolver_BiCGSTAB_PM(t,Coeff)
CASE(3)
 CALL LinearSolver_StabBiCGSTAB_P(t,Coeff)
CASE(4)
 CALL LinearSolver_GMRES_P(t,Coeff)
CASE(5)
 CALL LinearSolver_BiCGSTAB_LRP(t,Coeff)
CASE(6)
 CALL LinearSolver_BiCGSTAB_LP(t,Coeff)
CASE(7)
 CALL LinearSolver_BiCGSTABl(t,Coeff)
END SELECT

END SUBROUTINE LinearSolver

SUBROUTINE LinearSolver_CGS(t,Coeff)
!==================================================================================================================================
! Solves Linear system Ax=b using CGS
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U,Ut
USE MOD_LinearSolver_Vars,    ONLY:eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver,nInnerIter
USE MOD_LinearSolver_Vars,    ONLY:LinSolverRHS,ImplicitSource,nRestarts
USE MOD_LinearOperator,       ONLY:MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iter
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAl             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Q(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Tvec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,sigma,alpha,beta,Norm_R
REAL             :: AbortCrit
INTEGER          :: iterLinSolver,Restart
! preconditioner
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TvecQt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!==================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n+1 -( I - dt*DG_Operator ) U^n
!    = dt*Source^n+1  + dt*DG_Operator U^n 
!    = dt*Ut + dt*Source^n+1

! store here for later use
Un   = U ! here, n stands for U^n
Uold = U
Restart=0
nInnerIter = 0
! LinSolverRHS and X0 = U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
! compute  A*U^n
CALL VectorDotProduct(R0,R0,Norm_R0)
Norm_R0=SQRT(Norm_R0)

! Init
P=R0
R=R0
Tvec=R0
AbortCrit = Norm_R0*eps_LinearSolver

DO WHILE (Restart.LT.nRestarts) ! maximum number of trials with CGS
  DO iterLinSolver=1, maxIter_LinearSolver
    ! Preconditioner
    CALL Preconditioner(coeff,P,Pt)
    CALL MatrixVector(t,coeff,Pt,V)
    CALL VectorDotProduct(V,R0,alpha) ! sig,alpha turned compared to BiCGSTAB ! caution
    CALL VectorDotProduct(R,R0,sigma)
    alpha = sigma / alpha
    Q=Tvec - alpha*V
    ! Preconditioner
    TvecQt=Tvec+Q
    CALL PRECONDITIONER(coeff,TvecQt,TvecQt)
    CALL MatrixVector(t,coeff,TvecQt,V) ! we are using V because it is not needed again
    Un=Un+alpha*TvecQt
    ! R_j+1=R_j + alpha A(uj+qj)
    R = R - alpha*V
    CALL VectorDotProduct(R,R0,beta)
    beta = beta/sigma
    Tvec = R+beta*Q
    P    = Tvec+beta*(Q+beta*P)
    CALL VectorDotProduct(R,R,Norm_R)
    Norm_R=SQRT(Norm_R)
    ! test if success
    IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12))THEN
      U=Un
      nInnerIter=nInnerIter+iterLinSolver
      totalIterLinearSolver=totalIterLinearSolver+nInnerIter
#ifdef DLINANALYZE
      CALL CPU_TIME(tE)
      ! Debug Ausgabe, Anzahl der Iterationen...
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in CGS        : ',tE-tS
      SWRITE(UNIT_stdOut,'(A23,E16.8)')   ' Norm_R0            : ',Norm_R0
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Ratio Precond/DG   : ',tPrecond/tDG
#endif /* DLINANALYZE */
      RETURN
    ENDIF
  END DO ! iterLinSolver
  ! restart with new U
  ! LinSolverRHS and X0 = U
  !  U              = 0.5*(Uold+Un)
  ImplicitSource = 0.
  U             = Un
  ! LinSolverRHS and X0 = U
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
  ! compute  A*U^n
  CALL VectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
  P   = R0
  R   = R0
  Tvec= R0
  nInnerIter=nInnerIter+iterLinSolver
  Restart = Restart + 1
END DO ! while chance < 2 

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIROOT) CALL abort(&
__STAMP__ &
,'CGS NOT CONVERGED WITH RESTARTS AND CGS ITERATIONS:',Restart,REAL(nInnerIter+iterLinSolver))
END SUBROUTINE LinearSolver_CGS


SUBROUTINE LinearSolver_BiCGStab(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab 
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,        ONLY: U
USE MOD_LinearSolver_Vars,  ONLY: eps_LinearSolver,maxIter_LinearSolver!,epsTilde_LinearSolver
USE MOD_LinearSolver_Vars,  ONLY: LinSolverRHS,ImplicitSource
USE MOD_LinearOperator, ONLY: MatrixVector, MatrixVectorSource, VectorDotProduct
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
1002 CONTINUE 
! Save old time step
! Compute A*U^n
!CALL MatrixVector(t,Coeff)
! Compute first residual b-A*U^n 
! Equation System:
! U_n+1 = U_n + dt*(DG_Operator(U_n+1)+Sources_n+1)
! (I-dt*DG_Operator)U_n+1 = U_n+dt(Sources_n+1)
!          A          x   =        b
! ==>
! R0 = b - Ax_0
!    = U_n + dt*Sources_n - (I-dt*DG_Operator)*x_0
!    = U_n + dt*Sources_n - x_0 + dt*DG_Operator*x_0
!    = U_n + dt*Sourvce_n - U_n + dt*DG_Operator*U_n
! with x_0 = U_n and Ut from MatrixVectorSource (see above)
!    = U_n - Ut


Un=U
UOld=U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*DG_Operator(Un)

! Init
P=R0
R=R0

CALL VectorDotProduct(R,R,Norm_R0)
Norm_R0=SQRT(Norm_R0)

DO iter=1,maxIter_LinearSolver
  CALL MatrixVector(t,Coeff,P,V)
  CALL VectorDotProduct(V,R0,sigma)
  CALL VectorDotProduct(R,R0,alpha)
  alpha=alpha/sigma
  S = R - alpha*V
  CALL MatrixVector(t,Coeff,S,TVec)
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
    SWRITE(*,*)'iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    RETURN
  ENDIF
  IF((iter.GT.maxIter_LinearSolver).AND.(chance.LT.2)) THEN
    U=0.5*(UOld+Un)
    chance=chance+1
    ! restart
    ImplicitSource=0.
    LinSolverRHS  =U
    SWRITE(*,*)'SAUARSCH,iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0
    GOTO 1002
  END IF 
END DO !iter
SWRITE(*,*)'No Convergence: maxIter,last,initial residual',iter,Norm_R,Norm_R0
SWRITE(*,*)'flummi',Norm_R,eps_LinearSolver,Norm_R0
CALL abort(&
__STAMP__&
,' No convergence!')
END SUBROUTINE LinearSolver_BiCGSTAB


SUBROUTINE LinearSolver_BiCGStab_P(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab with right preconditioner P_r
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,             ONLY: U
USE MOD_TimeDisc_Vars,       ONLY: dt
USE MOD_LinearSolver_Vars,   ONLY: eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver
USE MOD_LinearSolver_Vars,   ONLY: LinSolverRHS,ImplicitSource
USE MOD_ApplyPreconditioner, ONLY:Preconditioner
USE MOD_LinearOperator,      ONLY: MatrixVector, MatrixVectorSource, VectorDotProduct
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,Norm_R,Norm_T2
INTEGER          :: iterLinSolver,chance
REAL             :: alpha,sigma,omega,beta
! preconditioner
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: St(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n -( I - dt*DG_Operator ) U^n
!    = dt*DG_Operator U^n 
!    = dt * Ut


! store here for later use
Un   = U ! here, n stands for U^n
Uold = U

chance=0
DO WHILE (chance.LT.2)  ! maximum of two trials with BiCGStab inner interation

  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut
  ! compute  A*U^n
  CALL VectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
  P  = R0
  R  = R0
  
  DO iterLinSolver = 1, maxIter_LinearSolver  ! two trials with half of iterations
    ! Preconditioner
    CALL Preconditioner(coeff,P,Pt)
    CALL MatrixVector(t,Coeff,Pt,V)

    CALL VectorDotProduct(V,R0,sigma)
    CALL VectorDotProduct(R,R0,alpha)

    alpha=alpha/sigma
    S = R - alpha*V

    ! Preconditioner
    CALL Preconditioner(coeff,S,St)
    CALL MatrixVector(t,Coeff,St,TVec)

    CALL VectorDotProduct(TVec,TVec,Norm_T2)
    CALL VectorDotProduct(TVec,S,omega)
    omega=omega/Norm_T2

    Un=Un+alpha*Pt+omega*St
    R=S-omega*TVec
    CALL VectorDotProduct(R,R0,beta)
    beta=beta/(omega*sigma)
    P=R+beta*(P-omega*V)
    CALL VectorDotProduct(R,R,Norm_R)
    Norm_R=SQRT(Norm_R)
    ! test if success
    IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
      U=Un
      totalIterLinearSolver=totalIterLinearSolver+iterLinSolver
      ! Debug Ausgabe, Anzahl der Iterationen...
      SWRITE(*,*)'Iter LinSolver: ',iterLinSolver
      SWRITE(*,*)'t             : ',t
      SWRITE(*,*)'Norm_R        : ',Norm_R
      SWRITE(*,*)'Norm_R0       : ',Norm_R0
      RETURN
    ENDIF
  END DO ! iterLinSolver
  ! restart with new U
  U    = 0.5*(Uold+Un)
  Un   = U
  Uold = U
  ! restart
  ImplicitSource=0.
  LinSolverRHS  =U
  chance = chance+1
  totalIterLinearSolver=totalIterLinearSolver+iterLinSolver
  SWRITE(*,*) 'No convergence during first half of iterations.'
  SWRITE(*,*) 'Iter,t,Norm_R,Norm_R0',iterLinSolver,t,Norm_R,Norm_R0

END DO ! while chance < 2 

SWRITE(*,*)'Norm_R        : ',Norm_R
SWRITE(*,*)'Norm_R0       : ',Norm_R0
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'BiCGSTAB NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',chance,REAL(iterLinSolver))

END SUBROUTINE LinearSolver_BiCGSTAB_P

SUBROUTINE LinearSolver_BiCGStab_PM(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab with right preconditioner P_r
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U
USE MOD_TimeDisc_Vars,        ONLY:dt
USE MOD_LinearSolver_Vars,    ONLY:eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver,nInnerIter
USE MOD_LinearSolver_Vars,    ONLY:LinSolverRHS,ImplicitSource,nRestarts
USE MOD_LinearOperator,       ONLY:MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,Norm_R,Norm_T2
INTEGER          :: iterLinSolver,Restart
REAL             :: alpha,sigma,omega,beta
REAL             :: AbortCrit
! preconditioner
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: St(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#ifdef DLINANALYZE
REAL             :: tS,tE,tStart,tend
REAL             :: tDG, tPrecond
#endif /* DLINANALYZE */
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n+1 -( I - dt*DG_Operator ) U^n
!    = dt*Source^n+1  + dt*DG_Operator U^n 
!    = dt*Ut + dt*Source^n+1

#ifdef DLINANALYZE
tPrecond=0.
tDG=0.
CALL CPU_TIME(tS)
#endif /* DLINANALYZE */

! store here for later use
Un   = U ! here, n stands for U^n
Uold = U
Restart=0
nInnerIter = 0
! LinSolverRHS and X0 = U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
! compute  A*U^n
CALL VectorDotProduct(R0,R0,Norm_R0)

! check this out
alpha=Norm_R0
Norm_R0=SQRT(Norm_R0)

P  = R0
R  = R0
AbortCrit = Norm_R0*eps_LinearSolver

DO WHILE (Restart.LT.nRestarts)  ! maximum of two trials with BiCGStab inner interation
  DO iterLinSolver = 1, maxIter_LinearSolver  ! two trials with half of iterations
    ! Preconditioner
#ifdef DLINANALYZE
        CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    CALL Preconditioner(coeff,P,Pt)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tPrecond=tPrecond+tend-tStart
    ! matrix vector
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    CALL MatrixVector(t,coeff,Pt,V)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tDG=tDG+tend-tStart
#endif /* DLINANALYZE */
    CALL VectorDotProduct(V,R0,sigma)
    !CALL VectorDotProduct(R,R0,alpha)

    alpha=alpha/sigma
    S = R - alpha*V

#ifdef DLINANALYZE
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    ! Preconditioner
    CALL Preconditioner(coeff,S,St)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tPrecond=tPrecond+tend-tStart
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    ! matrix vector
    CALL MatrixVector(t,coeff,St,TVec)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tDG=tDG+tend-tStart
#endif /* DLINANALYZE */

    CALL VectorDotProduct(TVec,TVec,Norm_T2)
    CALL VectorDotProduct(TVec,S,omega)
    omega=omega/Norm_T2

    Un=Un+alpha*Pt+omega*St
    R=S-omega*TVec
    CALL VectorDotProduct(R,R0,alpha)
    beta=alpha/(omega*sigma)
    P=R+beta*(P-omega*V)
    CALL VectorDotProduct(R,R,Norm_R)
    Norm_R=SQRT(Norm_R)
    ! test if success
    IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12)) THEN
      U=Un
      nInnerIter=nInnerIter+iterLinSolver
      totalIterLinearSolver=totalIterLinearSolver+nInnerIter
#ifdef DLINANALYZE
      CALL CPU_TIME(tE)
      ! Debug Ausgabe, Anzahl der Iterationen...
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in BiCGSTAB   : ',tE-tS
      SWRITE(UNIT_stdOut,'(A23,E16.8)')   ' Norm_R0            : ',Norm_R0
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Ratio Precond/DG   : ',tPrecond/tDG
#endif /* DLINANALYZE */
      RETURN
    ENDIF
  END DO ! iterLinSolver
  ! restart with new U
  ! LinSolverRHS and X0 = U
!  U              = 0.5*(Uold+Un)
  ImplicitSource = 0.
  U             = Un
  ! LinSolverRHS and X0 = U
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
  ! compute  A*U^n
  CALL VectorDotProduct(R0,R0,Norm_R0)
  alpha=Norm_R0
  Norm_R0=SQRT(Norm_R0)
  P  = R0
  R  = R0
  nInnerIter=nInnerIter+iterLinSolver
  Restart = Restart + 1
END DO ! while chance < 2 

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'BiCGSTAB NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',Restart,REAL(nInnerIter+iterLinSolver))

END SUBROUTINE LinearSolver_BiCGSTAB_PM

SUBROUTINE LinearSolver_StabBiCGSTAB(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using stabilized BiCGStab 
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,       ONLY:U
USE MOD_LinearSolver_Vars, ONLY:eps_LinearSolver,maxIter_LinearSolver!,epsTilde_LinearSolver
USE MOD_Equation_Vars, ONLY:eps0,c_corr
USE MOD_LinearSolver_Vars, ONLY:LinSolverRHS,ImplicitSource
USE MOD_LinearOperator, ONLY: MatrixVector, MatrixVectorSource, VectorDotProduct
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

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n -( I - dt*DG_Operator ) U^n
!    = dt*DG_Operator U^n 
!    = dt * Ut


! store here for later use
Un   = U ! here, n stands for U^n
Uold = U

chance=0
DO WHILE (chance.LT.2)  ! maximum of two trials with BiCGStab inner interation

  ! Compute residual vector
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*DG_Operator(Un)
  
  P=R0
  R=R0
  
  CALL VectorDotProduct(R,R,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
  
  DO iter=1,maxIter_LinearSolver
    CALL MatrixVector(t,Coeff,P,V)
    CALL VectorDotProduct(V,R0,sigma)
    CALL VectorDotProduct(R0,R0,Norm_R02)
    CALL VectorDotProduct(V,V,Norm_V) 
   ! IF((sigma.GT.1E-4*Norm_R02*Norm_V).AND.(iter.GT.10)) THEN
      CALL VectorDotProduct(R,R0,alpha)
      alpha=alpha/sigma
      S = R - alpha*V
      CALL VectorDotProduct(S,S,Norm_S)
      IF(Norm_S.GT.1.E-12) THEN
        CALL MatrixVector(t,Coeff,S,TVec)
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
        IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
          U=Un
          ! Debug Ausgabe, Anzahl der Iterationen...
          SWRITE(*,*)'Iter LinSolver: ',iter
          SWRITE(*,*)'t             : ',t
          SWRITE(*,*)'Norm_R        : ',Norm_R
          SWRITE(*,*)'Norm_R0       : ',Norm_R0
          RETURN
        ENDIF
      ELSE
        Un=Un+alpha*P
        R=S
        CALL VectorDotProduct(R,R,Norm_R)
        Norm_R=SQRT(Norm_R)
        IF((Norm_R.LE.eps_LinearSolver*Norm_R0).OR.(Norm_R.LT.1.E-12)) THEN
          U=Un
          ! Debug Ausgabe, Anzahl der Iterationen...
          SWRITE(*,*)'Iter LinSolver: ',iter
          SWRITE(*,*)'t             : ',t
          SWRITE(*,*)'Norm_R        : ',Norm_R
          SWRITE(*,*)'Norm_R0       : ',Norm_R0
          RETURN
        ENDIF
      END IF
  END DO !iter
    ! restart with new U
  U    = 0.5*(Uold+Un)
  Un   = U
  Uold = U
  ! restart
  ImplicitSource=0.
  LinSolverRHS  =U
  chance = chance+1
  SWRITE(*,*) 'No convergence during first half of iterations.'
  SWRITE(*,*) 'Iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0

END DO ! while

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'StabBiCGSTAB NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',chance)

END SUBROUTINE LinearSolver_StabBiCGSTAB


SUBROUTINE LinearSolver_StabBiCGSTAB_P(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using stabilized BiCGStab 
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Precond_Vars,            ONLY: PrecondType
USE MOD_LinearSolver_Vars,       ONLY: eps_LinearSolver,maxIter_LinearSolver,epsTilde_LinearSolver,totalIterLinearSolver
USE MOD_LinearSolver_Vars,       ONLY: LinSolverRHS,ImplicitSource
USE MOD_LinearOperator,          ONLY: MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,     ONLY:Preconditioner
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
REAL             :: AbortCrit,AbortCrit2
INTEGER          :: chance
! preconditioner
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: St(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: tS,tE,tStart,tend
REAL             :: tDG, tPrecond
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n -( I - dt*DG_Operator ) U^n
!    = dt*DG_Operator U^n 
!    = dt * Ut


! store here for later use
CALL CPU_TIME(tS)
Un   = U ! here, n stands for U^n
Uold = U
chance=0
DO WHILE (chance.LT.2)  ! maximum of two trials with BiCGStab inner interation
  ! init and get first error norm
  ! Compute A*U^n
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff * DG_Operator(Un)
  CALL VectorDotProduct(R0,R0,Norm_R0)
  alpha  = Norm_R0
  Norm_R0=SQRT(Norm_R0)
  P=R0
  R=R0
  AbortCrit  = Norm_R0*eps_LinearSolver
  AbortCrit2 = AbortCrit*AbortCrit
  
  DO iter=1,maxIter_LinearSolver
    ! Preconditioner
    CALL CPU_TIME(tStart)
    CALL Preconditioner(coeff,P,Pt)
    CALL CPU_TIME(tend)
    ! matrix vector
    CALL CPU_TIME(tStart)
    CALL MatrixVector(t,coeff,Pt,V)
    CALL CPU_TIME(tend)
    tDG=tDG+tend-tStart
    CALL VectorDotProduct(V,R0,sigma)

!    CALL VectorDotProduct(V,V,Norm_V) 
!    Norm_V = SQRT(Norm_V)
!    AbortCrit2 = epsTilde_LinearSolver*Norm_V*Norm_R0
!    IF((sigma.GT.AbortCrit2).AND.(iter.GT.10)) THEN
      alpha=alpha/sigma
      S = R - alpha*V
      CALL VectorDotProduct(S,S,Norm_S)
      !Norm_S = SQRT(Norm_S)
      !IF((Norm_S.GT.AbortCrit).OR.(Norm_S.GT.1e-12))THEN
      IF((Norm_S.GT.AbortCrit2).OR.(Norm_S.GT.1e-12))THEN
        ! Preconditioner
        CALL CPU_TIME(tStart)
        CALL Preconditioner(coeff,S,St)
        CALL CPU_TIME(tend)
        ! matrix vector
        CALL CPU_TIME(tStart)
        CALL MatrixVector(t,coeff,St,TVec)
        CALL CPU_TIME(tend)
        tDG=tDG+tend-tStart
        CALL VectorDotProduct(TVec,TVec,Norm_T2)
        CALL VectorDotProduct(TVec,S,omega)
        omega=omega/Norm_T2
        ! new un or target variable
        Un=Un+alpha*Pt+omega*St
        R=S-omega*TVec
        CALL VectorDotProduct(R,R0,alpha)
        beta=alpha/(omega*sigma)
        P=R+beta*(P-omega*V)
        CALL VectorDotProduct(R,R,Norm_R)
        Norm_R=SQRT(Norm_R)
        IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12)) THEN
          U=Un 
          totalIterLinearSolver=totalIterLinearSolver+iter
          CALL CPU_TIME(tE)
          ! Debug Ausgabe, Anzahl der Iterationen...
          SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',iter
          SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' T in STABBiCGSTAB  : ',tE-tS
          SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
          SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
          RETURN
        ENDIF ! Norm_R < AbortCrit
      ELSE ! Norm_S .LT. 1e-12
        Un=Un+alpha*Pt
        R=S
        CALL VectorDotProduct(R,R,Norm_R)
        Norm_R=SQRT(Norm_R)
        CALL VectorDotProduct(R,R0,alpha)
        IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12)) THEN
          U=Un
          totalIterLinearSolver=totalIterLinearSolver+iter
          ! Debug Ausgabe, Anzahl der Iterationen...
          CALL CPU_TIME(tE)
          ! Debug Ausgabe, Anzahl der Iterationen...
          SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',iter
          SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' T in STABBiCGSTAB  : ',tE-tS
          SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
          SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
          RETURN
        ENDIF ! Norm_R < AbortCrit
      END IF ! Norm_S
!    ELSE
!      SWRITE(*,*) 'Sigma crit killed us!!'
!      U = Un
!      totalIterLinearSolver=totalIterLinearSolver+iter
!    END IF ! Norm_Sigma
  END DO !iter
 ! restart with new U 
  Un   = 0.5*(Uold+Un)
  Uold = U
  LinSolverRHS = U
  ImplicitSource = 0.
  totalIterLinearSolver=totalIterLinearSolver+iter
  chance = chance+1
  SWRITE(*,*) 'No convergence during first half of iterations.'
  SWRITE(*,*) 'Iter,t,Norm_R,Norm_R0',iter,t,Norm_R,Norm_R0

END DO ! while

SWRITE(*,*)'Norm_R        : ',Norm_R
SWRITE(*,*)'Norm_R0       : ',Norm_R0
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'StabBiCGSTAB_P NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',chance,REAL(iter))

END SUBROUTINE LinearSolver_StabBiCGSTAB_P


SUBROUTINE LinearSolver_GMRES_P(t,coeff)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess   ! why not Un??
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY: U
USE MOD_LinearSolver_Vars,    ONLY: LinSolverRHS,ImplicitSource
USE MOD_LinearSolver_Vars,    ONLY: eps_LinearSolver,TotalIterLinearSolver
USE MOD_LinearSolver_Vars,    ONLY: nKDim,nRestarts,nInnerIter
USE MOD_LinearOperator,       ONLY: MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL              :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:nKDim)
REAL              :: V2P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL              :: W(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL              :: Z(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL              :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL              :: Gam(1:nKDim+1),C(1:nKDim),S(1:nKDim),H(1:nKDim+1,1:nKDim+1),Alp(1:nKDim+1)
REAL              :: Norm_R0,Resu,Temp,Bet
REAL              :: AbortCrit
INTEGER           :: Restart
INTEGER           :: m,nn,o
! preconditoner + Vt
REAL              :: Vt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:nKDim)
REAL              :: Vt2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
#ifdef DLINANALYZE
REAL              :: tS,tE, tS2,tE2,t1,t2
real              :: tstart,tend,tPrecond,tDG
#endif /* DLINANALYZE */
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n -( I - dt*DG_Operator ) U^n
!    = dt*DG_Operator U^n 
!    = dt * Ut

#ifdef DLINANALYZE
! time measurement
CALL CPU_TIME(tS)
! start GMRES
tPrecond=0.
tDG=0.
#endif /* DLINANALYZE */

Restart=0
nInnerIter=0
Un=U 

! compute starting residual 
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
CALL VectorDotProduct(R0,R0,Norm_R0)
Norm_R0=SQRT(Norm_R0)
! define relative abort criteria
AbortCrit=Norm_R0*eps_LinearSolver
! GMRES(m)  inner loop
V(:,:,:,:,:,1)=R0/Norm_R0
Gam(1)=Norm_R0

DO WHILE (Restart<nRestarts)
  DO m=1,nKDim
    nInnerIter=nInnerIter+1
    ! Preconditioner
#ifdef DLINANALYZE
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    CALL Preconditioner(coeff,V(:,:,:,:,:,m),Vt(:,:,:,:,:,m))
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tPrecond=tPrecond+tend-tStart
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    ! matrix vector
    CALL MatrixVector(t,coeff,Vt(:,:,:,:,:,m),W)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tDG=tDG+tend-tStart
#endif /* DLINANALYZE */
    ! Gram-Schmidt
    DO nn=1,m
      CALL VectorDotProduct(V(:,:,:,:,:,nn),W,H(nn,m))
      W=W-H(nn,m)*V(:,:,:,:,:,nn)
    END DO !nn
    CALL VectorDotProduct(W,W,Resu)
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
    IF ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) THEN !converge or max Krylov reached
      DO nn=m,1,-1
         Alp(nn)=Gam(nn) 
         DO o=nn+1,m
           Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
         END DO !o
         Alp(nn)=Alp(nn)/H(nn,nn)
      END DO !nn
      DO nn=1,m
        Un=Un+Alp(nn)*Vt(:,:,:,:,:,nn)
      END DO !nn
      IF (ABS(Gam(m+1)).LE.AbortCrit) THEN !converged
        totalIterLinearSolver=totalIterLinearSolver+nInnerIter
        U=Un
#ifdef DLINANALYZE
        CALL CPU_TIME(tE)
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' nRestarts          : ',Restart
        SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in GMRES      : ',tE-tS
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Gam(1)
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Gam(m+1)
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Ratio Precond/DG   : ',tPrecond/tDG
#endif /* DLINANALYZE */
        RETURN
      END IF  ! converged
    ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
      V(:,:,:,:,:,m+1)=W/H(m+1,m)
    END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
  END DO ! m 
  ! Restart needed
  Restart=Restart+1
  ! new settings for source
  U=Un
  ImplicitSource = 0.
  ! does not change LinSolverRHS
! start residuum berrechnen
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
  CALL VectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
  ! GMRES(m)  inner loop
  V(:,:,:,:,:,1)=R0/Norm_R0
  Gam(1)=Norm_R0
END DO ! Restart

IF(MPIRoot) CALL abort(&
__STAMP__ &
,'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(nInnerIter))

END SUBROUTINE LinearSolver_GMRES_P

SUBROUTINE LinearSolver_BiCGStab_LRP(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab with left and right preconditioners
! left preconditioner is right preconditioner
! alogrithm from Meister, p. 228 with direct residual control for left preconditoned algorithm
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U
USE MOD_TimeDisc_Vars,        ONLY:dt
USE MOD_LinearSolver_Vars,    ONLY:eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver,nInnerIter
USE MOD_LinearSolver_Vars,    ONLY:LinSolverRHS,ImplicitSource,nRestarts
USE MOD_LinearOperator,       ONLY:MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,Norm_R,Norm_T2,Norm_Rn
INTEGER          :: iterLinSolver,Restart
REAL             :: alpha,sigma,omega,beta
REAL             :: AbortCrit
! preconditioner
REAL             :: deltaX(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Rt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0t(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Vt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Tvect(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: St(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n+1 -( I - dt*DG_Operator ) U^n
!    = dt*Source^n+1  + dt*DG_Operator U^n 
!    = dt*Ut + dt*Source^n+1

! store here for later use
Un   = U ! here, n stands for U^n
deltaX = 0.
Uold = U
Restart=0
nInnerIter = 0
! LinSolverRHS and X0 = U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
!P = R0
R = R0
CALL VectorDotProduct(R0,R0,Norm_R0)
AbortCrit = eps_LinearSolver*SQRT(Norm_R0)

! left precondtioning of residuum
CALL Preconditioner(coeff,R0,R0t)

! compute preconditoned R0
Pt = R0t
Rt = R0t
CALL VectorDotProduct(R0t,R0t,Norm_R0)
Norm_R = Norm_R0


DO WHILE (Restart.LT.nRestarts)  ! maximum of two trials with BiCGStab inner interation
  DO iterLinSolver = 1, maxIter_LinearSolver  ! two trials with half of iterations
    ! Preconditioner
    ! right preconditioner before Maxtrix-Vector
    CALL Preconditioner(coeff,Pt,V)
    CALL MatrixVector(t,coeff,V,V) ! or V,Vt and V=Vt
    ! left preconditioner
    CALL Preconditioner(coeff,V,Vt)

    ! compute preconditioned alpha
    CALL VectorDotProduct(Rt,R0t,alpha)
    CALL VectorDotProduct(Vt,R0t,sigma)
    alpha=alpha/sigma
    S = R - alpha*V ! requires to save v 
    ! left preconditoner
    !CALL Preconditioner(coeff,S,St)
    ! or does not require the application of the preconditioner, Meister,p.227,EQ.5.9.3
    St = Rt - alpha*Vt

    ! next Matrix Vector
    ! right precondtioner
    CALL Preconditioner(coeff,St,Tvec) 
    CALL MatrixVector(t,coeff,Tvec,TVec) ! or Tvec,Tvect;  Tvec=TvecT
    ! left preconditioner
    CALL Preconditioner(coeff,Tvec,Tvect)

    ! compute omega
    CALL VectorDotProduct(TVect,TVect,Norm_T2)
    CALL VectorDotProduct(TVect,St,omega)
    omega=omega/Norm_T2

    ! change of x
    deltaX=deltaX+alpha*Pt+omega*St

    ! compute new residuum
    R =S -omega*Tvec ! requires to store Tvec and Tvect
    Rt=St-omega*Tvect

    CALL VectorDotProduct(Rt,R0t,Norm_RN)
    beta=alpha*Norm_RN/(omega*Norm_R)
    Pt = Rt+ beta*(Pt-omega*Vt)

    CALL VectorDotProduct(R,R,Norm_R)
    Norm_R=SQRT(Norm_R)
    ! test if success
    IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12)) THEN
      CALL Preconditioner(coeff,deltaX,deltaX)
      U=Un+deltaX
      nInnerIter=nInnerIter+iterLinSolver
      totalIterLinearSolver=totalIterLinearSolver+nInnerIter
#ifdef DLINANALYZE
      CALL CPU_TIME(tE)
      ! Debug Ausgabe, Anzahl der Iterationen...
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in BiCGSTAB   : ',tE-tS
      SWRITE(UNIT_stdOut,'(A23,E16.8)')   ' Norm_R0            : ',Norm_R0
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Ratio Precond/DG   : ',tPrecond/tDG
#endif /* DLINANALYZE */
      RETURN
    ENDIF
    Norm_R = Norm_RN
  END DO ! iterLinSolver
  CALL Preconditioner(coeff,deltaX,deltaX)
  U=Un+deltaX
  deltaX = 0.
  Uold = U

  ! restart with new U
  ImplicitSource = 0.
  ! LinSolverRHS and X0 = U
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
  R = R0
  ! left precondtioning of residuum
  CALL Preconditioner(coeff,R0,R0t)
  ! compute preconditoned R0
  Pt = R0t
  Rt = R0t
  CALL VectorDotProduct(R0t,R0t,Norm_R0)
  Norm_R = Norm_R0
  nInnerIter=nInnerIter+iterLinSolver
  Restart = Restart + 1
END DO ! while chance < 2 

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'BiCGSTAB NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',Restart,REAL(nInnerIter+iterLinSolver))

END SUBROUTINE LinearSolver_BiCGSTAB_LRP

SUBROUTINE LinearSolver_BiCGStab_LP(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab with left and right preconditioners
! left preconditioner is right preconditioner
! alogrithm from Meister, p. 228 with direct residual control for left preconditoned algorithm
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U
USE MOD_TimeDisc_Vars,        ONLY:dt
USE MOD_LinearSolver_Vars,    ONLY:eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver,nInnerIter
USE MOD_LinearSolver_Vars,    ONLY:LinSolverRHS,ImplicitSource,nRestarts
USE MOD_LinearOperator,       ONLY:MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: UOld(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: S(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: TVec(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Norm_R0,Norm_R,Norm_T2,Norm_Rn
INTEGER          :: iterLinSolver,Restart
REAL             :: alpha,sigma,omega,beta
REAL             :: AbortCrit
! preconditioner
!REAL             :: deltaX(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Rt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: R0t(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Vt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Tvect(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: St(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n+1 -( I - dt*DG_Operator ) U^n
!    = dt*Source^n+1  + dt*DG_Operator U^n 
!    = dt*Ut + dt*Source^n+1

! store here for later use
Un   = U ! here, n stands for U^n
Uold = U
Restart=0
nInnerIter = 0
! LinSolverRHS and X0 = U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
!P = R0
R = R0
CALL VectorDotProduct(R0,R0,Norm_R0)
AbortCrit = eps_LinearSolver*SQRT(Norm_R0)

! left precondtioning of residuum
CALL Preconditioner(coeff,R0,R0t)

! compute preconditoned R0
Pt = R0t
Rt = R0t
CALL VectorDotProduct(R0t,R0t,Norm_R0)
Norm_R = Norm_R0


DO WHILE (Restart.LT.nRestarts)  ! maximum of two trials with BiCGStab inner interation
  DO iterLinSolver = 1, maxIter_LinearSolver  ! two trials with half of iterations
    ! Preconditioner
    CALL MatrixVector(t,coeff,Pt,V) ! or V,Vt and V=Vt
    ! left preconditioner
    CALL Preconditioner(coeff,V,Vt)

    ! compute preconditioned alpha
    CALL VectorDotProduct(Rt,R0t,alpha)
    CALL VectorDotProduct(Vt,R0t,sigma)
    alpha=alpha/sigma
    S = R - alpha*V ! requires to save v 
    ! left preconditoner
    !CALL Preconditioner(coeff,S,St)
    ! or does not require the application of the preconditioner, Meister,p.227,EQ.5.9.3
    St = Rt - alpha*Vt

    ! next Matrix Vector
    CALL MatrixVector(t,coeff,St,TVec) ! or Tvec,Tvect;  Tvec=TvecT
    ! left preconditioner
    CALL Preconditioner(coeff,Tvec,Tvect)

    ! compute omega
    CALL VectorDotProduct(TVect,TVect,Norm_T2)
    CALL VectorDotProduct(TVect,St,omega)
    omega=omega/Norm_T2

    ! change of x
    Un=Un+alpha*Pt+omega*St

    ! compute new residuum
    R =S -omega*Tvec ! requires to store Tvec and Tvect
    Rt=St-omega*Tvect

    CALL VectorDotProduct(Rt,R0t,Norm_RN)
    beta=alpha*Norm_RN/(omega*Norm_R)
    Pt = Rt+ beta*(Pt-omega*Vt)

    CALL VectorDotProduct(R,R,Norm_R)
    Norm_R=SQRT(Norm_R)
    ! test if success
    IF((Norm_R.LE.AbortCrit).OR.(Norm_R.LT.1.E-12)) THEN
      U=Un
      nInnerIter=nInnerIter+iterLinSolver
      totalIterLinearSolver=totalIterLinearSolver+nInnerIter
#ifdef DLINANALYZE
      CALL CPU_TIME(tE)
      ! Debug Ausgabe, Anzahl der Iterationen...
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in BiCGSTAB   : ',tE-tS
      SWRITE(UNIT_stdOut,'(A23,E16.8)')   ' Norm_R0            : ',Norm_R0
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
      SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Ratio Precond/DG   : ',tPrecond/tDG
#endif /* DLINANALYZE */
      RETURN
    ENDIF
    Norm_R = Norm_RN
  END DO ! iterLinSolver
  U=Un
  Uold = U

  ! restart with new U
  ImplicitSource = 0.
  ! LinSolverRHS and X0 = U
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output
  R = R0
  ! left precondtioning of residuum
  CALL Preconditioner(coeff,R0,R0t)
  ! compute preconditoned R0
  Pt = R0t
  Rt = R0t
  CALL VectorDotProduct(R0t,R0t,Norm_R0)
  Norm_R = Norm_R0
  nInnerIter=nInnerIter+iterLinSolver
  Restart = Restart + 1
END DO ! while chance < 2 

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'BiCGSTAB NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',Restart,REAL(nInnerIter+iterLinSolver))

END SUBROUTINE LinearSolver_BiCGSTAB_LP

SUBROUTINE LinearSolver_BiCGSTABl(t,Coeff)
!===================================================================================================================================
! Solves Linear system Ax=b using BiCGStab(l) with right preconditioner P_r
! Matrix A = I - Coeff*R
! Attention: Vector x is U^n+1, initial guess set to U^n 
! Attention: Vector b is U^n 
! Sleijpen 1993: BiCGSTAB(l) for linear equations involving unsymmetric matrices with complex spectrum
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U
USE MOD_TimeDisc_Vars,        ONLY:dt
USE MOD_LinearSolver_Vars,    ONLY:eps_LinearSolver,maxIter_LinearSolver,totalIterLinearSolver,nInnerIter
USE MOD_LinearSolver_Vars,    ONLY:LinSolverRHS,ImplicitSource,nRestarts,ldim
USE MOD_LinearOperator,       ONLY:MatrixVector, MatrixVectorSource, VectorDotProduct
USE MOD_ApplyPreconditioner,  ONLY:Preconditioner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: deltaX(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: P(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,0:ldim)
REAL             :: R(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,0:ldim)
REAL             :: R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: sigma(0:ldim),tau(1:ldim,1:ldim)
REAL             :: phi(0:ldim),phis(0:ldim),phiss(0:ldim)
INTEGER          :: iterLinSolver,Restart
INTEGER          :: m,nn
REAL             :: alpha,omega,beta
REAL             :: Norm_R, Norm_R0, Norm_Abort
REAL             :: AbortCrit
! preconditioner
REAL             :: Pt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL             :: Rt(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!===================================================================================================================================

! U^n+1 = U^n + dt * DG_Operator U^n+1 + Sources^n+1
! (I - dt*DG_Operator) U^n+1 = U^n + dt*Sources^n+1
!       A                x   = b
! 
! Residuum
! for initial guess, x0 is set to U^n
! R0 = b - A x0
!    = U^n + dt*Sources^n+1 -( I - dt*DG_Operator ) U^n
!    = dt*Source^n+1  + dt*DG_Operator U^n 
!    = dt*Ut + dt*Source^n+1

! store here for later use
Un   = U ! here, n stands for U^n
Restart=0
nInnerIter = 0
deltaX=0.0d0
! LinSolverRHS and X0 = U
CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output

! compute  A*U^n
CALL VectorDotProduct(R0,R0,Norm_R0)
Norm_R0=SQRT(Norm_R0)
AbortCrit = Norm_R0*eps_LinearSolver

! starting direction accoring to old paper
P(:,:,:,:,:,0) = 0.
R(:,:,:,:,:,0) = R0
Norm_R0=1.0d0
alpha=0.0d0
omega=1.0d0


DO WHILE(Restart.LT.nRestarts)
  DO iterLinSolver=1,maxIter_LinearSolver
    Norm_R0= - omega* Norm_R0
    ! Bi-CG Part
    DO nn=0, ldim-1
      CALL VectorDotProduct(R(:,:,:,:,:,nn),R0,Norm_R)
      beta=alpha*Norm_R/Norm_R0
      Norm_R0=Norm_R
      DO m=0,nn
        P(:,:,:,:,:,m) = R(:,:,:,:,:,m) - beta * P(:,:,:,:,:,m)
      END DO ! m
      ! Preconditioner
      CALL Preconditioner(coeff,P(:,:,:,:,:,nn),Pt)
      ! matrix vector
      CALL MatrixVector(t,coeff,Pt,P(:,:,:,:,:,nn+1))
      CALL VectorDotProduct(P(:,:,:,:,:,nn+1),R0,phi(nn))
      alpha=Norm_R0/phi(nn)
      DO m=0,nn
        R(:,:,:,:,:,m) = R(:,:,:,:,:,m) - alpha * P(:,:,:,:,:,m+1)
      END DO ! m
      ! Preconditioner
      CALL Preconditioner(coeff,R(:,:,:,:,:,nn),Rt)
      ! matrix vector
      CALL MatrixVector(t,coeff,Rt,R(:,:,:,:,:,nn+1))
      deltaX=deltaX+alpha*P(:,:,:,:,:,0)
    END DO ! nn
    ! mod. G.-S.
    DO nn=1, ldim
      DO m=1,nn-1
        CALL VectorDotProduct(R(:,:,:,:,:,nn),R(:,:,:,:,:,m),tau(m,nn))
        tau(m,nn) = 1.0d0/sigma(m) * tau(m,nn)
        R(:,:,:,:,:,nn) = R(:,:,:,:,:,nn) - tau(m,nn) * R(:,:,:,:,:,m)
      END DO ! m
      CALL VectorDotProduct(R(:,:,:,:,:,nn),R(:,:,:,:,:,nn),sigma(nn))
      CALL VectorDotProduct(R(:,:,:,:,:,nn),R(:,:,:,:,:,0),phis(nn))
      phis(nn) = phis(nn)/sigma(nn)
    END DO ! nn
    phi(ldim) = phis(ldim)
    omega=phi(ldim)
    DO nn=ldim-1,1,-1
      phi(nn) = phis(nn)
      DO m=nn+1,ldim
        phi(nn) = phi(nn)-tau(nn,m)*phi(m)
      END DO ! m
    END DO ! nn
    DO nn=1,ldim-1
      phiss(nn) = phi(nn+1)
      DO m=nn+1,ldim-1
        phiss(nn)=phiss(nn) + tau(nn,m)*phi(m+1)
      END DO ! m
    END DO !  nn
    ! update
    deltaX = deltaX+phi(1)*R(:,:,:,:,:,0)
    R(:,:,:,:,:,0) = R(:,:,:,:,:,0) - phis(ldim)*R(:,:,:,:,:,ldim)
    P(:,:,:,:,:,0) = P(:,:,:,:,:,0) - phi (ldim)*P(:,:,:,:,:,ldim)
    DO nn=1,ldim-1
      P(:,:,:,:,:,0) = P(:,:,:,:,:,0) - phi (nn)*P(:,:,:,:,:,nn)
      deltaX = deltaX+phiss(nn)*R(:,:,:,:,:,nn)
      R(:,:,:,:,:,0) = R(:,:,:,:,:,0) - phis(nn)*R(:,:,:,:,:,nn)
    END DO ! nn
    CALL VectorDotProduct(R(:,:,:,:,:,0),R(:,:,:,:,:,0),Norm_Abort)
    Norm_Abort=SQRT(Norm_Abort)
    IF((Norm_Abort.LE.AbortCrit).OR.(Norm_Abort.LT.1.E-12)) THEN
      ! invert preconditioner
      CALL Preconditioner(coeff,deltaX,U)
      U=U+Un
      nInnerIter=nInnerIter+iterLinSolver*ldim
      totalIterLinearSolver=totalIterLinearSolver+nInnerIter
      !SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter LinSolver     : ',nInnerIter
      !SWRITE(UNIT_stdOut,'(A22,I5)')      ' Restarts           : ',Restart
      RETURN
    END IF
  END DO ! iterLinearSolver
  ! restart with new U
  ! LinSolverRHS and X0 = U
!  U              = 0.5*(Uold+Un)
  CALL Preconditioner(coeff,deltaX,U)
  U=U+Un
  ImplicitSource = 0.
  ! LinSolverRHS and X0 = U
  CALL MatrixVectorSource(t,Coeff,R0) ! coeff*Ut+Source^n+1 ! only output

  ! starting direction accoring to old paper
  P(:,:,:,:,:,0) = 0.
  R(:,:,:,:,:,0) = R0
  Norm_R0=1.
  alpha=0.
  omega=1.
  nInnerIter=nInnerIter+iterLinSolver*ldim
  Restart = Restart + 1
END DO ! while chance < 2 

SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Norm_R0
SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Norm_R
IF(MPIRoot) CALL abort(&
__STAMP__ &
,'BiCGSTAB(l) NOT CONVERGED WITH RESTARTS AND BiCGSTAB ITERATIONS:',Restart,REAL(nInnerIter+iterLinSolver))

END SUBROUTINE LinearSolver_BiCGSTABl

SUBROUTINE FinalizeLinearSolver()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_LinearSolver_Vars,ONLY:LinearSolverInitIsDone,ImplicitSource,LinSolverRHS
USE MOD_Predictor    ,ONLY:FinalizePredictor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

LinearSolverInitIsDone = .FALSE.
SDEALLOCATE(ImplicitSource)
SDEALLOCATE(LinSolverRHS)
CALL FinalizePredictor
!SDEALLOCATE(FieldSource)
END SUBROUTINE FinalizeLinearSolver


END MODULE MOD_LinearSolver
