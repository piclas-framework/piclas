#include "boltzplatz.h"

MODULE MOD_LinearOperator
!===================================================================================================================================
! Module containing matrix vector operations
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

#ifndef PP_HDG
INTERFACE MatrixVector
  MODULE PROCEDURE MatrixVector
END INTERFACE

INTERFACE MatrixVectorSource
  MODULE PROCEDURE MatrixVectorSource
END INTERFACE

INTERFACE VectorDotProduct
  MODULE PROCEDURE VectorDotProduct
END INTERFACE
#endif /*NOT HDG*/

#if defined(PARTICLES) && defined(IMPA)
INTERFACE PartVectorDotProduct
  MODULE PROCEDURE PartVectorDotProduct
END INTERFACE

INTERFACE PartMatrixVector
  MODULE PROCEDURE PartMatrixVector
END INTERFACE
#endif


#ifndef PP_HDG
PUBLIC:: MatrixVector, MatrixVectorSource, VectorDotProduct, ElementVectorDotProduct, DENSE_MATMUL
#endif /*NOT HDG*/
#if defined(PARTICLES) && defined(IMPA)
PUBLIC:: PartVectorDotProduct,PartMatrixVector
#endif
!===================================================================================================================================

CONTAINS

#ifndef PP_HDG
SUBROUTINE MatrixVector(t,Coeff,X,Y)
!===================================================================================================================================
! Computes Matrix Vector Product y=A*x for linear Equations only
! Matrix A = I - Coeff*R
! Attention: Vector x is always U 
! Attention: Result y is always stored in Ut
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY:U,Ut
USE MOD_DG,                 ONLY:DGTimeDerivative_weakForm
USE MOD_LinearSolver_Vars,  ONLY:mass
USE MOD_Equation_Vars,      ONLY:DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars,      ONLY:sdtCFLOne
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
REAL,INTENT(IN)  :: X(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Y(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: rTmp
REAL             :: locMass
INTEGER          :: i,j,k,iElem,iVar
!===================================================================================================================================
U=X
CALL DGTimeDerivative_weakForm(t,t,0,doSource=.FALSE.)
! y = (I-Coeff*R)*x = x - Coeff*R*x 
IF(DoParabolicDamping)THEN
  !rTmp=1.0-(fDamping-1.0)*dt*sdTCFLOne
  rTmp=1.0-(fDamping-1.0)*coeff*sdTCFLOne
  DO iElem=1,PP_nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          locMass=mass(1,i,j,k,iElem)
          DO iVar=1,6
            Y(iVar,i,j,k,iElem) = locMass*(     U(iVar,i,j,k,iElem) - Coeff*Ut(iVar,i,j,k,iElem))
          END DO ! iVar=1,6
          DO iVar=7,PP_nVar
            Y(iVar,i,j,k,iElem) = locMass*(rTmp*U(iVar,i,j,k,iElem) - Coeff*Ut(iVar,i,j,k,iElem))
          END DO ! iVar=7,PP_nVar
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO ! iElem=1,PP_nElems
ELSE
  Y = mass*(U - Coeff*Ut)
END IF
END SUBROUTINE MatrixVector


SUBROUTINE MatrixVectorSource(t,Coeff,Y)
!===================================================================================================================================
! Computes Matrix Vector Product y=A*x for linear Equations only
! Matrix A = I - Coeff*R
! Attention: Vector x is always U 
! Attention: Result y is always stored in Ut
! Attention: Is only Required for the calculation of the residuum
!            THEREFORE: coeff*DG_Operator(U) is the output pf this routine
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,           ONLY:U,Ut
USE MOD_DG,                ONLY:DGTimeDerivative_weakForm
USE MOD_Equation,          ONLY:CalcSource
USE MOD_Equation,          ONLY:DivCleaningDamping
USE MOD_LinearSolver_Vars, ONLY:ImplicitSource, LinSolverRHS,mass
USE MOD_Equation_Vars,     ONLY:DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars,     ONLY:sdtCFLOne
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Y(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: rTmp(1:8),locMass
INTEGER          :: i,j,k,iElem,iVar
!===================================================================================================================================

! y =  Coeff*Ut+source
! Residual = b - A x0 
!          = b - x0 + coeff*ut + coeff*Source^n+1

CALL DGTimeDerivative_weakForm(t,t,0,doSource=.FALSE.)
!Y = LinSolverRHS - X0 +coeff*ut
#ifndef PP_HDG
CALL CalcSource(t,1.0,ImplicitSource)
#endif

IF(DoParabolicDamping)THEN
  rTmp(1:6)=1.0
  rTmp( 7 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
  rTmp( 8 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
ELSE
  rTmp(1:8)=1.0
END IF

DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        locMass=mass(1,i,j,k,iElem)
        DO iVar=1,PP_nVar
          Y(iVar,i,j,k,iElem) = locMass*( LinSolverRHS(iVar,i,j,k,iElem)         &
                                         -rTmp(iVar)*U(iVar,i,j,k,iElem)         &
                                         +    coeff*Ut(iVar,i,j,k,iElem)         &
                                         +coeff*ImplicitSource(iVar,i,j,k,iElem) )
        END DO ! iVar=1,PP_nVar
      END DO ! i=0,PP_N
    END DO ! j=0,PP_N
  END DO ! k=0,PP_N
END DO ! iElem=1,PP_nElems

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
#ifdef MPI
REAL              :: ResuSend
#endif
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
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

END SUBROUTINE VectorDotProduct


SUBROUTINE ElementVectorDotProduct(a,b,resu)
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
REAL,INTENT(IN)   :: a(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: b(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar,i,j,k
!===================================================================================================================================

resu=0.
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      DO iVar=1,PP_nVar
        resu=resu + a(iVar,i,j,k)*b(iVar,i,j,k)
      END DO
    END DO
  END DO
END DO

END SUBROUTINE ElementVectorDotProduct

SUBROUTINE DENSE_MATMUL(n,A,x,y)
!===================================================================================================================================
! Computes a dense Matrix vector product
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: A(n,n), x(n)
INTEGER,INTENT(IN)    :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: y(n) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER               :: i,j
!===================================================================================================================================

DO i=1,n
 y(i)=0.
END DO

DO j=1,n
  DO i=1,n
    y(i) = y(i)+A(i,j)*x(j)
  END DO ! i
END DO ! j

END SUBROUTINE DENSE_MATMUL
#endif /*NOT HDG*/


#if defined(PARTICLES) && defined(IMPA)
SUBROUTINE PartVectorDotProduct(a,b,resu)
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
REAL,INTENT(IN)   :: a(1:6)
REAL,INTENT(IN)   :: b(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar
!===================================================================================================================================

resu=0.
DO iVar=1,6
  resu=resu + a(iVar)*b(iVar)
END DO

END SUBROUTINE PartVectorDotProduct


SUBROUTINE PartMatrixVector(t,Coeff,PartID,X,Y)
!===================================================================================================================================
! Computes Matrix Vector Product y=A*x for linear Equations only
! Matrix A = I - Coeff*R
! Attention: Vector x is always U 
! Attention: Result y is always stored in Ut
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,                 ONLY:Abort
USE MOD_LinearSolver_Vars,       ONLY:reps0,PartXK,R_PartXK
USE MOD_Equation_Vars,           ONLY:c2_inv
USE MOD_Particle_Vars,           ONLY:PartState, PartLorentzType
USE MOD_Part_RHS,                ONLY:SLOW_RELATIVISTIC_PUSH,FAST_RELATIVISTIC_PUSH &
                                     ,RELATIVISTIC_PUSH,NON_RELATIVISTIC_PUSH
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t,Coeff
INTEGER,INTENT(IN) :: PartID
REAL,INTENT(IN)    :: X(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Y(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: X_abs,epsFD
REAL               :: PartT(1:6)
REAL               :: LorentzFacInv
!REAL               :: FieldAtParticle(1:6)
!REAL               :: typ_v_abs, XK_V, sign_XK_V
!===================================================================================================================================

CALL PartVectorDotProduct(X,X,X_abs)
IF(X_abs.NE.0.)THEN
  EpsFD= rEps0/SQRT(X_abs)
ELSE
  EpsFD= rEps0*0.1
END IF

!CALL PartVectorDotProduct(PartState(PartID,1:6),ABS(X),typ_v_abs)
!CALL PartVectorDotProduct(PartXK(1:6,PartID),X,Xk_V)
!sign_XK_V=SIGN(1.,Xk_V)
!EpsFD= rEps0/X_abs*MAX(ABS(Xk_V),typ_v_abs)*SIGN_Xk_V

PartState(PartID,1:6) = PartXK(1:6,PartID)+EpsFD*X
! compute fields at particle position, if relaxation freez, therefore use fixed field and pt
!CALL InterpolateFieldToSingleParticle(PartID,FieldAtParticle)
!PartT(4:6)=Pt(PartID,1:3)
SELECT CASE(PartLorentzType)
CASE(0)
  PartT(4:6) = NON_RELATIVISTIC_PUSH(PartID,FieldAtParticle(PartID,1:6))
  LorentzFacInv = 1.0
CASE(1)
  PartT(4:6) = SLOW_RELATIVISTIC_PUSH(PartID,FieldAtParticle(PartID,1:6))
  LorentzFacInv = 1.0
CASE(3)
  PartT(4:6) = FAST_RELATIVISTIC_PUSH(PartID,FieldAtParticle(PartID,1:6))
  LorentzFacInv = 1.0
CASE(5)
  LorentzFacInv=1.0+DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))*c2_inv      
  LorentzFacInv=1.0/SQRT(LorentzFacInv)
  PartT(4:6) = RELATIVISTIC_PUSH(PartID,FieldAtParticle(PartID,1:6),LorentzFacInvIn=LorentzFacInv)
CASE DEFAULT
CALL abort(&
__STAMP__ &
,' Given PartLorentzType does not exist!',PartLorentzType)
END SELECT
PartT(1)=LorentzFacInv*PartState(PartID,4) ! funny, or PartXK
PartT(2)=LorentzFacInv*PartState(PartID,5) ! funny, or PartXK
PartT(3)=LorentzFacInv*PartState(PartID,6) ! funny, or PartXK
! or frozen version
Y = (X - (coeff/EpsFD)*(PartT - R_PartXk(:,PartID)))

! compiler warnings
IF(1.EQ.2)THEN
  PartT(1)=t
END IF

END SUBROUTINE PartMatrixVector
#endif


END MODULE MOD_LinearOperator
