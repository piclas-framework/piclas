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

#if !(USE_HDG)
INTERFACE MatrixVector
  MODULE PROCEDURE MatrixVector
END INTERFACE

INTERFACE MatrixVectorSource
  MODULE PROCEDURE MatrixVectorSource
END INTERFACE

INTERFACE VectorDotProduct
  MODULE PROCEDURE VectorDotProduct
END INTERFACE

#ifdef IMPA
INTERFACE EvalResidual
  MODULE PROCEDURE EvalResidual
END INTERFACE
#endif
#endif /*NOT HDG*/

#if defined(PARTICLES)
#if defined(IMPA) || defined(ROS)
INTERFACE PartVectorDotProduct
  MODULE PROCEDURE PartVectorDotProduct
END INTERFACE

INTERFACE PartMatrixVector
  MODULE PROCEDURE PartMatrixVector
END INTERFACE
#endif /*ROS OR IMPA*/
#endif /*PARTICLES*/

#if !(USE_HDG)
PUBLIC :: MatrixVector, MatrixVectorSource, VectorDotProduct, ElementVectorDotProduct, DENSE_MATMUL
#ifdef IMPA
PUBLIC :: EvalResidual
#endif /*IMPA*/
#endif /*NOT HDG*/
#if defined(PARTICLES)
#if defined(IMPA) || defined(ROS)
PUBLIC:: PartVectorDotProduct,PartMatrixVector
#endif /*ROS OR IMPA*/
#endif
!===================================================================================================================================

CONTAINS

#if !(USE_HDG)
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
#ifdef IMPA
  rTmp=1.0-(fDamping-1.0)*coeff*sdTCFLOne
#else
  ! use the inverse of coefficient because of inverse definition
  rTmp=(1.0-(fDamping-1.0)*1./coeff*sdTCFLOne)
#endif
  DO iElem=1,PP_nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          locMass=mass(1,i,j,k,iElem)
#ifdef IMPA
          DO iVar=1,6
            Y(iVar,i,j,k,iElem) = locMass*(     U(iVar,i,j,k,iElem) - Coeff*Ut(iVar,i,j,k,iElem))
          END DO ! iVar=1,6
          DO iVar=7,PP_nVar
            Y(iVar,i,j,k,iElem) = locMass*(rTmp*U(iVar,i,j,k,iElem) - Coeff*Ut(iVar,i,j,k,iElem))
          END DO ! iVar=7,PP_nVar
#else
          ! Rosenbrock, CAUTION: coeff = coeff_inv
          DO iVar=1,6
            Y(iVar,i,j,k,iElem) = locMass*(     coeff*U(iVar,i,j,k,iElem) - Ut(iVar,i,j,k,iElem))
          END DO ! iVar=1,6
          DO iVar=7,PP_nVar
            Y(iVar,i,j,k,iElem) = locMass*(rTmp*coeff*U(iVar,i,j,k,iElem) - Ut(iVar,i,j,k,iElem))
          END DO ! iVar=7,PP_nVar
#endif
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO ! iElem=1,PP_nElems
ELSE
#ifdef IMPA
  Y = mass*(U - Coeff*Ut)
#else
  Y = mass*(coeff* U - Ut)
#endif
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
USE MOD_LinearSolver_Vars, ONLY:LinSolverRHS,mass
USE MOD_Equation_Vars,     ONLY:DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars,     ONLY:sdtCFLOne
#if !(USE_HDG) && !defined(ROS)
USE MOD_LinearSolver_Vars, ONLY:ImplicitSource
#endif /*NO ROSENBROCK and no HDG*/
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
#if !(USE_HDG) && !defined(ROS)
CALL CalcSource(t,1.0,ImplicitSource)
#endif /*NO ROSENBROCK and no HDG*/

IF(DoParabolicDamping)THEN
  rTmp(1:6)=1.0
#ifdef IMPA
  rTmp( 7 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
  rTmp( 8 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
#else
  rTmp(1:6)=1.0
  ! Here, we have to use the inverse because the coeff is the inverse
  rTmp( 7 )=(1.0-(fDamping-1.0)/coeff*sdTCFLOne)
  rTmp( 8 )=(1.0-(fDamping-1.0)/coeff*sdTCFLOne)
#endif
ELSE
  rTmp(1:8)=1.0
END IF

DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        locMass=mass(1,i,j,k,iElem)
        DO iVar=1,PP_nVar
#if ROS
          ! matrix vector for rosenbrock-type RK.
          Y(iVar,i,j,k,iElem) = locMass*( LinSolverRHS(iVar,i,j,k,iElem)         &
                                         -rTmp(iVar)*coeff*U(iVar,i,j,k,iElem)   &
                                         +                 Ut(iVar,i,j,k,iElem)  )
#else
          ! non-rosenbrock RK
          Y(iVar,i,j,k,iElem) = locMass*( LinSolverRHS(iVar,i,j,k,iElem)         &
                                         -rTmp(iVar)*U(iVar,i,j,k,iElem)         &
                                         +    coeff*Ut(iVar,i,j,k,iElem)         &
                                         +coeff*ImplicitSource(iVar,i,j,k,iElem) )
#endif
        END DO ! iVar=1,PP_nVar
      END DO ! i=0,PP_N
    END DO ! j=0,PP_N
  END DO ! k=0,PP_N
END DO ! iElem=1,PP_nElems

END SUBROUTINE MatrixVectorSource

#ifdef IMPA
SUBROUTINE EvalResidual(t,Coeff,Norm_R0)
!===================================================================================================================================
! Compute Initial norm for linear solver by calling MatrixVectorSource and VectorDotProduct
! Matrix A = I - Coeff*R
! Attention: Vector x is always U
! Attention: Result y is always stored in Ut
! Attention: Is only Required for the calculation of the residuum
!            THEREFORE: coeff*DG_Operator(U) is the output pf this routine
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: t,Coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Norm_R0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: Y(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!===================================================================================================================================

! y =  Coeff*Ut+source
! Residual = b - A x0
!          = b - x0 + coeff*ut + coeff*Source^n+1
CALL MatrixVectorSource(t,Coeff,Y)
CALL VectorDotProduct(Y,Y,Norm_R0)
Norm_R0=SQRT(Norm_R0)

END SUBROUTINE EvalResidual
#endif


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
#if USE_MPI
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

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
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


#if defined(PARTICLES)
#if defined(IMPA) || defined(ROS)
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
USE MOD_Globals               ,ONLY: Abort
USE MOD_LinearSolver_Vars     ,ONLY: reps0,PartXK,R_PartXK
USE MOD_Globals_Vars          ,ONLY: c2_inv
USE MOD_Particle_Vars         ,ONLY: PartState, PartLorentzType
USE MOD_Part_RHS              ,ONLY: PartRHS
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle
USE MOD_PICInterpolation      ,ONLY: InterpolateFieldToSingleParticle
USE MOD_Particle_Vars         ,ONLY: PartState
#ifndef ROS
USE MOD_Particle_Vars         ,ONLY: PartDtFrac
#endif
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
!===================================================================================================================================

CALL PartVectorDotProduct(X,X,X_abs)
IF(X_abs.NE.0.)THEN
  EpsFD= rEps0/SQRT(X_abs)
ELSE
  EpsFD= rEps0*0.1
END IF

! CALL PartVectorDotProduct(X,X,X_abs)
! IF(X_abs.NE.0.)THEN
!   CALL PartVectorDotProduct(PartState(1:6,PartID),ABS(X),typ_v_abs)
!   CALL PartVectorDotProduct(PartXK(1:6,PartID),X,Xk_V)
!   sign_XK_V=SIGN(1.,Xk_V)
!   EpsFD= rEps0/X_abs*MAX(ABS(Xk_V),typ_v_abs)*SIGN_Xk_V
! ELSE
!   EpsFD= rEps0*0.1
! END IF

! ! Bassi gang
! CALL PartVectorDotProduct(PartXK(1:6,PartID),PartXK(1:6,PartID),Xk_v)
! CALL PartVectorDotProduct(X,X,X_abs)
! EpsFD=rEps0  * SQRT(1.+XK_V)/X_abs

PartState(1:6,PartID) = PartXK(1:6,PartID)+EpsFD*X
! compute fields at particle position, if relaxation freez, therefore use fixed field and pt
! CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),PEM%GlobalElemID(PartID))
! CALL InterpolateFieldToSingleParticle(PartID,FieldAtParticle(1:6,PartID))
!PartT(4:6)=Pt(1:3,PartID)
IF(PartLorentzType.EQ.5)THEN
  LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv)
  CALL PartRHS(PartID,FieldAtParticle(1:6,PartID),PartT(4:6),LorentzFacInv)
ELSE
  LorentzFacInv = 1.0
  CALL PartRHS(PartID,FieldAtParticle(1:6,PartID),PartT(4:6))
END IF ! PartLorentzType.EQ.5
PartT(1:3)=LorentzFacInv*PartState(4:6,PartID) ! funny, or PartXK
! or frozen version
#if ROS
!Y(1:6) = (Coeff*X(1:6) - (1./EpsFD)*(PartT(1:6) - R_PartXk(1:6,PartID)))
Y(1:6) = (X(1:6) - (coeff/EpsFD)*(PartT(1:6) - R_PartXk(1:6,PartID)))
#else
Y(1:6) = (X(1:6) - (PartDtFrac(PartID)*coeff/EpsFD)*(PartT(1:6) - R_PartXk(1:6,PartID)))
#endif

! compiler warnings
IF(1.EQ.2)THEN
  PartT(1)=t
END IF

END SUBROUTINE PartMatrixVector
#endif /*IMPA or ROS*/
#endif /*PARTICLES*/


END MODULE MOD_LinearOperator
