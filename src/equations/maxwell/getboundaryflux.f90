#include "boltzplatz.h"

MODULE MOD_GetBoundaryFlux
!===================================================================================================================================
! Contains FillBoundary (which depends on the considered equation)
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
INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

PUBLIC::GetBoundaryFlux
!===================================================================================================================================

CONTAINS



SUBROUTINE GetBoundaryFlux(F_Face,BCType,BCState,xGP_Face,normal,t,tDeriv,U_Face)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals,        ONLY:Abort,CROSS
USE MOD_PreProc
USE MOD_Riemann,        ONLY:Riemann
USE MOD_Equation,       ONLY:ExactFunc
USE MOD_Equation_vars,  ONLY:c,c_inv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: BCType,BCState
REAL,INTENT(INOUT)                   :: U_Face(PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: normal(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: t
INTEGER,INTENT(IN)                   :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: F_Face(PP_nVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q
REAL                                 :: n_loc(3),resul(PP_nVar),epsBC
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
SELECT CASE(BCType)
CASE(1) !Periodic already filled!
CASE(2) ! exact BC = Dirichlet BC !!
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      CALL ExactFunc(BCState,t,tDeriv,xGP_Face(:,p,q),U_Face_loc(:,p,q))
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(3) ! 1st order absorbing BC 
        ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  epsBC=1e-10
  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero. 
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,PP_N
    DO p=0,PP_N
      IF (SUM(abs(U_Face(4:6,p,q))).GT.epsBC)THEN
        U_Face_loc(7,p,q) = - U_Face(7,p,q) - c*(DOT_PRODUCT(U_Face(4:6,p,q),normal(1:3,p,q)))
        U_Face_loc(8,p,q) = - U_Face(8,p,q) - c_inv*(DOT_PRODUCT(U_Face(1:3,p,q),normal(1:3,p,q)))
      END IF ! sum(abs(B)) > epsBC
    END DO ! p
  END DO ! q

  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))
  
CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      resul=U_Face(:,p,q)
      n_loc=normal(:,p,q)
      U_Face_loc(1:3,p,q) = -resul(1:3) + 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
      U_Face_loc(4:6,p,q) =  resul(4:6) - 2*(DOT_PRODUCT(resul(4:6),n_loc))*n_loc
      U_Face_loc(  7,p,q) =  resul(  7)
      U_Face_loc(  8,p,q) = -resul(  8)
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))
 
CASE(5) ! 1st order absorbing BC 
        ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero. 
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,PP_N
    DO p=0,PP_N
      U_Face_loc(7,p,q) = - U_Face(7,p,q) - c*(DOT_PRODUCT(U_Face(4:6,p,q),normal(1:3,p,q)))
      U_Face_loc(8,p,q) = - U_Face(8,p,q) - c_inv*(DOT_PRODUCT(U_Face(1:3,p,q),normal(1:3,p,q)))
    END DO ! p
  END DO ! q

  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(6) ! 1st order absorbing BC + fix for low B field
        ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero. 
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,PP_N
    DO p=0,PP_N
      IF (DOT_PRODUCT(U_FACE(4:6,p,q),U_FACE(4:6,p,q))*c*10..GT.DOT_PRODUCT(U_FACE(1:3,p,q),U_FACE(1:3,p,q)))THEN
        U_Face_loc(7,p,q) = - U_Face(7,p,q) - c*(DOT_PRODUCT(U_Face(4:6,p,q),normal(1:3,p,q)))
        U_Face_loc(8,p,q) = - U_Face(8,p,q) - c_inv*(DOT_PRODUCT(U_Face(1:3,p,q),normal(1:3,p,q)))
      END IF ! sum(abs(B)) > epsBC
    END DO ! p
  END DO ! q

  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(10) ! symmetry BC (perfect MAGNETIC conductor, PMC)
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      resul=U_Face(:,p,q)
      n_loc=normal(:,p,q)
      U_Face_loc(1:3,p,q) =  resul(1:3) - 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
      U_Face_loc(4:6,p,q) = -resul(4:6) + 2*(DOT_PRODUCT(resul(4:6),n_loc))*n_loc
      U_Face_loc(  7,p,q) = -resul(  7)
      U_Face_loc(  8,p,q) =  resul(  8)
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in maxwell/getboundaryflux.f90!')
END SELECT ! BCType
END SUBROUTINE GetBoundaryFlux


END MODULE MOD_GetBoundaryFlux
