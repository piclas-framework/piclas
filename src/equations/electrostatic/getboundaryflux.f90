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
INTERFACE BR1_GetBoundaryFlux
  MODULE PROCEDURE BR1_GetBoundaryFlux
END INTERFACE

PUBLIC::GetBoundaryFlux,BR1_GetBoundaryFlux
!===================================================================================================================================

CONTAINS



SUBROUTINE GetBoundaryFlux(F_Face,BCType,BCState,xGP_Face,normal,tangent1,tangent2,t,tDeriv,U_Face)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Riemann,ONLY:Riemann
USE MOD_Equation,ONLY:ExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: BCType,BCState
REAL,INTENT(INOUT)                   :: U_Face(PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: normal(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: tangent1(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: tangent2(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: t
INTEGER,INTENT(IN)                   :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: F_Face(PP_nVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q
REAL                                 :: n_loc(3),resul(PP_nVar)
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
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))

CASE(3) ! 1st order absorbing BC
  U_Face_loc=0.
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))
  
CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      resul=U_Face(:,p,q)
      n_loc=normal(:,p,q)
      U_Face_loc(1:3,p,q) = resul(1:3) !-resul(1:3) + 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
      U_Face_loc(  4,p,q) = -resul(  4)
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in maxwell/getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE GetBoundaryFlux



SUBROUTINE BR1_GetBoundaryFlux(F_Face,BCType,BCState,xGP_Face,normal,tangent1,tangent2,t,tDeriv,U_Face,dir)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Equation,ONLY:ExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: BCType,BCState
REAL,INTENT(IN)                      :: U_Face(PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: normal(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: tangent1(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: tangent2(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: t
INTEGER,INTENT(IN)                   :: tDeriv,dir
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: F_Face(PP_nVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q,iVar
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N),ar,br,N_loc(1:3,1:3),Prim(8)
REAL                                 :: U_loc(PP_nVar)
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
  DO iVar=1,PP_nVar
  ! BR1 uses arithmetic mean value of states for the Riemann flux
    F_Face(iVar,:,:)=0.5*(U_Face(iVar,:,:)+U_Face_loc(iVar,:,:))*Normal(dir,:,:)
  END DO ! iVar=1,PP_nVar
CASE(4) ! ideal conductor
!      resul=U_Face(:,p,q)
!      n_loc=normal(:,p,q)
!      U_Face_loc(1:3,p,q) = -resul(1:3) + 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
!      U_Face_loc(4:6,p,q) =  resul(4:6) - 2*(DOT_PRODUCT(resul(4:6),n_loc))*n_loc
!      U_Face_loc(  7,p,q) =  resul(  7)
!      U_Face_loc(  8,p,q) = -resul(  8)
CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       ' unknown BC Type in maxwell/getboundaryflux.f90!',BCType,999.)
END SELECT ! BCType
END SUBROUTINE BR1_GetBoundaryFlux

END MODULE MOD_GetBoundaryFlux
