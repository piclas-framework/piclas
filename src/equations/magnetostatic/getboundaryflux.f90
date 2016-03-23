#include "boltzplatz.h"

MODULE MOD_GetBoundaryFlux
!===================================================================================================================================
! Fills the boundary part of the flux list (from nInnerSidess+1 to nSides)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE


! Public Part ----------------------------------------------------------------------------------------------------------------------

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
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Riemann,ONLY:Riemann
USE MOD_Equation,ONLY:ExactFunc
USE MOD_Equation_Vars,ONLY:IniExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: BCType,BCState
REAL,INTENT(IN)                   :: U_Face(PP_nVar,0:PP_N,0:PP_N)
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
REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
SELECT CASE(BCType)
CASE(1) !Periodic already filled!
CASE(2) ! exact BC = Dirichlet BC !!
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      CALL ExactFunc(IniExactFunc,t,tDeriv,xGP_Face(:,p,q),U_Face_loc(:,p,q))
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:), normal(:,:,:))
CASE(22) ! exact BC = Dirichlet BC !!
! SPECIAL BC: BCState specifies exactfunc to be used!!
DO q=0,PP_N
DO p=0,PP_N
  CALL ExactFunc(BCState,t,tDeriv,xGP_Face(:,p,q),U_Face_loc(:,p,q))
END DO ! p
END DO ! q
! Dirichlet means that we use the gradients from inside the grid cell
CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))
CASE(3) !Neumann Flux=0.
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE DEFAULT ! unknown BCType
  CALL abort(&
__STAMP__&
,'no BC defined in navierstokes/getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE GetBoundaryFlux

END MODULE MOD_GetBoundaryFlux
