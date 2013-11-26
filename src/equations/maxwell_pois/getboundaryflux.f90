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
INTERFACE FillFlux_BC_Pois
  MODULE PROCEDURE FillFlux_BC
END INTERFACE

PUBLIC::GetBoundaryFlux,FillFlux_BC_Pois
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
USE MOD_Globals,        ONLY:Abort
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
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))

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

  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))
  
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
  CALL Riemann(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),         &
               normal(:,:,:),tangent1(:,:,:),tangent2(:,:,:))

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in maxwell/getboundaryflux.f90!',999,999.)
END SELECT ! BCType
END SUBROUTINE GetBoundaryFlux




SUBROUTINE FillFlux_BC(t,tDeriv,Flux)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,   ONLY: Phi_Minus
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,SurfElem,BCFace_xGP,BC,BoundaryType
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides
USE MOD_Equation_Vars,   ONLY: IniExactFunc
USE MOD_GetBoundaryFlux_Pois,ONLY: GetBoundaryFlux_Pois
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
INTEGER,INTENT(IN) :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux(1:4,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            ::SideID,BCType,BCState,p,q
!===================================================================================================================================
! fill flux for boundary sides
DO SideID=1,nBCSides
  BCType=Boundarytype(BC(SideID),BC_TYPE)
  BCState=Boundarytype(BC(SideID),BC_STATE)
  IF (BCState.EQ.0)BCState=IniExactFunc
  CALL GetBoundaryFlux_Pois(Flux(:,:,:,SideID),BCType,BCState,BCFace_xGP(:,:,:,SideID),NormVec(:,:,:,SideID), &
                TangVec1(:,:,:,SideID), TangVec2(:,:,:,SideID),t,tDeriv,Phi_Minus(:,:,:,SideID) ,SideID)
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,SideID)
    END DO
  END DO
END DO! SideID
END SUBROUTINE FillFlux_BC

END MODULE MOD_GetBoundaryFlux
