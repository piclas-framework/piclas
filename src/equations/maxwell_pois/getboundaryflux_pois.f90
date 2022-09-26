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

MODULE MOD_GetBoundaryFlux_Pois
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
INTERFACE GetBoundaryFlux_Pois
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE


PUBLIC::GetBoundaryFlux_Pois
!===================================================================================================================================

CONTAINS



SUBROUTINE GetBoundaryFlux(F_Face,BCType,BCState,xGP_Face,normal,t,tDeriv,U_Face,iSide)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_PreProc
USE MOD_Riemann_Pois           ,ONLY: Riemann_Pois
USE MOD_Equation               ,ONLY: ExactFunc
USE MOD_Globals_Vars           ,ONLY: c,c_inv
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: BCType,BCState
REAL,INTENT(INOUT)                   :: U_Face(4,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: normal(3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                      :: t
INTEGER,INTENT(IN)                   :: tDeriv
INTEGER,INTENT(IN)                   :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: F_Face(4,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q
REAL                                 :: n_loc(3),resul(4)
REAL                                 :: U_Face_loc(4,0:PP_N,0:PP_N)
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
  CALL Riemann_Pois(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(3) ! 1st order absorbing BC
  U_Face_loc=0.

!  DO q=0,PP_N
!    DO p=0,PP_N
!      resul=U_Face(:,p,q)
!      n_loc=normal(:,p,q)
!!      U_Face_loc(1,p,q) = resul(1) - c*DOT_PRODUCT(resul(2:4),n_loc)
!      U_Face_loc(2:4,p,q) = -resul(2:4) + c_inv*U_Face(1,p,q)*normal(1:3,p,q)

!    END DO ! p
!  END DO

  CALL Riemann_Pois(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      resul=U_Face(:,p,q)
      n_loc=normal(:,p,q)
    ! U_Face_loc(1,p,q) = 2. * 1000. - resul(1)  !- c*DOT_PRODUCT(resul(2:4),n_loc)
      U_Face_loc(1,p,q) = 2. * (PartBound%Voltage(PartBound%MapToPartBC(BC(iSide))) - resul(1)  !+ c*DOT_PRODUCT(resul(2:4),n_loc)
      U_Face_loc(2:4,p,q) = + resul(2:4) !- 1./c*resul(1)*n_loc
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann_Pois(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE(10) ! symmetry BC
  ! Determine the exact BC state
  DO q=0,PP_N
    DO p=0,PP_N
      resul=U_Face(:,p,q)
      n_loc=normal(:,p,q)
      U_Face_loc(1,p,q) = resul(1)
      U_Face_loc(2:4,p,q) = - resul(2:4)
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  CALL Riemann_Pois(F_Face(:,:,:),U_Face(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))

CASE DEFAULT ! unknown BCType
  CALL abort(&
__STAMP__&
,'no BC defined in maxwell/getboundaryflux.f90!')
END SELECT ! BCType
END SUBROUTINE GetBoundaryFlux


END MODULE MOD_GetBoundaryFlux_Pois
