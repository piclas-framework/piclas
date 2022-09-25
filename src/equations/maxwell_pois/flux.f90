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
MODULE MOD_Flux
!===================================================================================================================================
! Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
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
INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D
END INTERFACE

PUBLIC::EvalFlux3D
!===================================================================================================================================

CONTAINS

SUBROUTINE EvalFlux3D(iElem,f,g,h)
!===================================================================================================================================
! Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c2
USE MOD_Equation_Vars ,ONLY: c_corr,c_corr_c2
USE MOD_DG_Vars       ,ONLY: U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Uin(8)
INTEGER             :: i,j,k
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      Uin=U(:,i,j,k,iElem)
      ! hier der physikalische Fluss ohne die Divergenzkorrektur!
      !A
      f(1,i,j,k) = Uin(8)*c_corr_c2          ! phi*chi*c^2
      f(2,i,j,k) = Uin(6)*c2                 ! B3*c^2
      f(3,i,j,k) =-Uin(5)*c2                 ! -B2*c^2
      f(4,i,j,k) = Uin(7)*c_corr             ! psi*c_corr
      f(5,i,j,k) =-Uin(3)                    ! -E3
      f(6,i,j,k) = Uin(2)                    ! E2
      f(7,i,j,k) = Uin(4)*c_corr_c2          ! B1*c_corr*c^2
      f(8,i,j,k) = Uin(1)*c_corr             ! E1*c_corr
      !B
      g(1,i,j,k) =-f(2,i,j,k)                ! -B3*c^2
      g(2,i,j,k) = f(1,i,j,k)                ! phi*c_corr*c^2
      g(3,i,j,k) = Uin(4)*c2                 ! B1*c^2
      g(4,i,j,k) = Uin(3)                    ! E3
      g(5,i,j,k) = f(4,i,j,k)                ! psi*c_corr
      g(6,i,j,k) =-Uin(1)                    ! -E1
      g(7,i,j,k) = Uin(5)*c_corr_c2          ! B2*c_corr*c^2
      g(8,i,j,k) = Uin(2)*c_corr             ! E2*c_corr
      !C
      h(1,i,j,k) =-f(3,i,j,k)                ! B2*c^2
      h(2,i,j,k) =-g(3,i,j,k)                ! -B1*c^2
      h(3,i,j,k) = f(1,i,j,k)                ! phi*c_corr*c^2
      h(4,i,j,k) =-Uin(2)                    ! -E2
      h(5,i,j,k) = Uin(1)                    ! E1
      h(6,i,j,k) = f(4,i,j,k)                ! psi*c_corr
      h(7,i,j,k) = Uin(6)*c_corr_c2          ! B3*c_corr*c^2
      h(8,i,j,k) = Uin(3)*c_corr             ! E3*c_corr
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3D

END MODULE MOD_Flux
