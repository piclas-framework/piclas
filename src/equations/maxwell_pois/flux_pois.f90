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
MODULE MOD_Flux_Pois
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
INTERFACE EvalFlux3D_Pois
  MODULE PROCEDURE EvalFlux3D
END INTERFACE

PUBLIC::EvalFlux3D_Pois
!===================================================================================================================================

CONTAINS

SUBROUTINE EvalFlux3D(iElem,f,g,h)
!===================================================================================================================================
! Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: c_corr,c_corr_c2,Phi
USE MOD_Globals_Vars  ,ONLY: c2

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(4,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Uin(4)
INTEGER             :: i,j,k
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      Uin=Phi(:,i,j,k,iElem)
      ! hier der physikalische Fluss ohne die Divergenzkorrektur!
      !A
      f(1,i,j,k) = Uin(2)*c_corr_c2 !Uin(2)*c_corr_c2          ! phi*chi*c^2
      f(2,i,j,k) = Uin(1)*c_corr             ! E1*c_corr
      f(3,i,j,k) = 0
      f(4,i,j,k) = 0
      !B
      g(1,i,j,k) = Uin(3)*c_corr_c2 !Uin(3)*c_corr_c2          ! phi*chi*c^2
      g(2,i,j,k) = 0
      g(3,i,j,k) = Uin(1)*c_corr             ! E2*c_corr
      g(4,i,j,k) = 0
      !C
      h(1,i,j,k) = Uin(4)*c_corr_c2  !Uin(4)*c_corr_c2         ! E3*c_corr
      h(2,i,j,k) = 0
      h(3,i,j,k) = 0
      h(4,i,j,k) = Uin(1)*c_corr             ! E3*c_corr


!      ! hier der physikalische Fluss ohne die Divergenzkorrektur!
!      !A
!      f(1,i,j,k) = -Uin(2)*c_corr_c2 !Uin(2)*c_corr_c2          ! phi*chi*c^2
!      f(2,i,j,k) = -Uin(1)*c_corr             ! E1*c_corr
!      f(3,i,j,k) = 0
!      f(4,i,j,k) = 0
!      !B
!      g(1,i,j,k) = -Uin(3)*c_corr_c2 !Uin(3)*c_corr_c2          ! phi*chi*c^2
!      g(2,i,j,k) = 0
!      g(3,i,j,k) = -Uin(1)*c_corr             ! E2*c_corr
!      g(4,i,j,k) = 0
!      !C
!      h(1,i,j,k) = -Uin(4)*c_corr_c2  !Uin(4)*c_corr_c2         ! E3*c_corr
!      h(2,i,j,k) = 0
!      h(3,i,j,k) = 0
!      h(4,i,j,k) = -Uin(1)*c_corr             ! E3*c_corr
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3D

END MODULE MOD_Flux_Pois
