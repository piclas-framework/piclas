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

INTERFACE EvalFlux3DDielectric
  MODULE PROCEDURE EvalFlux3DDielectric
END INTERFACE

PUBLIC::EvalFlux3D
PUBLIC::EvalFlux3DDielectric
!===================================================================================================================================

CONTAINS

SUBROUTINE EvalFlux3D(iElem,f,g,h)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: c_corr,c_corr_c2
USE MOD_Globals_Vars  ,ONLY: c2
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
        f(1,i,j,k) = Uin(8)*c_corr_c2          ! phi*chi*c^2      P1
        f(2,i,j,k) = Uin(6)*c2                 ! B3*c^2           Q1
        f(3,i,j,k) =-Uin(5)*c2                 ! -B2*c^2          R1
        f(4,i,j,k) = Uin(7)*c_corr             ! psi*c_corr       L1
        f(5,i,j,k) =-Uin(3)                    ! -E3              M1
        f(6,i,j,k) = Uin(2)                    ! E2               N1
        f(7,i,j,k) = Uin(4)*c_corr_c2          ! B1*c_corr*c^2    S1
        f(8,i,j,k) = Uin(1)*c_corr             ! E1*c_corr        T1
        !B
        g(1,i,j,k) =-f(2,i,j,k)                ! -B3*c^2          P2
        g(2,i,j,k) = f(1,i,j,k)                ! phi*c_corr*c^2   Q2
        g(3,i,j,k) = Uin(4)*c2                 ! B1*c^2           R2
        g(4,i,j,k) = Uin(3)                    ! E3               L2
        g(5,i,j,k) = f(4,i,j,k)                ! psi*c_corr       M2
        g(6,i,j,k) =-Uin(1)                    ! -E1              N2
        g(7,i,j,k) = Uin(5)*c_corr_c2          ! B2*c_corr*c^2    S2
        g(8,i,j,k) = Uin(2)*c_corr             ! E2*c_corr        T2
        !C
        h(1,i,j,k) =-f(3,i,j,k)                ! B2*c^2           P3
        h(2,i,j,k) =-g(3,i,j,k)                ! -B1*c^2          Q3
        h(3,i,j,k) = f(1,i,j,k)                ! phi*c_corr*c^2   R3
        h(4,i,j,k) =-Uin(2)                    ! -E2              L3
        h(5,i,j,k) = Uin(1)                    ! E1               M3
        h(6,i,j,k) = f(4,i,j,k)                ! psi*c_corr       N3
        h(7,i,j,k) = Uin(6)*c_corr_c2          ! B3*c_corr*c^2    S3
        h(8,i,j,k) = Uin(3)*c_corr             ! E3*c_corr        T3
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3D


SUBROUTINE EvalFlux3DDielectric(iElem,f,g,h)
!===================================================================================================================================
! calculate the 3D physical fluxes for each DOF i,j,k with additional dielectric material factors
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars   ,ONLY: c_corr,c_corr_c2
USE MOD_Globals_Vars    ,ONLY: c2
USE MOD_DG_Vars         ,ONLY: U
USE MOD_Dielectric_Vars ,ONLY: ElemToDielectric,DielectricConstant_inv
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
REAL                :: Uin(8),die
INTEGER             :: i,j,k
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      Uin=U(:,i,j,k,iElem)
      die=DielectricConstant_Inv(i,j,k,ElemToDielectric(iElem))
        ! hier der physikalische Fluss ohne die Divergenzkorrektur!
        !A
        f(1,i,j,k) = Uin(8)*c_corr_c2*die  ! phi*chi*c^2      P1
        f(2,i,j,k) = Uin(6)*c2 * die       ! B3*c^2           Q1 * 1/(EpsR*MuR)
        f(3,i,j,k) =-Uin(5)*c2 * die       ! -B2*c^2          R1 * 1/(EpsR*MuR)
        f(4,i,j,k) = Uin(7)*c_corr         ! psi*c_corr       L1
        f(5,i,j,k) =-Uin(3)                ! -E3              M1
        f(6,i,j,k) = Uin(2)                ! E2               N1
        f(7,i,j,k) = Uin(4)*c_corr_c2*die  ! B1*c_corr*c^2    S1
        f(8,i,j,k) = Uin(1)*c_corr         ! E1*c_corr        T1
        !B
        g(1,i,j,k) =-f(2,i,j,k)            ! -B3*c^2          P2
        g(2,i,j,k) = f(1,i,j,k)            ! phi*c_corr*c^2   Q2
        g(3,i,j,k) = Uin(4)*c2 * die       ! B1*c^2           R2 * 1/(EpsR*MuR)
        g(4,i,j,k) = Uin(3)                ! E3               L2
        g(5,i,j,k) = f(4,i,j,k)            ! psi*c_corr       M2
        g(6,i,j,k) =-Uin(1)                ! -E1              N2
        g(7,i,j,k) = Uin(5)*c_corr_c2*die  ! B2*c_corr*c^2    S2
        g(8,i,j,k) = Uin(2)*c_corr         ! E2*c_corr        T2
        !C
        h(1,i,j,k) =-f(3,i,j,k)            ! B2*c^2           P3
        h(2,i,j,k) =-g(3,i,j,k)            ! -B1*c^2          Q3
        h(3,i,j,k) = f(1,i,j,k)            ! phi*c_corr*c^2   R3
        h(4,i,j,k) =-Uin(2)                ! -E2              L3
        h(5,i,j,k) = Uin(1)                ! E1               M3
        h(6,i,j,k) = f(4,i,j,k)            ! psi*c_corr       N3
        h(7,i,j,k) = Uin(6)*c_corr_c2*die  ! B3*c_corr*c^2    S3
        h(8,i,j,k) = Uin(3)*c_corr         ! E3*c_corr        T3
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3DDielectric

END MODULE MOD_Flux
