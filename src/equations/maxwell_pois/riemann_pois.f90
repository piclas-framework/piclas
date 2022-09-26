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
MODULE MOD_Riemann_Pois
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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
INTERFACE Riemann_Pois
  MODULE PROCEDURE Riemann
END INTERFACE

PUBLIC::Riemann_Pois
!===================================================================================================================================

CONTAINS

SUBROUTINE Riemann(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: eta_c,c_corr,c_corr_c,c_corr_c2
USE MOD_Globals_Vars  ,ONLY: c,c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(4,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(4,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),A_p(4,4),A_n(4,4)
INTEGER                                          :: Count_1,Count_2
!REAL                                             :: D(3,3)                  ! auxiliary matrices used
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!===================================================================================================================================
! Gauss point i,j
DO Count_2=0,PP_N
  DO Count_1=0,PP_N
    n_loc(:)=nv(:,Count_1,Count_2)

!    A_p(1,1) = c_corr_c
!    A_p(2,2) = c_corr_c * n_loc(1)*n_loc(1)
!    A_p(2,3) = c_corr_c* n_loc(1)*n_loc(2)
!    A_p(2,4) = c_corr_c* n_loc(1)*n_loc(3)
!    A_p(3,2) = A_p(2,3)
!    A_p(3,3) = c_corr_c * n_loc(2)*n_loc(2)
!    A_p(3,4) = c_corr_c* n_loc(2)*n_loc(3)
!    A_p(4,2) = A_p(2,4)
!    A_p(4,3) = A_p(3,4)
!    A_p(4,4) = c_corr_c * n_loc(3)*n_loc(3)

!    !negative A-Matrix
!    A_n(1,1)=-A_p(1,1)
!    A_n(2:4,2:4)=-A_p(2:4,2:4)
!
!   ! !positive A-Matrix-Divergence-Correction-Term
!    A_p(1,2) = -c_corr_c2*n_loc(1)
!    A_p(1,3) = -c_corr_c2*n_loc(2)
!    A_p(1,4) = -c_corr_c2*n_loc(3)
!    A_p(2,1) = -c_corr*n_loc(1)
!    A_p(3,1) = -c_corr*n_loc(2)
!    A_p(4,1) = -c_corr*n_loc(3)
!    !negative A-Matrix-Divergence-Correction-Term
!    A_n(1,2:4) = A_p(1,2:4) !c_corr*c*c*n(1)
!    A_n(2:4,1)= A_p(2:4,1)


!--- for original version see below (easier to understand)

    A_p(1,1) = c_corr_c
    A_p(2,2) = c_corr_c * n_loc(1)*n_loc(1)
    A_p(2,3) = c_corr_c* n_loc(1)*n_loc(2)
    A_p(2,4) = c_corr_c* n_loc(1)*n_loc(3)
    A_p(3,2) = A_p(2,3)
    A_p(3,3) = c_corr_c * n_loc(2)*n_loc(2)
    A_p(3,4) = c_corr_c* n_loc(2)*n_loc(3)
    A_p(4,2) = A_p(2,4)
    A_p(4,3) = A_p(3,4)
    A_p(4,4) = c_corr_c * n_loc(3)*n_loc(3)

    !negative A-Matrix
    A_n(1,1)=-A_p(1,1)
    A_n(2:4,2:4)=-A_p(2:4,2:4)

   ! !positive A-Matrix-Divergence-Correction-Term
    A_p(1,2) = c_corr_c2*n_loc(1)
    A_p(1,3) = c_corr_c2*n_loc(2)
    A_p(1,4) = c_corr_c2*n_loc(3)
    A_p(2,1) = c_corr*n_loc(1)
    A_p(3,1) = c_corr*n_loc(2)
    A_p(4,1) = c_corr*n_loc(3)
    !negative A-Matrix-Divergence-Correction-Term
    A_n(1,2:4) = A_p(1,2:4) !c_corr*c*c*n(1)
    A_n(2:4,1)= A_p(2:4,1)

    ! Warum 0.5 -> Antwort im Taube/Dumbser-Paper. Im Munz/Schneider Paper fehlt das 0.5 lustigerweise.


    F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
  END DO
END DO
END SUBROUTINE Riemann


END MODULE MOD_Riemann_Pois
