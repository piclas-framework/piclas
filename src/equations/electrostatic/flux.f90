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
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars,ONLY:c2,c_corr,c_corr_c2
USE MOD_DG_Vars,ONLY:U
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
      Uin=U(:,i,j,k,iElem)
      ! hier der physikalische Fluss ohne die Divergenzkorrektur!
      !A
      f(1,i,j,k) = Uin(4)*c_corr_c2          ! phi*chi*c^2
      f(2,i,j,k) = 0   
      f(3,i,j,k) = 0 
      f(4,i,j,k) = Uin(1)*c_corr             ! E1*c_corr
      !B
      g(1,i,j,k) = 0
      g(2,i,j,k) = f(1,i,j,k)   
      g(3,i,j,k) = 0 
      g(4,i,j,k) = Uin(2)*c_corr             ! E2*c_corr
      !C                                    
      h(1,i,j,k) = 0
      h(2,i,j,k) = 0
      h(3,i,j,k) = f(1,i,j,k)
      h(4,i,j,k) = Uin(3)*c_corr             ! E3*c_corr
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3D

END MODULE MOD_Flux
