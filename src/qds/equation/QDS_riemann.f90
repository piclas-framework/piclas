#include "boltzplatz.h"

MODULE MOD_QDS_Riemann
!===================================================================================================================================
!> Contains the routines to
!> - determine the riemann flux for QDS DG method
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE RiemannQDS
  MODULE PROCEDURE RiemannQDS
END INTERFACE

PUBLIC::RiemannQDS
!===================================================================================================================================
CONTAINS
SUBROUTINE RiemannQDS(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_QDS_DG_Vars,     ONLY:QDSnVar,QDSMaxVelo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(QDSnVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(QDSnVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL                                             :: n_loc(3),A_p(4,4),A_n(4,4)
INTEGER                                          :: Count_1,Count_2, iVar
!REAL                                             :: D(3,3)                  ! auxiliary matrices used 
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!REAL                                            :: Lambda_L, Lambda_R
REAL                                            :: velocompL, velocompR,LambdaMax
!===================================================================================================================================
!Lax-Friedrich
DO iVar=0,7
  DO Count_2=0,PP_N; DO Count_1=0,PP_N 
    velocompL = U_L(2+iVar*5,Count_1,Count_2)*nv(1,Count_1,Count_2) + U_L(3+iVar*5,Count_1,Count_2)*nv(2,Count_1,Count_2) &
            + U_L(4+iVar*5,Count_1,Count_2)*nv(3,Count_1,Count_2)
    velocompR = U_R(2+iVar*5,Count_1,Count_2)*nv(1,Count_1,Count_2) + U_R(3+iVar*5,Count_1,Count_2)*nv(2,Count_1,Count_2) &
            + U_R(4+iVar*5,Count_1,Count_2)*nv(3,Count_1,Count_2)
    !IF (ABS(velocompL).GT.ABS(velocompR)) THEN
      !LambdaMax = ABS(velocompL)
    !ELSE
      !LambdaMax = ABS(velocompR)
    !END IF
!    LambdaMax = MERGE(velocompL, velocompR, ABS(velocompL).GT.ABS(velocompR))
    LambdaMax=QDSMaxVelo

!    Lambda_L = 0.5 * (LambdaMax + ABS(LambdaMax))
!    Lambda_R = 0.5 * (LambdaMax - ABS(LambdaMax))
!    F(1 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(1 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(1 + iVar*5,Count_1,Count_2)) 
!    F(2 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(2 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(2 + iVar*5,Count_1,Count_2))
!    F(3 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(3 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(3 + iVar*5,Count_1,Count_2))
!    F(4 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(4 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(4 + iVar*5,Count_1,Count_2))
!    F(5 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(5 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(5 + iVar*5,Count_1,Count_2)) 
    


     F(1 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(1 + iVar*5,Count_1,Count_2) + velocompR* U_R(1 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(1 + iVar*5,Count_1,Count_2)- U_L(1 + iVar*5,Count_1,Count_2))
     F(2 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(2 + iVar*5,Count_1,Count_2) + velocompR* U_R(2 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(2 + iVar*5,Count_1,Count_2)- U_L(2 + iVar*5,Count_1,Count_2))
     F(3 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(3 + iVar*5,Count_1,Count_2) + velocompR* U_R(3 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(3 + iVar*5,Count_1,Count_2)- U_L(3 + iVar*5,Count_1,Count_2))
     F(4 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(4 + iVar*5,Count_1,Count_2) + velocompR* U_R(4 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(4 + iVar*5,Count_1,Count_2)- U_L(4 + iVar*5,Count_1,Count_2))
     F(5 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(5 + iVar*5,Count_1,Count_2) + velocompR* U_R(5 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(5 + iVar*5,Count_1,Count_2)- U_L(5 + iVar*5,Count_1,Count_2))

!     F(1 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(1 + iVar*5,Count_1,Count_2) + velocompR* U_R(1 + iVar*5,Count_1,Count_2))
!     F(2 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(2 + iVar*5,Count_1,Count_2) + velocompR* U_R(2 + iVar*5,Count_1,Count_2))
!     F(3 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(3 + iVar*5,Count_1,Count_2) + velocompR* U_R(3 + iVar*5,Count_1,Count_2))
!     F(4 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(4 + iVar*5,Count_1,Count_2) + velocompR* U_R(4 + iVar*5,Count_1,Count_2))
!     F(5 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(5 + iVar*5,Count_1,Count_2) + velocompR* U_R(5 + iVar*5,Count_1,Count_2))

     
!    LambdaMax = MAX(ABS(velocompL), ABS(velocompR))
!    F(1 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(1 + iVar*5,Count_1,Count_2) +velocompR * U_R(1 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(1 + iVar*5,Count_1,Count_2) -  U_R(1 + iVar*5,Count_1,Count_2)))
!    F(2 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(2 + iVar*5,Count_1,Count_2) +velocompR * U_R(2 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(2 + iVar*5,Count_1,Count_2) -  U_R(2 + iVar*5,Count_1,Count_2)))
!    F(3 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(3 + iVar*5,Count_1,Count_2) +velocompR * U_R(3 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(3 + iVar*5,Count_1,Count_2) -  U_R(3 + iVar*5,Count_1,Count_2)))
!    F(4 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(4 + iVar*5,Count_1,Count_2) +velocompR * U_R(4 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(4 + iVar*5,Count_1,Count_2) -  U_R(4 + iVar*5,Count_1,Count_2)))
!    F(5 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(5 + iVar*5,Count_1,Count_2) +velocompR * U_R(5 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(5 + iVar*5,Count_1,Count_2) -  U_R(5 + iVar*5,Count_1,Count_2)))
  END DO; END DO
END DO


END SUBROUTINE RiemannQDS


END MODULE MOD_QDS_Riemann
