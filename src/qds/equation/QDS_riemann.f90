!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#if USE_QDS_DG
#include "piclas.h"

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
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
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
INTEGER                                          :: p,q, iVar,L
REAL                                             :: velocompL, velocompR,LambdaMax
!===================================================================================================================================
!Lax-Friedrich
DO iVar=0,7
  L=iVar*5
  DO q=0,PP_N; DO p=0,PP_N
    IF(U_L(1+L,p,q).GT.0.0)THEN
      velocompL = U_L(2+L,p,q)/U_L(1+L,p,q)*nv(1,p,q) + &
                  U_L(3+L,p,q)/U_L(1+L,p,q)*nv(2,p,q) + &
                  U_L(4+L,p,q)/U_L(1+L,p,q)*nv(3,p,q)
    ELSE
      velocompL = 0.0
    END IF
    IF(U_R(1+L,p,q).GT.0.0)THEN
      velocompR = U_R(2+L,p,q)/U_R(1+L,p,q)*nv(1,p,q) + &
                  U_R(3+L,p,q)/U_R(1+L,p,q)*nv(2,p,q) + &
                  U_R(4+L,p,q)/U_R(1+L,p,q)*nv(3,p,q)
    ELSE
      velocompR = 0.0
    END IF
    !IF (ABS(velocompL).GT.ABS(velocompR)) THEN
      !LambdaMax = ABS(velocompL)
    !ELSE
      !LambdaMax = ABS(velocompR)
    !END IF
    !LambdaMax = MERGE(ABS(velocompL),ABS(velocompR),ABS(velocompL).GT.ABS(velocompR))
    LambdaMax = MAX(ABS(velocompL),ABS(velocompR))
    !LambdaMax=QDSMaxVelo


!    Lambda_L = 0.5 * (LambdaMax + ABS(LambdaMax))
!    Lambda_R = 0.5 * (LambdaMax - ABS(LambdaMax))
!    F(1 + iVar*5,p,q) =  (Lambda_L * U_L(1 + iVar*5,p,q) + Lambda_R * U_R(1 + iVar*5,p,q))
!    F(2 + iVar*5,p,q) =  (Lambda_L * U_L(2 + iVar*5,p,q) + Lambda_R * U_R(2 + iVar*5,p,q))
!    F(3 + iVar*5,p,q) =  (Lambda_L * U_L(3 + iVar*5,p,q) + Lambda_R * U_R(3 + iVar*5,p,q))
!    F(4 + iVar*5,p,q) =  (Lambda_L * U_L(4 + iVar*5,p,q) + Lambda_R * U_R(4 + iVar*5,p,q))
!    F(5 + iVar*5,p,q) =  (Lambda_L * U_L(5 + iVar*5,p,q) + Lambda_R * U_R(5 + iVar*5,p,q))



     F(1+L,p,q) =   0.5*(velocompL* U_L(1+L,p,q) + velocompR* U_R(1+L,p,q)) &
                  - 0.5*LambdaMax *(U_R(1+L,p,q) -            U_L(1+L,p,q))

     F(2+L,p,q) =   0.5*(velocompL* U_L(2+L,p,q) + velocompR* U_R(2+L,p,q)) &
                  - 0.5*LambdaMax *(U_R(2+L,p,q) -            U_L(2+L,p,q))

     F(3+L,p,q) =   0.5*(velocompL* U_L(3+L,p,q) + velocompR* U_R(3+L,p,q)) &
                  - 0.5*LambdaMax *(U_R(3+L,p,q) -            U_L(3+L,p,q))

     F(4+L,p,q) =   0.5*(velocompL* U_L(4+L,p,q) + velocompR* U_R(4+L,p,q)) &
                  - 0.5*LambdaMax *(U_R(4+L,p,q) -            U_L(4+L,p,q))

     F(5+L,p,q) =   0.5*(velocompL* U_L(5+L,p,q) + velocompR* U_R(5+L,p,q)) &
                  - 0.5*LambdaMax *(U_R(5+L,p,q) -            U_L(5+L,p,q))

!     F(1 + iVar*5,p,q) = 0.5*(velocompL* U_L(1 + iVar*5,p,q) + velocompR* U_R(1 + iVar*5,p,q))
!     F(2 + iVar*5,p,q) = 0.5*(velocompL* U_L(2 + iVar*5,p,q) + velocompR* U_R(2 + iVar*5,p,q))
!     F(3 + iVar*5,p,q) = 0.5*(velocompL* U_L(3 + iVar*5,p,q) + velocompR* U_R(3 + iVar*5,p,q))
!     F(4 + iVar*5,p,q) = 0.5*(velocompL* U_L(4 + iVar*5,p,q) + velocompR* U_R(4 + iVar*5,p,q))
!     F(5 + iVar*5,p,q) = 0.5*(velocompL* U_L(5 + iVar*5,p,q) + velocompR* U_R(5 + iVar*5,p,q))


!    LambdaMax = MAX(ABS(velocompL), ABS(velocompR))
!    F(1 + iVar*5,p,q) =  0.5*(velocompL * U_L(1 + iVar*5,p,q) +velocompR * U_R(1 + iVar*5,p,q) &
!          + LambdaMax *(U_L(1 + iVar*5,p,q) -  U_R(1 + iVar*5,p,q)))
!    F(2 + iVar*5,p,q) =  0.5*(velocompL * U_L(2 + iVar*5,p,q) +velocompR * U_R(2 + iVar*5,p,q) &
!          + LambdaMax *(U_L(2 + iVar*5,p,q) -  U_R(2 + iVar*5,p,q)))
!    F(3 + iVar*5,p,q) =  0.5*(velocompL * U_L(3 + iVar*5,p,q) +velocompR * U_R(3 + iVar*5,p,q) &
!          + LambdaMax *(U_L(3 + iVar*5,p,q) -  U_R(3 + iVar*5,p,q)))
!    F(4 + iVar*5,p,q) =  0.5*(velocompL * U_L(4 + iVar*5,p,q) +velocompR * U_R(4 + iVar*5,p,q) &
!          + LambdaMax *(U_L(4 + iVar*5,p,q) -  U_R(4 + iVar*5,p,q)))
!    F(5 + iVar*5,p,q) =  0.5*(velocompL * U_L(5 + iVar*5,p,q) +velocompR * U_R(5 + iVar*5,p,q) &
!          + LambdaMax *(U_L(5 + iVar*5,p,q) -  U_R(5 + iVar*5,p,q)))
  END DO; END DO
END DO


END SUBROUTINE RiemannQDS


END MODULE MOD_QDS_Riemann
#endif /*USE_QDS_DG*/
