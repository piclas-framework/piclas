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

MODULE MOD_QDS_Flux
!===================================================================================================================================
!> Contains the routines to
!> - physical flux for the QDS-DG method
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
INTERFACE EvalFlux3DQDS
  MODULE PROCEDURE EvalFlux3DQDS
END INTERFACE

PUBLIC::EvalFlux3DQDS
!===================================================================================================================================
CONTAINS
SUBROUTINE EvalFlux3DQDS(iElem,f,g,h)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_QDS_DG_Vars,        ONLY:UQDS
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(QDSnVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                :: Uin(QDSnVar)
INTEGER             :: i,j,k,iVar,L
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      !Uin=UQDS(:,i,j,k,iElem)
      DO iVar=0,7
!        ! hier der physikalische Fluss ohne die Divergenzkorrektur!
!        !A
!        f(1+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(1+L,i,j,k,iElem)
!        f(2+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
!        f(3+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
!        f(4+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
!        f(5+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
!        !B
!        g(1+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(1+L,i,j,k,iElem)
!        g(2+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
!        g(3+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
!        g(4+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
!        g(5+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
!        !C
!        h(1+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(1+L,i,j,k,iElem)
!        h(2+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
!        h(3+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
!        h(4+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
!        h(5+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
        L=5*iVar
        f(1+L,i,j,k) = UQDS(2+L,i,j,k,iElem)
        IF (UQDS(1+L,i,j,k,iElem).GT.0.0) THEN
          f(2+L,i,j,k) = UQDS(2+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
          f(3+L,i,j,k) = UQDS(2+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
          f(4+L,i,j,k) = UQDS(2+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
        ELSE
          f(2+L,i,j,k) = 0.0
          f(3+L,i,j,k) = 0.0
          f(4+L,i,j,k) = 0.0
        END IF
        f(5+L,i,j,k) = UQDS(2+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
        !B
        g(1+L,i,j,k) = UQDS(3+L,i,j,k,iElem)
        IF (UQDS(1+L,i,j,k,iElem).GT.0.0) THEN
          g(2+L,i,j,k) = UQDS(3+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
          g(3+L,i,j,k) = UQDS(3+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
          g(4+L,i,j,k) = UQDS(3+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
        ELSE
          g(2+L,i,j,k) = 0.0
          g(3+L,i,j,k) = 0.0
          g(4+L,i,j,k) = 0.0
        END IF
        g(5+L,i,j,k) = UQDS(3+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
        !C
        h(1+L,i,j,k) = UQDS(4+L,i,j,k,iElem)
        IF (UQDS(1+L,i,j,k,iElem).GT.0.0) THEN
          h(2+L,i,j,k) = UQDS(4+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(2+L,i,j,k,iElem)
          h(3+L,i,j,k) = UQDS(4+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(3+L,i,j,k,iElem)
          h(4+L,i,j,k) = UQDS(4+L,i,j,k,iElem)/UQDS(1+L,i,j,k,iElem)*UQDS(4+L,i,j,k,iElem)
        ELSE
          h(2+L,i,j,k) = 0.0
          h(3+L,i,j,k) = 0.0
          h(4+L,i,j,k) = 0.0
        END IF
        h(5+L,i,j,k) = UQDS(4+L,i,j,k,iElem)*UQDS(5+L,i,j,k,iElem)
      END DO
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3DQDS


END MODULE MOD_QDS_Flux
#endif /*USE_QDS_DG*/
