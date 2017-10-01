#include "boltzplatz.h"

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
USE MOD_QDS_DG_Vars, ONLY:QDSnVar,UQDS
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
INTEGER             :: i,j,k,iVar
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      !Uin=UQDS(:,i,j,k,iElem)
      DO iVar=0,7
!        ! hier der physikalische Fluss ohne die Divergenzkorrektur!
!        !A
!        f(1+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
!        f(2+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
!        f(3+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
!        f(4+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
!        f(5+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
!        !B
!        g(1+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
!        g(2+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
!        g(3+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
!        g(4+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
!        g(5+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
!        !C
!        h(1+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
!        h(2+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
!        h(3+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
!        h(4+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
!        h(5+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
        f(1+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)
        IF (UQDS(1+5*iVar,i,j,k,iElem).GT.0.0) THEN
          f(2+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem)
          f(3+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem)
          f(4+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem)
        ELSE
          f(2+5*iVar,i,j,k) = 0.0
          f(3+5*iVar,i,j,k) = 0.0
          f(4+5*iVar,i,j,k) = 0.0
        END IF
        f(5+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem)
        !B
        g(1+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)
        IF (UQDS(1+5*iVar,i,j,k,iElem).GT.0.0) THEN
          g(2+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem)
          g(3+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem)
          g(4+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem)
        ELSE
          g(2+5*iVar,i,j,k) = 0.0
          g(3+5*iVar,i,j,k) = 0.0
          g(4+5*iVar,i,j,k) = 0.0
        END IF
        g(5+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem)
        !C
        h(1+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)
        IF (UQDS(1+5*iVar,i,j,k,iElem).GT.0.0) THEN
          h(2+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem)
          h(3+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem)
          h(4+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)/UQDS(1+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem)
        ELSE
          h(2+5*iVar,i,j,k) = 0.0
          h(3+5*iVar,i,j,k) = 0.0
          h(4+5*iVar,i,j,k) = 0.0
        END IF
        h(5+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem)
      END DO
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3DQDS


END MODULE MOD_QDS_Flux
