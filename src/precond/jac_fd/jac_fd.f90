#include "boltzplatz.h"

MODULE MOD_Jac_FD
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==100) || (PP_TimeDiscMethod==101) || (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==103) || (PP_TimeDiscMethod==104) || (PP_TimeDiscMethod==121)  || (PP_TimeDiscMethod==122)

INTERFACE InitJac_FD
  MODULE PROCEDURE InitJac_FD
END INTERFACE

INTERFACE Jac_FD_slow
  MODULE PROCEDURE Jac_FD_slow
END INTERFACE

!INTERFACE Jac_FD
!  MODULE PROCEDURE Jac_FD
!END INTERFACE

!INTERFACE R_xkEps_Zero
!  MODULE PROCEDURE R_xkEps_Zero
!END INTERFACE

!INTERFACE CalcSourceFD
!  MODULE PROCEDURE CalcSourceFD
!END INTERFACE

!INTERFACE FinalizeJac_FD
!  MODULE PROCEDURE FinalizeJac_FD
!END INTERFACE


PUBLIC::InitJac_FD
!PUBLIC::Jac_FD
PUBLIC::Jac_FD_slow
!PUBLIC::R_xkEps_Zero
!PUBLIC::FinalizeJac_FD
!===================================================================================================================================

CONTAINS


SUBROUTINE InitJac_FD
!===================================================================================================================================
! Allocate global variable R_xkEps
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Jac_FD_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF(PrecondFDInitisDone)THEN
  SWRITE(*,*) "Init Jacobian FD already called."
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT JACOBIAN FD...'

! the local R_xkEps
ALLOCATE(R_xkEps(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N))
ALLOCATE(XK     (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N))
ALLOCATE(RXK    (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N))

! Nullify R_xkEps
R_xkEps=0.

! get machine epsilon
rEps0=SQRT(EPSILON(0.0))
srEps0=1./rEps0

PrecondFDInitisDone=.TRUE.
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT JACOBIAN FD DONE!'
END SUBROUTINE InitJac_FD


SUBROUTINE Jac_FD_slow(t,tStage,tDeriv,iElem,dRdU)
!===================================================================================================================================
! Coputes the Finite Difference-Derivative dRdU for the Calculation of Preconditioner P 
! Attention: dRdU = 0 (in precond.f90)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Jac_FD_Vars,   ONLY: R_xkEps,Xk,sreps0,reps0,Rxk
USE MOD_LinearSolver_Vars, ONLY: nDOFelem
USE MOD_DG,            ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,       ONLY: U,Ut
#ifdef MPI
USE MOD_MPI_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: t, tStage
INTEGER,INTENT(IN)                   :: iElem,tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: dRdU(1:nDOFelem,1:nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                              :: iVar,r,s,ii,jj,kk,iiVar,i,j,k
#ifdef MPI
INTEGER                              :: iProc
#endif
!===================================================================================================================================

Xk=U(:,:,:,:,iElem)
CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.FALSE.)
Rxk=Ut(:,:,:,:,iElem) !linearization Ut of Xk for FD

! nullify
dRdU = 0.

#ifdef MPI
DO iProc=0,nProcessors-1
  s=1
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          IF (iProc.EQ.myRank) THEN
            IF(iElem.GT.PP_nElems) CYCLE
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k) + reps0
          END IF
          CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.FALSE.)
          IF (iProc.EQ.myRank) THEN
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k) 
          END IF
          IF (iProc.EQ.myRank) THEN
            r=1
            DO kk=0,PP_N
              DO jj=0,PP_N
                DO ii=0,PP_N
                  DO iiVar=1,PP_nVar
                    dRdU(r,s) = dRdU(r,s)+(Ut(iiVar,ii,jj,kk,iElem)-Rxk(iiVar,ii,jj,kk))*sreps0
                    r=r+1
                  END DO !iiVar
                END DO !ii
              END DO !jj
            END DO !kk
          END IF
          s=s+1
        END DO !PP_nVar
      END DO !i
    END DO !j
  END DO !k
END DO !iProc
#else
  s=1
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k) + reps0
            CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.FALSE.)
            U(iVar,i,j,k,iElem) = Xk(iVar,i,j,k) 
            r=1
            DO kk=0,PP_N
              DO jj=0,PP_N
                DO ii=0,PP_N
                  DO iiVar=1,PP_nVar
                    dRdU(r,s) = dRdU(r,s)+(Ut(iiVar,ii,jj,kk,iElem)-Rxk(iiVar,ii,jj,kk))*sreps0
                    r=r+1
                  END DO !iiVar
                END DO !ii
              END DO !jj
            END DO !kk
          s=s+1
        END DO !PP_nVar
      END DO !i
    END DO !j
  END DO !k
#endif

!cleanup
CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.TRUE.)

END SUBROUTINE Jac_FD_slow


SUBROUTINE FinalizeJac_FD()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Jac_FD_vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SDEALLOCATE(R_xkEps)
PrecondFdInitisDone = .FALSE.
END SUBROUTINE FinalizeJac_FD


#endif /*(PP_TimeDiscMethod==100) || (PP_TimeDiscMethod==101) || (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==103) || (PP_TimeDiscMethod==104) || (PP_TimeDiscMethod==121)  || (PP_TimeDiscMethod==122)*/
END MODULE MOD_Jac_FD
