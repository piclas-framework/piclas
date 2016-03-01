#include "boltzplatz.h"

MODULE MOD_CalcTimeStep
!===================================================================================================================================
! Low-Storage Runge-Kutta integration of degree 3 for one step.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!===================================================================================================================================

CONTAINS

FUNCTION CALCTIMESTEP()
!===================================================================================================================================
! Calculate the time step for the current update of U for the Euler-Equations
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars,ONLY:c,c_corr
USE MOD_TimeDisc_Vars,ONLY:CFLScale
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                         :: CalcTimeStep
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL                         :: Max_Lambda1,Max_Lambda2,Max_Lambda3
REAL                         :: TimeStepConv
!===================================================================================================================================
TimeStepConv=HUGE(1.)
DO iElem=1,PP_nElems
  Max_Lambda1=0.
  Max_Lambda2=0.
  Max_Lambda3=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        ! Convective Eigenvalues
! VERSION 1 & 2: -----------------------------
        Max_Lambda1=MAX(Max_Lambda1,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)))))
        Max_Lambda2=MAX(Max_Lambda2,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)))))
        Max_Lambda3=MAX(Max_Lambda3,sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)))))
! --------------------------------------------
! VERSION 3: ---------------------------------
!        Max_Lambda1=MAX(&
!                        Max_Lambda1,&
!                        sJ(i,j,k,iElem)*(MAX(1.,c_corr)*c &
!                        *(SQRT(SUM(Metrics_fTilde(:,i,j,k,iElem)*Metrics_fTilde(:,i,j,k,iElem)) &
!                              +SUM(Metrics_gTilde(:,i,j,k,iElem)*Metrics_gTilde(:,i,j,k,iElem)) &
!                              +SUM(Metrics_hTilde(:,i,j,k,iElem)*Metrics_hTilde(:,i,j,k,iElem)) ) ))&
!                        )
! --------------------------------------------
      END DO ! i
    END DO ! j
  END DO ! k
! VERSION 3: ---------------------------------
!   TimeStepConv=MIN(TimeStepConv,CFLScale*2./(Max_Lambda1))
! --------------------------------------------
! VERSION 2: quadratic superposition
  TimeStepConv=MIN(TimeStepConv,CFLScale*2./SQRT(Max_Lambda1**2+Max_Lambda2**2+Max_Lambda3**2))
! --------------------------------------------
! VERSION 1: linear superposition 
!  TimeStepConv=MIN(TimeStepConv,CFLScale*2./(Max_Lambda1+Max_Lambda2+Max_Lambda3))
! --------------------------------------------
  IF(TimeStepConv.NE.TimeStepConv)THEN
    ERRWRITE(*,'(A,I0,A,I0,A,I0,A,I0)')'Convective timestep NaN on proc ',myRank,' at global position (iElem): ',iElem
    ERRWRITE(*,*)'dt_conv=',TimeStepConv
    CALL abort(&
        __STAMP__&
        ,'Convective timestep NaN!',999,999.)
  END IF
END DO ! iElem
#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStepConv,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
CalcTimeStep=TimeStepConv
END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
