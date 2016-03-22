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
! Calculate the time step for the current update of U for the Linear Scalar Advection Equation du/dt + a du/dx = 0
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_TimeDisc_Vars,ONLY:CFLScale
USE MOD_PreProc
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
REAL                         :: TimeStep(2)
REAL                         :: Lambda_v1,Lambda_v2,Lambda_v3
REAL                         :: MaxLambda_v
!===================================================================================================================================

MaxLambda_v=1.0E-10  ! Viscous
Lambda_v1=1.0E-10
Lambda_v2=1.0E-10
Lambda_v3=1.0E-10
TimeStep=HUGE(1.)
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        Lambda_v1=MAX(Lambda_v1,(SUM((Metrics_fTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        Lambda_v2=MAX(Lambda_v2,(SUM((Metrics_gTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        Lambda_v3=MAX(Lambda_v3,(SUM((Metrics_hTilde(:,i,j,k,iElem)*sJ(i,j,k,iElem))**2)))
        maxLambda_v=MAX(maxLambda_v,(Lambda_v1+Lambda_v2+Lambda_v3))
      END DO ! i
    END DO ! j
  END DO ! k          
END DO ! iElem=1,PP_nElems
TimeStep(2)=MIN(TimeStep(2),4./maxLambda_v)
#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif
CalcTimeStep=MINVAL(TimeStep)
IF(CalcTimeStep.NE.CalcTimeStep)THEN
  SWRITE(*,*)' ******* Exit: Timestep NaN *******'
  CALL abort(&
  __STAMP__&
  ,'Flexi crashed!',999,999.)
END IF
END FUNCTION CalcTimeStep

END MODULE MOD_CalcTimeStep
