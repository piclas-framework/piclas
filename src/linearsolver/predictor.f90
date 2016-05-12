#include "boltzplatz.h"

MODULE MOD_Predictor
!===================================================================================================================================
! Contains routines to predict the solution at the new stage
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTEGER :: PredictorType


INTERFACE Predictor
  MODULE PROCEDURE Predictor
END INTERFACE

PUBLIC:: InitPredictor, FinalizePredictor,Predictor,StorePredictor
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPredictor()
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,          ONLY:GETINT
USE MOD_Linearsolver_vars,    ONLY:Upast
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
PredictorType = GETINT('Predictor','0')

SELECT CASE(PredictorType)
CASE(0,1,2,3)
  ! nothing to do
CASE(4) ! linear prediction
  ALLOCATE(Upast(1:PP_nVar,0:PP_nVar,0:PP_nVar,0:PP_nVar,1:PP_nElems,-1:0))
  upast=0.
CASE(5)
  ALLOCATE(Upast(1:PP_nVar,0:PP_nVar,0:PP_nVar,0:PP_nVar,1:PP_nElems,-2:0))
  upast=0.
CASE DEFAULT
  CALL abort(&
__STAMP__ &
,'PredictorType not implemented!',PredictorType,999.)
END SELECT

END SUBROUTINE InitPredictor


SUBROUTINE Predictor(iStage,dt,Un,FieldStage,FieldSource)
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,          ONLY: U
USE MOD_LinearSolver_Vars,ONLY: LinSolverRHS,Upast
#if (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==105) || (PP_TimeDiscMethod==122)
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bsO3,RK_bs
#endif
#if (PP_TimeDiscMethod==101)  || (PP_TimeDiscMethod==121)
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bs
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
INTEGER,INTENT(IN)           :: iStage
REAL,INTENT(IN)              :: dt
REAL,INTENT(IN)              :: FieldStage (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
REAL,INTENT(IN)              :: FieldSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems,1:5)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL               :: tphi
INTEGER            :: iCounter
!===================================================================================================================================

SELECT CASE(PredictorType)
  CASE(0)
    ! trival guess
    ! nothing to do, because old stage value is used as a prediction
    U=U
  CASE(1)
    ! use RHS as Predictor
    U=LinSolverRHS
  CASE(2)
    ! second order dense output
#if (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==101) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
    !tphi = 1.+RK_c(iStage)
    tphi = RK_c(iStage)
    U=Un
    DO iCounter = 1,iStage-1
      U = U + (RK_bs(iCounter,1)*tphi+RK_bs(iCounter,2)*tphi**2)*dt &
              *(FieldStage (:,:,:,:,:,iCounter) + FieldSource(:,:,:,:,:,iCounter))
    END DO
#else
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)
#endif
  CASE(3)
#if (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==105) || (PP_TimeDiscMethod==122)
    ! third order dense output
   ! tphi = 1.+RK_c(iStage)
    tphi = RK_c(iStage)
    U=Un
    DO iCounter = 1,iStage-1
      U = U + (RK_bsO3(iCounter,1)*tphi+RK_bsO3(iCounter,2)*tphi**2+RK_bsO3(iCounter,3)*tphi**3)*dt &
              *(FieldStage (:,:,:,:,:,iCounter) + FieldSource(:,:,:,:,:,iCounter))
    END DO
#else
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)
#endif
  CASE(4)
    U=2.*Upast(:,:,:,:,:,0)-Upast(:,:,:,:,:,-1)
  CASE DEFAULT
END SELECT

END SUBROUTINE Predictor

SUBROUTINE StorePredictor()
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,          ONLY: U
USE MOD_LinearSolver_Vars,ONLY: Upast
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SELECT CASE(PredictorType)
CASE(4)
  Upast(:,:,:,:,:,-1)=Upast(:,:,:,:,:,0)
  Upast(:,:,:,:,:,0) =U
CASE(5) 
  Upast(:,:,:,:,:,-2)=Upast(:,:,:,:,:,-1)
  Upast(:,:,:,:,:,-1)=Upast(:,:,:,:,:, 0)
  Upast(:,:,:,:,:, 0)=U
END SELECT

END SUBROUTINE StorePredictor

SUBROUTINE FinalizePredictor()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Linearsolver_vars,    ONLY:upast
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(Upast)
END SUBROUTINE FinalizePredictor

END MODULE MOD_Predictor
