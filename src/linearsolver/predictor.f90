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

#if defined(PARTICLES) && IMPA
INTERFACE PartPredictor
  MODULE PROCEDURE PartPredictor
END INTERFACE
#endif /*PARTICLES*/

PUBLIC:: InitPredictor, FinalizePredictor,Predictor,StorePredictor
#if defined(PARTICLES) && IMPA
PUBLIC:: PartPredictor
#endif /*PARTICLES*/

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
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bsO3,RK_bs,RK_b
#endif
#if (PP_TimeDiscMethod==101)  || (PP_TimeDiscMethod==121)
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bs,RK_b
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
#if (PP_TimeDiscMethod==101) || (PP_TimeDiscMethod==121)
REAL               :: tphi2
#endif
#if (PP_TimeDiscMethod==102) || (PP_TimeDiscMethod==122)
REAL               :: tphi2,tphi3
#endif
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
    tphi2= tphi*tphi
    U=(RK_bs(iStage-1,1)*tphi+RK_bs(iStage-1,2)*tphi2) & 
     *(FieldStage (:,:,:,:,:,iStage-1) + FieldSource(:,:,:,:,:,iStage-1))
    DO iCounter = 1,iStage-2
      U = U + (RK_bs(iCounter,1)*tphi+RK_bs(iCounter,2)*tphi2) &
              *(FieldStage (:,:,:,:,:,iCounter) + FieldSource(:,:,:,:,:,iCounter))
    END DO
    U=Un+dt*U
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
    tphi2= tphi*tphi
    tphi3= tphi*tphi2
    U=(RK_bsO3(iStage-1,1)*tphi+RK_bsO3(iStage-1,2)*tphi2+RK_bsO3(iStage-1,3)*tphi3) &
     *(FieldStage (:,:,:,:,:,iStage-1) + FieldSource(:,:,:,:,:,iStage-1))
    DO iCounter = 1,iStage-2
      U = U + (RK_bsO3(iCounter,1)*tphi+RK_bsO3(iCounter,2)*tphi2+RK_bsO3(iCounter,3)*tphi3) &
              *(FieldStage (:,:,:,:,:,iCounter) + FieldSource(:,:,:,:,:,iCounter))
    END DO
    U=Un+dt*U
#else
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)
#endif
  CASE(4)
    U=2.*Upast(:,:,:,:,:,0)-Upast(:,:,:,:,:,-1)
  CASE DEFAULT
END SELECT

! disable warnings
IF(1.EQ.2)THEN
  iCounter=iStage
  iCounter=INT(tPhi)
  tphi    =Un(1,1,1,1,1)
  tphi    =FieldSource(1,1,1,1,1,1)
  tphi    =FieldStage(1,1,1,1,1,1)
  tphi    =dt
END IF

END SUBROUTINE Predictor

#if defined(PARTICLES) && IMPA
SUBROUTINE PartPredictor(iStage,dt,PartID)
!===================================================================================================================================
! predicts the new Stage-value to decrease computational time
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars,    ONLY: PartState
USE MOD_Particle_Vars,    ONLY: PartStage, PartQ,PartStateN
#if (PP_TimeDiscMethod==122)
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bsO3,RK_bs,RK_b
#endif
#if (PP_TimeDiscMethod==121)
USE MOD_TimeDisc_Vars,    ONLY: RK_c,RK_bs,RK_b
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: iStage
REAL,INTENT(IN)              :: dt
INTEGER,INTENT(IN)           :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL               :: tphi
#if (PP_TimeDiscMethod==121) 
REAL               :: tphi2
#endif
#if (PP_TimeDiscMethod==122) 
REAL               :: tphi3,tphi2
#endif
INTEGER            :: iCounter
!===================================================================================================================================

SELECT CASE(PredictorType)
  CASE(0)
    ! trival guess
    ! nothing to do, because old stage value is used as a prediction
    PartState(PartID,:)=PartState(PartID,:)
  CASE(1)
    ! use RHS as Predictor
    PartState(PartID,:)=PartQ(:,PartID)
  CASE(2)
    ! second order dense output
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
    !tphi = 1.+RK_c(iStage)
    tphi = RK_c(iStage)
    tphi2= tphi*tphi 
    PartState(PartID,1:6)=(RK_bs(iStage-1,1)*tphi+RK_bs(iStage-1,2)*tphi2)*PartStage(PartID,1:6,iStage-1)
    DO iCounter=1,iStage-2
      PartState(PartID,1:6) = PartState(PartID,1:6) + &
                          (RK_bs(iCounter,1)*tphi+RK_bs(iCounter,2)*tphi2)*PartStage(PartID,1:6,iCounter)
    END DO
    PartState(PartID,1:6)=PartStateN(PartID,1:6)+dt*PartState(PartID,1:6)
#else
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)
#endif
  CASE(3)
#if (PP_TimeDiscMethod==122)
    ! third order dense output
   ! tphi = 1.+RK_c(iStage)
    tphi = RK_c(iStage)
    tphi2= tphi*tphi
    tphi3= tphi*tphi2
    PartState(PartID,1:6)=(RK_bsO3(iStage-1,1)*tphi+RK_bsO3(iStage-1,2)*tphi2+RK_bsO3(iStage-1,3)*tphi3) &
                         * PartStage(PartID,1:6,iStage-1)
    DO iCounter = 1,iStage-2
      PartState(PartID,1:6) = PartState(PartID,1:6) &
                            + (RK_bsO3(iCounter,1)*tphi+RK_bsO3(iCounter,2)*tphi2+RK_bsO3(iCounter,3)*tphi3) &
                            * (PartStage(PartID,1:6,iCounter) )
    END DO
    PartState(PartID,1:6)=PartStateN(PartID,1:6)+dt*PartState(PartID,1:6)
#else
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)
#endif
  CASE DEFAULT
   CALL abort(&
__STAMP__&
,'No Predictor for this timedisc!',999,999.)

END SELECT

END SUBROUTINE PartPredictor
#endif /*PARTICLES + IMPA*/


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
