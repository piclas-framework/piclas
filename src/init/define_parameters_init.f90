#include "boltzplatz.h"

MODULE MOD_Define_Parameters_Init
!===================================================================================================================================
! initialization of all defined parameters
!===================================================================================================================================

PUBLIC:: InitDefineParameters
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDefineParameters() 
!----------------------------------------------------------------------------------------------------------------------------------!
! Calls all parameter definition routines
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Restart         ,ONLY: DefineParametersRestart
USE MOD_Analyze         ,ONLY: DefineParametersAnalyze
USE MOD_RecordPoints    ,ONLY: DefineParametersRecordPoints
USE MOD_TimeDisc        ,ONLY: DefineParametersTimedisc
USE MOD_Boltzplatz_Init ,ONLY: DefineParametersBoltzplatz
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL DefineParametersBoltzplatz()
CALL DefineParametersRestart()
CALL DefineParametersTimedisc()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()
END SUBROUTINE InitDefineParameters

END MODULE MOD_Define_Parameters_Init
