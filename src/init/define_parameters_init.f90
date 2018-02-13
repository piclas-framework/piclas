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
USE MOD_MPI              ,ONLY: DefineParametersMPI
USE MOD_IO_HDF5          ,ONLY: DefineParametersIO
USE MOD_Interpolation    ,ONLY: DefineParametersInterpolation
#if defined(IMEX) || defined(IMPA)
USE MOD_LinearSolver     ,ONLY: DefineParametersLinearSolver
#endif
USE MOD_Output           ,ONLY: DefineParametersOutput
USE MOD_Restart          ,ONLY: DefineParametersRestart
USE MOD_LoadBalance      ,ONLY: DefineParametersLoadBalance
USE MOD_Analyze          ,ONLY: DefineParametersAnalyze
USE MOD_RecordPoints     ,ONLY: DefineParametersRecordPoints
USE MOD_TimeDisc         ,ONLY: DefineParametersTimedisc
USE MOD_Mesh             ,ONLY: DefineparametersMesh
USE MOD_Particle_Mesh    ,ONLY: DefineparametersParticleMesh
USE MOD_Equation         ,ONLY: DefineParametersEquation
#ifndef PP_HDG
USE MOD_PML              ,ONLY: DefineParametersPML
#endif /*PP_HDG*/
#ifdef PP_HDG
USE MOD_HDG              ,ONLY: DefineParametersHDG
#endif /*PP_HDG*/
USE MOD_TTMInit          ,ONLY: DefineParametersTTM
USE MOD_Dielectric       ,ONLY: DefineParametersDielectric
USE MOD_Filter           ,ONLY: DefineParametersFilter
USE MOD_Boltzplatz_Init  ,ONLY: DefineParametersBoltzplatz
USE MOD_ParticleInit     ,ONLY: DefineParametersParticles
USE MOD_Particle_Analyze ,ONLY: DefineParametersParticleAnalyze
USE MOD_PICInit          ,ONLY: DefineParametersPIC
USE MOD_Part_Emission    ,ONLY: DefineParametersParticleEmission
USE MOD_DSMC_Init        ,ONLY: DefineParametersDSMC
USE MOD_LD_Init          ,ONLY: DefineParametersLD
USE MOD_DSMC_SurfModelInit,ONLY: DefineParametersSurfModel
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL DefineParametersMPI()
CALL DefineParametersIO()
CALL DefineParametersInterpolation()
#if defined(IMEX) || defined(IMPA)
CALL DefineParametersLinearSolver()
#endif
CALL DefineParametersOutput()
CALL DefineParametersBoltzplatz()
CALL DefineParametersRestart()
CALL DefineParametersLoadBalance()
CALL DefineParametersTimedisc()
CALL DefineparametersMesh()
CALL DefineparametersParticleMesh()
CALL DefineParametersEquation()
#ifndef PP_HDG
CALL DefineParametersPML()
#endif /*PP_HDG*/
#ifdef PP_HDG
CALL DefineParametersHDG()
#endif /*PP_HDG*/
CALL DefineParametersTTM()
CALL DefineParametersDielectric()
CALL DefineParametersFilter()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()
CALL DefineParametersParticles()
CALL DefineParametersParticleAnalyze()
CALL DefineParametersPIC()
CALL DefineParametersParticleEmission()
CALL DefineParametersDSMC()
CALL DefineParametersLD()
CALL DefineParametersSurfModel()

END SUBROUTINE InitDefineParameters

END MODULE MOD_Define_Parameters_Init
