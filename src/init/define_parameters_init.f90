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
USE MOD_Output           ,ONLY: DefineParametersOutput
USE MOD_Restart          ,ONLY: DefineParametersRestart
USE MOD_LoadBalance      ,ONLY: DefineParametersLoadBalance
USE MOD_Analyze          ,ONLY: DefineParametersAnalyze
USE MOD_RecordPoints     ,ONLY: DefineParametersRecordPoints
USE MOD_TimeDisc         ,ONLY: DefineParametersTimedisc
USE MOD_Mesh             ,ONLY: DefineparametersMesh
USE MOD_Particle_Mesh    ,ONLY: DefineparametersParticleMesh
USE MOD_Particle_Surfaces,ONLY: DefineParametersParticleSurfaces
USE MOD_Equation         ,ONLY: DefineParametersEquation
USE MOD_PML              ,ONLY: DefineParametersPML
USE MOD_Dielectric       ,ONLY: DefineParametersDielectric
USE MOD_Filter           ,ONLY: DefineParametersFilter
USE MOD_Boltzplatz_Init  ,ONLY: DefineParametersBoltzplatz
USE MOD_ParticleInit     ,ONLY: DefineParametersParticles
!USE MOD_DSMC_Init       ,ONLY: DefineParametersDSMC
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
CALL DefineParametersOutput()
CALL DefineParametersBoltzplatz()
CALL DefineParametersRestart()
CALL DefineParametersLoadBalance()
CALL DefineParametersTimedisc()
CALL DefineparametersMesh()
CALL DefineparametersParticleMesh()
CALL DefineparametersParticleSurfaces()
CALL DefineParametersEquation()
CALL DefineParametersPML()
CALL DefineParametersDielectric()
CALL DefineParametersFilter()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()
CALL DefineParametersParticles()
!CALL DefineParametersDSMC()

END SUBROUTINE InitDefineParameters

END MODULE MOD_Define_Parameters_Init
