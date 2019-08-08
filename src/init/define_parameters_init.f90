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
#include "piclas.h"

MODULE MOD_Define_Parameters_Init
!===================================================================================================================================
! Initialization of all defined parameters
!===================================================================================================================================

PUBLIC:: InitDefineParameters
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDefineParameters()
!----------------------------------------------------------------------------------------------------------------------------------!
! Calls all parameter definition routines
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_ReadInTools                     ,ONLY: prms
USE MOD_MPI                             ,ONLY: DefineParametersMPI
USE MOD_IO_HDF5                         ,ONLY: DefineParametersIO
USE MOD_Interpolation                   ,ONLY: DefineParametersInterpolation
USE MOD_Output                          ,ONLY: DefineParametersOutput
USE MOD_Restart                         ,ONLY: DefineParametersRestart
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver                    ,ONLY: DefineParametersLinearSolver
#endif
USE MOD_LoadBalance                     ,ONLY: DefineParametersLoadBalance
USE MOD_Analyze                         ,ONLY: DefineParametersAnalyze
USE MOD_RecordPoints                    ,ONLY: DefineParametersRecordPoints
USE MOD_TimeDisc                        ,ONLY: DefineParametersTimedisc
USE MOD_Mesh                            ,ONLY: DefineparametersMesh
USE MOD_Equation                        ,ONLY: DefineParametersEquation
#if !(USE_HDG)
USE MOD_PML                             ,ONLY: DefineParametersPML
#endif /*USE_HDG*/
#if USE_QDS_DG
USE MOD_QDS                             ,ONLY: DefineParametersQDS
#endif
#if USE_HDG
USE MOD_HDG                             ,ONLY: DefineParametersHDG
#endif /*USE_HDG*/
USE MOD_Dielectric                      ,ONLY: DefineParametersDielectric
USE MOD_Filter                          ,ONLY: DefineParametersFilter
USE MOD_Piclas_Init                     ,ONLY: DefineParametersPiclas
#ifdef PARTICLES
USE MOD_ParticleInit                    ,ONLY: DefineParametersParticles
USE MOD_Particle_Mesh                   ,ONLY: DefineparametersParticleMesh
USE MOD_Particle_Analyze                ,ONLY: DefineParametersParticleAnalyze
USE MOD_TTMInit                         ,ONLY: DefineParametersTTM
USE MOD_PICInit                         ,ONLY: DefineParametersPIC
USE MOD_Part_Emission                   ,ONLY: DefineParametersParticleEmission
USE MOD_DSMC_Init                       ,ONLY: DefineParametersDSMC
USE MOD_LD_Init                         ,ONLY: DefineParametersLD
USE MOD_SurfaceModel_Init               ,ONLY: DefineParametersSurfModel
USE MOD_SurfaceModel_Analyze            ,ONLY: DefineParametersSurfModelAnalyze
USE MOD_BGK_Init                        ,ONLY: DefineParametersBGK
USE MOD_FPFlow_Init                     ,ONLY: DefineParametersFPFlow
USE MOD_Particle_Boundary_Porous        ,ONLY: DefineParametersPorousBC
USE MOD_Particle_VarTimeStep            ,ONLY: DefineParametersVaribleTimeStep
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! Insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' DEFINING PARAMETERS ...'
SWRITE(UNIT_stdOut,'(132("="))')

CALL DefineParametersMPI()
CALL DefineParametersIO()
CALL DefineParametersLoadBalance()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
#if defined(ROS) || defined(IMPA)
CALL DefineParametersLinearSolver()
#endif
CALL DefineParametersOutput()
CALL DefineParametersPiclas()
CALL DefineParametersTimedisc()
CALL DefineParametersMesh()
CALL DefineParametersEquation()
#if !(USE_HDG)
CALL DefineParametersPML()
#endif /*USE_HDG*/
#if USE_QDS_DG
CALL DefineParametersQDS()
#endif
#if USE_HDG
CALL DefineParametersHDG()
#endif /*USE_HDG*/
CALL DefineParametersDielectric()
CALL DefineParametersFilter()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()
#ifdef PARTICLES
CALL DefineParametersParticles()
CALL DefineParametersVaribleTimeStep()
CALL DefineParametersPorousBC()
CALL DefineParametersParticleMesh()
CALL DefineParametersParticleAnalyze()
CALL DefineParametersTTM()
CALL DefineParametersPIC()
CALL DefineParametersParticleEmission()
CALL DefineParametersDSMC()
CALL DefineParametersLD()
#if (PP_TimeDiscMethod==300)
CALL DefineParametersFPFlow()
#endif
#if (PP_TimeDiscMethod==400)
CALL DefineParametersBGK()
#endif
CALL DefineParametersSurfModel()
CALL DefineParametersSurfModelAnalyze()
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I0,A)') ' DEFINING PARAMETERS DONE! --> ',prms%count_entries(),' UNIQUE PARAMETERS DEFINED'
SWRITE(UNIT_stdOut,'(132("="))')


END SUBROUTINE InitDefineParameters

END MODULE MOD_Define_Parameters_Init
