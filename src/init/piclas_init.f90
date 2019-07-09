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

MODULE MOD_Piclas_Init
!===================================================================================================================================
! contains global init and finalize
!===================================================================================================================================

INTERFACE InitPiclas
  MODULE PROCEDURE InitPiclas
END INTERFACE

INTERFACE FinalizePiclas
  MODULE PROCEDURE FinalizePiclas
END INTERFACE

PUBLIC:: InitPiclas,FinalizePiclas
!===================================================================================================================================
!PUBLIC:: InitDefineParameters

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersPiclas()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Piclas Initialization")

#ifdef PARTICLES
CALL prms%CreateLogicalOption(  'UseDSMC'    , "Flag for using DSMC in Calculation", '.FALSE.')
CALL prms%CreateLogicalOption(  'UseLD'      , "Flag for using LD in Calculation", '.FALSE.')
#endif

END SUBROUTINE DefineParametersPiclas




SUBROUTINE InitPiclas(IsLoadBalance)
!----------------------------------------------------------------------------------------------------------------------------------!
! init Piclas data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,        ONLY:prms
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Restart_Vars,       ONLY:RestartInitIsDone
USE MOD_Restart,            ONLY:InitRestart
USE MOD_Restart_Vars,       ONLY:DoRestart
!USE MOD_ReadInTools,        ONLY:IgnoredStrings
USE MOD_Mesh,               ONLY:InitMesh
USE MOD_Equation,           ONLY:InitEquation
USE MOD_GetBoundaryFlux,    ONLY:InitBC
USE MOD_DG,                 ONLY:InitDG
USE MOD_Mortar,             ONLY:InitMortar
#ifndef PP_HDG
USE MOD_PML,                ONLY:InitPML
#endif /*PP_HDG*/
USE MOD_Dielectric,         ONLY:InitDielectric
USE MOD_Filter,             ONLY:InitFilter
USE MOD_Analyze,            ONLY:InitAnalyze
USE MOD_RecordPoints,       ONLY:InitRecordPoints
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver,       ONLY:InitLinearSolver
#endif /*ROS or IMPA*/
!#ifdef IMEX
!USE MOD_CSR,                ONLY:InitCSR
!#endif /*IMEX*/
USE MOD_Restart_Vars,       ONLY:N_Restart,InterpolateSolution,RestartNullifySolution
#ifdef MPI
USE MOD_MPI,                ONLY:InitMPIvars
#endif /*MPI*/
#ifdef PARTICLES
USE MOD_DSMC_Vars,          ONLY:UseDSMC
USE MOD_LD_Vars,            ONLY:UseLD
USE MOD_ParticleInit,       ONLY:InitParticles
USE MOD_TTMInit,            ONLY:InitTTM,InitIMD_TTM_Coupling
USE MOD_TTM_Vars,           ONLY:DoImportTTMFile
USE MOD_Particle_Surfaces,  ONLY:InitParticleSurfaces
USE MOD_Particle_Mesh,      ONLY:InitParticleMesh, InitElemBoundingBox
USE MOD_Particle_Analyze,   ONLY:InitParticleAnalyze
USE MOD_SurfaceModel_Analyze,ONLY:InitSurfModelAnalyze
USE MOD_Particle_MPI,       ONLY:InitParticleMPI
#if defined(IMPA) || defined(ROS)
USE MOD_ParticleSolver,     ONLY:InitPartSolver
#endif
#endif
#ifdef PP_HDG
USE MOD_HDG,                ONLY:InitHDG
#endif
USE MOD_Interfaces,         ONLY:InitInterfaces
#if USE_QDS_DG
USE MOD_QDS,                ONLY:InitQDS
#endif /*USE_QDS_DG*/
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETREALARRAY
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================


#ifdef PARTICLES
! DSMC handling:
useDSMC=GETLOGICAL('UseDSMC','.FALSE.')
useLD=GETLOGICAL('UseLD','.FALSE.')
IF(useLD) useDSMC=.TRUE.
#endif /*PARTICLES*/

! Initialization
!CALL InitInterpolation()
IF(IsLoadBalance)THEN
  DoRestart=.TRUE.
  RestartInitIsDone=.TRUE.
  InterpolationInitIsDone=.TRUE.
  RestartNullifySolution=.FALSE.
  !BuildNewMesh       =.FALSE. !not used anymore?
  !WriteNewMesh       =.FALSE. !not used anymore?
  InterpolateSolution=.FALSE.
  N_Restart=PP_N
  CALL InitMortar()
ELSE
  CALL InitMortar()
  CALL InitRestart()
END IF
CALL InitMesh()
#ifdef MPI
CALL InitMPIVars()
#endif /*MPI*/
#ifdef PARTICLES
!#ifdef MPI
CALL InitParticleMPI
CALL InitElemBoundingBox()
!#endif /*MPI*/
CALL InitParticleSurfaces()
!CALL InitParticleMesh()
#endif /*PARTICLES*/
CALL InitEquation()
CALL InitBC()
!#ifdef PARTICLES
!CALL InitParticles()
!#endif
#ifndef PP_HDG
CALL InitPML() ! Perfectly Matched Layer (PML): electromagnetic-wave-absorbing layer
#endif /*PP_HDG*/
CALL InitDielectric() ! Dielectric media
CALL InitDG()
CALL InitFilter()
!CALL InitTimeDisc()
#if defined(ROS) || defined(IMPA)
CALL InitLinearSolver()
#endif /*ROS /IMEX*/
!#if defined(IMEX)
!CALL InitCSR()
!#endif /*IMEX*/
#ifdef PARTICLES
CALL InitParticles()
#if defined(IMPA) || defined(ROS)
CALL InitPartSolver()
#endif
!CALL GetSideType
#endif
CALL InitAnalyze()
CALL InitRecordPoints()
#ifdef PARTICLES
CALL InitParticleAnalyze()
CALL InitSurfModelAnalyze()
#endif

#ifdef PP_HDG
CALL InitHDG()
#endif

#ifdef PARTICLES
  CALL InitTTM() ! FD grid based data from a Two-Temperature Model (TTM) from Molecular Dynamics (MD) Code IMD
IF(DoImportTTMFile)THEN
  CALL InitIMD_TTM_Coupling() ! use MD and TTM data to distribute the cell averaged charge to the atoms/ions
END IF
#endif /*PARTICLES*/

CALL InitInterfaces() ! set riemann solver identifier for face connectivity (vacuum, dielectric, PML ...)

#if USE_QDS_DG
CALL InitQDS()
#endif /*USE_QDS_DG*/

! do this last!
!CALL IgnoredStrings()
! write out parameters that are not used and remove multiple and unused, that are not needed to do restart if no parameter.ini is
! read in
IF (.NOT.IsLoadBalance) THEN
  CALL prms%WriteUnused()
  CALL prms%RemoveUnnecessary()
END IF

END SUBROUTINE InitPiclas


SUBROUTINE FinalizePiclas(IsLoadBalance)
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize Piclas data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ReadInTools               ,ONLY:prms
USE MOD_Restart,                   ONLY:FinalizeRestart
USE MOD_Interpolation,             ONLY:FinalizeInterpolation
USE MOD_Mesh,                      ONLY:FinalizeMesh
USE MOD_Equation,                  ONLY:FinalizeEquation
USE MOD_Interfaces,                ONLY:FinalizeInterfaces
#if USE_QDS_DG
USE MOD_QDS,                       ONLY:FinalizeQDS
#endif /*USE_QDS_DG*/
USE MOD_GetBoundaryFlux,           ONLY:FinalizeBC
USE MOD_DG,                        ONLY:FinalizeDG
USE MOD_Mortar,                    ONLY:FinalizeMortar
USE MOD_Dielectric,                ONLY:FinalizeDielectric
#ifndef PP_HDG
USE MOD_PML,                       ONLY:FinalizePML
#else
USE MOD_HDG,                       ONLY:FinalizeHDG
#endif /*PP_HDG*/
USE MOD_Filter,                    ONLY:FinalizeFilter
USE MOD_Analyze,                   ONLY:FinalizeAnalyze
USE MOD_RecordPoints,              ONLY:FinalizeRecordPoints
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver,              ONLY:FinalizeLinearSolver
!USE MOD_CSR,                       ONLY:FinalizeCSR
#endif /*IMEX*/
!USE MOD_TimeDisc,                  ONLY:FinalizeTimeDisc
#ifdef MPI
USE MOD_MPI,                       ONLY:FinalizeMPI
#endif /*MPI*/
#ifdef PARTICLES
USE MOD_Particle_Surfaces,         ONLY:FinalizeParticleSurfaces
USE MOD_InitializeBackgroundField, ONLY:FinalizeBackGroundField
USE MOD_Particle_Mesh,             ONLY:FinalizeParticleMesh
USE MOD_Particle_Analyze,          ONLY:FinalizeParticleAnalyze
USE MOD_PICDepo,                   ONLY:FinalizeDeposition
USE MOD_ParticleInit,              ONLY:FinalizeParticles
USE MOD_TTMInit,                   ONLY:FinalizeTTM
USE MOD_DSMC_Init,                 ONLY:FinalizeDSMC
#if (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Init,               ONLY:FinalizeFPFlow
#endif
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Init,                  ONLY:FinalizeBGK
#endif
USE MOD_SurfaceModel_Init,         ONLY:FinalizeSurfaceModel
USE MOD_Particle_Boundary_Sampling,ONLY:FinalizeParticleBoundarySampling
USE MOD_Particle_Vars,             ONLY:ParticlesInitIsDone
USE MOD_PIC_Vars,                  ONLY:PICInitIsDone
#ifdef MPI
USE MOD_Particle_MPI,              ONLY:FinalizeParticleMPI
USE MOD_Particle_MPI_Vars,         ONLY:ParticleMPIInitisdone
#endif /*MPI*/
#endif /*PARTICLES*/
USE MOD_IO_HDF5,                ONLY:ClearElemData,ElementOut
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CALL ClearElemData(ElementOut)
!Finalize
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
CALL FinalizeDG()
#if defined(IMPA) || defined(ROS)
!CALL FinalizeCSR()
CALL FinalizeLinearSolver()
#endif /*IMEX*/
#ifndef PP_HDG
CALL FinalizePML()
#else
CALL FinalizeDielectric()
CALL FinalizeHDG()
#endif /*PP_HDG*/
CALL FinalizeEquation()
CALL FinalizeBC()
IF(.NOT.IsLoadBalance) CALL FinalizeInterpolation()
!CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeMortar()
CALL FinalizeFilter()
#ifdef PARTICLES
CALL FinalizeSurfaceModel()
CALL FinalizeParticleBoundarySampling()
CALL FinalizeParticleSurfaces()
CALL FinalizeParticleMesh()
CALL FinalizeParticleAnalyze()
CALL FinalizeDeposition()
#ifdef MPI
CALL FinalizeParticleMPI()
#endif /*MPI*/
CALL FinalizeDSMC()
#if (PP_TimeDiscMethod==300)
CALL FinalizeFPFlow()
#endif
#if (PP_TimeDiscMethod==400)
CALL FinalizeBGK()
#endif
CALL FinalizeParticles()
CALL FinalizeBackGroundField()
#endif /*PARTICLES*/
#ifdef MPI
CALL FinalizeMPI()
#endif /*MPI*/

#ifdef PARTICLES
ParticlesInitIsDone = .FALSE.
PICInitIsDone = .FALSE.
#ifdef MPI
ParticleMPIInitIsDone=.FALSE.
#endif /*MPI*/

CALL FinalizeTTM() ! FD grid based data from a Two-Temperature Model (TTM) from Molecular Dynamics (MD) Code IMD
#endif /*PARTICLES*/

CALL FinalizeInterfaces()
#if USE_QDS_DG
CALL FinalizeQDS()
#endif /*USE_QDS_DG*/
CALL prms%finalize(IsLoadBalance)

END SUBROUTINE FinalizePiclas

END MODULE MOD_Piclas_Init
