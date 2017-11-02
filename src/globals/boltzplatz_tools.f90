#include "boltzplatz.h"

MODULE MOD_Boltzplatz_tools
!===================================================================================================================================
! contains global init and finalize
!===================================================================================================================================

INTERFACE InitBoltzplatz
  MODULE PROCEDURE InitBoltzplatz
END INTERFACE

INTERFACE FinalizeBoltzplatz
  MODULE PROCEDURE FinalizeBoltzplatz
END INTERFACE

PUBLIC:: InitBoltzplatz,FinalizeBoltzplatz
!===================================================================================================================================
CONTAINS

SUBROUTINE FinalizeBoltzplatz(IsLoadBalance) 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize Boltzplatz data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
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
#ifdef IMEX
USE MOD_LinearSolver,              ONLY:FinalizeLinearSolver
USE MOD_CSR,                       ONLY:FinalizeCSR
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
USE MOD_Particle_Vars,             ONLY:ParticlesInitIsDone
#ifdef MPI
USE MOD_Particle_MPI,              ONLY:FinalizeParticleMPI
USE MOD_Particle_MPI_Vars,         ONLY:ParticleMPIInitisdone
#endif /*MPI*/
#endif /*PARTICLES*/
USE MOD_ReadinTools,               ONLY: Readindone
#if (PP_TimeDiscMethod==104)
USE MOD_LinearSolver_Vars,ONLY:nNewton
#endif

!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

!Finalize
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
CALL FinalizeDG()
#ifdef IMEX
CALL FinalizeCSR()
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
CALL FinalizeParticleSurfaces()
CALL FinalizeParticleMesh()
CALL FinalizeParticleAnalyze()
CALL FinalizeDeposition() 
#ifdef MPI
CALL FinalizeParticleMPI()
#endif /*MPI*/
CALL FinalizeDSMC()
CALL FinalizeParticles()
CALL FinalizeBackGroundField()
#endif /*PARTICLES*/
#ifdef MPI
CALL FinalizeMPI()
#endif /*MPI*/

ReadInDone=.FALSE.
#ifdef PARTICLES
ParticlesInitIsDone=.FALSE.  
#ifdef MPI
ParticleMPIInitIsDone=.FALSE.
#endif /*MPI*/

CALL FinalizeTTM() ! FD grid based data from a Two-Temperature Model (TTM) from Molecular Dynamics (MD) Code IMD
#endif /*PARTICLES*/

CALL FinalizeInterfaces()
#if USE_QDS_DG
CALL FinalizeQDS()
#endif /*USE_QDS_DG*/

END SUBROUTINE FinalizeBoltzplatz


SUBROUTINE InitBoltzplatz(IsLoadBalance) 
!----------------------------------------------------------------------------------------------------------------------------------!
! init Boltzplatz data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone
USE MOD_Restart_Vars,       ONLY:RestartInitIsDone
USE MOD_Restart,            ONLY:InitRestart
USE MOD_Restart_Vars,       ONLY:DoRestart
USE MOD_ReadInTools,        ONLY:IgnoredStrings
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
#if defined(IMEX) || defined(IMPA)
USE MOD_LinearSolver,       ONLY:InitLinearSolver
#endif /*IMEX*/
#ifdef IMEX
USE MOD_CSR,                ONLY:InitCSR
#endif /*IMEX*/
USE MOD_Restart_Vars,       ONLY:N_Restart,InterpolateSolution
#ifdef MPI
USE MOD_MPI,                ONLY:InitMPIvars
#endif /*MPI*/
#ifdef PARTICLES
USE MOD_ParticleInit,       ONLY:InitParticles
USE MOD_TTMInit,            ONLY:InitTTM,InitIMD_TTM_Coupling
USE MOD_TTM_Vars,           ONLY:DoImportTTMFile
USE MOD_Particle_Surfaces,  ONLY:InitParticleSurfaces
USE MOD_Particle_Mesh,      ONLY:InitParticleMesh, InitElemBoundingBox
USE MOD_Particle_Analyze,   ONLY:InitParticleAnalyze
USE MOD_Particle_MPI,       ONLY:InitParticleMPI
#if defined(IMPA) || (PP_TimeDiscMethod==110)
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
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Initialization
!CALL InitInterpolation()
IF(IsLoadBalance)THEN
  DoRestart=.TRUE.
  RestartInitIsDone=.TRUE.
  InterpolationInitIsDone=.TRUE.
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
#if defined(IMEX) || defined(IMPA)
CALL InitLinearSolver()
#endif /*IMEX*/
#if defined(IMEX)
CALL InitCSR()
#endif /*IMEX*/
#ifdef PARTICLES
CALL InitParticles()
#if defined(IMPA) 
CALL InitPartSolver()
#endif
!CALL GetSideType
#endif
CALL InitAnalyze()
CALL InitRecordPoints()
#ifdef PARTICLES
CALL InitParticleAnalyze()
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
CALL IgnoredStrings()

END SUBROUTINE InitBoltzplatz


END MODULE MOD_Boltzplatz_tools
