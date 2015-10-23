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
USE MOD_DG,                        ONLY:FinalizeDG
USE MOD_PML,                       ONLY:FinalizePML
USE MOD_Filter,                    ONLY:FinalizeFilter
USE MOD_Output,                    ONLY:FinalizeOutput
USE MOD_Analyze,                   ONLY:FinalizeAnalyze
USE MOD_RecordPoints,              ONLY:FinalizeRecordPoints
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
USE MOD_DSMC_Init,                 ONLY:FinalizeDSMC
#ifdef MPI
USE MOD_Particle_MPI,              ONLY:FinalizeParticleMPI
#endif /*MPI*/
#endif /*PARTICLES*/
USE MOD_ReadinTools,  ONLY: Readindone
USE MOD_Particle_MPI_Vars, ONLY: ParticleMPIInitisdone
USE MOD_Particle_Vars, ONLY: ParticlesInitIsDone

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
CALL FinalizeOutput()
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
CALL FinalizeDG()
CALL FinalizePML()
CALL FinalizeEquation()
IF(.NOT.IsLoadBalance) CALL FinalizeInterpolation()
!CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
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
CALL FinalizeMPI()
#endif

ReadInDone=.FALSE.
ParticleMPIInitIsDone=.FALSE.
ParticlesInitIsDone=.FALSE.  

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
USE MOD_Restart,          ONLY:InitRestart
USE MOD_Restart_Vars,     ONLY:DoRestart
USE MOD_ReadInTools,      ONLY:IgnoredStrings
!USE MOD_Interpolation,    ONLY:InitInterpolation
USE MOD_Mesh,             ONLY:InitMesh
USE MOD_Equation,         ONLY:InitEquation
USE MOD_DG,               ONLY:InitDG
USE MOD_PML,              ONLY:InitPML
USE MOD_Filter,           ONLY:InitFilter
USE MOD_Output,           ONLY:InitOutput
USE MOD_Analyze,          ONLY:InitAnalyze
USE MOD_RecordPoints,     ONLY:InitRecordPoints
!USE MOD_TimeDisc,         ONLY:InitTimeDisc
USE MOD_Restart_Vars,     ONLY:N_Restart,InterpolateSolution
#ifdef MPI
USE MOD_MPI,              ONLY:InitMPIvars
#endif /*MPI*/
#ifdef PARTICLES
USE MOD_ParticleInit,     ONLY:InitParticles
USE MOD_Particle_Surfaces,ONLY:InitParticleSurfaces
USE MOD_Particle_Mesh,    ONLY:InitParticleMesh
USE MOD_Particle_Analyze, ONLY:InitParticleAnalyze
USE MOD_Particle_MPI,     ONLY:InitParticleMPI
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
ELSE
  CALL InitRestart()
END IF
CALL InitOutput()
CALL InitMesh()
#ifdef MPI
CALL InitMPIVars()
#endif /*MPI*/
#ifdef PARTICLES
!#ifdef MPI
CALL InitParticleMPI
!#endif /*MPI*/
CALL InitParticleSurfaces()
!CALL InitParticleMesh()
#endif /*PARTICLES*/
CALL InitEquation()
!1#ifdef PARTICLES
!CALL InitParticles()
!#endif
CALL InitPML()
CALL InitDG()
CALL InitFilter()
!CALL InitTimeDisc()
#ifdef PARTICLES
CALL InitParticles()
!CALL GetSideType
#endif
CALL InitAnalyze()
CALL InitRecordPoints()
#ifdef PARTICLES
CALL InitParticleAnalyze()
#endif
CALL IgnoredStrings()


END SUBROUTINE InitBoltzplatz


END MODULE MOD_Boltzplatz_tools
