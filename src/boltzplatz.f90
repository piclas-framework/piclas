#include "boltzplatz.h"

PROGRAM Boltzplatz
!===================================================================================================================================
! Control program of the Boltzplatz code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Boltzplatz_Tools,  ONLY:InitBoltzplatz,FinalizeBoltzplatz
USE MOD_Restart,           ONLY:Restart
USE MOD_Interpolation,     ONLY:InitInterpolation!,FinalizeInterpolation
USE MOD_IO_HDF5,           ONLY:InitIO
USE MOD_TimeDisc,          ONLY:InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI,               ONLY:InitMPI
USE MOD_RecordPoints_Vars, ONLY:RP_Data
#ifdef MPI
USE MOD_LoadBalance,       ONLY:InitLoadBalance,FinalizeLoadBalance
USE MOD_MPI,               ONLY:FinalizeMPI
#endif /*MPI*/
#ifdef IMEX
USE MOD_LinearSolver_Vars, ONLY:totalIterLinearSolver
#endif /*IMEX*/
!USE MOD_ReadInTools,      ONLY:IgnoredStrings
!
!USE MOD_Restart,          ONLY:InitRestart,Restart,FinalizeRestart
!USE MOD_Mesh,             ONLY:InitMesh,FinalizeMesh
!USE MOD_Equation,         ONLY:InitEquation,FinalizeEquation
!USE MOD_DG,               ONLY:InitDG,FinalizeDG
!USE MOD_PML,              ONLY:InitPML,FinalizePML
!USE MOD_Filter,           ONLY:InitFilter,FinalizeFilter
!USE MOD_Output,           ONLY:InitOutput,FinalizeOutput
!USE MOD_RecordPoints,     ONLY:InitRecordPoints,FinalizeRecordPoints
!USE MOD_TimeDisc,         ONLY:TimeDisc
!#ifdef MPI
!!USE MOD_MPI,              ONLY:InitMPIvars,FinalizeMPI
!!USE MOD_MPI,              ONLY:FinalizeMPI
!#endif /*MPI*/
!#ifdef PARTICLES
!USE MOD_ParticleInit,     ONLY:InitParticles
!USE MOD_Particle_Surfaces,ONLY:InitParticleSurfaces,FinalizeParticleSurfaces!, GetSideType
!USE MOD_InitializeBackgroundField, ONLY: FinalizeBackGroundField
!USE MOD_Particle_Mesh,    ONLY:InitParticleMesh,FinalizeParticleMesh
!USE MOD_Particle_Analyze, ONLY:InitParticleAnalyze,FinalizeParticleAnalyze
!USE MOD_Particle_MPI,     ONLY:InitParticleMPI
!USE MOD_PICDepo,          ONLY:FinalizeDeposition
!USE MOD_ParticleInit,     ONLY:FinalizeParticles
!USE MOD_DSMC_Init,        ONLY:FinalizeDSMC
!#ifdef MPI
!USE MOD_Particle_MPI,     ONLY:FinalizeParticleMPI
!#endif /*MPI*/
!#endif
!USE MOD_ReadinTools,  ONLY: Readindone
!USE MOD_Particle_MPI_Vars, ONLY: ParticleMPIInitisdone
!USE MOD_Particle_Vars, ONLY: ParticlesInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL    :: Time
!===================================================================================================================================

CALL InitMPI()
CALL InitIO()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')&
 "           ____            ___    __                    ___              __              "
SWRITE(UNIT_stdOut,'(A)')&
 "          /\\  _`\\         /\\_ \\  /\\ \\__                /\\_ \\            /\\ \\__           "
SWRITE(UNIT_stdOut,'(A)')&
 "          \\ \\ \\L\\ \\    ___\\//\\ \\ \\ \\ ,_\\  ____    _____\\//\\ \\       __  \\ \\ ,_\\  ____    "
SWRITE(UNIT_stdOut,'(A)')&
 "           \\ \\  _ <'  / __`\\\\ \\ \\ \\ \\ \\/ /\\_ ,`\\ /\\ '__`\\\\ \\ \\    /'__`\\ \\ \\ \\/ /\\_ ,`\\  "
SWRITE(UNIT_stdOut,'(A)')&
 "            \\ \\ \\L\\ \\/\\ \\L\\ \\\\_\\ \\_\\ \\ \\_\\/_/  /_\\ \\ \\L\\ \\\\_\\ \\_ /\\ \\L\\.\\_\\ \\ \\_\\/_/  /_ "
SWRITE(UNIT_stdOut,'(A)')&
 "             \\ \\____/\\ \\____//\\____\\\\ \\__\\ /\\____\\\\ \\ ,__//\\____\\\\ \\__/.\\_\\\\ \\__\\ /\\____\\"
SWRITE(UNIT_stdOut,'(A)')&
 "              \\/___/  \\/___/ \\/____/ \\/__/ \\/____/ \\ \\ \\/ \\/____/ \\/__/\\/_/ \\/__/ \\/____/"
SWRITE(UNIT_stdOut,'(A)')&
 "                                                    \\ \\_\\                                "
SWRITE(UNIT_stdOut,'(A)')&
 "                                                     \\/_/                                "
SWRITE(UNIT_stdOut,'(A)')&
 ' '
SWRITE(UNIT_stdOut,'(132("="))')
CALL InitGlobals()
#ifdef MPI
CALL InitLoadBalance()
#endif /*MPI*/
! call init routines
! Measure init duration
StartTime=BOLTZPLATZTIME()

! Initialization
CALL InitInterpolation()
CALL InitTimeDisc()

CALL InitBoltzplatz(IsLoadBalance=.FALSE.)
!CALL InitRestart()
!CALL InitOutput()
!CALL InitMesh()
!#ifdef MPI
!CALL InitMPIVars()
!#endif /*MPI*/
!#ifdef PARTICLES
!!#ifdef MPI
!CALL InitParticleMPI
!!#endif /*MPI*/
!CALL InitParticleSurfaces()
!!CALL InitParticleMesh()
!#endif /*PARTICLES*/
!CALL InitEquation()
!!1#ifdef PARTICLES
!!CALL InitParticles()
!!#endif
!CALL InitPML()
!CALL InitDG()
!CALL InitFilter()
!#ifdef PARTICLES
!CALL InitParticles()
!!CALL GetSideType
!#endif
!CALL InitAnalyze()
!CALL InitRecordPoints()
!#ifdef PARTICLES
!CALL InitParticleAnalyze()
!#endif
!CALL IgnoredStrings()

! RESTART
CALL Restart()

! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Run Simulation
CALL TimeDisc()

!Finalize
CALL FinalizeBoltzplatz(IsLoadBalance=.FALSE.)
!CALL FinalizeOutput()
!CALL FinalizeRecordPoints()
!CALL FinalizeAnalyze()
!CALL FinalizeDG()
!CALL FinalizePML()
!CALL FinalizeEquation()
!CALL FinalizeInterpolation()
!CALL FinalizeTimeDisc()
!CALL FinalizeRestart()
!CALL FinalizeMesh()
!CALL FinalizeFilter()
!#ifdef PARTICLES
!CALL FinalizeParticleSurfaces()
!CALL FinalizeParticleMesh()
!CALL FinalizeParticleAnalyze()
!CALL FinalizeDeposition() 
!#ifdef MPI
!CALL FinalizeParticleMPI()
!#endif /*MPI*/
!CALL FinalizeDSMC()
!CALL FinalizeParticles()
!CALL FinalizeBackGroundField()
!#endif
!#ifdef MPI
!CALL FinalizeMPI()
!#endif

CALL FinalizeTimeDisc()
! mssing arrays to deallocate
SDEALLOCATE(RP_Data)

!Measure simulation duration
Time=BOLTZPLATZTIME()

#ifdef MPI
!! and additional required for restart with load balance
!ReadInDone=.FALSE.
!ParticleMPIInitIsDone=.FALSE.
!ParticlesInitIsDone=.FALSE.
CALL FinalizeLoadBalance()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(&
  __STAMP__&
  ,'MPI finalize error',iError,999.)
#endif
#ifdef IMEX
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,I12)') ' Total iteration Linear Solver ',totalIterLinearSolver
#if (PP_TimeDiscMethod==104)
SWRITE(UNIT_stdOut,'(A,I12)') ' Total iteration Newton        ',nNewton
SWRITE(UNIT_stdOut,'(132("="))')
#endif
#endif /*IMEX*/
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)')  ' BOLTZPLATZ FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Boltzplatz
