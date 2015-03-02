#include "boltzplatz.h"

PROGRAM Boltzplatz
!===================================================================================================================================
! Control program of the Boltzplatz code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,      ONLY:IgnoredStrings
USE MOD_Restart,          ONLY:InitRestart,Restart,FinalizeRestart
USE MOD_Interpolation,    ONLY:InitInterpolation,FinalizeInterpolation
USE MOD_Mesh,             ONLY:InitMesh,FinalizeMesh
USE MOD_Equation,         ONLY:InitEquation,FinalizeEquation
USE MOD_DG,               ONLY:InitDG,FinalizeDG
USE MOD_PML,              ONLY:InitPML,FinalizePML
USE MOD_Filter,           ONLY:InitFilter,FinalizeFilter
USE MOD_Output,           ONLY:InitOutput,FinalizeOutput
USE MOD_IO_HDF5,          ONLY:InitIO
USE MOD_Analyze,          ONLY:InitAnalyze,FinalizeAnalyze
USE MOD_Particle_Analyze, ONLY:InitParticleAnalyze,FinalizeParticleAnalyze
USE MOD_RecordPoints,     ONLY:InitRecordPoints,FinalizeRecordPoints
USE MOD_TimeDisc,         ONLY:InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI,              ONLY:InitMPI
#ifdef MPI
USE MOD_MPI,              ONLY:InitMPIvars
#endif
#ifdef PARTICLES
USE MOD_ParticleInit,     ONLY:InitParticles
USE MOD_Particle_Surfaces,ONLY:InitParticleSurfaces,FinalizeParticleSurfaces, GetSideType
USE MOD_InitializeBackgroundField, ONLY: FinalizeBackGroundField
USE MOD_Particle_Mesh,    ONLY:InitParticleMesh
#ifdef MPI
USE MOD_Particle_MPI,     ONLY:InitParticleCommSize,InitParticleMPI
#endif
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL :: Time
!===================================================================================================================================
CALL InitMPI()
CALL InitIO()
CALL InitGlobals()
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
! Measure init duration
StartTime=BOLTZPLATZTIME()

! Initialization
CALL InitInterpolation()
CALL InitRestart()
CALL InitOutput()
CALL InitMesh()
#ifdef MPI
CALL InitMPIvars()
#endif
#ifdef PARTICLES
! missing
#ifdef MPI
CALL InitParticleMPI
#endif
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
CALL InitTimeDisc()
#ifdef PARTICLES
CALL InitParticles()
CALL GetSideType
#ifdef MPI
CALL InitParticleCommSize()
#endif
CALL InitParticleAnalyze()
#endif
CALL InitAnalyze()
CALL InitRecordPoints()
CALL IgnoredStrings()
CALL Restart()

! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Run Simulation
CALL TimeDisc()

!Finalize
CALL FinalizeOutput()
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
CALL FinalizeDG()
CALL FinalizePML()
CALL FinalizeEquation()
CALL FinalizeInterpolation()
CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeFilter()
#ifdef PARTICLES
CALL FinalizeParticleSurfaces()
CALL FinalizeBackGroundField()
#endif
! Measure simulation duration
Time=BOLTZPLATZTIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError,999.)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' BOLTZPLATZ FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
END PROGRAM Boltzplatz
