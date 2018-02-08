#include "boltzplatz.h"

PROGRAM Boltzplatz
!===================================================================================================================================
! Control program of the Boltzplatz code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Boltzplatz_Init   ,ONLY: InitBoltzplatz,FinalizeBoltzplatz,InitDefineParameters
USE MOD_Restart           ,ONLY: Restart
USE MOD_Interpolation     ,ONLY: InitInterpolation
USE MOD_IO_HDF5           ,ONLY: InitIO
USE MOD_TimeDisc          ,ONLY: InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI               ,ONLY: InitMPI
USE MOD_RecordPoints_Vars ,ONLY: RP_Data
USE MOD_Mesh_Vars         ,ONLY: DoSwapMesh
USE MOD_Mesh              ,ONLY: SwapMesh
#ifdef MPI
USE MOD_LoadBalance       ,ONLY: InitLoadBalance,FinalizeLoadBalance
USE MOD_MPI               ,ONLY: FinalizeMPI
#endif /*MPI*/
USE MOD_Output            ,ONLY: InitOutput

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL    :: Time
!===================================================================================================================================

CALL InitMPI()

CALL InitDefineParameters()

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

CALL InitOutput()
CALL InitIO()

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

! Do SwapMesh
IF(DoSwapMesh)THEN
  ! Measure init duration
  Time=BOLTZPLATZTIME()
  IF(MPIroot)THEN
    Call SwapMesh()
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' SWAPMESH DONE! BOLTZPLATZ DONE! [',Time-StartTime,' sec ]'
    SWRITE(UNIT_stdOut,'(132("="))')
    STOP
  ELSE
  CALL abort(&
  __STAMP__&
  ,'DO NOT CALL SWAPMESH WITH MORE THAN 1 Procs!',iError,999.)
  END IF
END IF

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
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)')  ' BOLTZPLATZ FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Boltzplatz
