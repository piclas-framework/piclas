#include "boltzplatz.h"

PROGRAM Boltzplatz
!===================================================================================================================================
! Control program of the Boltzplatz code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ParameterFile,ParameterDSMCFile
USE MOD_Commandline_Arguments
USE MOD_ReadInTools            ,ONLY: prms,PrintDefaultparameterFile,ExtractparameterFile
USE MOD_Boltzplatz_Init        ,ONLY: InitBoltzplatz,FinalizeBoltzplatz
USE MOD_Restart_Vars           ,ONLY: RestartFile
USE MOD_Restart                ,ONLY: Restart
USE MOD_Interpolation          ,ONLY: InitInterpolation
USE MOD_IO_HDF5                ,ONLY: InitIO
USE MOD_TimeDisc               ,ONLY: InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI                    ,ONLY: InitMPI
USE MOD_RecordPoints_Vars      ,ONLY: RP_Data
USE MOD_Mesh_Vars              ,ONLY: DoSwapMesh
USE MOD_Mesh                   ,ONLY: SwapMesh
#ifdef MPI
USE MOD_LoadBalance            ,ONLY: InitLoadBalance,FinalizeLoadBalance
USE MOD_MPI                    ,ONLY: FinalizeMPI
#endif /*MPI*/
USE MOD_Output                 ,ONLY: InitOutput
USE MOD_Define_Parameters_Init ,ONLY: InitDefineParameters
USE MOD_StringTools            ,ONLY:STRICMP, GetFileExtension
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: Time
LOGICAL                 :: userblockFound
!===================================================================================================================================

CALL InitMPI()

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

CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF ((nArgs.GT.3) .OR. ((nArgs.EQ.0).AND.(doPrintHelp.EQ.0)) ) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
    'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

CALL InitDefineParameters()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

ParameterFile = Args(1)
IF (nArgs.EQ.2) THEN
  ParameterDSMCFile = Args(2)
  IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
      'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "h5")) THEN
    RestartFile = ParameterDSMCFile
    ParameterDSMCFile = '' !'no file found'
  END IF
ELSE IF (nArgs.GT.2) THEN
  ParameterDSMCFile = Args(2)
  IF (STRICMP(GetFileExtension(ParameterDSMCFile), "h5").OR.STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
      'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
ELSE IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
  ! Print out error message containing valid syntax
  !CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
  !  'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  ParameterFile = ".boltzplatz.ini" 
  CALL ExtractParameterFile(Args(1), ParameterFile, userblockFound)
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(Args(1))//"'")
  END IF
  RestartFile = Args(1)
END IF

StartTime=BOLTZPLATZTIME()
CALL prms%read_options(ParameterFile)
! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' READIN INI DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
! Check if we want to read in DSMC.ini
IF(nArgs.GE.2)THEN
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "ini")) THEN
!StartTime=BOLTZPLATZTIME()
    CALL prms%read_options(ParameterDSMCFile,furtherini=.TRUE.)
! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' READIN FURTHER INI DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
  END IF
END IF

CALL InitOutput()
CALL InitIO()

CALL InitGlobals()
#ifdef MPI
CALL InitLoadBalance()
#endif /*MPI*/
! call init routines
! Measure init duration
!StartTime=BOLTZPLATZTIME()

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

