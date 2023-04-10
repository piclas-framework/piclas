!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"
#include "commit.h"

MODULE MOD_Piclas

IMPLICIT NONE
PRIVATE
SAVE

INTERFACE InitializePiclas
   MODULE PROCEDURE InitializePiclas
END INTERFACE


PUBLIC::InitializePiclas

CONTAINS
!===================================================================================================================================
! Control program of the Piclas code. Initialization of the computation
!===================================================================================================================================
SUBROUTINE InitializePiclas()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ParameterFile,ParameterDSMCFile
USE MOD_Globals_Vars           ,ONLY: InitializationWallTime,MajorVersion,MinorVersion,PatchVersion
USE MOD_Commandline_Arguments
USE MOD_Globals                ,ONLY: iError,Logging,MPIroot,StartTime,UNIT_stdOut,PiclasTime,doPrintHelp,abort
USE MOD_Globals                ,ONLY: SetStackSizeUnlimited,CollectiveStop,ReOpenLogFile
USE MOD_Globals_Init           ,ONLY: InitGlobals
USE MOD_Globals_Vars           ,ONLY: ParameterFile,ParameterDSMCFile,InitializationWallTime
USE MOD_ReadInTools            ,ONLY: prms,PrintDefaultparameterFile,ExtractparameterFile
USE MOD_Piclas_Init            ,ONLY: InitPiclas
USE MOD_Restart_Vars           ,ONLY: RestartFile
USE MOD_Restart                ,ONLY: Restart
USE MOD_Interpolation          ,ONLY: InitInterpolation
USE MOD_IO_HDF5                ,ONLY: InitIOHDF5
USE MOD_TimeDiscInit           ,ONLY: InitTime,InitTimeDisc
USE MOD_MPI                    ,ONLY: InitMPI
USE MOD_Mesh_Vars              ,ONLY: DoSwapMesh
USE MOD_Mesh                   ,ONLY: SwapMesh
#if USE_MPI
USE MOD_MPI_Shared!            ,ONLY: InitMPIShared
USE MOD_LoadBalance            ,ONLY: InitLoadBalance
#endif /*USE_MPI*/
USE MOD_Output                 ,ONLY: InitOutput
USE MOD_Define_Parameters_Init ,ONLY: InitDefineParameters
USE MOD_StringTools            ,ONLY: STRICMP, GetFileExtension
#ifdef PARTICLES
USE MOD_Particle_Vars          ,ONLY: DoInitialIonization
USE MOD_ParticleInit           ,ONLY: InitialIonization
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: SystemTime
LOGICAL                 :: userblockFound
!===================================================================================================================================
CALL SetStackSizeUnlimited()

CALL InitMPI()

#if defined(MEASURE_MPI_WAIT)
IF((MPIW8SIZEFIELD+MPIW8SIZEPART).EQ.0) CALL abort(__STAMP__,'PICLAS_MEASURE_MPI_WAIT=T not implemented for this system')
#endif /*defined(MEASURE_MPI_WAIT)*/

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(1))
#endif /*EXTRAE*/

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')"                                        _______ _________ _______  _        _______  _______ "
SWRITE(UNIT_stdOut,'(A)')"                                       (  ____ )\__   __/(  ____ \( \      (  ___  )(  ____ \"
SWRITE(UNIT_stdOut,'(A)')"                                       | (    )|   ) (   | (    \/| (      | (   ) || (    \/"
SWRITE(UNIT_stdOut,'(A)')"                                       | (____)|   | |   | |      | |      | (___) || (_____ "
SWRITE(UNIT_stdOut,'(A)')"                                       |  _____)   | |   | |      | |      |  ___  |(_____  )"
SWRITE(UNIT_stdOut,'(A)')"                                       | (         | |   | |      | |      | (   ) |      ) |"
SWRITE(UNIT_stdOut,'(A)')"                                       | )      ___) (___| (____/\| (____/\| )   ( |/\____) |"
SWRITE(UNIT_stdOut,'(A)')"                                       |/       \_______/(_______/(_______/|/     \|\_______)"
SWRITE(UNIT_stdOut,'(132(" "))')
SWRITE(UNIT_stdOut,'(A)')"piclas version "&
    //TRIM(int2strf(MajorVersion))//"."//TRIM(int2strf(MinorVersion))//"."//TRIM(int2strf(PatchVersion))&
    //" with commit "//TRIM(GIT_CURRENT_COMMIT)
SWRITE(UNIT_stdOut,'(132("="))')

CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF ((nArgs.GT.3) .OR. ((nArgs.EQ.0).AND.(doPrintHelp.EQ.0)) ) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: piclas parameter.ini [DSMC.ini] [restart.h5]'// &
    'or piclas --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
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
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: piclas parameter.ini [DSMC.ini] [restart.h5]'// &
      'or piclas --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "h5")) THEN
    RestartFile = ParameterDSMCFile
    ParameterDSMCFile = '' !'no file found'
  END IF
ELSE IF (nArgs.GT.2) THEN
  ParameterDSMCFile = Args(2)
  RestartFile = Args(3)
  IF (STRICMP(GetFileExtension(ParameterDSMCFile), "h5").OR.STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: piclas parameter.ini [DSMC.ini] [restart.h5]'// &
      'or piclas --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
ELSE IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
  ! Print out error message containing valid syntax
  !CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: piclas parameter.ini [DSMC.ini] [restart.h5]'// &
  !  'or piclas --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  ParameterFile = ".piclas.ini"
  CALL ExtractParameterFile(Args(1), ParameterFile, userblockFound)
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(Args(1))//"'")
  END IF
  RestartFile = Args(1)
END IF

GETTIME(StartTime)
CALL prms%read_options(ParameterFile)
! Measure init duration
GETTIME(SystemTime)
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A,I0,A)') ' READING INI DONE! [',SystemTime-StartTime,' sec ] NOW '&
,prms%count_setentries(),' PARAMETERS ARE SET'
SWRITE(UNIT_stdOut,'(132("="))')
! Check if we want to read in DSMC.ini
IF(nArgs.GE.2)THEN
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "ini")) THEN
    CALL prms%read_options(ParameterDSMCFile,furtherini=.TRUE.)
    ! Measure init duration
    GETTIME(SystemTime)
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A,F14.2,A,I0,A)') ' READING FURTHER INI DONE! [',SystemTime-StartTime,' sec ] NOW '&
    ,prms%count_setentries(),' PARAMETERS ARE SET'
    SWRITE(UNIT_stdOut,'(132("="))')
  END IF
END IF

CALL InitOutput()
CALL InitIOHDF5()

CALL InitGlobals()
#if USE_MPI
CALL InitLoadBalance()
CALL InitMPIShared()
#endif /*USE_MPI*/
! call init routines
! Measure init duration
!StartTime=PICLASTIME()

! Initialization
CALL InitInterpolation()
CALL InitTimeDisc()

CALL InitPiclas(IsLoadBalance=.FALSE.)

! Do SwapMesh
IF(DoSwapMesh)THEN
  ! Measure init duration
  GETTIME(SystemTime)
  IF(MPIroot)THEN
    Call SwapMesh()
    CALL DisplayMessageAndTime(SystemTime-StartTime, 'SWAPMESH DONE! PICLAS DONE', DisplayDespiteLB=.TRUE.)
    STOP
  ELSE
    CALL abort(__STAMP__,'DO NOT CALL SWAPMESH WITH MORE THAN 1 Procs!',iError,999.)
  END IF
END IF
LOGWRITE_BARRIER

! The beginning of time
CALL InitTime()

! RESTART
CALL Restart()

#ifdef PARTICLES
! Ionize the current particles
IF(DoInitialIonization) CALL InitialIonization()
#endif /*PARTICLES*/

! Measure init duration
GETTIME(SystemTime)
InitializationWallTime=SystemTime-StartTime
SWRITE(UNIT_stdOut,'(132("="))')
CALL DisplayMessageAndTime(InitializationWallTime, 'INITIALIZATION DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
END SUBROUTINE InitializePiclas


END MODULE MOD_Piclas
