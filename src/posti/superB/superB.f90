!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

PROGRAM SuperB_standalone
!===================================================================================================================================
!> Standalone version of SuperB for the calculation of magnetic fields
!===================================================================================================================================
! MODULES
USE MOD_Commandline_Arguments
!USE MOD_Globals               ,ONLY: doPrintHelp,iError,MPIRoot,StartTime,UNIT_stdOut,PiclasTime,SetStackSizeUnlimited
USE MOD_Globals!               ,ONLY: CollectiveStop
USE MOD_Globals_Init          ,ONLY: InitGlobals
USE MOD_SuperB_Init           ,ONLY: DefineParametersSuperB, FinalizeSuperB
USE MOD_SuperB                ,ONLY: SuperB
USE MOD_SuperB_Vars           ,ONLY: BGFieldTDep
USE MOD_Globals_Vars          ,ONLY: ParameterFile
USE MOD_ReadInTools           ,ONLY: prms,PrintDefaultparameterFile,ExtractparameterFile
USE MOD_Interpolation         ,ONLY: InitInterpolation
USE MOD_IO_HDF5               ,ONLY: InitIOHDF5
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_Equation              ,ONLY: InitEquation
USE MOD_Output                ,ONLY: InitOutput
USE MOD_Interpolation         ,ONLY: DefineParametersInterpolation
USE MOD_IO_HDF5               ,ONLY: DefineParametersIO
USE MOD_Output                ,ONLY: DefineParametersOutput
USE MOD_Mesh                  ,ONLY: DefineParametersMesh,FinalizeMesh
USE MOD_Equation              ,ONLY: DefineParametersEquation
USE MOD_Interpolation_Vars    ,ONLY: BGField,BGFieldAnalytic
USE MOD_Mesh                  ,ONLY: InitMesh
#if USE_MPI
USE MOD_MPI_Shared
#endif /*USE_MPI*/
USE MOD_Globals_Init          ,ONLY: DefineParametersGlobals
USE MOD_Mesh_ReadIn           ,ONLY: FinalizeMeshReadin
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: SystemTime
CHARACTER(32)           :: hilf
!===================================================================================================================================
! Initialize
!CALL InitializePiclas()

CALL SetStackSizeUnlimited()

CALL InitMPI()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') "                                                                             D  "
SWRITE(UNIT_stdOut,'(A)') "                                                                      :NNNNNNNN "
SWRITE(UNIT_stdOut,'(A)') "                                                                             D  "
SWRITE(UNIT_stdOut,'(A)') "                                                                                "
SWRITE(UNIT_stdOut,'(A)') "    MMMMMMMMI                                                         MMMMMMMM~ "
SWRITE(UNIT_stdOut,'(A)') "  MM,      D=                                                         MM     MM~"
SWRITE(UNIT_stdOut,'(A)') " ?MZ            MI     MM   MM NMMMM      ,MMMMN     MO +MMD         IM      MM "
SWRITE(UNIT_stdOut,'(A)') " ,MM           MM      M,   MMM   OMM   +MN    MM   MMMM             MM    IMM  "
SWRITE(UNIT_stdOut,'(A)') "   NMMMMMM     MM     MM   +M7     MM   M      MM   MM               MMMMMMMMM  "
SWRITE(UNIT_stdOut,'(A)') "         MM:  IM+     MN   MM      MM  MMMMMMMMMM  :M?              MM       MM "
SWRITE(UNIT_stdOut,'(A)') "         7M:  MM     $M,   MM     OM+  MM          MM               MM       MM "
SWRITE(UNIT_stdOut,'(A)') "M~      ?MM   MM    ?MM   IM     DM7   MM      ?   MM              ZMZ     ~MM  "
SWRITE(UNIT_stdOut,'(A)') "MMMMMMMMM     NMMMMM MM   MMMMMMMM      MMMMMMMN  ZM               MMMMMMMMM    "
SWRITE(UNIT_stdOut,'(A)') "                         M=                                                     "
SWRITE(UNIT_stdOut,'(A)') "                        MM                                                      "
SWRITE(UNIT_stdOut,'(A)') "                        ~~                                                      "
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')"superB version 1.0.0"
SWRITE(UNIT_stdOut,'(132("="))')

GETTIME(StartTime)

CALL ParseCommandlineArguments()


! Check if the number of arguments is correct
IF ((nArgs.GT.1) .OR. ((nArgs.EQ.0).AND.(doPrintHelp.EQ.0)) ) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: superB parameter.ini '// &
    'or superB --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

!CALL InitDefineParameters()

CALL DefineParametersIO()
CALL DefineParametersGlobals()
CALL DefineParametersInterpolation()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEquation()
CALL DefineParametersSuperB()
#if USE_MPI
CALL DefineParametersMPIShared()
#endif /*USE_MPI*/
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

ParameterFile = Args(1)

CALL prms%read_options(ParameterFile)
! Measure init duration
GETTIME(SystemTime)
SWRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT=hilf,FMT='(I0)') prms%count_setentries()
CALL DisplayMessageAndTime(SystemTime-StartTime, ' READING INI DONE! NOW '//TRIM(hilf)//' PARAMETERS ARE SET',&
DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')

CALL InitOutput()
CALL InitIOHDF5()

CALL InitGlobals()
#if USE_MPI
CALL InitMPIShared()
#endif /*USE_MPI*/

! Initialization
CALL InitInterpolation()
CALL InitEquation()

CALL InitMesh(3) ! 0: only read and build Elem_xGP,
                 ! 1: as 0 + build connectivity
                 ! 2: as 1 + calc metrics
                 ! 3: as 2 but skip InitParticleMesh
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate the background B-field via SuperB
!-----------------------------------------------------------------------------------------------------------------------------------
CALL SuperB()

! Deallocation of BGField
SDEALLOCATE(BGFieldTDep)
SDEALLOCATE(BGField)
SDEALLOCATE(BGFieldAnalytic)
! Finalize SuperB
CALL FinalizeSuperB()
CALL FinalizeMesh()
CALL FinalizeMeshReadin(2)

GETTIME(SystemTime)
SWRITE(UNIT_stdOut,'(132("="))')
CALL DisplayMessageAndTime(SystemTime-StartTime, ' SuperB finished!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')
! MPI
#if USE_MPI
! We also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

END PROGRAM SuperB_standalone
