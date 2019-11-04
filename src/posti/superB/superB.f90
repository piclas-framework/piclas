!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

PROGRAM SuperB_standalone
!===================================================================================================================================
!> Standalone version of SuperB for the calculation of magnetic fields
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_SuperB_Init           ,ONLY: InitializeSuperB
USE MOD_SuperB                ,ONLY: SuperB
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
USE MOD_Mesh                  ,ONLY: DefineParametersMesh
USE MOD_Equation              ,ONLY: DefineParametersEquation
USE MOD_SuperB_Init           ,ONLY: DefineParametersSuperB
USE MOD_StringTools           ,ONLY: STRICMP, GetFileExtension
USE MOD_Interpolation_Vars    ,ONLY: BGField
USE MOD_Mesh                  ,ONLY: InitMesh
#ifdef PARTICLES
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationType
#endif /*PARTICLES*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: SystemTime
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

CALL ParseCommandlineArguments()


! Check if the number of arguments is correct
IF ((nArgs.GT.1) .OR. ((nArgs.EQ.0).AND.(doPrintHelp.EQ.0)) ) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: superB parameter.ini '// &
    'or superB --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

!CALL InitDefineParameters()

CALL DefineParametersIO()
CALL DefineParametersInterpolation()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEquation()
CALL DefineParametersSuperB()
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

ParameterFile = Args(1)

StartTime=PICLASTIME()
CALL prms%read_options(ParameterFile)
! Measure init duration
SystemTime=PICLASTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A,I0,A)') ' READING INI DONE! [',SystemTime-StartTime,' sec ] NOW '&
,prms%count_setentries(),' PARAMETERS ARE SET'
SWRITE(UNIT_stdOut,'(132("="))')

CALL InitOutput()
CALL InitIOHDF5()

CALL InitGlobals()

! Initialization
CALL InitInterpolation()
CALL InitEquation()

CALL InitMesh(0)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate the background B-field via SuperB
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
InterpolationType = 'particle_position'
#endif /*PARTICLES*/
CALL SuperB()

! Deallocation of BGField
SDEALLOCATE(BGField)

SystemTime=PICLASTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A,I0,A)') ' SuperB finished! [',SystemTime-StartTime,' sec ] '
SWRITE(UNIT_stdOut,'(132("="))')
! MPI
#if USE_MPI
! We also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

END PROGRAM SuperB_standalone
