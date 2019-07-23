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

MODULE MOD_Output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER :: OUTPUTFORMAT_NONE         = 0
INTEGER,PARAMETER :: OUTPUTFORMAT_TECPLOT      = 1
INTEGER,PARAMETER :: OUTPUTFORMAT_TECPLOTASCII = 2
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW     = 3

! Output format for ASCII data files
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_CSV     = 0
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_TECPLOT = 1
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE
  FUNCTION get_userblock_size()
      INTEGER :: get_userblock_size
  END FUNCTION
END INTERFACE

INTERFACE
  FUNCTION get_inifile_size(filename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR,C_INT
      CHARACTER(KIND=C_CHAR) :: filename(*)
      INTEGER(KIND=C_INT)    :: get_inifile_size
  END FUNCTION get_inifile_size
END INTERFACE

INTERFACE
  SUBROUTINE insert_userblock(filename,filename2,inifilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: filename(*)
      CHARACTER(KIND=C_CHAR) :: filename2(*)
      CHARACTER(KIND=C_CHAR) :: inifilename(*)
  END SUBROUTINE insert_userblock
END INTERFACE

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,FinalizeOutput

PUBLIC::DefineParametersOutput
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersOutput()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Output")
!CALL prms%CreateIntOption(          'NVisu',       "Polynomial degree at which solution is sampled for visualization.")
!CALL prms%CreateIntOption(          'NOut',        "Polynomial degree at which solution is written. -1: NOut=N, >0: NOut", '-1')
CALL prms%CreateStringOption(       'ProjectName', "Name of the current simulation (mandatory).")
CALL prms%CreateLogicalOption(      'Logging',     "Write log files containing debug output.", '.FALSE.')
CALL prms%CreateLogicalOption(      'WriteErrorFiles',  "Write error files containing error output.", '.TRUE.')
CALL prms%CreateIntFromStringOption('OutputFormat',"File format for visualization: None, Tecplot, TecplotASCII, ParaView. "//&
                                                 " Note: Tecplot output is currently unavailable due to licensing issues.", 'None')
CALL addStrListEntry('OutputFormat','none',        OUTPUTFORMAT_NONE)
CALL addStrListEntry('OutputFormat','tecplot',     OUTPUTFORMAT_TECPLOT)
CALL addStrListEntry('OutputFormat','tecplotascii',OUTPUTFORMAT_TECPLOTASCII)
CALL addStrListEntry('OutputFormat','paraview',    OUTPUTFORMAT_PARAVIEW)
CALL prms%CreateIntFromStringOption('ASCIIOutputFormat',"File format for ASCII files, e.g. body forces: CSV, Tecplot."&
                                                       , 'CSV')
CALL addStrListEntry('ASCIIOutputFormat','csv',    ASCIIOUTPUTFORMAT_CSV)
CALL addStrListEntry('ASCIIOutputFormat','tecplot',ASCIIOUTPUTFORMAT_TECPLOT)
CALL prms%CreateLogicalOption(      'doPrintStatusLine','Print: percentage of time, ...', '.FALSE.')
CALL prms%CreateLogicalOption(      'WriteStateFiles','Write HDF5 state files. Disable this only for debugging issues. \n'// &
                                                      'NO SOLUTION WILL BE WRITTEN!', '.TRUE.')
END SUBROUTINE DefineParametersOutput

SUBROUTINE InitOutput()
!===================================================================================================================================
! Initialize all output (and analyze) variables.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars, ONLY: ParameterFile,ProjectName,ParameterDSMCFile
USE MOD_Preproc
USE MOD_ReadInTools,ONLY:GetStr,GetLogical,GETINT
USE MOD_Output_Vars,ONLY:OutputInitIsDone,OutputFormat
USE MOD_Output_Vars,ONLY:userblock_len, userblock_total_len, UserBlockTmpFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
LOGICAL                        :: WriteErrorFiles
LOGICAL                        :: LogIsOpen 
INTEGER                        :: inifile_len
!===================================================================================================================================
IF (OutputInitIsDone) THEN
  CALL abort(&
      __STAMP__&
      ,'InitOutput not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

ProjectName=GETSTR('ProjectName')
Logging    =GETLOGICAL('Logging','.FALSE.')

WriteErrorFiles=GETLOGICAL('WriteErrorFiles','.FALSE.')
IF(WriteErrorFiles)THEN
  ! Open file for error output
  WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'
  OPEN(UNIT=UNIT_errOut,  &
       FILE=ErrorFileName,&
       STATUS='REPLACE',  &
       ACTION='WRITE',    &
       IOSTAT=OpenStat)
ELSE
  UNIT_errOut=UNIT_stdOut
END IF

IF (MPIRoot) THEN
  ! read userblock length in bytes from data section of piclas-executable
  userblock_len = get_userblock_size()
  inifile_len = get_inifile_size(TRIM(ParameterFile)//C_NULL_CHAR)
  ! prepare userblock file
  CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(ParameterFile)//C_NULL_CHAR,TRIM(ParameterDSMCFile)//C_NULL_CHAR)
  INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)
END IF

OutputFormat = GETINT('OutputFormat','1')
! Open file for logging
IF(Logging)THEN
  INQUIRE(UNIT=UNIT_LogOut,OPENED=LogIsOpen)
  IF(.NOT.LogIsOpen)THEN
    WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
    OPEN(UNIT=UNIT_logOut,  &
         FILE=LogFile,      &
         STATUS='REPLACE',  &
         ACTION='WRITE',    &
         IOSTAT=OpenStat)
    CALL DATE_AND_TIME(StrDate,StrTime)
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,'(132("#"))')
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                        StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
  END IF !logIsOpen
END IF  ! Logging

OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput



SUBROUTINE InitOutputBasis(N_in,NVisu_in,xGP,wBary)
!===================================================================================================================================
! Initialize all output variables.
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars, ONLY:Vdm_GaussN_NVisu
USE MOD_Basis,       ONLY:InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: N_in,NVisu_in
REAL,INTENT(IN),DIMENSION(0:N_in)   :: xGP,wBary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:NVisu_in)  :: XiVisu
INTEGER                     :: i
!===================================================================================================================================
!equidistant visu points
DO i=0,NVisu_in
  XiVisu(i) = 2./REAL(NVisu_in) * REAL(i) - 1.
END DO
! Gauss/Gl -> Visu : computation -> visualization
ALLOCATE(Vdm_GaussN_NVisu(0:NVisu_in,0:N_in))
CALL InitializeVandermonde(N_in,NVisu_in,wBary,xGP,XiVisu,Vdm_GaussN_NVisu)
END SUBROUTINE InitOutputBasis



SUBROUTINE FinalizeOutput()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:Vdm_GaussN_NVisu,OutputInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NVisu)
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
