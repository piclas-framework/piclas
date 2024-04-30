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

INTERFACE PrintStatusLine
  MODULE PROCEDURE PrintStatusLine
END INTERFACE

INTERFACE PrintStatusLineRadiation
  MODULE PROCEDURE PrintStatusLineRadiation
END INTERFACE


PUBLIC:: InitOutput
PUBLIC:: PrintStatusLine
PUBLIC:: DefineParametersOutput
PUBLIC:: PrintStatusLineRadiation
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
CALL prms%CreateStringOption(  'ProjectName'        , "Name of the current simulation (mandatory).")
CALL prms%CreateLogicalOption( 'Logging'            , "Write log files containing debug output."     , '.FALSE.')
CALL prms%CreateLogicalOption( 'WriteErrorFiles'    , "Write error files containing error output."   , '.FALSE.')
CALL prms%CreateLogicalOption( 'doPrintStatusLine'  , 'Print: percentage of time                     , ...'       , '.FALSE.')
CALL prms%CreateLogicalOption( 'DoWriteStateToHDF5' , "Write state of calculation to hdf5-file."     , '.TRUE.')
END SUBROUTINE DefineParametersOutput

SUBROUTINE InitOutput()
!===================================================================================================================================
! Initialize all output (and analyze) variables.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars ,ONLY: ParameterFile,ProjectName,ParameterDSMCFile
USE MOD_Preproc
USE MOD_ReadInTools  ,ONLY: GetStr,GetLogical,GETINT
USE MOD_Output_Vars
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

ProjectName        = GETSTR('ProjectName')
Logging            = GETLOGICAL('Logging')
WriteErrorFiles    = GETLOGICAL('WriteErrorFiles')
doPrintStatusLine  = GETLOGICAL("doPrintStatusLine")
DoWriteStateToHDF5 = GETLOGICAL('DoWriteStateToHDF5')

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


!==================================================================================================================================
!> Displays the actual status of the simulation
!==================================================================================================================================
SUBROUTINE PrintStatusLine(t,dt,tStart,tEnd,mode)
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars , ONLY: doPrintStatusLine
USE MOD_TimeDisc_Vars,ONLY: time_start
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: t      ! < current simulation time
REAL,INTENT(IN)    :: dt     ! < current time step
REAL,INTENT(IN)    :: tStart ! < start time of simulation
REAL,INTENT(IN)    :: tEnd   ! < end time of simulation
INTEGER,INTENT(IN) :: mode   ! < mode: 1: only print when doPrintStatusLine=T, 2: print during dt_analyze
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: percent,time_remaining,mins,secs,hours,days
CHARACTER(LEN=60)  :: hilf
!==================================================================================================================================

IF(mode.EQ.1)THEN
  IF(.NOT.doPrintStatusLine) RETURN
  hilf = ACHAR(13) ! ACHAR(13) is carriage return
ELSEIF(mode.EQ.2)THEn
  IF(doPrintStatusLine) RETURN
  IF(ALMOSTEQUALRELATIVE(t,tEnd,1e-5)) RETURN
  hilf = ACHAR(10) ! ACHAR(10) is Line feed (newline)
ELSE
  CALL abort(__STAMP__,'wrong mode = ',IntInfoOpt=mode)
END IF ! mode.EQ.1

IF(MPIroot)THEN
#ifdef INTEL
  OPEN(UNIT_stdOut,CARRIAGECONTROL='fortran')
#endif
  percent = (t-tStart) / (tEnd-tStart)
  CALL CPU_TIME(time_remaining)
  time_remaining = time_remaining - time_start
  IF (percent.GT.0.0) time_remaining = time_remaining/percent - time_remaining
  percent = percent*100.
  secs = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  mins = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  hours = MOD(time_remaining,24.)
  time_remaining = time_remaining / 24
  !days = MOD(time_remaining,365.) ! Use this if years are also to be displayed
  days = time_remaining
#if USE_COFFEE
  WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,I6,A1,I0.2,A1,I0.2,A1,I0.2,A,A,A,A3,F6.2,A3,A1)',ADVANCE='NO') &
      '  Time = ', t,'  dt = ', dt, ' ', ' eta = ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),'     |',&
      REPEAT('☕',CEILING(percent/2)),REPEAT(' ',INT((100-percent)/2)),'| [',percent,'%] ', TRIM(hilf)
#else
  WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,I6,A1,I0.2,A1,I0.2,A1,I0.2,A,A,A1,A,A3,F6.2,A3,A1)',ADVANCE='NO') &
      '  Time = ', t,'  dt = ', dt, ' ', ' eta = ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),'     |',&
      REPEAT('=',MAX(CEILING(percent/2)-1,0)),'>',REPEAT(' ',INT((100-percent)/2)),'| [',percent,'%] ', TRIM(hilf)
#endif /*USE_COFFEE*/
#ifdef INTEL
  CLOSE(UNIT_stdOut)
#endif
END IF
END SUBROUTINE PrintStatusLine


!==================================================================================================================================
!> Displays the actual status of the simulation
!==================================================================================================================================
SUBROUTINE PrintStatusLineRadiation(t,tStart,tEnd,Phot,outputrank)
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: t      !< current simulation time
REAL,INTENT(IN) :: tStart !< start time of simulation
REAL,INTENT(IN) :: tEnd   !< end time of simulation
LOGICAL, INTENT(IN) :: Phot
INTEGER, INTENT(IN),OPTIONAL :: outputrank
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: percent,time_remaining,mins,secs,hours,days
INTEGER :: visRank
!==================================================================================================================================
IF (PRESENT(outputrank)) THEN
  visRank = outputrank
ELSE
  visRank = 0
END IF

IF(myRank.EQ.visRank)THEN
#ifdef INTEL
  OPEN(UNIT_stdOut,CARRIAGECONTROL='fortran')
#endif
  percent = (t-tStart) / (tend-tStart)
  CALL CPU_TIME(time_remaining)
  IF (percent.GT.0.0) time_remaining = time_remaining/percent - time_remaining
  percent =  percent*100.
  secs = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  mins = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  hours = MOD(time_remaining,24.)
  time_remaining = time_remaining / 24
  !days = MOD(time_remaining,365.) ! Use this if years are also to be displayed
  days = time_remaining
  IF (Phot) THEN
    WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,I6,A1,I0.2,A1,I0.2,A1,I0.2,A,A,A,A4,I3,A4,A1)',ADVANCE='NO') &
        '  Photon = ', t,'  TotalPhotons = ', tEnd, ' ', ' eta = ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),'     |',&
        REPEAT('☄️ ',CEILING(percent/4)),REPEAT('  ',INT((100-percent)/4)),'| [ ',NINT(percent),'% ] ',&
        ACHAR(13) ! ACHAR(13) is carriage return
  ELSE
    WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,A,I6,A1,I0.2,A1,I0.2,A1,I0.2,A,A,A,A4,I3,A4,A1)',ADVANCE='NO') &
        '  Elem = ', t,'  TotalElems = ', tEnd, ' ', ' eta = ',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),'     |',&
        REPEAT('☢ ',CEILING(percent/4)),REPEAT('  ',INT((100-percent)/4)),'| [ ',NINT(percent),'% ] ',&
        ACHAR(13) ! ACHAR(13) is carriage return
  END IF
#ifdef INTEL
  CLOSE(UNIT_stdOut)
#endif
END IF
END SUBROUTINE PrintStatusLineRadiation

END MODULE MOD_Output
