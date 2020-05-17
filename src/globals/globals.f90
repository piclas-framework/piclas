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

MODULE MOD_Globals
!===================================================================================================================================
!> Provides parameters, used globally (please use EXTREMELY carefully!)
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER  :: UNIT_stdOut=6
INTEGER,PARAMETER  :: UNIT_logOut=133
INTEGER            :: UNIT_errOut=999
LOGICAL            :: Logging
CHARACTER(LEN=255) :: LogFile
CHARACTER(LEN=255) :: ErrorFileName='NOT_SET'
INTEGER            :: iError
REAL               :: StartTime
INTEGER            :: myRank,myLocalRank,myLeaderRank,myWorkerRank
INTEGER            :: nProcessors,nLocalProcs,nLeaderProcs,nWorkerProcs
INTEGER            :: MPI_COMM_NODE    ! local node subgroup
INTEGER            :: MPI_COMM_LEADERS ! all node masters
INTEGER            :: MPI_COMM_WORKERS ! all non-master nodes
LOGICAL            :: MPIRoot,MPILocalRoot
#if USE_MPI
!#include "mpif.h"
INTEGER            :: MPIStatus(MPI_STATUS_SIZE)
#else
INTEGER,PARAMETER  :: MPI_COMM_WORLD=-1 ! DUMMY when compiling single (MPI=OFF)
#endif

INTEGER            :: doPrintHelp ! 0: no help, 1: help, 2: markdown-help

#ifdef INTKIND8
INTEGER, PARAMETER :: IK = SELECTED_INT_KIND(18)
#else
INTEGER, PARAMETER :: IK = SELECTED_INT_KIND(8)
#endif

INTERFACE InitGlobals
  MODULE PROCEDURE InitGlobals
END INTERFACE

INTERFACE ReOpenLogFile
  MODULE PROCEDURE ReOpenLogFile
END INTERFACE

INTERFACE Abort
  MODULE PROCEDURE AbortProg
END INTERFACE Abort

INTERFACE CollectiveStop
  MODULE PROCEDURE CollectiveStop
END INTERFACE CollectiveStop

INTERFACE PrintWarning
  MODULE PROCEDURE PrintWarning
END INTERFACE PrintWarning

INTERFACE FILEEXISTS
  MODULE PROCEDURE FILEEXISTS
END INTERFACE FILEEXISTS

INTERFACE INTSTAMP
  MODULE PROCEDURE INTSTAMP
END INTERFACE INTSTAMP

INTERFACE TIMESTAMP
  MODULE PROCEDURE TIMESTAMP
END INTERFACE

INTERFACE PICLASTIME
  MODULE PROCEDURE PICLASTIME
END INTERFACE

INTERFACE LOCALTIME
  MODULE PROCEDURE LOCALTIME
END INTERFACE

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE CreateErrFile
  MODULE PROCEDURE CreateErrFile
END INTERFACE CreateErrFile

INTERFACE CROSS
  MODULE PROCEDURE CROSS
END INTERFACE CROSS

INTERFACE
  SUBROUTINE setstacksizeunlimited() BIND(C)
  END SUBROUTINE setstacksizeunlimited
END INTERFACE

INTERFACE str2real
  MODULE PROCEDURE str2real
END INTERFACE

INTERFACE str2int
  MODULE PROCEDURE str2int
END INTERFACE

INTERFACE int2str
  MODULE PROCEDURE int2str
END INTERFACE

INTERFACE int2strf
  MODULE PROCEDURE int2strf
END INTERFACE

INTERFACE str2logical
  MODULE PROCEDURE str2logical
END INTERFACE

INTERFACE GetParameterFromFile
  MODULE PROCEDURE GetParameterFromFile
END INTERFACE

INTERFACE SphericalCoordinates
  MODULE PROCEDURE SphericalCoordinates
END INTERFACE

INTERFACE TransformVectorfieldSphericalCoordinates
  MODULE PROCEDURE TransformVectorfieldSphericalCoordinates
END INTERFACE

INTERFACE TransformVectorFromSphericalCoordinates
  MODULE PROCEDURE TransformVectorFromSphericalCoordinates
END INTERFACE

PUBLIC :: setstacksizeunlimited

!===================================================================================================================================
CONTAINS

SUBROUTINE InitGlobals()
!===================================================================================================================================
! Pre-compute required constants
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
LOGICAL                        :: LogIsOpen
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS ...'

PI=ACOS(-1.)
sPI = 1./PI

! get machine accuracy
epsMach=EPSILON(0.0)
TwoEpsMach=2.0d0*epsMach

! Open file for logging
IF(Logging)THEN
  INQUIRE(UNIT=UNIT_LogOut,OPENED=LogIsOpen)
  IF(.NOT.LogIsOpen)THEN
    WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
    OPEN(UNIT=UNIT_logOut,  &
         FILE=LogFile,      &
         STATUS='UNKNOWN',  &
         ACTION='WRITE',    &
         POSITION='APPEND', &
         IOSTAT=OpenStat)
    CALL DATE_AND_TIME(StrDate,StrTime)
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,'(132("#"))')
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                        StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
  END IF !logIsOpen
END IF  ! Logging

SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitGlobals


SUBROUTINE ReOpenLogFile()
!===================================================================================================================================
! re-open log file (used by preprocessor LOGWRITE_BARRIER) to be sure that all logwrites are written to file
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
LOGICAL                        :: LogIsOpen
!===================================================================================================================================
  INQUIRE(UNIT=UNIT_LogOut,OPENED=LogIsOpen)
  IF(logIsOpen)CLOSE(UNIT_logOut)
  OPEN(UNIT=UNIT_logOut, FILE=LogFile, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND', IOSTAT=OpenStat)
END SUBROUTINE ReOpenLogFile


! FUNCTION AlmostEqual(Num1,Num2) ! see piclas.h
! !===================================================================================================================================
! ! Bruce Dawson quote:
! ! "There is no silver bullet. You have to choose wisely."
! !    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
! !      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
! !      to your calculation. Maybe."
! !    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
! !      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
! !      An absolute epsilon could be used if you knew exactly what number you were comparing against."
! !    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
! !      Good luck and God speed."
! !===================================================================================================================================
! ! MODULES
! USE MOD_Globals_Vars,    ONLY:TwoEpsMach ! relative epsilon value: something like 4.???E-16 for double precision
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! REAL            :: Num1,Num2      ! Number
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! LOGICAL         :: ALMOSTEQUAL
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! !===================================================================================================================================
! IF(ABS(Num1-Num2).LE.MAX(ABS(Num1),ABS(Num2))*TwoEpsMach)THEN
!   ALMOSTEQUAL=.TRUE.
! ELSE
!   ALMOSTEQUAL=.FALSE.
! END IF
! END FUNCTION AlmostEqual


! FUNCTION ALMOSTEQUALRELATIVE(Num1,Num2,Tolerance) ! old name "AlmostEqualToTolerance", new is same as for flexi: see piclas.h
! !===================================================================================================================================
! ! Bruce Dawson quote:
! ! "There is no silver bullet. You have to choose wisely."
! !    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
! !      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
! !      to your calculation. Maybe."
! !    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
! !      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
! !      An absolute epsilon could be used if you knew exactly what number you were comparing against."
! !    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
! !      Good luck and God speed."
! !===================================================================================================================================
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! REAL            :: Num1,Num2
! REAL            :: Tolerance ! relative epsilon value as input
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! LOGICAL         :: ALMOSTEQUALRELATIVE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! !===================================================================================================================================
! IF(ABS(Num1-Num2).LE.MAX(ABS(Num1),ABS(Num2))*Tolerance)THEN
!    ALMOSTEQUALRELATIVE=.TRUE.
! ELSE
!   ALMOSTEQUALRELATIVE=.FALSE.
! END IF
! END FUNCTION ALMOSTEQUALRELATIVE


! FUNCTION AlmostZero(Num) ! see piclas.h
! !===================================================================================================================================
! ! Performe an almost zero check. But ...
! ! Bruce Dawson quote:
! ! "There is no silver bullet. You have to choose wisely."
! !    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
! !      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
! !      to your calculation. Maybe."
! !    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
! !      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
! !      An absolute epsilon could be used if you knew exactly what number you were comparing against."
! !    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
! !      Good luck and God speed."
! !===================================================================================================================================
! ! MODULES
! USE MOD_Globals_Vars,    ONLY:EpsMach
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! REAL            :: Num ! Number
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! LOGICAL         :: AlmostZero
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! !===================================================================================================================================
!
! AlmostZero=.FALSE.
! IF(ABS(Num).LE.EpsMach) AlmostZero=.TRUE.
!
! END FUNCTION AlmostZero

#if USE_MPI
SUBROUTINE AbortProg(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfoOpt,RealInfoOpt,SingleOpt)
#else
SUBROUTINE AbortProg(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfoOpt,RealInfoOpt)
#endif
!===================================================================================================================================
! Terminate program correctly if an error has occurred (important in MPI mode!).
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      ! Source file where error has occurred
INTEGER                           :: SourceLine      ! Line in source file
CHARACTER(LEN=*)                  :: CompDate        ! Compilation date
CHARACTER(LEN=*)                  :: CompTime        ! Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    ! Error message
INTEGER,OPTIONAL                  :: IntInfoOpt      ! Error info (integer)
REAL,OPTIONAL                     :: RealInfoOpt     ! Error info (real)
#if USE_MPI
LOGICAL,OPTIONAL                  :: SingleOpt       ! Only MPI-Root performs check
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   There is no way back!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: IntInfo         ! Error info (integer)
REAL                              :: RealInfo        ! Error info (real)
#if USE_MPI
INTEGER                           :: errOut          ! Output of MPI_ABORT
INTEGER                           :: signalout       ! Output errorcode
#endif /*USE_MPI*/
!===================================================================================================================================
IF(logging) CLOSE(UNIT_logOut)
#if USE_MPI
IF(PRESENT(SingleOpt))THEN
  IF(SingleOpt.AND.(.NOT.MPIRoot)) RETURN
END IF
#endif
IF(PRESENT(IntInfoOpt))THEN
  IntInfo=IntInfoOpt
ELSE
  IntInfo=999
END IF
IF(PRESENT(RealInfoOpt))THEN
  RealInfo=RealInfoOpt
ELSE
  RealInfo=999.
END IF
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,*)'_____________________________________________________________________________'
WRITE(UNIT_stdOut,*)'Program abort caused on Proc ',myRank,' in File : ',TRIM(SourceFile),' Line ',SourceLine
WRITE(UNIT_stdOut,*)'This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime)
WRITE(UNIT_stdOut,'(A10,A)',ADVANCE='NO')'Message: ',TRIM(ErrorMessage)
IF(PRESENT(IntInfoOpt)) WRITE(UNIT_stdOut,'(I8)',ADVANCE='NO')IntInfo
IF(PRESENT(RealInfoOpt)) WRITE(UNIT_stdOut,'(E16.8)')RealInfo
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,'(A,A,A)')'See ',TRIM(ErrorFileName),' for more details'
WRITE(UNIT_stdOut,*)
!CALL delete()
#if USE_MPI
signalout=2 ! MPI_ABORT requires an output error-code /=0
errOut = 1
CALL MPI_ABORT(MPI_COMM_WORLD,signalout,errOut)
#endif
STOP 2
END SUBROUTINE AbortProg


!==================================================================================================================================
!> print a warning to the command line (only MPI root)
!==================================================================================================================================
SUBROUTINE PrintWarning(msg)
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*) :: msg
!===================================================================================================================================
IF (myRank.EQ.0) THEN
  WRITE(UNIT_stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(UNIT_stdOut,*) 'WARNING:'
  WRITE(UNIT_stdOut,*) TRIM(msg)
  WRITE(UNIT_stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
END IF
END SUBROUTINE PrintWarning


!==================================================================================================================================
!> \brief Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!>
!> Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> Terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> This routine can only be used if ALL processes are guaranteed to generate the same error at the same time!
!> Prime use is to exit FLEXI without MPI errors and with a single error message if some parameters are not set in the init
!> routines or a file is not found.
!>
!> Criteria where CollectiveStop may be used:
!> 0. In case of doubt stick with Abort, which is always safe!
!> 1. A routine is BY DESIGN (!) called by all processes, i.e. does not permit to be called by single processes or subgroups.
!> 2. The criteria for the CollectiveStop must be identical among all processors.
!> 3. The routine is only used during the init phase.
!> 4. The error must not originate from MPI errors (e.g. during MPI init)
!> 5. The error must not originate from checking roundof errors (e.g. accuracy of interpolation matrices)
!>
!==================================================================================================================================
SUBROUTINE CollectiveStop(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      !< Source file where error has occurred
INTEGER                           :: SourceLine      !< Line in source file
CHARACTER(LEN=*)                  :: CompDate        !< Compilation date
CHARACTER(LEN=*)                  :: CompTime        !< Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    !< Error message
INTEGER,OPTIONAL                  :: IntInfo         !< Error info (integer)
REAL,OPTIONAL                     :: RealInfo        !< Error info (real)
!   There is no way back!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
!==================================================================================================================================
IntString = ""
RealString = ""

IF (PRESENT(IntInfo))  WRITE(IntString,"(A,I0)")  "\nIntInfo:  ", IntInfo
IF (PRESENT(RealInfo)) WRITE(RealString,"(A,F24.19)") "\nRealInfo: ", RealInfo

SWRITE(UNIT_stdOut,*) '_____________________________________________________________________________\n', &
                     'Program abort caused on Proc ',myRank, '\n', &
                     '  in File : ',TRIM(SourceFile),' Line ',SourceLine, '\n', &
                     '  This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime), '\n', &
                     'Message: ',TRIM(ErrorMessage), &
                     TRIM(IntString), TRIM(RealString)

CALL FLUSH(UNIT_stdOut)
#if USE_MPI
CALL MPI_FINALIZE(iError)
#endif
ERROR STOP 1
END SUBROUTINE CollectiveStop


SUBROUTINE CreateErrFile()
!===================================================================================================================================
! Open file for error output
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
LOGICAL                        :: isOpen
!===================================================================================================================================
INQUIRE(UNIT=UNIT_errOut,OPENED=isOpen)
IF(.NOT.isOpen)THEN
  OPEN(UNIT=UNIT_errOut,  &
       FILE=ErrorFileName,&
       STATUS='REPLACE',  &
       ACTION='WRITE',    &
       IOSTAT=OpenStat)
END IF
END SUBROUTINE CreateErrFile


!==================================================================================================================================
!> Convert a String to an Integer
!==================================================================================================================================
SUBROUTINE str2int(str,int_number,stat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=*),INTENT(IN) :: str
INTEGER,INTENT(OUT)         :: int_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  int_number
END SUBROUTINE str2int


!==================================================================================================================================
!> Convert a String to an Integer
!==================================================================================================================================
SUBROUTINE int2str(str,int_number,stat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=255),INTENT(OUT) :: str
INTEGER,INTENT(IN)             :: int_number
INTEGER,INTENT(OUT)            :: stat
!===================================================================================================================================
WRITE(str,'(I0)',IOSTAT=stat)  int_number
END SUBROUTINE int2str


!==================================================================================================================================
!> Convert an Integer to a String
!==================================================================================================================================
!SUBROUTINE int2strf(str,int_number,stat)
FUNCTION int2strf(int_number)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=3) :: int2strf
INTEGER,INTENT(IN) :: int_number
!===================================================================================================================================
WRITE(int2strf,'(I0)')  int_number
int2strf = TRIM(ADJUSTL(int2strf))
END FUNCTION


!==================================================================================================================================
!> Convert a String to a REAL
!==================================================================================================================================
SUBROUTINE str2real(str,real_number,stat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=*),INTENT(IN) :: str
REAL,INTENT(OUT)            :: real_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  real_number
END SUBROUTINE str2real


!==================================================================================================================================
!> Convert a String to a LOGICAL
!==================================================================================================================================
SUBROUTINE str2logical(str,logical_number,stat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(len=*),INTENT(IN) :: str
LOGICAL,INTENT(OUT)         :: logical_number
INTEGER,INTENT(OUT)         :: stat
!===================================================================================================================================
READ(str,*,IOSTAT=stat)  logical_number
END SUBROUTINE str2logical


!==================================================================================================================================
!> read compile flags from a specified file
!> example line in "configuration.cmake": SET(PICLAS_EQNSYSNAME "maxwell" CACHE STRING "Used equation system")
!> ParameterName: timestep
!> output: 0.1
!> Type of Msg: [G]et[P]arameter[F]rom[File] -> GPFF: not ordinary read-in tool
!==================================================================================================================================
SUBROUTINE GetParameterFromFile(FileName,ParameterName,output,DelimiterSymbolIN,CommentSymbolIN,DoDisplayInfo)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: FileName          !> e.g. './../laser.inp'
CHARACTER(LEN=*),INTENT(IN)          :: ParameterName     !> e.g. 'timestep'
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: DelimiterSymbolIN !> e.g. '=' (default is '=')
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: CommentSymbolIN   !> e.g. '#' (default is '!')
CHARACTER(LEN=*),INTENT(INOUT)       :: output            !> e.g. '0.1'
LOGICAL,OPTIONAL,INTENT(IN)          :: DoDisplayInfo     !> default is: TRUE
                                                          !> display DefMsg or errors if the parameter or the file is not found
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                              :: ExistFile         !> file exists=.true., file does not exist=.false.
INTEGER                              :: iSTATUS           !> status
CHARACTER(LEN=255)                   :: temp,temp2,temp3  !> temp variables for read in of file lines
CHARACTER(LEN=255)                   :: DelimiterSymbol   !> symbol for commenting out code, e.g., "#" or "!"
CHARACTER(LEN=255)                   :: CommentSymbol     !> symbol for commenting out code, e.g., "#" or "!"
INTEGER                              :: ioUnit            !> field handler unit and ??
INTEGER                              :: IndNum            !> Index Number
CHARACTER(LEN=8)                     :: DefMsg            !> additional flag like "DEFAULT" or "*CUSTOM"
!===================================================================================================================================
IF(PRESENT(DelimiterSymbolIN))THEN
  DelimiterSymbol=TRIM(ADJUSTL(DelimiterSymbolIN))
ELSE
  DelimiterSymbol='='
END IF
IF(PRESENT(CommentSymbolIN))THEN
  CommentSymbol=TRIM(ADJUSTL(CommentSymbolIN))
ELSE
  CommentSymbol='!'
END IF
output=''
! read from file
INQUIRE(File=TRIM(FileName),EXIST=ExistFile)
IF(ExistFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ')
  DO
    READ(ioUnit,'(A)',iostat=iSTATUS)temp
    IF(ADJUSTL(temp(1:LEN(TRIM(CommentSymbol)))).EQ.TRIM(CommentSymbol)) CYCLE  ! complete line is commented out
    IF(iSTATUS.EQ.-1)EXIT                           ! end of file is reached
    IF(LEN(trim(temp)).GT.1)THEN                    ! exclude empty lines
      IndNum=INDEX(temp,TRIM(ParameterName))        ! e.g. 'timestep'
      IF(IndNum.GT.0)THEN
        IF(IndNum-1.GT.0)THEN                       ! check if the parameter name is contained within a substring of another
          IF(temp(IndNum-1:IndNum-1).NE.' ')CYCLE   ! parameter, e.g., "timestep" within "fd_timestep" -> skip
        END IF
        temp2=TRIM(ADJUSTL(temp(IndNum+LEN(TRIM(ParameterName)):LEN(temp))))
        IF(DelimiterSymbol.NE.'')THEN               ! demiliting symbol must not be empty
          IndNum=INDEX(temp2,TRIM(DelimiterSymbol)) ! only use string FROM delimiting symbol +1
          IF(IndNum.GT.0)THEN
            temp3=TRIM(ADJUSTL(temp2(IndNum+1:LEN(temp2))))
            temp2=temp3
          END IF
        ELSE
          ! no nothing?
        END IF
        IndNum=INDEX(temp2,TRIM(CommentSymbol)) ! only use string UP TO commenting symbol
        IF(IndNum.EQ.0)IndNum=LEN(TRIM(temp2))+1
        output=TRIM(ADJUSTL(temp2(1:IndNum-1)))
        DefMsg='GPFF'
        SWRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(ParameterName),' | ', TRIM(output),' | ',TRIM(DefMsg),' | '
        EXIT ! found the parameter -> exit loop
      END IF
    END IF
  END DO
  CLOSE(ioUnit)
  IF(output.EQ.'')THEN
    IF(PRESENT(DoDisplayInfo))THEN
      IF(DoDisplayInfo)THEN
        SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: Parameter ['//TRIM(ParameterName)//'] not found.'
      END IF
    ELSE
      SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: Parameter ['//TRIM(ParameterName)//'] not found.'
    END IF
    output='ParameterName does not exist'
  END IF
ELSE
  IF(PRESENT(DoDisplayInfo))THEN
    IF(DoDisplayInfo)THEN
      SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: File ['//TRIM(FileName)//'] not found.'
    END IF
  ELSE
    SWRITE(UNIT_stdOut,'(A)') ' SUBROUTINE GetParameterFromFile: File ['//TRIM(FileName)//'] not found.'
  END IF
  output='file does not exist'
END IF
END SUBROUTINE GetParameterFromFile


!==================================================================================================================================
!> Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!==================================================================================================================================
FUNCTION FILEEXISTS(filename)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename
LOGICAL                     :: FILEEXISTS
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
INQUIRE(FILE=TRIM(filename), EXIST=FILEEXISTS)
END FUNCTION FILEEXISTS


FUNCTION INTSTAMP(Nam,Num)
!===================================================================================================================================
! Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)   :: Nam      ! Name
INTEGER            :: Num      ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=200) :: IntStamp ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(IntStamp,'(A,A5,I6.6)')TRIM(Nam),'_Proc',Num
END FUNCTION INTSTAMP



FUNCTION TIMESTAMP(Filename,Time)
!===================================================================================================================================
! Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals_Vars ,ONLY: TimeStampLenStr,TimeStampLenStr2
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)   :: Filename  ! (file)name
REAL               :: Time      ! time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255) :: TimeStamp ! the complete timestamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i         ! loop variable
!===================================================================================================================================
!WRITE(TimeStamp,'(F21.17)')Time
WRITE(TimeStamp,'(F'//TRIM(TimeStampLenStr)//'.'//TRIM(TimeStampLenStr2)//')')Time

! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP


#if USE_MPI
FUNCTION PICLASTIME(Comm)
#else
FUNCTION PICLASTIME()
#endif
!===================================================================================================================================
! Calculates current time (own function because of a laterMPI implementation)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_MPI
INTEGER, INTENT(IN),OPTIONAL    :: Comm
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: PiclasTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
PiclasTime=MPI_WTIME()
#else
CALL CPU_TIME(PiclasTime)
#endif
END FUNCTION PICLASTIME


FUNCTION LOCALTIME()
!===================================================================================================================================
! Calculates current local time (own function because of a laterMPI implementation)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: LocalTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
LocalTime=MPI_WTIME()
#else
CALL CPU_TIME(LocalTime)
#endif
END FUNCTION LOCALTIME


FUNCTION GETFREEUNIT()
!===================================================================================================================================
! Get unused file unit number
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetFreeUnit ! File unit number
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!===================================================================================================================================
GetFreeUnit=55
INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
IF(connected)THEN
  DO
    GetFreeUnit=GetFreeUnit+1
    INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
    IF(.NOT.connected)EXIT
  END DO
END IF
END FUNCTION GETFREEUNIT


PURE FUNCTION CROSS(v1,v2)
!===================================================================================================================================
! Computes the cross product of two 3-dimensional vectors: cross=v1 x v2
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)
REAL,INTENT(IN) :: v2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: CROSS(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS


PURE FUNCTION CROSSNORM(v1,v2)
!===================================================================================================================================
! Computes the cross product of to 3 dimensional vectors: cross=v1 x v2
! and normalizes the vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3),v2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: CROSSNORM(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: length
!===================================================================================================================================
CROSSNORM=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
length=SQRT(CROSSNORM(1)*CROSSNORM(1)+CROSSNORM(2)*CROSSNORM(2)+CROSSNORM(3)*CROSSNORM(3))
CROSSNORM=CROSSNORM/length
END FUNCTION CROSSNORM


PURE FUNCTION UNITVECTOR(v1)
!===================================================================================================================================
! compute  a unit vector from a given vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: UNITVECTOR(3)
REAL            :: invL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
invL=SQRT(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
IF(ABS(invL).GT.0.0)THEN
  invL=1./invL
  UNITVECTOR=v1*invL
ELSE
  UNITVECTOR = (/ 0. , 0. , 0./)
END IF ! ABS(invL).GT.0.0
END FUNCTION UNITVECTOR


PURE FUNCTION VECNORM(v1)
!===================================================================================================================================
! Computes the Euclidean norm (length) of a vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    ! Vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: VECNORM  ! Euclidean norm (length) of the vector v1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
VECNORM=SQRT(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
END FUNCTION VECNORM


PURE SUBROUTINE OrthoNormVec(v1,v2,v3)
!===================================================================================================================================
!> computes orthonormal basis from a given vector v1 (v1 must be normalized)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars, ONLY:EpsMach
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: v2(3), v3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(v1(3)).LT.100*EpsMach)THEN
  v2=(/-v1(2)-v1(3) , v1(1) , v1(1)       /)
ELSE
  v2=(/ v1(3)       , v1(3) ,-v1(1)-v1(2) /)
END IF
v2=UNITVECTOR(v2)
v3(:)=CROSSNORM(v1,v2)
END SUBROUTINE OrthoNormVec


PURE FUNCTION DOTPRODUCT(v1)
!===================================================================================================================================
! Computes the dot product of a vector with itself
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)       ! Input 3D vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: DOTPRODUCT  ! Dot product of v1 with itself
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DOTPRODUCT=v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
END FUNCTION DOTPRODUCT


PURE SUBROUTINE SphericalCoordinates(X,r,theta,phi)
!===================================================================================================================================
!> Computes the spherical coordinates of a Cartesian Vector X
!> r     : radial distance (Euclidean norm (length) of vector X)
!> theta : polar angle (angle between the positive Z-axis and the vector X) 
!>         0 <= theta <= Pi
!> phi   : azimuthal angle (angle between the projection of the vector onto the X-Y-plane and the positive X-axis) 
!>         0 <= phi < 2*Pi
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars, ONLY:Pi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: X(1:3) ! Vector in Cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: r,theta,phi ! radial distance, polar angle and azimuthal angle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Radial distance
r = VECNORM(X(1:3))

! Azimuthal angle
phi = ATAN2(X(2),X(1))
! If the angle comes out negative (this requires a negative Y value), add 2*Pi
IF(phi.LT.0.0) phi=phi+2*Pi

! Polar angle
IF(ABS(r).GT.0.0)THEN
  theta = ACOS(X(3)/r)
ELSE
  theta = 0.
END IF ! ABS(r).GT.0.0

END SUBROUTINE SphericalCoordinates


PURE SUBROUTINE TransformVectorfieldSphericalCoordinates(P,XHat,X)
!===================================================================================================================================
!> Transform a vector field component from spherical coordinates to Cartesian coordinates
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: P(1:3)    ! Position vector
REAL,INTENT(IN)  :: XHat(1:3) ! Vector in Spherical coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: X(1:3)    ! Resulting vector in Cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: r,theta,phi ! radial distance, polar angle and azimuthal angle
!===================================================================================================================================
! Get spherical coordinates
CALL SphericalCoordinates(P,r,theta,phi)

! Apply transformation matrix to vector in spherical coordinates to obtain the vector in Cartesian coordinates
CALL TransformVectorFromSphericalCoordinates(XHat,theta,phi,X)

END SUBROUTINE TransformVectorfieldSphericalCoordinates


PURE SUBROUTINE TransformVectorFromSphericalCoordinates(XHat,theta,phi,X)
!===================================================================================================================================
!> Transform a vector from spherical coordinates to Cartesian coordinates via supplied vector in spherical coordinates XHat,
!> azimuthal angle phi and polar angle theta
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: phi       !> azimuthal angle
REAL,INTENT(IN)  :: theta     !> polar angle
REAL,INTENT(IN)  :: XHat(1:3) !> Vector in Spherical coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: X(1:3)    ! Resulting vector in Cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Transformation matrix x vector in spherical coordinates
! X = M * XHat
X(1:3) = (/SIN(theta)*COS(phi)*XHat(1) + COS(theta)*COS(phi)*XHat(2) - SIN(phi)*XHat(3) ,&
           SIN(theta)*SIN(phi)*XHat(1) + COS(theta)*SIN(phi)*XHat(2) + COS(phi)*XHat(3) ,&
           COS(theta)*         XHat(1) - SIN(theta)*         XHat(2)                   /)

END SUBROUTINE TransformVectorFromSphericalCoordinates


END MODULE MOD_Globals
