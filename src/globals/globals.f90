#include "boltzplatz.h"

MODULE MOD_Globals
!===================================================================================================================================
! Provides parameters, used globally (please use EXTREMLY carefully!) 
!===================================================================================================================================
! MODULES
#ifdef MPI
USE mpi
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER ::UNIT_stdOut=6
INTEGER,PARAMETER ::UNIT_logOut=133
INTEGER           ::UNIT_errOut=999
LOGICAL           ::Logging
CHARACTER(LEN=255)::ErrorFileName='NOT_SET'
INTEGER           ::iError
REAL              ::StartTime
INTEGER           ::myRank
INTEGER           ::nProcessors
LOGICAL           ::MPIRoot
#ifdef MPI
!#include "mpif.h"
INTEGER           :: MPIStatus(MPI_STATUS_SIZE)
#endif

INTERFACE InitGlobals
  MODULE PROCEDURE InitGlobals
END INTERFACE

INTERFACE AlmostEqual
  MODULE PROCEDURE AlmostEqual
END INTERFACE

INTERFACE AlmostZero
  MODULE PROCEDURE AlmostZero
END INTERFACE

INTERFACE Abort
  MODULE PROCEDURE AbortProg
END INTERFACE Abort

INTERFACE INTSTAMP
  MODULE PROCEDURE INTSTAMP
END INTERFACE INTSTAMP

INTERFACE TIMESTAMP
  MODULE PROCEDURE TIMESTAMP
END INTERFACE

INTERFACE BOLTZPLATZTIME
  MODULE PROCEDURE BOLTZPLATZTIME
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
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS ...'

Pi=ACOS(-1.)
spi = 1./pi

! get machine accuracy
epsMach=EPSILON(0.0)
TwoEpsMach=2.0d0*epsMach

SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitGlobals


FUNCTION AlmostEqual(Num1,Num2)
!===================================================================================================================================
! Performe an almost equal check
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,    ONLY:TwoEpsMach
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL            :: Num1,Num2      ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: ALMOSTEQUAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ALMOSTEQUAL=.FALSE.
IF(ABS(Num1-Num2).LE.MAX(ABS(Num1),ABS(Num2))*TwoEpsMach) ALMOSTEQUAL=.TRUE.

END FUNCTION AlmostEqual


FUNCTION AlmostZero(Num)
!===================================================================================================================================
! Performe an almost zero check
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,    ONLY:EpsMach
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL            :: Num ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: AlmostZero
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

AlmostZero=.FALSE.
IF(ABS(Num).LE.EpsMach) AlmostZero=.TRUE.

END FUNCTION AlmostZero


SUBROUTINE AbortProg(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfoOpt,RealInfoOpt,SingleOpt)
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
LOGICAL,OPTIONAL                  :: SingleOpt       ! Only MPI-Root performs check
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   There is no way back!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: IntInfo         ! Error info (integer)
REAL                              :: RealInfo        ! Error info (real)
#ifdef MPI
INTEGER                           :: errOut          ! Output of MPI_ABORT
#endif /*MPI*/
!===================================================================================================================================
IF(PRESENT(SingleOpt))THEN
  IF(SingleOpt.AND.(.NOT.MPIRoot)) RETURN
END IF
IntInfo  = MERGE(IntInfoOpt ,999 ,PRESENT(IntInfoOpt) )
RealInfo = MERGE(RealInfoOpt,999.,PRESENT(RealInfoOpt))
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
#ifdef MPI
CALL MPI_ABORT(MPI_COMM_WORLD,iError,errOut)
#endif
STOP 0001
END SUBROUTINE AbortProg


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
USE MOD_Analyze_Vars,          ONLY: Analyze_dt
! IMPLICIT VARIABLE HANDLING
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
IF (Analyze_dt.LT.1E-10) THEN
  WRITE(TimeStamp,'(F15.14)')Time
ELSE
  WRITE(TimeStamp,'(F15.12)')Time
END IF
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP


#ifdef MPI
FUNCTION BOLTZPLATZTIME(Comm)
#else
FUNCTION BOLTZPLATZTIME()
#endif
!===================================================================================================================================
! Calculates current time (own function because of a laterMPI implementation)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef MPI
INTEGER, INTENT(IN),OPTIONAL    :: Comm
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: BoltzplatzTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
#ifdef MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
BoltzplatzTime=MPI_WTIME()
#else
CALL CPU_TIME(BoltzplatzTime)
#endif
END FUNCTION BOLTZPLATZTIME

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
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    ! 
REAL,INTENT(IN) :: v2(3)    ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: CROSS(3) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS

FUNCTION CROSSNORM(v1,v2)
!===================================================================================================================================
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
! and normalizes the vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    ! 
REAL,INTENT(IN) :: v2(3)    ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: CROSSNORM(3) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL            :: length
!===================================================================================================================================
CROSSNORM=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
length=SQRT(CROSSNORM(1)*CROSSNORM(1)+CROSSNORM(2)*CROSSNORM(2)+CROSSNORM(3)*CROSSNORM(3))
CROSSNORM=CROSSNORM/length
END FUNCTION CROSSNORM

END MODULE MOD_Globals
