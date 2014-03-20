#include "boltzplatz.h"

MODULE MOD_ReadInTools
!===================================================================================================================================
! MODULES
!===================================================================================================================================
USE MOD_Globals
USE MOD_ISO_VARYING_STRING
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! Public Part
PUBLIC::GETSTR
PUBLIC::CNTSTR
PUBLIC::GETINT
PUBLIC::GETREAL
PUBLIC::GETLOGICAL
PUBLIC::GETINTARRAY
PUBLIC::GETREALARRAY
PUBLIC::IgnoredStrings
!===================================================================================================================================

INTERFACE GETSTR
  MODULE PROCEDURE GETSTR
END INTERFACE

INTERFACE CNTSTR
  MODULE PROCEDURE CNTSTR
END INTERFACE

INTERFACE GETINT
  MODULE PROCEDURE GETINT
END INTERFACE

INTERFACE GETREAL
  MODULE PROCEDURE GETREAL
END INTERFACE

INTERFACE GETLOGICAL
  MODULE PROCEDURE GETLOGICAL
END INTERFACE

INTERFACE GETINTARRAY
  MODULE PROCEDURE GETINTARRAY
END INTERFACE

INTERFACE GETREALARRAY
  MODULE PROCEDURE GETREALARRAY
END INTERFACE

INTERFACE IgnoredStrings
  MODULE PROCEDURE IgnoredStrings
END INTERFACE

INTERFACE FillStrings
  MODULE PROCEDURE FillStrings
END INTERFACE

INTERFACE FindStr
  MODULE PROCEDURE FindStr
END INTERFACE

INTERFACE LowCase
  MODULE PROCEDURE LowCase
END INTERFACE

INTERFACE GetNewString
  MODULE PROCEDURE GetNewString
END INTERFACE

INTERFACE DeleteString
  MODULE PROCEDURE DeleteString
END INTERFACE

TYPE tString
  TYPE(Varying_String)::Str
  TYPE(tString),POINTER::NextStr,PrevStr
END TYPE tString

LOGICAL,PUBLIC::ReadInDone=.FALSE.
TYPE(tString),POINTER::FirstString

CONTAINS

FUNCTION GETSTR(Key,Proposal)
!===================================================================================================================================
! Read string named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)                   :: GetStr   ! String read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,GetStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,GetStr,DefMsg)
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(Key),' | ', TRIM(GetStr),' | ',TRIM(DefMsg),' | '
END FUNCTION GETSTR



FUNCTION CNTSTR(Key,Proposal)
!===================================================================================================================================
! Counts all occurances of string named "key" from inifile and store in "GETSTR". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Inifile was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: CntStr   ! Number of parameters named "Key" in inifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: IntProposal
CHARACTER(LEN=LEN(Key))              :: TmpKey
TYPE(tString),POINTER                :: Str1
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

CntStr=0
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)

! Search
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    CntStr=CntStr+1
  END IF ! (INDEX...
  ! Next string in list
  Str1=>Str1%NextStr
END DO
IF (CntStr.EQ.0) THEN
  IF (PRESENT(Proposal)) THEN
    READ(Proposal,'(I3)')IntProposal
    CntStr=IntProposal
  ELSE
    SWRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
    CALL abort(__STAMP__,'Code stopped!',999,999.)
  END IF
END IF
END FUNCTION CNTSTR



FUNCTION GETINT(Key,Proposal)
!===================================================================================================================================
! Read integer named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: GetInt  ! Integer read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetInt
SWRITE(UNIT_StdOut,'(a3,a30,a3,i33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetInt,' | ',TRIM(DefMsg),' | '
END FUNCTION GETINT



FUNCTION GETREAL(Key,Proposal)
!===================================================================================================================================
! Read real named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: GetReal  ! Real read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetReal
SWRITE(UNIT_StdOut,'(a3,a30,a3,e33.5,a3,a7,a3)')' | ',TRIM(Key),' | ', GetReal,' | ',TRIM(DefMsg),' | '
END FUNCTION GETREAL



FUNCTION GETLOGICAL(Key,Proposal)
!===================================================================================================================================
! Read logical named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key        ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal   ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: GetLogical ! Logical read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetLogical
SWRITE(UNIT_StdOut,'(a3,a30,a3,l33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetLogical,' | ',TRIM(DefMsg),' | '
END FUNCTION GETLOGICAL



FUNCTION GETINTARRAY(Key,nIntegers,Proposal)
!===================================================================================================================================
! Read array of "nIntegers" integer values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nIntegers        ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                   :: GetIntArray(nIntegers)      ! Integer array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)        :: HelpStr
CHARACTER(LEN=8)          :: DefMsg
INTEGER                   :: iInteger
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetIntArray
SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Integer array of size (',nIntegers,') | ',TRIM(DefMsg),' | '
DO iInteger=0,nIntegers-1
  IF ((iInteger.GT.0) .AND. (MOD(iInteger,8).EQ.0)) THEN
    SWRITE(UNIT_stdOut,*)
    SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  SWRITE(UNIT_stdOut,'(i5)',ADVANCE='NO')GetIntArray(iInteger+1)
END DO
SWRITE(UNIT_stdOut,*)
END FUNCTION GETINTARRAY



FUNCTION GETREALARRAY(Key,nReals,Proposal)
!===================================================================================================================================
! Read array of "nReals" real values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nReals           ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                      :: GetRealArray(nReals)        ! Real array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)        :: HelpStr
CHARACTER(LEN=8)          :: DefMsg
INTEGER                   :: iReal
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
CALL getPImultiplies(helpstr)
READ(HelpStr,*)GetRealArray
SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Real array of size (',nReals,') | ',TRIM(DefMsg),' | '
DO iReal=0,nReals-1
  IF ((iReal.GT.0) .AND. (MOD(iReal,8).EQ.0)) THEN
    SWRITE(UNIT_stdOut,*)
    SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  SWRITE(UNIT_stdOut,'(f5.2)',ADVANCE='NO')GetRealArray(iReal+1)
END DO
SWRITE(UNIT_stdOut,*)
END FUNCTION GETREALARRAY



SUBROUTINE IgnoredStrings()
!===================================================================================================================================
! Prints out remaining strings in list after read-in is complete
!===================================================================================================================================
! MODULES
USE MOD_ISO_VARYING_STRING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tString),POINTER                  :: Str1
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')" THE FOLLOWING INI-FILE PARAMETERS WERE IGNORED:"
Str1=>FirstString
DO WHILE(ASSOCIATED(Str1))
  SWRITE(UNIT_stdOut,'(A4,A)')" |- ",TRIM(CHAR(Str1%Str))
  Str1=>Str1%NextStr
END DO
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE IgnoredStrings



SUBROUTINE FillStrings(IniFile)
!===================================================================================================================================
! Read ini file and put each line in a string object. All string objects are connected to a list of string objects starting
! with "firstString"
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ISO_VARYING_STRING
USE MOD_DSMC_Vars,ONLY: UseDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL   :: IniFile                    ! Name of ini file to be read in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tString),POINTER                       :: Str1=>NULL(),Str2=>NULL()
CHARACTER(LEN=255)                          :: HelpStr,Str
CHARACTER(LEN=300)                          :: File, DSMCFile
TYPE(Varying_String)                        :: aStr,bStr,Separator
INTEGER                                     :: EOF
INTEGER                                     :: stat,iniUnit,nLines,i
CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE :: FileContent,FileContent2
CHARACTER(LEN=1)                            :: tmpChar=''
!===================================================================================================================================
! Check if we have read in ini file already
IF (ReadInDone) RETURN
! Get name of ini file
IF (PRESENT(IniFile)) THEN
  File = TRIM(IniFile)
ELSE
  CALL GETARG(1,File)
  !CALL GET_COMMAND_ARGUMENT(1,File)
END IF
SWRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(File),'":'

IF(MPIRoot)THEN
  iniUnit=GETFREEUNIT()
  OPEN(UNIT   = iniUnit,       &
       FILE   = File,          &
       STATUS = 'OLD',         &
       ACTION = 'READ',        &
       ACCESS = 'SEQUENTIAL',  &
       IOSTAT = stat)
  IF(stat.NE.0) THEN
    CALL abort(__STAMP__,"Could not open ini file.")
  ELSE
    nLines=0
    stat=0
    DO
      READ(iniunit,"(A)",IOSTAT=stat)tmpChar
      IF(stat.NE.0)EXIT
      nLines=nLines+1
    END DO
  END IF
END IF
#ifdef MPI
CALL MPI_BCAST(nLines,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#endif
ALLOCATE(FileContent(nLines))
IF (MPIRoot) THEN
  !read file
  REWIND(iniUnit)
  READ(iniUnit,'(A)')FileContent
  CLOSE(iniUnit)
END IF
#ifdef MPI
CALL MPI_BCAST(FileContent,LEN(FileContent)*nLines,MPI_CHARACTER,0,MPI_COMM_WORLD,iError)
#endif

NULLIFY(Str1,Str2)
DO i=1,nLines
  IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
  ! Read line from memory
  aStr=var_str(FileContent(i))
  Str=aStr
  ! Remove comments with "!"
  CALL Split(aStr,Str1%Str,"!")
  ! Remove comments with "#"
  CALL Split(Str1%Str,bStr,"#")
  Str1%Str=bStr
  ! Remove "%" sign from old ini files, i.e. mesh% disc% etc.
  CALL Split(Str1%Str,bStr,"%",Separator,Back=.false.)

  ! If we have a newtype ini file, take the other part
  IF(LEN(CHAR(Separator)).EQ.0) Str1%Str=bStr
  ! Remove blanks
  Str1%Str=Replace(Str1%Str," ","",Every=.true.)
  ! Replace brackets
  Str1%Str=Replace(Str1%Str,"(/"," ",Every=.true.)
  Str1%Str=Replace(Str1%Str,"/)"," ",Every=.true.)
  ! Replace commas
  Str1%Str=Replace(Str1%Str,","," ",Every=.true.)
  ! Lower case
  CALL LowCase(CHAR(Str1%Str),HelpStr)
  ! If we have a remainder (no comment only)
  IF(LEN_TRIM(HelpStr).GT.2) THEN
    Str1%Str=Var_Str(HelpStr)
    IF(.NOT.ASSOCIATED(Str2)) THEN
      FirstString=>Str1
    ELSE
      Str2%NextStr=>Str1
      Str1%PrevStr=>Str2
    END IF
    Str2=>Str1
    CALL GetNewString(Str1)
  END IF
END DO

IF (useDSMC) THEN
  IF(MPIRoot) THEN  
    CALL GETARG(2,DSMCFile)
    SWRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(DSMCFile),'":'
    iniUnit=GETFREEUNIT()
    OPEN(UNIT   = iniUnit,    &
         FILE   = DSMCFile,       &
         STATUS = 'OLD',      &
         ACTION = 'READ',     &
         ACCESS = 'SEQUENTIAL',&
         IOSTAT = stat)
    IF(stat.NE.0) THEN
      CALL abort(__STAMP__,"Could not open ini file.")
    ELSE
      nLines=0
      stat=0
      DO
        READ(iniunit,"(A)",IOSTAT=stat)tmpChar
        IF(stat.NE.0)EXIT
        nLines=nLines+1
      END DO
    END IF
  END IF
#ifdef MPI
  CALL MPI_BCAST(nLines,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#endif
  ALLOCATE(FileContent2(nLines))
  IF (MPIRoot) THEN
    !read file
    REWIND(iniUnit)
    READ(iniUnit,'(A)')FileContent2
    CLOSE(iniUnit)
  END IF
#ifdef MPI
  CALL MPI_BCAST(FileContent2,LEN(FileContent2)*nLines,MPI_CHARACTER,0,MPI_COMM_WORLD,iError)
#endif
  NULLIFY(Str1,Str2)
  DO i=1,nLines
    IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
    ! Read line from memory
    aStr=var_str(FileContent2(i))
    Str=aStr
    ! Remove comments with "!"
    CALL Split(aStr,Str1%Str,"!")
    ! Remove comments with "#"
    CALL Split(Str1%Str,bStr,"#")
    Str1%Str=bStr
    ! Remove "%" sign from old ini files, i.e. mesh% disc% etc.
    CALL Split(Str1%Str,bStr,"%",Separator,Back=.false.)
  
    ! If we have a newtype ini file, take the other part
    IF(LEN(CHAR(Separator)).EQ.0) Str1%Str=bStr
    ! Remove blanks
    Str1%Str=Replace(Str1%Str," ","",Every=.true.)
    ! Replace brackets
    Str1%Str=Replace(Str1%Str,"(/"," ",Every=.true.)
    Str1%Str=Replace(Str1%Str,"/)"," ",Every=.true.)
    ! Replace commas
    Str1%Str=Replace(Str1%Str,","," ",Every=.true.)
    ! Lower case
    CALL LowCase(CHAR(Str1%Str),HelpStr)
    ! If we have a remainder (no comment only)
    IF(LEN_TRIM(HelpStr).GT.2) THEN
      Str1%Str=Var_Str(HelpStr)
      IF(.NOT.ASSOCIATED(Str2)) THEN
        FirstString=>Str1
      ELSE
        Str2%NextStr=>Str1
        Str1%PrevStr=>Str2
      END IF
      Str2=>Str1
      CALL GetNewString(Str1)
    END IF
  END DO
END IF

END SUBROUTINE FillStrings



SUBROUTINE GetNewString(Str)
!===================================================================================================================================
! Create and initialize new string object.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tString),POINTER :: Str ! New string
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
NULLIFY(Str)
ALLOCATE(Str)
NULLIFY(Str%NextStr,Str%PrevStr)
END SUBROUTINE GetNewString



SUBROUTINE DeleteString(Str)
!===================================================================================================================================
! Remove string "Str" from list of strings witFirstString,h first element "DirstString" and delete string.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    TYPE(tString),POINTER :: Str         ! String to delete
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF (ASSOCIATED(Str%NextStr)) Str%NextStr%PrevStr=>Str%PrevStr
IF (ASSOCIATED(Str,FirstString)) THEN
  FirstString=>Str%NextStr
ELSE
  Str%PrevStr%NextStr=>Str%NextStr
END IF
DEALLOCATE(Str)
NULLIFY(Str)
END SUBROUTINE DeleteString



SUBROUTINE FindStr(Key,Str,DefMsg,Proposal)
!===================================================================================================================================
! Find parameter string containing keyword "Key" in list of strings starting with "FirstString" and return string "Str" without
! keyword. If keyword is not found in list of strings, return default values "Proposal" (error if not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key         ! Search for this keyword in ini file
CHARACTER(LEN=8),INTENT(INOUT)       :: DefMsg      ! Default message = keyword not found, return default parameters (if available)
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal    ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT)         :: Str         ! Parameter string without keyword
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=LEN(Key))              :: TmpKey 
TYPE(tString),POINTER                :: Str1
LOGICAL                              :: Found
!===================================================================================================================================
DefMsg='*CUSTOM'
! Convert to lower case
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)
Found=.FALSE.
Str1=>FirstString
DO WHILE(.NOT.Found)
  IF (.NOT.ASSOCIATED(Str1)) THEN
    IF (.NOT.PRESENT(Proposal)) THEN
      CALL abort(__STAMP__,  &
           'Inifile missing necessary keyword item : '//TRIM(TmpKey)//' - Code stopped!',999,999.)
    ELSE ! Return default value
      CALL LowCase(TRIM(Proposal),Str)
      IF (Str(1:1).NE.'@') THEN
        DefMsg='DEFAULT'
      END IF
      RETURN
    END IF ! (.NOT.PRESENT(Proposal))
  END IF ! (.NOT.ASSOCIATED(Str1))

  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    Found=.TRUE.
    Str1%Str=replace(Str1%Str,TRIM(TmpKey)//'=',"",Every=.TRUE.)
    Str=TRIM(CHAR(Str1%Str))
    ! Remove string from list
    CALL DeleteString(Str1)
  ELSE
    ! Next string in list
    Str1=>Str1%NextStr
  END IF

END DO
END SUBROUTINE FindStr



SUBROUTINE LowCase(Str1,Str2)
!===================================================================================================================================
! Transform upper case letters in "Str1" into lower case letters, result is "Str2"
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: Str1 ! Input string 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT) :: Str2 ! Output string, lower case letters only
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iLen,nLen,Upper
CHARACTER(LEN=*),PARAMETER   :: lc='abcdefghijklmnopqrstuvwxyz'
CHARACTER(LEN=*),PARAMETER   :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LOGICAL                      :: HasEq
!===================================================================================================================================
HasEq=.FALSE.
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
END DO
END SUBROUTINE LowCase

SUBROUTINE getPImultiplies(helpstr)
!===================================================================================================================================
! Searches for the occurence of 'PI','Pi','pi' and 'pI' in a helpstr and replaces
! it with the value of pi=3.1415... etc. and oes a multiplication.
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT OUTPUT VARIABLE
    CHARACTER(LEN=*) :: helpstr   ! Input character string
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
    type(varying_string)      :: separator
    type(varying_string)      :: astr,bstr,cstr,dstr
    CHARACTER(LEN=500)        :: dummystr
    REAL                      :: PI
    REAL                      :: adummy
    LOGICAL                   :: finished
!===================================================================================================================================
  !Initialiazation
  dstr=var_str("")
  PI=ACOS(-1.)
  finished=.false.
  !Replace all occurences of pi in the string by one symbol
  helpstr=replace(helpstr,"PI","@",every=.true.)  ! Insert value of Pi
  helpstr=replace(helpstr,"pi","@",every=.true.)  ! Insert value of Pi
  helpstr=replace(helpstr,"pI","@",every=.true.)  ! Insert value of Pi
  helpstr=replace(helpstr,"Pi","@",every=.true.)  ! Insert value of Pi
  astr=var_str(helpstr)
  ! loop over string
  DO WHILE(.NOT. finished)
    !split sting at "@-occurences"
    CALL split(astr,bstr,"@",separator,back=.false.) !bStr is string in front of @
    IF(len(char(separator)) .NE. 0)THEN 
      ! we have found something, bnow get the factor in front of @
      CALL split(bstr,cstr," ",separator,back=.true.)
      IF(LEN(char(cstr)) .EQ. 0)THEN
        !no factor
        adummy=1
      ELSE
        !extract factor 
        dummystr=trim(char(cstr))
        READ(dummystr,*)adummy
      ENDIF
      !do the multiplication and recombine the string into "dstr"
      adummy=PI*adummy
      WRITE(dummystr,'(2a,F15.12)')trim(char(dstr)),trim(char(bstr)),adummy
      dstr=var_str(dummystr)
    ELSE
      ! we did not find anything now recombine the remaining string into "dstr"
      WRITE(dummystr,'(2a)')trim(char(dstr)),trim(char(bstr))
      dstr=var_str(dummystr)
      finished=.true.
    END IF
  END DO
  helpstr=trim(char(dstr))

END SUBROUTINE getPImultiplies

END MODULE MOD_ReadInTools
