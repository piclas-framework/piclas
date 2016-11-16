#include "boltzplatz.h"

!==================================================================================================================================
!> Contains the routines to 
!> - perform the actual regressioncheck
!==================================================================================================================================
MODULE MOD_RegressionCheck_Run
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------


INTERFACE PerformRegressionCheck
  MODULE PROCEDURE PerformRegressionCheck
END INTERFACE

PUBLIC::PerformRegressionCheck
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which performs the actual regressioncheck. It triggers the builds and execute commands. Additionally, it performs
!> the checks for L2-error norms, h5-diff and runtime
!==================================================================================================================================
SUBROUTINE PerformRegressionCheck()
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Build,   ONLY: ReadConfiguration_boltzplatz,BuildConfiguration_boltzplatz,GetFlagFromFile
USE MOD_RegressionCheck_Compare, ONLY: CompareNorm,CompareDataSet,CompareRuntime,ReadNorm,IntegrateLine
USE MOD_RegressionCheck_Tools,   ONLY: CheckForExecutable,InitExample,CleanExample
USE MOD_RegressionCheck_Vars,    ONLY: nExamples,ExampleNames,Examples,EXECPATH,RuntimeOptionType
USE MOD_RegressionCheck_Vars,    ONLY: BuildTESTCASE,BuildTIMEDISCMETHOD,BuildDir,BuildSolver
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(len=255)             :: cwd                               ! current cworking directory CALL getcwd(cwd)
CHARACTER(LEN=500)             :: SYSCOMMAND,dummystr               !> string to fit the system command
CHARACTER(LEN=255)             :: TESTCASE                          !> compilation flag in FLEXI
CHARACTER(LEN=255)             :: TIMEDISCMETHOD                    !> compilation flag in PICLas
CHARACTER(LEN=255)             :: ReggieBuildExe                    !> cache boltzplatz executables when doing "BuildSolver":
                                                                    !> this means don't build the same cmake configuration twice
CHARACTER(LEN=255)             :: FileName                          !> path to a file or its name
CHARACTER(LEN=255)             :: FolderName                        !> path to a folder or its name
CHARACTER(LEN=255)             :: parameter_boltzplatz              !> input parameter file for boltzplatz depending on EQNSYS
CHARACTER(LEN=255)             :: parameter_boltzplatz2             !> input parameter file for boltzplatz depending on TIMEDISC
LOGICAL                        :: ExistFile                         !> file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS                           !> status
INTEGER                        :: iExample                          !> loop index for example
REAL,ALLOCATABLE               :: ReferenceNorm(:,:)                !> L2 and Linf norm of the executed example from a reference
                                                                    !> solution
CHARACTER(LEN=255),ALLOCATABLE :: BuildConfigurations(:,:)          !> ??
LOGICAL,ALLOCATABLE            :: BuildValid(:)                     !> ??
INTEGER,ALLOCATABLE            :: BuildCounter(:)                   !> ??
INTEGER,ALLOCATABLE            :: BuildIndex(:)                     !> ??
INTEGER                        :: ErrorStatus                       !> Error-code of regressioncheck
INTEGER                        :: N_compile_flags                   !> number of compile-flags
INTEGER                        :: iReggieBuild,nReggieBuilds ! field handler unit and ??
INTEGER                        :: IndNum
INTEGER                        :: iSubExample
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Performing tests ...'
ReggieBuildExe=''
!==================================================================================================================================
DO iExample = 1, nExamples ! loop level 1 of 3
!==================================================================================================================================
  ! currently: each parameter configuration is only built and tested for the "run_particle" example
  dummystr=TRIM(ADJUSTL(ExampleNames(iExample)))
  IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).EQ.RuntimeOptionType)THEN
    SWRITE(UNIT_stdOut,*) ''
  END IF
  SWRITE(UNIT_stdOut,'(A,2x,A50)',ADVANCE='no') ' Example-Name: ',  TRIM(ExampleNames(iExample))
  IF(dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType)))).NE.RuntimeOptionType)THEN
    SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...skipping'
    CYCLE
  ELSE
    SWRITE(UNIT_stdOut,'(A,2x,A)') '  ...running'
  END IF
!print*,"dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType))))",dummystr(1:LEN(TRIM(ADJUSTL(RuntimeOptionType))))
!print*,"RuntimeOptionType",RuntimeOptionType
!read*
  ! if "BuildSolver" is true, boltzplatz's complete valid compiler-flag parameter 
  ! combination is tested (specified in "configuration.boltzplatz", default example is "run_particle")
  IF(BuildSolver)THEN
    IF(ReggieBuildExe.EQ.'')THEN
      CALL ReadConfiguration_boltzplatz(iExample,nReggieBuilds,BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)
      BuildCounter=1 ! reset the counter between read and build (used for selecting the build configuration for compilation)
      SYSCOMMAND='rm -rf '//TRIM(BuildDir)//'build_reggie_bin > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      SYSCOMMAND= 'mkdir '//TRIM(BuildDir)//'build_reggie_bin'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
    ELSE
      !print*,ReggieBuildExe
      print*,"CALL ReadConfiguration_boltzplatz() has already been performed, skipping..."
    END IF
  ELSE ! pure run, no compilation
    nReggieBuilds=1
  END IF
!==================================================================================================================================
  DO iReggieBuild = 1, nReggieBuilds ! loop level 2 of 3: cycle the number of build configurations (no configuration = only 1 run)
!==================================================================================================================================
    IF(BuildSolver)THEN
      ! if build is not valid no boltzplatz binary has been built and the lopp can cycle here
      IF(.NOT.BuildValid(iReggieBuild))THEN ! invalid reggie build 
        WRITE (ReggieBuildExe, '(a, i4.4)') "invalid"
      ELSE
        WRITE (ReggieBuildExe, '(a, i4.4)') "boltzplatz", COUNT(BuildValid(1:iReggieBuild)) ! e.g. boltzplatz0001
      END IF
      ! check if build exists -> if it does, don't build a new executable with cmake
      FileName=TRIM(BuildDir)//'build_reggie_bin/'//ReggieBuildExe
      INQUIRE(File=FileName,EXIST=ExistFile)
      IF(ExistFile) THEN ! 1. build already exists (e.g. boltzplatz0001 located in ../build_reggie_bin/)
        EXECPATH=TRIM(FileName)
      ELSE ! 2. build does not exists -> create it
        CALL BuildConfiguration_boltzplatz(iReggieBuild,nReggieBuilds,&
                                      BuildCounter,BuildIndex,N_compile_flags,BuildConfigurations,BuildValid)
        IF(BuildValid(iReggieBuild))THEN ! only move boltzplatz exe if it has been created (only for valid builds)
          SYSCOMMAND='cd '//TRIM(BuildDir)//' && mv build_reggie/bin/boltzplatz build_reggie_bin/'//ReggieBuildExe
          CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
          EXECPATH=TRIM(BuildDir)//'build_reggie_bin/'//ReggieBuildExe
        END IF
      END IF
      ! if build is not valid no exe has been built and the lopp can cycle here
      IF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.NE.nReggieBuilds))THEN ! invalid reggie build but not last reggie build
        CYCLE
      ELSEIF((.NOT.BuildValid(iReggieBuild)).AND.(iReggieBuild.EQ.nReggieBuilds))THEN ! last reggie build -> exit 
                                                                                      ! ("cycle" would start infinite loop)
        EXIT
      END IF
      CALL CheckForExecutable(Mode=2) ! check if executable was created correctly
    END IF
  
    ! read the parameters for the current example (parameter_reggie.ini)
    CALL InitExample(TRIM(Examples(iExample)%PATH),LEN(TRIM(Examples(iExample)%PATH)),Examples(iExample))

    ! depending on the equation system -> get different Nvar 
    CALL GetNvar(iExample,iReggieBuild)

    ! check if executable is compiled with correct TESTCASE (e.g. for tylorgreenvortex)
    IF(BuildSolver)THEN
      TESTCASE=BuildTESTCASE(iReggieBuild)
      TIMEDISCMETHOD=BuildTIMEDISCMETHOD(iReggieBuild)
    ELSE
      FileName=EXECPATH(1:INDEX(EXECPATH,'/',BACK = .TRUE.))//'configuration.cmake'
      INQUIRE(File=FileName,EXIST=ExistFile)
      IF(ExistFile) THEN
        CALL  GetFlagFromFile(FileName,'BOLTZPLATZ_TESTCASE',TESTCASE)
        ! set default for, e.g., PICLas code (currently no testcases are implemented)
        IF((TRIM(TESTCASE).EQ.'flag does not exist'))TESTCASE='default'
        IF(TRIM(TESTCASE).EQ.'flag does not exist')CALL abort(&
          __STAMP__&
          ,'BOLTZPLATZ_TESTCASE flag not found in configuration.cmake!',999,999.)

        CALL GetFlagFromFile(FileName,'BOLTZPLATZ_TIMEDISCMETHOD',TIMEDISCMETHOD)
        IF(TRIM(TIMEDISCMETHOD).EQ.'flag does not exist')CALL abort(&
          __STAMP__&
          ,'BOLTZPLATZ_TIMEDISCMETHOD flag not found in configuration.cmake!',999,999.)
      ELSE
        SWRITE(UNIT_stdOut,'(A12,A)')     ' ERROR: ','no "configuration.cmake" found at the location of the boltzplatz binary.'
        SWRITE(UNIT_stdOut,'(A12,A)')  ' FileName: ', TRIM(FileName)
        SWRITE(UNIT_stdOut,'(A12,L)') ' ExistFile: ', ExistFile
        ERROR STOP '-1'
      END IF
    END IF

    ! check folder name and decide whether it can be executed with the current boltzplatz binary (e.g. testcases ...)
    !CALL getcwd(cwd)
    !print*,cwd
    !FolderName=cwd(INDEX(cwd,'/',BACK = .TRUE.):LEN(TRIM(cwd)))
    !print*,FolderName
    !IndNum=INDEX(FolderName, '/')
    !IF(IndNum.GT.0)FolderName=FolderName(2:LEN(FolderName))
   
    !print*,TRIM(ExampleNames(iExample))
    !print*,TRIM(TESTCASE)
    IndNum=INDEX(ExampleNames(iExample), 'TESTCASE') ! look for TESTCASE- flag in folder name of example
    IF(IndNum.GT.0)THEN ! folder name contains 'TESTCASE'
      FolderName=ExampleNames(iExample) ! e.g. run_TESTCASE-taylorgreenvortex/
      FolderName=FolderName(IndNum+9:LEN(FolderName)) ! e.g. taylorgreenvortex/
      IndNum=INDEX(FolderName, '_')                   ! e.g. taylorgreenvortex_full/
      IF(IndNum.GT.0)FolderName=FolderName(1:IndNum-1)! e.g. taylorgreenvortex
      IndNum=INDEX(FolderName, '/')                   ! e.g. taylorgreenvortex/
      IF(IndNum.GT.0)FolderName=FolderName(1:IndNum-1)! e.g. taylorgreenvortex
      IF(FolderName.NE.TESTCASE)THEN ! e.g. taylorgreenvortex .NE. phill
          ! TESTCASE folder and non-TESTCASE boltzplatz
          SWRITE(UNIT_stdOut,'(A,2x,A)') ' TESTCASE not found in boltzplatz binary ...skipping'
        CYCLE
      ELSE
        SWRITE(UNIT_stdOut,'(A,2x,A)') ' TESTCASE is correct ...running' ! TESTCASE folder and TESTCASE boltzplatz
      END IF
    ELSE ! folder name does not contain 'TESTCASE'
      FolderName='default' ! for non-TESTCASE setups, set the default settings
      IF(FolderName.NE.TESTCASE)THEN ! e.g. default .NE. phill
        ! non-TESTCASE folder, but TESTCASE boltzplatz
        SWRITE(UNIT_stdOut,'(A)') ' TESTCASE "default" not found in boltzplatz binary: TESTCASE=['//TRIM(TESTCASE)//'] ...skipping'
        CYCLE
      ELSE
        ! non-TESTCASE folder and non-TESTCASE boltzplatz
        SWRITE(UNIT_stdOut,'(A)') ' TESTCASE "default" is correct: TESTCASE=['//TRIM(TESTCASE)//'] ...running' 
      END IF
    END IF

    !Examples(iExample)%Nvar=1
    ! debug 
    print*,'EXECPATH:                ',TRIM(EXECPATH)
    print*,'EQNSYSNAME:              ',TRIM(Examples(iExample)%EQNSYSNAME)
    print*,'nVar:                    ',     Examples(iExample)%Nvar    
    print*,'PATH:                    ',TRIM(Examples(iExample)%PATH)
    print*,'EXEC:                    ',     Examples(iExample)%EXEC  
    print*,'Reference:               ',TRIM(Examples(iExample)%ReferenceFile)
    print*,'State:                   ',TRIM(Examples(iExample)%ReferenceStateFile)
    print*,'HDF5 dataset:            ',TRIM(Examples(iExample)%ReferenceDataSetName)
    print*,'Restart:                 ',TRIM(Examples(iExample)%RestartFileName)
    print*,'Example%SubExample:      ',TRIM(Examples(iExample)%SubExample)
    print*,'Example%SubExampleNumber:',     Examples(iExample)%SubExampleNumber
  !SWRITE(UNIT_stdOut,'(A,I2,A4,A)')"Example%SubExampleOption(",iSubExample,") = ",TRIM(Example%SubExampleOption(iSubExample))
    ! perform simulation
    parameter_boltzplatz2=''
    IF(TRIM(TIMEDISCMETHOD).EQ.'DSMC')THEN
      parameter_boltzplatz='parameter_boltzplatz_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'_DSMC.ini'
      parameter_boltzplatz2='parameter_DSMC.ini'
    ELSE
      parameter_boltzplatz='parameter_boltzplatz_'//TRIM(ADJUSTL(Examples(iExample)%EQNSYSNAME))//'.ini'
    END IF
    INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_boltzplatz),EXIST=ExistFile)
    IF(.NOT.ExistFile) THEN
      SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
      SWRITE(UNIT_stdOut,'(A,A)') ' parameter_boltzplatz:      ',TRIM(parameter_boltzplatz)
      ERROR STOP '-1'
    END IF
    IF(parameter_boltzplatz2.NE.'')THEN
      INQUIRE(File=TRIM(Examples(iExample)%PATH)//TRIM(parameter_boltzplatz2),EXIST=ExistFile)
      IF(.NOT.ExistFile) THEN
        SWRITE(UNIT_stdOut,'(A,A)') ' ERROR: no File found under ',TRIM(Examples(iExample)%PATH)
        SWRITE(UNIT_stdOut,'(A,A)') ' parameter_boltzplatz2:     ',TRIM(parameter_boltzplatz2)
        ERROR STOP '-1'
      END IF
    END IF
!==================================================================================================================================
    DO iSubExample = 1, MAX(1,Examples(iExample)%SubExampleNumber) ! loop level 3 of 3: SubExamples (e.g. different TimeDiscMethods)
!==================================================================================================================================
      IF(Examples(iExample)%SubExampleNumber.GT.0)THEN ! SubExample has been specified
        SWRITE(UNIT_stdOut,'(A)')" "
        SWRITE(UNIT_stdOut,'(A,I2,A,A)')" SubExampleOption(",iSubExample,")=",TRIM(Examples(iExample)%SubExampleOption(iSubExample))
        SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//& ! print the current SubExampleOption(iSubExample) to parameter_boltzplatz
   ' && sed -i -e "s/.*'//TRIM(Examples(iExample)%SubExample)//&
                  '=.*/'//TRIM(Examples(iExample)%SubExample)//&
                     '='//TRIM(Examples(iExample)%SubExampleOption(iSubExample))//&
                   '/" '//TRIM(parameter_boltzplatz)
        CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      END IF
    IF(Examples(iExample)%IntegrateLine)THEN ! delete pre-existing data files before running the code
      SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && rm '//TRIM(Examples(iExample)%IntegrateLineFile)//' > /dev/null 2>&1'
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS) ! delete, e.g., "TGVAnalysis.dat" or "Database.csv"
    END IF
      IF(Examples(iExample)%EXEC)THEN
        SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && mpirun -np 2 '//TRIM(EXECPATH)//' '//TRIM(parameter_boltzplatz)//' ' &
                    //TRIM(parameter_boltzplatz2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
      ELSE
        SYSCOMMAND='cd '//TRIM(Examples(iExample)%PATH)//' && '//TRIM(EXECPATH)//' '//TRIM(parameter_boltzplatz)//' ' &
                    //TRIM(parameter_boltzplatz2)//' '//TRIM(Examples(iExample)%RestartFileName)//' 1>std.out 2>err.out'
      END IF
      CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
      print*,' --------------------------------------', iSTATUS
      IF(iSTATUS.EQ.0)THEN ! Computation successfull
        SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' Successfull computation ...'
        CALL AddError('Successfull computation',iExample,iSubExample,ErrorStatus=0,ErrorCode=0)
      ELSE ! Computation failed
        SWRITE(UNIT_stdOut,'(A)')   ' Computation of example failed'
        SWRITE(UNIT_stdOut,'(A,A)') ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
        SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
        CALL AddError('Error while executing',iExample,iSubExample,ErrorStatus=1,ErrorCode=2)
        CYCLE ! when a computation fails no useful output is created for comparison, continue the subexample cycle
      END IF
      ! comparing results and writing error messages for the current case
      SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' Comparing results...'
      ! check error norms
      ALLOCATE(ReferenceNorm(Examples(iExample)%nVar,2))
      IF(Examples(iExample)%ReferenceFile.EQ.'')THEN
        ! constant value, should be zero no reference file given
        CALL CompareNorm(ErrorStatus,iExample)
      ELSE
        ! read in reference and compare to reference solution
        CALL ReadNorm(iExample,ReferenceNorm)
        CALL CompareNorm(ErrorStatus,iExample,ReferenceNorm)
      END IF
      DEALLOCATE(ReferenceNorm)
      IF(ErrorStatus.EQ.1)THEN
        SWRITE(UNIT_stdOut,'(A)')   ' Error-norm mismatched! Example failed! '
        SWRITE(UNIT_stdOut,'(A)')   ' For more information: '
        SWRITE(UNIT_stdOut,'(A,A)') ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
        SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
        CALL AddError('Mismatch of error norms',iExample,iSubExample,ErrorStatus=1,ErrorCode=3)
      END IF

      ! diff h5 file
      IF(Examples(iExample)%ReferenceStateFile.NE.'')THEN
        CALL CompareDataSet(iExample)
        IF(Examples(iExample)%ErrorStatus.EQ.3)THEN
          CALL AddError('Mismatch in HDF5-files. Datasets are unequal',iExample,iSubExample,ErrorStatus=3,ErrorCode=4)
          SWRITE(UNIT_stdOut,'(A)')  ' Mismatch in HDF5-files'
        END IF
      END IF

      ! Integrate over line
      IF(Examples(iExample)%IntegrateLine)THEN
        CALL IntegrateLine(ErrorStatus,iExample)
        IF(Examples(iExample)%ErrorStatus.EQ.5)THEN
        CALL AddError('Mismatch in LineIntegral',iExample,iSubExample,ErrorStatus=5,ErrorCode=5)
        END IF
      END IF

      IF(Examples(iExample)%ErrorStatus.EQ.0)THEN
        SWRITE(UNIT_stdOut,'(A)')  ' Example successfull! '
        CALL CleanExample(iExample)
      END IF
      !SWRITE(UNIT_stdOut,'(A)')  '-------------------------------------'
    END DO ! iSubExample = 1, MAX(1,SubExampleNumber) (for cases without specified SubExamples: SubExampleNumber=0)
  END DO ! iReggieBuild = 1, nReggieBuilds
END DO ! iExample=1,nExamples

END SUBROUTINE PerformRegressionCheck


!==================================================================================================================================
!> Add an Error entry to the list of pointers containing information regarding the compilation process, execution process, 
!> Error codes and example info
!==================================================================================================================================
SUBROUTINE AddError(Info,iExample,iSubExample,ErrorStatus,ErrorCode)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: ExampleNames,Examples,EXECPATH,firstError,aError
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(len=*),INTENT(IN) :: Info
INTEGER,INTENT(IN)          :: iExample,iSubExample,ErrorStatus,ErrorCode
!INTEGER         :: a
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!LOGICAL         :: ALMOSTEQUAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Examples(iExample)%ErrorStatus=ErrorStatus
IF(firstError%ErrorCode.EQ.-1)THEN ! first error pointer
  firstError%ErrorCode              =ErrorCode ! no error
  firstError%Example                =TRIM(ExampleNames(iExample))
  firstError%SubExampleOption       =TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  firstError%Info                   =TRIM(Info)
  firstError%Build                  =TRIM(EXECPATH(INDEX(EXECPATH,'/',BACK=.TRUE.)+1:LEN(EXECPATH)))
  ALLOCATE(aError)
  aError=>firstError
ELSE ! next error pointer
  ALLOCATE(aError%nextError)
  aError%nextError%ErrorCode        =ErrorCode ! no error
  aError%nextError%Example          =TRIM(ExampleNames(iExample))
  aError%nextError%SubExampleOption =TRIM(Examples(iExample)%SubExampleOption(iSubExample))
  aError%nextError%Info             =TRIM(Info)
  aError%nextError%Build            =TRIM(EXECPATH(INDEX(EXECPATH,'/',BACK=.TRUE.)+1:LEN(EXECPATH)))
  aError=>aError%nextError
END IF

END SUBROUTINE AddError


!==================================================================================================================================
!> depending on the equation system -> get different Nvar for the current example
!> currently supports: - navierstokes             ->    Examples(iExample)%Nvar=5
!>                     - linearscalaradvection    ->    Examples(iExample)%Nvar=1
!==================================================================================================================================
SUBROUTINE GetNvar(iExample,iReggieBuild)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples,EXECPATH,BuildEQNSYS,BuildSolver
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iReggieBuild
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: IndNum
CHARACTER(LEN=255)             :: FileName                          ! path to a file or its name
CHARACTER(LEN=255)             :: temp,temp2                        ! temp variables for read in of file lines
LOGICAL                        :: ExistFile                         ! file exists=.true., file does not exist=.false.
INTEGER                        :: iSTATUS                           ! status
INTEGER                        :: ioUnit                            ! IO channel
!===================================================================================================================================
ioUnit=GETFREEUNIT()
IndNum=INDEX(EXECPATH, '/')
IF(IndNum.GT.0)THEN
  IndNum=INDEX(EXECPATH,'/',BACK = .TRUE.) ! get path without boltzplatz binary
  FileName=EXECPATH(1:IndNum)//'configuration.cmake'
  INQUIRE(File=FileName,EXIST=ExistFile)
  IF(ExistFile) THEN
    OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS="OLD",IOSTAT=iSTATUS,ACTION='READ') 
    DO
      READ(ioUnit,'(A)',iostat=iSTATUS)temp
      IF(ADJUSTL(temp).EQ.'!') CYCLE
      IF(iSTATUS.EQ.-1)EXIT
      IF(LEN(trim(temp)).GT.1)THEN
        IndNum=INDEX(temp,'BOLTZPLATZ_EQNSYSNAME')
        IF(IndNum.GT.0)THEN
          temp2=TRIM(ADJUSTL(temp(IndNum+LEN('BOLTZPLATZ_EQNSYSNAME'):LEN(temp))))
          IndNum=INDEX(temp2, '"')
          IF(IndNum.GT.0)THEN
            temp2=temp2(IndNum+1:LEN(TRIM(temp2)))
            IF(INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum.GT.IndNum)THEN ! get boltzplatz exe path up to 2nd '"' in name
              IndNum=INDEX(temp2(IndNum+1:LEN(TRIM(temp2))), '"')+IndNum
            END IF
          END IF
          Examples(iExample)%EQNSYSNAME=temp2(1:IndNum-1)
          EXIT
        END IF
      END IF
    END DO
    CLOSE(ioUnit)
  ELSE ! could not find 'configuration.cmake' at location of boltzplatz binary
    IF(BuildSolver)THEN ! get EQNSYSNAME from cmake build configuration settings
      Examples(iExample)%EQNSYSNAME=BuildEQNSYS(iReggieBuild)
    ELSE ! stop 
      SWRITE(UNIT_stdOut,'(A12,A)')     ' ERROR: ','no "configuration.cmake" found at the location of the boltzplatz binary.'
      SWRITE(UNIT_stdOut,'(A12,A)')  ' FileName: ', TRIM(FileName)
      SWRITE(UNIT_stdOut,'(A12,L)') ' ExistFile: ', ExistFile
      ERROR STOP '-1'
    END IF
  END IF
END IF
SELECT CASE (TRIM(Examples(iExample)%EQNSYSNAME))
  CASE ('navierstokes')  
    Examples(iExample)%Nvar=5
  CASE ('linearscalaradvection')  
    Examples(iExample)%Nvar=1
  CASE ('maxwell')  
    Examples(iExample)%Nvar=8
  CASE DEFAULT
    Examples(iExample)%Nvar=-1
    SWRITE(UNIT_stdOut,'(A)')   ' ERROR: missing case select for this BOLTZPLATZ_EQNSYSNAME with appropriate Nvar. Fix it by'
    SWRITE(UNIT_stdOut,'(A)')   '        adding the correct line of code to ../regressioncheck/regressioncheck_run.f90'
    SWRITE(UNIT_stdOut,'(A,A)') '        Examples(iExample)%EQNSYSNAME=',Examples(iExample)%EQNSYSNAME
    ERROR STOP '77'
END SELECT
END SUBROUTINE GetNvar


END MODULE MOD_RegressionCheck_Run
