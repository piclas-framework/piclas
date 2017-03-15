#include "boltzplatz.h"

MODULE MOD_TTMInit
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitTTM
  MODULE PROCEDURE InitTTM
END INTERFACE

INTERFACE FinalizeTTM
  MODULE PROCEDURE FinalizeTTM
END INTERFACE

PUBLIC::InitTTM,FinalizeTTM
!===================================================================================================================================

CONTAINS

SUBROUTINE InitTTM()
!===================================================================================================================================
! IMD Two Temperature Model (TTM)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_TTM_Vars
#ifdef MPI
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif /* MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: ioUnit
INTEGER        :: IndNum
INTEGER        :: i,ix,iy,iz
INTEGER        :: io_error
INTEGER        :: TTMFileN
CHARACTER(255) :: StrTmp
REAL           :: tmp_array(1:15)
!===================================================================================================================================
IF(TTMInitIsDone)THEN
   SWRITE(*,*) "InitTTM already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TTM ...'


TTMFile=TRIM(GETSTR('TTMFile',''))
IF(TRIM(TTMFile).NE.'')THEN
  DoImportTTMFile=.TRUE.







  TTMGridFDdim=GETINTARRAY('TTMGridFDdim',3,'0 , 0 , 0')
  IF(ANY(TTMGridFDdim(1:3).LE.0))THEN
    DoImportTTMFile=.FALSE.
  ELSE
  !        1      2    3       4  5      6       7       8       9    10   11 12
  ! #x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z fd_k fd_g Z proc
  ALLOCATE( TTM_FD(1:11,TTMGridFDdim(1),TTMGridFDdim(2),TTMGridFDdim(3)) )


  print*,"Reading from file: ",TRIM(TTMFile)
#ifdef MPI
  IF(.NOT.PartMPI%MPIROOT)THEN
    !print*,"NOT root: myrank",myrank
    CALL abort(__STAMP__&
    ,'ERROR: Cannot SetParticlePosition in multi-core environment for SpaceIC=IMD!')
  ELSE
    !print*,"I am root: myrank",myrank
  END IF
#endif /*MPI*/
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(TTMFile),STATUS='OLD',ACTION='READ',IOSTAT=io_error)
  IF(io_error.NE.0)THEN
    CALL abort(&
    __STAMP__&
    ,'Cannot open specified File (ttm data) from '//TRIM(TTMFile))
  END IF
  IndNum=INDEX(TTMFile, '/',BACK = .TRUE.)
  IF(IndNum.GT.0)THEN
    StrTmp=TRIM(TTMFile(IndNum+1:LEN(TTMFile)))
    IndNum=INDEX(StrTmp,'.',BACK = .TRUE.)
    IF(IndNum.GT.0)THEN
      StrTmp=StrTmp(1:IndNum-1)
      IndNum=INDEX(StrTmp,'.')
      IF(IndNum.GT.0)THEN
        StrTmp=StrTmp(IndNum+1:LEN(StrTmp))
      END IF
    END IF
  END IF

  !read(Species(FractNbr)%Init(iInit)%IMDFile(7:11),*,iostat=io_error)  i
  read(StrTmp,*,iostat=io_error)  i
  TTMNumber = i
  print*,"IMD *.ttm file = ",StrTmp
  print*,"TTMNumber",TTMNumber

  TTMFileN=0
  DO i=1,1 ! header lines
  !print*,"i",i
    !READ(ioUnit,*)
    READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
    print*,i,' : ',TRIM(StrTmp)
  END DO
  DO !i=1,10!chunkSize
    READ(ioUnit,*,IOSTAT=io_error) tmp_array(1:15)
    IF(io_error>0)THEN
      SWRITE(*,*) "io_error=",io_error
      SWRITE(*,*) "A problem reading the TTM file occured in line: ",TTMFileN+2
      print*,'[',tmp_array(1:15),']'
      CALL abort(__STAMP__&
      ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
    ELSE IF(io_error<0)THEN
      SWRITE(*,*) "End of file reached: ",TTMFileN+1
      EXIT
    ELSE
    !print*,'[',tmp_array(1:15),']'
    END IF
    TTMFileN=TTMFileN+1

    !print*,'[',tmp_array(1:15),']'
    !read*

    ix=INT(tmp_array(1),4)+1
    iy=INT(tmp_array(2),4)+1
    iz=INT(tmp_array(3),4)+1

    IF( (ix.LE.TTMGridFDdim(1)) .AND. (iy.LE.TTMGridFDdim(2)) .AND. (iz.LE.TTMGridFDdim(3)) ) THEN
      !print*,ix,iy,iz
      TTM_FD(1:11,ix,iy,iz) = tmp_array(4:14)

    ELSE
      SWRITE(*,*) "A problem reading the TTM file occured in line: ",TTMFileN+1
      DoImportTTMFile=.FALSE.
      EXIT
    END IF





  END DO
  FD_Nelems=ix*iy*iz
  CLOSE(ioUnit)
  print*,"Read ",TTMFileN," lines from file (+1 header line)"






  ! the local DG solution in physical and reference space
  !ALLOCATE( TTM(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))








  END IF
ELSE
  DoImportTTMFile=.FALSE.
END IF


TTMInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TTM DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
stop
END SUBROUTINE InitTTM


SUBROUTINE FinalizeTTM() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize TTM variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_TTM_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SDEALLOCATE()
END SUBROUTINE FinalizeTTM


END MODULE MOD_TTMInit
