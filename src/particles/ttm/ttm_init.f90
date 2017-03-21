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
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_TTM_Vars
#ifdef MPI
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif /* MPI*/
USE MOD_Particle_Mesh_Vars,      ONLY:ElemBaryNGeo
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
INTEGER        :: i,j,k,ix,iy,iz
INTEGER        :: io_error
INTEGER        :: iLineTTM
INTEGER        :: iElem,iElemFD
CHARACTER(255) :: StrTmp
REAL           :: tmp_array(1:15)
REAL           :: fd_hx,fd_hy,fd_hz
!===================================================================================================================================
IF(TTMInitIsDone)THEN
   SWRITE(*,*) "InitTTM already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TTM ...'


TTMFile   =TRIM(GETSTR('TTMFile',''))
TTMLogFile=TRIM(GETSTR('TTMLogFile',''))
IF((TRIM(TTMFile).NE.'').AND.(TRIM(TTMLogFile).NE.''))THEN
  DoImportTTMFile=.TRUE.
  ! get FD grid array dimension in x, y and z direction
  TTMGridFDdim=GETINTARRAY('TTMGridFDdim',3,'0 , 0 , 0')
  ! get FD grid element bary center tolerance for checking distance to elem bary center distance of DG cells
  TTMElemBaryTolerance=GETREAL('TTMElemBaryTolerance','1e-6')
  ! get FD grid cell sizes
  CALL GetFDGridCellSize(fd_hx,fd_hy,fd_hz)
  IF(ANY(TTMGridFDdim(1:3).LE.0).OR.ANY((/fd_hx,fd_hy,fd_hz/).LE.0.0))THEN
    DoImportTTMFile=.FALSE.
  ELSE
    FD_nElems=TTMGridFDdim(1)*TTMGridFDdim(2)*TTMGridFDdim(3)
    ALLOCATE( TTM_FD(1:11,TTMGridFDdim(1),TTMGridFDdim(2),TTMGridFDdim(3)) )
    TTM_FD=0.0
    ALLOCATE( ElemBaryFD(3,FD_nElems) )
    ElemBaryFD=0.0
    ALLOCATE( ElemIndexFD(3,FD_nElems) )
    ElemIndexFD=0.0
    ALLOCATE( ElemIsDone(FD_nElems) )
    ElemIsDone=.FALSE.
    ! IMD TTM data format and target arrays for data input (from ASCII data file *.ttm)
    !        1      2    3       4  5      6       7       8       9    10   11       ->  TTM_FD(1:11,ix,iy,iz,iLineTTM) -> TTM_DG
    !  1 2 3 4      5    6       7  8      9       10      11      12   13   14 15    ->  tmp_array(1:15)
    ! #x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z fd_k fd_g Z  proc
    print*,"Reading from file: ",TRIM(TTMFile)
#ifdef MPI
    IF(.NOT.PartMPI%MPIROOT)THEN
      !print*,"NOT root: myrank",myrank
      CALL abort(&
      __STAMP__&
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
    ! get index of *.ttm data file, e.g., "5" (corresponds to specific time stamp calculated from timestep)
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

    iLineTTM=0
    DO i=1,1 ! header lines: currently 1 -> adjust?
    !print*,"i",i
      !READ(ioUnit,*)
      READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
      !print*,i,' : ',TRIM(StrTmp) ! print header line
    END DO
    DO
      READ(ioUnit,*,IOSTAT=io_error) tmp_array(1:15)
      IF(io_error>0)THEN
        SWRITE(*,*) "io_error=",io_error
        SWRITE(*,*) "A problem reading the TTM file occured in line: ",iLineTTM+2
        print*,'[',tmp_array(1:15),']'
        CALL abort(&
        __STAMP__&
        ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
      ELSE IF(io_error<0)THEN
        SWRITE(*,*) "End of file reached: ",iLineTTM+1
        EXIT
      ELSE ! io_error = 0
      !print*,'[',tmp_array(1:15),']'
      END IF
      iLineTTM=iLineTTM+1 ! *.ttm data file: line counter

      ix=INT(tmp_array(1),4)+1 ! FD grid index in x
      iy=INT(tmp_array(2),4)+1 ! FD grid index in y
      iz=INT(tmp_array(3),4)+1 ! FD grid index in z

      IF( (ix.LE.TTMGridFDdim(1)) .AND. (iy.LE.TTMGridFDdim(2)) .AND. (iz.LE.TTMGridFDdim(3)) ) THEN
        TTM_FD(1:11,ix,iy,iz) = tmp_array(4:14)
        ElemIndexFD(1,iLineTTM)=ix
        ElemIndexFD(2,iLineTTM)=iy
        ElemIndexFD(3,iLineTTM)=iz
!print*,"TTM_FD(1:11,ix,iy,iz)=",TTM_FD(1:11,ix,iy,iz)
!read*
        ElemBaryFD(1,iLineTTM)=(REAL(ix)-0.5)*fd_hx 
        ElemBaryFD(2,iLineTTM)=(REAL(iy)-0.5)*fd_hy 
        ElemBaryFD(3,iLineTTM)=(REAL(iz)-0.5)*fd_hz 
          !print*,ElemBaryFD(:,iLineTTM)
!read*
      ELSE
        SWRITE(*,*) "A problem reading the TTM file occured in line: ",iLineTTM+1
        DoImportTTMFile=.FALSE.
        EXIT
      END IF
    END DO
    IF(FD_nElems.NE.ix*iy*iz)THEN
      CALL abort(&
      __STAMP__&
      ,'Number of FD elements FD_nElems does not comply with the number of FD elements read from *.ttm file'//TRIM(TTMFile))
    END IF
    CLOSE(ioUnit)
    print*,"Read ",iLineTTM," lines from file (+1 header line)"

    IF(FD_nElems.EQ.PP_nElems)THEN ! same number of elements in FD grid and DG solution -> assume they coincide
      ! the local DG solution in physical and reference space
      ALLOCATE( TTM_DG(1:11,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
      TTM_DG=0.0

      DO iElem=1,PP_nElems
        DO iELemFD=1,FD_nElems
          !print*,iElem,ElemBaryNGeo(:,iElem)
          !print*,iElemFD,ElemBaryFD(:,iElemFD)
!read*
          IF(ElemIsDone(iElemFD))CYCLE
          
          IF(AlmostEqualToTolerance(ElemBaryNGeo(1,iElem),ElemBaryFD(1,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
          IF(AlmostEqualToTolerance(ElemBaryNGeo(2,iElem),ElemBaryFD(2,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
          IF(AlmostEqualToTolerance(ElemBaryNGeo(3,iElem),ElemBaryFD(3,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
          ElemIsDone(iElemFD)=.TRUE.

          ix=ElemIndexFD(1,iElemFD)
          iy=ElemIndexFD(2,iElemFD)
          iz=ElemIndexFD(3,iElemFD)

!print*,"TTM_FD(1:11,ix,iy,iz)=",TTM_FD(1:11,ix,iy,iz)
!read*

          DO i=0,PP_N ! set all DG DOF values equal to FD cell value
            DO j=0,PP_N ! set all DG DOF values equal to FD cell value
              DO k=0,PP_N ! set all DG DOF values equal to FD cell value
                TTM_DG(1:11,i,j,k,iElem)=TTM_FD(1:11,ix,iy,iz)
              END DO
            END DO
          END DO

          EXIT

        END DO
      END DO
      IF(ANY(ElemIsDone).EQV..FALSE.)THEN
print*,"NOT all found :("
      ELSE
print*,"ALL FOUND!!!! :-)"
      END IF
      SDEALLOCATE(TTM_FD)
      SDEALLOCATE(ElemBaryFD)
      SDEALLOCATE(ElemIndexFD)
      SDEALLOCATE(ElemIsDone)
  
    ELSE ! use swap-mesh or interpolate the FD grid to the DG solution
      SWRITE(*,*) "FD_nElems=",FD_nElems
      SWRITE(*,*) "PP_nElems=",PP_nElems
      CALL abort(&
      __STAMP__&
      ,'Error interpolating ttm data (FD grid) to DG solution. FD_nElems.NE.PP_nElems. Different elements not implemented')
     
    END IF

  END IF
ELSE
  DoImportTTMFile=.FALSE.
END IF


TTMInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TTM DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitTTM


SUBROUTINE GetFDGridCellSize(fd_hx,fd_hy,fd_hz)
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize TTM variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_TTM_Vars
USE MOD_Particle_Vars, ONLY:IMDLengthScale
!USE MOD_RegressionCheck_tools, ONLY: str2real
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: fd_hx,fd_hy,fd_hz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: ioUnit
INTEGER        :: IndNum,IndNum2
INTEGER        :: i
INTEGER        :: io_error
CHARACTER(255) :: StrTmp,StrTmp2
!===================================================================================================================================
! search for the following line and extract "fd_hx,fd_hy,fd_hz"
! fd_h.x:24.107143, fd_h.y:22.509134,fd_h.z:22.509134\nnVolume of one FD cell: 1.221415e+04 [cubic Angstroms] 
fd_hx = -1.0
fd_hy = -1.0
fd_hz = -1.0

OPEN(NEWUNIT=ioUnit,FILE=TRIM(TTMLogFile),STATUS='OLD',ACTION='READ',IOSTAT=io_error)
IF(io_error.NE.0)THEN
  CALL abort(&
  __STAMP__&
  ,'Cannot open specified File (ttm data) from '//TRIM(TTMLogFile))
END IF
print*,"Reading from file: ",TRIM(TTMLogFile)
i=1
DO !i=1,1 ! header lines: currently 1 -> adjust?
  READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
  IF(io_error>0)THEN
    SWRITE(*,*) "io_error=",io_error
    SWRITE(*,*) "A problem reading the TTM file occured in line: ",i+1
    SWRITE(*,*) StrTmp
    !print*,'[',tmp_array(1:15),']'
    CALL abort(&
    __STAMP__&
    ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
  ELSE IF(io_error<0)THEN
    SWRITE(*,*) "End of file reached: ",i
    EXIT
  ELSE ! io_error = 0
    i=i+1
    IndNum=INDEX(TRIM(StrTmp),'fd_h.x:')
    IF(IndNum.GT.0)THEN
      StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
      IndNum=INDEX(TRIM(StrTmp),',')
      IF(IndNum.GT.0)THEN
        CALL str2real(StrTmp(8:IndNum-1),fd_hx,io_error)
        IF(io_error.NE.0)THEN
          SWRITE(*,*) "io_error=",io_error,", failed to get fd_hx"
          fd_hx=-99.
          Exit
        END IF
        fd_hx=fd_hx*IMDLengthScale
print*,"fd_hx=",fd_hx
        StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
      END IF


      IndNum=INDEX(TRIM(StrTmp),'fd_h.y:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
        IndNum=INDEX(TRIM(StrTmp),',')
        IF(IndNum.GT.0)THEN
          CALL str2real(StrTmp(8:IndNum-1),fd_hy,io_error)
          IF(io_error.NE.0)THEN
            SWRITE(*,*) "io_error=",io_error,", failed to get fd_hy"
            fd_hy=-99.
            Exit
          END IF
          fd_hy=fd_hy*IMDLengthScale
print*,"fd_hy=",fd_hy
          StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
        END IF
      END IF


      IndNum=INDEX(TRIM(StrTmp),'fd_h.z:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
  
          IndNum=INDEX(TRIM(StrTmp),'\')
          IF(IndNum.GT.0)THEN
            CALL str2real(StrTmp(8:IndNum-1),fd_hz,io_error)
            IF(io_error.NE.0)THEN
              SWRITE(*,*) "io_error=",io_error,", failed to get fd_hz"
              fd_hz=-99.
              Exit
            END IF
            fd_hz=fd_hz*IMDLengthScale
print*,"fd_hz=",fd_hz
            StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
        END IF
      END IF

      EXIT
    END IF
    IF(i.GT.100)EXIT ! only read the first 100 lines
  END IF
END DO



CLOSE(ioUnit)
END SUBROUTINE GetFDGridCellSize


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
SDEALLOCATE(TTM_FD)
SDEALLOCATE(TTM_DG)
END SUBROUTINE FinalizeTTM


END MODULE MOD_TTMInit
