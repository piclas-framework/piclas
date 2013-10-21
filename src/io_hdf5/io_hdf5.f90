#include "boltzplatz.h"

MODULE MOD_io_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE HDF5
USE MOD_Globals,ONLY: iError
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTEGER(HID_T)           :: File_ID
INTEGER(HSIZE_T),POINTER :: HSize(:)
INTEGER                  :: nDims

INTERFACE OpenDataFile
  MODULE PROCEDURE OpenHDF5File
END INTERFACE

INTERFACE CloseDataFile
  MODULE PROCEDURE CloseHDF5File
END INTERFACE

!===================================================================================================================================

CONTAINS

SUBROUTINE OpenHDF5File(FileString,create,single,communicatorOpt)
!===================================================================================================================================
! Open HDF5 file and groups
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString
LOGICAL,INTENT(IN)            :: create 
LOGICAL,INTENT(IN),OPTIONAL   :: single
INTEGER,INTENT(IN),OPTIONAL   :: communicatorOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Plist_ID
INTEGER                        :: comm
#ifdef MPI
INTEGER                        :: info !for lustre file system
#endif /* MPI */
!===================================================================================================================================
LOGWRITE(*,'(A)')'  OPEN HDF5 FILE "',TRIM(FileString),'" ...'
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Setup file access property list with parallel I/O access (MPI) or with default property list.
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#ifdef MPI
IF(PRESENT(communicatorOpt))THEN
  comm=communicatorOpt
ELSE
  comm=MPI_COMM_WORLD
END IF
IF (PRESENT(single)) THEN
  IF (.NOT. single) THEN
    !normal case:
    !info=MPI_INFO_NULL
    
    !for lustre file system:
    ! Create info to be attached to HDF5 file 
    CALL MPI_info_create(info,iError)
    ! Disables ROMIO's data-sieving 
    CALL MPI_Info_set(info, "romio_ds_read", "disable",iError)
    CALL MPI_Info_set(info, "romio_ds_write", "disable",iError)
    ! Enable ROMIO's collective buffering 
    CALL MPI_Info_set(info, "romio_cb_read", "enable",iError)
    CALL MPI_Info_set(info, "romio_cb_write", "enable",iError)
    ! end lustre stuff
    
    CALL H5PSET_FAPL_MPIO_F(Plist_ID,comm, info, iError)
  END IF
ELSE
!normal case:
!info=MPI_INFO_NULL

!for lustre file system:
! Create info to be attached to HDF5 file 
CALL MPI_info_create(info,iError)
! Disables ROMIO's data-sieving 
CALL MPI_Info_set(info, "romio_ds_read", "disable",iError)
CALL MPI_Info_set(info, "romio_ds_write", "disable",iError)
! Enable ROMIO's collective buffering 
CALL MPI_Info_set(info, "romio_cb_read", "enable",iError)
CALL MPI_Info_set(info, "romio_cb_write", "enable",iError)
! end lustre stuff

CALL H5PSET_FAPL_MPIO_F(Plist_ID,comm, info, iError)
END IF
#endif /* MPI */

! Open the file collectively.
IF(create)THEN
 CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, access_prp = Plist_ID)
ELSE !read-only
 CALL H5FOPEN_F(TRIM(FileString), H5F_ACC_RDONLY_F, File_ID, iError, access_prp = Plist_ID)
 IF(iError.EQ.-1) THEN
    WRITE(*,*) "ERROR: Can't open file: ", TRIM(FileString)
    STOP
 END IF      
END IF
CALL H5PCLOSE_F(Plist_ID, iError)
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE OpenHDF5File



SUBROUTINE CloseHDF5File()
!===================================================================================================================================
! Close HDF5 file and groups
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,UNIT_logOut,Logging
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
LOGWRITE(*,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE CloseHDF5File

END MODULE MOD_io_HDF5
