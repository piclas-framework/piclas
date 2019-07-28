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
!==================================================================================================================================
!> Initializes HDF5 IO and sets HDF-MPI parameters, opens ans closes files.
!==================================================================================================================================
MODULE MOD_IO_HDF5
! MODULES
USE HDF5
USE MOD_Globals,ONLY: iError
IMPLICIT NONE

ABSTRACT INTERFACE
  SUBROUTINE EvalElemInt(ElemData)
  USE MOD_Mesh_Vars,ONLY:nElems
  REAL,INTENT(OUT) :: ElemData(nElems)
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: gatheredWrite       !< flag whether every process should output data or data should first be gathered
INTEGER(HID_T)           :: File_ID             !< file which is currently opened
INTEGER(HID_T)           :: Plist_File_ID       !< property list of file which is currently opened
INTEGER(HSIZE_T),POINTER :: HSize(:)            !< HDF5 array size (temporary variable)
INTEGER                  :: nDims               !<
INTEGER                  :: MPIInfo             !< hardware / storage specific / file system MPI parameters to pass to HDF5
                                                !< for optimized performance on specific systems

!> Type containing pointers to data to be written to HDF5 in an element-wise scalar fashion.
!> Alternatively a function pointer can be specified providing the desired data.
!> Only one of the pointers may be associated.
TYPE tElementOut
  CHARACTER(LEN=255)                    :: VarName                    !< variable name
  REAL,POINTER                          :: RealArray(:) => NULL()
  REAL,POINTER                          :: RealScalar   => NULL()
  INTEGER,POINTER                       :: IntArray(:)  => NULL()
  INTEGER(KIND=8),POINTER               :: LongIntArray(:) => NULL()
  INTEGER,POINTER                       :: IntScalar    => NULL()
  LOGICAL,POINTER                       :: LogArray(:)  => NULL()
  PROCEDURE(EvalElemInt),POINTER,NOPASS :: eval         => NULL()
  TYPE(tElementOut),POINTER             :: next         => NULL()     !< next list item
END TYPE

TYPE(tElementOut),POINTER    :: ElementOut   => NULL() !< linked list of output pointers

INTERFACE InitIO
  MODULE PROCEDURE InitIO_HDF5
END INTERFACE

INTERFACE OpenDataFile
  MODULE PROCEDURE OpenDataFile
END INTERFACE

INTERFACE CloseDataFile
  MODULE PROCEDURE CloseDataFile
END INTERFACE

INTERFACE AddToElemData
  MODULE PROCEDURE AddToElemData
END INTERFACE

INTERFACE ClearElemData
  MODULE PROCEDURE ClearElemData
END INTERFACE

PUBLIC::DefineParametersIO,InitIO,OpenDataFile,CloseDataFile
PUBLIC::AddToElemData

!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersIO()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("IO_HDF5")
CALL prms%CreateLogicalOption('gatheredWrite', "Set true to activate gathered HDF5 IO for parallel computations. "//&
                                               "Only local group masters will write data after gathering from local slaves.",&
                                               '.FALSE.')
END SUBROUTINE DefineParametersIO

SUBROUTINE InitIO_HDF5()
!===================================================================================================================================
! Initialize HDF5 IO
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,        ONLY:GETLOGICAL
#ifdef INTEL
USE IFPORT
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

gatheredWrite=.FALSE.
IF(nLeaderProcs.LT.nProcessors) gatheredWrite=GETLOGICAL('gatheredWrite','.FALSE.')
#ifdef MPI
  CALL MPI_Info_Create(MPIInfo, iError)

  !normal case:
  MPIInfo=MPI_INFO_NULL

  ! Large block IO extremely slow on Juqeen cluster (only available on IBM clusters)
  !CALL MPI_Info_set(MPIInfo, "IBM_largeblock_io", "true", ierror)
#ifdef LUSTRE
  CALL MPI_Info_Create(MPIInfo, iError)
  ! For lustre file system:
  ! Disables ROMIO's data-sieving
  CALL MPI_Info_set(MPIInfo, "romio_ds_read", "disable",iError)
  CALL MPI_Info_set(MPIInfo, "romio_ds_write","disable",iError)
  ! Enable ROMIO's collective buffering
  CALL MPI_Info_set(MPIInfo, "romio_cb_read", "enable", iError)
  CALL MPI_Info_set(MPIInfo, "romio_cb_write","enable", iError)
#endif
#endif /* MPI */
END SUBROUTINE InitIO_HDF5


!==================================================================================================================================
!> Open HDF5 file and groups
!==================================================================================================================================
SUBROUTINE OpenDataFile(FileString,create,single,readOnly,communicatorOpt,userblockSize)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString      !< filename to be opened
LOGICAL,INTENT(IN)           :: create          !< create file if it doesn't exist. Overwrited file if already present!
LOGICAL,INTENT(IN)           :: single          !< single=T : only one processor opens file, single=F : open/create collectively
LOGICAL,INTENT(IN)           :: readOnly        !< T : file is opened in read only mode, so file system timestamp remains unchanged
                                                !< F: file is open read/write mode
INTEGER,INTENT(IN),OPTIONAL  :: communicatorOpt !< only MPI and single=F: communicator to be used for collective access
INTEGER,INTENT(IN),OPTIONAL  :: userblockSize   !< size of the file to be prepended to HDF5 file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HSIZE_T)               :: userblockSize_loc, tmp, tmp2
!==================================================================================================================================
LOGWRITE(*,'(A)')'  OPEN HDF5 FILE "'//TRIM(FileString)//'" ...'

userblockSize_loc = 0
IF (PRESENT(userblockSize)) userblockSize_loc = userblockSize

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)

! Setup file access property list with parallel I/O access (MPI) or with default property list.
IF(create)THEN
  CALL H5PCREATE_F(H5P_FILE_CREATE_F, Plist_File_ID, iError)
  IF(iError.NE.0) CALL abort(__STAMP__,&
    'ERROR: Could not create file '//TRIM(FileString))
ELSE
  CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_File_ID, iError)
  IF(iError.NE.0) CALL abort(__STAMP__,&
    'ERROR: Could not open file '//TRIM(FileString))
END IF

#ifdef MPI
IF(.NOT.single)THEN
  IF(.NOT.PRESENT(communicatorOpt))CALL abort(__STAMP__,&
    'ERROR: communicatorOpt must be supplied in OpenDataFile when single=.FALSE.')
  CALL H5PSET_FAPL_MPIO_F(Plist_File_ID, communicatorOpt, MPIInfo, iError)
END IF
  IF(iError.NE.0) CALL abort(__STAMP__,&
    'ERROR: H5PSET_FAPL_MPIO_F failed in OpenDataFile')
#endif /* MPI */

! Open the file collectively.
IF(create)THEN
  IF (userblockSize_loc > 0) THEN
    tmp = userblockSize_loc/512
    IF (MOD(userblockSize_loc,512).GT.0) tmp = tmp+1
    tmp2 = 512*2**CEILING(LOG(REAL(tmp))/LOG(2.))
    CALL H5PSET_USERBLOCK_F(Plist_File_ID, tmp2, iError)
  END IF
  CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, creation_prp = Plist_File_ID)
ELSE !read-only ! and write (added later)
  IF(.NOT.FILEEXISTS(FileString)) CALL abort(__STAMP__,&
    'ERROR: Specified file '//TRIM(FileString)//' does not exist.')
  IF (readOnly) THEN
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDONLY_F,  File_ID, iError, access_prp = Plist_File_ID)
  ELSE
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDWR_F,  File_ID, iError, access_prp = Plist_File_ID)
  END IF
END IF
IF(iError.NE.0) CALL abort(__STAMP__,&
  'ERROR: Could not open or create file '//TRIM(FileString))

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE OpenDataFile



!==================================================================================================================================
!> Close HDF5 file and groups
!==================================================================================================================================
SUBROUTINE CloseDataFile()
! MODULES
USE MOD_Globals,ONLY:UNIT_logOut,Logging
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
LOGWRITE(*,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5PCLOSE_F(Plist_File_ID, iError)
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE CloseDataFile

!==================================================================================================================================
!> Set pointers to element-wise arrays or scalars which will be gathered and written out. Both real or integer data types
!> are supported. It is also possible to pass a function pointer which will be evaluated to calculate the data.
!==================================================================================================================================
SUBROUTINE AddToElemData(ElementOut_In,VarName,RealArray,RealScalar,IntArray,LongIntArray,IntScalar,LogArray,Eval)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElementOut),POINTER,INTENT(INOUT)    :: ElementOut_In        !< Pointer list of element-wise data that
                                                                   !< is written to the state file
CHARACTER(LEN=*),INTENT(IN)                :: VarName              !< Name of the current array/scalar
REAL,INTENT(IN),TARGET,OPTIONAL            :: RealArray(nElems)    !< Data is an array containing reals
REAL,INTENT(IN),TARGET,OPTIONAL            :: RealScalar           !< Data is a real scalar
INTEGER,INTENT(IN),TARGET,OPTIONAL         :: IntArray(nElems)     !< Data is an array containing integers
INTEGER,INTENT(IN),TARGET,OPTIONAL         :: IntScalar            !< Data is a integer scalar
INTEGER(KIND=8),INTENT(IN),TARGET,OPTIONAL :: LongIntArray(nElems) !< Data is a integer scalar
LOGICAL,INTENT(IN),TARGET,OPTIONAL         :: LogArray(nElems)     !< Data is a logical scalar
PROCEDURE(EvalElemInt),POINTER,OPTIONAL    :: Eval                 !< Data is evaluated using a function pointer
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElementOut),POINTER          :: eout
INTEGER                            :: nOpts
!==================================================================================================================================
IF(.NOT.ASSOCIATED(ElementOut_In))THEN
  ! list is empty, create first entry
  ALLOCATE(ElementOut_In)
  eout=>ElementOut_In
ELSE
  eout=>ElementOut_In
  ! loop until last entry
  DO WHILE(ASSOCIATED(eout%next))
    eout=>eout%next
  END DO
  ! insert new entry
  ALLOCATE(eout%next)
  eout=>eout%next
ENDIF

! set varname and data pointer
NULLIFY(eout%next)
eout%VarName=VarName
nOpts=0
IF(PRESENT(RealArray))THEN
  eout%RealArray  => RealArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(RealScalar))THEN
  eout%RealScalar => RealScalar
  nOpts=nOpts+1
ENDIF
IF(PRESENT(IntArray))THEN
  eout%IntArray   => IntArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(IntScalar))THEN
  eout%IntScalar  => IntScalar
  nOpts=nOpts+1
ENDIF
IF(PRESENT(LongIntArray))THEN
  eout%LongIntArray  => LongIntArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(LogArray))THEN
  eout%LogArray  => LogArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(eval))THEN
  eout%eval       => Eval
  nOpts=nOpts+1
ENDIF
IF(nOpts.NE.1) CALL Abort(__STAMP__,&
  'More then one optional argument passed to AddToElemData.')
END SUBROUTINE AddToElemData


!==================================================================================================================================
!> Deallocate all pointers to element-wise arrays or scalars which will be gathered and written out.
!> The linked list of the pointer is deallocated for each entry
!==================================================================================================================================
SUBROUTINE ClearElemData(ElementOut)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElementOut),POINTER,INTENT(INOUT)    :: ElementOut           !< Pointer list of element-wise data that
                                                                   !< is written to the state file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElementOut),POINTER          :: e,e2
!==================================================================================================================================
IF(.NOT.ASSOCIATED(ElementOut))THEN
  RETURN
ENDIF

e=>ElementOut
DO WHILE(ASSOCIATED(e))
  e%VarName = ''
  IF(ASSOCIATED(e%RealArray))    NULLIFY(e%RealArray    ) !=> NULL()
  IF(ASSOCIATED(e%RealScalar))   NULLIFY(e%RealScalar   ) !=> NULL()
  IF(ASSOCIATED(e%IntArray))     NULLIFY(e%IntArray     ) !=> NULL()
  IF(ASSOCIATED(e%IntScalar))    NULLIFY(e%IntScalar    ) !=> NULL()
  IF(ASSOCIATED(e%LongIntArray)) NULLIFY(e%LongIntArray ) !=> NULL()
  IF(ASSOCIATED(e%LogArray))     NULLIFY(e%LogArray     ) !=> NULL()
  IF(ASSOCIATED(e%eval))         NULLIFY(e%eval         ) !=> NULL()
  e2=>e%next
  ! deallocate stuff
  DEALLOCATE(e)
  e=> e2
END DO

!ElementOut   => NULL() !< linked list of output pointers
NULLIFY(ElementOut)

END SUBROUTINE ClearElemData

END MODULE MOD_io_HDF5
