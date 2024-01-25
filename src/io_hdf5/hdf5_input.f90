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

MODULE MOD_HDF5_Input
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ISVALIDHDF5FILE
  MODULE PROCEDURE ISVALIDHDF5FILE
END INTERFACE

INTERFACE ISVALIDMESHFILE
  MODULE PROCEDURE ISVALIDMESHFILE
END INTERFACE

INTERFACE GetHDF5NextFileName
  MODULE PROCEDURE GetHDF5NextFileName
END INTERFACE

INTERFACE DatasetExists
  MODULE PROCEDURE DatasetExists
END INTERFACE

INTERFACE GetDataSize
  MODULE PROCEDURE GetDataSize
END INTERFACE

INTERFACE GetAttributeSize
  MODULE PROCEDURE GetAttributeSize
END INTERFACE

INTERFACE GetDataProps
  MODULE PROCEDURE GetDataProps
END INTERFACE

! no interface because one argument is size of another argument-array -> nVal, Array(nVal)
!INTERFACE ReadArray
!  MODULE PROCEDURE ReadArray
!END INTERFACE

INTERFACE ReadAttribute
  MODULE PROCEDURE ReadAttribute
END INTERFACE

INTERFACE GetVarnames
  MODULE PROCEDURE GetVarnames
END INTERFACE

PUBLIC :: ISVALIDHDF5FILE,ISVALIDMESHFILE,GetDataProps,GetAttributeSize,GetHDF5NextFileName
PUBLIC :: ReadArray,ReadAttribute
PUBLIC :: File_ID,HSize,nDims        ! Variables that need to be public
PUBLIC :: OpenDataFile,CloseDataFile ! Subroutines that need to be public
PUBLIC :: DatasetExists
PUBLIC :: GetDataSize
PUBLIC :: GetVarnames
PUBLIC :: GetArrayAndName
!===================================================================================================================================

CONTAINS

FUNCTION ISVALIDHDF5FILE(FileName,FileVersionRealOpt,FileVersionIntOpt)
!===================================================================================================================================
! Subroutine to check if a file is a valid PICLas HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName
REAL,INTENT(IN),OPTIONAL       :: FileVersionRealOpt
INTEGER,INTENT(IN),OPTIONAL    :: FileVersionIntOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                        :: isValidHDF5File
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: FileVersionReal,FileVersionRealRef
INTEGER                        :: FileVersionInt,FileVersionIntRef
INTEGER(HID_T)                 :: Plist_ID
CHARACTER(LEN=255)             :: ProgramName
LOGICAL                        :: help
!===================================================================================================================================
IF(.NOT.FILEEXISTS(FileName))THEN
  SWRITE(UNIT_stdOut,'(A)')' ERROR: The file does not exit! FileName: '//TRIM(FileName)
  isValidHDF5File=.FALSE.
  RETURN
END IF ! .NOT.FILEEXISTS(FileName)
isValidHDF5File=.TRUE.
iError=0
FileVersionRealRef=-1.0
IF(PRESENT(FileVersionRealOpt)) FileVersionRealRef=FileVersionRealOpt
FileVersionIntRef=-1
IF(PRESENT(FileVersionRealOpt)) FileVersionIntRef=FileVersionIntOpt

! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
IF(iError.NE.0) CALL Abort(__STAMP__,'ERROR: COULD NOT OPEN FILE!')

! Open HDF5 file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
CALL H5PCLOSE_F(Plist_ID, iError)
IF(iError.EQ.0) THEN
  isValidHDF5File=.TRUE.
  ! Check program name -------------------------------------------------------------------------------------------------------------
  ! Open the attribute "Program" of root group
  CALL ReadAttribute(File_ID,'Program',1,StrScalar=ProgramName)
  help=.FALSE.
  IF(TRIM(ProgramName) .EQ. 'PICLas') help=.TRUE.
  IF(TRIM(ProgramName) .EQ. 'Boltzplatz') help=.TRUE.
  IF(TRIM(ProgramName) .EQ. 'Flexi') help=.TRUE.
  IF(.NOT.help) isValidHDF5File=.FALSE.

  ! Check file version -------------------------------------------------------------------------------------------------------------
  ! Open the attribute "File_Version" of root group
  IF(FileVersionRealRef.GT.0.0)THEN
    CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionReal)
    !IF(FileVersionReal .LT. FileVersionRealRef)THEN
    !  isValidHDF5File=.FALSE.
    !  SWRITE(UNIT_stdOut,'(A)')' ERROR: FILE VERSION TOO OLD! FileName: '//TRIM(FileName)
    !END IF
  END IF ! FileVersionRealRef.GT.0.0
  IF(FileVersionIntRef.GT.0)THEN
    CALL ReadAttribute(File_ID,'Piclas_VersionInt',1,IntScalar=FileVersionInt)
    !IF(FileVersionReal .LT. FileVersionIntRef)THEN
    !  isValidHDF5File=.FALSE.
    !  SWRITE(UNIT_stdOut,'(A)')' ERROR: FILE VERSION TOO OLD! FileName: '//TRIM(FileName)
    !END IF
  END IF ! FileVersionIntRef.GT.0.0
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  ! Close FORTRAN predefined datatypes
!!!  isValidHDF5File=.FALSE.
  CALL H5CLOSE_F(iError)
END IF
END FUNCTION ISVALIDHDF5FILE


!==================================================================================================================================
!> Subroutine to check if a file is a valid mesh file
!==================================================================================================================================
FUNCTION ISVALIDMESHFILE(MeshFileName)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName    !< name of mesh file to be checked
LOGICAL                        :: isValidMeshFile !< result: file is valid mesh file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: NGeoExists
INTEGER(HID_T)                 :: Plist_ID
!==================================================================================================================================
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Create property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if USE_MPI
! Setup file access property list with parallel I/O access (MPI)
CALL H5PSET_FAPL_MPIO_F(Plist_ID,MPI_COMM_PICLAS, MPIInfo, iError)
#endif /*USE_MPI*/

! Check if file exists
IF(.NOT.FILEEXISTS(MeshFileName)) THEN
  CALL abort(__STAMP__,'ERROR: Mesh file '//TRIM(MeshFileName)//' does not exist.')
  isValidMeshFile = .FALSE.
  RETURN
END IF

! Open HDF5 file
CALL H5FOPEN_F(TRIM(MeshFileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
IF(iError.EQ.0) THEN
  isValidMeshFile=.TRUE.

  ! Check NGeo attribute --------------------------------------------------------------------------------------------------------
  CALL DatasetExists(File_ID,'Ngeo',NGeoExists,attrib=.TRUE.)
  IF (.NOT.NGeoExists) isValidMeshFile = .FALSE.

  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  isValidMeshFile=.FALSE.
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
END IF

END FUNCTION ISVALIDMESHFILE

!==================================================================================================================================
!> Subroutine to determine HDF5 data size
!==================================================================================================================================
SUBROUTINE GetDataSize(Loc_ID,DSetName,nDims,IntSize,attrib)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName   !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID     !< ID of datase
LOGICAL,INTENT(IN),OPTIONAL          :: attrib     !< logical wether atrtibute or dataset
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)                  :: nDims      !< found data size dimensions
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: IntSize(:) !< found data size
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
LOGICAL                              :: attrib_loc
!===================================================================================================================================
IF (PRESENT(attrib)) THEN
  attrib_loc=attrib
ELSE
  attrib_loc=.FALSE.
END IF
IF(attrib_loc)THEN
  ! Open the dataset with default properties.
  CALL H5AOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
  ! Get the data space of the dataset.
  CALL H5AGET_SPACE_F(DSet_ID, FileSpace, iError)
  ! Get number of dimensions of data space
  CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
  ! Get size and max size of data space
  ALLOCATE(IntSize(nDims),SizeMax(nDims))
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, IntSize, SizeMax, iError)
  CALL H5SCLOSE_F(FileSpace, iError)
  CALL H5ACLOSE_F(DSet_ID, iError)
ELSE
  ! Open the dataset with default properties.
  CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
  ! Get the data space of the dataset.
  CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
  ! Get number of dimensions of data space
  CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
  ! Get size and max size of data space
  ALLOCATE(IntSize(nDims),SizeMax(nDims))
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, IntSize, SizeMax, iError)
  CALL H5SCLOSE_F(FileSpace, iError)
  CALL H5DCLOSE_F(DSet_ID, iError)
END IF
DEALLOCATE(SizeMax)
END SUBROUTINE GetDataSize


!==================================================================================================================================
!> Subroutine to determine HDF5 size of attribute
!==================================================================================================================================
SUBROUTINE GetAttributeSize(Loc_ID,AttribName,nDims,Size)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: AttribName !< name if attribute to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
INTEGER,INTENT(OUT)                  :: nDims    !< found data size dimensions
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  !< found data size
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: Attr_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
!==================================================================================================================================
! Open the dataset with default properties.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
! Get the data space of the dataset.
CALL H5AGET_SPACE_F(Attr_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5ACLOSE_F(Attr_ID, iError)
DEALLOCATE(SizeMax)
END SUBROUTINE GetAttributeSize


!==================================================================================================================================
!> @brief Subroutine to check whether a dataset in the HDF5 file exists
!>
!> We have no "h5dexists_f", so we use the error given by h5dopen_f.
!> this produces HDF5 error messages even if everything is okay, so we turn the error msgs off
!> during this operation.
!> auto error messages off
!==================================================================================================================================
SUBROUTINE DatasetExists(Loc_ID_in,DSetName,Exists,attrib,DSetName_attrib)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName        !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID_in       !< ID of dataset
LOGICAL,INTENT(IN),OPTIONAL          :: attrib          !< check dataset or attribute
CHARACTER(LEN=*),OPTIONAL            :: DSetName_attrib !< name if dataset to be checked
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: Exists          !< result: dataset exists
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: Loc_ID
LOGICAL                              :: attrib_loc
!==================================================================================================================================
Loc_ID=Loc_ID_in
IF (PRESENT(attrib)) THEN
  attrib_loc=attrib
  IF(PRESENT(DSetName_attrib))THEN
    ! Open dataset
    IF(TRIM(DSetName_attrib).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DSetName_attrib),Loc_ID, iError)
  END IF
ELSE
  attrib_loc=.FALSE.
END IF

! Check attribute or data set. Data sets can be checked by determining the existence of the corresponding link
IF(attrib_loc)THEN
  CALL H5AEXISTS_F(Loc_ID, TRIM(DSetName), Exists, iError)
ELSE
  CALL H5LEXISTS_F(Loc_ID, TRIM(DSetName), Exists, iError)
END IF

END SUBROUTINE DatasetExists


!==================================================================================================================================
!> Subroutine to determine HDF5 dataset properties
!==================================================================================================================================
SUBROUTINE GetDataProps(DatasetName,nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,nDimsOffset_opt)
! MODULES
USE MOD_Globals
USE MOD_ReadInTools        ,ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)             :: DatasetName   !< Name of Dataset that should be read
INTEGER,INTENT(OUT)                     :: nVar_HDF5     !< number of variables
INTEGER,INTENT(OUT)                     :: N_HDF5        !< polynomial degree
INTEGER,INTENT(OUT)                     :: nElems_HDF5   !< inumber of elements
CHARACTER(LEN=255),OPTIONAL,INTENT(OUT) :: NodeType_HDF5 !< nodetype string
INTEGER,INTENT(IN),OPTIONAL             :: nDimsOffset_opt !< optional shift when reading the values in each dimension of the array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: Rank,nDimsOffset_loc
INTEGER(HID_T)                          :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
!==================================================================================================================================
LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A,A)')' GET SIZE OF DATA IN HDF5 FILE... '

! Dimensional shift (optional) if arrays with rank > 5 are processed (e.g. DG_Solution from state files with an additional
! dimension that corresponds to time)
IF (PRESENT(nDimsOffset_opt)) THEN; nDimsOffset_loc = nDimsOffset_opt
ELSE                              ; nDimsOffset_loc = 0
END IF

! Read in attributes
! Open given dataset with default properties.
CALL H5DOPEN_F(File_ID, TRIM(DatasetName), Dset_ID, iError)

! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
CALL PrintOption('Rank of database','HDF5',IntOpt=Rank) ! 'HDF5.'
! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)
IF(PRESENT(NodeType_HDF5)) THEN
  ! Read in NodeType
  CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_HDF5)
END IF

! Display data
! nVar = first array index
nVar_HDF5 = INT(Dims(1),4)
CALL PrintOption('Number of variables nVar','HDF5',IntOpt=nVar_HDF5) ! 'HDF5.'
! N = index 2-4 of array, is expected to have the same value for each direction
IF (Rank.EQ.2) THEN
  N_HDF5 = 1
ELSE
  ! U(1:nVar,0:N,0:N,0:N,nElems) is the assumed array size,
  ! but possibly shifted if the array has more than 5 dimensions (-nDimsOffset_loc)
  N_HDF5 = INT(Dims(Rank-1-nDimsOffset_loc)) - 1
END IF
CALL PrintOption('Polynomial degree N','HDF5',IntOpt=N_HDF5) ! 'HDF5.'
IF(PRESENT(NodeType_HDF5)) THEN
  CALL PrintOption('Node type','HDF5',StrOpt=NodeType_HDF5) ! 'HDF5.'
END IF
! nElems = index Rank of array, possibly shifted if the array has more than 5 dimensions (-nDimsOffset_loc)
nElems_HDF5 = INT(Dims(Rank-nDimsOffset_loc),4)
CALL PrintOption('Number of Elements','HDF5',IntOpt=nElems_HDF5) ! 'HDF5.'

LBWRITE(UNIT_stdOut,'(A)')' DONE!'
LBWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE GetDataProps


SUBROUTINE GetVarnames(AttribName,VarNames,AttribExists)
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                :: AttribName
CHARACTER(LEN=255),ALLOCATABLE,INTENT(OUT) :: VarNames(:)
LOGICAL,INTENT(OUT)                        :: AttribExists
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dims, nVal
!===================================================================================================================================
SDEALLOCATE(VarNames)
CALL DatasetExists(File_ID,AttribName,AttribExists,attrib=.TRUE.)
IF (AttribExists) THEN
  ! get size of array
  CALL GetAttributeSize(File_ID,AttribName,dims,HSize)
  nVal=INT(HSize(1))
  DEALLOCATE(HSize)
  ALLOCATE(VarNames(nVal))

  ! read variable names
  CALL ReadAttribute(File_ID,TRIM(AttribName),nVal,StrArray=VarNames)
END IF
END SUBROUTINE GetVarnames

!===================================================================================================================================
!> High level wrapper to ReadArray and ReadAttrib. Check if array exists and directly
!> allocate, read array and attributes
!> Assume that the array to be read is of size (nVar,.,.,.,.,nElems) and that an associated
!> attribute containing the variable names exists
!===================================================================================================================================
SUBROUTINE GetArrayAndName(ArrayName,AttribName,nVal,Array,VarNames)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nElems,nGlobalElems,OffsetElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)     :: ArrayName   !< name of array to be read
CHARACTER(LEN=*),INTENT(IN)     :: AttribName  !< name of varnames to be read
INTEGER,INTENT(OUT)             :: nVal(15)    !< size of array
REAL,ALLOCATABLE,INTENT(OUT)    :: Array(:)    !< array to be read
CHARACTER(LEN=255),ALLOCATABLE,INTENT(OUT) :: VarNames(:) !< variable names
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL  :: found
INTEGER  :: dims
!===================================================================================================================================
nVal=-1
SDEALLOCATE(Array)
SDEALLOCATE(VarNames)

CALL DatasetExists(File_ID, TRIM(ArrayName), found)
IF (found) THEN
  ! get size of array
  CALL GetDataSize(File_ID,TRIM(ArrayName),dims,HSize)
  nVal(1:dims)=INT(HSize)
  IF(nVal(dims).NE.nGlobalElems) STOP 'Last array dimension != nElems !'
  nVal(dims)=nElems
  DEALLOCATE(HSize)
  ALLOCATE(Array(PRODUCT(nVal(1:dims))))
  ALLOCATE(VarNames(nVal(1)))

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVal       => INT(nVal(1:dims),IK)    ,&
        OffsetElem => INT(OffsetElem,IK) )
    ! read array
    CALL ReadArray(TRIM(ArrayName),dims,nVal,OffsetElem,dims,RealArray=Array)
  END ASSOCIATE

  ! read variable names
  CALL ReadAttribute(File_ID,TRIM(AttribName),nVal(1),StrArray=VarNames)
END IF

END SUBROUTINE GetArrayAndName


!==================================================================================================================================
!> Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!==================================================================================================================================
SUBROUTINE ReadArray(ArrayName,Rank,nVal,Offset_in,offset_dim,RealArray,IntegerArray,StrArray,IntegerArray_i4&
#if USE_MPI
  ,ReadNonOverlap_opt&
#endif /*USE_MPI*/
        )
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER                                                   :: Rank         !< number of dimensions of the array
INTEGER(KIND=IK)                                          :: offset_in    !< offset =0, start at beginning of the array
INTEGER                                                   :: offset_dim   !< which dimension is the offset (only one dim. possible)
INTEGER(KIND=IK)                                          :: nVal(Rank)   !< size of complete (local) array to write
CHARACTER(LEN=*),INTENT(IN)                               :: ArrayName    !< name of array to be read
REAL,&
    DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET  :: RealArray    !< only if real array shall be read
INTEGER(KIND=IK),&
    DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET  :: IntegerArray !< only if integer array shall be read of KIND=IK
INTEGER,&
    DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET  :: IntegerArray_i4 !< only if integer array shall be
!                                                                            !< read of KIND=4
CHARACTER(LEN=255),&
    DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET  :: StrArray        !< only if string shall be read
#if USE_MPI
LOGICAL,OPTIONAL,INTENT(IN)                               :: ReadNonOverlap_opt !< T: Set H5FD_MPIO_COLLECTIVE_F for collective I/O
!                                                                               !< F: Read individually
!                                                                               !< IMPORTANT: ReadNonOverlap_opt = T will break
!                                                                               !< when the data that is read has overlapping
!                                                                               !< elements, e.g.,
!                                                                               !<  proc 0 reads elements from index 1:10
!                                                                               !<  proc 1 reads elements from index 5:15
!                                                                               !<  ...
!                                                                               !< and will produce garbage without giving an error
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,Type_ID,MemSpace,FileSpace,PList_ID
INTEGER(HSIZE_T)               :: Offset(Rank),Dimsf(Rank)
TYPE(C_PTR)                    :: buf
#if USE_MPI
INTEGER(HID_T)                 :: driver
LOGICAL                        :: ReadNonOverlap
#endif /*USE_MPI*/
!==================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
Dimsf=nVal
LOGWRITE(*,*)'Dimsf,Offset=',Dimsf,Offset_in
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(File_ID, TRIM(ArrayName) , DSet_ID, iError)

IF(iError.NE.0) CALL Abort(__STAMP__,'Array '//TRIM(ArrayName)//' does not exist.')

! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
Offset(:)=0
Offset(offset_dim)=Offset_in
CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if USE_MPI
! Set property list to collective dataset read
!CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError) ! old
CALL H5PGET_DRIVER_F(Plist_File_ID, driver, iError) ! remove error "collective access for MPI-based drivers only"
IF(PRESENT(ReadNonOverlap_opt))THEN
  ReadNonOverlap=ReadNonOverlap_opt
ELSE
  ReadNonOverlap=.TRUE.
END IF ! PRESENT(ReadNonOverlap_opt)
IF((ReadNonOverlap).AND.(driver.EQ.H5FD_MPIO_F)) CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError)
#endif /*USE_MPI*/
CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)

! Read the data
IF(PRESENT(RealArray))       buf=C_LOC(RealArray)
IF(PRESENT(IntegerArray))    buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerArray_i4)) buf=C_LOC(IntegerArray_i4)
IF(PRESENT(StrArray))        buf=C_LOC(StrArray(1))
CALL H5DREAD_F(DSet_ID,Type_ID,buf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)

! Close the datatype, property list, dataspaces and dataset.
CALL H5TCLOSE_F(Type_ID, iError)
CALL H5PCLOSE_F(PList_ID,iError)
CALL H5SCLOSE_F(FileSpace,iError)! Close the file dataspace
CALL H5DCLOSE_F(DSet_ID, iError) ! Close the dataset
CALL H5SCLOSE_F(MemSpace,iError) ! Close the memory dataspace

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadArray



!==================================================================================================================================
!> Subroutine to read attributes from HDF5 file.
!==================================================================================================================================
SUBROUTINE ReadAttribute(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntScalar,&
                                 StrScalar,LogicalScalar,RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)                  :: Loc_ID_in         !< HDF5 file id of opened file
INTEGER           ,INTENT(IN)                  :: nVal              !< number of attributes in case an array is expected
CHARACTER(LEN=*)  ,INTENT(IN)                  :: AttribName        !< name of attribute to be read
CHARACTER(LEN=*)  ,INTENT(IN) ,OPTIONAL        :: DatasetName       !< dataset name in case attribute is located in a dataset
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealArray(nVal)   !< Array of real array attributes
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntArray(nVal)    !< Array for integer array for attributes
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealScalar        !< Scalar real attribute
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntScalar         !< Scalar integer attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrScalar         !< Scalar string attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrArray(nVal)    !< Array for character array attributes
LOGICAL           ,INTENT(OUT),OPTIONAL        :: LogicalScalar     !< Scalar logical attribute
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Attr_ID,Type_ID,Loc_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER,TARGET                 :: IntToLog
CHARACTER(LEN=255),TARGET      :: StrTmp(1)
TYPE(C_PTR)                    :: buf
!==================================================================================================================================

LOGWRITE(*,*)' READ ATTRIBUTE "',TRIM(AttribName),'" FROM HDF5 FILE...'
Dimsf(1)=nVal
Loc_ID=Loc_ID_in
IF(PRESENT(DatasetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
END IF

! Create the attribute for group Loc_ID.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)

IF(iError.NE.0) CALL abort(__STAMP__,&
    'Attribute ['//TRIM(AttribName)//'] does not exist or h5 file already opened by a differen program')

IF(PRESENT(RealArray))     RealArray=0.
IF(PRESENT(RealScalar))    RealScalar=0.
IF(PRESENT(IntArray))      IntArray=0
IF(PRESENT(IntScalar))     IntScalar=0
IF(PRESENT(LogicalScalar)) LogicalScalar=.FALSE.
IF(PRESENT(StrScalar))THEN
  StrScalar=''
  StrTmp(1)=''
END IF
IF(PRESENT(StrArray))      StrArray(:)=''

IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealScalar))    Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntArray))      Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntScalar))     Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar)) Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)

buf=C_NULL_PTR
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntArray))      buf=C_LOC(IntArray)
IF(PRESENT(IntScalar))     buf=C_LOC(IntScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(IntToLog)
IF(PRESENT(StrScalar))     buf=C_LOC(StrTmp(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))

! Read the attribute data.
CALL H5AREAD_F(Attr_ID, Type_ID, buf, iError)

IF(PRESENT(LogicalScalar)) LogicalScalar=(IntToLog.EQ.1)
IF(PRESENT(StrScalar))     StrScalar=StrTmp(1)
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)

! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadAttribute


!===================================================================================================================================
!> Subroutine to determine filename of next HDF5 file for FlushHDF5
!===================================================================================================================================
#if USE_MPI
SUBROUTINE GetHDF5NextFileName(FileName,NextFileName_HDF5,single)
#else
SUBROUTINE GetHDF5NextFileName(FileName,NextFileName_HDF5)
#endif
! MODULES
USE MOD_globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName          !< filename to check
#if USE_MPI
LOGICAL,INTENT(IN)             :: single            !< switch whether file is being accessed in parallel my MPI_COMM_PICLAS
#endif
CHARACTER(LEN=255),INTENT(OUT) :: NextFileName_HDF5 !< output: follow up file according to checked file opened
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: ReadError
INTEGER(HID_T)                 :: File_ID_loc,Plist_ID
LOGICAL                        :: Exists
!===================================================================================================================================
LOGWRITE(*,*)' GET NEXT FILE NAME FROM HDF5 FILE ', TRIM (FileName),' ...'
ReadError=0
NextFileName_HDF5=''
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Setup file access property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if USE_MPI
IF(.NOT.single)THEN
  ! Set property list to MPI IO
  CALL H5PSET_FAPL_MPIO_F(Plist_ID, MPI_COMM_PICLAS, MPI_INFO_NULL, iError)
END IF
#endif /*USE_MPI*/
! Open file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID_loc, iError,access_prp = Plist_ID)
ReadError=iError
CALL H5PCLOSE_F(Plist_ID, iError)
iError=ReadError
IF (iError.EQ.0) THEN
  ! Check if the attribute 'NextFile' exists (Create the attribute for group File_ID_loc.)
  CALL H5AEXISTS_F(File_ID_loc, 'NextFile', Exists, iError)
  IF (Exists) THEN
    ! Open the attribute "NextFile" of opened file
    CALL ReadAttribute(File_ID_loc,'NextFile',1,StrScalar=NextFileName_HDF5)
  ELSE
    NextFileName_HDF5=TRIM(FileName) ! fall back to old file name, which will be deleted anyway
  END IF ! NextFileExists
  ! Close the file.
  CALL H5FCLOSE_F(File_ID_loc, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
  iError=-1
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE GetHDF5NextFileName


END MODULE MOD_HDF5_Input
