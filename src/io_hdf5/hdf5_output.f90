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

MODULE MOD_HDF5_output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE FlushHDF5
  MODULE PROCEDURE FlushHDF5
END INTERFACE

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

INTERFACE GenerateFileSkeleton
  MODULE PROCEDURE GenerateFileSkeleton
END INTERFACE

INTERFACE GenerateNextFileInfo
  MODULE PROCEDURE GenerateNextFileInfo
END INTERFACE

INTERFACE WriteTimeAverage
  MODULE PROCEDURE WriteTimeAverage
END INTERFACE

INTERFACE
  SUBROUTINE copy_userblock(outfilename,infilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: outfilename(*)
      CHARACTER(KIND=C_CHAR) :: infilename(*)
  END SUBROUTINE copy_userblock
END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: FlushHDF5,WriteHDF5Header,GatheredWriteArray
PUBLIC :: WriteArrayToHDF5,WriteAttributeToHDF5,GenerateFileSkeleton
PUBLIC :: WriteTimeAverage,GenerateNextFileInfo, copy_userblock
#if USE_MPI && defined(PARTICLES)
PUBLIC :: DistributedWriteArray
#endif /*USE_MPI && defined(PARTICLES)*/
PUBLIC :: RemoveHDF5
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteTimeAverage(MeshFileName,OutputTime,PreviousTime,VarNamesAvg,VarNamesFluc,UAvg,UFluc,dtAvg,nVar_Avg,nVar_Fluc)
!==================================================================================================================================
!> Subroutine to write time averaged data and fluctuations HDF5 format
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars ,ONLY: ProjectName
USE MOD_Mesh_Vars    ,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVar_Avg                                     !< Number of averaged variables
INTEGER,INTENT(IN)             :: nVar_Fluc                                    !< Number of fluctuations
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName                                 !< Name of mesh file
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesAvg(nVar_Avg)                        !< Average variable names
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesFluc(nVar_Fluc)                      !< Fluctuations variable names
REAL,INTENT(IN)                :: OutputTime                                   !< Time of output
REAL,INTENT(IN),OPTIONAL       :: PreviousTime                                 !< Time of previous output
REAL,INTENT(IN),TARGET         :: UAvg(nVar_Avg,0:PP_N,0:PP_N,0:PP_N,nElems)   !< Averaged Solution
REAL,INTENT(IN),TARGET         :: UFluc(nVar_Fluc,0:PP_N,0:PP_N,0:PP_N,nElems) !< Fluctuations
REAL,INTENT(IN)                :: dtAvg                                        !< Timestep of averaging
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
REAL                           :: StartT,EndT
!==================================================================================================================================
IF((nVar_Avg.EQ.0).AND.(nVar_Fluc.EQ.0)) RETURN ! no time averaging

GETTIME(StartT)
SWRITE (UNIT_stdOut,'(A)',ADVANCE='NO') ' WRITE TIME AVERAGED STATE AND FLUCTUATIONS TO HDF5 FILE...'

! generate nextfile info in previous output file
IF(PRESENT(PreviousTime))THEN
  IF(MPIRoot .AND. PreviousTime.LT.OutputTime) CALL GenerateNextFileInfo('TimeAvg',OutputTime,PreviousTime)
END IF

! Write timeaverages ---------------------------------------------------------------------------------------------------------------
IF(nVar_Avg.GT.0)THEN
  ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_TimeAvg',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton('TimeAvg',nVar_Avg,VarNamesAvg,MeshFileName,OutputTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

  ! Reopen file and write DG solution
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        nElems          => INT(nElems,IK)          ,&
        nVar_Avg        => INT(nVar_Avg,IK)        ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName , create = .FALSE.  , &
                            DataSetName     = 'DG_Solution' , rank = 5          , &
                            nValGlobal      = (/nVar_Avg    , N+1_IK            , N+1_IK , N+1_IK , nGlobalElems/) , &
                            nVal            = (/nVar_Avg    , N+1_IK            , N+1_IK , N+1_IK , nElems/)       , &
                            offset          = (/0_IK        , 0_IK              , 0_IK   , 0_IK   , offsetElem/)   , &
                            collective      = .TRUE.        , RealArray = UAvg)
  END ASSOCIATE
END IF

! Write fluctuations ---------------------------------------------------------------------------------------------------------------
IF(nVar_Fluc.GT.0)THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Fluc',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton('Fluc',nVar_Fluc,VarNamesFluc,MeshFileName,OutputTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

  ! Reopen file and write DG solution
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        nElems          => INT(nElems,IK)          ,&
        nVar_Fluc       => INT(nVar_Fluc,IK)       ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName , create = .FALSE.                                                          , &
                            DataSetName     = 'DG_Solution' , rank = 5           , &
                            nValGlobal      = (/nVar_Fluc   , N+1_IK             , N+1_IK , N+1_IK , nGlobalElems/) , &
                            nVal            = (/nVar_Fluc   , N+1_IK             , N+1_IK , N+1_IK , nElems/)       , &
                            offset          = (/0_IK        , 0_IK               , 0_IK   , 0_IK   , offsetElem/)   , &
                            collective      = .TRUE.        , RealArray = UFluc)
  END ASSOCIATE
END IF

GETTIME(endT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteTimeAverage


SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FileNameIn,WriteUserblockIn,NIn,NodeType_in)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Output_Vars            ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
USE MOD_Interpolation_Vars     ,ONLY: NodeType
#ifdef INTEL
USE IFPORT                     ,ONLY: SYSTEM
#endif
#ifdef PARTICLES
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: UseBRElectronFluid
#endif /*USE_HDG*/
#endif /*PARTICLES*/
!USE MOD_PreProcFlags
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: TypeString
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FileNameIn
INTEGER,INTENT(IN)                   :: nVar
INTEGER,INTENT(IN),OPTIONAL          :: NIn
CHARACTER(LEN=255)                   :: StrVarNames(nVar)
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
LOGICAL,INTENT(IN),OPTIONAL          :: WriteUserblockIn
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: NodeType_in        !< Type of 1D points
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                               :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)                             :: Dimsf(5)
CHARACTER(LEN=255)                           :: FileName
#ifdef PARTICLES
CHARACTER(LEN=255), DIMENSION(1:3),PARAMETER :: TrackingString = (/'refmapping  ', 'tracing     ', 'triatracking'/)
#endif /*PARTICLES*/
LOGICAL                                      :: WriteUserblock
INTEGER                                      :: Nloc
!===================================================================================================================================
! Check if NIn is to be used
IF(PRESENT(NIn))THEN
  Nloc = NIn
ELSE
  Nloc = PP_N
END IF ! PRESENT(NIn)
! Create file
IF(PRESENT(FileNameIn))THEN
  FileName=FileNameIn
ELSE
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
END IF ! PRESENT(FileNameIn)
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)

! Write file header
CALL WriteHDF5Header(TRIM(TypeString),File_ID)

! Preallocate the data space for the dataset.
Dimsf=(/nVar,Nloc+1,Nloc+1,Nloc+1,nGlobalElems/)
CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,'DG_Solution', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=Nloc)
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
IF(PRESENT(NodeType_in))THEN
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType_in/))
ELSE
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
END IF ! PRESENT(NodeType_in)
CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)

CALL WriteAttributeToHDF5(File_ID,'NComputation',1,IntegerScalar=Nloc)

#ifdef PARTICLES
CALL WriteAttributeToHDF5(File_ID,'TrackingMethod',1,StrScalar=(/TRIM(TrackingString(TrackingMethod))/))
#if USE_HDG
IF(UseBRElectronFluid)THEN
  CALL WriteAttributeToHDF5(File_ID,'SimulationModel',1,StrScalar=(/'HDG-BR'/))
ELSE
  CALL WriteAttributeToHDF5(File_ID,'SimulationModel',1,StrScalar=(/'HDG'/))
END IF ! UseBRElectronFluid
#endif /*USE_HDG*/
#endif /*PARTICLES*/

CALL CloseDataFile()

! Add userblock to hdf5-file
IF(PRESENT(WriteUserblockIn))THEN
  WriteUserblock = WriteUserblockIn
ELSE
  WriteUserblock = .TRUE.
END IF ! PRESENT(WriteUserblockIn)
IF(WriteUserblock) CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)

END SUBROUTINE GenerateFileSkeleton


SUBROUTINE GenerateNextFileInfo(TypeString,OutputTime,PreviousTime)
!===================================================================================================================================
!> Subroutine that opens the previous written file on root processor and writes the necessary nextfile info
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
#ifdef INTEL
USE IFPORT                 ,ONLY: SYSTEM
#endif
!USE MOD_PreProcFlags
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN)                :: PreviousTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,MeshFile255
!===================================================================================================================================
! Set old file name
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),PreviousTime))//'.h5'
! The restart file name and the file name set here might differ (renamed restart file or changed project name).
! Therefore, the attribute is only written if the file exists.
IF(FILEEXISTS(Filename))THEN
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=(/MeshFile255/),Overwrite=.TRUE.)
  CALL CloseDataFile()
END IF ! FILEEXISTS(Filename)

END SUBROUTINE GenerateNextFileInfo


SUBROUTINE FlushHDF5(FlushTime_In)
!===================================================================================================================================
! Deletes all HDF5 output files, beginning from time Flushtime
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars     ,ONLY: ProjectName
USE MOD_HDF5_Input       ,ONLY: GetHDF5NextFileName
#if USE_LOADBALANCE
USE MOD_Loadbalance_Vars ,ONLY: DoLoadBalance,nLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Output_Vars      ,ONLY: DoWriteStateToHDF5
USE MOD_Restart_Vars     ,ONLY: DoRestart,FlushInitialState
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL :: FlushTime_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: stat,ioUnit
REAL                     :: FlushTime
CHARACTER(LEN=255)       :: InputFile,NextFile
!===================================================================================================================================
! Only MPI root does the flushing and only if DoWriteStateToHDF5 is true
IF((.NOT.MPIRoot).OR.(.NOT.DoWriteStateToHDF5)) RETURN

#if USE_LOADBALANCE
IF(DoLoadBalance.AND.nLoadBalance.GT.0) RETURN
#endif /*USE_LOADBALANCE*/

WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' DELETING OLD HDF5 FILES...'
IF (.NOT.PRESENT(FlushTime_In)) THEN
  FlushTime=0.0
ELSE
  FlushTime=FlushTime_In
END IF

! Delete state files
NextFile=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',FlushTime))//'.h5'

! If the original restart file is not to be deleted, skip this file and go to the next one
IF(DoRestart.AND.(.NOT.FlushInitialState))THEN
  ! Set next file name
#if USE_MPI
  CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
  CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
END IF ! .NOT.FlushInitialState

! Loop over all possible state files that can be deleted
DO
  InputFile=TRIM(NextFile)
  ! Set next file name
#if USE_MPI
  CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
  CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
  ! Delete File - only root
  stat=0
  OPEN ( NEWUNIT= ioUnit,         &
         FILE   = InputFile,      &
         STATUS = 'OLD',          &
         ACTION = 'WRITE',        &
         ACCESS = 'SEQUENTIAL',   &
         IOSTAT = stat          )
  IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
  IF(iError.NE.0) EXIT  ! iError is set in GetHDF5NextFileName !
END DO

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'DONE'

END SUBROUTINE FlushHDF5


SUBROUTINE RemoveHDF5(InputFile)
!===================================================================================================================================
! Deletes all HDF5 output files, beginning from time Flushtime
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: InputFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: stat,ioUnit
!===================================================================================================================================
! Only MPI root does the killing
IF(.NOT.MPIRoot) RETURN

WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' DELETING HDF5 FILE ['//TRIM(InputFile)//']...'

! Delete File - only root
stat=0
OPEN ( NEWUNIT= ioUnit,         &
       FILE   = TRIM(InputFile),&
       STATUS = 'OLD',          &
       ACTION = 'WRITE',        &
       ACCESS = 'SEQUENTIAL',   &
       IOSTAT = stat          )
IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
IF(iError.NE.0) WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '**** FAILED to remove ['//TRIM(InputFile)//'] with iError.NE.0 ****'

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'DONE'

END SUBROUTINE RemoveHDF5


SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: ProgramName,FileVersionReal,FileVersionInt,ProjectName,PiclasVersionStr
USE MOD_Globals_Vars ,ONLY: MajorVersion,MinorVersion,PatchVersion
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: FileType_in
INTEGER(HID_T),INTENT(IN)                :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                       :: tmp255
!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number

!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files

! First write program name
tmp255=TRIM(ProgramName)
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=(/tmp255/))
tmp255=TRIM(FileType_in)
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=(/tmp255/))
tmp255=TRIM(ProjectName)
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=(/tmp255/))
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersionReal)
WRITE(UNIT=PiclasVersionStr,FMT='(I0,A1,I0,A1,I0)') MajorVersion,".",MinorVersion,".",PatchVersion
tmp255=TRIM(PiclasVersionStr)
CALL WriteAttributeToHDF5(File_ID,'Piclas_Version',1,StrScalar=(/tmp255/))
CALL WriteAttributeToHDF5(File_ID,'Piclas_VersionInt',1,IntegerScalar=FileVersionInt)
END SUBROUTINE WriteHDF5Header


SUBROUTINE WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntegerArray,StrArray,IntegerArray_i4)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: DataSetName
INTEGER,INTENT(IN)                            :: rank                  ! < number of dimensions of the array
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)      ! < max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)            ! < size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)          ! < offset =0, start at beginning of the array
LOGICAL,INTENT(IN)                            :: collective            ! < use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL                   :: resizeDim(rank)       ! < specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL                   :: chunkSize(rank)       ! < specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(rank)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray_i4(rank) ! KIND=4
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray(rank)    ! KIND=IK
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(rank)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky
TYPE(C_PTR)                    :: buf
#if !defined(INTKIND8)
INTEGER(KIND=8)                :: Nbr8
INTEGER                        :: irank
#endif /*!defined(INTKIND8)*/
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

#if !defined(INTKIND8)
! Sanity check: Determine the total number of elements that are written to .h5 vs. maximum of INT4
IF(MPIRoot)THEN
  Nbr8 = 1
  DO irank = 1, rank
    Nbr8 = Nbr8 * INT(nValGlobal(irank),8)
  END DO ! i = 1, rank
  IF(Nbr8.GT.INT(HUGE(1_4),8))THEN
    WRITE (UNIT_stdOut,'(A,I0,A,I0,A1)',ADVANCE='NO') "WARNING: Number of entries in "//TRIM(DataSetName)//" ",Nbr8,&
        " is larger than ",HUGE(1_4)," "
  END IF ! Nbr8.GT.INT(HUGE(1_4),9)
END IF ! MPIRoot
#endif /*!defined(INTKIND8)*/

! specify chunk size if desired
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF
! make array extendable in case you want to append something
IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize)) CALL abort(__STAMP__,'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray))        Type_ID=H5T_NATIVE_DOUBLE
!IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntegerArray))     Type_ID=h5kind_to_type(IK,H5_INTEGER_KIND)
IF(PRESENT(IntegerArray_i4))  Type_ID=h5kind_to_type(4,H5_INTEGER_KIND)
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(Type_ID, SizeSet, iError)
END IF

Dimsf = nValGlobal ! we need the global array size
CALL H5ESET_AUTO_F(0,iError)
CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
IF(iError.NE.0)THEN ! does not exist
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), Type_ID, FileSpace, DSet_ID,iError,dsetparams)
  CALL H5SCLOSE_F(FileSpace, iError)
END IF
CALL H5ESET_AUTO_F(1,iError)
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
OffsetHDF = Offset
! Create the data space in the memory
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
ELSE
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SSELECT_NONE_F(FileSpace,iError)
ELSE
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, OffsetHDF, Dimsf, iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if USE_MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

!Write the dataset collectively.
IF(PRESENT(IntegerArray))    buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerArray_i4)) buf=C_LOC(IntegerArray_i4)
IF(PRESENT(RealArray))       buf=C_LOC(RealArray)
IF(PRESENT(StrArray))        buf=C_LOC(StrArray(1))
!IF(ANY(Dimsf.EQ.0)) buf =NULL()
IF(ANY(Dimsf.EQ.0)) THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,C_NULL_PTR,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
ELSE
  CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list, dataspaces and dataset.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArrayToHDF5


SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,&
                                RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                RealArray,IntegerArray,StrArray, &
                                Overwrite)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_HDF5_Input            ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)           :: Loc_ID_in
CHARACTER(LEN=*)  ,INTENT(IN)           :: AttribName
INTEGER           ,INTENT(IN)           :: nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray(nVal)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(nVal)
LOGICAL           ,INTENT(IN),OPTIONAL        :: LogicalScalar
LOGICAL           ,INTENT(IN),OPTIONAL        :: Overwrite
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: Loc_ID    ! Object identifier
INTEGER(HID_T)                 :: Type_ID   ! Attribute datatype identifier
INTEGER(HID_T)                 :: DataSpace ! Attribute dataspace identifier
INTEGER(HID_T)                 :: Attr_ID   ! Attribute identifier
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER,TARGET                 :: logtoint
TYPE(C_PTR)                    :: buf
LOGICAL                        :: AttribExists,Overwrite_loc
!===================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
IF(PRESENT(RealScalar))    Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerScalar)) Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar))THEN
  LogToInt=MERGE(1,0,LogicalScalar)
  Type_ID=H5T_NATIVE_INTEGER
END IF
IF(PRESENT(StrScalar).OR.PRESENT(StrArray))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  IF(PRESENT(StrScalar))THEN
    AttrLen=LEN(StrScalar(1))
  ELSE
    AttrLen=255
  END IF
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

! Check if attribute already exists
CALL DatasetExists(File_ID,TRIM(AttribName),AttribExists,attrib=.TRUE.)
IF(AttribExists)THEN
  IF(PRESENT(Overwrite))THEN
    Overwrite_loc = Overwrite
  ELSE
    Overwrite_loc = .FALSE.
  END IF
  IF(.NOT.Overwrite_loc) CALL abort(__STAMP__,'Attribute '//TRIM(AttribName)//' alreay exists in HDF5 File')
  ! Delete the old attribute only if it is re-writen below(otherwise the original info is lost)
  CALL H5ADELETE_F(Loc_ID, TRIM(AttribName), iError)
END IF ! AttribExists

! Create attribute
CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)

! Write the attribute data.
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntegerArray))  buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerScalar)) buf=C_LOC(IntegerScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(LogToInt)
IF(PRESENT(StrScalar))     buf=C_LOC(StrScalar(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))
CALL H5AWRITE_F(Attr_ID, Type_ID, buf, iError)

! Close datatype
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttributeToHDF5


SUBROUTINE GatheredWriteArray(FileName,create,DataSetName,rank,nValGlobal,nVal,offset,collective,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: FileName,DataSetName
LOGICAL,INTENT(IN)                            :: create,collective
INTEGER,INTENT(IN)                            :: rank
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)               ! max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)                     ! size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)                   ! offset =0, start at beginning of the array
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
REAL,              ALLOCATABLE          :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE          :: UStr(:)
INTEGER(KIND=IK),ALLOCATABLE            :: UInt(:)
INTEGER(KIND=IK)                        :: nValGather(rank),nDOFLocal
INTEGER(KIND=IK),DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
INTEGER(KIND=IK)                        :: i
!===================================================================================================================================
IF(gatheredWrite)THEN
  IF(ANY(offset(1:rank-1).NE.0)) &
    CALL abort(&
    __STAMP__&
    ,'Offset only allowed in last dimension for gathered IO.')

  ! Get last dim of each array on IO nodes
  nDOFLocal=PRODUCT(nVal)
  CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER_INT_KIND,nDOFPerNode,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_NODE,iError)

  ! Allocate big array and compute offsets of small arrs inside big
  offsetNode=0_IK
  IF(MPILocalRoot)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nLocalProcs
      offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    IF(PRESENT(RealArray))     ALLOCATE(UReal(PRODUCT(nValGather)))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( PRODUCT(nValGather)))
    IF(PRESENT(StrArray))      ALLOCATE(UStr( PRODUCT(nValGather)))
  ELSE
    IF(PRESENT(RealArray))     ALLOCATE(UReal(1))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( 1))
    IF(PRESENT(StrArray))      ALLOCATE(UStr( 1))
  ENDIF

  ! Associate construct for integer settings
  ASSOCIATE (&
        nDOFLocal    => INT(nDOFLocal)   ,&
        nDOFPerNode  => INT(nDOFPerNode) ,&
        offsetNode   => INT(offsetNode)  )
    ! Gather small arrays on IO nodes
    IF(PRESENT(RealArray)) CALL MPI_GATHERV(&
        RealArray , nDOFLocal     ,              MPI_DOUBLE_PRECISION , &
        UReal     , nDOFPerNode   , offsetNode , MPI_DOUBLE_PRECISION , &
        0         , MPI_COMM_NODE , iError)
    IF(PRESENT(IntegerArray))  CALL MPI_GATHERV(&
        IntegerArray  , nDOFLocal     ,              MPI_INTEGER_INT_KIND , &
        UInt          , nDOFPerNode   , offsetNode , MPI_INTEGER_INT_KIND , &
        0             , MPI_COMM_NODE , iError)
    !IF(PRESENT(StrArray))  CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE,&
    !                                        UReal,nDOFPerNode, offsetNode,MPI_DOUBLE,0,MPI_COMM_NODE,iError)
  END ASSOCIATE

  IF(MPILocalRoot)THEN
    ! Reopen file and write DG solution (only IO nodes)
    CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS)
    IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank                  , nValGlobal , nValGather , &
                                                 offset      , collective=collective , RealArray=UReal)
    IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank                  , nValGlobal , nValGather , &
                                                     offset      , collective=collective , IntegerArray =UInt)
    !IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nValGather,&
    !                                             offset,collective=collective,StrArr =UStr)
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(UReal)
  SDEALLOCATE(UInt)
  SDEALLOCATE(UStr)
ELSE
#endif
  CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                               offset      , collective , RealArray=RealArray)
  IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset          , collective , IntegerArray =IntegerArray)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                               offset      , collective , StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE GatheredWriteArray

#ifdef PARTICLES
#if USE_MPI
SUBROUTINE DistributedWriteArray(FileName,DataSetName,rank,nValGlobal,nVal,offset,collective,&
                                 offSetDim,communicator,RealArray,IntegerArray,StrArray,IntegerArray_i4)
!===================================================================================================================================
!> Write distributed data, that is not present in each proc of given communicator
!>   e.g. master surfaces that are not hosted by each proc
!> 1: check if every proc of given communicator has data
!> 2: if any proc has no data, split the communicator and write only with the new communicator
!> 3: else write with all procs of the given communicator
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: FileName,DataSetName
LOGICAL,INTENT(IN)                            :: collective
INTEGER,INTENT(IN)                            :: offSetDim,communicator
INTEGER,INTENT(IN)                            :: rank
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)               ! max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)                     ! size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)                   ! offset =0, start at beginning of the array
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray( PRODUCT(nVal))   ! KIND=IK
INTEGER         ,INTENT(IN),OPTIONAL,TARGET   :: IntegerArray_i4( PRODUCT(nVal))! KIND=4
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Color, OutPutCOMM,nOutPutProcs,MyOutputRank
LOGICAL                        :: DataOnProc, DoNotSplit
!===================================================================================================================================

DataOnProc=.FALSE.
! 1: check if every proc of given communicator has data
IF(nVal(offSetDim).GT.0) DataOnProc=.TRUE.
CALL MPI_ALLREDUCE(DataOnProc,DoNotSplit, 1, MPI_LOGICAL, MPI_LAND, COMMUNICATOR, IERROR)

IF(.NOT.DoNotSplit)THEN
! 2: if any proc has no data, split the communicator and write only with the new communicator
  color=MPI_UNDEFINED
  IF(DataOnProc) color=2001
  MyOutputRank=0

  CALL MPI_COMM_SPLIT(COMMUNICATOR, color, MyOutputRank, OutputCOMM,iError)
  IF(DataOnProc) THEN
    CALL MPI_COMM_SIZE(OutputCOMM, nOutPutProcs,iError)
    IF(nOutPutProcs.EQ.1)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal           , nVal , &
                                                   offset      , collective=.FALSE. , RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal                  , nVal , &
                                                   offset          , collective=.FALSE. , IntegerArray =IntegerArray)
      IF(PRESENT(IntegerArray_i4))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal                  , nVal , &
                                                   offset          , collective=.FALSE. , IntegerArray_i4 =IntegerArray_i4)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal          , nVal , &
                                                   offset      , collective=.FALSE. , StrArray =StrArray)
      CALL CloseDataFile()
    ELSE
      CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                                   offset      , collective , RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                                   offset          , collective , IntegerArray =IntegerArray)
      IF(PRESENT(IntegerArray_i4)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                                   offset          , collective , IntegerArray_i4 =IntegerArray_i4)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                                   offset      , collective , StrArray =StrArray)
      CALL CloseDataFile()
    END IF
    CALL MPI_BARRIER(OutputCOMM,IERROR)
    CALL MPI_COMM_FREE(OutputCOMM,iERROR)
  END IF
  ! MPI Barrier is required, that the other procs don't open the datafile while this procs are still writing
  CALL MPI_BARRIER(COMMUNICATOR,IERROR)
  OutputCOMM=MPI_UNDEFINED
ELSE
! 3: else write with all procs of the given communicator
  ! communicator_opt has to be the given communicator or else procs that are not in the given communicator might block the write out
  ! e.g. surface communicator contains only procs with physical surface and MPI_COMM_PICLAS contains every proc
  !      Consequently, MPI_COMM_PICLAS would block communication
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=communicator)
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                               offset      , collective , RealArray=RealArray)
  IF(PRESENT(IntegerArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset         , collective , IntegerArray =IntegerArray)
  IF(PRESENT(IntegerArray_i4)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset         , collective , IntegerArray_i4 =IntegerArray_i4)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                               offset      , collective , StrArray =StrArray)
  CALL CloseDataFile()
END IF

END SUBROUTINE DistributedWriteArray
#endif /*USE_MPI*/
#endif /*PARTICLES*/


END MODULE MOD_HDF5_output
