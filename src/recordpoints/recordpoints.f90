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

MODULE MOD_RecordPoints
!===================================================================================================================================
! Module contains the record points
! tracking of state variable at certain predefined points
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
INTERFACE InitRecordPoints
  MODULE PROCEDURE InitRecordPoints
END INTERFACE

INTERFACE RecordPoints
  MODULE PROCEDURE RecordPoints
END INTERFACE

INTERFACE WriteRPToHDF5
  MODULE PROCEDURE WriteRPToHDF5
END INTERFACE

INTERFACE FinalizeRecordPoints
  MODULE PROCEDURE FinalizeRecordPoints
END INTERFACE

PUBLIC::InitRecordPoints,RecordPoints,FinalizeRecordPoints,WriteRPToHDF5
!===================================================================================================================================
PUBLIC::DefineParametersRecordPoints

CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRecordPoints()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("RecordPoints")
CALL prms%CreateLogicalOption('RP_inUse',          "Set true to compute solution history at points defined in recordpoints file.",&
                                                   '.FALSE.')
CALL prms%CreateStringOption( 'RP_DefFile',        "File containing element-local parametric recordpoint coordinates and structure.")
CALL prms%CreateRealOption(   'RP_MaxMemory',      "Maximum memory in MiB to be used for storing recordpoint state history. "//&
                                                   "If memory is exceeded before regular IO level states are written to file.",&
                                                   '100.')
END SUBROUTINE DefineParametersRecordPoints

SUBROUTINE InitRecordPoints()
!===================================================================================================================================
! Read RP parameters from ini file and RP definitions from HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools         ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars  ,ONLY: InterpolationInitIsDone
USE MOD_RecordPoints_Vars   ,ONLY: RPDefFile,RP_inUse,RP_onProc,RecordpointsInitIsDone
USE MOD_RecordPoints_Vars   ,ONLY: RP_MaxBuffersize
USE MOD_RecordPoints_Vars   ,ONLY: nRP,nGlobalRP,lastSample,iSample,nSamples,RP_fileExists
#if USE_MPI
USE MOD_Recordpoints_Vars ,ONLY: RP_COMM
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: RP_maxMemory
INTEGER               :: maxRP
!===================================================================================================================================
! check if recordpoints are activated
RP_inUse=GETLOGICAL('RP_inUse','.FALSE.')
IF(.NOT.RP_inUse) RETURN
IF((.NOT.InterpolationInitIsDone) .OR. RecordPointsInitIsDone)THEN
   SWRITE(*,*) "InitRecordPoints not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINTS...'

nRP=0
iSample=0
nSamples=0
RPDefFile=GETSTR('RP_DefFile')                        ! Filename with RP coords
CALL ReadRPList(RPDefFile) ! RP_inUse is set to FALSE by ReadRPList if no RP is on proc.
maxRP=nGlobalRP
#if USE_MPI
  CALL InitRPCommunicator()
#endif /*USE_MPI*/

IF(RP_onProc)THEN
  RP_maxMemory=GETREAL('RP_MaxMemory','100.')         ! Max buffer (100MB)
  maxRP=nGlobalRP
#if USE_MPI
  CALL MPI_ALLREDUCE(nRP,maxRP,1,MPI_INTEGER,MPI_MAX,RP_COMM,iError)
#endif /*USE_MPI*/
  RP_MaxBufferSize = CEILING(RP_MaxMemory)*131072/(maxRP*(PP_nVar+1)) != size in bytes/(real*maxRP*nVar)
  SDEALLOCATE(lastSample)
  ALLOCATE(lastSample(0:PP_nVar,nRP))
END IF
RP_fileExists=.FALSE.

RecordPointsInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRecordPoints


#if USE_MPI
SUBROUTINE InitRPCommunicator()
!===================================================================================================================================
! Read RP parameters from ini file and RP definitions from HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc,myRPrank,RP_COMM,nRP_Procs
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color
!===================================================================================================================================

! create new RP communicator for RP output
color = MERGE(2,MPI_UNDEFINED,RP_onProc)
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,color,myRPrank,RP_COMM,iError)
IF (RP_onProc) THEN
  CALL MPI_COMM_RANK (RP_COMM,myRPrank,iError)
  CALL MPI_COMM_SIZE(RP_COMM, nRP_Procs,iError)
  IF(myRPrank.EQ.0) WRITE(*,*) 'RP COMM:',nRP_Procs,'procs'
END IF

END SUBROUTINE InitRPCommunicator
#endif /*USE_MPI*/


SUBROUTINE ReadRPList(FileString)
!===================================================================================================================================
! Read RP HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_Mesh_Vars             ,ONLY:MeshFile,nGlobalElems
USE MOD_Mesh_Vars             ,ONLY:OffsetElem
USE MOD_RecordPoints_Vars     ,ONLY:RP_onProc
USE MOD_RecordPoints_Vars     ,ONLY:OffsetRP,RP_ElemID,nRP,nGlobalRP,offsetRP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: MeshFile_RPList
INTEGER(8)                    :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP1,iRP_glob
INTEGER                       :: OffsetRPArray(2,PP_nElems)
REAL,ALLOCATABLE              :: xi_RP(:,:)
!===================================================================================================================================
IF(MPIRoot)THEN
  IF(.NOT.FILEEXISTS(FileString))  CALL abort(&
__STAMP__&
,'RPList from data file "'//TRIM(FileString)//'" does not exist',999,999.)
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

! compare mesh file names
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_RPList)
IF(TRIM(MeshFile_RPList).NE.TRIM(MeshFile)) THEN
  SWRITE(UNIT_stdOut,*) ' WARNING: MeshFileName from RPList differs from Mesh File specified in parameterfile!'
END IF

! Readin OffsetRP
CALL GetDataSize(File_ID,'OffsetRP',nDims,HSize)
nGlobalElems_RPList=HSize(2) !global number of elements
DEALLOCATE(HSize)
IF(nGlobalElems_RPList.NE.nGlobalElems) CALL abort(&
__STAMP__&
,'nGlobalElems from RPList differs from nGlobalElems from Mesh File!',999,999.)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      PP_nElems     => INT(PP_nElems,IK) ,&
      OffsetElem    => INT(OffsetElem,IK) )
  CALL ReadArray('OffsetRP',2,(/2_IK,PP_nElems/),OffsetElem,2,IntegerArray_i4=OffsetRPArray)
END ASSOCIATE

! Check if local domain contains any record points
! OffsetRP: first index: 1: offset in RP list for first RP on elem,
!                        2: offset in RP list for last RP on elem
! If these offsets are equal, no RP on elem.
nRP=OffsetRPArray(2,PP_nElems)-OffsetRPArray(1,1)
offsetRP = OffsetRPArray(1,1)
! Read in RP reference coordinates
CALL GetDataSize(File_ID,'xi_RP',nDims,HSize)
IF(HUGE(0).LT.HSize(2)) THEN
  CALL abort(&
__STAMP__&
,'Global number of record points exceeds INTEGER TYPE 4!',999,999.)
ELSE
  nGlobalRP=INT(HSize(2),4) !global number of RecordPoints
END IF
DEALLOCATE(HSize)
ALLOCATE(xi_RP(3,nRP))

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nRP      => INT(nRP,IK) ,&
      offsetRP => INT(offsetRP,IK) )
  CALL ReadArray('xi_RP',2,(/3_IK,nRP/),offsetRP,2,RealArray=xi_RP)
END ASSOCIATE

IF(nRP.LT.1) THEN
  RP_onProc=.FALSE.
ELSE
  RP_onProc=.TRUE.
  ! create mapping to elements
  ALLOCATE(RP_ElemID(nRP))
  DO iRP1=1,nRP
    iRP_glob=offsetRP+iRP1
    DO iElem=1,PP_nElems
      IF(iRP_glob .LE. OffsetRPArray(2,iElem) .AND. iRP_glob .GT. OffsetRPArray(1,iElem)) &
        RP_ElemID(iRP1)=iElem
    END DO
  END DO
END IF
CALL CloseDataFile()

IF(RP_onProc) CALL InitRPBasis(xi_RP)
DEALLOCATE(xi_RP)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' done.'
END SUBROUTINE ReadRPList



SUBROUTINE InitRPBasis(xi_RP)
!===================================================================================================================================
! Precalculate basis function values at recordpoint positions
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_RecordPoints_Vars     ,ONLY: nRP,L_xi_RP,L_eta_RP,L_zeta_RP
USE MOD_Interpolation_Vars    ,ONLY: xGP,wBary
USE MOD_Basis                 ,ONLY: LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: xi_RP(3,nRP)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iRP
!===================================================================================================================================
! build local basis for Recordpoints
ALLOCATE(L_xi_RP(0:PP_N,nRP), L_eta_RP(0:PP_N,nRP), L_zeta_RP(0:PP_N,nRP))
DO iRP=1,nRP
  CALL LagrangeInterpolationPolys(xi_RP(1,iRP),PP_N,xGP,wBary,L_xi_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(2,iRP),PP_N,xGP,wBary,L_eta_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(3,iRP),PP_N,xGP,wBary,L_zeta_RP(:,iRP))
END DO
END SUBROUTINE InitRPBasis


SUBROUTINE RecordPoints(t,Output)
!===================================================================================================================================
! Interpolate solution at time t to RecordPoint positions and fill output buffer
! The decision if an analysis is performed is done in PerformAnalysis.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars           ,ONLY: U
USE MOD_Timedisc_Vars     ,ONLY: dt
USE MOD_TimeDisc_Vars     ,ONLY: tAnalyze
USE MOD_Analyze_Vars      ,ONLY: Analyze_dt,FieldAnalyzeStep
USE MOD_RecordPoints_Vars ,ONLY: RP_Data,RP_ElemID
USE MOD_RecordPoints_Vars ,ONLY: RP_Buffersize,RP_MaxBuffersize,iSample
USE MOD_RecordPoints_Vars ,ONLY: l_xi_RP,l_eta_RP,l_zeta_RP,nRP
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars      ,ONLY: E
#endif /*PP_nVar==1*/
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: t
LOGICAL,INTENT(IN)             :: Output ! force sampling (e.g. first/last timestep)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iRP
#if USE_HDG
#if PP_nVar==1
INTEGER,PARAMETER       :: AddVar=3
#endif /*PP_nVar==1*/
#else
INTEGER,PARAMETER       :: AddVar=0
#endif /*USE_HDG*/
REAL                    :: U_RP(PP_nVar+AddVar,nRP)
REAL                    :: l_eta_zeta_RP
!-----------------------------------------------------------------------------------------------------------------------------------

! selection criterion for analysis is performed within PerformAnalysis

!IF(iter.EQ.0)THEN
!  ! Compute required buffersize from timestep and add 10% tolerance
!  RP_Buffersize = MIN(CEILING((1.05*Analyze_dt)/(dt*FieldAnalyzeStep))+1,RP_MaxBuffersize)
!  ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
!END IF
IF(.NOT.ALLOCATED(RP_Data))THEN
  ! Compute required buffersize from timestep and add 10% tolerance
  RP_Buffersize = MIN(CEILING((1.05*Analyze_dt)/(dt*FieldAnalyzeStep))+1,RP_MaxBuffersize)
  !IPWRITE(*,*) 'buffer',rp_buffersize,rp_maxbuffersize
  ALLOCATE(RP_Data(0:PP_nVar+AddVar,nRP,RP_Buffersize))
  RP_Data=0.
END IF

! evaluate state at RP
iSample=iSample+1
!IPWRITE(*,*) 'Sampling ...',iSample,size(U_RP),size(RP_Data),nRP,RP_BufferSize
U_RP=0.
DO iRP=1,nRP
  DO k=0,PP_N
    DO j=0,PP_N
      l_eta_zeta_RP=l_eta_RP(j,iRP)*l_zeta_RP(k,iRP)
      DO i=0,PP_N
#if USE_HDG
#if PP_nVar==1
        U_RP(:,iRP)=U_RP(:,iRP) + (/ U(:,i,j,k,RP_ElemID(iRP)), E(1:3,i,j,k,RP_ElemID(iRP)) /)*l_xi_RP(i,iRP)*l_eta_zeta_RP
#endif /*PP_nVar==1*/
#else
        U_RP(:,iRP)=U_RP(:,iRP) +    U(:,i,j,k,RP_ElemID(iRP))                                *l_xi_RP(i,iRP)*l_eta_zeta_RP
#endif /*USE_HDG*/
      END DO !i
    END DO !j
  END DO !k
END DO ! iRP
RP_Data(1:PP_nVar+AddVar , : , iSample)=U_RP
RP_Data(0                , : , iSample)=t

! dataset is full, write data and reset
!IF(iSample.EQ.RP_Buffersize) CALL WriteRPToHDF5(tWriteData,.FALSE.)
IF(iSample.EQ.RP_Buffersize) THEN
  SWRITE(UNIT_StdOut,*) ' BufferSize reached!'
  CALL WriteRPToHDF5(tAnalyze,.FALSE.)
END IF

! Suppress compiler warnings
RETURN
write(*,*) Output

END SUBROUTINE RecordPoints


SUBROUTINE WriteRPToHDF5(OutputTime,finalizeFile)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
! Information to time in HDF5-Format:
! file1: 0 :t1
! file2: t1:tend
! Hence, t1 and the fields are stored in both files
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE HDF5
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Equation_Vars     ,ONLY: StrVarNames
USE MOD_HDF5_Output       ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5
USE MOD_Globals_Vars      ,ONLY: ProjectName
USE MOD_Mesh_Vars         ,ONLY: MeshFile
USE MOD_Recordpoints_Vars ,ONLY: myRPrank,lastSample
USE MOD_Recordpoints_Vars ,ONLY: RPDefFile,RP_Data,iSample,nSamples
USE MOD_Recordpoints_Vars ,ONLY: offsetRP,nRP,nGlobalRP,lastSample
USE MOD_Recordpoints_Vars ,ONLY: RP_Buffersize,RP_Maxbuffersize,RP_fileExists,chunkSamples
#if USE_MPI
USE MOD_Recordpoints_Vars ,ONLY: RP_COMM,nRP_Procs
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime
LOGICAL,INTENT(IN)             :: finalizeFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileString,hilf
REAL                           :: startT,endT
#if USE_HDG
#if PP_nVar==1
INTEGER,PARAMETER       :: AddVar=3
#endif /*PP_nVar==1*/
#else
INTEGER,PARAMETER       :: AddVar=0
#endif /*USE_HDG*/
!===================================================================================================================================
WRITE(hilf,'(A)') ' WRITE RECORDPOINT DATA TO HDF5 FILE...'
SWRITE(UNIT_stdOut,'(A)')' '//TRIM(hilf)
GETTIME(startT)

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime))//'.h5'
IF(myRPrank.EQ.0)THEN
  CALL OpenDataFile(Filestring,create=.NOT.RP_fileExists,single=.TRUE.,readOnly=.FALSE.)
  IF(.NOT.RP_fileExists)THEN
    ! Create dataset attributes
    CALL WriteAttributeToHDF5(File_ID,'File_Type'  ,1,StrScalar=(/TRIM('RecordPoints_Data')/))
    CALL WriteAttributeToHDF5(File_ID,'MeshFile'   ,1,StrScalar=(/TRIM(MeshFile)/))
    CALL WriteAttributeToHDF5(File_ID,'ProjectName',1,StrScalar=(/TRIM(ProjectName)/))
    CALL WriteAttributeToHDF5(File_ID,'RPDefFile'  ,1,StrScalar=(/TRIM(RPDefFile)/))
    CALL WriteAttributeToHDF5(File_ID,'Time'       ,1,RealScalar=OutputTime)
    CALL WriteAttributeToHDF5(File_ID,'VarNames'   ,PP_nVar+AddVar,StrArray=StrVarNames)
  END IF
  CALL CloseDataFile()
END IF

#if USE_MPI
CALL MPI_BARRIER(RP_COMM,iError)
IF(nRP_Procs.EQ.1)THEN
  CALL OpenDataFile(Filestring,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
ELSE
  CALL OpenDataFile(Filestring,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=RP_COMM)
END IF
#else
CALL OpenDataFile(Filestring,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

IF(iSample.GT.0)THEN
  IF(.NOT.RP_fileExists) chunkSamples=iSample
  ! write buffer into file, we need two offset dimensions (one buffer, one processor)
  nSamples=nSamples+iSample

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        PP_nVarP1    => INT(PP_nVar+AddVar+1,IK)    ,&
        nSamples     => INT(nSamples,IK)     ,&
        nRP          => INT(nRP,IK)          ,&
        iSample      => INT(iSample,IK)      ,&
        offsetRP     => INT(offsetRP,IK)     )

#if USE_MPI
    IF(nRP_Procs.EQ.1)THEN
#endif
      CALL WriteArrayToHDF5(DataSetName = 'RP_Data'   , rank= 3       , &
                            nValGlobal  = (/PP_nVarP1 , INT(nGlobalRP , IK)                  , nSamples/) , &
                            nVal        = (/PP_nVarP1 , nRP           , iSample/)            , &
                            offset      = (/0_IK      , offsetRP      , nSamples-iSample/)   , &
                            resizeDim   = (/.FALSE.   , .FALSE.       , .TRUE./)             , &
                            chunkSize   = (/PP_nVar+1 , nGlobalRP     , chunkSamples      /) , &
                            RealArray   = RP_Data(:,:,1:iSample),&
                            collective  = .FALSE.)
#if USE_MPI
    ELSE
      CALL WriteArrayToHDF5(DataSetName = 'RP_Data'   , rank = 3      , &
                            nValGlobal  = (/PP_nVarP1 , INT(nGlobalRP , IK)                  , nSamples/) , &
                            nVal        = (/PP_nVarP1 , nRP           , iSample/)            , &
                            offset      = (/0_IK      , offsetRP      , nSamples-iSample/)   , &
                            resizeDim   = (/.FALSE.   , .FALSE.       , .TRUE./)             , &
                            chunkSize   = (/PP_nVar+1 , nGlobalRP     , chunkSamples      /) , &
                            RealArray   = RP_Data(:,:,1:iSample),&
                            collective  = .TRUE.)
    END IF
#endif
  END ASSOCIATE
  lastSample=RP_Data(:,:,iSample)
END IF
CALL CloseDataFile()
! Reset buffer
RP_Data=0.

iSample=0
RP_fileExists=.TRUE.
IF(finalizeFile)THEN
  IF(myRPrank.EQ.0)THEN
    WRITE(UNIT_stdOut,'(a,I4,a)')' RP Buffer  : ',nSamples,' samples.'
  END IF
  IF((nSamples.GT.RP_Buffersize).AND.(RP_Buffersize.LT.RP_Maxbuffersize))THEN
    ! Recompute required buffersize from timestep and add 10% tolerance
    RP_Buffersize=MIN(CEILING(1.2*nSamples),RP_MaxBuffersize)
    DEALLOCATE(RP_Data)
    ALLOCATE(RP_Data(0:PP_nVar+AddVar,nRP,RP_Buffersize))
  END IF
  RP_fileExists=.FALSE.
  ! last sample of previous file = first sample of next file
  iSample=1
  nSamples=0
  RP_Data(:,:,1)=lastSample
  SWRITE (UNIT_stdOut,'(A)') ' '
ELSE
  hilf=' '
END IF

GETTIME(endT)
CALL DisplayMessageAndTime(EndT-StartT, TRIM(hilf)//' DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteRPToHDF5


SUBROUTINE FinalizeRecordPoints()
!===================================================================================================================================
! Deallocate RP arrays
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars
USE MOD_LoadBalance_Vars, ONLY:DoLoadBalance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(DoLoadBalance)THEN
  IF(RP_onProc)THEN
    SDEALLOCATE(RP_Data)
#if USE_MPI
    CALL MPI_COMM_FREE(RP_Comm,iERROR)
#endif /*USE_MPI*/
  END IF
  nRP=0
  RP_onProc=.FALSE.
END IF
SDEALLOCATE(RP_ElemID)
SDEALLOCATE(L_xi_RP)
SDEALLOCATE(L_eta_RP)
SDEALLOCATE(L_zeta_RP)
SDEALLOCATE(LastSample)
RecordPointsInitIsDone = .FALSE.
END SUBROUTINE FinalizeRecordPoints


END MODULE MOD_RecordPoints
