#include "boltzplatz.h"

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

CONTAINS

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
USE MOD_RecordPoints_Vars   ,ONLY: RP_MaxBuffersize,RP_SamplingOffset
USE MOD_RecordPoints_Vars   ,ONLY: nRP,nGlobalRP,lastSample,iSample,nSamples
#ifdef MPI
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
#ifdef MPI
  CALL InitRPCommunicator()
#endif /*MPI*/

IF(RP_onProc)THEN
  RP_maxMemory=GETREAL('RP_MaxMemory','100.')         ! Max buffer (100MB)
  RP_SamplingOffset=GETINT('RP_SamplingOffset','1')   ! Sampling offset (iteration)
  maxRP=nGlobalRP
# ifdef MPI
  CALL MPI_ALLREDUCE(nRP,maxRP,1,MPI_INTEGER,MPI_MAX,RP_COMM,iError)
# endif /*MPI*/
  RP_MaxBufferSize = CEILING(RP_MaxMemory)*131072/(maxRP*(PP_nVar+1)) != size in bytes/(real*maxRP*nVar)
  SDEALLOCATE(lastSample)
  ALLOCATE(lastSample(0:PP_nVar,nRP))
END IF

RecordPointsInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRecordPoints


#ifdef MPI
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
INTEGER                   :: color,iProc
INTEGER                   :: noRPrank,RPrank
LOGICAL                   :: hasRP 
!===================================================================================================================================
color=MPI_UNDEFINED
IF(RP_onProc) color=2

! create ranks for RP communicator
IF(MPIRoot) THEN
  RPrank=-1
  noRPrank=-1
  myRPRank=0
  IF(RP_onProc) THEN
    RPrank=0
  ELSE 
    noRPrank=0
  END IF
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasRP,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    IF(hasRP) THEN
      RPrank=RPrank+1
      CALL MPI_SEND(RPrank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    ELSE
      noRPrank=noRPrank+1
      CALL MPI_SEND(noRPrank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    END IF
  END DO
ELSE
    CALL MPI_SEND(RP_onProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
    CALL MPI_RECV(myRPrank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
END IF

! create new RP communicator for RP output
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, myRPrank, RP_COMM,iError)
IF(RP_onProc) CALL MPI_COMM_SIZE(RP_COMM, nRP_Procs,iError)
IF(myRPrank.EQ.0 .AND. RP_onProc) WRITE(*,*) 'RP COMM:',nRP_Procs,'procs'

END SUBROUTINE InitRPCommunicator
#endif /*MPI*/


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
LOGICAL                       :: fileExists
CHARACTER(LEN=255)            :: MeshFile_RPList
INTEGER(8)                    :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP1,iRP_glob
INTEGER                       :: OffsetRPArray(2,PP_nElems)
REAL,ALLOCATABLE              :: xi_RP(:,:)
!===================================================================================================================================
IF(MPIRoot)THEN
  INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
  IF(.NOT.FileExists)  CALL abort(&
      __STAMP__, &
          'RPList from data file "'//TRIM(FileString)//'" does not exist',999,999.)
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'
! Open data file
#ifdef MPI
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
#else
CALL OpenDataFile(FileString,create=.FALSE.)
#endif

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
    __STAMP__, &
          'nGlobalElems from RPList differs from nGlobalElems from Mesh File!',999,999.)

CALL ReadArray('OffsetRP',2,(/2,PP_nElems/),OffsetElem,2,IntegerArray=OffsetRPArray)

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
      __STAMP__, &
  'Global number of record points exceeds INTEGER TYPE 4!',999,999.)
ELSE
  nGlobalRP=INT(HSize(2),4) !global number of RecordPoints
END IF
DEALLOCATE(HSize)
ALLOCATE(xi_RP(3,nRP)) 
CALL ReadArray('xi_RP',2,(/3,nRP/),offsetRP,2,RealArray=xi_RP)

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

SUBROUTINE RecordPoints(iter,t,forceSampling,Output)
!===================================================================================================================================
! Interpolate solution at time t to recordpoint positions and fill output buffer 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars          ,ONLY:U
USE MOD_Timedisc_Vars,    ONLY:dt
USE MOD_TimeDisc_Vars    ,ONLY:tAnalyze
USE MOD_Analyze_Vars     ,ONLY:Analyze_dt
USE MOD_RecordPoints_Vars,ONLY:RP_Data,RP_ElemID
USE MOD_RecordPoints_Vars,ONLY:RP_Buffersize,RP_MaxBuffersize,RP_SamplingOffset,iSample
USE MOD_RecordPoints_Vars,ONLY:l_xi_RP,l_eta_RP,l_zeta_RP,nRP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)     :: iter
REAL,INTENT(IN)                :: t
LOGICAL,INTENT(IN)             :: forceSampling,Output ! force sampling (e.g. first/last timestep)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iRP
REAL                    :: u_RP(PP_nVar,nRP)
REAL                    :: l_eta_zeta_RP 
!-----------------------------------------------------------------------------------------------------------------------------------
IF(.NOT.Output) THEN
  IF(MOD(iter,RP_SamplingOffset).NE.0 .AND. .NOT. forceSampling ) RETURN
END IF

!IF(iter.EQ.0)THEN
!  ! Compute required buffersize from timestep and add 10% tolerance
!  RP_Buffersize = MIN(CEILING((1.05*Analyze_dt)/(dt*RP_SamplingOffset))+1,RP_MaxBuffersize)
!  ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
!END IF
IF(.NOT.ALLOCATED(RP_Data))THEN
  ! Compute required buffersize from timestep and add 10% tolerance
  RP_Buffersize = MIN(CEILING((1.05*Analyze_dt)/(dt*RP_SamplingOffset))+1,RP_MaxBuffersize)
  !IPWRITE(*,*) 'buffer',rp_buffersize,rp_maxbuffersize
  ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
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
        U_RP(:,iRP)=U_RP(:,iRP) + U(:,i,j,k,RP_ElemID(iRP))*l_xi_RP(i,iRP)*l_eta_zeta_RP
      END DO !i
    END DO !j
  END DO !k
END DO ! iRP
RP_Data(1:PP_nVar,:,iSample)=U_RP
RP_Data(0,        :,iSample)=t

! dataset is full, write data and reset
!IF(iSample.EQ.RP_Buffersize) CALL WriteRPToHDF5(tWriteData,.FALSE.)
IF(iSample.EQ.RP_Buffersize) THEN
  SWRITE(UNIT_StdOut,*) ' BufferSize reached!'
  CALL WriteRPToHDF5(tAnalyze,.FALSE.)
END IF

END SUBROUTINE RecordPoints


SUBROUTINE WriteRPToHDF5(OutputTime,finalizeFile)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE HDF5
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Equation_Vars     ,ONLY: StrVarNames
USE MOD_HDF5_Output       ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5
USE MOD_Output_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars         ,ONLY: MeshFile
USE MOD_Recordpoints_Vars ,ONLY: myRPrank,lastSample
USE MOD_Recordpoints_Vars ,ONLY: RPDefFile,RP_Data,iSample,nSamples
USE MOD_Recordpoints_Vars ,ONLY: offsetRP,nRP,nGlobalRP,lastSample
USE MOD_Recordpoints_Vars ,ONLY: RP_Buffersize,RP_Maxbuffersize,RP_fileExists
#ifdef MPI
USE MOD_Recordpoints_Vars ,ONLY: RP_COMM
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime
LOGICAL,INTENT(IN)             :: finalizeFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileString
#ifdef MPI
REAL                           :: startT,endT
#endif
!===================================================================================================================================
IF(myRPrank.EQ.0) WRITE(UNIT_stdOut,'(a)')' WRITE RECORDPOINT DATA TO HDF5 FILE...'
#ifdef MPI
startT=MPI_WTIME()
#endif

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime))//'.h5'
#ifdef MPI
CALL OpenDataFile(Filestring,create=.NOT.RP_fileExists,single=.FALSE.,communicatorOpt=RP_COMM)
#else
CALL OpenDataFile(Filestring,create=.NOT.RP_fileExists)
#endif

IF(iSample.GT.0)THEN
  ! write buffer into file, we need two offset dimensions (one buffer, one processor)
  !CALL MPI_BARRIER(RP_COMM,iError)
  nSamples=nSamples+iSample
  !IPWRITE(*,*) 'blabla',iSample,nSamples,nRP,size(RP_Data),size(LastSample),RP_BufferSize
  !IPWRITE(*,*) 'values',RP_data(8,nRP,iSample)
  !CALL MPI_BARRIER(RP_COMM,iError)
  CALL WriteArrayToHDF5(DataSetName='RP_Data', rank=3,&
                        nValGlobal=(/PP_nVar+1,nGlobalRP,nSamples/),&
                        nVal=      (/PP_nVar+1,nRP      ,iSample/),&
                        offset=    (/0        ,offsetRP ,nSamples-iSample/),&
                        resizeDim= (/.FALSE.  ,.FALSE.  ,.TRUE./),&
                        chunkSize= (/PP_nVar+1,1        ,1/),&
                        RealArray=RP_Data(:,:,1:iSample),&
                        collective=.TRUE.)!, existing=RP_fileExists)
  lastSample=RP_Data(:,:,iSample)
END IF
! Reset buffer
RP_Data=0.

IF(finalizeFile)THEN
  IF(myRPrank.EQ.0) WRITE(UNIT_stdOut,'(a,I4,a)')' RP Buffer  : ',nSamples,' samples.'
  ! Create dataset attributes
  CALL WriteAttributeToHDF5(File_ID,'File_Type'  ,1,StrScalar=(/TRIM('RecordPoints_Data')/))
  CALL WriteAttributeToHDF5(File_ID,'MeshFile'   ,1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'ProjectName',1,StrScalar=(/TRIM(ProjectName)/))
  CALL WriteAttributeToHDF5(File_ID,'RPDefFile'  ,1,StrScalar=(/TRIM(RPDefFile)/))
  CALL WriteAttributeToHDF5(File_ID,'Time'       ,1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'VarNames'   ,PP_nVar,StrArray=StrVarNames)

  IF((nSamples.GT.RP_Buffersize).AND.(RP_Buffersize.LT.RP_Maxbuffersize))THEN
    ! Recompute required buffersize from timestep and add 10% tolerance
    RP_Buffersize=MIN(CEILING(1.2*nSamples),RP_MaxBuffersize)
    DEALLOCATE(RP_Data)
    ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
  END IF
  ! last sample = first sample
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
  nSamples=0
  iSample=0
  RP_fileExists=.FALSE.
#else
  RP_Data(:,:,1)=lastSample
  nSamples=0
  iSample=1
  RP_fileExists=.FALSE.
#endif
ELSE
  iSample=0
  RP_fileExists=.TRUE.
END IF
CALL CloseDataFile()

#ifdef MPI
endT=MPI_WTIME()
IF(myRPrank.EQ.0) WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE  [',EndT-StartT,'s]'
#else
IF(myRPrank.EQ.0) WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE'
#endif
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
#ifdef MPI
    CALL MPI_COMM_FREE(RP_Comm,iERROR)
#endif /*MPI*/
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
