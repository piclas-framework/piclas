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
USE MOD_ReadInTools         ,ONLY:CNTSTR,GETSTR,GETINT,GETREAL,GETREALARRAY,GETINTARRAY,GETLOGICAL
USE MOD_Interpolation_Vars  ,ONLY:InterpolationInitIsDone
USE MOD_RecordPoints_Vars   ,ONLY: RPDefFile  
USE MOD_RecordPoints_Vars   ,ONLY: RP_OutputInterval,RP_inUse,RP_onProc,RecordpointsInitIsDone,iSample
USE MOD_RecordPoints_Vars   ,ONLY: RP_SamplingOffset
USE MOD_RecordPoints_Vars   ,ONLY: getNewRPset
USE MOD_RecordPoints_Vars   ,ONLY: firstset,actualset,iSample_lastWrite 
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8) :: offset
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone) .OR. RecordPointsInitIsDone)THEN
   SWRITE(*,*) "InitRecordPoints not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINTS...'

! use RP_OutputInterval as flag for recordpoints module
RP_OutputInterval=GETINT('RP_OutputInterval','0')
RP_SamplingOffset=GETINT('RP_SamplingOffset','1')
IF (RP_OutputInterval.GT.0) THEN
  RP_inUse=.TRUE.
  ! get recordpoints file name
  RPDefFile=GETSTR('RP_DefFile')  
  CALL ReadRPList(RPDefFile)
  ! RP_inUse is set to FALSE by ReadRPList if no RP is on proc.
  IF(RP_onProc) THEN
    CALL InitRPBasis()
    ! initialize output buffer
    offset=0
    CALL getNewRPset(firstset,offset)
    actualset=>firstset
    iSample = 0
    iSample_lastWrite=0
  END IF
#ifdef MPI
  CALL InitRPCommunicator()
#endif /*MPI*/
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
!IF(myRPrank.EQ.0 .AND. .NOT.RP_onProc) WRITE(*,*) 'not RP COMM:',nRP_Procs,'procs'
!WRITE(*,*) 'myrank:',myRank,' RP_onProc',RP_onProc,' RP_COMM',RP_COMM,' color:',color,' myRPrank:', myRPrank

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
USE MOD_RecordPoints_Vars     ,ONLY:OffsetRP,xi_RP,RP_ElemID,nRP_loc,nRP_global,offsetRP_loc
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
INTEGER                       :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP1,iRP_glob
!===================================================================================================================================
IF(MPIRoot)THEN
  INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
  IF(.NOT.FileExists)  CALL abort(__STAMP__, &
          'RPList from data file "'//TRIM(FileString)//'" does not exist',999,999.)
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

! compare mesh file names
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_RPList)
IF(TRIM(MeshFile_RPList).NE.TRIM(MeshFile)) THEN
  SWRITE(UNIT_stdOut,*) ' WARNING: MeshFileName from RPList differs from Mesh File specified in parameterfile!'
END IF

! Readin OffsetRP 
CALL GetDataSize(File_ID,'OffsetRP',nDims,HSize)
nGlobalElems_RPList=HSize(2) !global number of elements
DEALLOCATE(HSize)
IF(nGlobalElems_RPList.NE.nGlobalElems) CALL abort(__STAMP__, &
          'nGlobalElems from RPList differs from nGlobalElems from Mesh File!',999,999.)

ALLOCATE(OffsetRP(2,PP_nElems)) 
CALL ReadArray('OffsetRP',2,(/2,PP_nElems/),OffsetElem,2,IntegerArray=OffsetRP)

! Check if local domain contains any record points
! OffsetRP: first index: 1: offset in RP list for first RP on elem,
!                        2: offset in RP list for last RP on elem
! If these offsets are equal, no RP on elem.
nRP_loc=OffsetRP(2,PP_nElems)-OffsetRP(1,1)
offsetRP_loc = OffsetRP(1,1)
! Read in RP reference coordinates
CALL GetDataSize(File_ID,'xi_RP',nDims,HSize)
nRP_global=HSize(2) !global number of RecordPoints
DEALLOCATE(HSize)
ALLOCATE(xi_RP(3,nRP_loc)) 
CALL ReadArray('xi_RP',2,(/3,nRP_loc/),offsetRP_loc,2,RealArray=xi_RP)

IF(nRP_loc.LT.1) THEN
  RP_onProc=.FALSE.
ELSE  
  RP_onProc=.TRUE.
  ! create mapping to elements
  ALLOCATE(RP_ElemID(nRP_loc))
  DO iRP1=1,nRP_loc
    iRP_glob=offsetRP_loc+iRP1
    DO iElem=1,PP_nElems
      IF(iRP_glob .LE. OffsetRP(2,iElem) .AND. iRP_glob .GT. OffsetRP(1,iElem)) &
        RP_ElemID(iRP1)=iElem
    END DO
  END DO
END IF
CALL CloseDataFile() 
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' done.'
END SUBROUTINE ReadRPList



SUBROUTINE InitRPBasis()
!===================================================================================================================================
! Precalculate basis function values at recordpoint positions 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_RecordPoints_Vars     ,ONLY:xi_RP,nRP_loc,L_xi_RP,L_eta_RP,L_zeta_RP
USE MOD_Interpolation_Vars    ,ONLY:xGP,wBary
USE MOD_Basis                 ,ONLY:LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iRP
!===================================================================================================================================
! build local basis for Recordpoints
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Setting up Record Point basis...'
ALLOCATE(L_xi_RP(0:PP_N,nRP_loc))
ALLOCATE(L_eta_RP(0:PP_N,nRP_loc))
ALLOCATE(L_zeta_RP(0:PP_N,nRP_loc))
DO iRP=1,nRP_loc 
  CALL LagrangeInterpolationPolys(xi_RP(1,iRP),PP_N,xGP,wBary,L_xi_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(2,iRP),PP_N,xGP,wBary,L_eta_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(3,iRP),PP_N,xGP,wBary,L_zeta_RP(:,iRP))
END DO
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' done.'
END SUBROUTINE InitRPBasis



SUBROUTINE RecordPoints(iter,t,forceSampling)
!===================================================================================================================================
! Interpolate solution at time t to recordpoint positions and fill output buffer 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RecordPoints_Vars,ONLY:RP_OutputInterval,RP_SamplingOffset,RP_ElemID,actualset,iSample
USE MOD_RecordPoints_Vars,ONLY:getNewRPset,RP_onProc
USE MOD_RecordPoints_Vars,ONLY:myRPrank
USE MOD_RecordPoints_Vars,ONLY:l_xi_RP,l_eta_RP,l_zeta_RP,nRP_loc
USE MOD_dg_vars          ,ONLY:U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8),INTENT(IN)     :: iter
REAL,INTENT(IN)                :: t
LOGICAL                        :: forceSampling ! force sampling independently from SamplingOffset
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iBuffer, i,j,k,iRP
REAL                    :: u_RP(1:PP_nVar,1:nRP_loc)
REAL                    :: l_eta_zeta_RP 
!-----------------------------------------------------------------------------------------------------------------------------------
IF(MOD(iter,RP_SamplingOffset).NE.0 .AND. .NOT. forceSampling) RETURN
iSample=iSample+1
iBuffer=iSample - actualset%Offset
IF(RP_onProc) THEN
  ! evaluate state at RP
  U_RP=0.  
  DO iRP=1,nRP_loc
    DO k=0,PP_N
      DO j=0,PP_N
        l_eta_zeta_RP=l_eta_RP(j,iRP)*l_zeta_RP(k,iRP)
        DO i=0,PP_N
          U_RP(:,iRP)=U_RP(:,iRP) + U(:,i,j,k,RP_ElemID(iRP))*l_xi_RP(i,iRP)*l_eta_zeta_RP
        END DO !i
      END DO !j
    END DO !k
  END DO ! iRP
  actualset%data(1:PP_nVar,:,iBuffer)=U_RP
  actualset%data(0,:,iBuffer)=t
END IF

! actualset is full, allocate new buffer
IF(iBuffer .EQ. RP_OutputInterval) THEN
  CALL GetNewRPset(actualset%nextset,iSample)
  actualset=>actualset%nextset
END IF
END SUBROUTINE RecordPoints



SUBROUTINE WriteRPToHDF5(iter,OutputTime)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE HDF5
USE MOD_io_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_HDF5_Output       ,ONLY: WriteAttribute,WriteToHDF5_multiD
USE MOD_Output_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars         ,ONLY: offsetElem,nGlobalElems
USE MOD_Mesh_Vars         ,ONLY: MeshFile
USE MOD_Recordpoints_Vars ,ONLY: RP_onProc,RP_COMM,myRPrank
USE MOD_Recordpoints_Vars ,ONLY: RPDefFile
USE MOD_Recordpoints_Vars ,ONLY: actualset,firstset,iSample_lastWrite,iSample
USE MOD_Recordpoints_Vars ,ONLY: offsetRP_loc,nRP_loc,nRP_global,RP_OutputInterval,RP_SamplingOffset
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8),INTENT(IN)     :: iter
REAL,INTENT(IN)                :: OutputTime
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVal
CHARACTER(LEN=255)             :: FileName,FileString,tmp255,StrVarNames(PP_nVar)
#ifdef MPI
REAL                           :: StartT,EndT
#endif
REAL                           :: RPBuffer(0:PP_nVar,nRP_loc)
INTEGER                        :: nSamples,iStart,iEnd,nBuffer,sample_offset
INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,counter,offset
LOGICAL                        :: existing
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)')' WRITE RECORDPOINT DATA TO HDF5 FILE...'
#ifdef MPI
StartT=MPI_WTIME()
#endif

IF(RP_onProc) THEN
  FileName=TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime)
  FileString=TRIM(FileName)//'.h5'
  CALL OpenDataFile(Filestring,.TRUE.,single=.FALSE.,communicatorOpt=RP_COMM) !create=.TRUE.
  
  tmp255=TRIM('RecordPoints_Data')
  CALL WriteAttribute(File_ID,'File_Type',1,StrScalar=tmp255)
  ! Create dataset attribute "MeshFile"
  tmp255=TRIM(MeshFile)
  CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=tmp255)
  ! Create dataset attribute "ProjectName"
  tmp255=TRIM(ProjectName)
  CALL WriteAttribute(File_ID,'ProjectName',1,StrScalar=tmp255)
  ! Create dataset attribute "RPDeffile"
  tmp255=TRIM(RPDefFile)
  CALL WriteAttribute(File_ID,'RPDefFile',1,StrScalar=tmp255)
  
  ! Create dataset attribute "Time"
  CALL WriteAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  ! Create dataset attribute "VarNames"
#if PP_nVar==8
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ'
StrVarNames(4)='MagneticFieldX'
StrVarNames(5)='MagneticFieldY'
StrVarNames(6)='MagneticFieldZ'
StrVarNames(7)='Phi'       
StrVarNames(8)='Psi'  
#endif        
  CALL WriteAttribute(File_ID,'VarNames',PP_nVar,StrArray=StrVarNames)
  nSamples = iSample - iSample_lastWrite  ! output size
  IF(myRPrank.EQ.0) WRITE(UNIT_stdOut,'(a,I4,a)')' RP Buffer  : ',nSamples,' samples.'
  
  ! write each buffer directly to the hdf5 file
  ! for this we need seperate write commands and two offset dimensions (one buffer, one processor)
  actualset=>firstset
  sample_offset=0
  existing=.FALSE.
  DO WHILE(ASSOCIATED(actualset))
    iEnd=MIN(sample_offset+RP_OutputInterval,nSamples)
    IF(iEnd .LE. nSamples) THEN 
      nBuffer=iEnd-sample_offset 
      dimsf(1:3)   = (/PP_nVar+1,nRP_global,nSamples/)
      counter(1:3) = (/PP_nVar+1,nRP_loc,nBuffer/)
      offset(1:3)  = (/0,offsetRP_loc,sample_offset/)
      CALL WritetoHDF5_multiD('RP_Data',                         &
                               3,                                &   ! rank: no. dimensions of array
                              dimsf,                             &   ! dimsf: global array dims
                              counter,                           &   ! counter: local array dims
                              offset,                            &   ! offset 
                              actualset%data(:,:,1:nBuffer),     &   ! RealArray
                              existing,                          &   ! does this dataset already exist?!
                              .TRUE.)                                ! Flag to determine whether this
                                                                     ! proc has data to write
      existing=.TRUE.
      ! save last sample temporarily
      IF(iEnd.EQ.nSamples) THEN
        RPBuffer(:,:)=actualset%data(:,:,nBuffer)
        EXIT  ! IMPORTANT!
      END IF
    END IF
    sample_offset=iEnd
    actualset%data=0.
    actualset%Offset=0
    actualset=>actualset%nextset
  END DO  
  CALL CloseDataFile()
  
  ! Reset buffer
  actualset=>firstset
  ! but keep last sample : first and last sample of each rp data file overlap with neighbors..
  actualset%data(:,:,1) = RPBuffer(:,:)
  iSample_lastWrite = iSample
  iSample=iSample+1 
  actualset%Offset=iSample-1
END IF

#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
EndT=MPI_WTIME()
SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE  [',EndT-StartT,'s]'
#else
SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE'
#endif
END SUBROUTINE WriteRPToHDF5




SUBROUTINE FinalizeRecordPoints()
!===================================================================================================================================
! Deallocate RP arrays 
!===================================================================================================================================
! MODULES
USE MOD_RecordPoints_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(RP_ElemID)
SDEALLOCATE(xi_RP)
SDEALLOCATE(L_xi_RP)
SDEALLOCATE(L_eta_RP)
SDEALLOCATE(L_zeta_RP)
SDEALLOCATE(u_RP)
RecordPointsInitIsDone = .FALSE.
END SUBROUTINE FinalizeRecordPoints



END MODULE MOD_RecordPoints
