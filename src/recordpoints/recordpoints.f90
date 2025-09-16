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
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
#if defined(LSERK) || USE_HDG || defined(discrete_velocity)
PUBLIC:: InitRecordPoints
PUBLIC:: RecordPoints
PUBLIC:: EvalRecordPoints
PUBLIC:: WriteRP
PUBLIC:: FinalizeRecordPoints
PUBLIC:: DefineParametersRecordPoints
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRecordPoints()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
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

!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRecordPoints()
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
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: RP_maxMemory
INTEGER               :: maxRP
!==================================================================================================================================
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
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitRecordPoints


#if USE_MPI
!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRPCommunicator()
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc,myRPrank,RP_COMM,nRP_Procs
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color
!==================================================================================================================================

!--- Split communicator from MPI_COMM_FLEXI
color = MERGE(2,MPI_UNDEFINED,RP_onProc)
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS,color,myRPrank,RP_COMM,iError)

! ignore comm if proc not on RP_COMM
IF (RP_onProc) THEN
  ! Find my rank on the RP communicator, comm size and proc name
  CALL MPI_COMM_RANK (RP_COMM,myRPrank,iError)
  CALL MPI_COMM_SIZE(RP_COMM, nRP_Procs,iError)
END IF

END SUBROUTINE InitRPCommunicator
#endif /*USE_MPI*/


!==================================================================================================================================
!> Read Recordpoint coordinates from HDF5 file
!==================================================================================================================================
SUBROUTINE ReadRPList(FileString)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_Mesh_Vars             ,ONLY:MeshFile,nGlobalElems
USE MOD_Mesh_Vars             ,ONLY:OffsetElem,nElems
USE MOD_RecordPoints_Vars     ,ONLY:RP_onProc
USE MOD_RecordPoints_Vars     ,ONLY:OffsetRP,RP_ElemID,nRP,nGlobalRP,offsetRP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileString !< name of hdf5 file for readin of recordpoints
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: MeshFile_RPList
INTEGER(8)                    :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP,iRP_glob
INTEGER,ALLOCATABLE           :: OffsetRPArray(:,:)
REAL,ALLOCATABLE              :: xi_RP(:,:)
!==================================================================================================================================
IF(MPIRoot)THEN
    IF (.NOT.FILEEXISTS(FileString)) &
      CALL Abort(__STAMP__,'RPList from data file "'//TRIM(FileString)//'" does not exist')
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

! compare mesh file names
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_RPList)
IF(TRIM(MeshFile_RPList).NE.TRIM(MeshFile)) THEN
    ! Print empty line to break the ADVANCE=NO
    SWRITE(UNIT_stdOut,'(/,A,A,A)') ' WARNING: MeshFileName ',TRIM(MeshFile_RPList), &
                                    ' from RPList differs from Mesh File specified in parameterfile!'
END IF

! Readin OffsetRP
CALL GetDataSize(File_ID,'OffsetRP',nDims,HSize)
nGlobalElems_RPList=HSize(2) !global number of elements
DEALLOCATE(HSize)

IF(nGlobalElems_RPList.NE.nGlobalElems) &
    CALL CollectiveStop(__STAMP__,'nGlobalElems from RPList differs from nGlobalElems from Mesh File!')

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      PP_nElems     => INT(PP_nElems,IK) ,&
      OffsetElem    => INT(OffsetElem,IK) )
  ALLOCATE(OffsetRPArray(2,offsetElem+1:offsetElem+nElems))
  CALL ReadArray('OffsetRP',2,(/2_IK,PP_nElems/),OffsetElem,2,IntegerArray_i4=OffsetRPArray)
END ASSOCIATE

! Check if local domain contains any record points
! OffsetRP: first index: 1: offset in RP list for first RP on elem,
!                        2: offset in RP list for last RP on elem
! If these offsets are equal, no RP on elem.
  nRP      = OffsetRPArray(2,offsetElem+nElems) - OffsetRPArray(1,offsetElem+1)
  offsetRP = OffsetRPArray(1,offsetElem+1)

! Read in RP reference coordinates
CALL GetDataSize(File_ID,'xi_RP',nDims,HSize)
  CHECKSAFEINT(HSize(2),4)
  nGlobalRP=INT(HSize(2),4) !global number of RecordPoints
DEALLOCATE(HSize)
ALLOCATE(xi_RP(3,nRP))

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nRP      => INT(nRP,IK) ,&
      offsetRP => INT(offsetRP,IK) )
  CALL ReadArray('xi_RP',2,(/3_IK,nRP/),offsetRP,2,RealArray=xi_RP)
END ASSOCIATE

RP_onProc = MERGE(.TRUE.,.FALSE.,nRP.GT.0)

IF (RP_onProc) THEN
  ! create mapping to elements
  ALLOCATE(RP_ElemID(nRP))
  DO iRP = 1,nRP
    iRP_glob = offsetRP+iRP
    DO iElem = 1,nElems
      IF (iRP_glob.LE.OffsetRPArray(2,offsetElem+iElem) .AND. iRP_glob.GT.OffsetRPArray(1,offsetElem+iElem)) &
        RP_ElemID(iRP) = iElem
    END DO
  END DO
END IF
CALL CloseDataFile()

IF(RP_onProc) CALL InitRPBasis(xi_RP)
DEALLOCATE(xi_RP)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' done.'
END SUBROUTINE ReadRPList


!==================================================================================================================================
!> Precompute Lagrange basis function values at recordpoints
!==================================================================================================================================
SUBROUTINE InitRPBasis(xi_RP)
! MODULES
USE MOD_PreProc
USE MOD_RecordPoints_Vars     ,ONLY: nRP,L_xi_RP,L_eta_RP,L_zeta_RP,RP_ElemID
USE MOD_Interpolation_Vars    ,ONLY: N_Inter,NMax
USE MOD_Basis                 ,ONLY: LagrangeInterpolationPolys
#if !defined(discrete_velocity)
USE MOD_DG_Vars               ,ONLY: N_DG_Mapping
#endif /*!defined(discrete_velocity)*/
USE MOD_Mesh_Vars             ,ONLY: offSetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: xi_RP(3,nRP)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iRP,Nloc,ElemID
!==================================================================================================================================
! build local basis for Recordpoints. Allocate with NMax and fill with zeros where Nloc<NMax
ALLOCATE(L_xi_RP(0:NMax,nRP), L_eta_RP(0:NMax,nRP), L_zeta_RP(0:NMax,nRP))
L_xi_RP   = 0.0
L_eta_RP  = 0.0
L_zeta_RP = 0.0
DO iRP=1,nRP
  ElemID = RP_ElemID(iRP)
#if !defined(discrete_velocity)
  Nloc   = N_DG_Mapping(2,ElemID+offSetElem)
#else
  Nloc = PP_N
#endif /*!defined(discrete_velocity)*/
  CALL LagrangeInterpolationPolys(xi_RP(1,iRP),Nloc,N_Inter(Nloc)%xGP,N_Inter(Nloc)%wBary,L_xi_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(2,iRP),Nloc,N_Inter(Nloc)%xGP,N_Inter(Nloc)%wBary,L_eta_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(3,iRP),Nloc,N_Inter(Nloc)%xGP,N_Inter(Nloc)%wBary,L_zeta_RP(:,iRP))
END DO

END SUBROUTINE InitRPBasis


!==================================================================================================================================
!> Evaluate solution at current time t at recordpoint positions and fill output buffer
!==================================================================================================================================
SUBROUTINE RecordPoints(t,forceSampling)
!SUBROUTINE RecordPoints(nVar,StrVarNames,iter,t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_TimeDisc_Vars    ,ONLY: iter
USE MOD_Analyze_Vars     ,ONLY: Analyze_dt,FieldAnalyzeStep
USE MOD_RecordPoints_Vars,ONLY: RP_Data
USE MOD_RecordPoints_Vars,ONLY: RP_Buffersize,RP_MaxBufferSize,RP_SamplingOffset,iSample
USE MOD_RecordPoints_Vars,ONLY: nRP
USE MOD_Timedisc_Vars    ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: t                       !< current time t
LOGICAL,INTENT(IN)             :: forceSampling           !< force sampling (e.g. at first/last timestep of computation)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_HDG && !(USE_FV)
#if PP_nVar==1
INTEGER,PARAMETER       :: AddVar=3
#endif /*PP_nVar==1*/
#else
INTEGER,PARAMETER       :: AddVar=0
#endif /*USE_HDG*/
INTEGER                 :: nVar
REAL, ALLOCATABLE       :: U_RP(:,:)
!----------------------------------------------------------------------------------------------------------------------------------
RP_SamplingOffset = 1 ! Set to one for now
IF(MOD(iter,INT(RP_SamplingOffset,8)).NE.0 .AND. .NOT. forceSampling) RETURN
#ifdef discrete_velocity
nVar = PP_nVar_FV
#else
nVar = PP_nVar+AddVar
#endif /*discrete_velocity*/
ALLOCATE(U_RP(nVar,nRP))
IF(.NOT.ALLOCATED(RP_Data))THEN
  ! Compute required buffersize from timestep and add 20% tolerance
  ! +1 is added to ensure a minimum buffersize of 2
  RP_Buffersize = MIN(CEILING((1.2*Analyze_dt)/(dt*FieldAnalyzeStep))+1,RP_MaxBufferSize)
  ALLOCATE(RP_Data(0:nVar,nRP,RP_Buffersize))
END IF

! evaluate state at RPs
CALL EvalRecordPoints(t,nVar,U_RP)

! Increment counter and fill buffer
iSample=iSample+1
RP_Data(1:nVar,:,iSample) = U_RP
RP_Data(0,     :,iSample) = t

! dataset is full, write data and reset
IF(iSample.EQ.RP_Buffersize) CALL WriteRP(Analyze_dt,.FALSE.)

END SUBROUTINE RecordPoints


!==================================================================================================================================
!> Evaluate solution at current time t at recordpoint positions
!==================================================================================================================================
SUBROUTINE EvalRecordPoints(t,nVar,U_RP)
! MODULES
USE MOD_Globals
USE MOD_Preproc
#ifdef discrete_velocity
USE MOD_FV_Vars                ,ONLY: U_FV
USE MOD_Timedisc_Vars          ,ONLY: dt
USE MOD_Equation_Vars_FV       ,ONLY: DVMMethod, DVMnSpecies, DVMSpecData, DVMColl, DVMnMacro
USE MOD_DistFunc               ,ONLY: MacroValuesFromDistribution, TargetDistribution, MoleculeRelaxEnergy
USE MOD_Mesh_Vars_FV           ,ONLY: Elem_xGP_FV
#else
USE MOD_DG_Vars           ,ONLY: U_N,N_DG_Mapping
#endif /*DVM*/
USE MOD_RecordPoints_Vars ,ONLY: RP_ElemID
USE MOD_RecordPoints_Vars ,ONLY: L_xi_RP,L_eta_RP,L_zeta_RP,nRP
USE MOD_Mesh_Vars         ,ONLY: offSetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)         :: t                       !< current time t
INTEGER,INTENT(IN)      :: nVar
REAL,INTENT(INOUT)      :: U_RP(nVar,nRP)          !< State at recordpoints
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iRP,ElemID,Nloc
REAL                    :: L_eta_zeta_RP
#ifdef discrete_velocity
REAL                    :: tau, prefac, MacroVal(DVMnMacro,DVMnSpecies+1), fTarget(PP_nVar_FV), rho, Pr
INTEGER                 :: iSpec, vFirstID, vLastID
REAL                    :: Erot(DVMnSpecies+1), ErelaxTrans, ErelaxRot(DVMnSpecies)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
U_RP=0.
DO iRP=1,nRP
  ElemID = RP_ElemID(iRP)
#ifndef discrete_velocity
  Nloc = N_DG_Mapping(2,ElemID+offSetElem)
  DO k = 0,Nloc
    DO j = 0,Nloc
      L_eta_zeta_RP=L_eta_RP(j,iRP)*L_zeta_RP(k,iRP)
      DO i=0,Nloc
#endif /*discrete_velocity*/
#ifdef discrete_velocity
        IF (t.GT.0..AND.DVMColl.AND.DVMMethod.GT.0) THEN
          CALL MacroValuesFromDistribution(MacroVal,U_FV(:,RP_ElemID(iRP)),dt,tau,1,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot)
          CALL MoleculeRelaxEnergy(ErelaxTrans,ErelaxRot,MacroVal(5,DVMnSpecies+1),Erot(1:DVMnSpecies),Pr)
          SELECT CASE(DVMMethod)
            CASE(1)
              prefac = tau*(1.-EXP(-dt/tau))/dt
            CASE(2)
              prefac = 2.*tau/(2.*tau+dt)
          END SELECT
          vFirstID = 1
          vLastID = 0
          DO iSpec=1,DVMnSpecies
            vLastID = vLastID + DVMSpecData(iSpec)%nVar
            CALL TargetDistribution(MacroVal(:,DVMnSpecies+1),fTarget(vFirstID:vLastID),iSpec,MacroVal(1,iSpec),rho,Pr,ErelaxTrans,ErelaxRot(iSpec))
            vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
          END DO
          U_RP(:,iRP)=U_FV(:,RP_ElemID(iRP))*prefac + fTarget(:)*(1.-prefac)
        ELSE
          U_RP(:,iRP)=U_FV(:,RP_ElemID(iRP))
          U_RP(1,iRP)=Elem_xGP_FV(1,RP_ElemID(iRP)) !when RP position is needed
        END IF
#else /*NOT discrete_velocity*/
#if USE_HDG
#if (PP_nVar==1) && !(USE_FV)
        U_RP(:,iRP)=U_RP(:,iRP) + (/ U_N(ElemID)%U(:,i,j,k), U_N(ElemID)%E(1:3,i,j,k) /)*L_xi_RP(i,iRP)*L_eta_zeta_RP
#endif /*PP_nVar==1*/
#else /*NOT USE_HDG*/
        U_RP(:,iRP)=U_RP(:,iRP) +    U_N(ElemID)%U(:,i,j,k)                             *L_xi_RP(i,iRP)*L_eta_zeta_RP
#endif /*USE_HDG*/
#endif /*discrete_velocity*/
#ifndef discrete_velocity
      END DO !i
    END DO !j
  END DO !k
#endif /*discrete_velocity*/
END DO ! iRP

! Surpress compiler warning
RETURN
L_eta_zeta_RP = t
END SUBROUTINE EvalRecordPoints


!==================================================================================================================================
!> Writes the time history of the solution at the recordpoints to an HDF5 file
!==================================================================================================================================
SUBROUTINE WriteRP(OutputTime,finalizeFile)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE HDF5
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
#if USE_FV
#if USE_HDG
USE MOD_Equation_Vars     ,ONLY: StrVarNames
#endif
USE MOD_Equation_Vars_FV       ,ONLY: StrVarNames_FV
#else
USE MOD_Equation_Vars          ,ONLY: StrVarNames
#endif
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
#ifdef discrete_velocity
USE MOD_Equation_Vars_FV  ,ONLY: DVMSpecData, DVMnSpecies
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime
LOGICAL,INTENT(IN)             :: finalizeFile
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileString,hilf,hilf2
CHARACTER(LEN=255)             :: tmp255
REAL                           :: startT,endT
#if USE_HDG
#if PP_nVar==1
INTEGER,PARAMETER       :: AddVar=3
#endif /*PP_nVar==1*/
#else
INTEGER,PARAMETER       :: AddVar=0
#endif /*USE_HDG*/

#ifdef discrete_velocity
! REAL,DIMENSION(PP_nVar_FV) :: Weights, VeloX, VeloY, VeloZ
! INTEGER                    :: iVel, jVel, kVel, upos
INTEGER                      :: iSpec
#endif
!===================================================================================================================================
WRITE(hilf,'(A)') ' WRITE RECORDPOINT DATA TO HDF5 FILE...'
SWRITE(UNIT_stdOut,'(A)')' '//TRIM(hilf)
GETTIME(startT)

#if USE_FV
#ifdef discrete_velocity
ASSOCIATE (PP_nVar_loc     => PP_nVar_FV)
#elif USE_HDG
ASSOCIATE (StrVarNames_loc => StrVarNames, &
  PP_nVar_loc     => PP_nVar+AddVar)
#else
ASSOCIATE (StrVarNames_loc => StrVarNames_FV, &
           PP_nVar_loc     => PP_nVar_FV)
#endif /*discrete_velocity*/
#else
ASSOCIATE (StrVarNames_loc => StrVarNames, &
           PP_nVar_loc     => PP_nVar+AddVar)
#endif

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime))//'.h5'

! init file or just update time
#if USE_MPI
IF(myRPrank.EQ.0)THEN
#endif /* USE_MPI */
  CALL OpenDataFile(Filestring,create=.NOT.RP_fileExists,single=.TRUE.,readOnly=.FALSE.)
  IF(.NOT.RP_fileExists)THEN
    ! Create dataset attributes
    CALL WriteAttributeToHDF5(File_ID,'File_Type'  ,1,StrScalar=(/CHARACTER(LEN=255)::'RecordPoints_Data'/))
    tmp255=TRIM(MeshFile)
    CALL WriteAttributeToHDF5(File_ID,'MeshFile'   ,1,StrScalar=(/tmp255/))
    tmp255=TRIM(ProjectName)
    CALL WriteAttributeToHDF5(File_ID,'ProjectName',1,StrScalar=(/tmp255/))
    tmp255=TRIM(RPDefFile)
    CALL WriteAttributeToHDF5(File_ID,'RPDefFile'  ,1,StrScalar=(/tmp255/))
    CALL WriteAttributeToHDF5(File_ID,'Time'       ,1,RealScalar=OutputTime)
#ifndef discrete_velocity
    CALL WriteAttributeToHDF5(File_ID,'VarNames'   ,PP_nVar_loc,StrArray=StrVarNames_loc)
#endif
  END IF
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif /* USE_MPI */


#if USE_MPI
CALL MPI_BARRIER(RP_COMM,iError)
IF(nRP_Procs.EQ.1)THEN
  CALL OpenDataFile(Filestring,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
ELSE
  CALL OpenDataFile(Filestring,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=RP_COMM)
END IF
#else
CALL OpenDataFile(Filestring,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif /* USE_MPI */

#ifdef discrete_velocity
  ! Associate construct for integer KIND=8 possibility
DO iSpec = 1,DVMnSpecies
  ASSOCIATE (nVel    => INT(DVMSpecData(iSpec)%nVelos(1),IK))
    WRITE(UNIT=hilf2,FMT='(I0.3)') iSpec

    CALL WriteArrayToHDF5(DataSetName = 'Weights'//TRIM(hilf2), rank= 2, &
    nValGlobal  = (/nVel,3_IK/)     , &
    nVal        = (/nVel,3_IK/)     , &
    offset      = (/0_IK,0_IK/)          , &
    RealArray   = DVMSpecData(iSpec)%Weights      , &
    collective  = .FALSE.)

    CALL WriteArrayToHDF5(DataSetName = 'Velos'//TRIM(hilf2), rank= 2, &
    nValGlobal  = (/nVel,3_IK/)     , &
    nVal        = (/nVel,3_IK/)     , &
    offset      = (/0_IK,0_IK/)          , &
    RealArray   = DVMSpecData(iSpec)%Velos        , &
    collective  = .FALSE.)

  END ASSOCIATE
END DO
#endif /*DVM*/

IF(iSample.GT.0)THEN
  IF(.NOT.RP_fileExists) chunkSamples=iSample
  ! write buffer into file, we need two offset dimensions (one buffer, one processor)
  nSamples=nSamples+iSample

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        PP_nVarP1    => INT(PP_nVar_loc+1,IK)    ,&
        nSamples     => INT(nSamples,IK)     ,&
        nRP          => INT(nRP,IK)          ,&
        iSample      => INT(iSample,IK)      ,&
        offsetRP     => INT(offsetRP,IK)     )

#if USE_MPI
    IF(nRP_Procs.EQ.1)THEN
#endif
      CALL WriteArrayToHDF5(DataSetName ='RP_Data'                                                         ,&
                            rank        = 3                                                               , &
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
    RP_Buffersize=MIN(CEILING(1.1*nSamples)+1,RP_MaxBufferSize)
    DEALLOCATE(RP_Data)
    ALLOCATE(RP_Data(0:PP_nVar_loc,nRP,RP_Buffersize))
  END IF

  RP_fileExists=.FALSE.
  iSample=1
  nSamples=0
  ! last sample of previous file is first sample of next file
  RP_Data(:,:,1)=lastSample
  SWRITE (UNIT_stdOut,'(A)') ' '
ELSE
  hilf=' '
END IF
END ASSOCIATE

GETTIME(endT)
CALL DisplayMessageAndTime(EndT-StartT, TRIM(hilf)//' DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteRP


!==================================================================================================================================
!> Deallocate recordpoint arrays
!==================================================================================================================================
SUBROUTINE FinalizeRecordPoints()
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars
USE MOD_LoadBalance_Vars, ONLY:DoLoadBalance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
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
SDEALLOCATE(lastSample)

#if USE_MPI
! Free MPI communicator
IF(RP_COMM.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(RP_COMM, iError)
#endif /* USE_MPI */

RecordPointsInitIsDone = .FALSE.

END SUBROUTINE FinalizeRecordPoints
#endif /*defined(LSERK) || USE_HDG || defined(discrete_velocity)*/

END MODULE MOD_RecordPoints